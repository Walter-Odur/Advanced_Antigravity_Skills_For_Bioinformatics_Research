"""REINVENT 4 setup, config generation, scoring profiles and analysis.

Generated configs are parsed with ``tomllib`` and checked against the
schema taken from the REINVENT 4 checkout itself
(``configs/SCORING.md`` and ``configs/scoring_components_example.toml``),
so a config cannot drift away from what the installed version accepts.

The analysis tests run on a CSV in REINVENT's real output format:
``Agent, Prior, Target, Score, SMILES, SMILES_state`` plus one column per
scoring component, with no step column, which is what
``write_summary()`` actually emits.
"""

from __future__ import annotations

import tomllib

import pytest

from agskills.errors import InvalidInputError, UsageError
from agskills.reinvent.analyze import analyze_results
from agskills.reinvent.config import RUN_MODES, build_config, toml_value
from agskills.reinvent.priors import (
    GENERATORS,
    check_setup,
    find_reinvent_dir,
    resolve_prior,
)
from agskills.reinvent.runner import prepare_seeds, preflight
from agskills.reinvent.scoring import (
    AGGREGATION_TYPES,
    SCORING_PROFILES,
    build_scoring_section,
    parse_component_spec,
)

from panel import EGFR_ACTIVES


def config_of(**kwargs) -> dict:
    """Build a config and parse it, so every test proves it is valid TOML."""
    kwargs.setdefault("device", "cpu")
    return tomllib.loads(build_config(**kwargs).toml)


# ---------------------------------------------------------------------------
# TOML rendering
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("value,expected", [
    (True, "true"), (False, "false"), (42, "42"), (0.0001, "0.0001"),
    (128, "128"), ("priors/reinvent.prior", '"priors/reinvent.prior"'),
    ([1, 2], "[1, 2]"), ([], "[]"),
])
def test_toml_values_render_correctly(value, expected):
    assert toml_value(value) == expected


def test_windows_paths_are_escaped():
    """A backslash in a path must not become an escape sequence."""
    rendered = toml_value(r"E:\models\reinvent.prior")
    assert tomllib.loads(f"path = {rendered}")["path"] == \
        r"E:\models\reinvent.prior"


def test_floats_avoid_scientific_notation():
    """A user has to read and edit these files."""
    assert "e-" not in toml_value(0.0001)


def test_unsupported_type_is_refused():
    with pytest.raises(UsageError):
        toml_value({"a": 1})


# ---------------------------------------------------------------------------
# Config generation: every mode produces valid TOML
# ---------------------------------------------------------------------------


def test_staged_learning_matches_the_reinvent_schema():
    parsed = config_of(mode="staged_learning", generator="reinvent",
                       scoring_profile="drug-like", num_steps=300,
                       batch_size=128)
    assert parsed["run_type"] == "staged_learning"
    assert parsed["device"] == "cpu"
    # Stages are an array of tables, which is what REINVENT requires.
    assert isinstance(parsed["stage"], list)
    assert len(parsed["stage"]) == 1
    stage = parsed["stage"][0]
    assert stage["max_steps"] == 300
    assert parsed["parameters"]["batch_size"] == 128
    assert parsed["learning_strategy"]["type"] == "dap"
    assert parsed["learning_strategy"]["sigma"] == 128
    # The agent starts from the prior.
    assert parsed["parameters"]["agent_file"] == \
        parsed["parameters"]["prior_file"]


def test_scoring_section_nests_as_reinvent_expects():
    """``[[stage.scoring.component.<Name>.endpoint]]`` is the required shape."""
    parsed = config_of(mode="staged_learning", generator="reinvent",
                       scoring_profile="drug-like")
    scoring = parsed["stage"][0]["scoring"]
    assert scoring["type"] in AGGREGATION_TYPES
    components = scoring["component"]
    assert isinstance(components, list) and components

    for entry in components:
        assert len(entry) == 1, "one component per [[component]] table"
        name = next(iter(entry))
        endpoints = entry[name]["endpoint"]
        assert isinstance(endpoints, list) and endpoints
        endpoint = endpoints[0]
        assert "name" in endpoint and "weight" in endpoint


def test_component_names_are_the_ones_reinvent_registers():
    """Names come from the checkout's own documentation, not from memory."""
    parsed = config_of(mode="staged_learning", generator="reinvent",
                       scoring_profile="drug-like")
    names = {next(iter(entry))
             for entry in parsed["stage"][0]["scoring"]["component"]}
    # Verified against configs/SCORING.md and
    # configs/scoring_components_example.toml in REINVENT 4.
    known = {"QED", "MolecularWeight", "SlogP", "TPSA", "NumRotBond",
             "HBondDonors", "HBondAcceptors", "SAScore", "custom_alerts",
             "NumAromaticRings", "GraphLength", "Csp3", "NumHeavyAtoms"}
    assert names <= known, f"unrecognised component name(s): {names - known}"


def test_transfer_learning_config():
    parsed = config_of(mode="transfer_learning", generator="reinvent",
                       smiles_file="actives.smi",
                       validation_smiles_file="val.smi", num_epochs=50)
    parameters = parsed["parameters"]
    assert parsed["run_type"] == "transfer_learning"
    assert parameters["num_epochs"] == 50
    assert parameters["smiles_file"] == "actives.smi"
    assert parameters["validation_smiles_file"] == "val.smi"
    assert parameters["output_model_file"].endswith(".model")
    # num_refs is 0 so a large training set does not crawl.
    assert parameters["num_refs"] == 0


def test_mol2mol_transfer_learning_adds_pair_thresholds():
    """Mol2Mol learns from molecule pairs and needs similarity bounds."""
    parsed = config_of(mode="transfer_learning", generator="mol2mol",
                       prior="scaffold_generic", smiles_file="actives.smi",
                       num_epochs=10)
    pairs = parsed["parameters"]["pairs"]
    assert pairs["type"] == "tanimoto"
    assert pairs["lower_threshold"] < pairs["upper_threshold"]


def test_sampling_config():
    parsed = config_of(mode="sampling", generator="reinvent",
                       num_smiles=500)
    assert parsed["run_type"] == "sampling"
    assert parsed["parameters"]["num_smiles"] == 500
    assert parsed["parameters"]["unique_molecules"] is True


def test_scoring_config():
    parsed = config_of(mode="scoring", generator="reinvent",
                       smiles_file="compounds.smi",
                       scoring_profile="anti-tb")
    assert parsed["run_type"] == "scoring"
    assert parsed["parameters"]["smiles_file"] == "compounds.smi"
    assert parsed["scoring"]["component"]


@pytest.mark.parametrize("mode", sorted(RUN_MODES))
def test_every_run_mode_produces_parseable_toml(mode):
    extra = {
        "staged_learning": {"scoring_profile": "drug-like"},
        "transfer_learning": {"smiles_file": "a.smi"},
        "sampling": {},
        "scoring": {"smiles_file": "a.smi", "scoring_profile": "drug-like"},
    }[mode]
    parsed = config_of(mode=mode, generator="reinvent", **extra)
    assert parsed["run_type"] == mode


@pytest.mark.parametrize("generator", sorted(GENERATORS))
def test_every_generator_produces_parseable_toml(generator):
    kwargs = {"mode": "staged_learning", "generator": generator,
              "scoring_profile": "drug-like"}
    if GENERATORS[generator]["needs_seeds"]:
        kwargs["smiles_file"] = "seeds.smi"
    parsed = config_of(**kwargs)
    assert parsed["run_type"] == "staged_learning"


# ---------------------------------------------------------------------------
# Config guardrails
# ---------------------------------------------------------------------------


def test_a_seeded_generator_without_seeds_is_refused():
    """LibInvent decorates a scaffold; without one there is nothing to do."""
    for generator in ("libinvent", "linkinvent", "mol2mol", "pepinvent"):
        with pytest.raises(UsageError) as excinfo:
            build_config(mode="staged_learning", generator=generator,
                         scoring_profile="drug-like", device="cpu")
        message = str(excinfo.value)
        assert "seed" in message.lower()
        assert "prepare-seeds" in message


def test_the_unseeded_generator_needs_no_seeds():
    parsed = config_of(mode="staged_learning", generator="reinvent",
                       scoring_profile="drug-like")
    assert "smiles_file" not in parsed["parameters"]


def test_optimisation_without_an_objective_is_refused():
    with pytest.raises(UsageError) as excinfo:
        build_config(mode="staged_learning", generator="reinvent",
                     device="cpu")
    assert "--scoring-profile" in str(excinfo.value)


def test_transfer_learning_without_training_data_is_refused():
    with pytest.raises(UsageError) as excinfo:
        build_config(mode="transfer_learning", generator="reinvent",
                     device="cpu")
    assert "--smiles-file" in str(excinfo.value)


def test_unknown_mode_and_generator_list_the_options():
    with pytest.raises(UsageError) as excinfo:
        build_config(mode="bogus", generator="reinvent", device="cpu")
    assert "staged_learning" in str(excinfo.value)

    with pytest.raises(UsageError) as excinfo:
        build_config(mode="sampling", generator="bogus", device="cpu")
    assert "reinvent" in str(excinfo.value)


def test_missing_prior_is_warned_about_not_silently_accepted():
    result = build_config(mode="sampling", generator="reinvent", device="cpu")
    assert any("not found" in warning for warning in result.warnings)
    assert any("check-setup" in warning for warning in result.warnings)


def test_cpu_with_many_steps_warns_about_runtime():
    result = build_config(mode="staged_learning", generator="reinvent",
                          scoring_profile="drug-like", num_steps=500,
                          device="cpu")
    assert any("HPC" in w or "GPU" in w for w in result.warnings)


def test_transfer_learning_without_validation_warns_about_overfitting():
    result = build_config(mode="transfer_learning", generator="reinvent",
                          smiles_file="a.smi", device="cpu")
    assert any("over-fitting" in w for w in result.warnings)


def test_a_profile_with_a_caveat_surfaces_it():
    """The anti-TB profile is physicochemical only, and must say so."""
    result = build_config(mode="staged_learning", generator="reinvent",
                          scoring_profile="anti-tb", device="cpu")
    assert any("not an activity model" in w for w in result.warnings)


def test_diversity_filter_is_on_by_default():
    """It is what stops the agent collapsing onto a single scaffold."""
    parsed = config_of(mode="staged_learning", generator="reinvent",
                       scoring_profile="drug-like")
    assert parsed["diversity_filter"]["type"] == "IdenticalMurckoScaffold"

    without = config_of(mode="staged_learning", generator="reinvent",
                        scoring_profile="drug-like", diversity_filter=False)
    assert "diversity_filter" not in without


def test_config_metadata_reports_expected_outputs_and_next_steps():
    result = build_config(mode="staged_learning", generator="reinvent",
                          scoring_profile="drug-like", num_steps=100,
                          batch_size=64, device="cpu", csv_prefix="run")
    assert "run_1.csv" in result.expected_outputs
    assert result.next_steps
    assert any("analyze-results" in step for step in result.next_steps)
    # The molecule budget is stated, so the cost is visible up front.
    assert any("6,400" in step for step in result.next_steps)


# ---------------------------------------------------------------------------
# Scoring profiles
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("profile", sorted(SCORING_PROFILES))
def test_every_profile_is_well_formed(profile):
    section = build_scoring_section(profile)
    assert section["n_components"] >= 3
    assert section["type"] in AGGREGATION_TYPES
    for entry in section["components"]:
        assert entry["component"]
        endpoint = entry["endpoint"]
        assert endpoint["name"]
        assert endpoint["weight"] > 0


@pytest.mark.parametrize("profile", sorted(SCORING_PROFILES))
def test_numeric_components_carry_a_transform(profile):
    """A raw property has to be mapped into [0, 1] to be aggregated.

    Without a transform the aggregate is meaningless, and a hard step would
    give the optimiser no gradient to follow.
    """
    plain_in_zero_one = {"QED", "custom_alerts"}
    for entry in build_scoring_section(profile)["components"]:
        if entry["component"] in plain_in_zero_one:
            continue
        transform = entry["endpoint"].get("transform")
        assert transform, f"{profile}/{entry['component']} has no transform"
        assert transform["type"] in (
            "double_sigmoid", "reverse_sigmoid", "sigmoid", "step",
            "left_step", "right_step", "value_mapping", "exponential_decay",
        )


def test_profiles_differ_in_the_property_space_they_target():
    """Profiles must encode genuinely different chemistry.

    The anti-TB profile deliberately permits higher lipophilicity than the
    general drug-like one, because the mycobacterial envelope is lipid-rich
    and known actives sit there - bedaquiline is 555 Da at clogP above 7.
    """
    def window(profile, component):
        for entry in SCORING_PROFILES[profile]["components"]:
            if entry["component"] == component:
                transform = entry["endpoint"].get("transform", {})
                return transform.get("low"), transform.get("high")
        return None

    drug_like = window("drug-like", "SlogP")
    anti_tb = window("anti-tb", "SlogP")
    assert anti_tb[1] > drug_like[1], (
        "anti-TB should tolerate more lipophilicity than drug-like"
    )

    assert window("fragment-like", "MolecularWeight")[1] < \
        window("drug-like", "MolecularWeight")[1]
    assert window("cns", "TPSA")[1] < window("drug-like", "TPSA")[1]


def test_geometric_mean_is_the_default_aggregation():
    """One near-zero component should drag the total down.

    With an arithmetic mean a compound can ignore an objective entirely and
    compensate elsewhere, which is rarely what multi-parameter optimisation
    is for.
    """
    for profile in SCORING_PROFILES.values():
        assert profile["aggregation"] == "geometric_mean"


def test_alert_component_carries_real_smarts():
    section = build_scoring_section("drug-like")
    alerts = next(e for e in section["components"]
                  if e["component"] == "custom_alerts")
    smarts = alerts["endpoint"]["params"]["smarts"]
    assert len(smarts) >= 10
    assert all(isinstance(pattern, str) and pattern for pattern in smarts)


# ---------------------------------------------------------------------------
# Custom components
# ---------------------------------------------------------------------------


def test_component_spec_parsing():
    parsed = parse_component_spec("SlogP:low=1:high=4:weight=2")
    assert parsed["component"] == "SlogP"
    assert parsed["endpoint"]["weight"] == 2.0
    transform = parsed["endpoint"]["transform"]
    assert transform["type"] == "double_sigmoid"
    assert transform["low"] == 1.0 and transform["high"] == 4.0


def test_bare_component_spec_needs_no_transform():
    parsed = parse_component_spec("QED")
    assert parsed["component"] == "QED"
    assert parsed["endpoint"]["name"] == "QED"
    assert "transform" not in parsed["endpoint"]


def test_component_spec_accepts_an_explicit_name_and_transform():
    parsed = parse_component_spec(
        "SAScore:name=Synthesis:transform=reverse_sigmoid:low=1:high=5:k=0.4")
    assert parsed["endpoint"]["name"] == "Synthesis"
    assert parsed["endpoint"]["transform"]["type"] == "reverse_sigmoid"
    assert parsed["endpoint"]["transform"]["k"] == 0.4


@pytest.mark.parametrize("spec", ["", "SlogP:low", "SlogP:low=abc",
                                  "SlogP:weight=heavy"])
def test_malformed_component_spec_is_refused(spec):
    with pytest.raises(UsageError):
        parse_component_spec(spec)


def test_custom_components_extend_a_profile():
    section = build_scoring_section("drug-like",
                                     components=["NumAromaticRings:low=1:high=3"])
    names = [e["component"] for e in section["components"]]
    assert "NumAromaticRings" in names
    assert len(names) == len(SCORING_PROFILES["drug-like"]["components"]) + 1


def test_components_alone_need_no_profile():
    section = build_scoring_section(components=["QED", "SlogP:low=1:high=4"])
    assert section["n_components"] == 2


def test_no_objective_at_all_is_refused():
    with pytest.raises(UsageError):
        build_scoring_section()


def test_unknown_aggregation_is_refused():
    with pytest.raises(UsageError):
        build_scoring_section("drug-like", aggregation="median")


def test_unknown_profile_lists_the_options():
    with pytest.raises(UsageError) as excinfo:
        build_scoring_section("made-up")
    for profile in SCORING_PROFILES:
        assert profile in str(excinfo.value)


# ---------------------------------------------------------------------------
# Prior and installation discovery
# ---------------------------------------------------------------------------


def test_reinvent_directory_is_discovered_not_hardcoded(tmp_path,
                                                         monkeypatch):
    """No absolute path may be baked in.

    The documentation this replaces hardcoded
    ``E:\\ANTIGRAVITY_WORKSHOP\\REINVENT4`` in every example, which fails on
    anyone else's machine.
    """
    fake = tmp_path / "REINVENT4"
    (fake / "reinvent").mkdir(parents=True)

    # Explicit argument wins.
    assert find_reinvent_dir(str(fake)) == fake.resolve()

    # Then the environment variable.
    monkeypatch.setenv("REINVENT_DIR", str(fake))
    assert find_reinvent_dir() == fake.resolve()


def test_absent_reinvent_directory_returns_none_not_a_guess(tmp_path,
                                                             monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert find_reinvent_dir("/definitely/not/here") is None


def test_prior_resolution_for_each_generator():
    for name, spec in GENERATORS.items():
        filename, _ = resolve_prior(name)
        assert filename == spec["prior"]


def test_mol2mol_variants_resolve_by_short_name():
    for variant in ("scaffold_generic", "similarity", "mmp"):
        filename, _ = resolve_prior("mol2mol", variant)
        assert variant in filename
        assert filename.endswith(".prior")
    # A leading dot is tolerated, matching the old documented style.
    assert resolve_prior("mol2mol", ".scaffold_generic")[0] == \
        resolve_prior("mol2mol", "scaffold_generic")[0]


def test_unknown_mol2mol_variant_lists_the_real_ones():
    with pytest.raises(UsageError) as excinfo:
        resolve_prior("mol2mol", "nonsense")
    assert "scaffold_generic" in str(excinfo.value)


def test_unknown_generator_is_refused():
    with pytest.raises(UsageError):
        resolve_prior("not_a_generator")


def test_missing_prior_file_names_the_download(tmp_path):
    from agskills.errors import ResourceNotFoundError
    with pytest.raises(ResourceNotFoundError) as excinfo:
        resolve_prior("reinvent", priors_dir=tmp_path, require_exists=True)
    message = str(excinfo.value)
    assert "zenodo" in message.lower()
    assert "REINVENT_PRIOR_BASE" in message


def test_check_setup_reports_problems_without_raising():
    """The point of the command is to describe a broken environment."""
    report = check_setup("/nonexistent/path")
    assert report["ready"] is False
    assert report["problems"]
    assert report["advice"]
    assert report["summary"]
    # It still reports what it can.
    assert "torch" in report
    assert report["recommended_device"] in ("cpu", "cuda:0")


def test_check_setup_finds_the_local_checkout_when_present(monkeypatch,
                                                            tmp_path):
    fake = tmp_path / "REINVENT4"
    (fake / "reinvent").mkdir(parents=True)
    (fake / "priors").mkdir()
    monkeypatch.setenv("REINVENT_DIR", str(fake))
    report = check_setup()
    assert report["reinvent_dir"] == str(fake.resolve())
    assert report["priors_dir"] == str((fake / "priors").resolve())
    # No priors are present, so that must be reported as a problem.
    assert report["priors_missing"]


# ---------------------------------------------------------------------------
# Seed preparation
# ---------------------------------------------------------------------------


def test_seed_preparation_deduplicates_and_accounts_for_everything(tmp_path):
    from agskills.chem.smiles import CompoundRecord
    records = [
        CompoundRecord("a", EGFR_ACTIVES[0]),
        CompoundRecord("a_again", EGFR_ACTIVES[0]),
        CompoundRecord("b", EGFR_ACTIVES[1]),
        CompoundRecord("broken", "NOT_A_MOLECULE"),
    ]
    output = tmp_path / "seeds.smi"
    report = prepare_seeds(records, str(output))

    assert report["input_count"] == 4
    assert report["written"] == 2
    assert report["rejected"] == 2
    assert report["written"] + report["rejected"] == report["input_count"]

    reasons = {r["name"]: r["reason"] for r in report["rejected_detail"]}
    assert "duplicate" in reasons["a_again"]
    assert "invalid" in reasons["broken"]

    lines = output.read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 2


def test_seed_file_has_unix_line_endings(tmp_path):
    from agskills.chem.smiles import CompoundRecord
    output = tmp_path / "seeds.smi"
    prepare_seeds([CompoundRecord("a", "CCO")], str(output))
    assert b"\r" not in output.read_bytes()


def test_seed_preparation_warns_on_a_format_mismatch(tmp_path):
    """LibInvent needs attachment points; plain molecules will not work."""
    from agskills.chem.smiles import CompoundRecord
    report = prepare_seeds([CompoundRecord("a", "CCO")],
                           str(tmp_path / "s.smi"), generator="libinvent")
    assert any("attachment" in w for w in report["warnings"])

    report = prepare_seeds([CompoundRecord("a", "CCO")],
                           str(tmp_path / "s2.smi"), generator="linkinvent")
    assert any("warhead" in w for w in report["warnings"])


def test_seed_preparation_with_nothing_valid_raises(tmp_path):
    from agskills.chem.smiles import CompoundRecord
    with pytest.raises(InvalidInputError):
        prepare_seeds([CompoundRecord("x", "@@@")], str(tmp_path / "s.smi"))


def test_seed_names_can_be_included(tmp_path):
    from agskills.chem.smiles import CompoundRecord
    output = tmp_path / "named.smi"
    prepare_seeds([CompoundRecord("ethanol", "CCO")], str(output),
                  include_names=True)
    assert "ethanol" in output.read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------


def test_preflight_rejects_unparseable_toml(tmp_path):
    bad = tmp_path / "broken.toml"
    bad.write_text("this is not = = toml\n", encoding="utf-8")
    with pytest.raises(InvalidInputError) as excinfo:
        preflight(bad)
    assert "generate-config" in str(excinfo.value)


def test_preflight_requires_a_run_type(tmp_path):
    config = tmp_path / "no_type.toml"
    config.write_text('device = "cpu"\n', encoding="utf-8")
    with pytest.raises(InvalidInputError) as excinfo:
        preflight(config)
    assert "run_type" in str(excinfo.value)


def test_preflight_names_every_missing_file(tmp_path):
    """A missing prior must be reported before the run starts, not during."""
    from agskills.errors import ResourceNotFoundError
    config = tmp_path / "run.toml"
    config.write_text(
        'run_type = "sampling"\ndevice = "cpu"\n\n[parameters]\n'
        'model_file = "does_not_exist.prior"\n',
        encoding="utf-8")
    with pytest.raises(ResourceNotFoundError) as excinfo:
        preflight(config)
    message = str(excinfo.value)
    assert "does_not_exist.prior" in message
    assert "check-setup" in message


def test_preflight_accepts_a_config_whose_files_exist(tmp_path):
    model = tmp_path / "model.prior"
    model.write_bytes(b"not a real model, but it exists")
    config = tmp_path / "run.toml"
    config.write_text(
        f'run_type = "sampling"\ndevice = "cpu"\n\n[parameters]\n'
        f'model_file = "{model.name}"\n',
        encoding="utf-8")
    report = preflight(config, device="cpu")
    assert report["run_type"] == "sampling"
    assert len(report["files_present"]) == 1
    assert report["files_missing"] == []


# ---------------------------------------------------------------------------
# Result analysis
# ---------------------------------------------------------------------------


def test_analysis_reads_the_real_reinvent_column_format(reinvent_csv):
    result = analyze_results(reinvent_csv, top_n=10)
    assert result["smiles_column"] == "SMILES"
    assert result["score_column"] == "Score"
    # The real format has no step column.
    assert result["step_column"] is None
    assert "Agent" in result["columns"]
    assert "SMILES_state" in result["columns"]


def test_analysis_collapses_duplicate_samples(reinvent_csv):
    """An RL run re-samples the same molecule many times.

    A "top 20" taken straight from the CSV is often the same few structures
    repeated, which overstates how much was discovered.
    """
    result = analyze_results(reinvent_csv, top_n=20)
    assert result["rows_read"] == 96
    assert result["unique_molecules"] == 13
    assert result["duplicate_samples"] == 96 - 13

    smiles = [m["smiles"] for m in result["top_molecules"]]
    assert len(smiles) == len(set(smiles)), "the ranking contains duplicates"
    assert all(m["times_sampled"] >= 1 for m in result["top_molecules"])
    assert sum(m["times_sampled"] for m in result["top_molecules"]) > \
        len(result["top_molecules"])


def test_analysis_ranks_by_score_descending(reinvent_csv):
    result = analyze_results(reinvent_csv, top_n=13)
    scores = [m["score"] for m in result["top_molecules"]]
    assert scores == sorted(scores, reverse=True)


def test_analysis_reports_scaffold_diversity(reinvent_csv):
    """Scaffold collapse is the main failure mode of an RL run."""
    result = analyze_results(reinvent_csv, top_n=5)
    assert result["distinct_scaffolds"] >= 1
    assert 0.0 < result["scaffold_diversity_ratio"] <= 1.0
    assert result["most_common_scaffolds"]
    assert result["interpretation"]["diversity"]
    assert "scaffold" in result["interpretation"]["diversity"].lower()


def test_analysis_attaches_descriptors_to_the_top_molecules(reinvent_csv):
    result = analyze_results(reinvent_csv, top_n=5)
    for molecule in result["top_molecules"]:
        properties = molecule["properties"]
        assert properties["mw"] > 0
        assert "qed" in properties
        assert molecule["murcko_scaffold"] is not None


def test_analysis_can_skip_descriptors(reinvent_csv):
    result = analyze_results(reinvent_csv, top_n=3, compute_properties=False)
    assert all("properties" not in m for m in result["top_molecules"])


def test_analysis_reports_score_statistics(reinvent_csv):
    statistics = analyze_results(reinvent_csv, top_n=5)["score_statistics"]
    assert statistics["n"] == 96
    assert 0.0 <= statistics["min"] <= statistics["mean"] <= statistics["max"]
    assert statistics["stdev"] > 0


def test_min_score_filters_before_ranking(reinvent_csv):
    everything = analyze_results(reinvent_csv, top_n=20)
    filtered = analyze_results(reinvent_csv, top_n=20, min_score=0.8)
    assert filtered["rows_below_min_score"] > 0
    assert filtered["unique_molecules"] <= everything["unique_molecules"]
    assert all(m["score"] >= 0.8 for m in filtered["top_molecules"])


def test_top_n_limits_the_ranking(reinvent_csv):
    assert len(analyze_results(reinvent_csv, top_n=3)["top_molecules"]) == 3


def test_sorting_by_a_component_column(reinvent_csv):
    result = analyze_results(reinvent_csv, top_n=5, sort_by="QED drug-likeness")
    assert result["score_column"] == "QED drug-likeness"


def test_unknown_sort_column_lists_the_real_ones(reinvent_csv):
    with pytest.raises(InvalidInputError) as excinfo:
        analyze_results(reinvent_csv, sort_by="NotAColumn")
    assert "SMILES" in str(excinfo.value)


def test_csv_without_a_structure_column_is_refused(tmp_path):
    csv_file = tmp_path / "wrong.csv"
    csv_file.write_text("Score,Agent\n0.5,10\n", encoding="utf-8")
    with pytest.raises(InvalidInputError) as excinfo:
        analyze_results(csv_file)
    assert "no SMILES column" in str(excinfo.value)


def test_empty_csv_is_refused(tmp_path):
    csv_file = tmp_path / "header_only.csv"
    csv_file.write_text("SMILES,Score\n", encoding="utf-8")
    with pytest.raises(InvalidInputError):
        analyze_results(csv_file)


def test_invalid_top_n_is_refused(reinvent_csv):
    with pytest.raises(InvalidInputError):
        analyze_results(reinvent_csv, top_n=0)


def test_analysis_counts_rows_without_a_structure(tmp_path):
    csv_file = tmp_path / "gappy.csv"
    csv_file.write_text(
        "SMILES,Score\nCCO,0.5\n,0.9\nCCC,0.7\n", encoding="utf-8")
    result = analyze_results(csv_file, top_n=5)
    assert result["rows_read"] == 3
    assert result["rows_without_structure"] == 1
    assert result["unique_molecules"] == 2


def test_analysis_output_serialises(reinvent_csv):
    import json
    json.dumps(analyze_results(reinvent_csv, top_n=5))
