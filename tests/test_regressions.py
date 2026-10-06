"""One test per defect found in the implementation this replaces.

Each test names the original behaviour, so if a future change reintroduces
it the failure says which known problem came back. The defects are grouped
by how they failed:

* **Crashes** - the command could not run at all.
* **Silent wrong answers** - the worst class: the command succeeded and
  reported something untrue.
* **Scientific errors** - the computation was not what the documentation
  or the literature says.
* **Portability** - worked on the author's machine only.
* **Documentation drift** - the docs described something the code did not do.
"""

from __future__ import annotations

import json

import pytest

from agskills.chem.admet import TIERS, screen_compounds, tier_summary
from agskills.chem.descriptors import compute_properties
from agskills.chem.rules import lipinski
from agskills.chem.smiles import CompoundRecord, load_compounds, parse_smiles
from agskills.errors import UsageError


def properties_for(smiles: str):
    return compute_properties(parse_smiles(smiles))


# ===========================================================================
# Crashes
# ===========================================================================


def test_docking_box_arguments_do_not_crash():
    """``--size`` was defined while ``args.box_size`` was read.

    Original: ``p_dock.add_argument("--size", ...)`` but
    ``box_size = tuple(args.box_size)``, so every ``--center`` invocation
    raised ``AttributeError: 'Namespace' object has no attribute
    'box_size'`` before any work was done.
    """
    from agskills.cli.compound_screening import build_parser
    parser = build_parser()
    args = parser.parse_args([
        "dock", "--receptor", "r.pdbqt", "--smiles", "CCO",
        "--center", "1", "2", "3", "--box-size", "20", "20", "20",
        "--output", "out.json",
    ])
    assert hasattr(args, "box_size"), "the parsed name must match the reader"
    assert args.box_size == [20.0, 20.0, 20.0]
    assert args.center == [1.0, 2.0, 3.0]


def test_a_missing_affinity_is_not_formatted_as_a_float():
    """``print(f"{score:.1f}")`` ran unconditionally.

    Original: when Vina produced no parseable score, ``score`` was ``None``
    and the status line raised ``TypeError: unsupported format string passed
    to NoneType.__format__`` - after the docking work was already finished.
    """
    from agskills.docking.vina import DockingResult
    result = DockingResult(name="x", smiles="CCO", status="no_pose")
    data = result.as_dict()
    json.dumps(data)
    assert data["best_affinity_kcal_mol"] is None


def test_the_vina_parser_accepts_an_integer_rmsd():
    """Real Vina prints mode 1's RMSD columns as a bare ``0``.

    A pattern insisting on a decimal point drops the best pose - the one
    that matters most. Found by testing against authentic v1.2.5 output.
    """
    from agskills.docking.vina import parse_vina_log
    poses = parse_vina_log("   1       -9.428          0          0\n")
    assert len(poses) == 1
    assert poses[0].affinity_kcal_mol == pytest.approx(-9.428)


# ===========================================================================
# Silent wrong answers
# ===========================================================================


def test_compound_names_cannot_silently_truncate_the_input():
    """``list(zip(args.smiles, names))`` discarded the excess.

    Original: three SMILES with one name produced a one-compound run that
    reported success. Two thirds of the input vanished without a word.
    """
    with pytest.raises(UsageError):
        load_compounds(smiles=["CCO", "CCC", "CCCC"], names=["one"])


def test_screening_counts_balance_with_invalid_input():
    """Passes and failures were counted with incompatible defaults.

    Original: ``sum(1 for c if c.get("overall_pass", False))`` for passes
    and ``sum(1 for c if not c.get("overall_pass", True))`` for failures.
    A compound with an unparseable SMILES has no ``overall_pass`` key, so
    it was counted in neither. Verified on the real code, three inputs gave
    ``total=3, passed=2, failed=0``.
    """
    records = [
        CompoundRecord("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        CompoundRecord("ethanol", "CCO"),
        CompoundRecord("bogus", "NOT_A_SMILES"),
    ]
    report = screen_compounds(records, strictness="relaxed")
    assert report["total_compounds"] == 3
    assert report["passed"] + report["failed"] + report["errored"] == 3
    assert report["errored"] == 1


def test_a_zinc_download_url_cannot_serve_a_web_page():
    """The generated URL returned HTML with HTTP 200.

    Original: ``zinc20.docking.org/tranches/{code}/download.smi`` serves a
    12 KB HTML page, not compound data. The generated script tested only
    that a file existed and counted its lines, so every tranche appeared to
    succeed and the pipeline continued with web pages saved as ``.smi``.
    """
    from agskills.libraries.zinc import generate_zinc_script, select_tranches, tranche_url

    url = tranche_url("CDAA")
    assert "files.docking.org" in url
    assert "zinc20.docking.org" not in url

    # And the script defends against the failure mode regardless.
    script = generate_zinc_script(select_tranches("fragment-like"))
    assert "doctype html" in script.lower()
    assert "not compound data" in script


def test_the_zinc_logp_letter_table_is_not_shifted():
    """The original logP mapping was off by one bin.

    Original table: ``D`` at logP 0-1, ``E`` at 1-2. Measured from live
    tranches with RDKit: ``D`` holds 1-2, ``E`` holds 2-2.5. Even with a
    working URL the old table downloaded the wrong chemistry.
    """
    from agskills.libraries.zinc import LOGP_BINS
    bins = {letter: (low, high) for letter, low, high in LOGP_BINS}
    assert bins["D"] == (1, 2), "D must be logP 1-2, not 0-1"
    assert bins["E"] == (2, 2.5), "E must be logP 2-2.5, not 1-2"


def test_the_requested_md_temperature_reaches_the_mdp_file():
    """``--temperature`` was accepted and then discarded.

    Original: the flag was interpolated into a comment in the generated
    Python workflow, but ``generate_mdp()`` took only
    ``(mdp_type, production_ns)``. Requesting 310 K produced
    ``ref_t = 300 300``, and grompp reads the MDP - so the simulation ran
    at 300 K while the user believed otherwise. Confirmed by running the
    original code.
    """
    from agskills.md.mdp import generate_mdp, parse_mdp
    for stage in ("nvt", "npt", "production"):
        values = parse_mdp(generate_mdp(stage, temperature=310.0))
        assert all(float(v) == 310.0 for v in values["ref-t"].split()), stage
    assert float(parse_mdp(generate_mdp("nvt", temperature=310.0))[
        "gen-temp"]) == 310.0


def test_binding_site_selection_prefers_the_inhibitor_over_the_cofactor(
        pdb_6hez):
    """Ranking heterogens by size picks the cofactor.

    In PDB 6HEZ (*M. tuberculosis* DprE1) the FAD cofactor has 53 atoms
    and the inhibitor 0SK has 28. A size-only ranking centres the docking
    box on the flavin site, so every affinity afterwards is for the wrong
    cavity. Observed: the box moved from (6.71, -9.95, 39.04) to
    (13.87, -20.62, 37.13) once cofactors were de-prioritised.
    """
    from agskills.targets.site import list_ligands, site_from_hetatm
    ligands = list_ligands(pdb_6hez)
    assert ligands[0]["resname"] == "0SK"
    fad = next(e for e in ligands if e["resname"] == "FAD")
    assert fad["n_atoms"] > ligands[0]["n_atoms"]
    assert site_from_hetatm(pdb_6hez).reference.startswith("0SK")


def test_multichain_plddt_is_not_overwritten(pdb_6hez):
    """pLDDT was keyed on residue number alone.

    Original: ``plddt_by_residue: dict[int, float]``, so residue 42 of
    chain B overwrote residue 42 of chain A and a dimer reported the
    confidence of roughly half its residues.
    """
    from agskills.targets.quality import extract_plddt
    residues = extract_plddt(pdb_6hez.read_text(encoding="utf-8"))
    keys = [(r.chain, r.resseq, r.icode) for r in residues]
    assert len(keys) == len(set(keys))
    assert len({r.chain for r in residues}) >= 2


def test_the_cocrystallised_ligand_is_actually_removed(pdb_1m17, tmp_path):
    """The cleaner removed only water, despite documenting more.

    Original: ``accept_residue`` tested ``hetflag == "W"``. BioPython flags
    water as ``"W"`` and every other heterogen as ``"H_<resname>"``, so
    ligands, cofactors and additives all survived into the docking
    receptor.
    """
    from agskills.targets.clean import clean_structure
    from agskills.targets.site import parse_atoms

    output = tmp_path / "clean.pdb"
    clean_structure(pdb_1m17, output, chain="A")
    remaining = {a.resname for a in
                 parse_atoms(output.read_text(encoding="utf-8"))
                 if a.is_hetatm}
    assert "AQ4" not in remaining, "erlotinib survived cleaning"
    assert "HOH" not in remaining


def test_a_skipped_check_is_not_reported_as_a_pass():
    """Zero checks performed yielded ``overall_pass: True``.

    Original: ``overall_pass = len(flags) == 0``. When a prediction source
    returned no recognised keys, no check ran, ``flags`` was empty, and the
    compound was reported as passing - a false positive on no evidence.
    """
    report = screen_compounds([CompoundRecord("x", "@@@")],
                              strictness="strict")
    compound = report["compounds"][0]
    assert compound["overall_pass"] is False
    assert compound["status"] == "error"


# ===========================================================================
# Scientific errors
# ===========================================================================


def test_molecular_weight_is_the_average_mass_everywhere():
    """One path used MolWt, another ExactMolWt, both labelled "mw".

    Original: the docking filter used ``Descriptors.MolWt`` while the ADMET
    fallback used ``Descriptors.ExactMolWt``. Verified by running the
    original code: aspirin was reported as ``mw: 180.04``, the
    monoisotopic mass, when its molecular weight is 180.16 - and that value
    was then compared against Lipinski's 500 Da limit.
    """
    aspirin = properties_for("CC(=O)Oc1ccccc1C(=O)O")
    assert aspirin.mw == pytest.approx(180.16, abs=0.02)
    assert aspirin.exact_mw == pytest.approx(180.042, abs=0.01)

    # Both tiers report the same, average, mass.
    for strictness in ("relaxed", "strict"):
        report = screen_compounds(
            [CompoundRecord("aspirin", "CC(=O)Oc1ccccc1C(=O)O")], strictness)
        assert report["compounds"][0]["properties"]["mw"] == \
            pytest.approx(180.16, abs=0.02)


def test_lipinski_violations_are_absolute_not_tier_relative():
    """The tier's limits were substituted into the rule itself.

    Original: in relaxed mode ``mw_max`` became 700, and violations were
    counted against that, so a 650 Da compound was reported with
    ``lipinski_violations: 0``. That is false: 650 > 500 is a rule-of-five
    violation whatever the screening policy.
    """
    heavy = properties_for(
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)NCCCCNC(=O)c1ccccc1")
    assert heavy.mw > 500
    result = lipinski(heavy)
    assert result.n_violations >= 1
    assert any("> 500" in violation for violation in result.violations)

    # Reported identically at both tiers, including the limit named.
    for strictness in ("relaxed", "strict"):
        report = screen_compounds([CompoundRecord("heavy", heavy.smiles)],
                                  strictness)
        check = next(c for c in report["compounds"][0]["checks"]
                     if c["check"] == "lipinski")
        assert check["value"] >= 1, strictness
        assert "500" in check["detail"], strictness


def test_lipinski_permits_one_violation():
    """Two definitions of the rule coexisted in one codebase.

    Original: the docking filter used ``violations <= max_violations`` with
    1 for strict, while the ADMET fallback used ``violations == 0``. The
    1997 paper predicts poor absorption at two or more violations, so up
    to one passes.
    """
    assert lipinski(properties_for("CCO")).max_allowed_violations == 1


def test_charmm_and_amber_non_bonded_settings_are_not_interchangeable():
    """AMBER settings were emitted under a CHARMM36 default.

    Original: the default force field was ``charmm36`` while every MDP
    template carried ``rvdw = 1.0`` with ``DispCorr = EnerPres``, which is
    the AMBER convention. CHARMM36 requires a force-switched
    Lennard-Jones potential between 1.0 and 1.2 nm and no dispersion
    correction; adding one double-counts a tail the parameters already
    absorb. Cross-checking the original md.mdp against charmm36 produces
    four distinct errors.
    """
    from agskills.md.mdp import generate_mdp, parse_mdp, validate_mdp_text

    charmm = parse_mdp(generate_mdp("production", force_field="charmm36"))
    assert charmm["vdw-modifier"] == "force-switch"
    assert float(charmm["rvdw"]) == 1.2
    assert "dispcorr" not in charmm

    amber_style = (
        "integrator = md\nnsteps = 1000\ndt = 0.002\ncoulombtype = PME\n"
        "rcoulomb = 1.0\nrvdw = 1.0\nDispCorr = EnerPres\n"
        "tcoupl = v-rescale\ntc-grps = System\ntau_t = 0.1\nref_t = 300\n"
    )
    report = validate_mdp_text(amber_style, force_field="charmm36")
    assert not report["valid"]
    assert len(report["errors"]) >= 3


def test_position_restraints_with_pressure_coupling_scale_the_reference():
    """``refcoord_scaling`` was missing from the NPT stage.

    Without it grompp warns and the restraint reference coordinates are not
    scaled with the box. Validating the original npt.mdp reports exactly
    this.
    """
    from agskills.md.mdp import generate_mdp, parse_mdp
    npt = parse_mdp(generate_mdp("npt"))
    assert "-DPOSRES" in npt["define"]
    assert npt["pcoupl"] != "no"
    assert npt["refcoord-scaling"] == "com"


def test_pains_screening_reports_every_alert_not_just_the_first():
    """Two different PAINS configurations were used on two code paths.

    Original: ``PAINS_A|PAINS_B|PAINS_C`` with ``GetMatches`` in the
    docking filter, but the combined ``PAINS`` catalogue with
    ``GetFirstMatch`` in the ADMET fallback. The second form can only ever
    report one alert, so the same molecule had one or four alerts depending
    on which subcommand ran.
    """
    from agskills.chem.alerts import screen_alerts
    # Catechol, a classic PAINS motif.
    mol = parse_smiles("Oc1ccccc1O")
    report = screen_alerts(mol, ("pains", "brenk"))
    assert report.n_alerts >= 1
    assert len(report.hits) == report.n_alerts
    assert set(report.catalogs_screened) == {"pains", "brenk"}


def test_structure_confidence_is_not_applied_to_crystallographic_b_factors(
        pdb_1m17):
    """A pLDDT verdict on a crystal structure is meaningless.

    1M17 has a mean B-factor near 30, which the original banding would
    call "very low confidence - not recommended for docking". It is a
    2.6 A crystal structure and perfectly usable; for B-factors lower is
    better and the 0-100 scale does not apply.
    """
    from agskills.targets.quality import assess_structure
    report = assess_structure(pdb_1m17)
    assert report.verdict == "NOT_APPLICABLE"
    assert "LOWER is better" in report.b_factor_caveat


def test_a_toxicity_risk_flag_does_not_reject_lead_like_chemistry():
    """Gating triage on Pfizer 3/75 empties the seed set.

    The criterion flags logP above 3 with TPSA below 75, which describes
    much of lead-like space: ibuprofen and the whole 4-anilinoquinazoline
    kinase series trip it. It was published as a signal to look harder at
    toxicity, not as an exclusion filter. Found while testing: the relaxed
    tier was rejecting every real EGFR active.
    """
    from panel import EGFR_ACTIVES
    records = [CompoundRecord(f"egfr_{i}", s)
               for i, s in enumerate(EGFR_ACTIVES)]
    report = screen_compounds(records, strictness="relaxed")
    assert report["passed"] == len(EGFR_ACTIVES)
    # The flag is still reported.
    assert any(c["advisories"] for c in report["compounds"])


# ===========================================================================
# Portability
# ===========================================================================


def test_no_absolute_developer_path_appears_in_the_package():
    """The documentation hardcoded one machine's directory layout.

    Original: every REINVENT example command contained
    ``E:\\ANTIGRAVITY_WORKSHOP\\REINVENT4``, which fails on anyone else's
    computer.
    """
    import pathlib
    root = pathlib.Path(__file__).resolve().parents[1]
    offenders = []
    for path in list((root / "src").rglob("*.py")) + \
            list((root / ".agents").rglob("*.md")) + \
            list((root / ".agents").rglob("*.py")):
        text = path.read_text(encoding="utf-8", errors="replace")
        for needle in ("E:\\ANTIGRAVITY", "E:/ANTIGRAVITY",
                       "C:\\Users\\Walter", "/home/walter"):
            if needle.lower() in text.lower():
                offenders.append(f"{path.relative_to(root)} contains {needle}")
    assert not offenders, offenders


def test_the_cache_directory_is_not_a_hardcoded_posix_path():
    """``robots_cache_dir="/tmp"`` does not exist on Windows.

    Original: every cache write failed silently, so ``robots.txt`` was
    re-fetched on every client construction - an extra network round trip
    per command.
    """
    import os

    from agskills.http import cache_dir
    directory = cache_dir()
    assert os.path.isdir(directory)
    if os.name == "nt":
        assert not directory.startswith("/tmp")


def test_the_default_user_agent_is_not_empty():
    """An empty User-Agent earns a 403 from several endpoints.

    Original: ``DEFAULT_USER_AGENT = os.environ.get(..., "")``.
    """
    from agskills.http import DEFAULT_USER_AGENT
    assert DEFAULT_USER_AGENT.strip()


def test_the_reinvent_checkout_is_discovered_not_assumed(tmp_path,
                                                          monkeypatch):
    from agskills.reinvent.priors import find_reinvent_dir
    fake = tmp_path / "REINVENT4"
    (fake / "reinvent").mkdir(parents=True)
    monkeypatch.setenv("REINVENT_DIR", str(fake))
    assert find_reinvent_dir() == fake.resolve()


def test_generated_scripts_use_unix_line_endings(tmp_path):
    """A CRLF shell script fails on a Linux cluster before it runs."""
    from agskills.hpc import SlurmResources, generate_slurm_script
    from agskills.io_utils import write_text
    script = generate_slurm_script(
        "reinvent", SlurmResources(mem="64G", time="24:00:00"),
        {"config": "run.toml"})
    path = write_text(script, tmp_path / "submit.sh")
    assert b"\r" not in path.read_bytes()


def test_reports_are_written_as_utf8_not_the_locale_codepage(tmp_path):
    """``Path.write_text`` without an encoding is cp1252 on Windows.

    Original: the structure-prediction quality report was written with a
    bare ``write_text``, so a report containing an en dash raised
    ``UnicodeEncodeError``.
    """
    from agskills.io_utils import write_text
    content = "mean pLDDT 85.2 - PROCEED; box 18 A wide\n"
    path = write_text(content, tmp_path / "report.md")
    assert path.read_text(encoding="utf-8") == content


def test_the_package_source_is_ascii_only():
    """Non-ASCII in a user-facing string is mangled by a cp1252 console.

    Observed before the fix: "high - the box is centred" printed as
    "high ? the box is centred".
    """
    import pathlib
    root = pathlib.Path(__file__).resolve().parents[1] / "src"
    offenders = []
    for path in sorted(root.rglob("*.py")):
        for number, line in enumerate(
                path.read_text(encoding="utf-8").splitlines(), 1):
            if any(ord(char) > 127 for char in line):
                offenders.append(f"{path.name}:{number}")
    assert not offenders, offenders


# ===========================================================================
# Documentation drift
# ===========================================================================


def test_the_documented_check_count_is_the_performed_check_count():
    """The docs claimed 14 and 34 checks; the code ran three.

    Original: ``admet-filter`` performed exactly three checks (Lipinski,
    PAINS, synthetic accessibility) while the skill documentation and the
    wizard's pipeline diagram both advertised "14 checks" and
    "34 checks". The counts are now derived from the registry, so they
    cannot drift.
    """
    summary = tier_summary()
    for name, tier in TIERS.items():
        assert summary[name] == len(tier.checks)
        report = screen_compounds([CompoundRecord("x", "CCO")], name)
        performed = report["compounds"][0]
        assert (performed["checks_passed"] + performed["checks_failed"]
                + performed["checks_skipped"]) == summary[name]
        assert report["checks_per_compound"] == summary[name]
    # And the real numbers are substantially more than three.
    assert summary["relaxed"] >= 10
    assert summary["strict"] >= 20


def test_the_skill_documentation_states_the_real_check_counts():
    """The SKILL.md files must agree with the registry."""
    import pathlib
    import re

    root = pathlib.Path(__file__).resolve().parents[1] / ".agents" / "skills"
    summary = tier_summary()
    checked = 0
    # The word boundary matters: without it this also matches the "4 check"
    # inside "REINVENT 4 checkout".
    pattern = re.compile(r"(\d+)\s+checks?\b")
    for path in root.rglob("SKILL.md"):
        text = path.read_text(encoding="utf-8")
        for match in pattern.finditer(text):
            count = int(match.group(1))
            assert count in summary.values(), (
                f"{path.relative_to(root)} claims {count} checks; the "
                f"registry performs {summary}"
            )
            checked += 1
    assert checked > 0, "no check counts found in the skill documentation"


def test_the_binding_site_file_the_docs_reference_can_be_produced():
    """``dock --site binding_site.json`` had no producer.

    Original: the docking subcommand documented a ``--site`` input and the
    SLURM generator defaulted to ``binding_site.json``, but no subcommand
    in any skill could create that file - and the alternative,
    ``--center``, crashed.
    """
    from agskills.cli.target_preparation import build_parser
    parser = build_parser()
    subcommands = [a for a in parser._actions
                   if hasattr(a, "choices") and a.choices
                   and hasattr(next(iter(a.choices.values())), "_actions")]
    assert "define-site" in subcommands[0].choices


def test_library_queries_go_through_the_rate_limited_client():
    """The docs promised rate limiting that the code bypassed.

    Original: the skill's Core Rules said the wrapper "enforces rate
    limits and handles retries", but ``compound_libraries.py`` and
    ``predict_structure.py`` used raw ``urllib.request`` with a hand-rolled
    ``time.sleep(0.5)`` - no shared budget, no retry, no backoff.
    """
    import inspect

    from agskills.libraries import chembl, coconut
    from agskills.targets import esmfold, fetch

    for module in (chembl, coconut, esmfold, fetch):
        source = inspect.getsource(module)
        assert "HttpClient" in source, module.__name__
        assert "urllib.request.urlopen" not in source, (
            f"{module.__name__} bypasses the rate-limited client"
        )
