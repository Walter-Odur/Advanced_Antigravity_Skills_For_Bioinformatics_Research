"""Command-line interfaces for all five skills.

Every subcommand is exercised through its real ``main(argv)`` entry point,
so the argument parser, the handler, the error translation and the output
contract are all covered together.

The shared contract under test:

* ``--output`` is required and receives machine-readable results.
* A deliberate failure produces a one-line message plus a hint on stderr
  and a documented exit code - never a traceback, because the exit codes
  are part of the interface.
* Nothing calls ``sys.exit`` from inside the science layer, which is what
  made the previous implementation untestable without catching
  ``SystemExit``.
"""

from __future__ import annotations

import json

import pytest

from agskills.cli import (
    compound_screening,
    compound_synthesis,
    md_simulation,
    target_preparation,
    wizard,
)
from agskills.errors import (
    InvalidInputError,
    ResourceNotFoundError,
    UsageError,
)

ALL_CLIS = {
    "target-preparation": target_preparation,
    "compound-screening": compound_screening,
    "compound-synthesis": compound_synthesis,
    "md-simulation": md_simulation,
    "wizard": wizard,
}


def run(module, *argv: str) -> int:
    return module.main(list(argv))


# ---------------------------------------------------------------------------
# Shared contract
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(ALL_CLIS))
def test_help_exits_zero(name, capsys):
    with pytest.raises(SystemExit) as excinfo:
        ALL_CLIS[name].main(["--help"])
    assert excinfo.value.code == 0
    assert capsys.readouterr().out


@pytest.mark.parametrize("name", sorted(ALL_CLIS))
def test_version_is_reported(name, capsys):
    with pytest.raises(SystemExit) as excinfo:
        ALL_CLIS[name].main(["--version"])
    assert excinfo.value.code == 0
    assert "2.0.0" in capsys.readouterr().out


@pytest.mark.parametrize("name", sorted(ALL_CLIS))
def test_no_subcommand_is_a_usage_error(name):
    with pytest.raises(SystemExit) as excinfo:
        ALL_CLIS[name].main([])
    assert excinfo.value.code == UsageError.exit_code


@pytest.mark.parametrize("name", sorted(ALL_CLIS))
def test_unknown_subcommand_is_a_usage_error(name):
    with pytest.raises(SystemExit) as excinfo:
        ALL_CLIS[name].main(["not-a-subcommand"])
    assert excinfo.value.code == UsageError.exit_code


@pytest.mark.parametrize("name", sorted(ALL_CLIS))
def test_every_subcommand_requires_output(name):
    """``--output`` is the documented contract for all of them."""
    module = ALL_CLIS[name]
    parser = module.build_parser()
    subparsers = [action for action in parser._actions
                  if hasattr(action, "choices") and action.choices
                  and hasattr(next(iter(action.choices.values())), "_actions")]
    assert subparsers, f"{name} exposes no subcommands"
    for subcommand, subparser in subparsers[0].choices.items():
        output = [a for a in subparser._actions if a.dest == "output"]
        assert output, f"{name} {subcommand} has no --output"
        assert output[0].required, f"{name} {subcommand}: --output optional"


@pytest.mark.parametrize("name", sorted(ALL_CLIS))
def test_every_subcommand_supports_quiet_and_verbose(name):
    parser = ALL_CLIS[name].build_parser()
    subparsers = [a for a in parser._actions
                  if hasattr(a, "choices") and a.choices
                  and hasattr(next(iter(a.choices.values())), "_actions")]
    for subcommand, subparser in subparsers[0].choices.items():
        dests = {a.dest for a in subparser._actions}
        assert "quiet" in dests, f"{name} {subcommand}"
        assert "verbose" in dests, f"{name} {subcommand}"


# ---------------------------------------------------------------------------
# compound-screening
# ---------------------------------------------------------------------------


def test_admet_filter_writes_a_report(tmp_path, capsys):
    output = tmp_path / "admet.json"
    code = run(compound_screening, "admet-filter",
               "--smiles", "CC(=O)Oc1ccccc1C(=O)O", "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
               "--names", "aspirin", "caffeine",
               "--strictness", "relaxed", "--output", str(output))
    assert code == 0
    data = json.loads(output.read_text(encoding="utf-8"))
    assert data["total_compounds"] == 2
    assert data["passed"] + data["failed"] + data["errored"] == 2
    assert "Tier 1" in capsys.readouterr().out


def test_admet_filter_rejects_mismatched_names(tmp_path, capsys):
    """Silent truncation is the defect; a clear error is the fix."""
    code = run(compound_screening, "admet-filter",
               "--smiles", "CCO", "CCC", "CCCC",
               "--names", "only_one",
               "--output", str(tmp_path / "x.json"))
    assert code == UsageError.exit_code
    captured = capsys.readouterr()
    assert "error:" in captured.err
    assert "one to one" in captured.err
    assert "hint:" in captured.err
    assert not (tmp_path / "x.json").exists()


def test_admet_filter_rejects_an_unknown_tier(tmp_path, capsys):
    with pytest.raises(SystemExit) as excinfo:
        compound_screening.main(["admet-filter", "--smiles", "CCO",
                                 "--strictness", "medium",
                                 "--output", str(tmp_path / "x.json")])
    assert excinfo.value.code == UsageError.exit_code


def test_admet_filter_reads_a_csv(tmp_path):
    csv_file = tmp_path / "compounds.csv"
    csv_file.write_text("SMILES,Name\nCCO,ethanol\nCCC,propane\n",
                        encoding="utf-8")
    output = tmp_path / "out.json"
    assert run(compound_screening, "admet-filter", "--csv", str(csv_file),
               "--output", str(output)) == 0
    assert json.loads(output.read_text(encoding="utf-8"))[
        "total_compounds"] == 2


def test_admet_filter_passing_only_filters_the_output(tmp_path):
    output = tmp_path / "out.json"
    run(compound_screening, "admet-filter",
        "--smiles", "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC", "C",
        "--names", "good", "methane",
        "--strictness", "relaxed", "--passing-only", "--output", str(output))
    data = json.loads(output.read_text(encoding="utf-8"))
    assert all(c["overall_pass"] for c in data["compounds"])
    assert data["total_compounds"] == 2, "counts describe the whole input"


def test_quiet_suppresses_the_summary(tmp_path, capsys):
    run(compound_screening, "admet-filter", "--smiles", "CCO",
        "--quiet", "--output", str(tmp_path / "out.json"))
    assert capsys.readouterr().out == ""


def test_describe_tiers_documents_the_check_registry(tmp_path):
    output = tmp_path / "tiers.json"
    assert run(compound_screening, "describe-tiers",
               "--output", str(output)) == 0
    data = json.loads(output.read_text(encoding="utf-8"))
    assert set(data["check_counts"]) == {"relaxed", "strict"}
    for name, tier in data["tiers"].items():
        assert tier["n_checks"] == data["check_counts"][name]
        assert len(tier["checks"]) == tier["n_checks"]


def test_query_zinc_writes_a_script_and_its_metadata(tmp_path):
    script = tmp_path / "download.sh"
    assert run(compound_screening, "query-zinc", "--subset", "fragment-like",
               "--output", str(script)) == 0
    assert script.is_file()
    assert b"\r" not in script.read_bytes()
    metadata = script.with_suffix(".json")
    assert metadata.is_file()
    data = json.loads(metadata.read_text(encoding="utf-8"))
    assert data["n_tranches"] >= 1
    assert data["tranches"][0]["url"].startswith("https://files.docking.org")


def test_dock_requires_a_box(tmp_path, capsys):
    receptor = tmp_path / "r.pdbqt"
    receptor.write_text("ATOM\n", encoding="utf-8")
    with pytest.raises(SystemExit) as excinfo:
        compound_screening.main(["dock", "--receptor", str(receptor),
                                 "--smiles", "CCO",
                                 "--output", str(tmp_path / "d.json")])
    assert excinfo.value.code == UsageError.exit_code


def test_dock_rejects_a_malformed_site_file(tmp_path, capsys):
    receptor = tmp_path / "r.pdbqt"
    receptor.write_text("ATOM\n", encoding="utf-8")
    site = tmp_path / "site.json"
    site.write_text('{"not": "a site"}', encoding="utf-8")
    code = run(compound_screening, "dock", "--receptor", str(receptor),
               "--smiles", "CCO", "--site", str(site),
               "--output", str(tmp_path / "d.json"))
    assert code == UsageError.exit_code
    assert "define-site" in capsys.readouterr().err


def test_dock_accepts_center_and_box_size_without_crashing(tmp_path, capsys):
    """The old parser defined ``--size`` but read ``args.box_size``.

    Any ``--center`` invocation therefore died with AttributeError before
    doing any work. Here the command must get far enough to report a real,
    actionable problem instead.
    """
    receptor = tmp_path / "r.pdbqt"
    receptor.write_text("ATOM\n", encoding="utf-8")
    code = run(compound_screening, "dock", "--receptor", str(receptor),
               "--smiles", "CCO",
               "--center", "13.87", "-20.62", "37.13",
               "--box-size", "18", "18", "15",
               "--output", str(tmp_path / "d.json"))
    captured = capsys.readouterr()
    assert "AttributeError" not in captured.err
    assert "box_size" not in captured.err
    # Without Vina installed this is a dependency error, which is correct.
    assert code in (0, 4)
    if code == 4:
        assert "Vina" in captured.err


def test_admet_predict_is_honest_about_the_online_service(tmp_path):
    """No public ADMETlab REST endpoint is reachable.

    Reporting that plainly, while still returning the deterministic offline
    result, is better than presenting a failure as a result.
    """
    output = tmp_path / "admet.json"
    assert run(compound_screening, "admet-predict", "--smiles", "CCO",
               "--output", str(output)) == 0
    data = json.loads(output.read_text(encoding="utf-8"))
    online = data["online_prediction"]
    assert online["attempted"] is False
    assert "404" in online["reason"]
    assert online["how_to_add_it"]
    assert data["engine"].startswith("RDKit")


# ---------------------------------------------------------------------------
# target-preparation
# ---------------------------------------------------------------------------


def test_assess_structure_on_a_real_alphafold_model(alphafold_model, tmp_path,
                                                    capsys):
    output = tmp_path / "quality.json"
    assert run(target_preparation, "assess-structure",
               "--pdb", str(alphafold_model), "--output", str(output)) == 0
    data = json.loads(output.read_text(encoding="utf-8"))
    assert data["n_residues"] == 125
    assert data["verdict"] in ("PROCEED", "REFINE", "STOP")
    assert "VERDICT" in capsys.readouterr().out


def test_assess_structure_refuses_a_crystal_b_factor_verdict(pdb_1m17,
                                                              tmp_path):
    output = tmp_path / "quality.json"
    run(target_preparation, "assess-structure", "--pdb", str(pdb_1m17),
        "--output", str(output))
    data = json.loads(output.read_text(encoding="utf-8"))
    assert data["verdict"] == "NOT_APPLICABLE"


def test_define_site_from_a_real_ligand(pdb_1m17, tmp_path, capsys):
    output = tmp_path / "site.json"
    assert run(target_preparation, "define-site", "--pdb", str(pdb_1m17),
               "--from-ligand", "--output", str(output)) == 0
    site = json.loads(output.read_text(encoding="utf-8"))
    for key in ("center_x", "center_y", "center_z", "size_x", "size_y",
                "size_z"):
        assert key in site
    assert site["reference"] == "AQ4:A:999"
    assert "dock" in site["usage"]
    assert "AQ4" in capsys.readouterr().out


def test_define_site_picks_the_inhibitor_over_the_cofactor(pdb_6hez,
                                                           tmp_path):
    output = tmp_path / "site.json"
    run(target_preparation, "define-site", "--pdb", str(pdb_6hez),
        "--from-ligand", "--output", str(output))
    site = json.loads(output.read_text(encoding="utf-8"))
    assert site["reference"].startswith("0SK")


def test_define_site_requires_exactly_one_selector(pdb_1m17, tmp_path,
                                                    capsys):
    """Two selectors is a usage error reported through the return code.

    Argparse-level mistakes (an unknown flag) raise ``SystemExit``;
    handler-level ones return the documented exit code from ``main()``.
    Both are part of the contract, and neither produces a traceback.
    """
    code = run(target_preparation, "define-site", "--pdb", str(pdb_1m17),
               "--from-ligand", "--whole-structure",
               "--output", str(tmp_path / "s.json"))
    assert code == UsageError.exit_code
    captured = capsys.readouterr()
    assert "exactly one" in captured.err
    assert "Traceback" not in captured.err
    assert not (tmp_path / "s.json").exists()


def test_define_site_from_residues(pdb_1m17, tmp_path):
    output = tmp_path / "site.json"
    assert run(target_preparation, "define-site", "--pdb", str(pdb_1m17),
               "--residues", "A:766", "A:769",
               "--output", str(output)) == 0
    assert json.loads(output.read_text(encoding="utf-8"))["method"] == \
        "pocket residue centroid"


def test_list_ligands_labels_cofactors(pdb_6hez, tmp_path):
    output = tmp_path / "ligands.json"
    assert run(target_preparation, "list-ligands", "--pdb", str(pdb_6hez),
               "--output", str(output)) == 0
    data = json.loads(output.read_text(encoding="utf-8"))
    roles = {entry["resname"]: entry["is_cofactor"]
             for entry in data["ligands"]}
    assert roles["FAD"] is True
    assert roles["0SK"] is False


def test_validate_sequence_accepts_a_real_sequence(tmp_path):
    sequence = ("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVK"
                "VKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQ")
    output = tmp_path / "seq.json"
    assert run(target_preparation, "validate-sequence",
               "--sequence", sequence, "--output", str(output)) == 0
    data = json.loads(output.read_text(encoding="utf-8"))
    assert data["valid"] is True
    assert data["length"] == len(sequence)


def test_validate_sequence_rejects_nonsense(tmp_path, capsys):
    code = run(target_preparation, "validate-sequence",
               "--sequence", "ACGTACGTACGTACGTACGTXXXX",
               "--output", str(tmp_path / "s.json"))
    assert code == InvalidInputError.exit_code
    assert "error:" in capsys.readouterr().err


def test_validate_sequence_rejects_an_overlong_one(tmp_path, capsys):
    code = run(target_preparation, "validate-sequence",
               "--sequence", "A" * 2000, "--output",
               str(tmp_path / "s.json"))
    assert code == InvalidInputError.exit_code
    assert "generate-af2-script" in capsys.readouterr().err


def test_generate_af2_script_is_executable_text(tmp_path):
    output = tmp_path / "submit_af2.sh"
    assert run(target_preparation, "generate-af2-script",
               "--fasta", "target.fasta", "--method", "colabfold",
               "--output", str(output)) == 0
    text = output.read_text(encoding="utf-8")
    assert text.startswith("#!/bin/bash")
    assert "colabfold_batch" in text
    assert b"\r" not in output.read_bytes()


def test_prepare_receptor_requires_a_source(tmp_path):
    with pytest.raises(SystemExit):
        target_preparation.main(["prepare-receptor",
                                 "--output", str(tmp_path / "r.json")])


def test_prepare_receptor_rejects_a_malformed_pdb_id(tmp_path, capsys):
    code = run(target_preparation, "prepare-receptor", "--pdb-id", "NOTANID",
               "--output", str(tmp_path / "r.json"))
    assert code == InvalidInputError.exit_code
    assert "rcsb.org" in capsys.readouterr().err


def test_prepare_receptor_rejects_a_gene_name_as_an_accession(tmp_path,
                                                              capsys):
    code = run(target_preparation, "prepare-receptor", "--uniprot-id", "EGFR",
               "--output", str(tmp_path / "r.json"))
    assert code == InvalidInputError.exit_code
    assert "uniprot.org" in capsys.readouterr().err


def test_prepare_receptor_from_a_local_file(pdb_1m17, tmp_path):
    output = tmp_path / "receptor.json"
    assert run(target_preparation, "prepare-receptor",
               "--input-file", str(pdb_1m17), "--chain", "A",
               "--output", str(output)) == 0
    data = json.loads(output.read_text(encoding="utf-8"))
    assert data["cleaning"]["atoms_kept"] > 1000
    assert "AQ4" in {e["resname"] for e in
                     data["cleaning"]["heterogens_removed"]}
    assert data["binding_site"]["reference"] == "AQ4:A:999"


# ---------------------------------------------------------------------------
# md-simulation
# ---------------------------------------------------------------------------


def test_gmxapi_setup_writes_the_workflow_and_mdp_files(tmp_path, capsys):
    output = tmp_path / "md" / "run_md.py"
    assert run(md_simulation, "gmxapi-setup", "--pdb", "protein_clean.pdb",
               "--ff", "charmm36", "--temperature", "310",
               "--production-ns", "50", "--output", str(output)) == 0
    assert output.is_file()
    for name in ("ions.mdp", "em.mdp", "nvt.mdp", "npt.mdp", "md.mdp",
                 "setup_md.sh"):
        assert (output.parent / name).is_file(), name
    # The requested temperature must reach the file grompp reads.
    assert "ref_t" in (output.parent / "md.mdp").read_text(encoding="utf-8")
    assert "310" in (output.parent / "md.mdp").read_text(encoding="utf-8")
    assert "310" in capsys.readouterr().out


def test_gmxapi_setup_zero_production_omits_the_production_input(tmp_path):
    output = tmp_path / "md" / "run_md.py"
    run(md_simulation, "gmxapi-setup", "--pdb", "p.pdb",
        "--production-ns", "0", "--output", str(output))
    assert not (output.parent / "md.mdp").exists()
    assert (output.parent / "npt.mdp").is_file()


def test_gmxapi_validate_accepts_generated_files(tmp_path):
    workflow = tmp_path / "md" / "run_md.py"
    run(md_simulation, "gmxapi-setup", "--pdb", "p.pdb", "--ff", "charmm36",
        "--production-ns", "10", "--output", str(workflow))
    mdps = [str(p) for p in sorted(workflow.parent.glob("*.mdp"))]
    report = tmp_path / "validation.json"
    assert run(md_simulation, "gmxapi-validate", "--mdp", *mdps,
               "--ff", "charmm36", "--output", str(report)) == 0
    data = json.loads(report.read_text(encoding="utf-8"))
    assert data["all_valid"] is True
    assert data["n_invalid"] == 0


def test_gmxapi_validate_detects_a_force_field_mismatch(tmp_path, capsys):
    """AMBER settings under a CHARMM36 label must be caught."""
    workflow = tmp_path / "md" / "run_md.py"
    run(md_simulation, "gmxapi-setup", "--pdb", "p.pdb",
        "--ff", "amber99sb-ildn", "--production-ns", "10",
        "--output", str(workflow))
    report = tmp_path / "validation.json"
    code = run(md_simulation, "gmxapi-validate",
               "--mdp", str(workflow.parent / "md.mdp"),
               "--ff", "charmm36", "--output", str(report))
    assert code == 0  # the report is the deliverable
    data = json.loads(report.read_text(encoding="utf-8"))
    assert data["all_valid"] is False
    errors = " ".join(data["reports"][0]["errors"])
    assert "force-switch" in errors
    assert "DispCorr" in errors
    assert "ERROR" in capsys.readouterr().out


def test_generate_hpc_script_for_each_job_type(tmp_path):
    from agskills.hpc import JOB_TYPES
    for job_type in JOB_TYPES:
        output = tmp_path / f"{job_type}.sh"
        argv = ["generate-hpc-script", "--job-type", job_type,
                "--output", str(output)]
        if job_type == "reinvent":
            argv += ["--config", "run.toml"]
        if job_type == "vina-screen":
            argv += ["--array", "1-50", "--ngpu", "0"]
        if job_type in ("colabfold", "alphafold2"):
            argv += ["--fasta", "t.fasta"]
        assert run(md_simulation, *argv) == 0, job_type
        assert output.read_text(encoding="utf-8").startswith("#!/bin/bash")


def test_generate_hpc_script_refuses_a_bad_walltime(tmp_path, capsys):
    code = run(md_simulation, "generate-hpc-script", "--job-type",
               "gromacs-md", "--time", "24h",
               "--output", str(tmp_path / "x.sh"))
    assert code == UsageError.exit_code
    assert "HH:MM:SS" in capsys.readouterr().err


def test_generate_hpc_script_refuses_vina_without_an_array(tmp_path, capsys):
    code = run(md_simulation, "generate-hpc-script", "--job-type",
               "vina-screen", "--ngpu", "0",
               "--output", str(tmp_path / "x.sh"))
    assert code == UsageError.exit_code
    assert "--array" in capsys.readouterr().err


def test_unknown_force_field_is_rejected_by_the_parser(tmp_path):
    with pytest.raises(SystemExit) as excinfo:
        md_simulation.main(["gmxapi-setup", "--pdb", "p.pdb", "--ff",
                            "made_up", "--output", str(tmp_path / "x.py")])
    assert excinfo.value.code == UsageError.exit_code


# ---------------------------------------------------------------------------
# compound-synthesis
# ---------------------------------------------------------------------------


def test_check_setup_reports_without_failing(tmp_path, capsys):
    output = tmp_path / "setup.json"
    assert run(compound_synthesis, "check-setup",
               "--output", str(output)) == 0
    data = json.loads(output.read_text(encoding="utf-8"))
    assert "ready" in data
    assert "summary" in data
    assert "REINVENT 4 setup" in capsys.readouterr().out


def test_list_profiles_documents_every_profile(tmp_path):
    output = tmp_path / "profiles.json"
    assert run(compound_synthesis, "list-profiles",
               "--output", str(output)) == 0
    data = json.loads(output.read_text(encoding="utf-8"))
    assert "drug-like" in data["profiles"]
    assert "anti-tb" in data["profiles"]
    assert data["profiles"]["anti-tb"]["caveat"]
    assert set(data["run_modes"]) >= {"staged_learning", "sampling"}


def test_generate_config_writes_toml_and_metadata(tmp_path):
    output = tmp_path / "run.toml"
    assert run(compound_synthesis, "generate-config",
               "--mode", "staged_learning", "--generator", "reinvent",
               "--scoring-profile", "drug-like", "--num-steps", "100",
               "--device", "cpu", "--output", str(output)) == 0
    import tomllib
    parsed = tomllib.loads(output.read_text(encoding="utf-8"))
    assert parsed["run_type"] == "staged_learning"
    metadata = tmp_path / "run.meta.json"
    assert metadata.is_file()
    data = json.loads(metadata.read_text(encoding="utf-8"))
    assert data["settings"]["max_steps"] == 100
    assert data["next_steps"]


def test_generate_config_refuses_a_seeded_generator_without_seeds(tmp_path,
                                                                   capsys):
    code = run(compound_synthesis, "generate-config",
               "--mode", "staged_learning", "--generator", "libinvent",
               "--scoring-profile", "drug-like", "--device", "cpu",
               "--output", str(tmp_path / "run.toml"))
    assert code == UsageError.exit_code
    assert "prepare-seeds" in capsys.readouterr().err


def test_prepare_seeds_writes_a_smi_file_and_a_report(tmp_path, capsys):
    output = tmp_path / "seeds.smi"
    report = tmp_path / "seeds.json"
    assert run(compound_synthesis, "prepare-seeds",
               "--smiles", "CCO", "OCC", "NOT_A_MOLECULE",
               "--names", "a", "b", "c",
               "--report", str(report), "--output", str(output)) == 0
    assert output.read_text(encoding="utf-8").strip() == "CCO"
    data = json.loads(report.read_text(encoding="utf-8"))
    assert data["input_count"] == 3
    assert data["written"] == 1
    assert data["rejected"] == 2
    assert "rejected 2" in capsys.readouterr().out


def test_analyze_results_on_the_real_csv_format(reinvent_csv, tmp_path,
                                                 capsys):
    output = tmp_path / "analysis.json"
    assert run(compound_synthesis, "analyze-results",
               "--csv", str(reinvent_csv), "--top-n", "5",
               "--output", str(output)) == 0
    data = json.loads(output.read_text(encoding="utf-8"))
    assert data["unique_molecules"] == 13
    assert len(data["top_molecules"]) == 5
    out = capsys.readouterr().out
    assert "unique molecules" in out
    assert "scaffolds" in out


def test_analyze_results_rejects_a_missing_file(tmp_path, capsys):
    code = run(compound_synthesis, "analyze-results",
               "--csv", str(tmp_path / "absent.csv"),
               "--output", str(tmp_path / "a.json"))
    assert code == ResourceNotFoundError.exit_code
    assert "error:" in capsys.readouterr().err


def test_run_dry_run_reports_without_executing(tmp_path, capsys):
    """A dry run must surface a missing prior rather than starting work."""
    config = tmp_path / "run.toml"
    config.write_text(
        'run_type = "sampling"\ndevice = "cpu"\n\n[parameters]\n'
        'model_file = "missing.prior"\n', encoding="utf-8")
    code = run(compound_synthesis, "run", "--config", str(config),
               "--dry-run", "--output", str(tmp_path / "report.json"))
    assert code == ResourceNotFoundError.exit_code
    assert "missing.prior" in capsys.readouterr().err


def test_generate_hpc_script_for_reinvent(tmp_path):
    output = tmp_path / "submit.sh"
    assert run(compound_synthesis, "generate-hpc-script",
               "--config", "run.toml", "--output", str(output)) == 0
    text = output.read_text(encoding="utf-8")
    assert "reinvent" in text
    assert "#SBATCH" in text


# ---------------------------------------------------------------------------
# wizard
# ---------------------------------------------------------------------------


def test_wizard_describe_documents_the_stages(tmp_path, capsys):
    output = tmp_path / "describe.json"
    assert run(wizard, "describe", "--output", str(output)) == 0
    data = json.loads(output.read_text(encoding="utf-8"))
    assert set(data["stages"]) >= {"target", "sourcing", "triage", "docking",
                                   "generative", "selection", "dynamics"}
    assert data["sub_skills"]
    assert data["do_not_use_when"]
    assert "drug-discovery-wizard" in capsys.readouterr().out


def test_wizard_plan_reports_readiness_per_stage(tmp_path, capsys):
    output = tmp_path / "plan.json"
    assert run(wizard, "plan", "--pdb-id", "6HEZ",
               "--stages", "target", "triage", "dynamics",
               "--output", str(output)) == 0
    data = json.loads(output.read_text(encoding="utf-8"))
    assert [s["stage"] for s in data["stages"]] == ["target", "triage",
                                                     "dynamics"]
    for stage in data["stages"]:
        assert "ready" in stage
        assert stage["produces"]
    assert data["two_tier_rationale"]
    assert "Pipeline plan" in capsys.readouterr().out


def test_wizard_plan_rejects_an_unknown_stage(tmp_path):
    with pytest.raises(SystemExit):
        wizard.main(["plan", "--stages", "not-a-stage",
                     "--output", str(tmp_path / "p.json")])


def test_wizard_run_executes_local_stages_on_a_real_structure(pdb_6hez,
                                                               tmp_path,
                                                               capsys):
    """The orchestrator must run the target stage end to end."""
    output = tmp_path / "pipeline.json"
    code = run(wizard, "run", "--input-file", str(pdb_6hez), "--chain", "A",
               "--stages", "target", "--workdir", str(tmp_path / "work"),
               "--output", str(output))
    assert code == 0
    data = json.loads(output.read_text(encoding="utf-8"))
    assert data["stages"]["target"]["status"] == "ok"
    assert data["artefacts"]["cleaned_pdb"]
    assert data["artefacts"]["binding_site"]
    assert data["skill"] == "drug-discovery-wizard"
    assert "drug-discovery-wizard" in data["note"]


def test_wizard_run_prepares_cluster_work_rather_than_pretending(tmp_path):
    """GPU stages must produce submission scripts, not fake results."""
    output = tmp_path / "pipeline.json"
    assert run(wizard, "run", "--stages", "generative", "dynamics",
               "--input-file", "protein_clean.pdb",
               "--workdir", str(tmp_path / "work"), "--device", "cpu",
               "--output", str(output)) == 0
    data = json.loads(output.read_text(encoding="utf-8"))
    assert data["stages"]["generative"]["status"] == "prepared"
    assert data["stages"]["dynamics"]["status"] == "prepared"
    for key in ("reinvent_config", "reinvent_script", "md_script"):
        assert key in data["artefacts"], key
    from pathlib import Path
    assert Path(data["artefacts"]["md_script"]).is_file()


def test_wizard_run_skips_stages_whose_inputs_are_absent(tmp_path):
    """A skipped stage must say why, not fail the pipeline."""
    output = tmp_path / "pipeline.json"
    assert run(wizard, "run", "--stages", "triage", "docking",
               "--workdir", str(tmp_path / "work"),
               "--output", str(output)) == 0
    data = json.loads(output.read_text(encoding="utf-8"))
    assert data["stages"]["triage"]["status"] == "skipped"
    assert data["stages"]["triage"]["reason"]
    assert data["stages"]["docking"]["status"] == "skipped"


def test_wizard_mentions_itself_in_its_output(tmp_path):
    """The skill contract requires the wizard to identify itself."""
    output = tmp_path / "pipeline.json"
    run(wizard, "run", "--stages", "generative",
        "--workdir", str(tmp_path / "w"), "--device", "cpu",
        "--output", str(output))
    data = json.loads(output.read_text(encoding="utf-8"))
    assert "drug-discovery-wizard" in json.dumps(data)
