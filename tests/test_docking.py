"""AutoDock Vina docking.

Vina is not installable via pip on Windows (the ``vina`` module needs Boost
headers), so the binary-dependent path is marked ``requires_vina`` and
skipped by default. Everything that can be verified without the binary is
verified here deterministically:

* the score parser, against authentic Vina v1.2.5 stdout including its
  banner and progress bar;
* conformer generation, which must be seeded so a run is reproducible;
* the box and argument handling that previously crashed outright;
* filename sanitisation, which previously allowed a compound name to
  escape the output directory;
* every failure path, each of which must be reported per compound rather
  than aborting the batch.
"""

from __future__ import annotations

import pytest

from agskills.chem.smiles import CompoundRecord
from agskills.docking.vina import (
    dock_compounds,
    find_vina,
    parse_vina_log,
    prepare_ligand_3d,
)
from agskills.errors import MissingDependencyError, ResourceNotFoundError, UsageError
from agskills.targets.site import BindingSite, site_from_hetatm

from panel import EGFR_ACTIVES


def a_site() -> BindingSite:
    return BindingSite(
        center_x=13.87, center_y=-20.62, center_z=37.13,
        size_x=18.0, size_y=18.0, size_z=15.0,
        method="test fixture",
    )


# ---------------------------------------------------------------------------
# Score parsing, against real Vina output
# ---------------------------------------------------------------------------


def test_parses_all_poses_from_real_vina_output(vina_output):
    poses = parse_vina_log(vina_output)
    assert len(poses) == 9
    assert poses[0].mode == 1
    assert poses[0].affinity_kcal_mol == pytest.approx(-9.428)
    assert poses[0].rmsd_lb == pytest.approx(0.0)
    assert poses[-1].affinity_kcal_mol == pytest.approx(-7.221)


def test_poses_come_back_in_mode_order(vina_output):
    poses = parse_vina_log(vina_output)
    assert [p.mode for p in poses] == list(range(1, 10))
    # Vina orders by affinity, so mode 1 is the best.
    affinities = [p.affinity_kcal_mol for p in poses]
    assert affinities == sorted(affinities)


def test_banner_and_progress_bar_are_not_mistaken_for_scores(vina_output):
    """The old parser used a loose regex that matched non-result lines.

    Vina's output contains a citation block, a grid description and a
    progress bar of asterisks; none of it is a score.
    """
    assert "AutoDock Vina v1.2.5" in vina_output
    assert "*****" in vina_output
    assert "0%   10   20" in vina_output
    poses = parse_vina_log(vina_output)
    assert len(poses) == 9
    # No pose may carry an implausible affinity scraped from the header.
    assert all(-20.0 < p.affinity_kcal_mol < 5.0 for p in poses)


def test_output_with_no_result_table_yields_no_poses():
    text = "AutoDock Vina v1.2.5\nSomething went wrong.\n"
    assert parse_vina_log(text) == []


def test_parser_requires_all_four_columns():
    """A truncated row is not a score."""
    assert parse_vina_log("   1       -9.428\n") == []
    assert len(parse_vina_log("   1       -9.428          0          0\n")) == 1


def test_parser_accepts_positive_affinities():
    """A badly clashing pose can score positive; it must still parse."""
    poses = parse_vina_log("   1        2.500      0.000      0.000\n")
    assert poses[0].affinity_kcal_mol == pytest.approx(2.5)


# ---------------------------------------------------------------------------
# Conformer generation
# ---------------------------------------------------------------------------


def test_ligand_embedding_produces_3d_coordinates():
    """Vina needs coordinates; a SMILES string has none."""
    mol = prepare_ligand_3d(EGFR_ACTIVES[0], name="egfr")
    assert mol.GetNumConformers() == 1
    conformer = mol.GetConformer()
    positions = [conformer.GetAtomPosition(i)
                 for i in range(mol.GetNumAtoms())]
    assert any(abs(p.z) > 0.01 for p in positions), "the molecule is flat"


def test_hydrogens_are_added_before_embedding():
    mol = prepare_ligand_3d("CCO")
    assert mol.GetNumAtoms() > 3, "explicit hydrogens should be present"


def test_embedding_is_seeded_and_therefore_reproducible():
    """Unseeded embedding gives different coordinates on every run.

    Different coordinates mean different docking scores, so an unseeded
    run is not reproducible. The previous implementation left ETKDG
    unseeded.
    """
    def coordinates(seed: int):
        mol = prepare_ligand_3d(EGFR_ACTIVES[0], seed=seed)
        conformer = mol.GetConformer()
        return [(round(conformer.GetAtomPosition(i).x, 6),
                 round(conformer.GetAtomPosition(i).y, 6),
                 round(conformer.GetAtomPosition(i).z, 6))
                for i in range(mol.GetNumAtoms())]

    assert coordinates(42) == coordinates(42)
    assert coordinates(42) != coordinates(99)


def test_unparseable_ligand_is_refused_by_name():
    with pytest.raises(UsageError) as excinfo:
        prepare_ligand_3d("NOT_A_MOLECULE", name="bad_entry")
    assert "bad_entry" in str(excinfo.value)


@pytest.mark.parametrize("smiles", EGFR_ACTIVES)
def test_real_drug_like_molecules_all_embed(smiles):
    assert prepare_ligand_3d(smiles).GetNumConformers() == 1


def test_a_macrocycle_still_embeds():
    """Strained rings need the random-coordinate fallback."""
    macrocycle = "C1CCCCCCCCCCCCCCCCC1"
    assert prepare_ligand_3d(macrocycle).GetNumConformers() == 1


# ---------------------------------------------------------------------------
# Box and argument handling: the crash that made docking unusable
# ---------------------------------------------------------------------------


def test_binding_site_exposes_centre_and_size_as_tuples():
    """The old code read ``args.box_size`` while the parser defined
    ``--size``, so any ``--center`` invocation raised AttributeError."""
    site = a_site()
    assert site.center == (13.87, -20.62, 37.13)
    assert site.size == (18.0, 18.0, 15.0)
    assert site.volume == pytest.approx(18.0 * 18.0 * 15.0)


def test_site_from_a_real_structure_feeds_docking_directly(pdb_1m17):
    """The object produced by site detection is the one docking consumes."""
    site = site_from_hetatm(pdb_1m17)
    assert len(site.center) == 3
    assert len(site.size) == 3
    assert all(value > 0 for value in site.size)


def test_missing_receptor_is_reported_before_any_work(tmp_path):
    with pytest.raises(ResourceNotFoundError):
        dock_compounds([CompoundRecord("a", "CCO")],
                       tmp_path / "absent.pdbqt", a_site(),
                       output_dir=tmp_path)


def test_a_pdb_receptor_is_rejected_with_guidance(tmp_path):
    """Vina needs PDBQT; handing it a PDB fails deep inside the binary."""
    receptor = tmp_path / "receptor.pdb"
    receptor.write_text("ATOM\n", encoding="utf-8")
    with pytest.raises(UsageError) as excinfo:
        dock_compounds([CompoundRecord("a", "CCO")], receptor, a_site(),
                       output_dir=tmp_path)
    message = str(excinfo.value)
    assert "PDBQT" in message
    assert "prepare-receptor" in message


def test_an_empty_compound_list_is_refused(tmp_path):
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("ATOM\n", encoding="utf-8")
    with pytest.raises(UsageError):
        dock_compounds([], receptor, a_site(), output_dir=tmp_path)


def test_absent_vina_is_reported_with_install_advice(tmp_path, monkeypatch):
    """A missing binary must say how to get it, not fail obscurely."""
    monkeypatch.setattr("agskills.docking.vina.find_vina",
                        lambda: ("none", ""))
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("ATOM\n", encoding="utf-8")
    with pytest.raises(MissingDependencyError) as excinfo:
        dock_compounds([CompoundRecord("a", "CCO")], receptor, a_site(),
                       output_dir=tmp_path)
    message = str(excinfo.value)
    assert "Vina" in message
    assert "conda" in message or "github" in message


def test_find_vina_reports_a_kind_and_detail():
    kind, detail = find_vina()
    assert kind in ("python", "cli", "none")
    if kind == "none":
        assert detail == ""


# ---------------------------------------------------------------------------
# Filename safety
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name,expected", [
    ("ZINC00001/trans", "ZINC00001_trans"),
    ("compound:1", "compound_1"),
    ("a b c", "a_b_c"),
    # Leading dots are stripped too, so a name cannot produce a hidden
    # file or a parent-directory reference.
    ("../escape", "_escape"),
    (".hidden", "hidden"),
    ("CON", "CON_"),
    ("", "compound"),
    ("   ", "compound"),
])
def test_compound_names_are_made_filesystem_safe(name, expected):
    """Docking writes one pose file per compound, named after it.

    Catalogue names routinely contain ``/`` and ``:``. Interpolating them
    into a path either escapes the output directory or fails outright on
    Windows.
    """
    from agskills.io_utils import safe_stem
    assert safe_stem(name) == expected


def test_safe_stem_is_length_limited():
    from agskills.io_utils import safe_stem
    assert len(safe_stem("x" * 500)) <= 64


# ---------------------------------------------------------------------------
# Reporting contract
# ---------------------------------------------------------------------------


def test_docking_result_serialises_and_explains_itself():
    """A None affinity must not be formatted as a float.

    The previous code ran ``print(f"{score:.1f} kcal/mol")`` unconditionally
    and raised TypeError whenever no affinity could be parsed - after the
    docking work was already done.
    """
    from agskills.docking.vina import DockingResult
    import json

    failed = DockingResult(name="x", smiles="CCO", status="no_pose",
                           error="Vina returned no scored pose")
    data = failed.as_dict()
    json.dumps(data)
    assert data["best_affinity_kcal_mol"] is None
    assert data["error"]


def test_interpretation_warns_against_overreading_scores(tmp_path,
                                                          monkeypatch):
    """A Vina score is triage, not a binding free energy.

    Reporting a ranking without that caveat invites over-interpretation of
    sub-kcal differences.
    """
    # Build a report without needing the binary by stubbing the engine.
    from agskills.docking import vina as vina_module

    monkeypatch.setattr(vina_module, "find_vina", lambda: ("cli", "/fake/vina"))

    def fake_dock(*args, **kwargs):
        raise RuntimeError("stubbed")

    monkeypatch.setattr(vina_module, "_dock_one_cli", fake_dock)
    monkeypatch.setattr(vina_module, "ligand_to_pdbqt",
                        lambda mol, name="ligand": "FAKE PDBQT")

    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("ATOM\n", encoding="utf-8")
    report = dock_compounds([CompoundRecord("a", "CCO")], receptor, a_site(),
                            output_dir=tmp_path)

    assert "kcal/mol" in report["interpretation"]
    assert "1 kcal/mol" in report["interpretation"]
    assert "not a binding free energy" in report["interpretation"]
    assert report["citation"]
    # The failure was recorded per compound, not raised.
    assert report["failed"] == 1
    assert report["docked"] == 0
    assert report["ranked_results"][0]["status"] == "docking_failed"


def test_report_records_the_parameters_used(tmp_path, monkeypatch):
    """A docking run must be reproducible from its own report."""
    from agskills.docking import vina as vina_module
    monkeypatch.setattr(vina_module, "find_vina", lambda: ("cli", "/fake/vina"))
    monkeypatch.setattr(vina_module, "ligand_to_pdbqt",
                        lambda mol, name="ligand": "FAKE")
    monkeypatch.setattr(vina_module, "_dock_one_cli",
                        lambda *a, **k: (_ for _ in ()).throw(
                            RuntimeError("stub")))
    receptor = tmp_path / "r.pdbqt"
    receptor.write_text("ATOM\n", encoding="utf-8")
    report = dock_compounds([CompoundRecord("a", "CCO")], receptor, a_site(),
                            output_dir=tmp_path, exhaustiveness=64,
                            n_poses=5, seed=7)
    parameters = report["parameters"]
    assert parameters["exhaustiveness"] == 64
    assert parameters["n_poses"] == 5
    assert parameters["seed"] == 7
    assert report["binding_site"]["center_x"] == pytest.approx(13.87)


def test_every_compound_appears_in_the_report(tmp_path, monkeypatch):
    """A batch must account for all its inputs, including failures."""
    from agskills.docking import vina as vina_module
    monkeypatch.setattr(vina_module, "find_vina", lambda: ("cli", "/fake/vina"))
    monkeypatch.setattr(vina_module, "ligand_to_pdbqt",
                        lambda mol, name="ligand": "FAKE")
    monkeypatch.setattr(vina_module, "_dock_one_cli",
                        lambda *a, **k: (_ for _ in ()).throw(
                            RuntimeError("stub")))
    receptor = tmp_path / "r.pdbqt"
    receptor.write_text("ATOM\n", encoding="utf-8")

    records = [CompoundRecord(f"c{i}", s) for i, s in enumerate(EGFR_ACTIVES)]
    records.append(CompoundRecord("broken", "NOT_A_MOLECULE"))
    report = dock_compounds(records, receptor, a_site(), output_dir=tmp_path)

    assert report["total_compounds"] == len(records)
    assert report["docked"] + report["failed"] == len(records)
    assert len(report["ranked_results"]) == len(records)
    # The unparseable one is reported as a preparation failure, not a crash.
    broken = next(r for r in report["ranked_results"]
                  if r["name"] == "broken")
    assert broken["status"] == "preparation_failed"


# ---------------------------------------------------------------------------
# Real docking, opt-in
# ---------------------------------------------------------------------------


@pytest.mark.requires_vina
@pytest.mark.slow
def test_real_docking_against_the_dprE1_receptor(tmp_path, pdb_6hez):
    """End-to-end docking into the real 6HEZ inhibitor site.

    Skipped unless a Vina binary is present. When it runs, it checks that
    affinities are physically plausible and that the known inhibitor
    chemotype scores better than a trivial solvent molecule.
    """
    from agskills.targets.clean import clean_structure
    from agskills.targets.pdbqt import receptor_to_pdbqt

    cleaned = tmp_path / "6HEZ_clean.pdb"
    clean_structure(pdb_6hez, cleaned, chain="A")
    site = site_from_hetatm(pdb_6hez)
    conversion = receptor_to_pdbqt(cleaned, site=site)
    if not conversion.ok:
        pytest.skip(f"no PDBQT converter: {conversion.message}")

    records = [
        CompoundRecord("ethanol", "CCO"),
        CompoundRecord("inhibitor_like", EGFR_ACTIVES[0]),
    ]
    report = dock_compounds(records, conversion.output_path, site,
                            output_dir=tmp_path / "poses",
                            exhaustiveness=8, n_poses=3)

    assert report["docked"] == 2, report["ranked_results"]
    for entry in report["ranked_results"]:
        affinity = entry["best_affinity_kcal_mol"]
        assert -20.0 < affinity < 0.0, f"{entry['name']}: {affinity}"

    scores = {e["name"]: e["best_affinity_kcal_mol"]
              for e in report["ranked_results"]}
    assert scores["inhibitor_like"] < scores["ethanol"], (
        "a drug-sized ligand should outscore ethanol"
    )
