"""PDBQT conversion: discovery, the two-pass strategy, and the in-box check.

The branch this replaces imported ``PDBQTWriterLegacy`` and then, on
success, set a note reading "Open Babel not available. PDBQT conversion
requires obabel" and produced no PDBQT at all. The import result was never
used.

The more interesting behaviour is new. Meeko matches every residue against
a chemical template, and residues with missing side-chain atoms fail that
match - 6HEZ has 16 such residues in chain A alone, so strict conversion
simply fails on a routine crystal structure. The retry is reported, and so
is the part that actually matters: whether any approximated residue lies
**inside** the docking box.
"""

from __future__ import annotations

import os
import sys

import pytest

from agskills.errors import MissingDependencyError, ResourceNotFoundError
from agskills.targets.clean import clean_structure
from agskills.targets.pdbqt import (
    ConversionResult,
    _parse_bad_residues,
    _residues_in_box,
    available_converters,
    find_executable,
    ligand_to_pdbqt,
    receptor_to_pdbqt,
)
from agskills.targets.site import BindingSite, site_from_hetatm

# Real mk_prepare_receptor output for 6HEZ chain A.
MEEKO_FAILURE = """\
Error: Creation of data structure for receptor failed.

Details:
- Template matching failed for: ['A:7', 'A:18', 'A:37', 'A:41', 'A:254', \
'A:259', 'A:263', 'A:266', 'A:267', 'A:294', 'A:299', 'A:317', 'A:334', \
'A:349', 'A:393', 'A:421']
These residues can be ignored with option --delete_bad_res or \
--delete_bad_res_from_box_radius.
"""


# ---------------------------------------------------------------------------
# Executable discovery
# ---------------------------------------------------------------------------


def test_finds_a_console_script_inside_this_environment():
    """``shutil.which`` consults PATH only.

    A console script installed into a virtual environment sits next to
    ``sys.executable``, and is frequently absent from PATH - which is how
    most tooling invokes the interpreter. Missing it made the converter
    report "none available" on a machine where Meeko was installed and
    working.
    """
    # The interpreter running these tests is itself discoverable this way.
    scripts = os.path.dirname(sys.executable)
    assert os.path.isdir(scripts)
    found = find_executable("pytest")
    if found is None:
        pytest.skip("pytest is not installed as a console script here")
    assert os.path.isfile(found)


def test_an_absent_executable_returns_none_rather_than_raising():
    assert find_executable("definitely-not-a-real-executable-xyz") is None


def test_available_converters_answers_without_raising():
    converters = available_converters()
    assert isinstance(converters, list)
    assert all(isinstance(name, str) for name in converters)


# ---------------------------------------------------------------------------
# Parsing Meeko's complaint
# ---------------------------------------------------------------------------


def test_unmatched_residues_are_extracted_from_real_output():
    residues = _parse_bad_residues(MEEKO_FAILURE)
    assert len(residues) == 16
    assert residues[0] == "A:7"
    assert "A:317" in residues
    # Order is preserved and duplicates removed.
    assert len(residues) == len(set(residues))


def test_clean_output_yields_no_unmatched_residues():
    assert _parse_bad_residues("Files written:\nreceptor.pdbqt\n") == []


def test_repeated_reports_are_deduplicated():
    residues = _parse_bad_residues(MEEKO_FAILURE + MEEKO_FAILURE)
    assert len(residues) == 16


# ---------------------------------------------------------------------------
# The in-box check: the part that matters scientifically
# ---------------------------------------------------------------------------


def test_identifies_which_approximated_residues_lie_inside_the_box(pdb_6hez,
                                                                    tmp_path):
    """A defect inside the docking box is the one that invalidates a score.

    An approximated side chain elsewhere in the protein is harmless; one
    lining the pocket makes affinities for that site unreliable. For 6HEZ
    chain A, exactly one of the 16 unmatched residues falls inside the
    inhibitor box.
    """
    cleaned = tmp_path / "6HEZ_clean.pdb"
    clean_structure(pdb_6hez, cleaned, chain="A")
    site = site_from_hetatm(pdb_6hez)

    residues = _parse_bad_residues(MEEKO_FAILURE)
    inside = _residues_in_box(residues, cleaned, site)

    assert inside, "expected at least one affected residue inside the box"
    assert set(inside) <= set(residues)
    assert len(inside) < len(residues), (
        "not every affected residue should be in the box, or the check is "
        "not discriminating"
    )
    assert "A:317" in inside


def test_no_box_means_no_in_box_claim(pdb_6hez, tmp_path):
    cleaned = tmp_path / "clean.pdb"
    clean_structure(pdb_6hez, cleaned, chain="A")
    assert _residues_in_box(_parse_bad_residues(MEEKO_FAILURE), cleaned,
                            None) == []


def test_an_empty_residue_list_is_handled(pdb_1m17):
    site = site_from_hetatm(pdb_1m17)
    assert _residues_in_box([], pdb_1m17, site) == []


def test_malformed_residue_selectors_are_skipped(pdb_1m17):
    site = site_from_hetatm(pdb_1m17)
    assert _residues_in_box(["not-a-residue", "A:abc"], pdb_1m17, site) == []


def test_a_residue_far_from_the_box_is_not_reported(pdb_1m17):
    """The check must be geometric, not a pass-through."""
    tiny_box = BindingSite(
        center_x=9999.0, center_y=9999.0, center_z=9999.0,
        size_x=5.0, size_y=5.0, size_z=5.0, method="deliberately elsewhere",
    )
    assert _residues_in_box(["A:766", "A:769"], pdb_1m17, tiny_box) == []


# ---------------------------------------------------------------------------
# Conversion behaviour
# ---------------------------------------------------------------------------


def test_a_missing_input_raises(tmp_path):
    with pytest.raises(ResourceNotFoundError):
        receptor_to_pdbqt(tmp_path / "absent.pdb")


def test_no_converter_reports_advice_rather_than_succeeding(tmp_path,
                                                             monkeypatch):
    """A failure must not be presented as a result.

    The branch this replaces returned a success note while producing no
    PDBQT file.
    """
    monkeypatch.setattr("agskills.targets.pdbqt.find_executable",
                        lambda name: None)
    source = tmp_path / "receptor.pdb"
    source.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00"
        "           C\n", encoding="utf-8")

    result = receptor_to_pdbqt(source)
    assert result.ok is False
    assert result.output_path is None
    assert result.converter == "none"
    assert "Auto-installation" in result.message
    assert str(source) in result.message


def test_a_successful_conversion_is_reported_with_its_converter(tmp_path,
                                                                 monkeypatch):
    """Simulate a converter that works, without needing one installed."""
    source = tmp_path / "receptor.pdb"
    source.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00"
        "           C\n", encoding="utf-8")
    destination = source.with_suffix(".pdbqt")

    monkeypatch.setattr(
        "agskills.targets.pdbqt.find_executable",
        lambda name: "/fake/mk_prepare_receptor"
        if name == "mk_prepare_receptor" else None)

    class _Completed:
        returncode = 0
        stdout = "Files written:\nreceptor.pdbqt\n"
        stderr = ""

    def fake_run(command, timeout=900.0):
        destination.write_text("ATOM      1  C   ALA A   1   0.0 0.0 0.0\n",
                               encoding="utf-8")
        return _Completed()

    monkeypatch.setattr("agskills.targets.pdbqt._run", fake_run)

    result = receptor_to_pdbqt(source)
    assert result.ok is True
    assert result.converter == "meeko"
    assert result.strict is True
    assert result.compromised_residues == []
    assert destination.is_file()


def test_a_template_failure_triggers_a_reported_permissive_retry(tmp_path,
                                                                  monkeypatch):
    """Strict then permissive, with the compromise itemised."""
    source = tmp_path / "receptor.pdb"
    source.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00"
        "           C\n", encoding="utf-8")
    destination = source.with_suffix(".pdbqt")

    monkeypatch.setattr(
        "agskills.targets.pdbqt.find_executable",
        lambda name: "/fake/mk_prepare_receptor"
        if name == "mk_prepare_receptor" else None)

    attempts: list[list[str]] = []

    class _Completed:
        returncode = 0
        stderr = ""

        def __init__(self, stdout):
            self.stdout = stdout

    def fake_run(command, timeout=900.0):
        attempts.append(command)
        if "-a" in command:
            destination.write_text("ATOM      1  C   ALA A   1   0 0 0\n",
                                   encoding="utf-8")
            return _Completed("Ignored due to allow_bad_res.\n")
        return _Completed(MEEKO_FAILURE)

    monkeypatch.setattr("agskills.targets.pdbqt._run", fake_run)

    result = receptor_to_pdbqt(source)
    assert result.ok is True
    assert result.strict is False, "a permissive pass must say so"
    assert len(result.compromised_residues) == 16
    assert len(attempts) == 2, "strict first, then permissive"
    assert "-a" not in attempts[0]
    assert "-a" in attempts[1]
    assert any("did not match a chemical template" in w
               for w in result.warnings)
    assert any("missing from the crystal structure" in w
               for w in result.warnings)


def test_permissive_retry_can_be_refused(tmp_path, monkeypatch):
    """A caller may require a strictly-correct receptor."""
    source = tmp_path / "receptor.pdb"
    source.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00"
        "           C\n", encoding="utf-8")

    monkeypatch.setattr(
        "agskills.targets.pdbqt.find_executable",
        lambda name: "/fake/mk_prepare_receptor"
        if name == "mk_prepare_receptor" else None)

    class _Completed:
        returncode = 0
        stdout = MEEKO_FAILURE
        stderr = ""

    monkeypatch.setattr("agskills.targets.pdbqt._run",
                        lambda command, timeout=900.0: _Completed())

    result = receptor_to_pdbqt(source, allow_incomplete_residues=False)
    assert result.ok is False
    assert any("permissive mode was disabled" in w for w in result.warnings)


def test_conversion_result_serialises():
    import json
    result = ConversionResult(ok=False, output_path=None, converter="none",
                              message="nothing available")
    json.dumps(result.as_dict())


def test_ligand_conversion_without_meeko_names_the_install(monkeypatch):
    """A missing optional dependency must explain how to get it."""
    import builtins
    real_import = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name == "meeko":
            raise ImportError("no meeko")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    # Also block auto-install so the error surfaces.
    monkeypatch.setattr("agskills.targets.pdbqt.ensure_meeko", lambda: False)
    with pytest.raises(MissingDependencyError) as excinfo:
        ligand_to_pdbqt(None, name="test")
    message = str(excinfo.value)
    assert "meeko" in message
    assert "Auto-install" in message


# ---------------------------------------------------------------------------
# Real conversion, when a converter is present
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_real_receptor_conversion_on_6hez(pdb_6hez, tmp_path):
    """End-to-end conversion of a real structure, if a converter exists."""
    if not available_converters():
        pytest.skip("no PDBQT converter is installed")

    cleaned = tmp_path / "6HEZ_clean.pdb"
    clean_structure(pdb_6hez, cleaned, chain="A")
    site = site_from_hetatm(pdb_6hez)
    result = receptor_to_pdbqt(cleaned, site=site)

    if not result.ok:
        pytest.skip(f"conversion unavailable: {result.message}")

    from pathlib import Path
    output = Path(result.output_path)
    assert output.is_file()
    text = output.read_text(encoding="utf-8")
    # A PDBQT carries a partial charge and an AutoDock atom type per atom.
    atom_lines = [line for line in text.splitlines()
                  if line.startswith(("ATOM", "HETATM"))]
    assert len(atom_lines) > 1000
    assert any(line.rstrip().endswith(("A", "C", "N", "OA", "NA", "SA", "HD"))
               for line in atom_lines)

    # 6HEZ has incomplete side chains, so this is the permissive path, and
    # the report must say which residues were affected and where.
    if not result.strict:
        assert result.compromised_residues
        assert any("template" in w for w in result.warnings)


@pytest.mark.slow
def test_real_ligand_conversion(tmp_path):
    """Convert a real drug-like molecule to a ligand PDBQT."""
    try:
        import meeko  # noqa: F401
    except ImportError:
        pytest.skip("meeko is not installed")

    from agskills.docking.vina import prepare_ligand_3d
    mol = prepare_ligand_3d("COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC", name="egfr")
    pdbqt = ligand_to_pdbqt(mol, name="egfr")
    assert "ROOT" in pdbqt
    assert "TORSDOF" in pdbqt
    atom_lines = [line for line in pdbqt.splitlines()
                  if line.startswith("ATOM")]
    assert len(atom_lines) > 20
