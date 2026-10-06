"""Receptor cleaning on real structures.

The defect under test: the previous filter removed only residues flagged
``"W"`` by BioPython, which is water. Every other heterogen - the
co-crystallised ligand, cofactors, cryoprotectants, buffer ions - was
retained, despite documentation claiming they were removed. Docking into a
receptor that still contains its own ligand searches an occupied box.
"""

from __future__ import annotations

import pytest

from agskills.errors import InvalidInputError
from agskills.targets.clean import (
    BIOLOGICAL_IONS,
    COMMON_ADDITIVES,
    SOLVENT_RESNAMES,
    clean_structure,
)
from agskills.targets.site import parse_atoms


def _resnames(path):
    return {a.resname for a in parse_atoms(path.read_text(encoding="utf-8"))
            if a.is_hetatm}


# ---------------------------------------------------------------------------
# The core defect
# ---------------------------------------------------------------------------


def test_the_cocrystallised_ligand_is_removed(pdb_1m17, tmp_path):
    """Erlotinib must not survive into the docking receptor.

    This is the regression test. 1M17 contains erlotinib as residue AQ4;
    leaving it in place means Vina searches a box that is already occupied.
    """
    assert "AQ4" in _resnames(pdb_1m17), "fixture should contain the ligand"

    output = tmp_path / "clean.pdb"
    report = clean_structure(pdb_1m17, output, chain="A")

    assert "AQ4" not in _resnames(output)
    removed = {entry["resname"] for entry in report.heterogens_removed}
    assert "AQ4" in removed


def test_waters_are_removed_and_counted(pdb_1m17, tmp_path):
    output = tmp_path / "clean.pdb"
    report = clean_structure(pdb_1m17, output, chain="A")
    assert report.waters_removed > 0
    assert not (_resnames(output) & SOLVENT_RESNAMES)


def test_every_removal_states_a_reason(pdb_6hez, tmp_path):
    """A user must be able to see what was discarded and why."""
    output = tmp_path / "clean.pdb"
    report = clean_structure(pdb_6hez, output, chain="A")
    assert report.heterogens_removed
    for entry in report.heterogens_removed:
        assert entry["reason"], entry
        assert entry["count"] >= 1
    reasons = {entry["reason"] for entry in report.heterogens_removed}
    assert any("water" in r or "heterogen" in r or "additive" in r
               for r in reasons)


# ---------------------------------------------------------------------------
# Chain selection
# ---------------------------------------------------------------------------


def test_chain_selection_keeps_only_that_chain(pdb_6hez, tmp_path):
    output = tmp_path / "chain_a.pdb"
    report = clean_structure(pdb_6hez, output, chain="A")
    assert report.chains_kept == ["A"]
    assert len(report.chains_in_input) >= 2
    chains = {a.chain for a in parse_atoms(output.read_text(encoding="utf-8"))}
    assert chains == {"A"}


def test_omitting_the_chain_keeps_them_all(pdb_6hez, tmp_path):
    output = tmp_path / "all.pdb"
    report = clean_structure(pdb_6hez, output)
    assert len(report.chains_kept) >= 2


def test_requesting_an_absent_chain_lists_the_real_ones(pdb_1m17, tmp_path):
    with pytest.raises(InvalidInputError) as excinfo:
        clean_structure(pdb_1m17, tmp_path / "x.pdb", chain="Z")
    message = str(excinfo.value)
    assert "Z" in message
    assert "A" in message
    assert "--chain" in message


def test_selecting_one_chain_reduces_the_atom_count(pdb_6hez, tmp_path):
    one = clean_structure(pdb_6hez, tmp_path / "a.pdb", chain="A")
    both = clean_structure(pdb_6hez, tmp_path / "ab.pdb")
    assert one.atoms_kept < both.atoms_kept


# ---------------------------------------------------------------------------
# Configurable retention
# ---------------------------------------------------------------------------


def test_a_cofactor_can_be_retained_deliberately(pdb_6hez, tmp_path):
    """Keeping FAD is a legitimate choice for a flavoenzyme.

    What matters is that it is a *choice*, reported either way, rather than
    an accident of the filter.
    """
    without = tmp_path / "no_fad.pdb"
    clean_structure(pdb_6hez, without, chain="A")
    assert "FAD" not in _resnames(without)

    with_fad = tmp_path / "with_fad.pdb"
    report = clean_structure(pdb_6hez, with_fad, chain="A",
                             keep_ligands=["FAD"])
    assert "FAD" in _resnames(with_fad)
    kept = {entry["resname"] for entry in report.heterogens_kept}
    assert "FAD" in kept
    reasons = {entry["resname"]: entry["reason"]
               for entry in report.heterogens_kept}
    assert "--keep-ligand" in reasons["FAD"]


def test_waters_can_be_retained(pdb_1m17, tmp_path):
    output = tmp_path / "wet.pdb"
    report = clean_structure(pdb_1m17, output, chain="A", keep_waters=True)
    assert report.waters_removed == 0
    assert _resnames(output) & SOLVENT_RESNAMES


def test_metals_are_kept_by_default_but_can_be_removed(tmp_path):
    """A catalytic zinc changes the chemistry of its site.

    Removing it silently would alter the pocket a user is docking into, so
    metals are retained unless explicitly dropped.
    """
    source = tmp_path / "metalloprotein.pdb"
    source.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00"
        "           N\n"
        "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00 20.00"
        "           C\n"
        "ATOM      3  C   ALA A   1       2.000   0.000   0.000  1.00 20.00"
        "           C\n"
        "HETATM    4 ZN    ZN A 100       5.000   0.000   0.000  1.00 20.00"
        "          ZN\n"
        "HETATM    5  O   HOH A 200      10.000   0.000   0.000  1.00 20.00"
        "           O\n",
        encoding="utf-8")

    kept = tmp_path / "with_zn.pdb"
    report = clean_structure(source, kept)
    assert "ZN" in _resnames(kept)
    assert "ZN" in {e["resname"] for e in report.heterogens_kept}

    dropped = tmp_path / "no_zn.pdb"
    report = clean_structure(source, dropped, keep_metals=False)
    assert "ZN" not in _resnames(dropped)
    assert "ZN" in BIOLOGICAL_IONS


def test_buffer_ions_are_removed_even_when_metals_are_kept(tmp_path):
    """Sodium and chloride are counter-ions, not catalytic metals."""
    source = tmp_path / "salty.pdb"
    source.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00"
        "           C\n"
        "HETATM    2 NA    NA A 100       5.000   0.000   0.000  1.00 20.00"
        "          NA\n"
        "HETATM    3 CL    CL A 101       6.000   0.000   0.000  1.00 20.00"
        "          CL\n"
        "HETATM    4  C1  GOL A 102       7.000   0.000   0.000  1.00 20.00"
        "           C\n",
        encoding="utf-8")
    output = tmp_path / "clean.pdb"
    clean_structure(source, output, keep_metals=True)
    remaining = _resnames(output)
    assert not remaining & {"NA", "CL", "GOL"}
    assert {"NA", "CL", "GOL"} <= COMMON_ADDITIVES


def test_hydrogens_can_be_stripped(tmp_path):
    """Docking tools add their own hydrogens consistently."""
    source = tmp_path / "with_h.pdb"
    source.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00"
        "           N\n"
        "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00 20.00"
        "           C\n"
        "ATOM      3  HA  ALA A   1       1.500   0.500   0.000  1.00 20.00"
        "           H\n",
        encoding="utf-8")
    output = tmp_path / "no_h.pdb"
    report = clean_structure(source, output, remove_hydrogens=True)
    assert report.hydrogens_removed == 1
    assert report.atoms_kept == 2


# ---------------------------------------------------------------------------
# Structural edge cases
# ---------------------------------------------------------------------------


def test_only_the_primary_altloc_is_kept(pdb_6hez, tmp_path):
    """A receptor must have one unambiguous conformation.

    Keeping both alternate locations gives two atoms at nearly the same
    coordinates, which most docking tools either reject or treat as a
    clash.
    """
    output = tmp_path / "clean.pdb"
    report = clean_structure(pdb_6hez, output, chain="A")
    if report.altloc_atoms_dropped:
        assert report.altloc_atoms_dropped > 0
    # No (chain, resseq, atom name) may appear twice in the output.
    seen = set()
    duplicates = []
    for atom in parse_atoms(output.read_text(encoding="utf-8")):
        key = (atom.chain, atom.resseq, atom.icode, atom.name)
        if key in seen:
            duplicates.append(key)
        seen.add(key)
    assert not duplicates, f"duplicate atoms retained: {duplicates[:5]}"


def test_atom_and_residue_counts_are_measured_from_the_output(pdb_1m17,
                                                               tmp_path):
    """Counts must describe the file written, not the filter's intent."""
    output = tmp_path / "clean.pdb"
    report = clean_structure(pdb_1m17, output, chain="A")
    actual = len(parse_atoms(output.read_text(encoding="utf-8")))
    assert report.atoms_kept == actual
    assert report.residues_kept > 0
    assert report.residues_kept < report.atoms_kept


def test_cleaning_everything_away_raises_rather_than_writing_an_empty_file(
        tmp_path):
    """An empty receptor would fail confusingly several steps later."""
    source = tmp_path / "ligand_only.pdb"
    source.write_text(
        "HETATM    1  C1  GOL A 100       0.000   0.000   0.000  1.00 20.00"
        "           C\n"
        "HETATM    2  O   HOH A 200       5.000   0.000   0.000  1.00 20.00"
        "           O\n",
        encoding="utf-8")
    with pytest.raises(InvalidInputError) as excinfo:
        clean_structure(source, tmp_path / "empty.pdb")
    assert "--keep-ligand" in str(excinfo.value)


def test_multiple_models_are_reported_and_only_one_kept(tmp_path):
    source = tmp_path / "ensemble.pdb"
    source.write_text(
        "MODEL        1\n"
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00"
        "           C\n"
        "ATOM      2  CA  ALA A   2       1.000   0.000   0.000  1.00 20.00"
        "           C\n"
        "ENDMDL\n"
        "MODEL        2\n"
        "ATOM      3  CA  ALA A   1       9.000   0.000   0.000  1.00 20.00"
        "           C\n"
        "ATOM      4  CA  ALA A   2      10.000   0.000   0.000  1.00 20.00"
        "           C\n"
        "ENDMDL\n",
        encoding="utf-8")
    output = tmp_path / "clean.pdb"
    report = clean_structure(source, output)
    assert any("model" in warning.lower() for warning in report.warnings)
    assert report.atoms_kept == 2


def test_report_serialises(pdb_1m17, tmp_path):
    import json
    report = clean_structure(pdb_1m17, tmp_path / "clean.pdb", chain="A")
    json.dumps(report.as_dict())


def test_missing_input_raises(tmp_path):
    from agskills.errors import ResourceNotFoundError
    with pytest.raises(ResourceNotFoundError):
        clean_structure(tmp_path / "nope.pdb", tmp_path / "out.pdb")


def test_cleaning_is_deterministic(pdb_1m17, tmp_path):
    first = clean_structure(pdb_1m17, tmp_path / "a.pdb", chain="A")
    second = clean_structure(pdb_1m17, tmp_path / "b.pdb", chain="A")
    assert first.atoms_kept == second.atoms_kept
    assert (tmp_path / "a.pdb").read_text(encoding="utf-8") == \
        (tmp_path / "b.pdb").read_text(encoding="utf-8")
