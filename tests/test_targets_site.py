"""Binding-site definition on real holo structures.

Two scientific properties dominate here.

**The box must land on the ligand of interest.** For 1M17 (EGFR with
erlotinib) there is one drug-like ligand and the answer is unambiguous. For
6HEZ (DprE1) there is both an FAD cofactor of 53 atoms and an inhibitor of
28, so ranking heterogens by size picks the cofactor and aims the box at
the flavin site rather than the inhibitor pocket.

**The whole feature is new.** The code this replaces documented a
``--site binding_site.json`` input to docking, and defaulted the SLURM
generator to that filename, but no subcommand in any skill produced such a
file; the alternative, ``--center``, raised ``AttributeError`` because the
parser defined ``--size`` while the function read ``args.box_size``.
"""

from __future__ import annotations

import pytest

from agskills.errors import InvalidInputError
from agskills.targets.site import (
    COFACTORS,
    MIN_BOX_A,
    BindingSite,
    is_cofactor,
    list_ligands,
    parse_atoms,
    site_from_hetatm,
    site_from_residues,
    site_from_structure,
)


# ---------------------------------------------------------------------------
# Ligand discovery
# ---------------------------------------------------------------------------


def test_finds_erlotinib_in_1m17(pdb_1m17):
    """1M17 has exactly one drug-like ligand: erlotinib, residue AQ4."""
    ligands = list_ligands(pdb_1m17)
    assert len(ligands) == 1
    ligand = ligands[0]
    assert ligand["resname"] == "AQ4"
    assert ligand["chain"] == "A"
    assert ligand["resseq"] == 999
    assert ligand["n_atoms"] == 29
    assert ligand["is_cofactor"] is False
    assert ligand["role"] == "drug-like ligand"


def test_solvent_and_additives_are_not_offered_as_ligands(pdb_1m17):
    """Water and buffer components never mark a druggable pocket."""
    ligands = list_ligands(pdb_1m17)
    names = {entry["resname"] for entry in ligands}
    assert "HOH" not in names
    assert not names & {"GOL", "EDO", "SO4", "CL", "NA", "DMS"}


def test_drug_like_ligand_outranks_the_cofactor_in_6hez(pdb_6hez):
    """The inhibitor must be ranked above FAD despite being smaller.

    This is the regression test for the cofactor defect. 6HEZ holds FAD
    (53 atoms) and inhibitor 0SK (28 atoms). Size-only ranking returns FAD
    and centres the box on the flavin site, so every affinity computed
    afterwards would be for the wrong cavity.
    """
    ligands = list_ligands(pdb_6hez)
    names = [entry["resname"] for entry in ligands]
    assert "FAD" in names and "0SK" in names

    first = ligands[0]
    assert first["resname"] == "0SK", (
        f"ranked {first['resname']} first; expected the inhibitor 0SK"
    )
    assert first["is_cofactor"] is False

    # FAD is still reported, correctly labelled, just not preferred.
    fad = next(entry for entry in ligands if entry["resname"] == "FAD")
    assert fad["is_cofactor"] is True
    assert "cofactor" in fad["role"]
    assert fad["n_atoms"] > first["n_atoms"], (
        "the point of the test is that the cofactor is the larger one"
    )


def test_cofactor_table_covers_the_common_prosthetic_groups():
    for name in ("FAD", "FMN", "NAD", "NAP", "HEM", "ATP", "ADP", "SAM",
                 "PLP", "COA", "NAG", "GSH"):
        assert is_cofactor(name), name
    for name in ("0SK", "AQ4", "STI", "LIG"):
        assert not is_cofactor(name), name
    assert "FAD" in COFACTORS


def test_is_cofactor_is_case_and_whitespace_insensitive():
    assert is_cofactor("fad") and is_cofactor(" FAD ")


# ---------------------------------------------------------------------------
# Site from a ligand
# ---------------------------------------------------------------------------


def test_site_from_erlotinib_is_centred_on_the_ligand(pdb_1m17):
    """The box centre must equal the ligand's own centroid."""
    site = site_from_hetatm(pdb_1m17)
    assert site.reference == "AQ4:A:999"
    assert site.n_reference_atoms == 29
    assert "high" in site.confidence

    # Compute the centroid independently from the coordinates.
    atoms = [a for a in parse_atoms(pdb_1m17.read_text(encoding="utf-8"))
             if a.is_hetatm and a.resname == "AQ4" and a.element != "H"]
    assert len(atoms) == 29
    expected = (sum(a.x for a in atoms) / len(atoms),
                sum(a.y for a in atoms) / len(atoms),
                sum(a.z for a in atoms) / len(atoms))
    assert site.center_x == pytest.approx(expected[0], abs=0.01)
    assert site.center_y == pytest.approx(expected[1], abs=0.01)
    assert site.center_z == pytest.approx(expected[2], abs=0.01)


def test_site_in_6hez_defaults_to_the_inhibitor_not_the_cofactor(pdb_6hez):
    """The automatic choice must be the inhibitor, and must say so."""
    site = site_from_hetatm(pdb_6hez)
    assert site.reference.startswith("0SK")
    assert site.selected_ligand is not None
    assert site.selected_ligand["is_cofactor"] is False
    # The alternatives, including FAD, are surfaced so the user can switch.
    labels = {entry["label"] for entry in site.alternative_ligands}
    assert any(label.startswith("FAD") for label in labels)
    assert any("--ligand-resname" in warning for warning in site.warnings)


def test_the_box_encloses_the_reference_ligand_with_padding(pdb_1m17):
    """Every ligand atom must sit inside the box, with room to spare."""
    padding = 4.0
    site = site_from_hetatm(pdb_1m17, padding=padding)
    atoms = [a for a in parse_atoms(pdb_1m17.read_text(encoding="utf-8"))
             if a.is_hetatm and a.resname == "AQ4" and a.element != "H"]
    half = (site.size_x / 2, site.size_y / 2, site.size_z / 2)
    for atom in atoms:
        assert abs(atom.x - site.center_x) <= half[0]
        assert abs(atom.y - site.center_y) <= half[1]
        assert abs(atom.z - site.center_z) <= half[2]

    span_x = max(a.x for a in atoms) - min(a.x for a in atoms)
    assert site.size_x == pytest.approx(
        max(MIN_BOX_A, span_x + 2 * padding), abs=0.01)


def test_padding_widens_the_box(pdb_1m17):
    narrow = site_from_hetatm(pdb_1m17, padding=2.0)
    wide = site_from_hetatm(pdb_1m17, padding=8.0)
    assert wide.volume > narrow.volume
    assert wide.center == pytest.approx(narrow.center, abs=1e-6)


def test_a_small_ligand_box_is_clamped_to_a_usable_minimum(tmp_path):
    """A box smaller than a drug cannot hold a pose in any orientation."""
    pdb = tmp_path / "tiny.pdb"
    pdb.write_text(
        "HETATM    1  C1  LIG A 501       0.000   0.000   0.000  1.00 20.00"
        "           C\n"
        "HETATM    2  C2  LIG A 501       1.000   0.000   0.000  1.00 20.00"
        "           C\n",
        encoding="utf-8")
    site = site_from_hetatm(pdb, padding=1.0)
    assert min(site.size) >= MIN_BOX_A
    assert any("minimum" in warning for warning in site.warnings)


def test_explicit_ligand_selection_overrides_the_ranking(pdb_6hez):
    """A user who wants the cofactor site must be able to ask for it."""
    site = site_from_hetatm(pdb_6hez, resname="FAD", chain="A")
    assert site.reference.startswith("FAD")
    inhibitor = site_from_hetatm(pdb_6hez, resname="0SK", chain="A")
    # The two sites are genuinely different pockets.
    separation = sum((a - b) ** 2 for a, b
                     in zip(site.center, inhibitor.center)) ** 0.5
    assert separation > 5.0, (
        f"FAD and 0SK centres are only {separation:.1f} A apart"
    )


def test_selecting_a_single_chain_narrows_the_box(pdb_6hez):
    """Averaging a ligand across two chains would centre the box on nothing.

    6HEZ is a dimer with 0SK bound in both chains. A box spanning both
    copies would be centred in the space between them.
    """
    chain_a = site_from_hetatm(pdb_6hez, resname="0SK", chain="A")
    chain_b = site_from_hetatm(pdb_6hez, resname="0SK", chain="B")
    assert chain_a.n_reference_atoms == 28
    assert chain_b.n_reference_atoms == 28
    separation = sum((a - b) ** 2 for a, b
                     in zip(chain_a.center, chain_b.center)) ** 0.5
    assert separation > 10.0
    # The automatic choice must pick one copy, not both.
    automatic = site_from_hetatm(pdb_6hez)
    assert automatic.n_reference_atoms == 28


def test_structure_without_a_ligand_raises_with_alternatives(tmp_path):
    pdb = tmp_path / "apo.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00"
        "           C\n"
        "HETATM    2  O   HOH A 100       5.000   0.000   0.000  1.00 20.00"
        "           O\n",
        encoding="utf-8")
    with pytest.raises(InvalidInputError) as excinfo:
        site_from_hetatm(pdb)
    message = str(excinfo.value)
    assert "no co-crystallised ligand" in message
    assert "--residues" in message


def test_unmatched_ligand_selector_lists_the_candidates(pdb_1m17):
    with pytest.raises(InvalidInputError) as excinfo:
        site_from_hetatm(pdb_1m17, resname="XYZ")
    assert "AQ4" in str(excinfo.value)


# ---------------------------------------------------------------------------
# Site from residues
# ---------------------------------------------------------------------------


def test_site_from_pocket_residues(pdb_1m17):
    """Residues lining the erlotinib site must produce a nearby box.

    Thr766, Met769 and Leu820 are the gatekeeper and hinge residues of the
    EGFR ATP pocket, so a box built from them should sit close to the box
    built from erlotinib itself.
    """
    from_residues = site_from_residues(pdb_1m17, ["A:766", "A:769", "A:820"])
    from_ligand = site_from_hetatm(pdb_1m17)
    separation = sum((a - b) ** 2 for a, b
                     in zip(from_residues.center, from_ligand.center)) ** 0.5
    assert separation < 10.0, (
        f"pocket residues gave a centre {separation:.1f} A from the ligand"
    )
    assert "medium" in from_residues.confidence


def test_bare_residue_numbers_match_any_chain(pdb_6hez):
    site = site_from_residues(pdb_6hez, ["100"])
    assert site.n_reference_atoms > 0


def test_unparseable_residue_selector_raises(pdb_1m17):
    with pytest.raises(InvalidInputError) as excinfo:
        site_from_residues(pdb_1m17, ["A:not_a_number"])
    assert "CHAIN:RESSEQ" in str(excinfo.value)


def test_no_matching_residues_raises_with_a_numbering_hint(pdb_1m17):
    with pytest.raises(InvalidInputError) as excinfo:
        site_from_residues(pdb_1m17, ["Z:99999"])
    assert "numbering" in str(excinfo.value)


def test_partial_residue_match_is_reported(pdb_1m17):
    site = site_from_residues(pdb_1m17, ["A:766", "Z:99999"])
    assert any("found" in warning for warning in site.warnings)


def test_empty_residue_list_raises(pdb_1m17):
    with pytest.raises(InvalidInputError):
        site_from_residues(pdb_1m17, [])


# ---------------------------------------------------------------------------
# Whole-structure (blind) docking
# ---------------------------------------------------------------------------


def test_blind_box_encloses_the_protein_and_is_labelled_low_confidence(
        pdb_1m17):
    site = site_from_structure(pdb_1m17)
    assert "low" in site.confidence
    assert "blind" in site.method
    assert any("sparsely" in warning for warning in site.warnings)

    ligand_site = site_from_hetatm(pdb_1m17)
    assert site.volume > ligand_site.volume * 5, (
        "a whole-protein box should be far larger than a pocket box"
    )


def test_blind_box_warns_about_exhaustiveness(pdb_1m17):
    site = site_from_structure(pdb_1m17)
    assert any("exhaustiveness" in warning for warning in site.warnings)


# ---------------------------------------------------------------------------
# Serialisation: the contract docking consumes
# ---------------------------------------------------------------------------


def test_serialised_site_uses_the_flat_keys_docking_reads(pdb_1m17):
    """``dock --site`` reads these exact keys, so they are an interface."""
    import json
    data = site_from_hetatm(pdb_1m17).as_dict()
    json.dumps(data)
    for key in ("center_x", "center_y", "center_z",
                "size_x", "size_y", "size_z"):
        assert key in data, key
        assert isinstance(data[key], float)
    assert data["volume_A3"] > 0
    assert data["method"]
    assert data["confidence"]


def test_site_round_trips_through_json(pdb_1m17, tmp_path):
    """A written site must reconstruct exactly, since docking reloads it."""
    import json
    original = site_from_hetatm(pdb_1m17)
    path = tmp_path / "site.json"
    path.write_text(json.dumps(original.as_dict()), encoding="utf-8")
    data = json.loads(path.read_text(encoding="utf-8"))
    restored = BindingSite(
        center_x=data["center_x"], center_y=data["center_y"],
        center_z=data["center_z"], size_x=data["size_x"],
        size_y=data["size_y"], size_z=data["size_z"],
        method=data["method"],
    )
    assert restored.center == pytest.approx(
        tuple(round(v, 3) for v in original.center))


def test_site_detection_is_deterministic(pdb_6hez):
    """The same structure must always give the same box."""
    first = site_from_hetatm(pdb_6hez).as_dict()
    second = site_from_hetatm(pdb_6hez).as_dict()
    assert first == second


# ---------------------------------------------------------------------------
# Atom parsing
# ---------------------------------------------------------------------------


def test_parse_atoms_reads_real_coordinates(pdb_1m17):
    atoms = parse_atoms(pdb_1m17.read_text(encoding="utf-8"))
    assert len(atoms) > 2000
    assert any(a.is_hetatm for a in atoms)
    assert any(not a.is_hetatm for a in atoms)
    first = next(a for a in atoms if not a.is_hetatm)
    assert first.chain == "A"
    assert first.resname.isalpha()


def test_parse_atoms_tolerates_truncated_lines():
    """Some PDB files omit the element column."""
    text = ("ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 "
            "20.00\n")
    atoms = parse_atoms(text)
    assert len(atoms) == 1
    assert atoms[0].x == 1.0 and atoms[0].z == 3.0
    assert atoms[0].element == ""
