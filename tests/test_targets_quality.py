"""Structure confidence assessment, on real predicted and experimental files.

Two corrections are under test:

* pLDDT was keyed on residue number alone, so in a multi-chain file residue
  42 of chain B overwrote residue 42 of chain A.
* The "median" was ``sorted(values)[len(values) // 2]``, which is the upper
  central value for an even count.

A third behaviour is new: refusing to apply a pLDDT verdict to a
crystallographic B-factor column, where low is good and the 0-100 scale
does not apply.
"""

from __future__ import annotations

import statistics

import pytest

from agskills.errors import InvalidInputError
from agskills.targets.quality import (
    PlddtReport,
    assess_structure,
    extract_plddt,
    find_low_confidence_regions,
)


# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------


def test_extracts_plddt_from_a_real_alphafold_model(alphafold_model):
    """A real AlphaFold model must yield one value per residue, in order."""
    residues = extract_plddt(alphafold_model.read_text(encoding="utf-8"))
    assert len(residues) == 125, "AF-P9WJA5-F1 is a 125-residue protein"
    assert all(0.0 <= r.plddt <= 100.0 for r in residues)
    # Residues must come back in sequence order, which is what makes
    # "contiguous low-confidence region" meaningful.
    numbers = [r.resseq for r in residues]
    assert numbers == sorted(numbers)
    assert numbers == list(range(1, 126))
    assert {r.chain for r in residues} == {"A"}


def test_multichain_residues_are_not_conflated(pdb_6hez):
    """Chains A and B both have residue 100, and both must be reported.

    Keying on residue number alone silently discarded one of them, so a
    dimer reported the confidence of roughly half its residues.
    """
    residues = extract_plddt(pdb_6hez.read_text(encoding="utf-8"))
    chains = {r.chain for r in residues}
    assert len(chains) >= 2, f"expected several chains, got {chains}"

    # The same residue number in two chains must appear twice.
    by_number: dict[int, set[str]] = {}
    for residue in residues:
        by_number.setdefault(residue.resseq, set()).add(residue.chain)
    shared = [n for n, c in by_number.items() if len(c) >= 2]
    assert shared, "no residue number is shared between chains"

    # Every (chain, resseq) pair is unique: nothing was overwritten.
    keys = [(r.chain, r.resseq, r.icode) for r in residues]
    assert len(keys) == len(set(keys))


def test_only_the_first_model_is_read():
    """An NMR ensemble must not have its models averaged together."""
    text = (
        "MODEL        1\n"
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C\n"
        "ENDMDL\n"
        "MODEL        2\n"
        "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00 10.00           C\n"
        "ENDMDL\n"
    )
    residues = extract_plddt(text)
    assert len(residues) == 1
    assert residues[0].plddt == 90.0


def test_malformed_atom_lines_are_skipped_not_fatal():
    text = (
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C\n"
        "ATOM      2  CA  ALA A   X       0.000   0.000   0.000  1.00 ABCDE           C\n"
        "ATOM      3  CA  ALA A   3       0.000   0.000   0.000  1.00 80.00           C\n"
    )
    residues = extract_plddt(text)
    assert [r.resseq for r in residues] == [1, 3]


def test_only_ca_atoms_contribute():
    """pLDDT is a per-residue quantity, taken from the CA atom."""
    text = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N\n"
        "ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C\n"
        "ATOM      3  CB  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n"
    )
    residues = extract_plddt(text)
    assert len(residues) == 1
    assert residues[0].plddt == 90.0


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------


def test_median_is_a_true_median(tmp_path):
    """With an even count the median is the mean of the two central values.

    ``sorted(values)[len // 2]`` returns the upper one, which is wrong.
    """
    lines = []
    for index, value in enumerate([10.0, 20.0, 30.0, 40.0], start=1):
        lines.append(
            f"ATOM  {index:5d}  CA  ALA A{index:4d}       0.000   0.000   "
            f"0.000  1.00{value:6.2f}           C"
        )
    pdb = tmp_path / "even.pdb"
    pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")

    report = assess_structure(pdb)
    assert report.median_plddt == pytest.approx(25.0)
    assert report.median_plddt == pytest.approx(
        statistics.median([10.0, 20.0, 30.0, 40.0]))


def test_statistics_match_a_manual_computation(alphafold_model):
    report = assess_structure(alphafold_model)
    residues = extract_plddt(alphafold_model.read_text(encoding="utf-8"))
    values = [r.plddt for r in residues]
    assert report.n_residues == len(values)
    assert report.mean_plddt == pytest.approx(statistics.fmean(values), abs=0.01)
    assert report.min_plddt == pytest.approx(min(values))
    assert report.max_plddt == pytest.approx(max(values))


def test_bands_partition_the_residues(alphafold_model):
    """The four confidence bands must sum to the residue count."""
    report = assess_structure(alphafold_model)
    total = (report.residues_very_high + report.residues_confident
             + report.residues_low + report.residues_very_low)
    assert total == report.n_residues


def test_per_chain_breakdown_is_reported(pdb_6hez):
    report = assess_structure(pdb_6hez)
    assert len(report.per_chain) == len(report.chains)
    assert sum(c["n_residues"] for c in report.per_chain.values()) == \
        report.n_residues


# ---------------------------------------------------------------------------
# Low-confidence regions
# ---------------------------------------------------------------------------


def _residues(values, chain="A", start=1):
    from agskills.targets.quality import ResidueConfidence
    return [
        ResidueConfidence(chain=chain, resseq=start + i, icode="",
                          resname="ALA", plddt=v)
        for i, v in enumerate(values)
    ]


def test_contiguous_low_confidence_run_is_one_region():
    regions = find_low_confidence_regions(
        _residues([90, 90, 30, 30, 30, 90, 90]))
    assert len(regions) == 1
    assert regions[0]["start_residue"] == 3
    assert regions[0]["end_residue"] == 5
    assert regions[0]["length"] == 3
    assert regions[0]["mean_plddt"] == pytest.approx(30.0)


def test_two_separated_runs_are_two_regions():
    regions = find_low_confidence_regions(
        _residues([30, 30, 90, 90, 30, 30, 30]))
    assert [r["length"] for r in regions] == [2, 3]


def test_a_numbering_gap_breaks_a_region():
    """A region must be contiguous in sequence, not merely consecutive rows.

    Crystal structures have unmodelled stretches, so residue 50 can be
    followed by residue 80. Treating those as adjacent would report a
    region that does not exist.
    """
    residues = _residues([30, 30], start=10) + _residues([30, 30], start=80)
    regions = find_low_confidence_regions(residues)
    assert len(regions) == 2
    assert regions[0]["start_residue"] == 10
    assert regions[1]["start_residue"] == 80


def test_a_chain_change_breaks_a_region():
    residues = _residues([30, 30], chain="A") + _residues([30, 30], chain="B")
    regions = find_low_confidence_regions(residues)
    assert len(regions) == 2
    assert {r["chain"] for r in regions} == {"A", "B"}


def test_a_region_at_the_end_is_closed():
    regions = find_low_confidence_regions(_residues([90, 90, 30, 30]))
    assert len(regions) == 1
    assert regions[0]["end_residue"] == 4


def test_no_low_confidence_residues_gives_no_regions():
    assert find_low_confidence_regions(_residues([90, 95, 88])) == []


# ---------------------------------------------------------------------------
# Verdicts
# ---------------------------------------------------------------------------


def test_real_alphafold_model_gets_a_plddt_verdict(alphafold_model):
    report = assess_structure(alphafold_model)
    assert report.verdict in ("PROCEED", "REFINE", "STOP")
    assert "pLDDT" in report.b_factor_caveat
    assert report.reason
    assert "AlphaFold DB convention" in report.interpretation


def test_crystal_structure_is_not_given_a_plddt_verdict(pdb_1m17):
    """A B-factor column is not pLDDT, and must not be judged as one.

    For PDB 1M17 the mean B-factor is around 30, which a pLDDT scale would
    call "very low - not recommended for docking". That verdict would be
    meaningless: 1M17 is a 2.6 A crystal structure and perfectly usable.
    Lower B-factors are better, the scale is unbounded, and the whole
    pLDDT banding does not apply.
    """
    report = assess_structure(pdb_1m17)
    assert report.verdict == "NOT_APPLICABLE"
    assert "experimental structure" in report.b_factor_caveat
    assert "LOWER is better" in report.b_factor_caveat
    assert any("resolution" in action for action in report.actions)


def test_high_confidence_model_proceeds(tmp_path):
    lines = [
        f"ATOM  {i:5d}  CA  ALA A{i:4d}       0.000   0.000   0.000  1.00 "
        f"{95.0:6.2f}           C"
        for i in range(1, 101)
    ]
    pdb = tmp_path / "good.pdb"
    pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")
    report = assess_structure(pdb)
    assert report.verdict == "PROCEED"
    assert report.residues_very_high == 100


def test_very_low_confidence_model_is_stopped(tmp_path):
    lines = [
        f"ATOM  {i:5d}  CA  ALA A{i:4d}       0.000   0.000   0.000  1.00 "
        f"{25.0:6.2f}           C"
        for i in range(1, 101)
    ]
    pdb = tmp_path / "bad.pdb"
    pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")
    report = assess_structure(pdb)
    assert report.verdict == "STOP"
    assert any("AlphaFold2" in a or "colabfold" in a for a in report.actions)
    assert any("ligand-based" in a for a in report.actions)


def test_mixed_confidence_model_is_refined(tmp_path):
    """A model with a long disordered loop needs work, not rejection."""
    values = [90.0] * 60 + [30.0] * 15 + [90.0] * 25
    lines = [
        f"ATOM  {i:5d}  CA  ALA A{i:4d}       0.000   0.000   0.000  1.00 "
        f"{v:6.2f}           C"
        for i, v in enumerate(values, start=1)
    ]
    pdb = tmp_path / "mixed.pdb"
    pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")
    report = assess_structure(pdb)
    assert report.verdict == "REFINE"
    assert report.low_confidence_regions
    assert any("disordered" in a.lower() or "truncate" in a.lower()
               for a in report.actions)


def test_a_file_without_ca_atoms_raises_with_guidance(tmp_path):
    pdb = tmp_path / "ligand_only.pdb"
    pdb.write_text(
        "HETATM    1  C1  LIG A 501       0.000   0.000   0.000  1.00 20.00  "
        "         C\n", encoding="utf-8")
    with pytest.raises(InvalidInputError) as excinfo:
        assess_structure(pdb)
    assert "no CA atoms" in str(excinfo.value)


def test_report_serialises_and_rounds(alphafold_model):
    import json
    report = assess_structure(alphafold_model)
    assert isinstance(report, PlddtReport)
    data = report.as_dict()
    json.dumps(data)
    assert isinstance(data["mean_plddt"], float)
    assert round(data["mean_plddt"], 2) == data["mean_plddt"]


# ---------------------------------------------------------------------------
# B-factor scale detection
# ---------------------------------------------------------------------------


def _write_pdb(tmp_path, values, name="model.pdb"):
    lines = [
        f"ATOM  {i:5d}  CA  ALA A{i:4d}       0.000   0.000   0.000  1.00"
        f"{v:6.2f}           C"
        for i, v in enumerate(values, start=1)
    ]
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_esmfold_fractional_plddt_is_detected_and_rescaled(tmp_path):
    """The ESMFold API writes pLDDT as a fraction, not a percentage.

    Verified against the live endpoint: a 120-residue fold returned
    B-factors from 0.40 to 0.94, mean 0.84. That is an 84 pLDDT model. Read
    as a 0-100 score it becomes 0.84, which the banding calls "very low
    confidence, not recommended for docking" - the opposite of the truth,
    and how the code this replaces behaved.
    """
    from agskills.targets.quality import detect_scale, extract_plddt

    fractional = [0.40, 0.55, 0.78, 0.84, 0.88, 0.91, 0.94, 0.86, 0.80, 0.75]
    path = _write_pdb(tmp_path, fractional)
    assert detect_scale(extract_plddt(path.read_text(encoding="utf-8"))) == \
        "plddt_0_1"

    report = assess_structure(path)
    assert report.scale == "plddt_0_1"
    # Rescaled, so the mean is on the familiar 0-100 scale.
    assert report.mean_plddt == pytest.approx(77.1, abs=0.5)
    assert report.verdict in ("PROCEED", "REFINE")
    assert report.verdict != "STOP"
    assert "rescaled by 100" in report.b_factor_caveat


def test_a_percentage_plddt_model_is_detected_as_such(alphafold_model):
    from agskills.targets.quality import detect_scale, extract_plddt
    residues = extract_plddt(alphafold_model.read_text(encoding="utf-8"))
    assert detect_scale(residues) == "plddt_0_100"
    assert assess_structure(alphafold_model).scale == "plddt_0_100"


def test_crystallographic_b_factors_are_detected_as_such(pdb_1m17):
    from agskills.targets.quality import detect_scale, extract_plddt
    residues = extract_plddt(pdb_1m17.read_text(encoding="utf-8"))
    assert detect_scale(residues) == "bfactor"
    assert assess_structure(pdb_1m17).scale == "bfactor"


def test_high_b_factors_are_never_mistaken_for_confidence(tmp_path):
    """A B-factor above 100 is unambiguous."""
    from agskills.targets.quality import detect_scale, extract_plddt
    path = _write_pdb(tmp_path, [15.0, 45.0, 120.0, 85.0, 30.0])
    assert detect_scale(extract_plddt(path.read_text(encoding="utf-8"))) == \
        "bfactor"


def test_low_complexity_protein_is_not_misread_as_fractional(tmp_path):
    """A well-ordered crystal structure has low B-factors, not sub-1 ones."""
    from agskills.targets.quality import detect_scale, extract_plddt
    path = _write_pdb(tmp_path, [8.0, 9.5, 11.2, 7.8, 10.1, 12.4])
    assert detect_scale(extract_plddt(path.read_text(encoding="utf-8"))) == \
        "bfactor"
