"""SMILES loading, validation and de-duplication.

The defect that drove most of this module is silent data loss:
``list(zip(args.smiles, names))`` paired compounds with names, so three
SMILES and one name produced a one-compound run that reported success.
"""

from __future__ import annotations

import pytest

from agskills.chem.smiles import (
    CompoundRecord,
    dedupe,
    load_compounds,
    load_smiles_file,
    parse_smiles,
    standardize_smiles,
)
from agskills.errors import InvalidInputError, UsageError

from panel import DRUG_SMILES, INVALID_SMILES


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name,smiles", sorted(DRUG_SMILES.items()))
def test_real_drugs_parse(name, smiles):
    assert parse_smiles(smiles) is not None, name


@pytest.mark.parametrize("smiles", INVALID_SMILES)
def test_invalid_smiles_return_none_without_raising(smiles):
    """Malformed input must be data, not an exception.

    A batch screen has to report a bad row and carry on, so parsing never
    raises and never writes to stderr.
    """
    assert parse_smiles(smiles) is None


def test_parsing_does_not_pollute_stderr(capfd):
    """RDKit's C++ parse errors must be suppressed.

    RDKit writes them straight to the process stderr from C++, which would
    interleave with a command's own status output.
    """
    for smiles in INVALID_SMILES:
        parse_smiles(smiles)
    captured = capfd.readouterr()
    assert "SMILES Parse Error" not in captured.err
    assert "SMILES Parse Error" not in captured.out


def test_standardize_canonicalises_equivalent_spellings():
    assert standardize_smiles("C1=CC=CC=C1") == "c1ccccc1"
    assert standardize_smiles("OCC") == standardize_smiles("CCO")
    assert standardize_smiles("not a molecule") is None


def test_standardize_preserves_stereochemistry():
    """Canonicalisation must not silently discard stereo information.

    Losing a stereocentre would make two different enantiomers collapse to
    one entry during de-duplication.
    """
    l_alanine = standardize_smiles("C[C@@H](N)C(=O)O")
    d_alanine = standardize_smiles("C[C@H](N)C(=O)O")
    assert l_alanine != d_alanine
    assert "@" in l_alanine


# ---------------------------------------------------------------------------
# Loading: the silent-truncation defect
# ---------------------------------------------------------------------------


def test_mismatched_names_raise_rather_than_truncating():
    """Three SMILES with one name must be an error, not a one-compound run.

    This is the regression test for the ``zip`` defect: the old code
    returned a single compound and reported success, so two thirds of the
    input vanished without a word.
    """
    with pytest.raises(UsageError) as excinfo:
        load_compounds(smiles=["CCO", "CCC", "CCCC"], names=["only_one"])
    message = str(excinfo.value)
    assert "3" in message and "1" in message
    assert "one to one" in message


def test_too_many_names_also_raises():
    with pytest.raises(UsageError):
        load_compounds(smiles=["CCO"], names=["a", "b"])


def test_matching_names_are_paired_in_order():
    records = load_compounds(smiles=["CCO", "CCC"], names=["ethanol", "propane"])
    assert [r.name for r in records] == ["ethanol", "propane"]
    assert [r.smiles for r in records] == ["CCO", "CCC"]


def test_names_default_to_one_based_numbering():
    """Generated names must be consistent across the package.

    The old code numbered from 0 in some places and 1 in others, so the
    same compound could be ``compound_0`` or ``compound_1`` depending on
    which subcommand produced the file.
    """
    records = load_compounds(smiles=["CCO", "CCC", "CCCC"])
    assert [r.name for r in records] == ["compound_1", "compound_2",
                                         "compound_3"]


def test_no_source_raises_with_guidance():
    with pytest.raises(UsageError) as excinfo:
        load_compounds()
    assert "--smiles" in str(excinfo.value)


def test_multiple_sources_are_rejected(tmp_path):
    csv_file = tmp_path / "c.csv"
    csv_file.write_text("SMILES\nCCO\n", encoding="utf-8")
    with pytest.raises(UsageError, match="exactly one source"):
        load_compounds(smiles=["CCO"], csv_path=str(csv_file))


def test_empty_smiles_list_is_rejected():
    with pytest.raises(UsageError):
        load_compounds(smiles=[])


# ---------------------------------------------------------------------------
# CSV loading
# ---------------------------------------------------------------------------


def test_csv_column_matching_is_case_insensitive(tmp_path):
    """Real CSVs spell the column every possible way."""
    csv_file = tmp_path / "compounds.csv"
    csv_file.write_text("smiles,name\nCCO,ethanol\n", encoding="utf-8")
    records = load_compounds(csv_path=str(csv_file))
    assert len(records) == 1
    assert records[0].name == "ethanol"


def test_csv_missing_column_names_the_columns_present(tmp_path):
    csv_file = tmp_path / "compounds.csv"
    csv_file.write_text("structure,label\nCCO,ethanol\n", encoding="utf-8")
    with pytest.raises(InvalidInputError) as excinfo:
        load_compounds(csv_path=str(csv_file))
    message = str(excinfo.value)
    assert "structure" in message and "label" in message
    assert "--smiles-col" in message


def test_csv_handles_a_utf8_bom(tmp_path):
    """A CSV exported from Excel begins with a byte-order mark.

    Without ``utf-8-sig`` the first column name becomes ``\\ufeffSMILES``
    and the lookup fails with a confusing message.
    """
    csv_file = tmp_path / "excel.csv"
    csv_file.write_bytes("SMILES,Name\nCCO,ethanol\n".encode("utf-8-sig"))
    records = load_compounds(csv_path=str(csv_file))
    assert len(records) == 1


def test_csv_skips_blank_structure_rows_but_keeps_the_rest(tmp_path):
    csv_file = tmp_path / "gappy.csv"
    csv_file.write_text(
        "SMILES,Name\nCCO,ethanol\n,missing\nCCC,propane\n", encoding="utf-8")
    records = load_compounds(csv_path=str(csv_file))
    assert [r.name for r in records] == ["ethanol", "propane"]


def test_csv_with_no_usable_rows_raises(tmp_path):
    csv_file = tmp_path / "empty.csv"
    csv_file.write_text("SMILES,Name\n,\n,\n", encoding="utf-8")
    with pytest.raises(InvalidInputError, match="No usable rows"):
        load_compounds(csv_path=str(csv_file))


def test_csv_records_their_provenance(tmp_path):
    """Each record must say where it came from, for error reporting."""
    csv_file = tmp_path / "compounds.csv"
    csv_file.write_text("SMILES,Name\nCCO,ethanol\n", encoding="utf-8")
    records = load_compounds(csv_path=str(csv_file))
    assert records[0].source == "compounds.csv:1"


def test_csv_without_a_name_column_generates_names(tmp_path):
    csv_file = tmp_path / "anon.csv"
    csv_file.write_text("SMILES\nCCO\nCCC\n", encoding="utf-8")
    records = load_compounds(csv_path=str(csv_file))
    assert [r.name for r in records] == ["compound_1", "compound_2"]


# ---------------------------------------------------------------------------
# SMILES file loading
# ---------------------------------------------------------------------------


def test_smi_file_reads_structure_and_optional_name(tmp_path):
    smi = tmp_path / "library.smi"
    smi.write_text("CCO ethanol\nCCC\nc1ccccc1\tbenzene\n", encoding="utf-8")
    records = load_smiles_file(str(smi))
    assert [r.smiles for r in records] == ["CCO", "CCC", "c1ccccc1"]
    assert records[0].name == "ethanol"
    assert records[2].name == "benzene"


def test_smi_file_skips_a_zinc_header_row(tmp_path):
    """ZINC tranche files begin with ``smiles zinc_id``.

    Treating that header as a structure puts an unparseable entry at the
    top of every downloaded library.
    """
    smi = tmp_path / "zinc.smi"
    smi.write_text("smiles zinc_id\nCCO ZINC000001\nCCC ZINC000002\n",
                   encoding="utf-8")
    records = load_smiles_file(str(smi))
    assert len(records) == 2
    assert all(r.smiles != "smiles" for r in records)
    assert records[0].name == "ZINC000001"


def test_smi_file_ignores_comments_and_blank_lines(tmp_path):
    smi = tmp_path / "commented.smi"
    smi.write_text("# a comment\n\nCCO ethanol\n\n# another\nCCC propane\n",
                   encoding="utf-8")
    records = load_smiles_file(str(smi))
    assert len(records) == 2


def test_empty_smi_file_raises(tmp_path):
    smi = tmp_path / "blank.smi"
    smi.write_text("# nothing but comments\n", encoding="utf-8")
    with pytest.raises(InvalidInputError, match="No SMILES"):
        load_smiles_file(str(smi))


# ---------------------------------------------------------------------------
# De-duplication
# ---------------------------------------------------------------------------


def test_dedupe_collapses_equivalent_spellings_and_explains_each_drop():
    """Every rejected input must be reported with its reason.

    Shrinking a dataset silently is the failure mode this guards against: a
    user who supplied 100 structures needs to know that 40 were duplicates
    rather than discovering a 60-compound result.
    """
    records = [
        CompoundRecord("benzene_lower", "c1ccccc1"),
        CompoundRecord("benzene_kekule", "C1=CC=CC=C1"),
        CompoundRecord("ethanol", "CCO"),
        CompoundRecord("ethanol_reversed", "OCC"),
        CompoundRecord("broken", "NOT_A_MOLECULE"),
    ]
    kept, rejected = dedupe(records)

    assert len(kept) == 2
    assert len(rejected) == 3
    assert len(kept) + len(rejected) == len(records)

    reasons = {r["name"]: r["reason"] for r in rejected}
    assert reasons["benzene_kekule"] == "duplicate of benzene_lower"
    assert reasons["ethanol_reversed"] == "duplicate of ethanol"
    assert reasons["broken"] == "invalid SMILES"


def test_dedupe_returns_canonical_structures():
    kept, _ = dedupe([CompoundRecord("b", "C1=CC=CC=C1")])
    assert kept[0].smiles == "c1ccccc1"


def test_dedupe_keeps_the_first_occurrence_name():
    kept, _ = dedupe([
        CompoundRecord("first", "CCO"),
        CompoundRecord("second", "CCO"),
    ])
    assert len(kept) == 1
    assert kept[0].name == "first"


def test_dedupe_does_not_merge_stereoisomers():
    """Enantiomers are different compounds and must both survive."""
    kept, rejected = dedupe([
        CompoundRecord("l_ala", "C[C@@H](N)C(=O)O"),
        CompoundRecord("d_ala", "C[C@H](N)C(=O)O"),
    ])
    assert len(kept) == 2, rejected


def test_dedupe_of_empty_input_is_empty():
    kept, rejected = dedupe([])
    assert kept == [] and rejected == []


# ---------------------------------------------------------------------------
# CompoundRecord
# ---------------------------------------------------------------------------


def test_record_rejects_an_empty_structure():
    with pytest.raises(InvalidInputError, match="empty SMILES"):
        CompoundRecord(name="x", smiles="")
    with pytest.raises(InvalidInputError):
        CompoundRecord(name="x", smiles="   ")


def test_record_is_immutable():
    """Records are passed through a pipeline; accidental mutation would be
    hard to trace."""
    record = CompoundRecord(name="x", smiles="CCO")
    with pytest.raises(Exception):
        record.smiles = "CCC"  # type: ignore[misc]
