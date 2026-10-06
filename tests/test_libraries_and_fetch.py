"""Database clients, structure fetching and sequence validation.

Input validation is tested offline; the live behaviour is marked
``network`` and opted into.

The AlphaFold fix is the most consequential thing here. The previous code
built the download URL by hardcoding the model version,
``AF-{accession}-F1-model_v4.pdb``. The AlphaFold Database has since moved
to v6 and withdrew the v4 files, so **every** ``--uniprot-id`` fetch
returned HTTP 404. Verified against the live service on 2026-10-06: v4 and
v5 both 404, v6 succeeds. Resolving the URL through the prediction API
keeps working across future version bumps.
"""

from __future__ import annotations

import pytest

from agskills.errors import InvalidInputError, RemoteServiceError
from agskills.libraries.chembl import ASSAY_TYPES, query_activities, search_targets
from agskills.libraries.coconut import SOURCE_COLLECTIONS, search_coconut
from agskills.targets.esmfold import (
    ESMFOLD_MAX_RESIDUES,
    read_fasta,
    validate_sequence,
)
from agskills.targets.fetch import (
    fetch_alphafold,
    fetch_pdb,
    resolve_alphafold_url,
    validate_pdb_id,
    validate_uniprot_id,
)

# A real sequence: the N-terminal domain of E. coli maltose-binding protein.
REAL_SEQUENCE = (
    "MKIEEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDII"
    "FWAHDRFGGYAQSGLLAEITPDKAFQDKLYPFTWDAVRYNGKLIAYPIAVEALSLIYNKDL"
    "LPNPPKTWEEIPALDKELKAKGKSALMFNLQEPYFTWPLIAADGGYAFKYENGKYDIKDVG"
)


# ---------------------------------------------------------------------------
# Identifier validation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("pdb_id,expected", [
    ("6HEZ", "6HEZ"), ("1m17", "1M17"), ("  4hhb  ", "4HHB"), ("1abc", "1ABC"),
])
def test_valid_pdb_ids_are_normalised(pdb_id, expected):
    assert validate_pdb_id(pdb_id) == expected


@pytest.mark.parametrize("pdb_id", [
    "", "ABCD", "6HE", "6HEZZ", "protein.pdb", "EGFR", "6-EZ", "hemoglobin",
])
def test_invalid_pdb_ids_are_refused_with_guidance(pdb_id):
    with pytest.raises(InvalidInputError) as excinfo:
        validate_pdb_id(pdb_id)
    assert "rcsb.org" in str(excinfo.value)


@pytest.mark.parametrize("accession,expected", [
    ("P00533", "P00533"), ("p9wja5", "P9WJA5"), ("  Q9Y6K9 ", "Q9Y6K9"),
    ("A0A0B4J2D5", "A0A0B4J2D5"),
])
def test_valid_uniprot_accessions_are_normalised(accession, expected):
    assert validate_uniprot_id(accession) == expected


@pytest.mark.parametrize("accession", ["", "EGFR", "P005", "1P00533", "DprE1"])
def test_a_gene_name_is_not_an_accession(accession):
    """The commonest user error: supplying a gene symbol."""
    with pytest.raises(InvalidInputError) as excinfo:
        validate_uniprot_id(accession)
    assert "uniprot.org" in str(excinfo.value)


# ---------------------------------------------------------------------------
# ChEMBL input validation
# ---------------------------------------------------------------------------


def test_blank_target_search_is_refused():
    with pytest.raises(InvalidInputError):
        search_targets("")
    with pytest.raises(InvalidInputError):
        search_targets("EGFR", limit=0)


@pytest.mark.parametrize("target_id", ["", "EGFR", "CHEMBL", "203",
                                       "CHEMBLXYZ", "chembl-203"])
def test_malformed_chembl_target_ids_are_refused(target_id):
    with pytest.raises(InvalidInputError) as excinfo:
        query_activities(target_id)
    assert "search-chembl-target" in str(excinfo.value)


def test_a_lowercase_target_id_is_accepted_after_normalisation():
    """``chembl203`` is a reasonable thing for a user to type."""
    with pytest.raises(InvalidInputError):
        query_activities("chembl")       # still malformed
    # The valid form normalises rather than being rejected on case.
    assert validate_pdb_id("1m17") == "1M17"


def test_unknown_assay_type_lists_the_codes():
    with pytest.raises(InvalidInputError) as excinfo:
        query_activities("CHEMBL203", assay_type="Z")
    message = str(excinfo.value)
    for code in ASSAY_TYPES:
        assert code in message


def test_assay_types_are_documented():
    assert ASSAY_TYPES["B"].startswith("binding")
    assert "functional" in ASSAY_TYPES["F"]


# ---------------------------------------------------------------------------
# COCONUT input validation
# ---------------------------------------------------------------------------


def test_coconut_requires_exactly_one_search_mode():
    with pytest.raises(InvalidInputError) as excinfo:
        search_coconut()
    assert "exactly one" in str(excinfo.value)

    with pytest.raises(InvalidInputError):
        search_coconut("manzamine", smiles="CCO")


def test_coconut_rejects_nonsensical_paging():
    with pytest.raises(InvalidInputError):
        search_coconut("x", limit=0)
    with pytest.raises(InvalidInputError):
        search_coconut("x", page=0)


def test_coconut_documents_its_source_collections():
    """The African and marine collections are the reason to use COCONUT."""
    for collection in ("AfroDB", "ANPDB", "NANPDB", "SANCDB", "CMNPD"):
        assert collection in SOURCE_COLLECTIONS


# ---------------------------------------------------------------------------
# Sequence validation
# ---------------------------------------------------------------------------


def test_a_real_sequence_validates():
    info = validate_sequence(REAL_SEQUENCE)
    assert info.length == len(REAL_SEQUENCE)
    assert info.sequence == REAL_SEQUENCE


def test_whitespace_and_case_are_normalised():
    info = validate_sequence("  mkie egkl\nviwi ngdk  ")
    assert info.sequence == "MKIEEGKLVIWINGDK"


def test_a_trailing_stop_codon_marker_is_tolerated():
    """Translated nucleotide input often carries a trailing asterisk."""
    assert validate_sequence(REAL_SEQUENCE + "*").length == len(REAL_SEQUENCE)


@pytest.mark.parametrize("sequence,expected", [
    ("", "empty"),
    ("MKT", "at least"),
    ("A" * 500, "at most"),
])
def test_sequences_of_the_wrong_length_are_refused(sequence, expected):
    with pytest.raises(InvalidInputError) as excinfo:
        validate_sequence(sequence)
    assert expected in str(excinfo.value)


def test_the_length_limit_is_the_one_the_endpoint_enforces():
    """A declared constant that nothing reads is worse than none.

    The previous module set ``ESMFOLD_MAX_LENGTH = 400`` and never used it,
    while its validator accepted up to 2500 residues. Sequences in between
    passed validation and then failed at the API.
    """
    assert ESMFOLD_MAX_RESIDUES == 400
    validate_sequence("A" * ESMFOLD_MAX_RESIDUES)
    with pytest.raises(InvalidInputError):
        validate_sequence("A" * (ESMFOLD_MAX_RESIDUES + 1))


def test_an_overlong_sequence_is_pointed_at_the_cluster_route():
    with pytest.raises(InvalidInputError) as excinfo:
        validate_sequence("A" * 1000)
    assert "generate-af2-script" in str(excinfo.value)


@pytest.mark.parametrize("code,note", [
    ("X", "any residue"), ("B", "aspartate"), ("Z", "glutamate"),
    ("J", "leucine"), ("U", "selenocysteine"),
])
def test_ambiguity_codes_are_explained_not_merely_rejected(code, note):
    with pytest.raises(InvalidInputError) as excinfo:
        validate_sequence("MKIEEGKLVIWING" + code + "DKGYNG")
    message = str(excinfo.value)
    assert note in message
    assert "ambiguity" in message


@pytest.mark.parametrize("sequence,label", [
    ("ATGGCGCATTGCAGCTAGCTAGCTAGCATCG", "DNA"),
    ("AUGGCGCAUUGCAGCUAGCUAGCUAGCAUCG", "RNA"),
    ("ACGTACGTACGTNNNNACGTACGTACGT", "DNA"),
])
def test_a_nucleotide_sequence_is_caught_by_composition(sequence, label):
    """A DNA sequence cannot be caught by checking characters.

    A, C, G and T are alanine, cysteine, glycine and threonine, and N is
    asparagine, so "ACGTACGT..." is a valid - if improbable - protein
    sequence. Folding DNA as protein produces confident-looking nonsense,
    so the composition is checked instead.
    """
    with pytest.raises(InvalidInputError) as excinfo:
        validate_sequence(sequence)
    message = str(excinfo.value)
    assert label in message
    assert "Translate" in message


def test_a_real_protein_is_not_mistaken_for_a_nucleotide_sequence():
    """The composition check must not fire on genuine protein."""
    assert validate_sequence(REAL_SEQUENCE).length == len(REAL_SEQUENCE)
    # A short peptide that happens to use only ACGTN is left alone, since
    # the check needs enough residues to be confident.
    assert validate_sequence("ACGTNACGTN").length == 10


# ---------------------------------------------------------------------------
# FASTA reading
# ---------------------------------------------------------------------------


def test_reads_a_single_record_fasta(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(f">sp|P0AEX9|MALE_ECOLI Maltose-binding protein\n"
                     f"{REAL_SEQUENCE[:60]}\n{REAL_SEQUENCE[60:]}\n",
                     encoding="utf-8")
    info = read_fasta(fasta)
    assert info.sequence == REAL_SEQUENCE
    assert "MALE_ECOLI" in info.header
    assert info.warnings == []


def test_a_multi_record_fasta_reads_the_first_and_says_so(tmp_path):
    """Silently folding only the first of several sequences is confusing."""
    fasta = tmp_path / "multi.fasta"
    fasta.write_text(">first\nMKIEEGKLVIWINGDK\n>second\nAAAAAAAAAAAA\n",
                     encoding="utf-8")
    info = read_fasta(fasta)
    assert info.sequence == "MKIEEGKLVIWINGDK"
    assert any("more than one record" in w for w in info.warnings)


def test_an_empty_fasta_is_refused(tmp_path):
    fasta = tmp_path / "empty.fasta"
    fasta.write_text(">header only\n", encoding="utf-8")
    with pytest.raises(InvalidInputError) as excinfo:
        read_fasta(fasta)
    assert "no sequence" in str(excinfo.value)


def test_a_missing_fasta_raises():
    from agskills.errors import ResourceNotFoundError
    with pytest.raises(ResourceNotFoundError):
        read_fasta("definitely_absent.fasta")


# ---------------------------------------------------------------------------
# Guarding against HTML served with HTTP 200
# ---------------------------------------------------------------------------


def test_an_html_body_is_not_accepted_as_a_structure():
    """A server returning an error page with HTTP 200 must be detected.

    This is exactly how the old ZINC URL corrupted compound libraries, and
    a PDB download is equally exposed to it.
    """
    from agskills.targets.fetch import _looks_like_pdb
    assert not _looks_like_pdb("<!DOCTYPE html><html><body>404</body></html>")
    assert not _looks_like_pdb("<html>Not Found</html>")
    assert _looks_like_pdb("HEADER    OXIDOREDUCTASE\nATOM      1  N ...")
    assert _looks_like_pdb(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000")


def test_an_html_body_is_not_accepted_as_a_prediction():
    from agskills.targets.esmfold import _looks_like_pdb
    assert not _looks_like_pdb("<!DOCTYPE html>")
    assert _looks_like_pdb(
        "ATOM      1  N   MET A   1      -8.901   4.127  -0.555")


# ---------------------------------------------------------------------------
# Live behaviour, opt-in
# ---------------------------------------------------------------------------


@pytest.mark.network
def test_fetches_a_real_pdb_entry(tmp_path):
    download = fetch_pdb("1M17", tmp_path)
    assert download.identifier == "1M17"
    assert download.source == "RCSB PDB"
    assert download.path.is_file()
    text = download.path.read_text(encoding="utf-8")
    assert text.startswith("HEADER")
    assert "AQ4" in text, "1M17 should contain erlotinib"


@pytest.mark.network
def test_a_nonexistent_pdb_entry_is_a_clear_error(tmp_path):
    with pytest.raises(RemoteServiceError) as excinfo:
        fetch_pdb("9ZZZ", tmp_path)
    assert "404" in str(excinfo.value) or "not exist" in str(excinfo.value)


@pytest.mark.network
def test_alphafold_url_is_resolved_not_guessed():
    """The hardcoded v4 URL is dead; the API reports the current version."""
    entry = resolve_alphafold_url("P9WJA5")
    assert entry["pdbUrl"].startswith("https://alphafold.ebi.ac.uk/files/")
    assert int(entry["latestVersion"]) >= 4
    assert entry["entryId"] == "AF-P9WJA5-F1"


@pytest.mark.network
def test_the_hardcoded_v4_url_really_is_dead():
    """Documents why version resolution is necessary rather than tidy.

    If this ever starts passing, the AlphaFold Database has restored the
    v4 files and the note in the module docstring should be revisited.
    """
    import urllib.error
    import urllib.request
    request = urllib.request.Request(
        "https://alphafold.ebi.ac.uk/files/AF-P9WJA5-F1-model_v4.pdb",
        headers={"User-Agent": "agskills-tests/2.0"})
    with pytest.raises(urllib.error.HTTPError) as excinfo:
        urllib.request.urlopen(request, timeout=30)
    assert excinfo.value.code == 404


@pytest.mark.network
def test_fetches_a_real_alphafold_model(tmp_path):
    download = fetch_alphafold("P9WJA5", tmp_path)
    assert download.source == "AlphaFold Database"
    assert download.path.is_file()
    assert download.metadata["model_version"]
    # A predicted model carries no ligands, and the caveat must say so.
    assert "no ligands" in download.metadata["caveat"]
    from agskills.targets.quality import extract_plddt
    residues = extract_plddt(download.path.read_text(encoding="utf-8"))
    assert len(residues) == 125
    assert all(0 <= r.plddt <= 100 for r in residues)


@pytest.mark.network
def test_an_accession_without_a_model_is_reported_clearly():
    with pytest.raises(RemoteServiceError) as excinfo:
        resolve_alphafold_url("A0A000A000")
    assert "predict" in str(excinfo.value).lower()


@pytest.mark.network
def test_searches_real_chembl_targets():
    """DprE1 is a validated anti-tubercular target."""
    result = search_targets("DprE1", limit=5)
    assert result["total_returned"] >= 1
    organisms = " ".join(t["organism"] or "" for t in result["targets"])
    assert "Mycobacterium" in organisms
    for target in result["targets"]:
        assert target["target_chembl_id"].startswith("CHEMBL")


@pytest.mark.network
def test_retrieves_real_chembl_activities():
    """EGFR actives at pChEMBL 8 or better are 4-anilinoquinazolines."""
    result = query_activities("CHEMBL203", pchembl_min=8.0, limit=10)
    assert result["returned"] >= 5
    assert result["distinct_molecules_found"] >= result["returned"]

    for compound in result["compounds"]:
        assert compound["pchembl_median"] >= 8.0
        assert compound["n_measurements"] >= 1
        from agskills.chem.smiles import parse_smiles
        assert parse_smiles(compound["canonical_smiles"]) is not None

    # Ranked on the median, descending.
    medians = [c["pchembl_median"] for c in result["compounds"]]
    assert medians == sorted(medians, reverse=True)
    # And aggregated: no molecule appears twice.
    ids = [c["molecule_chembl_id"] for c in result["compounds"]]
    assert len(ids) == len(set(ids))


@pytest.mark.network
def test_chembl_aggregation_collapses_repeated_measurements():
    """ChEMBL holds many activities per molecule.

    Returning them unaggregated means a "top 20" can be one compound
    twenty times.
    """
    result = query_activities("CHEMBL203", pchembl_min=7.0, limit=50)
    assert result["activities_examined"] > result["distinct_molecules_found"]
    assert any(c["n_measurements"] > 1 for c in result["compounds"])


@pytest.mark.network
def test_searches_real_coconut_natural_products():
    result = search_coconut("antimycobacterial", limit=10)
    assert result["source"] == "COCONUT 2.0"
    assert result["returned"] >= 1
    from agskills.chem.smiles import parse_smiles
    for compound in result["compounds"]:
        if compound["smiles"]:
            assert parse_smiles(compound["smiles"]) is not None


@pytest.mark.network
@pytest.mark.slow
def test_predicts_a_real_structure_with_esmfold(tmp_path):
    """Fold a real sequence and assess the model, end to end."""
    from agskills.targets.esmfold import predict_structure
    output = tmp_path / "predicted.pdb"
    result = predict_structure(REAL_SEQUENCE[:120], output)
    assert output.is_file()
    assert result["sequence_length"] == 120
    assert result["method"].startswith("ESMFold")
    quality = result["quality"]
    assert quality["n_residues"] == 120
    assert quality["verdict"] in ("PROCEED", "REFINE", "STOP")
    assert 0 <= quality["mean_plddt"] <= 100
    # The single-sequence limitation must be stated.
    assert "no multiple sequence alignment" in result["limitation"]
