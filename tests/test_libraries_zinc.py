"""ZINC tranche selection and download-script generation.

The defect this module exists to fix is a silent-corruption failure. The
previous URL template, ``zinc20.docking.org/tranches/{code}/download.smi``,
does not serve compound data: it returns **HTTP 200 with a 12 KB HTML
page**. Because the generated script tested only that a file existed and
then counted its lines, every tranche appeared to succeed and the pipeline
continued with web pages saved as ``.smi`` files.

The tranche letter mapping was also wrong in a second way. Its logP table
placed ``D`` at 0-1 and ``E`` at 1-2, whereas sampling the live tranches
shows ``D`` holds 1-2 and ``E`` holds 2-2.5: the table was shifted by one
bin, so even with a working URL it would have fetched the wrong chemistry.
Both mappings here were derived empirically by downloading tranches across
the grid and characterising their contents with RDKit.
"""

from __future__ import annotations

import subprocess

import pytest

from agskills.errors import UsageError
from agskills.io_utils import write_text
from agskills.libraries.zinc import (
    LOGP_BINS,
    PURCHASABILITY_CODES,
    REACTIVITY_CODES,
    SIZE_BINS,
    SUBSETS,
    ZINC_FILES_BASE,
    generate_zinc_script,
    select_tranches,
    tranche_url,
)


# ---------------------------------------------------------------------------
# The URL, which is the whole point
# ---------------------------------------------------------------------------


def test_tranche_url_uses_the_verified_bulk_endpoint():
    """The URL must be the one that actually serves SMILES.

    Confirmed against the live service: this path returns a 4.2 MB file
    beginning with the header ``smiles zinc_id``.
    """
    assert tranche_url("CDAA") == \
        "https://files.docking.org/2D/CD/CDAA.smi"
    assert ZINC_FILES_BASE == "https://files.docking.org"


def test_tranche_url_does_not_use_the_broken_template():
    """The old path served an HTML error page with HTTP 200."""
    url = tranche_url("CDAA")
    assert "zinc20.docking.org" not in url
    assert "/tranches/" not in url
    assert "download.smi" not in url


def test_the_subdirectory_is_the_first_two_letters():
    assert tranche_url("ABCD") == f"{ZINC_FILES_BASE}/2D/AB/ABCD.smi"
    assert tranche_url("JJAA") == f"{ZINC_FILES_BASE}/2D/JJ/JJAA.smi"


def test_tranche_url_honours_the_format():
    assert tranche_url("CDAA", "sdf").endswith("/CDAA.sdf")


@pytest.mark.parametrize("code", ["CD", "CDAAA", "", "CD1A", "cd-a"])
def test_malformed_tranche_codes_are_refused(code):
    with pytest.raises(UsageError):
        tranche_url(code)


def test_tranche_codes_are_case_normalised():
    assert tranche_url("cdaa") == tranche_url("CDAA")


# ---------------------------------------------------------------------------
# The letter mapping, derived empirically
# ---------------------------------------------------------------------------


def test_logp_bins_are_in_ascending_contiguous_order():
    """A gap or an overlap would silently drop or duplicate chemistry."""
    previous_high = None
    for letter, low, high in LOGP_BINS:
        assert low < high, letter
        if previous_high is not None:
            assert low == previous_high, (
                f"logP bin {letter} starts at {low}, previous ended at "
                f"{previous_high}"
            )
        previous_high = high
    letters = [b[0] for b in LOGP_BINS]
    assert letters == sorted(letters)


def test_size_bins_are_in_ascending_contiguous_order():
    previous_high = None
    for letter, low, high in SIZE_BINS:
        assert low < high, letter
        if previous_high is not None:
            assert low == previous_high, letter
        previous_high = high


def test_logp_letters_match_the_measured_bins():
    """These edges were measured, not assumed.

    Sampling tranche ``E?AA`` across the logP letters with RDKit gave:
    ``D`` 1.0-2.0, ``E`` 2.0-2.5, ``G`` 3.0-3.5, ``H`` 3.5-4.0,
    ``I`` 4.0-4.5, ``J`` 4.5-5.0. The previous table had ``D`` at 0-1 and
    ``E`` at 1-2, one bin out.
    """
    bins = {letter: (low, high) for letter, low, high in LOGP_BINS}
    assert bins["D"] == (1, 2)
    assert bins["E"] == (2, 2.5)
    assert bins["G"] == (3, 3.5)
    assert bins["H"] == (3.5, 4)
    assert bins["I"] == (4, 4.5)
    assert bins["J"] == (4.5, 5)


def test_size_letters_match_the_measured_mean_weights():
    """Measured mean molecular weight per size letter must fall in its bin.

    Sampling tranche ``?EAA`` gave means of A 174, B 227, C 279, D 313,
    E 338, F 362, G 386, H 410, I 436, J 475, K 613 Da.
    """
    measured = {"A": 174.2, "B": 227.2, "C": 278.6, "D": 312.5, "E": 337.7,
                "F": 362.1, "G": 386.4, "H": 410.1, "I": 435.6, "J": 474.6,
                "K": 612.5}
    bins = {letter: (low, high) for letter, low, high in SIZE_BINS}
    for letter, mean in measured.items():
        low, high = bins[letter]
        assert low <= mean <= high, (
            f"size letter {letter}: measured mean {mean} Da is outside its "
            f"declared bin {low}-{high}"
        )


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------


def test_drug_like_subset_selects_a_sensible_grid():
    selection = select_tranches("drug-like")
    assert selection.mw_range == (250.0, 500.0)
    assert selection.logp_range == (-1.0, 5.0)
    # Every selected code is four letters and resolves to a real URL shape.
    assert selection.codes
    for code in selection.codes:
        assert len(code) == 4
        assert tranche_url(code).startswith(ZINC_FILES_BASE)
    assert len(selection.codes) == len(selection.size_letters) * \
        len(selection.logp_letters)


def test_selection_covers_the_requested_window():
    """Every bin overlapping the window must be included.

    Dropping a boundary bin loses the compounds at the edge of the user's
    range without saying so.
    """
    selection = select_tranches("drug-like", mw_range=(300.0, 350.0))
    # 300-325 is D and 325-350 is E.
    assert "D" in selection.size_letters
    assert "E" in selection.size_letters
    # And nothing far outside is included.
    assert "A" not in selection.size_letters
    assert "K" not in selection.size_letters


def test_fragment_subset_selects_small_tranches():
    selection = select_tranches("fragment-like")
    assert selection.mw_range[1] <= 250
    assert "A" in selection.size_letters


def test_narrower_windows_select_fewer_tranches():
    wide = select_tranches("all-purchasable")
    narrow = select_tranches("fragment-like")
    assert len(wide.codes) > len(narrow.codes)


def test_reactivity_and_purchasability_are_the_third_and_fourth_letters():
    selection = select_tranches("drug-like", reactivity="C",
                                purchasability="B")
    for code in selection.codes:
        assert code[2] == "C"
        assert code[3] == "B"


def test_unknown_subset_lists_the_options():
    with pytest.raises(UsageError) as excinfo:
        select_tranches("everything")
    message = str(excinfo.value)
    for name in SUBSETS:
        assert name in message


def test_unknown_filter_letters_explain_themselves():
    with pytest.raises(UsageError) as excinfo:
        select_tranches("drug-like", reactivity="Z")
    message = str(excinfo.value)
    assert "anodyne" in message

    with pytest.raises(UsageError) as excinfo:
        select_tranches("drug-like", purchasability="Z")
    assert "in-stock" in str(excinfo.value)


def test_unsupported_format_is_refused():
    with pytest.raises(UsageError):
        select_tranches("drug-like", fmt="xyz")


def test_inverted_range_is_refused():
    with pytest.raises(UsageError):
        select_tranches("drug-like", mw_range=(500.0, 250.0))


def test_a_window_outside_zinc_coverage_is_refused_with_guidance():
    with pytest.raises(UsageError) as excinfo:
        select_tranches("drug-like", logp_range=(500.0, 600.0))
    assert "Widen" in str(excinfo.value)


def test_selection_serialises_with_its_provenance():
    import json
    data = select_tranches("lead-like").as_dict()
    json.dumps(data)
    assert data["n_tranches"] == len(data["tranches"])
    assert data["tranches"][0]["url"].startswith(ZINC_FILES_BASE)
    # The filter letters are explained, not just coded.
    assert "anodyne" in data["reactivity"]
    assert "in-stock" in data["purchasability"]


def test_filter_code_tables_are_documented():
    assert "A" in REACTIVITY_CODES and "A" in PURCHASABILITY_CODES
    for description in list(REACTIVITY_CODES.values()) + \
            list(PURCHASABILITY_CODES.values()):
        assert len(description) > 5


# ---------------------------------------------------------------------------
# The generated script
# ---------------------------------------------------------------------------


def test_generated_script_is_valid_bash(tmp_path, bash):
    selection = select_tranches("fragment-like")
    path = write_text(generate_zinc_script(selection), tmp_path / "dl.sh")
    from conftest import bash_check
    result = bash_check(bash, path)
    assert result.returncode == 0, result.stderr


def test_generated_script_rejects_html_responses():
    """The script must detect the exact failure the old URL produced.

    A tranche that does not exist is served as an HTML page with HTTP 200,
    so checking the exit status alone is not enough; the content has to be
    inspected.
    """
    script = generate_zinc_script(select_tranches("fragment-like"))
    assert "doctype html" in script.lower()
    assert "not compound data" in script
    # And the partial file is discarded rather than kept.
    assert ".part" in script
    assert "rm -f" in script


def test_generated_script_resumes_rather_than_restarting():
    """Tranches are large; an interrupted download must not start over."""
    script = generate_zinc_script(select_tranches("fragment-like"))
    assert "already downloaded" in script
    assert "[skip]" in script


def test_generated_script_counts_compounds_and_reports_failures():
    script = generate_zinc_script(select_tranches("fragment-like"))
    assert "OK_COUNT" in script and "FAIL_COUNT" in script
    assert "Compounds total" in script
    # A run where nothing downloaded must exit non-zero.
    assert 'if [ "$OK_COUNT" -eq 0 ]; then' in script
    assert "exit 1" in script


def test_generated_script_subtracts_the_header_row_from_the_count():
    """ZINC files begin with ``smiles zinc_id``, which is not a compound."""
    script = generate_zinc_script(select_tranches("fragment-like"))
    assert "- 1" in script


def test_generated_script_documents_the_selection():
    selection = select_tranches("lead-like")
    script = generate_zinc_script(selection)
    assert "# Subset:        lead-like" in script
    assert "MW window" in script
    assert "logP window" in script
    assert "Reactivity" in script
    assert f"# Tranches:      {len(selection.codes)}" in script
    assert "disk space" in script


def test_generated_script_contains_every_selected_url():
    selection = select_tranches("fragment-like")
    script = generate_zinc_script(selection)
    for code, url in selection.urls:
        assert url in script, code


def test_generated_script_requires_curl():
    script = generate_zinc_script(select_tranches("fragment-like"))
    assert "command -v curl" in script


def test_generated_script_suggests_the_next_pipeline_step():
    script = generate_zinc_script(select_tranches("fragment-like"))
    assert "admet-filter" in script
    assert "relaxed" in script


def test_download_directory_is_configurable():
    script = generate_zinc_script(select_tranches("fragment-like"),
                                  out_dir="my_library")
    assert 'OUT_DIR="my_library"' in script


def test_generation_is_deterministic():
    selection = select_tranches("drug-like")
    assert generate_zinc_script(selection) == generate_zinc_script(selection)


# ---------------------------------------------------------------------------
# Live checks, opt-in
# ---------------------------------------------------------------------------


@pytest.mark.network
def test_selected_tranche_urls_are_reachable():
    """A sample of generated URLs must serve data, not HTML.

    This is the test that would have caught the original defect.
    """
    from agskills.libraries.zinc import verify_tranches
    selection = select_tranches("fragment-like")
    result = verify_tranches(selection, limit=4)
    assert result["checked"] == 4
    assert result["reachable"] == 4, result["unreachable"]


@pytest.mark.network
def test_a_real_tranche_serves_smiles_not_a_web_page():
    """Fetch the start of a real tranche and confirm it parses as SMILES."""
    import urllib.request

    from agskills.chem.smiles import parse_smiles

    url = tranche_url("CDAA")
    request = urllib.request.Request(
        url, headers={"User-Agent": "agskills-tests/2.0",
                      "Range": "bytes=0-20000"})
    with urllib.request.urlopen(request, timeout=60) as response:
        text = response.read().decode("utf-8", errors="replace")

    assert not text.lstrip().lower().startswith("<!doctype")
    lines = text.splitlines()
    assert lines[0].split()[:2] == ["smiles", "zinc_id"]
    parsed = [parse_smiles(line.split()[0]) for line in lines[1:40]
              if line.split()]
    valid = [m for m in parsed if m is not None]
    assert len(valid) >= 30, "a real tranche should be mostly valid SMILES"
