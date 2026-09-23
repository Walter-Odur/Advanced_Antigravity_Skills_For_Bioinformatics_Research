"""Compound Screening — pytest test suite.

Tests ADMET filtering (two-tier), compound library queries, and docking setup.

Run:
    pytest test_compound_screening.py -v
    pytest test_compound_screening.py -v -k "not network"
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Paths — derived dynamically for portability
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent
SKILL_DIR = ROOT / ".agents" / "skills" / "compound-screening"
API_SCRIPT = SKILL_DIR / "scripts" / "compound_screening_api.py"
PYTHON = sys.executable


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run_api(*args: str, expect_rc: int = 0, timeout: int = 60):
    cmd = [PYTHON, str(API_SCRIPT)] + list(args)
    result = subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(ROOT), timeout=timeout,
    )
    assert result.returncode == expect_rc, (
        f"Expected rc={expect_rc}, got {result.returncode}\n"
        f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )
    return result


def load_json(p: Path) -> dict:
    with open(p) as f:
        return json.load(f)


@pytest.fixture(scope="session")
def output_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("screening_tests")


# ---------------------------------------------------------------------------
# Script Existence
# ---------------------------------------------------------------------------

class TestScriptExists:
    def test_api_script_exists(self):
        assert API_SCRIPT.exists()

    def test_admetlab_client_exists(self):
        assert (SKILL_DIR / "scripts" / "admetlab_client.py").exists()

    def test_compound_libraries_exists(self):
        assert (SKILL_DIR / "scripts" / "compound_libraries.py").exists()

    def test_help(self):
        r = run_api("--help")
        assert "admet-filter" in r.stdout
        assert "query-coconut" in r.stdout
        assert "query-zinc" in r.stdout
        assert "dock" in r.stdout


# ---------------------------------------------------------------------------
# ADMET Filtering
# ---------------------------------------------------------------------------

class TestADMETLipinski:
    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out = output_dir / "lipinski.json"
        run_api(
            "admet-filter",
            "--smiles", "CC(=O)OC1=CC=CC=C1C(=O)O",
            "--names", "aspirin",
            "--output", str(self.out),
        )

    def _get_aspirin(self):
        data = load_json(self.out)
        return data["compounds"][0]

    def test_file_exists(self):
        assert self.out.exists()

    def test_aspirin_passes(self):
        assert self._get_aspirin()["lipinski_pass"] is True


class TestADMETPAINS:
    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out = output_dir / "pains.json"
        run_api(
            "admet-filter",
            "--smiles", "O=C1NC(=S)SC1",
            "--names", "rhodanine",
            "--output", str(self.out),
        )

    def test_rhodanine_flagged(self):
        data = load_json(self.out)
        rhodanine = data["compounds"][0]
        assert rhodanine["pains_pass"] is False


class TestADMETTwoTier:
    """Two-tier filtering: relaxed passes more than strict."""

    COMPOUNDS = {
        "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "rhodanine": "O=C1NC(=S)SC1",
        "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    }

    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        smiles = list(self.COMPOUNDS.values())
        names = list(self.COMPOUNDS.keys())

        self.strict_out = output_dir / "strict.json"
        run_api(
            "admet-filter",
            "--smiles", *smiles, "--names", *names,
            "--strictness", "strict",
            "--output", str(self.strict_out),
        )

        self.relaxed_out = output_dir / "relaxed.json"
        run_api(
            "admet-filter",
            "--smiles", *smiles, "--names", *names,
            "--strictness", "relaxed",
            "--output", str(self.relaxed_out),
        )

    def test_relaxed_passes_more(self):
        strict = load_json(self.strict_out)
        relaxed = load_json(self.relaxed_out)
        assert relaxed["passed"] >= strict["passed"]

    def test_strict_blocks_rhodanine(self):
        data = load_json(self.strict_out)
        rhodanine = next(c for c in data["compounds"] if c["name"] == "rhodanine")
        assert rhodanine["pains_pass"] is False

    def test_relaxed_passes_rhodanine(self):
        data = load_json(self.relaxed_out)
        rhodanine = next(c for c in data["compounds"] if c["name"] == "rhodanine")
        assert rhodanine["pains_pass"] is True

    def test_strictness_in_output(self):
        assert load_json(self.strict_out)["strictness"] == "strict"
        assert load_json(self.relaxed_out)["strictness"] == "relaxed"

    def test_tier_label_in_output(self):
        assert "strict" in load_json(self.strict_out)["tier"].lower()
        assert "relaxed" in load_json(self.relaxed_out)["tier"].lower()


# ---------------------------------------------------------------------------
# Compound Libraries
# ---------------------------------------------------------------------------

class TestZINCDownloadScript:
    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out = output_dir / "zinc.sh"
        run_api("query-zinc", "--subset", "drug-like", "--output", str(self.out))

    def test_script_generated(self):
        assert self.out.exists()
        text = self.out.read_text()
        assert "wget" in text or "curl" in text


class TestZINCLeadLike:
    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out = output_dir / "zinc_lead.sh"
        run_api("query-zinc", "--subset", "lead-like", "--output", str(self.out))

    def test_lead_like(self):
        assert self.out.exists()


@pytest.mark.network
class TestCOCONUTSearch:
    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out = output_dir / "coconut.json"
        run_api(
            "query-coconut", "--query", "manzamine",
            "--limit", "5", "--output", str(self.out),
            timeout=30,
        )

    def test_search_results(self):
        data = load_json(self.out)
        assert data.get("total_results", data.get("total_found", 0)) > 0 or \
               len(data.get("compounds", data.get("results", []))) > 0
