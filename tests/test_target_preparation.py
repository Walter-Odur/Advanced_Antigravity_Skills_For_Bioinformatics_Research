"""Target Preparation — pytest test suite.

Tests structure prediction, receptor preparation, and quality assessment.

Run:
    pytest test_target_preparation.py -v
    pytest test_target_preparation.py -v -k "not network"
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

# test file lives at <ROOT>/tests/test_target_preparation.py
ROOT = Path(__file__).resolve().parent.parent  # up from tests/ to workspace root
SKILL_DIR = ROOT / ".agents" / "skills" / "target-preparation"
API_SCRIPT = SKILL_DIR / "scripts" / "target_preparation_api.py"
PYTHON = sys.executable

# Optional: DeepMind science-skills (not required for this skill to work)
SCIENCE_SKILLS = ROOT / "science-skills" / "skills"
PDB_SCRIPTS = SCIENCE_SKILLS / "pdb_database" / "scripts"
AF_SCRIPTS = SCIENCE_SKILLS / "alphafold_database_fetch_and_analyze" / "scripts"


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
    return tmp_path_factory.mktemp("target_prep_tests")


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestScriptExists:
    def test_api_script_exists(self):
        assert API_SCRIPT.exists()

    def test_predict_structure_exists(self):
        assert (SKILL_DIR / "scripts" / "predict_structure.py").exists()

    def test_help(self):
        r = run_api("--help")
        assert "prepare-receptor" in r.stdout
        assert "predict-structure" in r.stdout
        assert "assess-structure" in r.stdout


@pytest.mark.skipif(
    not SCIENCE_SKILLS.exists(),
    reason="science-skills not installed (optional dependency)",
)
class TestDeepMindScripts:
    """These tests are optional — only run if science-skills is installed."""

    def test_pdb_scripts_exist(self):
        assert PDB_SCRIPTS.exists()
        assert any(PDB_SCRIPTS.glob("*.py"))

    def test_alphafold_scripts_exist(self):
        assert AF_SCRIPTS.exists()
        assert any(AF_SCRIPTS.glob("*.py"))


class TestStructurePrediction:
    def test_sequence_validation(self):
        from importlib.util import spec_from_file_location, module_from_spec
        spec = spec_from_file_location(
            "predict_structure",
            str(SKILL_DIR / "scripts" / "predict_structure.py"),
        )
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        valid_result = mod.validate_sequence("MKTLLILAVVAAALA")
        assert valid_result[0] is True
        invalid_result = mod.validate_sequence("INVALID123!")
        assert invalid_result[0] is False

    def test_fasta_reading(self, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">test\nMKTLLI\nLAVVAA\n")
        from importlib.util import spec_from_file_location, module_from_spec
        spec = spec_from_file_location(
            "predict_structure",
            str(SKILL_DIR / "scripts" / "predict_structure.py"),
        )
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
        header, seq = mod.read_fasta(str(fasta))
        assert seq == "MKTLLILAVVAA"


class TestAF2ScriptGeneration:
    @pytest.fixture(autouse=True)
    def _run(self, output_dir, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">test\nMKTLLILAVVAAALA\n")
        self.colabfold_out = output_dir / "colabfold.sh"
        run_api(
            "generate-af2-script",
            "--fasta", str(fasta), "--method", "colabfold",
            "--gpu", "a100", "--mem", "64G",
            "--output", str(self.colabfold_out),
        )

    def test_colabfold_script(self):
        text = self.colabfold_out.read_text()
        assert "colabfold" in text.lower()
        assert "SBATCH" in text


class TestQualityAssessment:
    @pytest.fixture(autouse=True)
    def _run(self, output_dir, tmp_path):
        pdb = tmp_path / "test.pdb"
        pdb.write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 85.00\n"
            "ATOM      2  CA  ALA A   2       3.800   0.000   0.000  1.00 45.00\n"
            "END\n"
        )
        self.out = output_dir / "quality.json"
        run_api("assess-structure", "--pdb", str(pdb), "--output", str(self.out))

    def test_assess_from_pdb(self):
        data = load_json(self.out)
        assert "mean_plddt" in data or "average_plddt" in data or "plddt" in data
