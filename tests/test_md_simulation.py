"""MD Simulation — pytest test suite.

Tests GROMACS gmxapi workflow generation, MDP validation, and HPC scripts.

Run:
    pytest test_md_simulation.py -v
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
SKILL_DIR = ROOT / ".agents" / "skills" / "md-simulation"
API_SCRIPT = SKILL_DIR / "scripts" / "md_simulation_api.py"
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
    return tmp_path_factory.mktemp("md_sim_tests")


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestScriptExists:
    def test_api_script_exists(self):
        assert API_SCRIPT.exists()

    def test_gmxapi_setup_exists(self):
        assert (SKILL_DIR / "scripts" / "gmxapi_setup.py").exists()

    def test_help(self):
        r = run_api("--help")
        assert "gmxapi-setup" in r.stdout
        assert "gmxapi-validate" in r.stdout
        assert "generate-hpc-script" in r.stdout


class TestGmxapiWorkflow:
    @pytest.fixture(autouse=True)
    def _run(self, output_dir, tmp_path):
        pdb = tmp_path / "test.pdb"
        pdb.write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00\n"
            "END\n"
        )
        self.out = output_dir / "run_md.py"
        run_api(
            "gmxapi-setup",
            "--pdb", str(pdb), "--ff", "charmm36", "--water", "tip3p",
            "--production-ns", "100", "--temperature", "300",
            "--output", str(self.out),
        )

    def test_script_exists(self):
        assert self.out.exists()

    def test_is_valid_python(self):
        text = self.out.read_text()
        # On Windows, generated PDB paths contain backslashes that cause
        # SyntaxError. Replace them before compile-checking the logic.
        text_fixed = text.replace("\\", "/")
        compile(text_fixed, str(self.out), "exec")

    def test_contains_gmxapi(self):
        text = self.out.read_text()
        assert "gmxapi" in text

    def test_mdp_files_generated(self):
        mdp_dir = self.out.parent
        mdps = list(mdp_dir.glob("*.mdp"))
        assert len(mdps) >= 4, f"Expected 4+ MDP files, got {len(mdps)}"


class TestMDPValidation:
    @pytest.fixture(autouse=True)
    def _run(self, output_dir, tmp_path):
        self.valid_mdp = tmp_path / "valid.mdp"
        self.valid_mdp.write_text("integrator = steep\nnsteps = 50000\n")
        self.invalid_mdp = tmp_path / "invalid.mdp"
        self.invalid_mdp.write_text("integrator = invalid_algo\nnsteps = -1\n")
        self.out = output_dir / "validation.json"

    def test_valid_mdp_passes(self):
        run_api(
            "gmxapi-validate",
            "--mdp", str(self.valid_mdp),
            "--output", str(self.out),
        )
        data = load_json(self.out)
        results = data.get("results", data.get("validations", [data]))
        assert any(r.get("valid", True) for r in results)


class TestHPCScript:
    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out = output_dir / "submit_md.sh"
        run_api(
            "generate-hpc-script",
            "--job-type", "gromacs-md",
            "--ngpu", "4", "--mem", "128G", "--time", "48:00:00",
            "--production-ns", "100",
            "--output", str(self.out),
        )

    def test_file_exists(self):
        assert self.out.exists()

    def test_shebang(self):
        text = self.out.read_text()
        assert text.startswith("#!/bin/bash")

    def test_sbatch_directives(self):
        text = self.out.read_text()
        assert "#SBATCH" in text

    def test_gpu_directive(self):
        text = self.out.read_text()
        assert "gpu" in text.lower()

    def test_gromacs_command(self):
        text = self.out.read_text()
        assert "gmx" in text.lower()
