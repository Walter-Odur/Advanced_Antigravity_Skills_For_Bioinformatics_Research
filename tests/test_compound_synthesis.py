"""Compound Synthesis Skill — pytest test suite.

Tests all 8 subcommands of compound_synthesis.py with minimal parameters
for fast execution (<30 seconds total). No PyTorch or REINVENT4 install
required — tests config generation, seed prep, analysis, and HPC scripts.

Run:
    pytest test_compound_synthesis.py -v
    pytest test_compound_synthesis.py -v -k "not network"  # skip network tests
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Paths — derived dynamically for portability
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent
SKILL_DIR = ROOT / ".agents" / "skills" / "compound-synthesis"
SCRIPT = SKILL_DIR / "scripts" / "compound_synthesis.py"
REINVENT_DIR = ROOT / "REINVENT4"
PYTHON = sys.executable
OUT_DIR = ROOT / "compound_synthesis_tests"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run_cli(*args: str, expect_rc: int = 0) -> subprocess.CompletedProcess:
    """Run compound_synthesis.py with given args and assert return code."""
    cmd = [PYTHON, str(SCRIPT)] + list(args)
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(ROOT),
        timeout=60,
    )
    assert result.returncode == expect_rc, (
        f"Expected rc={expect_rc}, got {result.returncode}\n"
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )
    return result


def load_json(path: str | Path) -> dict:
    """Load and return JSON from file."""
    with open(path) as f:
        return json.load(f)


def load_toml(path: str | Path) -> dict:
    """Load and return TOML from file."""
    import tomllib
    with open(path, "rb") as f:
        return tomllib.load(f)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session", autouse=True)
def output_dir():
    """Create the test output directory."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return OUT_DIR


@pytest.fixture(scope="session")
def seeds_file(output_dir):
    """Create a validated seeds file used by multiple tests."""
    output = output_dir / "seeds.smi"
    report = output_dir / "seeds_report.json"
    run_cli(
        "prepare-seeds",
        "--smiles", "CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O",
        "--report", str(report),
        "--output", str(output),
    )
    return output


@pytest.fixture(scope="session")
def mock_csv(output_dir):
    """Create a mock REINVENT results CSV."""
    csv_path = output_dir / "mock_results.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["step", "SMILES", "total_score", "Qed", "SlogP"])
        w.writerow([1, "CCO", 0.45, 0.32, 0.80])
        w.writerow([1, "c1ccccc1", 0.62, 0.44, 0.70])
        w.writerow([2, "CC(=O)O", 0.71, 0.55, 0.90])
        w.writerow([2, "c1ccc(O)cc1", 0.83, 0.68, 0.85])
        w.writerow([3, "CC(C)CC", 0.39, 0.28, 0.60])
    return csv_path


# ===========================================================================
# T01 — check-setup
# ===========================================================================

class TestCheckSetup:
    """Test the check-setup subcommand."""

    def test_runs_and_produces_json(self, output_dir):
        out = output_dir / "t01_setup.json"
        result = run_cli(
            "check-setup",
            "--reinvent-dir", str(REINVENT_DIR),
            "--output", str(out),
        )
        assert out.exists()
        data = load_json(out)
        assert "python_version" in data
        assert "torch_available" in data
        assert "rdkit_available" in data
        assert "gpu_available" in data
        assert "priors_found" in data
        assert "priors_missing" in data
        assert "ready" in data

    def test_reports_rdkit(self, output_dir):
        data = load_json(output_dir / "t01_setup.json")
        assert data["rdkit_available"] is True
        assert data["rdkit_version"] != ""

    def test_reports_python_version(self, output_dir):
        data = load_json(output_dir / "t01_setup.json")
        assert data["python_version"].startswith("3.")


# ===========================================================================
# T02 — generate-config: staged_learning
# ===========================================================================

class TestGenerateConfigRL:
    """Test staged_learning config generation."""

    @pytest.fixture(autouse=True)
    def _generate(self, output_dir):
        self.out = output_dir / "t02_rl.toml"
        run_cli(
            "generate-config",
            "--mode", "staged_learning",
            "--generator", "reinvent",
            "--prior", ".reinvent",
            "--scoring-profile", "drug-like",
            "--num-steps", "5",
            "--batch-size", "16",
            "--device", "cpu",
            "--output", str(self.out),
        )

    def test_file_exists(self):
        assert self.out.exists()

    def test_valid_toml(self):
        data = load_toml(self.out)
        assert data["run_type"] == "staged_learning"

    def test_device_cpu(self):
        data = load_toml(self.out)
        assert data["device"] == "cpu"

    def test_batch_size(self):
        data = load_toml(self.out)
        assert data["parameters"]["batch_size"] == 16

    def test_max_steps(self):
        data = load_toml(self.out)
        assert data["stage"][0]["max_steps"] == 5

    def test_scoring_components_count(self):
        data = load_toml(self.out)
        components = data["stage"][0]["scoring"]["component"]
        assert len(components) == 6  # drug-like profile has 6 components

    def test_paths_forward_slashed(self):
        data = load_toml(self.out)
        prior = data["parameters"]["prior_file"]
        assert "\\" not in prior, f"Backslash found in path: {prior}"

    def test_aggregation_geometric_mean(self):
        data = load_toml(self.out)
        assert data["stage"][0]["scoring"]["type"] == "geometric_mean"

    def test_learning_strategy_dap(self):
        data = load_toml(self.out)
        assert data["learning_strategy"]["type"] == "dap"
        assert data["learning_strategy"]["sigma"] == 128


# ===========================================================================
# T03 — generate-config: transfer_learning
# ===========================================================================

class TestGenerateConfigTL:
    """Test transfer_learning config generation."""

    @pytest.fixture(autouse=True)
    def _generate(self, output_dir, seeds_file):
        self.out = output_dir / "t03_tl.toml"
        run_cli(
            "generate-config",
            "--mode", "transfer_learning",
            "--prior", ".reinvent",
            "--smiles-file", str(seeds_file),
            "--num-epochs", "3",
            "--batch-size", "16",
            "--device", "cpu",
            "--output", str(self.out),
        )

    def test_valid_toml(self):
        data = load_toml(self.out)
        assert data["run_type"] == "transfer_learning"

    def test_epochs(self):
        data = load_toml(self.out)
        assert data["parameters"]["num_epochs"] == 3

    def test_batch_size(self):
        data = load_toml(self.out)
        assert data["parameters"]["batch_size"] == 16

    def test_smiles_file_set(self):
        data = load_toml(self.out)
        assert "smiles_file" in data["parameters"]
        assert data["parameters"]["smiles_file"] != ""

    def test_paths_forward_slashed(self):
        data = load_toml(self.out)
        for key in ("input_model_file", "smiles_file", "output_model_file"):
            val = data["parameters"].get(key, "")
            assert "\\" not in val, f"Backslash in {key}: {val}"


# ===========================================================================
# T04 — generate-config: sampling
# ===========================================================================

class TestGenerateConfigSampling:
    """Test sampling config generation."""

    @pytest.fixture(autouse=True)
    def _generate(self, output_dir):
        self.out = output_dir / "t04_sampling.toml"
        run_cli(
            "generate-config",
            "--mode", "sampling",
            "--prior", ".reinvent",
            "--num-smiles", "50",
            "--device", "cpu",
            "--output", str(self.out),
        )

    def test_valid_toml(self):
        data = load_toml(self.out)
        assert data["run_type"] == "sampling"

    def test_num_smiles(self):
        data = load_toml(self.out)
        assert data["parameters"]["num_smiles"] == 50


# ===========================================================================
# T05 — prepare-seeds
# ===========================================================================

class TestPrepareSeeds:
    """Test seed SMILES validation and deduplication."""

    def test_valid_count(self, output_dir):
        out = output_dir / "t05_seeds.smi"
        report_path = output_dir / "t05_report.json"
        run_cli(
            "prepare-seeds",
            "--smiles", "CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O",
            "INVALID_SMILES", "CCO",
            "--report", str(report_path),
            "--output", str(out),
        )
        report = load_json(report_path)
        assert report["valid"] == 3
        assert report["invalid"] == 1
        assert report["duplicates_removed"] == 1

    def test_output_file_line_count(self, output_dir):
        out = output_dir / "t05_seeds.smi"
        if not out.exists():
            run_cli(
                "prepare-seeds",
                "--smiles", "CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O",
                "INVALID_SMILES", "CCO",
                "--output", str(out),
            )
        lines = [l.strip() for l in out.read_text().splitlines() if l.strip()]
        assert len(lines) == 3

    def test_invalid_smiles_reported(self, output_dir):
        report_path = output_dir / "t05_report.json"
        if report_path.exists():
            report = load_json(report_path)
            assert "INVALID_SMILES" in report["invalid_smiles"]

    def test_canonical_smiles(self, output_dir):
        out = output_dir / "t05_seeds.smi"
        if out.exists():
            lines = [l.strip() for l in out.read_text().splitlines() if l.strip()]
            # Benzene should be canonical
            assert "c1ccccc1" in lines

    def test_empty_input_fails(self, output_dir):
        """No SMILES provided should exit with error."""
        out = output_dir / "t05_empty.smi"
        run_cli(
            "prepare-seeds",
            "--output", str(out),
            expect_rc=1,
        )


# ===========================================================================
# T06 — generate-config: scoring
# ===========================================================================

class TestGenerateConfigScoring:
    """Test scoring config generation."""

    @pytest.fixture(autouse=True)
    def _generate(self, output_dir, seeds_file):
        self.out = output_dir / "t06_scoring.toml"
        run_cli(
            "generate-config",
            "--mode", "scoring",
            "--smiles-file", str(seeds_file),
            "--scoring-profile", "lead-like",
            "--output", str(self.out),
        )

    def test_valid_toml(self):
        data = load_toml(self.out)
        assert data["run_type"] == "scoring"

    def test_lead_like_components(self):
        data = load_toml(self.out)
        components = data["scoring"]["component"]
        assert len(components) == 7  # lead-like has 7 components

    def test_smiles_file_set(self):
        data = load_toml(self.out)
        assert "smiles_file" in data["parameters"]


# ===========================================================================
# T07 — generate-config: custom components
# ===========================================================================

class TestCustomComponents:
    """Test custom component specification."""

    @pytest.fixture(autouse=True)
    def _generate(self, output_dir):
        self.out = output_dir / "t07_custom.toml"
        run_cli(
            "generate-config",
            "--mode", "staged_learning",
            "--prior", ".reinvent",
            "--component", "Qed:weight=1.0",
            "--component", "MolecularWeight:weight=0.5,transform=double_sigmoid,low=300,high=500",
            "--component", "SAScore:weight=0.4,transform=reverse_sigmoid,high=4,low=1,k=0.5",
            "--num-steps", "5",
            "--batch-size", "16",
            "--device", "cpu",
            "--output", str(self.out),
        )

    def test_valid_toml(self):
        data = load_toml(self.out)
        assert data["run_type"] == "staged_learning"

    def test_three_components(self):
        data = load_toml(self.out)
        components = data["stage"][0]["scoring"]["component"]
        assert len(components) == 3

    def test_mw_transform(self):
        """MolecularWeight should have double_sigmoid transform."""
        content = self.out.read_text()
        assert "double_sigmoid" in content
        assert "300" in content


# ===========================================================================
# T08 — generate-config: anti-TB with seeds
# ===========================================================================

class TestAntiTBWithSeeds:
    """Test anti-TB profile with inline seeds."""

    SEED = "Cc1ccnc(NC(=O)c2ccsc2NC(=O)c2ccn3cnnc3c2)n1"

    @pytest.fixture(autouse=True)
    def _generate(self, output_dir):
        self.out = output_dir / "t08_antitb.toml"
        run_cli(
            "generate-config",
            "--mode", "staged_learning",
            "--prior", ".reinvent",
            "--scoring-profile", "anti-tb",
            "--seeds", self.SEED,
            "--num-steps", "5",
            "--batch-size", "16",
            "--device", "cpu",
            "--output", str(self.out),
        )

    def test_valid_toml(self):
        data = load_toml(self.out)
        assert data["run_type"] == "staged_learning"

    def test_anti_tb_components(self):
        data = load_toml(self.out)
        components = data["stage"][0]["scoring"]["component"]
        assert len(components) == 5  # anti-tb profile

    def test_seeds_file_created(self, output_dir):
        seeds_path = output_dir / "t08_antitb_seeds.smi"
        assert seeds_path.exists()
        content = seeds_path.read_text().strip()
        assert self.SEED in content

    def test_smiles_file_referenced(self):
        data = load_toml(self.out)
        assert "smiles_file" in data["parameters"]


# ===========================================================================
# T09 — generate-hpc-script
# ===========================================================================

class TestHPCScript:
    """Test SLURM HPC script generation."""

    @pytest.fixture(autouse=True)
    def _generate(self, output_dir):
        # Need a config file first
        config = output_dir / "t02_rl.toml"
        if not config.exists():
            run_cli(
                "generate-config",
                "--mode", "staged_learning",
                "--prior", ".reinvent",
                "--scoring-profile", "drug-like",
                "--num-steps", "5",
                "--batch-size", "16",
                "--device", "cpu",
                "--output", str(config),
            )
        self.out = output_dir / "t09_hpc.sh"
        run_cli(
            "generate-hpc-script",
            "--config", str(config),
            "--gpu", "a100",
            "--ngpu", "1",
            "--mem", "64G",
            "--time", "24:00:00",
            "--conda-env", "reinvent4",
            "--job-name", "hydra_test",
            "--output", str(self.out),
        )

    def test_file_exists(self):
        assert self.out.exists()

    def test_shebang(self):
        content = self.out.read_text()
        assert content.startswith("#!/bin/bash")

    def test_gpu_directive(self):
        content = self.out.read_text()
        assert "#SBATCH --gres=gpu:a100:1" in content

    def test_memory_directive(self):
        content = self.out.read_text()
        assert "#SBATCH --mem=64G" in content

    def test_time_directive(self):
        content = self.out.read_text()
        assert "#SBATCH --time=24:00:00" in content

    def test_conda_activate(self):
        content = self.out.read_text()
        assert "conda activate reinvent4" in content

    def test_job_name(self):
        content = self.out.read_text()
        assert "hydra_test" in content

    def test_reinvent_command(self):
        content = self.out.read_text()
        assert "reinvent" in content


# ===========================================================================
# T10 — smart device detection
# ===========================================================================

class TestSmartDeviceDetection:
    """Test auto device detection when --device is omitted."""

    @pytest.fixture(autouse=True)
    def _generate(self, output_dir):
        self.out = output_dir / "t10_smart.toml"
        self.result = run_cli(
            "generate-config",
            "--mode", "staged_learning",
            "--prior", ".reinvent",
            "--scoring-profile", "drug-like",
            "--num-steps", "300",
            "--batch-size", "128",
            "--output", str(self.out),
        )

    def test_valid_toml(self):
        data = load_toml(self.out)
        assert data["run_type"] == "staged_learning"

    def test_device_auto_detected(self):
        assert "Auto-selected device:" in self.result.stdout

    def test_device_is_valid(self):
        data = load_toml(self.out)
        assert data["device"] in ("cpu", "cuda:0")

    def test_hpc_recommended_on_cpu(self):
        """If no GPU, HPC should be recommended for 300 steps."""
        data = load_toml(self.out)
        if data["device"] == "cpu":
            assert "[HPC RECOMMENDED]" in self.result.stdout

    def test_hpc_script_auto_generated_on_cpu(self, output_dir):
        """If HPC recommended, script should be auto-generated."""
        data = load_toml(self.out)
        if data["device"] == "cpu":
            hpc = output_dir / "t10_smart_hpc.sh"
            assert hpc.exists()


# ===========================================================================
# T11 — analyze-results
# ===========================================================================

class TestAnalyzeResults:
    """Test CSV result analysis."""

    @pytest.fixture(autouse=True)
    def _analyze(self, output_dir, mock_csv):
        self.out = output_dir / "t11_analysis.json"
        run_cli(
            "analyze-results",
            "--csv", str(mock_csv),
            "--top-n", "3",
            "--sort-by", "total_score",
            "--output", str(self.out),
        )

    def test_file_exists(self):
        assert self.out.exists()

    def test_total_compounds(self):
        data = load_json(self.out)
        assert data["stats"]["total_compounds"] == 5

    def test_unique_molecules(self):
        data = load_json(self.out)
        assert data["stats"]["unique_molecules"] == 5

    def test_max_score(self):
        data = load_json(self.out)
        assert abs(data["stats"]["max_score"] - 0.83) < 0.01

    def test_mean_score(self):
        data = load_json(self.out)
        assert abs(data["stats"]["mean_score"] - 0.60) < 0.01

    def test_top_n_count(self):
        data = load_json(self.out)
        assert len(data["top_compounds"]) == 3

    def test_top_compound_is_phenol(self):
        data = load_json(self.out)
        top = data["top_compounds"][0]
        assert top["SMILES"] == "c1ccc(O)cc1"

    def test_sorted_descending(self):
        data = load_json(self.out)
        scores = [float(c["total_score"]) for c in data["top_compounds"]]
        assert scores == sorted(scores, reverse=True)


# ===========================================================================
# T12 — download-priors (network)
# ===========================================================================

@pytest.mark.network
class TestDownloadPriors:
    """Test prior model download (requires network)."""

    def test_produces_json(self, output_dir):
        out = output_dir / "t12_download.json"
        run_cli(
            "download-priors",
            "--reinvent-dir", str(REINVENT_DIR),
            "--models", "reinvent",
            "--output", str(out),
        )
        assert out.exists()
        data = load_json(out)
        assert "priors_dir" in data
        assert "downloaded" in data
        assert "skipped" in data
        assert "errors" in data


# ===========================================================================
# T13 — Scoring profiles completeness
# ===========================================================================

class TestScoringProfiles:
    """Test that all scoring profiles generate valid configs."""

    @pytest.mark.parametrize("profile,expected_count", [
        ("drug-like", 6),
        ("lead-like", 7),
        ("fragment-like", 5),
        ("kinase-inhibitor", 6),
        ("anti-tb", 5),
    ])
    def test_profile_generates_valid_toml(self, output_dir, profile, expected_count):
        out = output_dir / f"t13_{profile}.toml"
        run_cli(
            "generate-config",
            "--mode", "staged_learning",
            "--prior", ".reinvent",
            "--scoring-profile", profile,
            "--num-steps", "5",
            "--batch-size", "16",
            "--device", "cpu",
            "--output", str(out),
        )
        data = load_toml(out)
        assert data["run_type"] == "staged_learning"
        components = data["stage"][0]["scoring"]["component"]
        assert len(components) == expected_count, (
            f"Profile '{profile}' expected {expected_count} components, got {len(components)}"
        )


# ===========================================================================
# T14 — Error handling
# ===========================================================================

class TestErrorHandling:
    """Test graceful error handling for bad inputs."""

    def test_missing_smiles_file_for_tl(self, output_dir):
        """transfer_learning without --smiles-file should fail."""
        out = output_dir / "t14_tl_error.toml"
        run_cli(
            "generate-config",
            "--mode", "transfer_learning",
            "--prior", ".reinvent",
            "--device", "cpu",
            "--output", str(out),
            expect_rc=1,
        )

    def test_missing_smiles_file_for_scoring(self, output_dir):
        """scoring without --smiles-file should fail."""
        out = output_dir / "t14_scoring_error.toml"
        run_cli(
            "generate-config",
            "--mode", "scoring",
            "--device", "cpu",
            "--output", str(out),
            expect_rc=1,
        )

    def test_nonexistent_csv_for_analyze(self, output_dir):
        """analyze-results with missing CSV should fail."""
        out = output_dir / "t14_analyze_error.json"
        run_cli(
            "analyze-results",
            "--csv", "nonexistent_file.csv",
            "--output", str(out),
            expect_rc=1,
        )


# ===========================================================================
# T15 — Run subcommand (structure test)
# ===========================================================================

class TestRunSubcommand:
    """Test the run subcommand structure (without actual REINVENT4)."""

    def test_missing_config_fails(self, output_dir):
        """run with nonexistent config should fail."""
        out = output_dir / "t15_run_error.json"
        run_cli(
            "run",
            "--config", "nonexistent_config.toml",
            "--output", str(out),
            expect_rc=1,
        )
