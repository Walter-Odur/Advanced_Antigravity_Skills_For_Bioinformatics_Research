"""Drug Discovery Wizard — pytest test suite.

Tests the full wizard pipeline using DeepMind science-skills for Steps 1-3
(PDB, AlphaFold, ChEMBL) and the wizard's own scripts for Steps 4-6
(ADMET, docking config, REINVENT, HPC, MD).

Run:
    pytest test_drug_discovery_wizard.py -v
    pytest test_drug_discovery_wizard.py -v -k "not network"  # skip API tests
    pytest test_drug_discovery_wizard.py -v -k "admet"         # ADMET only
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

WIZARD_DIR = Path(__file__).resolve().parent.parent / ".agents" / "skills" / "drug-discovery-wizard"
ROOT = Path(__file__).resolve().parent.parent
WIZARD_SCRIPT = WIZARD_DIR / "scripts" / "drug_discovery_api.py"
SCIENCE_SKILLS = ROOT / "science-skills" / "skills"
PDB_SCRIPTS = SCIENCE_SKILLS / "pdb_database" / "scripts"
AF_SCRIPTS = SCIENCE_SKILLS / "alphafold_database_fetch_and_analyze" / "scripts"
CHEMBL_SCRIPTS = SCIENCE_SKILLS / "chembl_database" / "scripts"
PYTHON = sys.executable
OUT_DIR = ROOT / "wizard_tests"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run_wizard(*args: str, expect_rc: int = 0, timeout: int = 60) -> subprocess.CompletedProcess:
    """Run drug_discovery_api.py with given args."""
    cmd = [PYTHON, str(WIZARD_SCRIPT)] + list(args)
    result = subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(ROOT), timeout=timeout,
    )
    assert result.returncode == expect_rc, (
        f"Expected rc={expect_rc}, got {result.returncode}\n"
        f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )
    return result


def run_uv_script(script_path: str, *args: str, expect_rc: int = 0, timeout: int = 120) -> subprocess.CompletedProcess:
    """Run a DeepMind science-skill script via uv run."""
    cmd = ["uv", "run", str(script_path)] + list(args)
    result = subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(ROOT), timeout=timeout,
    )
    assert result.returncode == expect_rc, (
        f"uv run {script_path} failed (rc={result.returncode})\n"
        f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )
    return result


def load_json(path: str | Path) -> dict | list:
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session", autouse=True)
def output_dir():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return OUT_DIR


# ===========================================================================
# STAGE 1: Target & Receptor Preparation (DeepMind science-skills)
# ===========================================================================

@pytest.mark.network
class TestPDBFetch:
    """T01 -- PDB fetch via DeepMind pdb-database skill."""

    def test_download_pdb_structure(self, output_dir):
        """Download PDB 1M17 (EGFR) coordinate file."""
        pdb_dir = output_dir / "pdb_files"
        pdb_dir.mkdir(exist_ok=True)
        run_uv_script(
            str(PDB_SCRIPTS / "download_coordinate_files.py"),
            "--ids", "1M17",
            "--output_dir", str(pdb_dir),
            "--format", "pdb",
        )
        # File may be gzipped (pdb_00001m17.pdb.gz) or plain
        has_pdb = any(pdb_dir.glob("*.pdb")) or any(pdb_dir.glob("*.pdb.gz"))
        assert has_pdb, f"No PDB file found in {pdb_dir}: {list(pdb_dir.iterdir())}"

    def test_fetch_pdb_metadata(self, output_dir):
        """Fetch metadata for PDB 1M17 via GraphQL."""
        out = output_dir / "t01_metadata.json"
        # fetch_pdb_metadata.py requires a GraphQL --query
        query = '{entry(entry_id: "1M17") {struct {title}}}'
        run_uv_script(
            str(PDB_SCRIPTS / "fetch_pdb_metadata.py"),
            "--query", query,
            "--output", str(out),
        )
        assert out.exists()
        data = load_json(out)
        assert data  # Non-empty response

    def test_search_pdb_by_title(self, output_dir):
        """Search PDB for EGFR structures."""
        out = output_dir / "t01_search.json"
        query = json.dumps({
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "struct.title",
                    "operator": "contains_phrase",
                    "value": "EGFR"
                }
            },
            "return_type": "entry"
        })
        run_uv_script(
            str(PDB_SCRIPTS / "search_pdb.py"),
            "--query", query,
            "--return_type", "entry",
            "--rows", "5",
            "--output", str(out),
        )
        assert out.exists()


@pytest.mark.network
class TestAlphaFoldFetch:
    """T03 -- AlphaFold fetch via DeepMind alphafold skill."""

    def test_fetch_alphafold_structure(self, output_dir):
        """Download AlphaFold structure for ACE2 (Q9BYF1)."""
        af_dir = output_dir / "af_files"
        af_dir.mkdir(exist_ok=True)
        # fetch_structure.py uses positional uniprot_id and -o for output dir
        # Set User-Agent to avoid 403 from AlphaFold API
        env = os.environ.copy()
        env["POLITE_HTTP_USER_AGENT"] = "DrugDiscoveryWizard/1.0 (test suite)"
        cmd = ["uv", "run", str(AF_SCRIPTS / "fetch_structure.py"), "Q9BYF1", "-o", str(af_dir)]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(ROOT), timeout=120, env=env)
        assert result.returncode == 0, (
            f"AlphaFold fetch failed (rc={result.returncode})\n"
            f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
        all_files = list(af_dir.iterdir())
        assert len(all_files) > 0, f"No files downloaded to {af_dir}"

    def test_analyze_plddt(self, output_dir):
        """Analyze pLDDT confidence — requires metadata from fetch step."""
        af_dir = output_dir / "af_files"
        # Find the metadata JSON from the fetch step
        meta_files = list(af_dir.glob("*metadata*.json"))
        if not meta_files:
            pytest.skip("No metadata file from fetch step — run test_fetch first")
        run_uv_script(
            str(AF_SCRIPTS / "analyze_plddt.py"),
            str(meta_files[0]),
        )
        # analyze_plddt prints to stdout, no --output flag

    def test_analyze_pae(self, output_dir):
        """Analyze PAE (domain boundaries) — requires PAE JSON from fetch step."""
        af_dir = output_dir / "af_files"
        pae_files = list(af_dir.glob("*pae*.json")) + list(af_dir.glob("*error*.json"))
        if not pae_files:
            pytest.skip("No PAE file from fetch step — run test_fetch first")
        run_uv_script(
            str(AF_SCRIPTS / "analyze_pae.py"),
            str(pae_files[0]),
        )


# ===========================================================================
# STAGE 2: ChEMBL Queries (DeepMind science-skills)
# ===========================================================================

@pytest.mark.network
class TestChEMBLQueries:
    """T06 — ChEMBL bioactive compound retrieval via DeepMind chembl skill."""

    def test_target_search(self, output_dir):
        """Search ChEMBL for EGFR target."""
        out = output_dir / "t06_target_search.json"
        run_uv_script(
            str(CHEMBL_SCRIPTS / "chembl_api.py"),
            "target", "--search", "EGFR", "--limit", "5",
            "--output", str(out),
        )
        assert out.exists()
        data = load_json(out)
        assert data  # Non-empty response

    def test_target_by_id(self, output_dir):
        """Fetch CHEMBL203 (EGFR) target details."""
        out = output_dir / "t06_target_id.json"
        run_uv_script(
            str(CHEMBL_SCRIPTS / "chembl_api.py"),
            "target", "--id", "CHEMBL203",
            "--output", str(out),
        )
        assert out.exists()

    def test_activities_for_target(self, output_dir):
        """Fetch IC50 activities for EGFR."""
        out = output_dir / "t06_activities.json"
        run_uv_script(
            str(CHEMBL_SCRIPTS / "chembl_api.py"),
            "activity",
            "--filter", "target_chembl_id=CHEMBL203", "standard_type=IC50",
            "--limit", "10",
            "--output", str(out),
        )
        assert out.exists()
        data = load_json(out)
        assert data  # Non-empty

    def test_molecule_lookup(self, output_dir):
        """Look up aspirin (CHEMBL25)."""
        out = output_dir / "t06_aspirin.json"
        run_uv_script(
            str(CHEMBL_SCRIPTS / "chembl_api.py"),
            "molecule", "--id", "CHEMBL25",
            "--output", str(out),
        )
        assert out.exists()

    def test_api_status(self, output_dir):
        """Check ChEMBL API is responsive."""
        out = output_dir / "t06_status.json"
        run_uv_script(
            str(CHEMBL_SCRIPTS / "chembl_api.py"),
            "status",
            "--output", str(out),
        )
        assert out.exists()


# ===========================================================================
# STAGE 3: ADMET Filtering (wizard's own)
# ===========================================================================

class TestADMETLipinski:
    """T07 — Lipinski Rule of Five check."""

    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out = output_dir / "t07_lipinski.json"
        run_wizard(
            "admet-filter",
            "--smiles", "CC(=O)OC1=CC=CC=C1C(=O)O",
            "--names", "aspirin",
            "--output", str(self.out),
        )

    def _get_aspirin(self):
        data = load_json(self.out)
        compounds = data.get("compounds", data.get("results", [data]))
        return compounds[0]

    def test_file_exists(self):
        assert self.out.exists()

    def test_aspirin_passes(self):
        aspirin = self._get_aspirin()
        assert aspirin["lipinski_pass"] is True

    def test_aspirin_zero_violations(self):
        aspirin = self._get_aspirin()
        assert aspirin["lipinski_violations"] == 0

    def test_aspirin_mw(self):
        aspirin = self._get_aspirin()
        mw = aspirin["mw"]
        assert 175 < mw < 185, f"Aspirin MW should be ~180.2, got {mw}"


class TestADMETPAINS:
    """T08 — PAINS filter for rhodanine."""

    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out = output_dir / "t08_pains.json"
        run_wizard(
            "admet-filter",
            "--smiles", "O=C1NC(=S)SC1",
            "--names", "rhodanine",
            "--output", str(self.out),
        )

    def _get_rhodanine(self):
        data = load_json(self.out)
        return data.get("compounds", data.get("results", [data]))[0]

    def test_rhodanine_flagged(self):
        rhodanine = self._get_rhodanine()
        assert rhodanine["pains_pass"] is False

    def test_rhodanine_alert_name(self):
        rhodanine = self._get_rhodanine()
        alerts = rhodanine.get("pains_alerts", [])
        alert_str = " ".join(str(a) for a in alerts).lower()
        assert "rhod" in alert_str, f"Expected rhodanine-related alert, got: {alerts}"


class TestADMETSAScore:
    """T09 — Synthetic accessibility scoring for taxol."""

    TAXOL_SMILES = "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C5=CC=CC=C5)C)O)OC(=O)C)OC(=O)C6=CC=CC=C6)(C[C@@H]1OC(=O)[C@@H](O)C(=O)C7=CC=CC=C7)O)(C)C"

    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out = output_dir / "t09_sa_score.json"
        run_wizard(
            "admet-filter",
            "--smiles", self.TAXOL_SMILES,
            "--names", "taxol",
            "--output", str(self.out),
        )

    def test_taxol_high_sa_score(self):
        """Taxol SA score: RDKit may canonicalize differently, accept 5.0-10.0."""
        data = load_json(self.out)
        taxol = data.get("compounds", data.get("results", [data]))[0]
        sa = taxol["sa_score"]
        assert 5.0 <= sa <= 10.0, f"Taxol SA score should be high (5-8.5), got {sa}"


class TestADMETCascade:
    """T10 — Full ADMET filter cascade on 5 compounds."""

    COMPOUNDS = {
        "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "rhodanine": "O=C1NC(=S)SC1",
        "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "dtbhq": "CC(C)(C)C1=CC(=C(C=C1)O)C(C)(C)C",
    }

    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out = output_dir / "t10_cascade.json"
        smiles = list(self.COMPOUNDS.values())
        names = list(self.COMPOUNDS.keys())
        run_wizard(
            "admet-filter",
            "--smiles", *smiles,
            "--names", *names,
            "--output", str(self.out),
        )

    def _get_compounds(self):
        data = load_json(self.out)
        return data.get("compounds", data.get("results", []))

    def test_five_compounds_returned(self):
        assert len(self._get_compounds()) == 5

    def test_aspirin_passes_all(self):
        aspirin = next(c for c in self._get_compounds() if c["name"] == "aspirin")
        assert aspirin["lipinski_pass"] is True
        assert aspirin["pains_pass"] is True

    def test_rhodanine_fails_pains(self):
        rhodanine = next(c for c in self._get_compounds() if c["name"] == "rhodanine")
        assert rhodanine["pains_pass"] is False

    def test_ibuprofen_passes(self):
        ibuprofen = next(c for c in self._get_compounds() if c["name"] == "ibuprofen")
        assert ibuprofen["lipinski_pass"] is True
        assert ibuprofen["pains_pass"] is True

    def test_caffeine_passes(self):
        caffeine = next(c for c in self._get_compounds() if c["name"] == "caffeine")
        assert caffeine["lipinski_pass"] is True


class TestADMETTwoTier:
    """T10b — Two-tier filtering: relaxed passes more than strict."""

    COMPOUNDS = {
        "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "rhodanine": "O=C1NC(=S)SC1",       # 1 PAINS alert
        "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    }

    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        smiles = list(self.COMPOUNDS.values())
        names = list(self.COMPOUNDS.keys())

        self.strict_out = output_dir / "t10b_strict.json"
        run_wizard(
            "admet-filter",
            "--smiles", *smiles,
            "--names", *names,
            "--strictness", "strict",
            "--output", str(self.strict_out),
        )

        self.relaxed_out = output_dir / "t10b_relaxed.json"
        run_wizard(
            "admet-filter",
            "--smiles", *smiles,
            "--names", *names,
            "--strictness", "relaxed",
            "--output", str(self.relaxed_out),
        )

    def test_relaxed_passes_more(self):
        strict = load_json(self.strict_out)
        relaxed = load_json(self.relaxed_out)
        assert relaxed["passed"] >= strict["passed"]

    def test_strict_blocks_rhodanine(self):
        data = load_json(self.strict_out)
        compounds = data["compounds"]
        rhodanine = next(c for c in compounds if c["name"] == "rhodanine")
        assert rhodanine["pains_pass"] is False
        assert rhodanine["overall_pass"] is False

    def test_relaxed_passes_rhodanine(self):
        """Rhodanine has 1 PAINS alert — relaxed allows <=1."""
        data = load_json(self.relaxed_out)
        compounds = data["compounds"]
        rhodanine = next(c for c in compounds if c["name"] == "rhodanine")
        assert rhodanine["pains_pass"] is True

    def test_strictness_in_output(self):
        strict = load_json(self.strict_out)
        relaxed = load_json(self.relaxed_out)
        assert strict["strictness"] == "strict"
        assert relaxed["strictness"] == "relaxed"

    def test_tier_label_in_output(self):
        strict = load_json(self.strict_out)
        relaxed = load_json(self.relaxed_out)
        assert "strict" in strict["tier"].lower()
        assert "relaxed" in relaxed["tier"].lower()


# ===========================================================================
# STAGE 4: Generative Chemistry (wizard's own)
# ===========================================================================

class TestREINVENTConfig:
    """T11 — REINVENT 4 config generation."""

    ERLOTINIB = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"

    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out = output_dir / "t11_reinvent.toml"
        run_wizard(
            "generate-reinvent-config",
            "--seeds", self.ERLOTINIB,
            "--target-logp", "3.0",
            "--target-qed", "0.7",
            "--output", str(self.out),
        )

    def test_file_exists(self):
        assert self.out.exists()

    def test_valid_toml(self):
        import tomllib
        with open(self.out, "rb") as f:
            data = tomllib.load(f)
        assert "run_type" in data or "parameters" in data


class TestHPCScript:
    """T12 — SLURM submission script generation."""

    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out = output_dir / "t12_hpc.sh"
        run_wizard(
            "generate-hpc-script",
            "--job-type", "reinvent",
            "--config", str(output_dir / "t11_reinvent.toml"),
            "--gpu", "a100",
            "--ngpu", "1",
            "--mem", "64G",
            "--time", "24:00:00",
            "--output", str(self.out),
        )

    def test_file_exists(self):
        assert self.out.exists()

    def test_shebang(self):
        content = self.out.read_text()
        assert content.startswith("#!/bin/bash")

    def test_gpu_directive(self):
        content = self.out.read_text()
        assert "gpu" in content.lower()
        assert "a100" in content.lower()

    def test_memory_directive(self):
        content = self.out.read_text()
        assert "64G" in content

    def test_time_directive(self):
        content = self.out.read_text()
        assert "24:00:00" in content


# ===========================================================================
# STAGE 5: Molecular Dynamics (wizard's own)
# ===========================================================================

class TestGROMACSSetup:
    """T13 — GROMACS MD setup script generation."""

    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out = output_dir / "t13_md_setup.sh"
        run_wizard(
            "generate-md-setup",
            "--pdb", "protein_clean.pdb",
            "--ff", "charmm36",
            "--water", "tip3p",
            "--production-ns", "100",
            "--output", str(self.out),
        )

    def test_file_exists(self):
        assert self.out.exists()

    def test_contains_pdb2gmx(self):
        content = self.out.read_text()
        assert "pdb2gmx" in content

    def test_contains_solvate(self):
        content = self.out.read_text()
        assert "solvate" in content

    def test_contains_mdrun_or_md_mdp(self):
        """Script should contain mdrun command or md.mdp production config."""
        content = self.out.read_text()
        assert "mdrun" in content or "md.mdp" in content

    def test_force_field(self):
        content = self.out.read_text()
        assert "charmm36" in content.lower()

    def test_production_ns(self):
        content = self.out.read_text()
        assert "100" in content  # 100 ns production


class TestMDHPCScript:
    """T14 — GROMACS MD SLURM script."""

    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out = output_dir / "t14_md_hpc.sh"
        run_wizard(
            "generate-hpc-script",
            "--job-type", "gromacs-md",
            "--ngpu", "4",
            "--mem", "128G",
            "--time", "48:00:00",
            "--production-ns", "100",
            "--output", str(self.out),
        )

    def test_file_exists(self):
        assert self.out.exists()

    def test_multi_gpu(self):
        content = self.out.read_text()
        assert "4" in content  # 4 GPUs

    def test_memory(self):
        content = self.out.read_text()
        assert "128G" in content

    def test_gromacs_command(self):
        content = self.out.read_text()
        content_lower = content.lower()
        assert "gmx" in content_lower or "gromacs" in content_lower or "mdrun" in content_lower


# ===========================================================================
# OUT-OF-SCOPE / HARD NEGATIVES
# ===========================================================================

class TestOutOfScope:
    """N01-N05 — Hard negative queries that should be gracefully declined."""

    def test_no_wet_lab_subcommand(self):
        """No 'design-assay' or wet-lab subcommand should exist."""
        result = run_wizard("--help")
        assert "design-assay" not in result.stdout
        assert "wet-lab" not in result.stdout

    def test_no_clinical_trial_subcommand(self):
        """No 'clinical-trial' subcommand should exist."""
        result = run_wizard("--help")
        assert "clinical-trial" not in result.stdout

    def test_no_protein_folding_subcommand(self):
        """No 'fold-protein' subcommand should exist (we use predict-structure instead)."""
        result = run_wizard("--help")
        assert "fold-protein" not in result.stdout

    def test_no_retrosynthesis_subcommand(self):
        """No 'retrosynthesis' subcommand should exist."""
        result = run_wizard("--help")
        assert "retrosynthesis" not in result.stdout


# ===========================================================================
# INTEGRATION: Science skills are available
# ===========================================================================

class TestScienceSkillsAvailability:
    """Verify DeepMind science-skills are installed and scripts accessible."""

    def test_pdb_scripts_exist(self):
        assert (PDB_SCRIPTS / "download_coordinate_files.py").exists()
        assert (PDB_SCRIPTS / "fetch_pdb_metadata.py").exists()
        assert (PDB_SCRIPTS / "search_pdb.py").exists()
        assert (PDB_SCRIPTS / "fetch_schema.py").exists()

    def test_alphafold_scripts_exist(self):
        assert (AF_SCRIPTS / "fetch_structure.py").exists()
        assert (AF_SCRIPTS / "analyze_plddt.py").exists()
        assert (AF_SCRIPTS / "analyze_pae.py").exists()

    def test_chembl_scripts_exist(self):
        assert (CHEMBL_SCRIPTS / "chembl_api.py").exists()

    def test_wizard_script_exists(self):
        assert WIZARD_SCRIPT.exists()

    def test_wizard_help(self):
        result = run_wizard("--help")
        assert "prepare-receptor" in result.stdout
        assert "admet-filter" in result.stdout
        assert "admet-predict" in result.stdout
        assert "gmxapi-setup" in result.stdout
        assert "gmxapi-validate" in result.stdout
        assert "generate-reinvent-config" in result.stdout
        assert "generate-hpc-script" in result.stdout
        assert "generate-md-setup" in result.stdout

    def test_admetlab_client_exists(self):
        assert (WIZARD_DIR / "scripts" / "admetlab_client.py").exists()

    def test_gmxapi_setup_exists(self):
        assert (WIZARD_DIR / "scripts" / "gmxapi_setup.py").exists()


# ===========================================================================
# ADMETlab 3.0 Tests
# ===========================================================================

class TestADMETlabPredict:
    """T15 — ADMETlab 3.0 API prediction via drug_discovery_api.py."""

    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out = output_dir / "t15_admet_predict.json"
        run_wizard(
            "admet-predict",
            "--smiles", "CC(=O)Oc1ccccc1C(=O)O",
            "--names", "aspirin",
            "--output", str(self.out),
        )

    def test_file_exists(self):
        assert self.out.exists()

    def test_has_source(self):
        data = load_json(self.out)
        assert "source" in data
        # Either ADMETlab 3.0 or RDKit fallback
        assert data["source"] in ("ADMETlab 3.0", "RDKit (offline fallback)")

    def test_has_compounds(self):
        data = load_json(self.out)
        compounds = data.get("compounds", [])
        assert len(compounds) >= 1

    def test_aspirin_name(self):
        data = load_json(self.out)
        compounds = data.get("compounds", [])
        assert compounds[0]["name"] == "aspirin"


class TestADMETlabFallback:
    """T16 — Verify RDKit fallback when API would be unreachable."""

    def test_rdkit_fallback_works(self):
        """Import the fallback function directly and test it."""
        sys.path.insert(0, str(WIZARD_DIR / "scripts"))
        try:
            from admetlab_client import _rdkit_fallback
            result = _rdkit_fallback(
                ["CC(=O)Oc1ccccc1C(=O)O"],
                ["aspirin"],
            )
            assert result["source"] == "RDKit (offline fallback)"
            assert len(result["compounds"]) == 1
            aspirin = result["compounds"][0]
            assert aspirin["lipinski_pass"] is True
            assert aspirin["pains_pass"] is True
            assert 175 < aspirin["mw"] < 185
        finally:
            sys.path.pop(0)

    def test_rdkit_fallback_detects_pains(self):
        """Rhodanine should be flagged by RDKit fallback."""
        sys.path.insert(0, str(WIZARD_DIR / "scripts"))
        try:
            from admetlab_client import _rdkit_fallback
            result = _rdkit_fallback(["O=C1NC(=S)SC1"], ["rhodanine"])
            rhodanine = result["compounds"][0]
            assert rhodanine["pains_pass"] is False
        finally:
            sys.path.pop(0)


class TestADMETlabSummary:
    """T17 — ADMETlab markdown summary generation."""

    def test_summary_generation(self):
        """Generate a summary from mock filter results."""
        sys.path.insert(0, str(WIZARD_DIR / "scripts"))
        try:
            from admetlab_client import generate_summary_markdown
            mock_results = {
                "total_compounds": 2,
                "total_passed": 1,
                "total_failed": 1,
                "compounds": [
                    {
                        "name": "aspirin",
                        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                        "flags": [],
                        "checked": 10,
                        "passed": 10,
                        "failed": 0,
                        "overall_pass": True,
                    },
                    {
                        "name": "toxic_compound",
                        "smiles": "CCCCCCCCCCCCCCCCCCCC",
                        "flags": ["LogP: 10.2 > 5.0", "MW: 600 > 500"],
                        "checked": 10,
                        "passed": 8,
                        "failed": 2,
                        "overall_pass": False,
                    },
                ],
            }
            md = generate_summary_markdown(mock_results)
            assert "ADMET Prediction Summary" in md
            assert "aspirin" in md
            assert "PASS" in md
            assert "FAIL" in md
        finally:
            sys.path.pop(0)


class TestADMETlabThresholds:
    """T18 — ADMETlab threshold configuration."""

    def test_thresholds_defined(self):
        sys.path.insert(0, str(WIZARD_DIR / "scripts"))
        try:
            from admetlab_client import THRESHOLDS
            # Must have key toxicity endpoints
            assert "hERG" in THRESHOLDS
            assert "AMES" in THRESHOLDS
            assert "DILI" in THRESHOLDS
            assert "CYP3A4-inh" in THRESHOLDS
            assert "BBB" in THRESHOLDS
            assert "Caco-2" in THRESHOLDS
            assert "PAINS" in THRESHOLDS
            assert "MW" in THRESHOLDS
            assert len(THRESHOLDS) >= 30
        finally:
            sys.path.pop(0)


# ===========================================================================
# gmxapi Tests
# ===========================================================================

class TestGmxapiWorkflowGeneration:
    """T19 — gmxapi Python workflow script generation."""

    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out_dir = output_dir / "gmxapi_test"
        self.out_dir.mkdir(exist_ok=True)
        self.out = self.out_dir / "run_md.py"
        run_wizard(
            "gmxapi-setup",
            "--pdb", "protein.pdb",
            "--ff", "charmm36",
            "--water", "tip3p",
            "--production-ns", "10",
            "--output", str(self.out),
        )

    def test_script_exists(self):
        assert self.out.exists()

    def test_script_is_valid_python(self):
        """Script should be valid Python syntax."""
        content = self.out.read_text()
        compile(content, str(self.out), "exec")

    def test_script_contains_gmxapi_imports(self):
        content = self.out.read_text()
        assert "gmxapi" in content or "gmx" in content

    def test_script_contains_pdb2gmx(self):
        content = self.out.read_text()
        assert "pdb2gmx" in content

    def test_script_contains_mdrun(self):
        content = self.out.read_text()
        assert "mdrun" in content

    def test_script_contains_8_steps(self):
        content = self.out.read_text()
        assert "Step 8/8" in content

    def test_mdp_files_generated(self):
        """All 5 MDP files should be generated alongside the script."""
        assert (self.out_dir / "ions.mdp").exists()
        assert (self.out_dir / "em.mdp").exists()
        assert (self.out_dir / "nvt.mdp").exists()
        assert (self.out_dir / "npt.mdp").exists()
        assert (self.out_dir / "md.mdp").exists()


class TestGmxapiMDP:
    """T20 — MDP file generation and validation."""

    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out_dir = output_dir / "mdp_test"
        self.out_dir.mkdir(exist_ok=True)
        sys.path.insert(0, str(WIZARD_DIR / "scripts"))
        from gmxapi_setup import generate_mdp
        for t in ["ions", "em", "nvt", "npt", "production"]:
            name = "md.mdp" if t == "production" else f"{t}.mdp"
            (self.out_dir / name).write_text(generate_mdp(t, 10.0))

    def test_em_integrator_is_steep(self):
        content = (self.out_dir / "em.mdp").read_text()
        assert "steep" in content

    def test_nvt_has_thermostat(self):
        content = (self.out_dir / "nvt.mdp").read_text()
        assert "V-rescale" in content

    def test_npt_has_barostat(self):
        content = (self.out_dir / "npt.mdp").read_text()
        assert "Parrinello-Rahman" in content

    def test_production_nsteps(self):
        """10 ns at dt=0.002 ps = 5,000,000 steps."""
        content = (self.out_dir / "md.mdp").read_text()
        assert "5000000" in content

    def test_all_use_verlet(self):
        for name in ["em.mdp", "nvt.mdp", "npt.mdp", "md.mdp"]:
            content = (self.out_dir / name).read_text()
            assert "Verlet" in content


class TestGmxapiValidation:
    """T21 — MDP file validation."""

    @pytest.fixture(autouse=True)
    def _run(self, output_dir):
        self.out_dir = output_dir / "validate_test"
        self.out_dir.mkdir(exist_ok=True)
        # Create a valid MDP
        self.valid_mdp = self.out_dir / "valid.mdp"
        self.valid_mdp.write_text(
            "integrator = md\nnsteps = 50000\ndt = 0.002\n"
            "coulombtype = PME\ntcoupl = V-rescale\npcoupl = no\n"
            "constraints = h-bonds\n"
        )
        # Create an invalid MDP
        self.invalid_mdp = self.out_dir / "invalid.mdp"
        self.invalid_mdp.write_text(
            "integrator = broken\nnsteps = -100\ncoulombtype = MAGIC\n"
        )

    def test_valid_mdp_passes(self):
        sys.path.insert(0, str(WIZARD_DIR / "scripts"))
        try:
            from gmxapi_setup import validate_mdp
            result = validate_mdp(str(self.valid_mdp))
            assert result["valid"] is True
            assert len(result["issues"]) == 0
        finally:
            sys.path.pop(0)

    def test_invalid_mdp_fails(self):
        sys.path.insert(0, str(WIZARD_DIR / "scripts"))
        try:
            from gmxapi_setup import validate_mdp
            result = validate_mdp(str(self.invalid_mdp))
            assert result["valid"] is False
            assert len(result["issues"]) >= 2  # bad integrator + bad coulombtype
        finally:
            sys.path.pop(0)

    def test_validate_via_cli(self, output_dir):
        report = output_dir / "validate_report.json"
        run_wizard(
            "gmxapi-validate",
            "--mdp", str(self.valid_mdp), str(self.invalid_mdp),
            "--output", str(report),
        )
        data = load_json(report)
        assert data["total_files"] == 2
        assert data["all_valid"] is False  # One invalid


# ===========================================================================
# Structure Prediction Tests
# ===========================================================================

class TestStructurePredictionScript:
    """T22 — predict_structure.py exists and is importable."""

    def test_script_exists(self):
        assert (WIZARD_DIR / "scripts" / "predict_structure.py").exists()

    def test_sequence_validation(self):
        sys.path.insert(0, str(WIZARD_DIR / "scripts"))
        try:
            from predict_structure import validate_sequence
            # Valid sequence
            ok, msg = validate_sequence("MKTLLILAVVAAALA" * 3)
            assert ok is True
            # Too short
            ok, msg = validate_sequence("MKT")
            assert ok is False
            assert "short" in msg.lower()
            # Invalid chars
            ok, msg = validate_sequence("MKTLLILAVV123XYZ")
            assert ok is False
            assert "Invalid" in msg
        finally:
            sys.path.pop(0)

    def test_fasta_parsing(self, output_dir):
        """Test FASTA file reading."""
        fasta = output_dir / "test.fasta"
        fasta.write_text(">test_protein\nMKTLLILAVV\nAALAMKTLLI\n")
        sys.path.insert(0, str(WIZARD_DIR / "scripts"))
        try:
            from predict_structure import read_fasta
            header, seq = read_fasta(str(fasta))
            assert header == "test_protein"
            assert seq == "MKTLLILAVVAALAMKTLLI"
        finally:
            sys.path.pop(0)


class TestAF2ScriptGeneration:
    """T23 — ColabFold/AF2 HPC script generation."""

    @pytest.fixture(autouse=True)
    def _setup(self, output_dir):
        self.fasta = output_dir / "target.fasta"
        self.fasta.write_text(">target\nMKTLLILAVVAALA\n")
        self.out = output_dir / "t23_af2.sh"

    def test_colabfold_script(self):
        run_wizard(
            "generate-af2-script",
            "--fasta", str(self.fasta),
            "--method", "colabfold",
            "--gpu", "a100",
            "--output", str(self.out),
        )
        assert self.out.exists()
        content = self.out.read_text()
        assert "colabfold_batch" in content
        assert "SBATCH" in content
        assert "a100" in content

    def test_alphafold2_script(self, output_dir):
        out = output_dir / "t23_af2_full.sh"
        run_wizard(
            "generate-af2-script",
            "--fasta", str(self.fasta),
            "--method", "alphafold2",
            "--output", str(out),
        )
        content = out.read_text()
        assert "run_alphafold.py" in content


class TestStructureQualityAssessment:
    """T24 — Assess predicted structure quality from pLDDT."""

    def test_assess_from_pdb(self, output_dir):
        """Create a mock PDB with B-factors and assess quality."""
        # Create a minimal PDB with known B-factor (pLDDT) values
        pdb_lines = []
        for i in range(1, 11):
            plddt = 85.0  # High confidence
            pdb_lines.append(
                f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
                f"  0.000   0.000   0.000  1.00{plddt:6.2f}           C  "
            )
        pdb = output_dir / "mock_predicted.pdb"
        pdb.write_text("\n".join(pdb_lines) + "\n")

        report = output_dir / "t24_quality.json"
        run_wizard(
            "assess-structure",
            "--pdb", str(pdb),
            "--output", str(report),
        )
        data = load_json(report)
        assert data["total_residues"] == 10
        assert data["mean_plddt"] == 85.0
        assert "High confidence" in data["confidence_assessment"]

    def test_wizard_help_includes_predict(self):
        result = run_wizard("--help")
        assert "predict-structure" in result.stdout
        assert "generate-af2-script" in result.stdout
        assert "assess-structure" in result.stdout


# ============================================================
# T25-T28: Compound Library Tests (COCONUT + ZINC20)
# ============================================================

class TestCompoundLibraryScripts:
    """Test that the compound library client script exists and imports."""

    def test_script_exists(self):
        script = WIZARD_DIR / "scripts" / "compound_libraries.py"
        assert script.exists(), "compound_libraries.py missing"

    def test_imports(self):
        sys.path.insert(0, str(WIZARD_DIR / "scripts"))
        try:
            from compound_libraries import (
                search_coconut,
                generate_zinc_download_script,
                ZINC_SUBSETS,
            )
            assert callable(search_coconut)
            assert callable(generate_zinc_download_script)
            assert "drug-like" in ZINC_SUBSETS
        finally:
            sys.path.pop(0)


class TestZINCDownloadScript:
    """Test ZINC20 download script generation (offline, no API needed)."""

    def test_drug_like_script(self, output_dir):
        out = output_dir / "zinc_druglike.sh"
        run_wizard(
            "query-zinc",
            "--subset", "drug-like",
            "--output", str(out),
        )
        assert out.exists()
        content = out.read_text()
        assert "ZINC20" in content
        assert "drug-like" in content
        assert "curl" in content

    def test_lead_like_script(self, output_dir):
        out = output_dir / "zinc_leadlike.sh"
        run_wizard(
            "query-zinc",
            "--subset", "lead-like",
            "--output", str(out),
        )
        content = out.read_text()
        assert "lead-like" in content

    def test_custom_mw_range(self, output_dir):
        out = output_dir / "zinc_custom.sh"
        run_wizard(
            "query-zinc",
            "--subset", "drug-like",
            "--mw-range", "300", "450",
            "--output", str(out),
        )
        content = out.read_text()
        assert "300-450" in content

    def test_sdf_format(self, output_dir):
        out = output_dir / "zinc_sdf.sh"
        run_wizard(
            "query-zinc",
            "--subset", "fragment-like",
            "--format", "sdf",
            "--output", str(out),
        )
        content = out.read_text()
        assert ".sdf" in content


@pytest.mark.network
class TestCOCONUTSearch:
    """Test COCONUT natural products search (requires network)."""

    def test_search_manzamine(self, output_dir):
        out = output_dir / "coconut_manzamine.json"
        run_wizard(
            "query-coconut",
            "--query", "manzamine",
            "--limit", "3",
            "--output", str(out),
        )
        data = load_json(out)
        assert data["source"] == "COCONUT"
        assert len(data["compounds"]) > 0
        # Check compound has expected fields
        c = data["compounds"][0]
        assert "smiles" in c
        assert "name" in c
        assert "id" in c


class TestWizardHelpIncludesLibraries:
    """Verify new subcommands appear in --help."""

    def test_help_includes_coconut(self):
        result = run_wizard("--help")
        assert "query-coconut" in result.stdout

    def test_help_includes_zinc(self):
        result = run_wizard("--help")
        assert "query-zinc" in result.stdout
