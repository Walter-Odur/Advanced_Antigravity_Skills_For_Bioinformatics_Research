"""Shared fixtures.

The default test run is hermetic: no network, no external binaries. Tests
that need either are marked and deselected by the ``addopts`` in
``pyproject.toml``, so a fresh checkout on another machine runs the whole
default suite with nothing installed beyond the package and its
dependencies.

Auto-installation of missing tools (pip/conda/apt) is disabled during
tests via the ``AGSKILLS_NO_AUTO_INSTALL`` environment variable.

Real data rather than synthetic stand-ins is used wherever the thing under
test is scientific:

``1M17.pdb``
    EGFR kinase with erlotinib bound. One chain, one drug-like ligand, so
    binding-site detection has an unambiguous right answer.
``6HEZ.pdb``
    *M. tuberculosis* DprE1 with both an FAD cofactor and an inhibitor
    (0SK), across two chains. This is the structure that exposes the
    cofactor-versus-ligand problem and the incomplete-side-chain problem.
``AF-P9WJA5-F1.pdb``
    A real AlphaFold model, so the B-factor column holds genuine pLDDT
    rather than crystallographic B-factors.
``pubchem_reference.json``
    Independently-computed physicochemical properties for ten marketed
    drugs, used to check the descriptor layer against an outside source
    instead of against itself.
``reinvent_staged_learning_1.csv``
    A staged-learning CSV in REINVENT 4's real column format
    (``Agent, Prior, Target, Score, SMILES, SMILES_state`` plus component
    columns), with deliberate duplicate sampling and a rising score trend.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

# Disable auto-installation of missing tools during tests.
os.environ["AGSKILLS_NO_AUTO_INSTALL"] = "1"

DATA_DIR = Path(__file__).parent / "data"


# ---------------------------------------------------------------------------
# Data fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def data_dir() -> Path:
    """The directory holding the real test data."""
    assert DATA_DIR.is_dir(), f"missing test data directory: {DATA_DIR}"
    return DATA_DIR


@pytest.fixture(scope="session")
def pdb_1m17(data_dir: Path) -> Path:
    """EGFR kinase domain with erlotinib (ligand AQ4) bound."""
    return data_dir / "1M17.pdb"


@pytest.fixture(scope="session")
def pdb_6hez(data_dir: Path) -> Path:
    """M. tuberculosis DprE1 with FAD and inhibitor 0SK, two chains."""
    return data_dir / "6HEZ.pdb"


@pytest.fixture(scope="session")
def alphafold_model(data_dir: Path) -> Path:
    """A real AlphaFold model, with genuine pLDDT in the B-factor column."""
    return data_dir / "AF-P9WJA5-F1.pdb"


@pytest.fixture(scope="session")
def reinvent_csv(data_dir: Path) -> Path:
    """A staged-learning CSV in REINVENT 4's real output format."""
    return data_dir / "reinvent_staged_learning_1.csv"


@pytest.fixture(scope="session")
def vina_output(data_dir: Path) -> str:
    """Authentic AutoDock Vina v1.2.5 stdout, including its banner."""
    return (data_dir / "vina_output.txt").read_text(encoding="utf-8")


@pytest.fixture(scope="session")
def reference_properties(data_dir: Path) -> dict:
    """PubChem-computed properties for ten marketed drugs.

    Keyed by drug name. Each entry carries ``SMILES``,
    ``MolecularWeight``, ``XLogP``, ``TPSA``, ``HBondDonorCount``,
    ``HBondAcceptorCount`` and ``RotatableBondCount``.
    """
    return json.loads((data_dir / "pubchem_reference.json").read_text(
        encoding="utf-8"))


# ---------------------------------------------------------------------------
# Chemistry fixtures
# ---------------------------------------------------------------------------


#: Marketed drugs with properties worth asserting on individually.
DRUG_SMILES: dict[str, str] = {
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "imatinib": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
    "bedaquiline": ("CCN(CC)CC[C@@H](O)[C@@H](c1ccc2ccccc2n1)c1cc2ccccc2cc1OC"),
    "rifampicin": ("C[C@H]1/C=C/C=C(\\C)C(=O)Nc2c(O)c3c(c(O)c2/C=N/N2CCN(C)CC2)"
                   "C(=O)[C@](C)(O3)O/C=C/[C@@H](OC)[C@H](C)[C@@H](OC(C)=O)"
                   "[C@H](C)[C@H](O)[C@H](C)C1"),
    "penicillin_g": "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",
    "ethanol": "CCO",
    "benzene": "c1ccccc1",
}

#: Structures that must fail to parse, each for a different reason.
INVALID_SMILES: list[str] = [
    "",
    "   ",
    "NOT_A_MOLECULE",
    "C(C(C",            # unbalanced parentheses
    "c1ccccc",          # unclosed ring
    "[Xx]",             # not an element
    "C1CC",             # unclosed ring bond
]


@pytest.fixture
def drug_smiles() -> dict[str, str]:
    return dict(DRUG_SMILES)


@pytest.fixture
def aspirin():
    """A parsed aspirin molecule."""
    from agskills.chem.smiles import parse_smiles
    mol = parse_smiles(DRUG_SMILES["aspirin"])
    assert mol is not None
    return mol


@pytest.fixture
def records_from_drugs():
    """Build CompoundRecords for a named subset of the drug panel."""
    from agskills.chem.smiles import CompoundRecord

    def _build(*names: str):
        chosen = names or tuple(DRUG_SMILES)
        return [CompoundRecord(name=n, smiles=DRUG_SMILES[n]) for n in chosen]

    return _build


# ---------------------------------------------------------------------------
# Environment isolation
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def isolated_cache(tmp_path, monkeypatch):
    """Point the HTTP cache at a per-test directory.

    Without this, the rate-limiter and robots.txt caches are shared with the
    developer's real runs, so one test could make another sleep.
    """
    cache = tmp_path / "cache"
    cache.mkdir()
    monkeypatch.setenv("AGSKILLS_CACHE_DIR", str(cache))
    monkeypatch.setenv("AGSKILLS_USER_AGENT",
                       "agskills-tests/2.0 (+automated test suite)")
    return cache


@pytest.fixture(autouse=True)
def no_reinvent_env(monkeypatch, request):
    """Clear REINVENT environment variables unless a test sets them.

    Several tests assert on discovery behaviour, which a developer's own
    ``REINVENT_DIR`` would otherwise change.
    """
    if "keep_reinvent_env" in request.keywords:
        return
    monkeypatch.delenv("REINVENT_DIR", raising=False)
    monkeypatch.delenv("REINVENT_PRIOR_BASE", raising=False)


# ---------------------------------------------------------------------------
# Capability detection for marked tests
# ---------------------------------------------------------------------------


def _have_vina() -> bool:
    from agskills.docking.vina import find_vina
    return find_vina()[0] != "none"


def _have_gromacs() -> bool:
    return bool(shutil.which("gmx") or shutil.which("gmx_mpi"))


def _have_reinvent() -> bool:
    import importlib.util
    return importlib.util.find_spec("reinvent") is not None


def pytest_collection_modifyitems(config, items):
    """Skip capability-marked tests when the capability is absent.

    The markers are deselected by default, so this only matters when a run
    opts into them explicitly (``-m requires_vina``). Skipping rather than
    failing means the opt-in run is still informative on a machine that
    lacks the tool.
    """
    checks = {
        "requires_vina": (_have_vina, "AutoDock Vina is not installed"),
        "requires_gromacs": (_have_gromacs, "GROMACS is not installed"),
        "requires_reinvent": (_have_reinvent,
                              "the reinvent package is not importable"),
    }
    cache: dict[str, bool] = {}
    for item in items:
        for marker, (probe, reason) in checks.items():
            if marker in item.keywords:
                if marker not in cache:
                    cache[marker] = probe()
                if not cache[marker]:
                    item.add_marker(pytest.mark.skip(reason=reason))


@pytest.fixture
def bash() -> str:
    """Path to a bash interpreter, for syntax-checking generated scripts.

    On Windows we use WSL's ``bash.exe``.  The generated scripts target
    Linux clusters, so WSL is the right environment to validate them in.
    """
    if sys.platform == "win32":
        wsl = shutil.which("wsl") or shutil.which("wsl.exe")
        if not wsl:
            pytest.skip("WSL is not installed; cannot syntax-check bash scripts")
        return wsl  # caller uses _wsl_bash_check() below
    found = shutil.which("bash")
    if not found:
        pytest.skip("bash is not available to syntax-check generated scripts")
    return found


def _win_to_wsl_path(win_path: str | Path) -> str:
    """Convert ``C:\\Users\\...\\file.sh`` to ``/mnt/c/Users/.../file.sh``."""
    p = str(win_path).replace("\\", "/")
    # "E:/foo/bar" -> "/mnt/e/foo/bar"
    if len(p) >= 2 and p[1] == ":":
        drive = p[0].lower()
        p = f"/mnt/{drive}{p[2:]}"
    return p


def bash_check(bash_exe: str, script_path: Path) -> subprocess.CompletedProcess:
    """Run ``bash -n`` on a script, handling Windows/WSL path conversion."""
    if sys.platform == "win32":
        wsl_path = _win_to_wsl_path(script_path)
        return subprocess.run(
            [bash_exe, "bash", "-n", wsl_path],
            capture_output=True, text=True, check=False,
        )
    return subprocess.run(
        [bash_exe, "-n", str(script_path)],
        capture_output=True, text=True, check=False,
    )

