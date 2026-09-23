"""REINVENT4 runner internals.

Provides helper functions for interacting with the REINVENT4 codebase:
- Installation checking
- Prior model resolution and download
- Config generation (TOML)
- Subprocess execution
- Result parsing

This module is imported by compound_synthesis.py and should NOT be called
directly by the agent.
"""

from __future__ import annotations

import csv
import io
import json
import os
import platform
import shutil
import subprocess
import sys
import urllib.request
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ZENODO_RECORD = "15641296"
ZENODO_BASE = f"https://zenodo.org/api/records/{ZENODO_RECORD}"

PRIOR_REGISTRY = {
    ".reinvent": "reinvent.prior",
    ".libinvent": "libinvent.prior",
    ".linkinvent": "linkinvent.prior",
    ".m2m_high": "mol2mol_high_similarity.prior",
    ".m2m_medium": "mol2mol_medium_similarity.prior",
    ".m2m_mmp": "mol2mol_mmp.prior",
    ".m2m_scaffold": "mol2mol_scaffold.prior",
    ".m2m_scaffold_generic": "mol2mol_scaffold_generic.prior",
    ".pepinvent": "pepinvent.prior",
}

GENERATOR_PRIOR_MAP = {
    "reinvent": ".reinvent",
    "libinvent": ".libinvent",
    "linkinvent": ".linkinvent",
    "mol2mol": ".m2m_scaffold_generic",
    "pepinvent": ".pepinvent",
}

# SMARTS for common PAINS-like unwanted substructures
PAINS_SMARTS = [
    "[*;r{8-17}]",
    "[#8][#8]",
    "[#6;+]",
    "[#16][#16]",
    "[#7;!n][S;!$(S(=O)=O)]",
    "[#7;!n][#7;!n]",
    "C#C",
    "C(=[O,S])[O,S]",
    "[#7;!n][C;!$(C(=[O,N])[N,O])][#16;!s]",
    "[#7;!n][C;!$(C(=[O,N])[N,O])][#7;!n]",
    "[#7;!n][C;!$(C(=[O,N])[N,O])][#8;!o]",
    "[#8;!o][C;!$(C(=[O,N])[N,O])][#16;!s]",
    "[#8;!o][C;!$(C(=[O,N])[N,O])][#8;!o]",
    "[#16;!s][C;!$(C(=[O,N])[N,O])][#16;!s]",
]


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class SetupReport:
    """Result of check-setup."""
    python_version: str = ""
    torch_available: bool = False
    torch_version: str = ""
    rdkit_available: bool = False
    rdkit_version: str = ""
    reinvent_importable: bool = False
    reinvent_version: str = ""
    reinvent_cli: str = ""
    gpu_available: bool = False
    gpu_name: str = ""
    gpu_memory_mb: int = 0
    priors_found: list = field(default_factory=list)
    priors_missing: list = field(default_factory=list)
    ready: bool = False
    errors: list = field(default_factory=list)


@dataclass
class SeedReport:
    """Result of prepare-seeds."""
    total_input: int = 0
    valid: int = 0
    invalid: int = 0
    duplicates_removed: int = 0
    invalid_smiles: list = field(default_factory=list)
    output_file: str = ""


# ---------------------------------------------------------------------------
# Installation checking
# ---------------------------------------------------------------------------

def check_reinvent_installation(reinvent_dir: str) -> SetupReport:
    """Validate REINVENT4 installation and dependencies.

    Args:
        reinvent_dir: Path to the REINVENT4 repository.

    Returns:
        SetupReport with all check results.
    """
    report = SetupReport()
    report.python_version = platform.python_version()

    # Check torch
    try:
        import torch
        report.torch_available = True
        report.torch_version = torch.__version__
        if torch.cuda.is_available():
            report.gpu_available = True
            report.gpu_name = torch.cuda.get_device_name(0)
            report.gpu_memory_mb = torch.cuda.get_device_properties(0).total_mem // (1024 * 1024)
    except ImportError:
        report.errors.append("PyTorch not installed. Run: pip install torch==2.12.0")

    # Check RDKit
    try:
        import rdkit
        report.rdkit_available = True
        report.rdkit_version = rdkit.__version__
    except ImportError:
        report.errors.append("RDKit not installed. Run: pip install rdkit")

    # Check REINVENT importable
    rdir = Path(reinvent_dir)
    if rdir.exists():
        # Add to path temporarily
        str_rdir = str(rdir)
        if str_rdir not in sys.path:
            sys.path.insert(0, str_rdir)
        try:
            from reinvent import version
            report.reinvent_importable = True
            report.reinvent_version = version.__version__
        except ImportError as e:
            report.errors.append(f"Cannot import reinvent: {e}")

        # Check for reinvent CLI
        reinvent_cli = shutil.which("reinvent")
        if reinvent_cli:
            report.reinvent_cli = reinvent_cli
        else:
            # Try in the Python Scripts dir
            scripts_dir = Path(sys.executable).parent / ("Scripts" if os.name == "nt" else "")
            candidate = scripts_dir / ("reinvent.exe" if os.name == "nt" else "reinvent")
            if candidate.exists():
                report.reinvent_cli = str(candidate)
    else:
        report.errors.append(f"REINVENT4 directory not found: {reinvent_dir}")

    # Check prior models
    priors_dir = rdir / "priors"
    if not priors_dir.exists():
        env_base = os.environ.get("REINVENT_PRIOR_BASE", "")
        if env_base:
            priors_dir = Path(env_base)

    for key, filename in PRIOR_REGISTRY.items():
        fpath = priors_dir / filename
        if fpath.exists():
            report.priors_found.append(filename)
        else:
            report.priors_missing.append(filename)

    # Overall readiness
    report.ready = (
        report.torch_available
        and report.rdkit_available
        and report.reinvent_importable
        and len(report.priors_found) > 0
    )

    return report


# ---------------------------------------------------------------------------
# Prior model resolution
# ---------------------------------------------------------------------------

def resolve_prior(prior_key: str, reinvent_dir: str) -> str:
    """Resolve a prior key (e.g. '.reinvent') to an absolute path.

    If the key starts with '.', look up in the registry.
    Otherwise, treat as a direct file path.

    Args:
        prior_key: Prior model key or file path.
        reinvent_dir: Path to the REINVENT4 repository.

    Returns:
        Absolute path to the prior model file.
    """
    if prior_key.startswith("."):
        filename = PRIOR_REGISTRY.get(prior_key)
        if not filename:
            raise ValueError(
                f"Unknown prior key '{prior_key}'. "
                f"Available: {', '.join(PRIOR_REGISTRY.keys())}"
            )
        # Check standard locations
        candidates = [
            Path(reinvent_dir) / "priors" / filename,
            Path(os.environ.get("REINVENT_PRIOR_BASE", "")) / filename,
        ]
        for c in candidates:
            if c.exists():
                return str(c.resolve())

        # Auto-download from Zenodo if not found
        target_path = Path(reinvent_dir) / "priors" / filename
        model_short = prior_key.lstrip(".")
        print(f"  Prior model '{filename}' not found locally. Auto-downloading from Zenodo...")
        result = download_priors(reinvent_dir, [model_short])
        if result["downloaded"]:
            print(f"  Downloaded: {filename}")
            return str(target_path.resolve())
        elif result["errors"]:
            print(f"  WARNING: Download failed: {result['errors']}")
            print(f"  Using path anyway (manual download required): {target_path}")

        return str(target_path.resolve())
    else:
        return str(Path(prior_key).resolve())


# ---------------------------------------------------------------------------
# Smart device detection
# ---------------------------------------------------------------------------

# Minimum VRAM (MB) thresholds for different workloads
_VRAM_THRESHOLDS = {
    "staged_learning": 4000,     # 4 GB minimum for RL
    "transfer_learning": 3000,   # 3 GB minimum for TL
    "sampling": 2000,            # 2 GB minimum for sampling
    "scoring": 0,                # CPU is fine for scoring
}


def detect_device(
    run_mode: str = "staged_learning",
    num_steps: int = 100,
    batch_size: int = 128,
) -> dict:
    """Auto-detect the best device and whether HPC is recommended.

    Checks for GPU availability and VRAM. If GPU is available and has
    enough memory, uses cuda. Otherwise falls back to CPU. For large
    jobs, recommends HPC submission.

    Args:
        run_mode: REINVENT4 run mode.
        num_steps: Number of RL steps (for HPC recommendation).
        batch_size: Batch size (for VRAM estimation).

    Returns:
        Dict with:
            device: "cuda:0" or "cpu"
            gpu_available: bool
            gpu_name: str
            gpu_memory_mb: int
            recommend_hpc: bool
            hpc_reason: str (why HPC is recommended, if applicable)
    """
    result = {
        "device": "cpu",
        "gpu_available": False,
        "gpu_name": "",
        "gpu_memory_mb": 0,
        "recommend_hpc": False,
        "hpc_reason": "",
    }

    # Check for GPU
    try:
        import torch
        if torch.cuda.is_available():
            result["gpu_available"] = True
            result["gpu_name"] = torch.cuda.get_device_name(0)
            result["gpu_memory_mb"] = (
                torch.cuda.get_device_properties(0).total_mem // (1024 * 1024)
            )

            # Check if VRAM is sufficient for the workload
            min_vram = _VRAM_THRESHOLDS.get(run_mode, 2000)
            # Scale VRAM requirement by batch size
            estimated_vram = min_vram + (batch_size * 2)  # rough estimate

            if result["gpu_memory_mb"] >= estimated_vram:
                result["device"] = "cuda:0"
            else:
                result["device"] = "cpu"
                result["recommend_hpc"] = True
                result["hpc_reason"] = (
                    f"GPU {result['gpu_name']} has {result['gpu_memory_mb']} MB VRAM, "
                    f"but estimated requirement is ~{estimated_vram} MB for "
                    f"{run_mode} with batch_size={batch_size}. "
                    f"Recommend HPC with a larger GPU (A100/V100)."
                )
        else:
            # No GPU available
            result["device"] = "cpu"
    except ImportError:
        # torch not installed
        result["device"] = "cpu"

    # Check if the job is too large for local execution regardless of GPU
    if run_mode == "staged_learning":
        if num_steps > 200 and batch_size >= 128:
            if not result["gpu_available"]:
                result["recommend_hpc"] = True
                result["hpc_reason"] = (
                    f"No GPU detected. RL with {num_steps} steps and "
                    f"batch_size={batch_size} on CPU will be very slow "
                    f"(estimated {num_steps * 30 // 60} min+). "
                    f"Recommend HPC submission."
                )
            elif result["gpu_memory_mb"] < 8000 and num_steps > 500:
                result["recommend_hpc"] = True
                result["hpc_reason"] = (
                    f"Large job ({num_steps} steps). GPU has "
                    f"{result['gpu_memory_mb']} MB VRAM. "
                    f"Recommend HPC with A100 (40/80 GB) for production runs."
                )

    elif run_mode == "transfer_learning":
        if num_steps > 100:  # num_steps = epochs for TL
            if not result["gpu_available"]:
                result["recommend_hpc"] = True
                result["hpc_reason"] = (
                    f"No GPU detected. Transfer learning with {num_steps} epochs "
                    f"on CPU will be slow. Recommend HPC submission."
                )

    return result


# ---------------------------------------------------------------------------
# Seed SMILES preparation
# ---------------------------------------------------------------------------

def prepare_seeds(
    smiles_list: list[str],
    output_file: str,
) -> SeedReport:
    """Validate, canonicalize, and deduplicate seed SMILES.

    Args:
        smiles_list: List of SMILES strings.
        output_file: Path to write validated SMILES (one per line).

    Returns:
        SeedReport with validation results.
    """
    report = SeedReport(total_input=len(smiles_list))

    try:
        from rdkit import Chem
    except ImportError:
        # Fallback: write as-is without validation
        os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
        with open(output_file, "w") as f:
            for smi in smiles_list:
                f.write(smi.strip() + "\n")
        report.valid = len(smiles_list)
        report.output_file = output_file
        return report

    canonical = []
    seen = set()

    for smi in smiles_list:
        smi = smi.strip()
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            report.invalid += 1
            report.invalid_smiles.append(smi)
            continue
        can = Chem.MolToSmiles(mol)
        if can in seen:
            report.duplicates_removed += 1
            continue
        seen.add(can)
        canonical.append(can)

    report.valid = len(canonical)
    report.output_file = output_file

    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    with open(output_file, "w") as f:
        for smi in canonical:
            f.write(smi + "\n")

    return report


# ---------------------------------------------------------------------------
# Scoring profiles
# ---------------------------------------------------------------------------

def _component_block(
    name: str,
    weight: float = 1.0,
    transform: dict | None = None,
    params: dict | None = None,
) -> dict:
    """Build a single scoring component dictionary."""
    endpoint = {"name": name, "weight": weight}
    if transform:
        endpoint["transform"] = transform
    if params:
        endpoint["params"] = params
    return {name: {"endpoint": [endpoint]}}


def get_scoring_profile(profile_name: str) -> list[dict]:
    """Return a list of scoring component dicts for a named profile.

    Args:
        profile_name: One of drug-like, lead-like, fragment-like,
            kinase-inhibitor, anti-tb.

    Returns:
        List of component dictionaries for TOML generation.
    """
    profiles = {
        "drug-like": [
            _component_block("Qed", weight=1.0),
            _component_block(
                "MolecularWeight", weight=0.5,
                transform={"type": "double_sigmoid", "high": 500.0, "low": 200.0,
                            "coef_div": 500.0, "coef_si": 20.0, "coef_se": 20.0},
            ),
            _component_block(
                "SlogP", weight=0.5,
                transform={"type": "reverse_sigmoid", "high": 5.0, "low": 1.0, "k": 0.5},
            ),
            _component_block(
                "TPSA", weight=0.3,
                transform={"type": "double_sigmoid", "high": 140.0, "low": 20.0,
                            "coef_div": 140.0, "coef_si": 20.0, "coef_se": 20.0},
            ),
            _component_block(
                "NumRotBond", weight=0.2,
                transform={"type": "reverse_sigmoid", "high": 10.0, "low": 3.0, "k": 0.5},
            ),
            _component_block(
                "custom_alerts", weight=0.79,
                params={"smarts": PAINS_SMARTS},
            ),
        ],
        "lead-like": [
            _component_block("Qed", weight=1.0),
            _component_block(
                "MolecularWeight", weight=0.6,
                transform={"type": "double_sigmoid", "high": 350.0, "low": 200.0,
                            "coef_div": 350.0, "coef_si": 20.0, "coef_se": 20.0},
            ),
            _component_block(
                "SlogP", weight=0.5,
                transform={"type": "double_sigmoid", "high": 3.0, "low": 0.0,
                            "coef_div": 100.0, "coef_si": 20.0, "coef_se": 20.0},
            ),
            _component_block(
                "HBondDonors", weight=0.3,
                transform={"type": "reverse_sigmoid", "high": 3.0, "low": 0.0, "k": 0.5},
            ),
            _component_block(
                "HBondAcceptors", weight=0.3,
                transform={"type": "reverse_sigmoid", "high": 6.0, "low": 0.0, "k": 0.5},
            ),
            _component_block(
                "NumRotBond", weight=0.3,
                transform={"type": "reverse_sigmoid", "high": 5.0, "low": 0.0, "k": 0.5},
            ),
            _component_block(
                "custom_alerts", weight=0.79,
                params={"smarts": PAINS_SMARTS},
            ),
        ],
        "fragment-like": [
            _component_block(
                "MolecularWeight", weight=1.0,
                transform={"type": "double_sigmoid", "high": 250.0, "low": 100.0,
                            "coef_div": 250.0, "coef_si": 20.0, "coef_se": 20.0},
            ),
            _component_block(
                "SlogP", weight=0.5,
                transform={"type": "double_sigmoid", "high": 3.0, "low": -1.0,
                            "coef_div": 100.0, "coef_si": 20.0, "coef_se": 20.0},
            ),
            _component_block(
                "HBondDonors", weight=0.3,
                transform={"type": "reverse_sigmoid", "high": 3.0, "low": 0.0, "k": 0.5},
            ),
            _component_block(
                "HBondAcceptors", weight=0.3,
                transform={"type": "reverse_sigmoid", "high": 3.0, "low": 0.0, "k": 0.5},
            ),
            _component_block(
                "NumRotBond", weight=0.3,
                transform={"type": "reverse_sigmoid", "high": 3.0, "low": 0.0, "k": 0.5},
            ),
        ],
        "kinase-inhibitor": [
            _component_block("Qed", weight=1.0),
            _component_block(
                "MolecularWeight", weight=0.5,
                transform={"type": "double_sigmoid", "high": 550.0, "low": 300.0,
                            "coef_div": 500.0, "coef_si": 20.0, "coef_se": 20.0},
            ),
            _component_block(
                "SlogP", weight=0.5,
                transform={"type": "double_sigmoid", "high": 5.0, "low": 1.0,
                            "coef_div": 100.0, "coef_si": 20.0, "coef_se": 20.0},
            ),
            _component_block(
                "HBondAcceptors", weight=0.3,
                transform={"type": "double_sigmoid", "high": 8.0, "low": 3.0,
                            "coef_div": 100.0, "coef_si": 20.0, "coef_se": 20.0},
            ),
            _component_block("SAScore", weight=0.4,
                transform={"type": "reverse_sigmoid", "high": 5.0, "low": 1.0, "k": 0.5},
            ),
            _component_block(
                "custom_alerts", weight=0.79,
                params={"smarts": PAINS_SMARTS},
            ),
        ],
        "anti-tb": [
            _component_block("Qed", weight=1.0),
            _component_block(
                "MolecularWeight", weight=0.5,
                transform={"type": "double_sigmoid", "high": 600.0, "low": 250.0,
                            "coef_div": 500.0, "coef_si": 20.0, "coef_se": 20.0},
            ),
            _component_block(
                "SlogP", weight=0.5,
                transform={"type": "double_sigmoid", "high": 5.0, "low": 1.0,
                            "coef_div": 100.0, "coef_si": 20.0, "coef_se": 20.0},
            ),
            _component_block("SAScore", weight=0.6,
                transform={"type": "reverse_sigmoid", "high": 4.0, "low": 1.0, "k": 0.5},
            ),
            _component_block(
                "custom_alerts", weight=0.79,
                params={"smarts": PAINS_SMARTS},
            ),
        ],
    }

    if profile_name not in profiles:
        raise ValueError(
            f"Unknown scoring profile '{profile_name}'. "
            f"Available: {', '.join(profiles.keys())}"
        )
    return profiles[profile_name]


def parse_custom_component(spec: str) -> dict:
    """Parse a custom component spec string.

    Format: "ComponentName:key1=val1,key2=val2"

    Args:
        spec: Component specification string.

    Returns:
        Component dictionary for TOML generation.
    """
    if ":" not in spec:
        return _component_block(spec)

    name, params_str = spec.split(":", 1)
    kwargs = {}
    transform = {}

    for pair in params_str.split(","):
        k, v = pair.split("=", 1)
        k = k.strip()
        v = v.strip()

        # Try numeric conversion
        try:
            v = float(v)
            if v == int(v):
                v = int(v)
        except ValueError:
            pass

        if k == "weight":
            kwargs["weight"] = float(v)
        elif k == "transform":
            transform["type"] = v
        elif k in ("high", "low", "k", "coef_div", "coef_si", "coef_se"):
            transform[k] = float(v)
        elif k in ("smiles", "smarts"):
            # These are list params
            kwargs.setdefault("params", {})[k] = [v]
        elif k == "radius":
            kwargs.setdefault("params", {})[k] = [int(v)]
        else:
            kwargs.setdefault("params", {})[k] = v

    comp = _component_block(
        name,
        weight=kwargs.get("weight", 1.0),
        transform=transform if transform else None,
        params=kwargs.get("params"),
    )
    return comp


# ---------------------------------------------------------------------------
# TOML config generation
# ---------------------------------------------------------------------------

def _toml_path(p: str) -> str:
    """Convert a file path to TOML-safe forward slashes."""
    return p.replace("\\", "/")


def _toml_value(v: Any) -> str:
    """Convert a Python value to a TOML-compatible string."""
    if isinstance(v, bool):
        return "true" if v else "false"
    elif isinstance(v, str):
        # Escape backslashes for TOML string safety
        safe = v.replace("\\", "/")
        return f'"{safe}"'
    elif isinstance(v, (int, float)):
        return str(v)
    elif isinstance(v, list):
        items = ", ".join(_toml_value(i) for i in v)
        return f"[{items}]"
    elif isinstance(v, dict):
        items = ", ".join(f"{k} = {_toml_value(val)}" for k, val in v.items())
        return f"{{{items}}}"
    return str(v)


def _write_toml_section(lines: list[str], prefix: str, data: dict, indent: str = "") -> None:
    """Recursively write a TOML section."""
    for key, value in data.items():
        if isinstance(value, dict):
            lines.append(f"{indent}[{prefix}.{key}]" if prefix else f"{indent}[{key}]")
            _write_toml_section(lines, f"{prefix}.{key}" if prefix else key, value, indent)
        else:
            lines.append(f"{indent}{key} = {_toml_value(value)}")


def generate_staged_learning_toml(
    *,
    prior_file: str,
    agent_file: str | None = None,
    device: str = "cpu",
    batch_size: int = 128,
    num_steps: int = 300,
    max_score: float = 0.9,
    min_steps: int = 50,
    sigma: int = 128,
    learning_rate: float = 0.0001,
    scoring_components: list[dict],
    aggregation: str = "geometric_mean",
    diversity_filter_type: str = "IdenticalMurckoScaffold",
    bucket_size: int = 25,
    smiles_file: str | None = None,
    tb_logdir: str = "tb_logs",
    summary_csv_prefix: str = "staged_learning",
    chkpt_file: str = "agent.chkpt",
) -> str:
    """Generate a complete staged_learning TOML config string.

    Args:
        prior_file: Path to prior model.
        agent_file: Path to agent model (defaults to prior_file).
        device: Torch device string.
        batch_size: Batch size per RL step.
        num_steps: Maximum RL steps.
        max_score: Early termination score threshold.
        min_steps: Minimum steps before early termination check.
        sigma: DAP sigma parameter.
        learning_rate: Optimizer learning rate.
        scoring_components: List of component dicts from get_scoring_profile.
        aggregation: Score aggregation function.
        diversity_filter_type: DF type.
        bucket_size: DF bucket size.
        smiles_file: Optional seed SMILES file for non-Reinvent generators.
        tb_logdir: TensorBoard log directory.
        summary_csv_prefix: CSV output prefix.
        chkpt_file: Checkpoint file name.

    Returns:
        TOML config as a string.
    """
    if agent_file is None:
        agent_file = prior_file

    lines = [
        "# REINVENT4 Staged Learning Configuration",
        "# Generated by Compound Synthesis Skill",
        "",
        'run_type = "staged_learning"',
        f'device = "{device}"',
        f'tb_logdir = "{tb_logdir}"',
        "",
        "[parameters]",
        "",
        f'summary_csv_prefix = "{summary_csv_prefix}"',
        "use_checkpoint = false",
        "purge_memories = false",
        f'prior_file = "{_toml_path(prior_file)}"',
        f'agent_file = "{_toml_path(agent_file)}"',
    ]

    if smiles_file:
        lines.append(f'smiles_file = "{_toml_path(smiles_file)}"')

    lines.extend([
        f"batch_size = {batch_size}",
        "unique_sequences = true",
        "randomize_smiles = true",
        "",
        "[learning_strategy]",
        "",
        'type = "dap"',
        f"sigma = {sigma}",
        f"rate = {learning_rate}",
        "",
        "[diversity_filter]",
        "",
        f'type = "{diversity_filter_type}"',
        f"bucket_size = {bucket_size}",
        "minscore = 0.4",
        "minsimilarity = 0.4",
        "",
        "# --- Stage 1 ---",
        "[[stage]]",
        "",
        f'chkpt_file = "{_toml_path(chkpt_file)}"',
        'termination = "simple"',
        f"max_score = {max_score}",
        f"min_steps = {min_steps}",
        f"max_steps = {num_steps}",
        "",
        "[stage.scoring]",
        f'type = "{aggregation}"',
        "",
    ])

    # Write scoring components
    for comp in scoring_components:
        for comp_name, comp_data in comp.items():
            lines.append("[[stage.scoring.component]]")
            lines.append(f"[stage.scoring.component.{comp_name}]")
            lines.append("")

            for endpoint in comp_data.get("endpoint", []):
                lines.append(f"[[stage.scoring.component.{comp_name}.endpoint]]")
                lines.append(f'name = "{endpoint.get("name", comp_name)}"')
                lines.append(f'weight = {endpoint.get("weight", 1.0)}')

                if "params" in endpoint:
                    for pk, pv in endpoint["params"].items():
                        lines.append(f"params.{pk} = {_toml_value(pv)}")

                if "transform" in endpoint:
                    t = endpoint["transform"]
                    for tk, tv in t.items():
                        lines.append(f"transform.{tk} = {_toml_value(tv)}")

                lines.append("")

    return "\n".join(lines) + "\n"


def generate_transfer_learning_toml(
    *,
    input_model_file: str,
    output_model_file: str,
    smiles_file: str,
    device: str = "cpu",
    num_epochs: int = 50,
    batch_size: int = 64,
    sample_batch_size: int = 100,
    save_every_n_epochs: int = 5,
    num_refs: int = 0,
    tb_logdir: str = "tb_TL",
    validation_smiles_file: str | None = None,
) -> str:
    """Generate a transfer_learning TOML config string."""
    lines = [
        "# REINVENT4 Transfer Learning Configuration",
        "# Generated by Compound Synthesis Skill",
        "",
        'run_type = "transfer_learning"',
        f'device = "{device}"',
        f'tb_logdir = "{tb_logdir}"',
        "",
        "[parameters]",
        "",
        f"num_epochs = {num_epochs}",
        f"save_every_n_epochs = {save_every_n_epochs}",
        f"batch_size = {batch_size}",
        f"num_refs = {num_refs}",
        f"sample_batch_size = {sample_batch_size}",
        f'input_model_file = "{_toml_path(input_model_file)}"',
        f'smiles_file = "{_toml_path(smiles_file)}"',
        f'output_model_file = "{_toml_path(output_model_file)}"',
    ]

    if validation_smiles_file:
        lines.append(f'validation_smiles_file = "{_toml_path(validation_smiles_file)}"')

    lines.append("")
    return "\n".join(lines) + "\n"


def generate_sampling_toml(
    *,
    model_file: str,
    device: str = "cpu",
    num_smiles: int = 1000,
    output_file: str = "samples.csv",
    unique_molecules: bool = True,
    randomize_smiles: bool = True,
    smiles_file: str | None = None,
) -> str:
    """Generate a sampling TOML config string."""
    lines = [
        "# REINVENT4 Sampling Configuration",
        "# Generated by Compound Synthesis Skill",
        "",
        'run_type = "sampling"',
        f'device = "{device}"',
        "",
        "[parameters]",
        "",
        f'model_file = "{_toml_path(model_file)}"',
    ]

    if smiles_file:
        lines.append(f'smiles_file = "{_toml_path(smiles_file)}"')

    lines.extend([
        f'output_file = "{output_file}"',
        f"num_smiles = {num_smiles}",
        f"unique_molecules = {_toml_value(unique_molecules)}",
        f"randomize_smiles = {_toml_value(randomize_smiles)}",
        "",
    ])

    return "\n".join(lines) + "\n"


def generate_scoring_toml(
    *,
    smiles_file: str,
    scoring_components: list[dict],
    aggregation: str = "geometric_mean",
    output_csv: str = "scoring_results.csv",
) -> str:
    """Generate a scoring-only TOML config string."""
    lines = [
        "# REINVENT4 Scoring Configuration",
        "# Generated by Compound Synthesis Skill",
        "",
        'run_type = "scoring"',
        "",
        "[parameters]",
        f'smiles_file = "{_toml_path(smiles_file)}"',
        f'output_csv = "{_toml_path(output_csv)}"',
        "",
        "[scoring]",
        f'type = "{aggregation}"',
        "",
    ]

    for comp in scoring_components:
        for comp_name, comp_data in comp.items():
            lines.append("[[scoring.component]]")
            lines.append(f"[scoring.component.{comp_name}]")
            lines.append("")
            for endpoint in comp_data.get("endpoint", []):
                lines.append(f"[[scoring.component.{comp_name}.endpoint]]")
                lines.append(f'name = "{endpoint.get("name", comp_name)}"')
                lines.append(f'weight = {endpoint.get("weight", 1.0)}')
                if "params" in endpoint:
                    for pk, pv in endpoint["params"].items():
                        lines.append(f"params.{pk} = {_toml_value(pv)}")
                if "transform" in endpoint:
                    t = endpoint["transform"]
                    for tk, tv in t.items():
                        lines.append(f"transform.{tk} = {_toml_value(tv)}")
                lines.append("")

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# HPC script generation
# ---------------------------------------------------------------------------

def generate_hpc_script(
    *,
    config_file: str,
    gpu: str = "a100",
    ngpu: int = 1,
    mem: str = "64G",
    time: str = "24:00:00",
    conda_env: str = "reinvent4",
    job_name: str = "reinvent4",
    partition: str = "gpu",
) -> str:
    """Generate a SLURM submission script for REINVENT4."""
    return f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --gres=gpu:{gpu}:{ngpu}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --cpus-per-task=8
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# --- Environment setup ---
module purge
module load anaconda3 cuda/12.6

conda activate {conda_env}

# --- Run REINVENT4 ---
echo "Starting REINVENT4 at $(date)"
echo "Config: {config_file}"
echo "Device: cuda:0"

reinvent -l reinvent_${{SLURM_JOB_ID}}.log -d cuda:0 {config_file}

echo "Finished REINVENT4 at $(date)"
"""


# ---------------------------------------------------------------------------
# Run REINVENT
# ---------------------------------------------------------------------------

def run_reinvent(
    config_path: str,
    reinvent_dir: str,
    log_file: str | None = None,
    device: str | None = None,
) -> dict:
    """Execute REINVENT4 as a subprocess.

    Args:
        config_path: Path to the TOML config file.
        reinvent_dir: Path to the REINVENT4 repository.
        log_file: Path to write the log file.
        device: Override device (e.g. "cpu", "cuda:0").

    Returns:
        Dict with return_code, stdout_tail, stderr_tail, log_file.
    """
    env = os.environ.copy()
    rdir = str(Path(reinvent_dir).resolve())

    # Add REINVENT4 to PYTHONPATH
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = f"{rdir}{os.pathsep}{existing}" if existing else rdir

    # Build command
    reinvent_cli = shutil.which("reinvent")
    if not reinvent_cli:
        # Try running as module
        cmd = [sys.executable, "-m", "reinvent.Reinvent"]
    else:
        cmd = [reinvent_cli]

    if log_file:
        cmd.extend(["-l", log_file])
    if device:
        cmd.extend(["-d", device])
    cmd.append(str(Path(config_path).resolve()))

    result = subprocess.run(
        cmd,
        env=env,
        capture_output=True,
        text=True,
        cwd=os.path.dirname(os.path.abspath(config_path)) or ".",
    )

    stdout_lines = result.stdout.strip().split("\n") if result.stdout else []
    stderr_lines = result.stderr.strip().split("\n") if result.stderr else []

    return {
        "return_code": result.returncode,
        "stdout_tail": stdout_lines[-20:] if len(stdout_lines) > 20 else stdout_lines,
        "stderr_tail": stderr_lines[-20:] if len(stderr_lines) > 20 else stderr_lines,
        "log_file": log_file,
        "success": result.returncode == 0,
    }


# ---------------------------------------------------------------------------
# Result analysis
# ---------------------------------------------------------------------------

def analyze_results(
    csv_path: str,
    top_n: int = 20,
    sort_by: str = "total_score",
) -> dict:
    """Parse REINVENT4 CSV output and extract top compounds.

    Args:
        csv_path: Path to the REINVENT4 output CSV file.
        top_n: Number of top compounds to return.
        sort_by: Column name to sort by (descending).

    Returns:
        Dict with summary stats and top compounds.
    """
    if not os.path.exists(csv_path):
        return {"error": f"CSV file not found: {csv_path}"}

    compounds = []
    headers = []

    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames or []
        for row in reader:
            compounds.append(row)

    if not compounds:
        return {"error": "No data in CSV file", "headers": headers}

    # Parse numeric columns
    for comp in compounds:
        for key in comp:
            try:
                comp[key] = float(comp[key])
            except (ValueError, TypeError):
                pass

    # Sort by the specified column
    if sort_by in headers:
        compounds.sort(key=lambda x: float(x.get(sort_by, 0)), reverse=True)

    top = compounds[:top_n]

    # Compute stats
    scores = [float(c.get(sort_by, 0)) for c in compounds if c.get(sort_by) is not None]
    stats = {}
    if scores:
        stats = {
            "total_compounds": len(compounds),
            "mean_score": sum(scores) / len(scores),
            "max_score": max(scores),
            "min_score": min(scores),
            "median_score": sorted(scores)[len(scores) // 2],
        }

    # Extract unique SMILES
    smiles_col = None
    for h in headers:
        if "smiles" in h.lower():
            smiles_col = h
            break

    unique_smiles = set()
    if smiles_col:
        unique_smiles = {c[smiles_col] for c in compounds if isinstance(c.get(smiles_col), str)}
        stats["unique_molecules"] = len(unique_smiles)

    return {
        "csv_file": csv_path,
        "headers": headers,
        "sort_by": sort_by,
        "stats": stats,
        "top_compounds": top,
    }


# ---------------------------------------------------------------------------
# Prior model download
# ---------------------------------------------------------------------------

def download_priors(
    reinvent_dir: str,
    models: list[str] | None = None,
) -> dict:
    """Download prior models from Zenodo.

    Args:
        reinvent_dir: Path to the REINVENT4 repository.
        models: List of model names (e.g. ['reinvent', 'libinvent']).
            If None, downloads all standard models.

    Returns:
        Dict with download results.
    """
    priors_dir = Path(reinvent_dir) / "priors"
    priors_dir.mkdir(parents=True, exist_ok=True)

    if models is None:
        models = ["reinvent", "libinvent", "linkinvent"]

    # Map short names to filenames
    name_to_file = {
        "reinvent": "reinvent.prior",
        "libinvent": "libinvent.prior",
        "linkinvent": "linkinvent.prior",
        "mol2mol_scaffold": "mol2mol_scaffold.prior",
        "mol2mol_scaffold_generic": "mol2mol_scaffold_generic.prior",
        "mol2mol_high": "mol2mol_high_similarity.prior",
        "mol2mol_medium": "mol2mol_medium_similarity.prior",
        "mol2mol_mmp": "mol2mol_mmp.prior",
        "pepinvent": "pepinvent.prior",
    }

    results = {"priors_dir": str(priors_dir), "downloaded": [], "skipped": [], "errors": []}

    # Get Zenodo file list
    try:
        resp = urllib.request.urlopen(ZENODO_BASE, timeout=30)
        record = json.loads(resp.read())
        zenodo_files = {f["key"]: f["links"]["self"] for f in record.get("files", [])}
    except Exception as e:
        results["errors"].append(f"Failed to fetch Zenodo record: {e}")
        return results

    for model_name in models:
        filename = name_to_file.get(model_name, f"{model_name}.prior")
        target = priors_dir / filename

        if target.exists():
            results["skipped"].append({"model": model_name, "file": str(target), "reason": "already exists"})
            continue

        url = zenodo_files.get(filename)
        if not url:
            results["errors"].append(f"Model '{filename}' not found on Zenodo record {ZENODO_RECORD}")
            continue

        try:
            print(f"  Downloading {filename}...", end=" ", flush=True)
            urllib.request.urlretrieve(url, str(target))
            size_mb = target.stat().st_size / (1024 * 1024)
            print(f"{size_mb:.1f} MB")
            results["downloaded"].append({"model": model_name, "file": str(target), "size_mb": round(size_mb, 1)})
        except Exception as e:
            results["errors"].append(f"Failed to download {filename}: {e}")

    return results
