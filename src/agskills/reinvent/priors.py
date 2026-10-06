"""Locate a REINVENT 4 installation and its prior model files.

The skill documentation this replaces hardcoded an absolute path,
``E:\\ANTIGRAVITY_WORKSHOP\\REINVENT4``, in every example command. On anyone
else's machine every one of those commands fails. Discovery is now, in
order: an explicit argument, the ``REINVENT_DIR`` environment variable, the
installed ``reinvent`` package's location, then a search of likely
directories near the working tree.
"""

from __future__ import annotations

import importlib.util
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..errors import ResourceNotFoundError, UsageError

__all__ = [
    "GENERATORS",
    "PRIOR_FILES",
    "find_reinvent_dir",
    "find_priors_dir",
    "resolve_prior",
    "check_setup",
]

#: The five REINVENT generators, their prior files and what each is for.
GENERATORS: dict[str, dict[str, Any]] = {
    "reinvent": {
        "prior": "reinvent.prior",
        "architecture": "RNN",
        "needs_seeds": False,
        "use_for": "Unconstrained de novo design. No input structure needed.",
    },
    "libinvent": {
        "prior": "libinvent.prior",
        "architecture": "Transformer",
        "needs_seeds": True,
        "seed_format": "one scaffold per line with attachment points, "
                       "written as '*'",
        "use_for": "R-group decoration of a fixed scaffold.",
    },
    "linkinvent": {
        "prior": "linkinvent.prior",
        "architecture": "Transformer",
        "needs_seeds": True,
        "seed_format": "two warheads per line separated by '|'",
        "use_for": "Designing a linker between two fragments.",
    },
    "mol2mol": {
        "prior": "mol2mol_similarity.prior",
        "architecture": "Transformer",
        "needs_seeds": True,
        "seed_format": "one compound per line",
        "use_for": "Optimising or hopping from a known compound.",
        "variants": {
            "similarity": "mol2mol_similarity.prior",
            "medium_similarity": "mol2mol_medium_similarity.prior",
            "high_similarity": "mol2mol_high_similarity.prior",
            "mmp": "mol2mol_mmp.prior",
            "scaffold": "mol2mol_scaffold.prior",
            "scaffold_generic": "mol2mol_scaffold_generic.prior",
        },
    },
    "pepinvent": {
        "prior": "pepinvent.prior",
        "architecture": "Transformer",
        "needs_seeds": True,
        "seed_format": "one peptide per line with masked positions",
        "use_for": "Peptide design.",
    },
}

#: Every prior filename the package may look for.
PRIOR_FILES: tuple[str, ...] = tuple(sorted({
    spec["prior"] for spec in GENERATORS.values()
} | {
    name for spec in GENERATORS.values()
    for name in (spec.get("variants") or {}).values()
}))

_PRIOR_DOI = "https://doi.org/10.5281/zenodo.15641296"


def find_reinvent_dir(explicit: str | None = None) -> Path | None:
    """Find a REINVENT 4 source checkout.

    Resolution order:

    1. *explicit* (a ``--reinvent-dir`` argument).
    2. ``$REINVENT_DIR``.
    3. The directory containing an importable ``reinvent`` package.
    4. A ``REINVENT4`` directory beside the current working tree.
    """
    candidates: list[Path] = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    env = os.environ.get("REINVENT_DIR")
    if env:
        candidates.append(Path(env).expanduser())

    spec = importlib.util.find_spec("reinvent")
    if spec and spec.origin:
        # .../REINVENT4/reinvent/__init__.py -> .../REINVENT4
        candidates.append(Path(spec.origin).parent.parent)

    cwd = Path.cwd()
    for base in (cwd, *cwd.parents[:3]):
        candidates.append(base / "REINVENT4")
        candidates.append(base / "reinvent4")

    for candidate in candidates:
        try:
            if candidate.is_dir() and (candidate / "reinvent").is_dir():
                return candidate.resolve()
        except OSError:
            continue
    return None


def find_priors_dir(reinvent_dir: Path | None = None) -> Path | None:
    """Find the directory holding the ``.prior`` model files.

    ``$REINVENT_PRIOR_BASE`` wins, so priors can live on fast local storage
    separate from the source checkout.
    """
    env = os.environ.get("REINVENT_PRIOR_BASE")
    if env:
        candidate = Path(env).expanduser()
        if candidate.is_dir():
            return candidate.resolve()
    if reinvent_dir:
        for name in ("priors", "models"):
            candidate = reinvent_dir / name
            if candidate.is_dir():
                return candidate.resolve()
    return None


def resolve_prior(generator: str, prior: str | None = None, *,
                  priors_dir: Path | None = None,
                  require_exists: bool = False) -> tuple[str, Path | None]:
    """Resolve a generator and optional prior name to a prior file path.

    Args:
        generator: A key from :data:`GENERATORS`.
        prior: A prior filename, a Mol2Mol variant name such as
            ``scaffold_generic``, or ``None`` for the generator's default.
        priors_dir: Where to look for the file.
        require_exists: Raise if the file is absent.

    Returns:
        ``(filename, resolved_path_or_None)``.

    Raises:
        UsageError: On an unknown generator or Mol2Mol variant.
        ResourceNotFoundError: If *require_exists* and the file is missing.
    """
    key = (generator or "").strip().lower()
    if key not in GENERATORS:
        raise UsageError(
            f"Unknown generator {generator!r}.",
            hint="Choose one of: " + ", ".join(
                f"{name} ({spec['use_for']})" for name, spec in GENERATORS.items()
            ),
        )
    spec = GENERATORS[key]
    variants = spec.get("variants") or {}

    if prior is None:
        filename = spec["prior"]
    else:
        wanted = prior.strip()
        # Accept a bare variant name, with or without a leading dot, as well
        # as a full filename.
        stripped = wanted.lstrip(".")
        if stripped in variants:
            filename = variants[stripped]
        elif wanted.endswith(".prior") or "/" in wanted or "\\" in wanted:
            filename = wanted
        elif variants:
            raise UsageError(
                f"{prior!r} is not a Mol2Mol variant.",
                hint="Variants: " + ", ".join(sorted(variants)),
            )
        else:
            filename = stripped if stripped.endswith(".prior") else f"{stripped}.prior"

    path: Path | None = None
    as_given = Path(filename).expanduser()
    if as_given.is_file():
        path = as_given.resolve()
    elif priors_dir:
        candidate = priors_dir / Path(filename).name
        if candidate.is_file():
            path = candidate.resolve()

    if path is None and require_exists:
        raise ResourceNotFoundError(
            f"Prior model file not found: {filename}",
            hint=f"Download the REINVENT 4 priors from {_PRIOR_DOI} and put "
                 "them in <REINVENT4>/priors, or point REINVENT_PRIOR_BASE "
                 "at the directory holding them.",
        )
    return filename, path


def check_setup(reinvent_dir: str | None = None) -> dict[str, Any]:
    """Report on the REINVENT 4 installation, priors and compute device.

    Never raises for a missing component: the point of this command is to
    tell the user exactly what is absent and how to obtain it.
    """
    report: dict[str, Any] = {"ready": False, "problems": [], "advice": []}

    directory = find_reinvent_dir(reinvent_dir)
    report["reinvent_dir"] = str(directory) if directory else None
    if directory is None:
        report["problems"].append("No REINVENT 4 checkout found.")
        report["advice"].append(
            "Clone it and install: git clone "
            "https://github.com/MolecularAI/REINVENT4 && cd REINVENT4 && "
            "pip install -e . ; then pass --reinvent-dir or set REINVENT_DIR."
        )

    # Is the package importable, and is the CLI on PATH?
    spec = importlib.util.find_spec("reinvent")
    report["reinvent_importable"] = spec is not None
    executable = shutil.which("reinvent")
    report["reinvent_executable"] = executable
    if spec is None:
        report["problems"].append("The 'reinvent' package is not importable.")
        report["advice"].append(
            "Install it in this interpreter: pip install -e <REINVENT4>"
        )
    if executable is None:
        report["problems"].append("The 'reinvent' command is not on PATH.")

    version = None
    if spec is not None:
        try:
            from reinvent.version import __version__ as version  # type: ignore
        except Exception:
            try:
                result = subprocess.run(
                    [sys.executable, "-c",
                     "import reinvent.version as v; print(v.__version__)"],
                    capture_output=True, text=True, timeout=60, check=False,
                )
                version = (result.stdout or "").strip() or None
            except Exception:
                version = None
    report["reinvent_version"] = version

    # Priors.
    priors_dir = find_priors_dir(directory)
    report["priors_dir"] = str(priors_dir) if priors_dir else None
    found: dict[str, str] = {}
    missing: list[str] = []
    for filename in PRIOR_FILES:
        if priors_dir and (priors_dir / filename).is_file():
            size = (priors_dir / filename).stat().st_size
            found[filename] = f"{size / 1e6:.1f} MB"
        else:
            missing.append(filename)
    report["priors_found"] = found
    report["priors_missing"] = missing
    if missing:
        report["problems"].append(
            f"{len(missing)} of {len(PRIOR_FILES)} prior model files are missing."
        )
        report["advice"].append(
            f"Download the priors from {_PRIOR_DOI} into "
            f"{priors_dir or '<REINVENT4>/priors'}, or set "
            "REINVENT_PRIOR_BASE to where they live."
        )

    # Compute device.
    torch_info: dict[str, Any] = {"available": False}
    try:
        import torch  # noqa: PLC0415
        torch_info = {
            "available": True,
            "version": torch.__version__,
            "cuda_available": bool(torch.cuda.is_available()),
            "cuda_device_count": torch.cuda.device_count()
            if torch.cuda.is_available() else 0,
        }
        if torch.cuda.is_available():
            torch_info["devices"] = [
                torch.cuda.get_device_name(i)
                for i in range(torch.cuda.device_count())
            ]
            torch_info["recommended_device"] = "cuda:0"
        else:
            torch_info["recommended_device"] = "cpu"
            report["advice"].append(
                "No CUDA device is visible. Reinforcement learning on CPU is "
                "very slow; generate an HPC script and submit the run "
                "instead."
            )
    except ImportError:
        report["problems"].append("PyTorch is not installed.")
        report["advice"].append(
            "REINVENT needs PyTorch. Install the build matching your CUDA "
            "version from https://pytorch.org."
        )
    report["torch"] = torch_info
    report["recommended_device"] = torch_info.get("recommended_device", "cpu")

    try:
        import rdkit  # noqa: PLC0415
        report["rdkit_version"] = rdkit.__version__
    except ImportError:
        report["problems"].append("RDKit is not installed.")

    report["ready"] = not report["problems"]
    report["summary"] = (
        "REINVENT 4 is ready to run."
        if report["ready"] else
        f"{len(report['problems'])} problem(s) must be resolved before "
        "REINVENT can run. Config generation and seed preparation work "
        "regardless."
    )
    return report
