"""Execute REINVENT 4 and prepare its inputs.

Running a generative-chemistry job is long and expensive, so this module
checks everything it can *before* starting: that the config parses, that
every file the config references exists, and that the requested device is
available. The previous implementation invoked ``reinvent`` and left the
user to discover a missing prior file from a traceback several minutes in.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
import tomllib
from pathlib import Path
from typing import Any

from ..chem.smiles import CompoundRecord, dedupe, load_compounds
from ..errors import (
    InvalidInputError,
    MissingDependencyError,
    ResourceNotFoundError,
    UsageError,
)
from ..io_utils import require_file, write_text
from .priors import find_reinvent_dir

__all__ = ["prepare_seeds", "preflight", "run_reinvent"]


def prepare_seeds(records: list[CompoundRecord], output_path: str, *,
                  generator: str = "reinvent",
                  include_names: bool = False) -> dict[str, Any]:
    """Validate and de-duplicate seed structures, then write a ``.smi`` file.

    Every input is accounted for in the report: kept, rejected as invalid,
    or rejected as a duplicate, with the reason.

    Raises:
        InvalidInputError: If no structure survives validation.
    """
    kept, rejected = dedupe(records)
    if not kept:
        raise InvalidInputError(
            f"None of the {len(records)} input structure(s) were usable.",
            hint="Check the SMILES are valid; the report lists each failure.",
        )

    if include_names:
        body = "\n".join(f"{r.smiles}\t{r.name}" for r in kept)
    else:
        body = "\n".join(r.smiles for r in kept)
    path = write_text(body + "\n", output_path)

    warnings: list[str] = []
    if generator == "libinvent" and not any("*" in r.smiles for r in kept):
        warnings.append(
            "LibInvent expects scaffolds with attachment points marked '*', "
            "but none of the seeds contain one."
        )
    if generator == "linkinvent" and not any("|" in r.smiles for r in kept):
        warnings.append(
            "LinkInvent expects two warheads per line separated by '|', but "
            "no seed line contains one."
        )
    if generator == "transfer_learning" and len(kept) < 50:
        warnings.append(
            f"Only {len(kept)} molecules. Transfer learning normally needs "
            "several hundred to meaningfully shift a prior."
        )

    return {
        "output_file": str(path),
        "generator": generator,
        "input_count": len(records),
        "written": len(kept),
        "rejected": len(rejected),
        "rejected_detail": rejected,
        "warnings": warnings,
        "note": (
            "Structures were canonicalised with RDKit before de-duplication, "
            "so differently-written forms of the same molecule collapsed to "
            "one entry."
        ),
    }


def _referenced_files(config: dict[str, Any]) -> list[tuple[str, str]]:
    """Collect ``(key, path)`` pairs for every file a config references."""
    keys = (
        "prior_file", "agent_file", "smiles_file", "input_model_file",
        "model_file", "validation_smiles_file",
    )
    found: list[tuple[str, str]] = []

    def walk(node: Any, prefix: str = "") -> None:
        if isinstance(node, dict):
            for key, value in node.items():
                if key in keys and isinstance(value, str):
                    found.append((f"{prefix}{key}", value))
                else:
                    walk(value, f"{prefix}{key}.")
        elif isinstance(node, list):
            for index, item in enumerate(node):
                walk(item, f"{prefix}{index}.")

    walk(config)
    return found


def preflight(config_path: str | Path, *, device: str | None = None,
              reinvent_dir: str | None = None) -> dict[str, Any]:
    """Check a config and its environment before committing to a run.

    Raises:
        InvalidInputError: If the TOML cannot be parsed or has no run_type.
        ResourceNotFoundError: If a file the config references is missing.
    """
    path = require_file(config_path)
    try:
        config = tomllib.loads(path.read_text(encoding="utf-8"))
    except tomllib.TOMLDecodeError as exc:
        raise InvalidInputError(
            f"{path} is not valid TOML: {exc}",
            hint="Regenerate it with 'generate-config' rather than editing "
                 "by hand.",
        ) from exc

    run_type = config.get("run_type")
    if not run_type:
        raise InvalidInputError(
            f"{path} does not set run_type.",
            hint="A REINVENT config must declare run_type at the top level.",
        )

    missing: list[dict[str, str]] = []
    present: list[dict[str, str]] = []
    for key, reference in _referenced_files(config):
        candidate = Path(reference)
        if not candidate.is_absolute():
            # Resolve relative to the config, then to the working directory,
            # then to the REINVENT checkout - all three are conventional.
            bases = [path.parent, Path.cwd()]
            directory = find_reinvent_dir(reinvent_dir)
            if directory:
                bases.append(directory)
            for base in bases:
                if (base / reference).is_file():
                    candidate = (base / reference).resolve()
                    break
        if candidate.is_file():
            present.append({"key": key, "path": str(candidate),
                            "size_mb": f"{candidate.stat().st_size / 1e6:.1f}"})
        else:
            missing.append({"key": key, "reference": reference})

    requested_device = device or config.get("device", "cuda:0")
    device_ok = True
    device_note = ""
    if str(requested_device).startswith("cuda"):
        try:
            import torch
            if torch.cuda.is_available():
                device_note = (f"{torch.cuda.device_count()} CUDA device(s): "
                               f"{torch.cuda.get_device_name(0)}")
            else:
                device_ok = False
                device_note = ("No CUDA device is visible. Pass --device cpu, "
                               "or submit this to a GPU node.")
        except ImportError:
            device_ok = False
            device_note = "PyTorch is not installed in this interpreter."
    else:
        device_note = f"Using {requested_device}."

    executable = shutil.which("reinvent")
    report = {
        "config": str(path),
        "run_type": run_type,
        "device_requested": requested_device,
        "device_ok": device_ok,
        "device_note": device_note,
        "reinvent_executable": executable,
        "files_present": present,
        "files_missing": missing,
        "ready": bool(executable) and not missing and device_ok,
    }
    if missing:
        raise ResourceNotFoundError(
            "The config references "
            f"{len(missing)} file(s) that do not exist: "
            + ", ".join(f"{m['key']} = {m['reference']}" for m in missing),
            hint="Run 'check-setup' to see whether the prior models are "
                 "installed, and 'prepare-seeds' to create a SMILES file.",
        )
    return report


def run_reinvent(config_path: str | Path, *, log_path: str | None = None,
                 device: str | None = None, reinvent_dir: str | None = None,
                 timeout: float | None = None,
                 dry_run: bool = False) -> dict[str, Any]:
    """Run REINVENT 4 against a config.

    Args:
        config_path: The TOML config.
        log_path: Where REINVENT should write its log.
        device: Override the config's device.
        reinvent_dir: The REINVENT checkout, if not discoverable.
        timeout: Wall-clock limit in seconds. ``None`` means no limit.
        dry_run: Perform the preflight checks and return without running.

    Raises:
        MissingDependencyError: If the ``reinvent`` command is unavailable.
        ResourceNotFoundError: If a referenced file is missing.
    """
    checks = preflight(config_path, device=device, reinvent_dir=reinvent_dir)
    path = Path(config_path).resolve()

    executable = checks["reinvent_executable"]
    if not executable:
        raise MissingDependencyError(
            "the 'reinvent' command",
            install="pip install -e <REINVENT4 checkout>, then re-run "
                    "'check-setup' to confirm",
        )

    if dry_run:
        return {"dry_run": True, "would_run": [executable, str(path)],
                "preflight": checks}

    command = [executable]
    if log_path:
        command += ["--log-filename", str(log_path)]
    if device:
        command += ["--device", device]
    command.append(str(path))

    started = time.monotonic()
    # REINVENT resolves relative paths in the config against the working
    # directory, so run from the config's own directory for predictability.
    try:
        proc = subprocess.run(
            command, cwd=str(path.parent), capture_output=True, text=True,
            timeout=timeout, check=False,
        )
    except subprocess.TimeoutExpired:
        return {
            "status": "timeout",
            "command": command,
            "elapsed_seconds": round(time.monotonic() - started, 1),
            "message": f"REINVENT exceeded the {timeout:.0f}s limit. Partial "
                       "CSV output and the checkpoint file may still be "
                       "usable.",
        }
    elapsed = time.monotonic() - started

    produced = sorted(
        str(p) for p in path.parent.iterdir()
        if p.is_file() and p.stat().st_mtime >= started - 1
        and p.suffix in (".csv", ".chkpt", ".model", ".json", ".log")
    )

    result = {
        "status": "ok" if proc.returncode == 0 else "failed",
        "return_code": proc.returncode,
        "command": command,
        "working_directory": str(path.parent),
        "elapsed_seconds": round(elapsed, 1),
        "preflight": checks,
        "output_files": produced,
        "stdout_tail": (proc.stdout or "")[-4000:],
        "stderr_tail": (proc.stderr or "")[-4000:],
    }
    if proc.returncode != 0:
        result["message"] = (
            f"REINVENT exited with status {proc.returncode}. The stderr tail "
            "above usually names the cause; a missing prior or an "
            "unparseable scoring component are the common ones."
        )
    else:
        csvs = [f for f in produced if f.endswith(".csv")]
        result["message"] = (
            f"Completed in {elapsed / 60:.1f} min. "
            + (f"Analyse {Path(csvs[0]).name} with 'analyze-results'."
               if csvs else "No CSV output was produced; check the log.")
        )
    return result
