"""Automatic dependency installation.

When a pip-installable package is missing at runtime, these helpers install
it silently and return ``True`` on success.  This removes the single most
common friction point when a skill is used on a fresh machine: the user
had the core package installed but not the optional extras, and the skill
just stopped with an install hint instead of doing the work.

Only **pure pip packages** are handled.  System-level tools (GROMACS,
Open Babel via conda, REINVENT 4) cannot be pip-installed and are left
to the existing error-and-hint path.

Every function is idempotent: if the package is already present, the
import succeeds immediately and ``pip`` is never invoked.
"""

from __future__ import annotations

import importlib
import subprocess
import sys
from typing import Sequence

__all__ = [
    "ensure_packages",
    "ensure_meeko",
    "ensure_vina",
    "ensure_pandas",
    "ensure_biopython",
]


def ensure_packages(packages: Sequence[str], *, import_name: str) -> bool:
    """Import *import_name*; if that fails, pip-install *packages* and retry.

    Args:
        packages: Pip specifiers, e.g. ``["meeko>=0.5", "scipy>=1.10"]``.
        import_name: The top-level module to ``import`` (e.g. ``"meeko"``).

    Returns:
        ``True`` if the import now succeeds (whether or not an install was
        needed).  ``False`` if the install failed.
    """
    try:
        importlib.import_module(import_name)
        return True
    except ImportError:
        pass
    try:
        proc = subprocess.run(
            [sys.executable, "-m", "pip", "install", "--quiet", *packages],
            capture_output=True, text=True, timeout=300, check=False,
        )
        if proc.returncode != 0:
            return False
        # Force re-evaluation of the module search path so the freshly
        # installed package is visible without restarting the interpreter.
        importlib.invalidate_caches()
        importlib.import_module(import_name)
        return True
    except Exception:
        return False


# ----- Convenience wrappers for each dependency group ---------------------

def ensure_meeko() -> bool:
    """Install ``meeko``, ``scipy`` and ``numpy`` if absent."""
    return ensure_packages(
        ["meeko>=0.5", "scipy>=1.10", "numpy>=1.24"],
        import_name="meeko",
    )


def ensure_vina() -> bool:
    """Install the ``vina`` Python bindings if absent."""
    return ensure_packages(["vina"], import_name="vina")


def ensure_pandas() -> bool:
    """Install ``pandas`` and ``numpy`` if absent."""
    return ensure_packages(
        ["pandas>=2.0", "numpy>=1.24"],
        import_name="pandas",
    )


def ensure_biopython() -> bool:
    """Install ``biopython`` if absent."""
    return ensure_packages(["biopython>=1.81"], import_name="Bio")
