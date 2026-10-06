"""Automatic dependency installation.

When a required tool is missing at runtime, these helpers install it
automatically using whatever method works: pip, conda/mamba, system
package managers (apt, brew, winget), or git clone.  The user never
sees an install prompt -- the skill just works.

The resolution order for each tool:

1. **Check** -- is the tool already available?  If so, return immediately.
2. **pip** -- for pure Python packages (meeko, vina, biopython, pandas).
3. **conda / mamba** -- for packages that need compiled libraries
   (openbabel, gromacs).  Mamba is preferred over conda when available.
4. **System package manager** -- apt-get (Debian/Ubuntu), brew (macOS),
   or winget/choco (Windows) as a last resort.
5. **git clone + pip install** -- for projects like REINVENT 4 that
   require a source checkout.

Every function is idempotent: if the tool is already present, nothing
is installed.
"""

from __future__ import annotations

import importlib
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Sequence

__all__ = [
    "ensure_packages",
    "ensure_meeko",
    "ensure_vina",
    "ensure_pandas",
    "ensure_biopython",
    "ensure_openbabel",
    "ensure_gromacs",
    "ensure_reinvent",
]

# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def _install_disabled() -> bool:
    """True when auto-install has been suppressed via the environment.

    Set ``AGSKILLS_NO_AUTO_INSTALL=1`` in tests or CI to prevent the
    module from running pip/conda/apt during collection.
    """
    return os.environ.get("AGSKILLS_NO_AUTO_INSTALL", "") == "1"

def _run_quiet(command: list[str], *, timeout: float = 600) -> bool:
    """Run a command silently.  Return True if it exits 0."""
    try:
        proc = subprocess.run(
            command, capture_output=True, text=True,
            timeout=timeout, check=False,
        )
        return proc.returncode == 0
    except Exception:
        return False


def _pip_install(packages: Sequence[str], *, timeout: float = 300) -> bool:
    """Run pip install for the given specifiers."""
    return _run_quiet(
        [sys.executable, "-m", "pip", "install", "--quiet", *packages],
        timeout=timeout,
    )


def _find_conda() -> str | None:
    """Find mamba or conda, preferring mamba for speed."""
    for name in ("mamba", "micromamba", "conda"):
        found = shutil.which(name)
        if found:
            return found
    return None


def _conda_install(packages: Sequence[str], *,
                   channel: str = "conda-forge",
                   timeout: float = 600) -> bool:
    """Install packages via conda/mamba.  Returns True on success."""
    conda = _find_conda()
    if not conda:
        return False
    cmd = [conda, "install", "-y", "-c", channel, *packages]
    return _run_quiet(cmd, timeout=timeout)


def _system_install(packages_by_manager: dict[str, list[str]],
                    *, timeout: float = 300) -> bool:
    """Try the system package manager.

    Args:
        packages_by_manager: Mapping of manager name to package list, e.g.
            ``{"apt": ["openbabel"], "brew": ["open-babel"],
              "choco": ["openbabel"]}``.
    """
    system = platform.system().lower()

    if system == "linux":
        # Try apt-get (Debian/Ubuntu), then dnf (Fedora/RHEL)
        for mgr, flag in [("apt-get", "install"), ("dnf", "install"),
                          ("yum", "install"), ("pacman", "-S")]:
            exe = shutil.which(mgr)
            if exe and mgr.replace("-", "") in packages_by_manager:
                key = mgr.replace("-", "")
                pkgs = packages_by_manager.get(key, [])
                if not pkgs:
                    continue
                # apt needs update first
                if mgr == "apt-get":
                    _run_quiet(["sudo", exe, "update", "-qq"], timeout=120)
                cmd = ["sudo", exe, flag, "-y", *pkgs]
                if _run_quiet(cmd, timeout=timeout):
                    return True
            # Also try with the actual key name (apt-get -> apt)
            if mgr in packages_by_manager:
                pkgs = packages_by_manager[mgr]
                if not pkgs:
                    continue
                if mgr == "apt-get":
                    _run_quiet(["sudo", exe, "update", "-qq"], timeout=120)
                cmd = ["sudo", exe, flag, "-y", *pkgs]
                if _run_quiet(cmd, timeout=timeout):
                    return True

    elif system == "darwin":
        brew = shutil.which("brew")
        if brew and "brew" in packages_by_manager:
            pkgs = packages_by_manager["brew"]
            if pkgs and _run_quiet([brew, "install", *pkgs], timeout=timeout):
                return True

    elif system == "windows":
        # Try winget, then chocolatey
        for mgr, subcmd in [("winget", ["install", "--accept-source-agreements",
                                         "--accept-package-agreements"]),
                            ("choco", ["install", "-y"])]:
            exe = shutil.which(mgr)
            if exe and mgr in packages_by_manager:
                pkgs = packages_by_manager[mgr]
                if pkgs and _run_quiet([exe, *subcmd, *pkgs], timeout=timeout):
                    return True

    return False


def _git_clone_and_install(repo_url: str, *, target_dir: Path | None = None,
                           editable: bool = True,
                           timeout: float = 600) -> bool:
    """Clone a git repository and pip-install it."""
    git = shutil.which("git")
    if not git:
        return False

    if target_dir is None:
        # Install into a standard location
        base = Path.home() / ".agskills" / "packages"
        repo_name = repo_url.rstrip("/").rstrip(".git").rsplit("/", 1)[-1]
        target_dir = base / repo_name

    target_dir.parent.mkdir(parents=True, exist_ok=True)

    if target_dir.exists():
        # Already cloned -- pull latest
        _run_quiet([git, "-C", str(target_dir), "pull"], timeout=120)
    else:
        if not _run_quiet([git, "clone", repo_url, str(target_dir)],
                          timeout=timeout):
            return False

    # pip install from the checkout
    flag = "-e" if editable else "."
    install_arg = str(target_dir) if editable else str(target_dir)
    cmd = [sys.executable, "-m", "pip", "install", "--quiet"]
    if editable:
        cmd += ["-e", str(target_dir)]
    else:
        cmd += [str(target_dir)]
    return _run_quiet(cmd, timeout=300)


# ---------------------------------------------------------------------------
# Import-based checks (Python packages)
# ---------------------------------------------------------------------------

def ensure_packages(packages: Sequence[str], *, import_name: str) -> bool:
    """Import *import_name*; if that fails, pip-install *packages* and retry.

    Returns True if the import now succeeds (whether or not an install was
    needed).  False if the install failed.
    """
    try:
        importlib.import_module(import_name)
        return True
    except ImportError:
        pass
    if _install_disabled():
        return False
    if not _pip_install(packages):
        return False
    importlib.invalidate_caches()
    try:
        importlib.import_module(import_name)
        return True
    except ImportError:
        return False


# ---------------------------------------------------------------------------
# Per-dependency convenience functions
# ---------------------------------------------------------------------------

def ensure_meeko() -> bool:
    """Install meeko, scipy, and numpy if absent."""
    return ensure_packages(
        ["meeko>=0.5", "scipy>=1.10", "numpy>=1.24"],
        import_name="meeko",
    )


def ensure_vina() -> bool:
    """Install the vina Python bindings if absent."""
    return ensure_packages(["vina"], import_name="vina")


def ensure_pandas() -> bool:
    """Install pandas and numpy if absent."""
    return ensure_packages(
        ["pandas>=2.0", "numpy>=1.24"],
        import_name="pandas",
    )


def ensure_biopython() -> bool:
    """Install biopython if absent."""
    return ensure_packages(["biopython>=1.81"], import_name="Bio")


def ensure_openbabel() -> bool:
    """Install Open Babel (obabel binary) if absent.

    Resolution order:
    1. Already on PATH  -> done.
    2. conda/mamba       -> ``conda install -c conda-forge openbabel``
    3. System pkg mgr    -> apt/brew/choco
    """
    if shutil.which("obabel"):
        return True

    if _install_disabled():
        return False

    # Try conda first -- most reliable cross-platform route
    if _conda_install(["openbabel"]):
        # conda may put it in a place not on PATH yet; refresh
        if shutil.which("obabel"):
            return True

    # Try system package managers
    if _system_install({
        "apt-get": ["openbabel"],
        "apt": ["openbabel"],
        "dnf": ["openbabel"],
        "yum": ["openbabel"],
        "pacman": ["openbabel"],
        "brew": ["open-babel"],
        "choco": ["openbabel"],
        "winget": ["openbabel"],
    }):
        if shutil.which("obabel"):
            return True

    return False


def ensure_gromacs() -> bool:
    """Install GROMACS (gmx binary) if absent.

    Resolution order:
    1. Already on PATH   -> done.
    2. conda/mamba        -> ``conda install -c conda-forge gromacs``
    3. System pkg mgr     -> apt/brew
    """
    if shutil.which("gmx") or shutil.which("gmx_mpi"):
        return True

    if _install_disabled():
        return False

    # conda -- best cross-platform route
    if _conda_install(["gromacs"]):
        if shutil.which("gmx") or shutil.which("gmx_mpi"):
            return True

    # System package managers
    if _system_install({
        "apt-get": ["gromacs"],
        "apt": ["gromacs"],
        "dnf": ["gromacs"],
        "yum": ["gromacs"],
        "pacman": ["gromacs"],
        "brew": ["gromacs"],
    }):
        if shutil.which("gmx") or shutil.which("gmx_mpi"):
            return True

    return False


def ensure_reinvent() -> bool:
    """Install REINVENT 4 if absent.

    Resolution order:
    1. Already importable -> done.
    2. git clone + pip install -e .
    """
    try:
        importlib.import_module("reinvent")
        return True
    except ImportError:
        pass

    # Also check if the CLI is on PATH
    if shutil.which("reinvent"):
        return True

    if _install_disabled():
        return False

    return _git_clone_and_install(
        "https://github.com/MolecularAI/REINVENT4.git",
    )
