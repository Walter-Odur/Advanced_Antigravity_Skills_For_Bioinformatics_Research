"""Filesystem helpers used by every command.

Three properties matter here and are not negotiable, because the original
scripts got each of them wrong at least once:

1. **Always UTF-8.** ``Path.write_text`` without an encoding uses the
   locale codepage, which is cp1252 on a default Windows install. Reports
   containing an en dash then die with ``UnicodeEncodeError``.
2. **Always LF for generated scripts.** A SLURM script or ``.mdp`` file
   written with CRLF fails on a Linux cluster (``bash\\r: bad interpreter``).
3. **Atomic writes.** A crash halfway through serialisation must not leave
   a truncated JSON file that a downstream step then half-parses.
"""

from __future__ import annotations

import json
import os
import re
import tempfile
from pathlib import Path
from typing import Any

from .errors import ResourceNotFoundError

__all__ = [
    "write_json",
    "write_text",
    "read_json",
    "read_lines",
    "ensure_parent",
    "safe_stem",
    "require_file",
]

# Characters that are illegal in a Windows filename, plus path separators.
# Compound names arrive from CSV files and SMILES catalogues, so they can
# contain anything at all.
_UNSAFE_FILENAME_CHARS = re.compile(r'[<>:"/\\|?*\x00-\x1f]')
_WINDOWS_RESERVED = {
    "CON", "PRN", "AUX", "NUL",
    *(f"COM{i}" for i in range(1, 10)),
    *(f"LPT{i}" for i in range(1, 10)),
}


def ensure_parent(path: str | os.PathLike[str]) -> Path:
    """Create the parent directory of *path* and return *path* as a ``Path``."""
    p = Path(path)
    parent = p.parent
    if str(parent) not in ("", "."):
        parent.mkdir(parents=True, exist_ok=True)
    return p


def _atomic_write(path: Path, data: str, newline: str) -> None:
    """Write *data* to *path* atomically, replacing any existing file."""
    parent = path.parent if str(path.parent) not in ("", ".") else Path(".")
    fd, tmp_name = tempfile.mkstemp(dir=parent, prefix=f".{path.name}.", suffix=".tmp")
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline=newline) as fh:
            fh.write(data)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def write_json(data: Any, path: str | os.PathLike[str], *, indent: int = 2) -> Path:
    """Serialise *data* to *path* as UTF-8 JSON, atomically.

    ``default=str`` is applied so that values the science layer may hand us
    (``Path``, ``numpy`` scalars, ``Decimal``) serialise instead of raising.
    """
    p = ensure_parent(path)
    payload = json.dumps(data, indent=indent, default=str, ensure_ascii=False)
    _atomic_write(p, payload + "\n", newline="\n")
    return p


def write_text(
    text: str,
    path: str | os.PathLike[str],
    *,
    executable: bool = False,
) -> Path:
    """Write *text* to *path* as UTF-8 with LF line endings, atomically.

    Every generated artefact (SLURM scripts, ``.mdp`` files, ``.smi`` files,
    TOML configs, shell scripts) goes through here, so all of them are
    cluster-safe regardless of the platform that produced them.

    Args:
        text: File contents.
        path: Destination path.
        executable: Set the owner/group/other execute bits. A no-op on
            Windows, where execute permission is not a file mode bit.
    """
    p = ensure_parent(path)
    # Normalise any CRLF the caller's template may have introduced.
    normalised = text.replace("\r\n", "\n").replace("\r", "\n")
    _atomic_write(p, normalised, newline="\n")
    if executable:
        try:
            mode = p.stat().st_mode
            p.chmod(mode | 0o111)
        except OSError:
            # Windows, or a filesystem without permission bits.
            pass
    return p


def read_json(path: str | os.PathLike[str]) -> Any:
    """Read UTF-8 JSON from *path*."""
    p = require_file(path)
    return json.loads(p.read_text(encoding="utf-8"))


def read_lines(path: str | os.PathLike[str]) -> list[str]:
    """Read *path* and return its non-empty, non-comment lines, stripped.

    Lines beginning with ``#`` are treated as comments. Used for ``.smi``
    and plain SMILES/ID list files.
    """
    p = require_file(path)
    out: list[str] = []
    for raw in p.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            out.append(line)
    return out


def require_file(path: str | os.PathLike[str]) -> Path:
    """Return *path* as a ``Path``, raising if it is not an existing file."""
    p = Path(path)
    if not p.exists():
        raise ResourceNotFoundError(
            f"File not found: {p}",
            hint="Check the path, or run the step that produces it first.",
        )
    if not p.is_file():
        raise ResourceNotFoundError(f"Not a regular file: {p}")
    return p


def safe_stem(name: str, *, fallback: str = "compound", max_len: int = 64) -> str:
    """Turn an arbitrary compound name into a safe filename stem.

    Docking writes one PDBQT per compound, named after the compound. Names
    coming from ChEMBL, COCONUT or a user CSV routinely contain ``/``,
    ``:``, quotes and spaces, which would either escape the output directory
    or fail outright on Windows.

    >>> safe_stem("ZINC00001/trans-isomer")
    'ZINC00001_trans-isomer'
    >>> safe_stem("  ")
    'compound'
    >>> safe_stem("NUL")
    'NUL_'
    """
    cleaned = _UNSAFE_FILENAME_CHARS.sub("_", name).strip().strip(".")
    cleaned = re.sub(r"\s+", "_", cleaned)
    if not cleaned:
        return fallback
    if cleaned.upper() in _WINDOWS_RESERVED:
        cleaned += "_"
    return cleaned[:max_len]
