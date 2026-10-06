"""PDBQT conversion for AutoDock Vina.

Two things the code this replaces got wrong.

**The Meeko branch never converted anything.** It imported
``PDBQTWriterLegacy`` and then, on success, set a note reading "Open Babel
not available. PDBQT conversion requires obabel", returning no PDBQT. The
import result was never used.

**Real crystal structures were not handled.** Meeko matches every residue
against a chemical template, and residues with missing side-chain atoms,
alternate locations or unusual connectivity fail that match. This is not an
edge case: PDB 6HEZ (*M. tuberculosis* DprE1) has 16 such residues in
chain A alone, and strict conversion simply fails on it. Rather than
reporting failure, this module retries in Meeko's permissive mode and
**reports exactly which residues were compromised**, including whether any
of them sit inside the docking box - which is the part that would actually
invalidate a result.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .._deps import ensure_meeko
from ..errors import MissingDependencyError
from ..io_utils import ensure_parent, require_file, write_text

__all__ = ["ConversionResult", "receptor_to_pdbqt", "ligand_to_pdbqt",
           "available_converters", "find_executable"]

_INSTALL_HINT = (
    "Install one of: 'pip install meeko scipy gemmi' (Python, "
    "cross-platform), or 'conda install -c conda-forge openbabel' (provides "
    "the obabel binary)."
)

#: Meeko reports unmatched residues as a Python-style list of "CHAIN:RESSEQ".
_BAD_RES_RE = re.compile(r"Template matching failed for:\s*\[([^\]]*)\]")
_RES_TOKEN_RE = re.compile(r"'([^']+)'")


# _auto_install_meeko removed -- using centralized _deps.ensure_meeko()


@dataclass
class ConversionResult:
    """Outcome of a receptor or ligand PDBQT conversion."""

    ok: bool
    output_path: str | None
    converter: str
    message: str = ""
    strict: bool = True
    compromised_residues: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        from dataclasses import asdict
        return asdict(self)


def find_executable(name: str) -> str | None:
    """Locate an executable, including inside the running interpreter's env.

    ``shutil.which`` consults ``PATH`` only. A console script installed into
    a virtual environment lives next to ``sys.executable`` (in ``Scripts``
    on Windows, ``bin`` on POSIX) and is frequently absent from ``PATH`` -
    for instance when the interpreter was invoked by its full path, which
    is how most tooling runs it.
    """
    found = shutil.which(name)
    if found:
        return found
    candidates = [Path(sys.executable).parent]
    base = Path(sys.prefix)
    candidates += [base / "Scripts", base / "bin"]
    suffixes = (".exe", ".bat", ".cmd", "") if os.name == "nt" else ("",)
    for directory in candidates:
        for suffix in suffixes:
            candidate = directory / f"{name}{suffix}"
            if candidate.is_file():
                return str(candidate)
    return None


def available_converters() -> list[str]:
    """Which PDBQT converters this machine can actually use."""
    found: list[str] = []
    if find_executable("mk_prepare_receptor"):
        found.append("meeko")
    else:
        try:
            import meeko  # noqa: F401,PLC0415
            found.append("meeko-python")
        except Exception:
            pass
    if find_executable("obabel"):
        found.append("obabel")
    if find_executable("prepare_receptor") or \
            find_executable("prepare_receptor4.py"):
        found.append("adfr")
    return found


def _parse_bad_residues(output: str) -> list[str]:
    """Extract the residues Meeko could not template-match."""
    residues: list[str] = []
    for match in _BAD_RES_RE.finditer(output):
        residues.extend(_RES_TOKEN_RE.findall(match.group(1)))
    # Preserve order while removing duplicates.
    return list(dict.fromkeys(residues))


def _residues_in_box(residues: list[str], pdb_path: Path,
                     site: Any) -> list[str]:
    """Which of *residues* have an atom inside the docking box."""
    if site is None or not residues:
        return []
    from .site import parse_atoms

    wanted = set()
    for token in residues:
        chain, _, number = token.partition(":")
        try:
            wanted.add((chain.strip(), int(number)))
        except ValueError:
            continue
    if not wanted:
        return []

    half = (site.size_x / 2.0, site.size_y / 2.0, site.size_z / 2.0)
    centre = (site.center_x, site.center_y, site.center_z)
    inside: set[str] = set()
    for atom in parse_atoms(pdb_path.read_text(encoding="utf-8",
                                               errors="replace")):
        key = (atom.chain, atom.resseq)
        if key not in wanted:
            continue
        if (abs(atom.x - centre[0]) <= half[0]
                and abs(atom.y - centre[1]) <= half[1]
                and abs(atom.z - centre[2]) <= half[2]):
            inside.add(f"{atom.chain}:{atom.resseq}")
    return sorted(inside)


def _run(command: list[str], timeout: float = 900.0) -> subprocess.CompletedProcess:
    return subprocess.run(command, capture_output=True, text=True,
                          timeout=timeout, check=False)


def receptor_to_pdbqt(pdb_path: str | Path,
                      output_path: str | Path | None = None,
                      *, site: Any = None,
                      allow_incomplete_residues: bool = True,
                      default_altloc: str = "A") -> ConversionResult:
    """Convert a cleaned receptor PDB to PDBQT.

    Tries Meeko's receptor preparation first, strictly, then permissively if
    strict mode rejects residues; then Open Babel; then the ADFR Suite.

    Args:
        pdb_path: A cleaned receptor PDB.
        output_path: Destination ``.pdbqt``. Defaults to *pdb_path* with the
            suffix replaced.
        site: An optional :class:`~agskills.targets.site.BindingSite`. When
            given, any residue that had to be approximated is checked
            against the docking box, because a defect inside the box is the
            one that matters.
        allow_incomplete_residues: Retry permissively when strict template
            matching fails. With this off, an incomplete structure is
            reported as a failure.
        default_altloc: Alternate location to take when a residue has
            several.

    Returns ``ok=False`` with actionable advice rather than raising, since a
    caller may legitimately want the cleaned PDB even when no converter is
    installed.
    """
    source = require_file(pdb_path)
    destination = (Path(output_path) if output_path
                   else source.with_suffix(".pdbqt"))
    ensure_parent(destination)
    warnings: list[str] = []

    # --- Meeko ------------------------------------------------------------
    meeko_exe = find_executable("mk_prepare_receptor")
    if meeko_exe:
        stem = destination.with_suffix("")
        base_command = [meeko_exe, "--read_pdb", str(source),
                        "-o", str(stem), "-p",
                        "--default_altloc", default_altloc]
        try:
            proc = _run(base_command)
            combined = (proc.stdout or "") + (proc.stderr or "")
            if destination.is_file() and destination.stat().st_size > 0:
                return ConversionResult(
                    True, str(destination), "meeko",
                    "Receptor PDBQT written by mk_prepare_receptor "
                    "(all residues matched their templates).",
                    strict=True, warnings=warnings,
                )

            bad = _parse_bad_residues(combined)
            if bad and allow_incomplete_residues:
                # Meeko's own recommended route for batch processing. It
                # pads the unmatched residues from the template rather than
                # deleting them, so the chain stays intact.
                proc = _run(base_command + ["-a"])
                combined2 = (proc.stdout or "") + (proc.stderr or "")
                if destination.is_file() and destination.stat().st_size > 0:
                    in_box = _residues_in_box(bad, source, site)
                    warnings.append(
                        f"{len(bad)} residue(s) did not match a chemical "
                        "template, most often because side-chain atoms are "
                        "missing from the crystal structure. Meeko padded "
                        f"them from its templates: {', '.join(bad[:12])}"
                        + (f" (+{len(bad) - 12} more)" if len(bad) > 12 else "")
                    )
                    if in_box:
                        warnings.append(
                            "IMPORTANT: "
                            f"{len(in_box)} of those residue(s) lie inside "
                            f"the docking box ({', '.join(in_box)}). Their "
                            "geometry is approximated, so affinities for "
                            "this site are less reliable. Consider "
                            "rebuilding those side chains, or choosing a "
                            "structure with better occupancy."
                        )
                    elif site is not None:
                        warnings.append(
                            "None of the affected residues lie inside the "
                            "docking box, so the site itself is unaffected."
                        )
                    return ConversionResult(
                        True, str(destination), "meeko",
                        "Receptor PDBQT written by mk_prepare_receptor in "
                        "permissive mode.",
                        strict=False, compromised_residues=bad,
                        warnings=warnings,
                    )
                warnings.append(
                    "Permissive mode also failed: "
                    f"{combined2.strip()[-400:]}"
                )
            elif bad:
                warnings.append(
                    f"{len(bad)} residue(s) failed template matching and "
                    "permissive mode was disabled: "
                    f"{', '.join(bad[:12])}"
                )
            else:
                warnings.append(
                    f"mk_prepare_receptor produced no output: "
                    f"{combined.strip()[-400:]}"
                )
        except subprocess.TimeoutExpired:
            warnings.append("mk_prepare_receptor timed out.")
        except Exception as exc:
            warnings.append(f"mk_prepare_receptor failed: {exc}")

    # --- Open Babel -------------------------------------------------------
    obabel = find_executable("obabel")
    if obabel:
        try:
            proc = _run([obabel, str(source), "-O", str(destination),
                         "-xr", "-p", "7.4"])
            if destination.is_file() and destination.stat().st_size > 0:
                warnings.append(
                    "Open Babel assigns Gasteiger charges and does not "
                    "validate residue chemistry, so it is more permissive "
                    "than Meeko but less careful."
                )
                return ConversionResult(
                    True, str(destination), "obabel",
                    "Receptor PDBQT written by Open Babel (rigid, pH 7.4).",
                    strict=False, warnings=warnings,
                )
            warnings.append(
                f"obabel exited {proc.returncode}: "
                f"{(proc.stderr or '').strip()[:300]}"
            )
        except Exception as exc:
            warnings.append(f"obabel failed: {exc}")

    # --- ADFR Suite -------------------------------------------------------
    adfr = (find_executable("prepare_receptor")
            or find_executable("prepare_receptor4.py"))
    if adfr:
        try:
            proc = _run([adfr, "-r", str(source), "-o", str(destination)])
            if destination.is_file() and destination.stat().st_size > 0:
                return ConversionResult(
                    True, str(destination), "adfr",
                    "Receptor PDBQT written by the ADFR Suite.",
                    strict=False, warnings=warnings,
                )
            warnings.append(f"prepare_receptor exited {proc.returncode}")
        except Exception as exc:
            warnings.append(f"prepare_receptor failed: {exc}")

    # --- Auto-install Meeko if nothing was found ---------------------------
    if not meeko_exe and not obabel and not adfr:
        installed = ensure_meeko()
        if installed:
            meeko_exe = find_executable("mk_prepare_receptor")
            if meeko_exe:
                stem = destination.with_suffix("")
                base_command = [meeko_exe, "--read_pdb", str(source),
                                "-o", str(stem), "-p",
                                "--default_altloc", default_altloc]
                try:
                    proc = _run(base_command)
                    combined = (proc.stdout or "") + (proc.stderr or "")
                    if destination.is_file() and destination.stat().st_size > 0:
                        warnings.append(
                            "meeko was auto-installed to provide PDBQT "
                            "conversion."
                        )
                        return ConversionResult(
                            True, str(destination), "meeko",
                            "Receptor PDBQT written by mk_prepare_receptor "
                            "(auto-installed; all residues matched).",
                            strict=True, warnings=warnings,
                        )
                    bad = _parse_bad_residues(combined)
                    if bad and allow_incomplete_residues:
                        proc = _run(base_command + ["-a"])
                        if destination.is_file() and \
                                destination.stat().st_size > 0:
                            in_box = _residues_in_box(bad, source, site)
                            warnings.append(
                                "meeko was auto-installed to provide PDBQT "
                                "conversion."
                            )
                            warnings.append(
                                f"{len(bad)} residue(s) did not match a "
                                "chemical template. Meeko padded them: "
                                f"{', '.join(bad[:12])}"
                                + (f" (+{len(bad) - 12} more)"
                                   if len(bad) > 12 else "")
                            )
                            if in_box:
                                warnings.append(
                                    "IMPORTANT: "
                                    f"{len(in_box)} of those residue(s) lie "
                                    "inside the docking box "
                                    f"({', '.join(in_box)})."
                                )
                            return ConversionResult(
                                True, str(destination), "meeko",
                                "Receptor PDBQT written by mk_prepare_receptor"
                                " (auto-installed; permissive mode).",
                                strict=False, compromised_residues=bad,
                                warnings=warnings,
                            )
                except Exception as exc:
                    warnings.append(
                        f"mk_prepare_receptor failed after auto-install: {exc}"
                    )

    return ConversionResult(
        ok=False, output_path=None, converter="none",
        message=(
            "No PDBQT converter succeeded, so the receptor remains a cleaned "
            f"PDB at {source}. " + _INSTALL_HINT
        ),
        warnings=warnings,
    )


def ligand_to_pdbqt(mol, *, name: str = "ligand") -> str:
    """Convert an embedded, hydrogen-added RDKit molecule to a PDBQT string.

    Args:
        mol: An RDKit molecule that already has 3D coordinates and explicit
            hydrogens. Embedding is the caller's job, so that a failure can
            be reported per compound rather than aborting a batch.

    Raises:
        MissingDependencyError: If Meeko is unavailable.
        RuntimeError: If Meeko rejects the molecule.
    """
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
    except ImportError:
        if not ensure_meeko():
            raise MissingDependencyError(
                "meeko (required to write ligand PDBQT files)",
                install="pip install meeko scipy gemmi",
            )
        from meeko import MoleculePreparation, PDBQTWriterLegacy

    preparator = MoleculePreparation()
    setups = preparator.prepare(mol)
    if not setups:
        raise RuntimeError(f"Meeko produced no setup for {name}")
    pdbqt, ok, error = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise RuntimeError(f"Meeko could not write PDBQT for {name}: {error}")
    return pdbqt
