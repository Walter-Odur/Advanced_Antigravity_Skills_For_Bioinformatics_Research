"""AutoDock Vina docking.

Defects corrected here, all of which made the previous ``dock`` subcommand
unusable:

* **The command crashed on startup.** The parser defined ``--size`` but the
  implementation read ``args.box_size``, so any invocation using
  ``--center`` died with ``AttributeError: 'Namespace' object has no
  attribute 'box_size'``. The alternative input, ``--site``, pointed at a
  JSON file that no subcommand in any skill could produce.
* **A crash on the reporting line.** ``print(f"{score:.1f} kcal/mol")`` runs
  unconditionally, but ``score`` is ``None`` whenever no affinity could be
  parsed from Vina's output, raising ``TypeError: unsupported format string
  passed to NoneType.__format__`` after the docking work was already done.
* **Compound names used unsanitised as filenames.** Output PDBQT paths were
  built by string interpolation of the compound name, so a name containing
  ``/`` or ``:`` - routine in catalogue identifiers - wrote outside the
  output directory or failed outright on Windows.
* **Fragile score parsing.** Affinities were scraped with a regular
  expression that also matched Vina's banner and progress output.

Reference:
    Eberhardt J, Santos-Martins D, Tillack AF, Forli S.
    *J Chem Inf Model* 2021;61:3891-3898.
    Trott O, Olson AJ. *J Comput Chem* 2010;31:455-461.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from .._deps import ensure_vina
from ..chem.smiles import CompoundRecord, parse_smiles
from ..errors import MissingDependencyError, ResourceNotFoundError, UsageError
from ..io_utils import ensure_parent, require_file, safe_stem, write_text
from ..targets.pdbqt import ligand_to_pdbqt
from ..targets.site import BindingSite

__all__ = [
    "DockingResult",
    "find_vina",
    "prepare_ligand_3d",
    "dock_compounds",
    "parse_vina_log",
]

#: Vina's result table rows look like::
#:
#:        1       -9.428          0          0
#:        2       -9.106      1.892      2.771
#:
#: i.e. mode, affinity (kcal/mol), RMSD lower bound, RMSD upper bound.
#:
#: Note that the RMSD columns of mode 1 are printed as a bare ``0``, not
#: ``0.000``: a pattern that insists on a decimal point silently drops the
#: best pose, which is the one that matters most. The affinity may also be
#: an integer, and can be positive for a badly clashing pose.
_VINA_ROW = re.compile(
    r"^\s*(?P<mode>\d+)\s+"
    r"(?P<affinity>[-+]?\d+(?:\.\d+)?)\s+"
    r"(?P<rmsd_lb>\d+(?:\.\d+)?)\s+"
    r"(?P<rmsd_ub>\d+(?:\.\d+)?)\s*$"
)


@dataclass
class PoseScore:
    """One docked pose."""

    mode: int
    affinity_kcal_mol: float
    rmsd_lb: float
    rmsd_ub: float


@dataclass
class DockingResult:
    """Docking outcome for one compound."""

    name: str
    smiles: str
    status: str
    best_affinity_kcal_mol: float | None = None
    poses: list[PoseScore] = field(default_factory=list)
    output_pdbqt: str | None = None
    elapsed_seconds: float | None = None
    error: str = ""

    def as_dict(self) -> dict[str, Any]:
        data = asdict(self)
        if isinstance(data.get("elapsed_seconds"), float):
            data["elapsed_seconds"] = round(data["elapsed_seconds"], 2)
        return data


def find_vina() -> tuple[str, str]:
    """Locate a Vina implementation.

    Returns ``(kind, detail)`` where *kind* is ``"python"``, ``"cli"`` or
    ``"none"``.
    """
    try:
        import vina  # noqa: F401,PLC0415
        return "python", getattr(vina, "__version__", "unknown")
    except ImportError:
        pass
    for candidate in ("vina", "vina.exe", "vina_1.2.5", "AutoDockVina"):
        found = shutil.which(candidate)
        if found:
            return "cli", found
    # Nothing found -- try auto-installing the Python bindings.
    if ensure_vina():
        import vina  # noqa: F811,PLC0415
        return "python", getattr(vina, "__version__", "unknown")
    return "none", ""


def prepare_ligand_3d(smiles: str, *, name: str = "ligand",
                      max_iters: int = 400, seed: int = 0xF00D):
    """Build a minimised 3D conformer for docking.

    Vina needs 3D coordinates; a SMILES string has none. ETKDGv3 generates
    the conformer and MMFF (falling back to UFF) relaxes it.

    Args:
        smiles: Input SMILES.
        name: Used in error messages.
        max_iters: Force-field iteration cap.
        seed: Fixed so a docking run is reproducible. The previous code
            left embedding unseeded, so the same input gave different
            coordinates - and so different scores - on every run.

    Returns:
        An RDKit molecule with hydrogens and one conformer.

    Raises:
        UsageError: If the SMILES cannot be parsed or embedded.
    """
    from rdkit.Chem import AllChem

    mol = parse_smiles(smiles)
    if mol is None:
        raise UsageError(f"{name}: RDKit could not parse the SMILES {smiles!r}")

    from rdkit import Chem

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        # Retry with random coordinates, which rescues most macrocycles and
        # highly constrained ring systems.
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise UsageError(
                f"{name}: 3D embedding failed.",
                hint="Highly strained or macrocyclic structures sometimes "
                     "need manual conformer generation.",
            )
    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, maxIters=max_iters)
        else:
            AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
    except Exception:
        # An unminimised conformer still docks; geometry just starts rougher.
        pass
    return mol


def parse_vina_log(text: str) -> list[PoseScore]:
    """Extract the pose table from Vina's stdout.

    Only lines matching the exact four-column result row are accepted, so
    the banner, warnings and progress bar cannot be mistaken for scores.

    >>> poses = parse_vina_log("   1      -9.5      0.000      0.000\\n")
    >>> poses[0].affinity_kcal_mol
    -9.5
    """
    poses: list[PoseScore] = []
    for line in text.splitlines():
        match = _VINA_ROW.match(line)
        if not match:
            continue
        poses.append(PoseScore(
            mode=int(match.group("mode")),
            affinity_kcal_mol=float(match.group("affinity")),
            rmsd_lb=float(match.group("rmsd_lb")),
            rmsd_ub=float(match.group("rmsd_ub")),
        ))
    # Mode 1 must be the best; if the table was partially captured, sort.
    poses.sort(key=lambda p: p.mode)
    return poses


def _dock_one_cli(executable: str, receptor: Path, ligand_pdbqt: str,
                  site: BindingSite, out_pdbqt: Path, *, exhaustiveness: int,
                  n_poses: int, cpu: int, seed: int,
                  timeout: float) -> tuple[list[PoseScore], str]:
    """Dock one ligand with the Vina command-line binary."""
    ligand_file = out_pdbqt.with_name(out_pdbqt.stem + "_ligand.pdbqt")
    write_text(ligand_pdbqt, ligand_file)
    command = [
        executable,
        "--receptor", str(receptor),
        "--ligand", str(ligand_file),
        "--center_x", f"{site.center_x:.3f}",
        "--center_y", f"{site.center_y:.3f}",
        "--center_z", f"{site.center_z:.3f}",
        "--size_x", f"{site.size_x:.3f}",
        "--size_y", f"{site.size_y:.3f}",
        "--size_z", f"{site.size_z:.3f}",
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(n_poses),
        "--cpu", str(cpu),
        "--seed", str(seed),
        "--out", str(out_pdbqt),
    ]
    proc = subprocess.run(command, capture_output=True, text=True,
                          timeout=timeout, check=False)
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError(f"vina exited {proc.returncode}: {detail[:400]}")
    return parse_vina_log(proc.stdout), proc.stdout


def dock_compounds(records: list[CompoundRecord], receptor: str | Path,
                   site: BindingSite, *, output_dir: str | Path,
                   exhaustiveness: int = 32, n_poses: int = 9,
                   cpu: int = 0, seed: int = 42,
                   timeout: float = 600.0) -> dict[str, Any]:
    """Dock a set of compounds into a receptor.

    Args:
        records: Compounds to dock.
        receptor: Prepared receptor PDBQT.
        site: Search box, from :mod:`agskills.targets.site`.
        output_dir: Where pose files are written.
        exhaustiveness: Vina search effort. 8 is Vina's default; 32 is a
            reasonable screening value; raise it for a large box.
        n_poses: Poses to keep per compound.
        cpu: Vina worker threads. 0 lets Vina decide.
        seed: Fixed RNG seed, so a run is reproducible.
        timeout: Per-compound wall-clock limit in seconds.

    Raises:
        ResourceNotFoundError: If the receptor file is missing.
        MissingDependencyError: If neither Vina nor Meeko is available.
    """
    receptor_path = require_file(receptor)
    if receptor_path.suffix.lower() != ".pdbqt":
        raise UsageError(
            f"The receptor must be a PDBQT file, got {receptor_path.name}.",
            hint="Produce one with 'prepare-receptor', which converts a "
                 "cleaned PDB to PDBQT.",
        )
    if not records:
        raise UsageError("No compounds to dock.")

    kind, detail = find_vina()
    if kind == "none":
        raise MissingDependencyError(
            "AutoDock Vina",
            install=("conda install -c conda-forge vina, or download a "
                     "binary from "
                     "https://github.com/ccsb-scripps/AutoDock-Vina/releases"),
        )

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    results: list[DockingResult] = []
    used_stems: set[str] = set()

    vina_module = None
    if kind == "python":
        from vina import Vina  # noqa: PLC0415
        vina_module = Vina

    for index, record in enumerate(records, start=1):
        started = time.monotonic()
        # Names come from catalogues and CSVs, so they are sanitised and
        # de-duplicated before touching the filesystem.
        stem = safe_stem(record.name, fallback=f"compound_{index}")
        candidate = stem
        suffix = 2
        while candidate in used_stems:
            candidate = f"{stem}_{suffix}"
            suffix += 1
        used_stems.add(candidate)
        out_pdbqt = out_dir / f"{candidate}_docked.pdbqt"

        try:
            mol = prepare_ligand_3d(record.smiles, name=record.name, seed=seed)
            ligand_pdbqt = ligand_to_pdbqt(mol, name=record.name)
        except (UsageError, MissingDependencyError, RuntimeError) as exc:
            results.append(DockingResult(
                name=record.name, smiles=record.smiles, status="preparation_failed",
                error=str(exc), elapsed_seconds=time.monotonic() - started,
            ))
            continue

        try:
            if kind == "python":
                engine = vina_module(sf_name="vina", cpu=cpu, seed=seed,
                                     verbosity=0)
                engine.set_receptor(str(receptor_path))
                engine.set_ligand_from_string(ligand_pdbqt)
                engine.compute_vina_maps(center=list(site.center),
                                         box_size=list(site.size))
                engine.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
                energies = engine.energies(n_poses=n_poses)
                poses = [
                    PoseScore(mode=i + 1, affinity_kcal_mol=float(row[0]),
                              rmsd_lb=float(row[1]) if len(row) > 1 else 0.0,
                              rmsd_ub=float(row[2]) if len(row) > 2 else 0.0)
                    for i, row in enumerate(energies)
                ]
                engine.write_poses(str(out_pdbqt), n_poses=n_poses,
                                   overwrite=True)
            else:
                poses, _ = _dock_one_cli(
                    detail, receptor_path, ligand_pdbqt, site, out_pdbqt,
                    exhaustiveness=exhaustiveness, n_poses=n_poses, cpu=cpu or 0,
                    seed=seed, timeout=timeout,
                )
        except subprocess.TimeoutExpired:
            results.append(DockingResult(
                name=record.name, smiles=record.smiles, status="timeout",
                error=f"exceeded {timeout:.0f}s",
                elapsed_seconds=time.monotonic() - started,
            ))
            continue
        except Exception as exc:
            results.append(DockingResult(
                name=record.name, smiles=record.smiles, status="docking_failed",
                error=f"{type(exc).__name__}: {exc}",
                elapsed_seconds=time.monotonic() - started,
            ))
            continue

        if not poses:
            results.append(DockingResult(
                name=record.name, smiles=record.smiles, status="no_pose",
                error="Vina returned no scored pose",
                output_pdbqt=str(out_pdbqt) if out_pdbqt.exists() else None,
                elapsed_seconds=time.monotonic() - started,
            ))
            continue

        results.append(DockingResult(
            name=record.name, smiles=record.smiles, status="ok",
            best_affinity_kcal_mol=poses[0].affinity_kcal_mol,
            poses=poses,
            output_pdbqt=str(out_pdbqt) if out_pdbqt.exists() else None,
            elapsed_seconds=time.monotonic() - started,
        ))

    scored = [r for r in results if r.best_affinity_kcal_mol is not None]
    unscored = [r for r in results if r.best_affinity_kcal_mol is None]
    scored.sort(key=lambda r: r.best_affinity_kcal_mol)  # most negative first

    return {
        "engine": f"AutoDock Vina ({kind}: {detail})",
        "receptor": str(receptor_path),
        "binding_site": site.as_dict(),
        "parameters": {
            "exhaustiveness": exhaustiveness,
            "n_poses": n_poses,
            "seed": seed,
            "cpu": cpu or "auto",
        },
        "total_compounds": len(results),
        "docked": len(scored),
        "failed": len(unscored),
        "best_affinity_kcal_mol": (scored[0].best_affinity_kcal_mol
                                   if scored else None),
        "ranked_results": [r.as_dict() for r in scored + unscored],
        "interpretation": (
            "Affinity is in kcal/mol and more negative is stronger. Vina "
            "scores are a rough guide: differences under about 1 kcal/mol "
            "are not meaningful, and the function is a scoring, not a "
            "binding free energy. Treat the ranking as triage and confirm "
            "the top poses visually and by MD."
        ),
        "citation": "Eberhardt J et al. J Chem Inf Model 2021;61:3891-3898",
    }
