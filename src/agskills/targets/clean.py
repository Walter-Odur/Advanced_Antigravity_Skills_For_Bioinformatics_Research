"""Receptor cleaning: chain selection, solvent and heterogen removal.

The code this replaces documented that it "removes water, ligands and
non-standard residues", but its residue filter only tested
``residue.id[0] == "W"``. BioPython flags a water as ``"W"`` and any other
heterogen as ``"H_<resname>"``, so co-crystallised ligands, buffer
components, cryoprotectants and ions were all retained. Docking into a
receptor that still contains its original ligand gives a box that is
already occupied, and Vina then reports scores for a site that is not free.

What counts as "removable" is a scientific judgement, so it is explicit and
configurable here, and whatever is removed is itemised in the report.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from ..errors import InvalidInputError, MissingDependencyError
from ..io_utils import ensure_parent, require_file, write_text

__all__ = [
    "CleanReport",
    "clean_structure",
    "SOLVENT_RESNAMES",
    "COMMON_ADDITIVES",
    "BIOLOGICAL_IONS",
]

#: Water, in the spellings found across the PDB archive.
SOLVENT_RESNAMES = frozenset({"HOH", "DOD", "WAT", "H2O", "TIP", "SOL"})

#: Crystallisation and cryoprotection additives that are artefacts of the
#: experiment rather than part of the biological assembly.
COMMON_ADDITIVES = frozenset({
    "GOL",  # glycerol
    "EDO",  # ethylene glycol
    "PEG", "PG4", "PGE", "1PE", "P6G", "2PE",  # polyethylene glycols
    "DMS",  # dimethyl sulfoxide
    "MPD",  # 2-methyl-2,4-pentanediol
    "SO4", "PO4", "ACT", "ACY",  # buffer anions / acetate
    "CL", "BR", "IOD",  # halides
    "NA", "K", "CS", "RB",  # alkali counter-ions
    "FMT", "TRS", "EPE", "MES", "IMD", "BME", "TLA", "CIT", "MLI",
    "NO3", "AZI", "SCN", "FLC", "UNX", "UNL",
})

#: Ions and metals that are frequently catalytic or structural, and whose
#: removal would change the chemistry of the site.
BIOLOGICAL_IONS = frozenset({
    "ZN", "MG", "MN", "FE", "FE2", "FES", "CA", "CU", "CU1", "NI", "CO",
    "MO", "W", "V", "CD", "HG",
})


@dataclass
class CleanReport:
    """What cleaning produced and what it discarded."""

    input_path: str
    output_path: str
    chains_in_input: list[str]
    chains_kept: list[str]
    residues_kept: int
    atoms_kept: int
    waters_removed: int
    heterogens_removed: list[dict[str, Any]] = field(default_factory=list)
    heterogens_kept: list[dict[str, Any]] = field(default_factory=list)
    altloc_atoms_dropped: int = 0
    hydrogens_removed: int = 0
    insertion_codes_present: bool = False
    warnings: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        from dataclasses import asdict
        return asdict(self)


def _biopython():
    try:
        from Bio.PDB import PDBIO, PDBParser, Select  # noqa: PLC0415
        return PDBParser, PDBIO, Select
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise MissingDependencyError(
            "biopython",
            install="run the skill with 'uv run', which resolves the "
                    "dependencies declared in the script's PEP 723 header, "
                    "or install it directly with 'pip install biopython'",
        ) from exc


#: Header records carried from the source into the cleaned file. Between
#: them these say what the structure *is*, which cannot be recovered from
#: coordinates alone.
#:
#: ``EXPDTA`` and ``REMARK   2 RESOLUTION`` mark an experimental structure;
#: ``TITLE`` is where the AlphaFold DB and ESMFold announce a prediction.
#: Without them :func:`agskills.targets.quality.detect_scale` has only
#: magnitude to go on, and a well-refined crystal structure (mean B-factor
#: 20-50 A^2) is indistinguishable from a poor predicted model (mean pLDDT
#: 20-50) - the case that made a 2.0 A structure of HIV-1 protease report
#: "not reliable enough for docking".
#:
#: ``CRYST1`` is kept because downstream tools read it: GROMACS ``pdb2gmx``
#: and ``editconf`` both use the unit cell when one is present.
_PROVENANCE_RECORDS = (
    "HEADER",
    "TITLE ",
    "EXPDTA",
    "CRYST1",
    "REMARK   2",
)

#: Ceiling on carried header lines, so a pathological entry cannot bloat
#: the cleaned file. The records above run to a handful of lines normally.
_MAX_PROVENANCE_LINES = 40


def _carry_provenance(source: Path, destination: Path) -> None:
    """Prepend the source's provenance records to a cleaned PDB.

    ``PDBIO.save`` writes coordinates only, so the cleaned file would
    otherwise be anonymous. See :data:`_PROVENANCE_RECORDS` for why that
    matters.
    """
    carried: list[str] = []
    for line in source.read_text(encoding="utf-8",
                                 errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM", "MODEL")):
            break
        if line.startswith(_PROVENANCE_RECORDS):
            carried.append(line.rstrip())
            if len(carried) >= _MAX_PROVENANCE_LINES:
                break
    if not carried:
        return

    # The PDB format fixes the order of these records, and every REMARK
    # belongs before CRYST1. Appending the cleaning note to the end would
    # put it after the unit cell, which is a malformed file even though
    # tolerant parsers accept it.
    note = ("REMARK 500 CLEANED BY AGSKILLS: SOLVENT AND NON-BIOLOGICAL "
            "HETEROGENS REMOVED")
    insert_at = next((i for i, line in enumerate(carried)
                      if line.startswith("CRYST1")), len(carried))
    carried.insert(insert_at, note)

    body = destination.read_text(encoding="utf-8", errors="replace")
    write_text("\n".join(carried) + "\n" + body, destination)


def clean_structure(
    input_path: str | Path,
    output_path: str | Path,
    *,
    chain: str | None = None,
    keep_ligands: Iterable[str] | None = None,
    keep_metals: bool = True,
    keep_waters: bool = False,
    remove_hydrogens: bool = False,
    model_index: int = 0,
) -> CleanReport:
    """Clean a PDB file for docking and write the result.

    Args:
        input_path: Source PDB file.
        output_path: Destination for the cleaned PDB.
        chain: Keep only this chain id. ``None`` keeps every chain.
        keep_ligands: Residue names to retain even though they are
            heterogens, e.g. a cofactor you intend to model explicitly.
        keep_metals: Retain the ions in :data:`BIOLOGICAL_IONS`, which are
            often catalytic. Buffer counter-ions are removed regardless.
        keep_waters: Retain water molecules. Off by default; structural
            waters have to be chosen deliberately, not kept wholesale.
        remove_hydrogens: Strip hydrogens so the docking tool can add its
            own consistently.
        model_index: Which MODEL to keep from an ensemble.

    Raises:
        InvalidInputError: If the file has no chains, or the requested chain
            is absent, or nothing survives cleaning.
    """
    PDBParser, PDBIO, Select = _biopython()
    source = require_file(input_path)
    keep_ligand_set = {r.strip().upper() for r in (keep_ligands or ()) if r.strip()}

    import warnings as _warnings
    from Bio import BiopythonWarning

    parser = PDBParser(QUIET=True)
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore", BiopythonWarning)
        structure = parser.get_structure("receptor", str(source))

    models = list(structure)
    if not models:
        raise InvalidInputError(f"{source} contains no models")
    if model_index >= len(models):
        raise InvalidInputError(
            f"{source} has {len(models)} model(s); model_index "
            f"{model_index} is out of range"
        )
    model = models[model_index]

    chains_in_input = [c.id.strip() or "_" for c in model.get_chains()]
    if not chains_in_input:
        raise InvalidInputError(f"{source} contains no chains")

    if chain is not None:
        wanted = chain.strip()
        if wanted not in [c.id for c in model.get_chains()]:
            raise InvalidInputError(
                f"Chain {wanted!r} is not in {source.name}. "
                f"Chains present: {', '.join(chains_in_input)}",
                hint="Omit --chain to keep every chain.",
            )

    report = CleanReport(
        input_path=str(source), output_path=str(output_path),
        chains_in_input=chains_in_input, chains_kept=[],
        residues_kept=0, atoms_kept=0, waters_removed=0,
    )
    if len(models) > 1:
        report.warnings.append(
            f"{len(models)} models present; kept model {model_index + 1} only."
        )

    removed_counter: dict[str, int] = {}
    kept_counter: dict[str, int] = {}

    def _classify_heterogen(resname: str) -> tuple[bool, str]:
        """Return ``(keep, why)`` for a heterogen residue name."""
        name = resname.strip().upper()
        if name in keep_ligand_set:
            return True, "explicitly requested with --keep-ligand"
        if name in SOLVENT_RESNAMES:
            return keep_waters, "water"
        if name in BIOLOGICAL_IONS:
            return keep_metals, "metal or catalytic ion"
        if name in COMMON_ADDITIVES:
            return False, "crystallisation additive or buffer component"
        return False, "heterogen (ligand, cofactor or modified residue)"

    class _Cleaner(Select):
        def accept_model(self, m):  # noqa: N802
            return m.id == model.id

        def accept_chain(self, c):  # noqa: N802
            return True if chain is None else c.id == chain.strip()

        def accept_residue(self, residue):  # noqa: N802
            hetflag = residue.id[0]
            if residue.id[2] not in (" ", ""):
                report.insertion_codes_present = True
            if hetflag == " ":
                return True  # standard polymer residue
            resname = residue.get_resname()
            keep, why = _classify_heterogen(resname)
            counter = kept_counter if keep else removed_counter
            counter[f"{resname}|{why}"] = counter.get(f"{resname}|{why}", 0) + 1
            if not keep and resname.strip().upper() in SOLVENT_RESNAMES:
                report.waters_removed += 1
            return keep

        def accept_atom(self, atom):  # noqa: N802
            # Keep only the primary alternate location so the receptor has
            # one unambiguous conformation.
            if atom.is_disordered() and atom.get_altloc() not in (" ", "A"):
                report.altloc_atoms_dropped += 1
                return False
            if remove_hydrogens and atom.element == "H":
                report.hydrogens_removed += 1
                return False
            return True

    destination = ensure_parent(output_path)
    io = PDBIO()
    io.set_structure(structure)
    selector = _Cleaner()
    io.save(str(destination), selector)
    _carry_provenance(source, destination)

    def _unpack(counter: dict[str, int]) -> list[dict[str, Any]]:
        out = []
        for key, count in sorted(counter.items()):
            resname, why = key.split("|", 1)
            out.append({"resname": resname, "count": count, "reason": why})
        return out

    report.heterogens_removed = _unpack(removed_counter)
    report.heterogens_kept = _unpack(kept_counter)

    # Count what actually landed in the output rather than trusting the filter.
    kept_atoms = 0
    kept_residues: set[tuple[str, int, str]] = set()
    kept_chains: set[str] = set()
    for line in destination.read_text(encoding="utf-8").splitlines():
        if line[:6] in ("ATOM  ", "HETATM"):
            kept_atoms += 1
            try:
                kept_chains.add(line[21])
                kept_residues.add((line[21], int(line[22:26]), line[26]))
            except (ValueError, IndexError):
                continue
    report.atoms_kept = kept_atoms
    report.residues_kept = len(kept_residues)
    report.chains_kept = sorted(c.strip() or "_" for c in kept_chains)

    if kept_atoms == 0:
        raise InvalidInputError(
            f"Cleaning {source.name} removed every atom.",
            hint="If the file contains only heterogens, pass --keep-ligand "
                 "with the residue names you need.",
        )
    if report.insertion_codes_present:
        report.warnings.append(
            "Residue insertion codes are present. Some force fields and "
            "docking tools renumber or reject them; check the output."
        )
    return report
