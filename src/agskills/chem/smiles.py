"""SMILES parsing, validation and compound-list loading.

Centralises three things the original scripts repeated inconsistently:

* **Silent truncation.** ``list(zip(args.smiles, names))`` was used to pair
  compounds with names. Passing three SMILES and one name silently dropped
  two compounds - the run reported success on one third of the input.
  :func:`load_compounds` raises instead.
* **RDKit console noise.** RDKit writes parse errors straight to stderr,
  which corrupts machine-readable output. Parsing is wrapped so failures
  become structured data.
* **Name defaults.** Some call sites numbered compounds from 0, others
  from 1. Names are now always ``compound_1 ... compound_N``.
"""

from __future__ import annotations

import contextlib
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

from ..errors import InvalidInputError, UsageError
from ..io_utils import read_lines, require_file

__all__ = [
    "CompoundRecord",
    "parse_smiles",
    "standardize_smiles",
    "load_compounds",
    "load_smiles_file",
    "dedupe",
    "rdkit_quiet",
]


def _rdkit():
    """Import RDKit, raising a MissingDependencyError with install advice."""
    try:
        from rdkit import Chem  # noqa: PLC0415
        return Chem
    except ImportError as exc:  # pragma: no cover - environment dependent
        from ..errors import MissingDependencyError
        raise MissingDependencyError("rdkit", install="pip install rdkit") from exc


@contextlib.contextmanager
def rdkit_quiet() -> Iterator[None]:
    """Suppress RDKit's C++ logging for the duration of the block.

    RDKit prints ``SMILES Parse Error`` directly to the process stderr from
    C++. When a command's contract is "stdout/stderr is a status line, the
    JSON file is the result", that noise is still undesirable, and it makes
    test output unreadable when deliberately feeding in bad input.
    """
    try:
        from rdkit import RDLogger
    except ImportError:  # pragma: no cover
        yield
        return
    RDLogger.DisableLog("rdApp.*")
    try:
        yield
    finally:
        RDLogger.EnableLog("rdApp.*")


@dataclass(frozen=True)
class CompoundRecord:
    """One input compound: a name, the SMILES as supplied, and its origin."""

    name: str
    smiles: str
    source: str = "input"

    def __post_init__(self) -> None:
        if not self.smiles or not self.smiles.strip():
            raise InvalidInputError(f"Compound {self.name!r} has an empty SMILES")


def parse_smiles(smiles: str, *, sanitize: bool = True):
    """Parse *smiles* and return an RDKit ``Mol``, or ``None`` if invalid.

    Never raises for malformed input and never writes to stderr, so callers
    can report invalid compounds as data.

    >>> parse_smiles("CCO") is not None
    True
    >>> parse_smiles("not a molecule") is None
    True
    >>> parse_smiles("") is None
    True
    """
    Chem = _rdkit()
    if not smiles or not smiles.strip():
        return None
    with rdkit_quiet():
        try:
            return Chem.MolFromSmiles(smiles.strip(), sanitize=sanitize)
        except Exception:
            return None


def standardize_smiles(smiles: str) -> str | None:
    """Return the RDKit canonical SMILES for *smiles*, or ``None``.

    Canonicalisation is what makes de-duplication meaningful: ``C1=CC=CC=C1``
    and ``c1ccccc1`` are the same molecule and must collapse to one entry.

    >>> standardize_smiles("C1=CC=CC=C1")
    'c1ccccc1'
    >>> standardize_smiles("OCC") == standardize_smiles("CCO")
    True
    >>> standardize_smiles("Xx") is None
    True
    """
    Chem = _rdkit()
    mol = parse_smiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def load_compounds(
    smiles: Iterable[str] | None = None,
    names: Iterable[str] | None = None,
    *,
    csv_path: str | None = None,
    smiles_col: str = "SMILES",
    name_col: str = "Name",
    smiles_file: str | None = None,
) -> list[CompoundRecord]:
    """Build a compound list from inline SMILES, a CSV, or a SMILES file.

    Exactly one source must be supplied.

    Raises:
        UsageError: If no source, or more than one source, is given, or if
            ``names`` is provided with a length that does not match
            ``smiles``.
        InvalidInputError: If a CSV lacks the requested SMILES column.
    """
    sources = [s for s in (smiles, csv_path, smiles_file) if s]
    if not sources:
        raise UsageError(
            "No compounds supplied.",
            hint="Pass --smiles, --csv or --smiles-file.",
        )
    if len(sources) > 1:
        raise UsageError(
            "Supply compounds from exactly one source.",
            hint="--smiles, --csv and --smiles-file are mutually exclusive.",
        )

    if csv_path:
        return _load_from_csv(csv_path, smiles_col, name_col)
    if smiles_file:
        return load_smiles_file(smiles_file)

    smiles_list = [s for s in (smiles or [])]
    if not smiles_list:
        raise UsageError("--smiles was given but contained no values.")

    if names is None:
        name_list = [f"compound_{i}" for i in range(1, len(smiles_list) + 1)]
    else:
        name_list = list(names)
        if len(name_list) != len(smiles_list):
            # The original code used zip() here, which silently discarded
            # the excess and reported success on a truncated run.
            raise UsageError(
                f"Got {len(smiles_list)} SMILES but {len(name_list)} names; "
                "they must correspond one to one.",
                hint="Omit --names to have them generated automatically.",
            )
    return [
        CompoundRecord(name=n, smiles=s, source="argument")
        for n, s in zip(name_list, smiles_list)
    ]


def _load_from_csv(csv_path: str, smiles_col: str, name_col: str) -> list[CompoundRecord]:
    path = require_file(csv_path)
    records: list[CompoundRecord] = []
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise InvalidInputError(f"{path} is empty or has no header row")
        # Match the header case-insensitively; REINVENT writes "SMILES",
        # ChEMBL exports "canonical_smiles", user files vary.
        lookup = {(c or "").strip().lower(): c for c in reader.fieldnames}
        smi_key = lookup.get(smiles_col.strip().lower())
        if smi_key is None:
            raise InvalidInputError(
                f"{path} has no column {smiles_col!r}. "
                f"Columns present: {', '.join(reader.fieldnames)}",
                hint="Select the right column with --smiles-col.",
            )
        name_key = lookup.get(name_col.strip().lower())
        for index, row in enumerate(reader, start=1):
            smi = (row.get(smi_key) or "").strip()
            if not smi:
                continue
            raw_name = (row.get(name_key) or "").strip() if name_key else ""
            records.append(
                CompoundRecord(
                    name=raw_name or f"compound_{index}",
                    smiles=smi,
                    source=f"{path.name}:{index}",
                )
            )
    if not records:
        raise InvalidInputError(
            f"No usable rows in {path} (column {smiles_col!r} was empty "
            "in every row)"
        )
    return records


def load_smiles_file(smiles_file: str) -> list[CompoundRecord]:
    """Load a ``.smi`` file: one SMILES per line, optional whitespace-separated name.

    This is the SMILES format REINVENT and ZINC both use.
    """
    path = require_file(smiles_file)
    records: list[CompoundRecord] = []
    for index, line in enumerate(read_lines(path), start=1):
        parts = line.split()
        if not parts:
            continue
        smi = parts[0]
        # Skip a header row such as "smiles zinc_id" from a ZINC tranche.
        if index == 1 and smi.lower() in ("smiles", "smi", "canonical_smiles"):
            continue
        name = parts[1] if len(parts) > 1 else f"compound_{index}"
        records.append(CompoundRecord(name=name, smiles=smi,
                                      source=f"{path.name}:{index}"))
    if not records:
        raise InvalidInputError(f"No SMILES found in {path}")
    return records


def dedupe(records: Iterable[CompoundRecord]) -> tuple[list[CompoundRecord], list[dict]]:
    """Drop duplicate and invalid structures, canonicalising as we go.

    Returns ``(kept, rejected)``. Each rejected entry records why, so a
    caller can report exactly what happened to every input row rather than
    quietly shrinking the dataset.
    """
    kept: list[CompoundRecord] = []
    rejected: list[dict] = []
    seen: dict[str, str] = {}

    for record in records:
        canonical = standardize_smiles(record.smiles)
        if canonical is None:
            rejected.append({
                "name": record.name,
                "smiles": record.smiles,
                "source": record.source,
                "reason": "invalid SMILES",
            })
            continue
        if canonical in seen:
            rejected.append({
                "name": record.name,
                "smiles": record.smiles,
                "source": record.source,
                "reason": f"duplicate of {seen[canonical]}",
                "canonical_smiles": canonical,
            })
            continue
        seen[canonical] = record.name
        kept.append(CompoundRecord(name=record.name, smiles=canonical,
                                   source=record.source))
    return kept, rejected
