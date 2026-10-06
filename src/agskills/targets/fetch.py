"""Fetch experimental and predicted protein structures.

The AlphaFold fix here matters: the code this replaces built the download
URL by hardcoding the model version, ``AF-{acc}-F1-model_v4.pdb``. The
AlphaFold Database has since advanced to v6 and the v4 files were withdrawn,
so **every** ``--uniprot-id`` fetch returned HTTP 404. Verified against the
live service on 2026-10-06:

* ``AF-P9WJA5-F1-model_v4.pdb`` -> 404
* ``AF-P9WJA5-F1-model_v5.pdb`` -> 404
* ``AF-P9WJA5-F1-model_v6.pdb`` -> 200

:func:`resolve_alphafold_url` therefore asks the prediction API which file
is current rather than guessing, which keeps working across future version
bumps.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..errors import InvalidInputError, RemoteServiceError
from ..http import HttpClient
from ..io_utils import ensure_parent, write_text

__all__ = [
    "StructureDownload",
    "fetch_pdb",
    "fetch_alphafold",
    "resolve_alphafold_url",
    "validate_pdb_id",
    "validate_uniprot_id",
]

RCSB_FILES = "https://files.rcsb.org"
RCSB_DATA = "https://data.rcsb.org"
ALPHAFOLD = "https://alphafold.ebi.ac.uk"

# A PDB entry id is 4 characters: a digit then three alphanumerics.
_PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
# UniProt accession format per the official regular expression.
_UNIPROT_RE = re.compile(
    r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
)


@dataclass
class StructureDownload:
    """A structure file written to disk, with its provenance."""

    path: Path
    identifier: str
    source: str
    url: str
    n_bytes: int
    metadata: dict[str, Any]

    def as_dict(self) -> dict[str, Any]:
        return {
            "path": str(self.path),
            "identifier": self.identifier,
            "source": self.source,
            "url": self.url,
            "bytes": self.n_bytes,
            **({"metadata": self.metadata} if self.metadata else {}),
        }


def validate_pdb_id(pdb_id: str) -> str:
    """Normalise and validate a PDB entry id.

    >>> validate_pdb_id("6hez")
    '6HEZ'
    >>> validate_pdb_id("protein.pdb")
    Traceback (most recent call last):
    agskills.errors.InvalidInputError: ...
    """
    candidate = (pdb_id or "").strip().upper()
    if not _PDB_ID_RE.match(candidate):
        raise InvalidInputError(
            f"{pdb_id!r} is not a PDB entry id.",
            hint="A PDB id is four characters starting with a digit, e.g. 6HEZ. "
                 "Search https://www.rcsb.org if you only have a protein name.",
        )
    return candidate


def validate_uniprot_id(accession: str) -> str:
    """Normalise and validate a UniProt accession.

    >>> validate_uniprot_id("p00533")
    'P00533'
    >>> validate_uniprot_id("EGFR")
    Traceback (most recent call last):
    agskills.errors.InvalidInputError: ...
    """
    candidate = (accession or "").strip().upper()
    if not _UNIPROT_RE.match(candidate):
        raise InvalidInputError(
            f"{accession!r} is not a UniProt accession.",
            hint="Accessions look like P00533 or A0A0B4J2D5. A gene name such "
                 "as EGFR is not an accession; look it up at "
                 "https://www.uniprot.org.",
        )
    return candidate


def _looks_like_pdb(text: str) -> bool:
    """Whether *text* is plausibly a PDB coordinate file.

    Guards against a server returning an HTML error page with HTTP 200,
    which is exactly how the ZINC download step used to silently save web
    pages in place of compound files.
    """
    head = text.lstrip()[:2048].upper()
    if head.startswith("<!DOCTYPE") or head.startswith("<HTML"):
        return False
    return any(token in head for token in ("HEADER", "ATOM  ", "HETATM", "CRYST1",
                                           "MODEL ", "TITLE ", "REMARK"))


def fetch_pdb(pdb_id: str, out_dir: str | Path, *,
              client: HttpClient | None = None) -> StructureDownload:
    """Download an experimental structure from the RCSB PDB.

    Args:
        pdb_id: Four-character PDB entry id.
        out_dir: Directory to write ``<ID>.pdb`` into.
        client: Injectable HTTP client, for tests.

    Raises:
        InvalidInputError: If *pdb_id* is malformed.
        RemoteServiceError: If the download fails or is not a PDB file.
    """
    entry = validate_pdb_id(pdb_id)
    http = client or HttpClient(RCSB_FILES, qps=3.0)
    url = f"{RCSB_FILES}/download/{entry}.pdb"
    text = http.fetch_text(f"/download/{entry}.pdb")
    if not _looks_like_pdb(text):
        raise RemoteServiceError(
            f"{url} did not return a PDB file "
            f"(got {len(text)} bytes starting {text[:80]!r}).",
            hint="Some entries are distributed only as mmCIF. Download the "
                 "CIF from RCSB and convert it, or use --input-file.",
        )
    path = ensure_parent(Path(out_dir) / f"{entry}.pdb")
    write_text(text, path)
    return StructureDownload(
        path=path, identifier=entry, source="RCSB PDB", url=url,
        n_bytes=len(text.encode("utf-8")), metadata={},
    )


def resolve_alphafold_url(accession: str, *,
                          client: HttpClient | None = None) -> dict[str, Any]:
    """Ask the AlphaFold API for the current model file for *accession*.

    Returns the prediction entry, including ``pdbUrl`` and ``latestVersion``.

    Raises:
        InvalidInputError: If *accession* is malformed.
        RemoteServiceError: If the protein has no AlphaFold model.
    """
    acc = validate_uniprot_id(accession)
    http = client or HttpClient(ALPHAFOLD, qps=3.0)
    payload = http.fetch_json(f"/api/prediction/{acc}")
    if not isinstance(payload, list) or not payload:
        raise RemoteServiceError(
            f"The AlphaFold Database has no model for {acc}.",
            hint="Predict the structure instead: 'predict-structure "
                 "--sequence ...' for a quick ESMFold model, or "
                 "'generate-af2-script' for an MSA-based cluster run.",
        )
    entry = payload[0]
    if not entry.get("pdbUrl"):
        raise RemoteServiceError(
            f"The AlphaFold entry for {acc} does not list a PDB download URL.",
            hint=f"Inspect {ALPHAFOLD}/entry/{acc} directly.",
        )
    return entry


def fetch_alphafold(accession: str, out_dir: str | Path, *,
                    client: HttpClient | None = None) -> StructureDownload:
    """Download the current AlphaFold model for a UniProt accession.

    The model version is resolved from the API, never hardcoded.
    """
    acc = validate_uniprot_id(accession)
    entry = resolve_alphafold_url(acc, client=client)
    url = entry["pdbUrl"]
    files = client or HttpClient(ALPHAFOLD, qps=3.0)
    text = files.fetch_text(url)
    if not _looks_like_pdb(text):
        raise RemoteServiceError(f"{url} did not return a PDB file.")

    version = entry.get("latestVersion", "unknown")
    path = ensure_parent(Path(out_dir) / f"AF-{acc}-F1-v{version}.pdb")
    write_text(text, path)
    return StructureDownload(
        path=path,
        identifier=acc,
        source="AlphaFold Database",
        url=url,
        n_bytes=len(text.encode("utf-8")),
        metadata={
            "model_version": version,
            "entry_id": entry.get("entryId"),
            "uniprot_description": entry.get("uniprotDescription"),
            "organism": entry.get("organismScientificName"),
            "sequence_length": entry.get("uniprotEnd"),
            "caveat": (
                "A predicted model carries no ligands, cofactors, metals or "
                "waters, and no crystallographic binding site. Define the "
                "docking box from a homologous holo structure or a pocket "
                "detection run before docking."
            ),
        },
    )
