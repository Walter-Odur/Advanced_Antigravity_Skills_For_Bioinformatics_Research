"""Structure prediction from sequence via the ESMFold API.

Corrections over the code this replaces:

* The length limit was inconsistent: a module constant declared
  ``ESMFOLD_MAX_LENGTH = 400`` but was never referenced, while the validator
  accepted up to 2500 residues. Sequences between those bounds passed
  validation and then failed at the API with an opaque error.
* The response was written to disk without checking it was a PDB file, so
  an HTML error page served with HTTP 200 became ``predicted.pdb``.
* The quality report was written with ``Path.write_text`` and no encoding,
  which raises ``UnicodeEncodeError`` on a default Windows console codepage
  as soon as the report contains a non-ASCII character.

Reference:
    Lin Z et al. *Science* 2023;379:1123-1130.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..errors import InvalidInputError, RemoteServiceError
from ..http import HttpClient
from ..io_utils import ensure_parent, require_file, write_text

__all__ = [
    "ESMFOLD_URL",
    "ESMFOLD_MAX_RESIDUES",
    "SequenceInfo",
    "read_fasta",
    "validate_sequence",
    "predict_structure",
]

ESMFOLD_URL = "https://api.esmatlas.com"
_FOLD_PATH = "/foldSequence/v1/pdb/"

#: The public ESMFold endpoint rejects sequences longer than this. Enforced
#: locally so the user gets a clear message and a route forward instead of a
#: server-side failure after a long wait.
ESMFOLD_MAX_RESIDUES = 400
MIN_RESIDUES = 10

#: The 20 standard amino acids. Ambiguity codes are rejected explicitly
#: because ESMFold cannot fold them.
STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AA = {
    "B": "aspartate or asparagine", "Z": "glutamate or glutamine",
    "J": "leucine or isoleucine", "X": "any residue",
    "U": "selenocysteine", "O": "pyrrolysine",
}


@dataclass
class SequenceInfo:
    """A validated amino-acid sequence."""

    sequence: str
    length: int
    header: str = ""
    source: str = ""
    warnings: list[str] = field(default_factory=list)


def read_fasta(fasta_path: str | Path) -> SequenceInfo:
    """Read the first record from a FASTA file.

    Raises:
        InvalidInputError: If the file holds no sequence.
    """
    path = require_file(fasta_path)
    header = ""
    chunks: list[str] = []
    records = 0
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith((">", ";")):
            records += 1
            if records > 1:
                break
            header = line[1:].strip()
        elif records <= 1:
            chunks.append("".join(line.split()).upper())

    sequence = "".join(chunks)
    if not sequence:
        raise InvalidInputError(
            f"{path} contains no sequence data.",
            hint="A FASTA file needs a '>' header line followed by residues.",
        )
    info = SequenceInfo(sequence=sequence, length=len(sequence), header=header,
                        source=str(path))
    if records > 1:
        info.warnings.append(
            f"{path.name} holds more than one record; only the first "
            f"({header or 'unnamed'}) was read."
        )
    return info


def validate_sequence(sequence: str, *, max_residues: int = ESMFOLD_MAX_RESIDUES
                      ) -> SequenceInfo:
    """Validate a sequence for ESMFold and return it normalised.

    Raises:
        InvalidInputError: If the sequence is empty, too short, too long, or
            contains characters ESMFold cannot fold.
    """
    cleaned = "".join((sequence or "").split()).upper()
    # Tolerate a trailing stop codon marker from translated nucleotide input.
    cleaned = cleaned.rstrip("*")

    if not cleaned:
        raise InvalidInputError("The sequence is empty.")
    if len(cleaned) < MIN_RESIDUES:
        raise InvalidInputError(
            f"The sequence has {len(cleaned)} residues; at least "
            f"{MIN_RESIDUES} are needed for a meaningful fold."
        )
    if len(cleaned) > max_residues:
        raise InvalidInputError(
            f"The sequence has {len(cleaned)} residues but the public ESMFold "
            f"endpoint accepts at most {max_residues}.",
            hint="Use 'generate-af2-script' to run ColabFold or AlphaFold2 on "
                 "a cluster, or fold a single domain rather than the full "
                 "chain.",
        )

    # A nucleotide sequence cannot be caught by checking characters: A, C,
    # G and T are alanine, cysteine, glycine and threonine, and N is
    # asparagine, so "ACGTACGT..." is a perfectly valid - if improbable -
    # protein sequence. Folding DNA as protein produces confident-looking
    # nonsense, so the composition is checked instead.
    # At least three distinct nucleotide letters are required, so a
    # low-complexity protein such as poly-alanine is not misread as DNA.
    alphabet = set(cleaned)
    for label, letters in (("DNA", set("ACGTN")), ("RNA", set("ACGUN"))):
        if alphabet <= letters and len(alphabet) >= 3 and len(cleaned) >= 20:
            raise InvalidInputError(
                f"This looks like a {label} sequence, not a protein: it uses "
                f"only the letters {''.join(sorted(alphabet))}.",
                hint=("Translate it to protein before folding. If it really "
                      "is a protein composed only of these residues, fold it "
                      "through the library API, which does not apply this "
                      "check."),
            )

    unknown = sorted(set(cleaned) - STANDARD_AA)
    if unknown:
        ambiguous = [c for c in unknown if c in AMBIGUOUS_AA]
        detail = ", ".join(
            f"{c} ({AMBIGUOUS_AA[c]})" if c in AMBIGUOUS_AA else repr(c)
            for c in unknown
        )
        raise InvalidInputError(
            f"The sequence contains residues ESMFold cannot fold: {detail}",
            hint=("Replace ambiguity codes with a definite residue."
                  if ambiguous else
                  "Check that this is a protein sequence, not DNA or RNA."),
        )
    return SequenceInfo(sequence=cleaned, length=len(cleaned), source="argument")


def _looks_like_pdb(text: str) -> bool:
    head = text.lstrip()[:2048].upper()
    if head.startswith("<!DOCTYPE") or head.startswith("<HTML"):
        return False
    return "ATOM  " in head or head.startswith("HEADER") or "MODEL " in head


def predict_structure(sequence: str, output_path: str | Path, *,
                      client: HttpClient | None = None,
                      timeout: float = 300.0) -> dict[str, Any]:
    """Fold *sequence* with ESMFold and write the PDB to *output_path*.

    Returns metadata including the prediction's pLDDT assessment.

    Raises:
        InvalidInputError: If the sequence is not foldable.
        RemoteServiceError: If the API fails or returns a non-PDB body.
    """
    info = validate_sequence(sequence)
    http = client or HttpClient(ESMFOLD_URL, qps=0.5, timeout=timeout,
                                deadline=timeout * 2)
    body = http.fetch_text(
        _FOLD_PATH,
        method="POST",
        data=info.sequence.encode("ascii"),
        headers={"Content-Type": "text/plain"},
        timeout=timeout,
    )
    if not _looks_like_pdb(body):
        raise RemoteServiceError(
            "The ESMFold endpoint did not return a PDB structure "
            f"(received {len(body)} bytes starting {body[:120]!r}).",
            hint="The service is intermittently unavailable. Retry, or use "
                 "'generate-af2-script' for an offline cluster run.",
        )

    destination = ensure_parent(output_path)
    write_text(body, destination)

    # Assess the model we just wrote, so the caller gets the verdict in the
    # same step rather than having to remember to run assess-structure.
    from .quality import assess_structure

    report = assess_structure(destination)
    return {
        "method": "ESMFold v1 (single-sequence, no MSA)",
        "api": f"{ESMFOLD_URL}{_FOLD_PATH}",
        "sequence_length": info.length,
        "output_pdb": str(destination),
        "pdb_bytes": len(body.encode("utf-8")),
        "quality": report.as_dict(),
        "citation": "Lin Z et al. Science 2023;379:1123-1130",
        "limitation": (
            "ESMFold predicts from a single sequence with no multiple "
            "sequence alignment. It is fast but typically less accurate than "
            "AlphaFold2 on multi-domain proteins and on targets with shallow "
            "evolutionary coverage."
        ),
        "warnings": info.warnings,
    }
