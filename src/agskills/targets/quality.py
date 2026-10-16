"""Structure confidence assessment from the B-factor / pLDDT column.

Two corrections over the code this replaces:

* **Chain awareness.** pLDDT was collected into ``dict[int, float]`` keyed on
  residue sequence number alone. In any multi-chain file residue 42 of chain
  B overwrites residue 42 of chain A, so a dimer silently reported the
  confidence of roughly half its residues. Residues are now keyed on
  ``(chain, resseq, icode)`` and sorted, so contiguous low-confidence regions
  are genuinely contiguous.
* **A real median.** ``sorted(values)[len(values) // 2]`` is the upper of the
  two central values for an even count, not the median. ``statistics.median``
  is used instead.

pLDDT bands follow the AlphaFold DB convention (Jumper et al., *Nature*
2021;596:583-589; Varadi et al., *Nucleic Acids Res* 2022;50:D439-D444):
very high >= 90, confident 70-90, low 50-70, very low < 50.
"""

from __future__ import annotations

import statistics
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Literal

from ..errors import InvalidInputError
from ..io_utils import require_file

__all__ = [
    "ResidueConfidence",
    "PlddtReport",
    "Scale",
    "detect_scale",
    "provenance",
    "extract_plddt",
    "assess_structure",
    "find_low_confidence_regions",
]


@dataclass(frozen=True)
class ResidueConfidence:
    """One residue's CA B-factor, used as pLDDT for predicted models."""

    chain: str
    resseq: int
    icode: str
    resname: str
    plddt: float

    @property
    def key(self) -> tuple[str, int, str]:
        return (self.chain, self.resseq, self.icode)


@dataclass
class PlddtReport:
    """Confidence assessment for one structure file."""

    path: str
    n_residues: int
    chains: list[str]
    mean_plddt: float | None
    median_plddt: float | None
    min_plddt: float | None
    max_plddt: float | None
    residues_very_high: int
    residues_confident: int
    residues_low: int
    residues_very_low: int
    low_confidence_regions: list[dict[str, Any]]
    verdict: str
    reason: str
    scale: str = "plddt_0_100"
    actions: list[str] = field(default_factory=list)
    interpretation: str = ""
    b_factor_caveat: str = ""
    per_chain: dict[str, dict[str, Any]] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        data = asdict(self)
        for key in ("mean_plddt", "median_plddt", "min_plddt", "max_plddt"):
            if isinstance(data[key], float):
                data[key] = round(data[key], 2)
        return data

    def summary(self) -> dict[str, Any]:
        """The verdict and the numbers behind it, without the per-region detail.

        Embedded by other commands' reports. The full form is roughly five
        times larger, mostly ``low_confidence_regions`` and ``per_chain``,
        and an agent reading a receptor report wants the verdict rather than
        a residue-by-residue breakdown. ``assess-structure`` produces the
        full report when it is actually needed.
        """
        return {
            "verdict": self.verdict,
            "reason": self.reason,
            "scale": self.scale,
            "n_residues": self.n_residues,
            "chains": self.chains,
            "mean_plddt": (round(self.mean_plddt, 2)
                           if self.mean_plddt is not None else None),
            "residues_below_50": self.residues_very_low,
            "n_low_confidence_regions": len(self.low_confidence_regions),
            "actions": self.actions,
            "full_report_command":
                "ag-target-preparation assess-structure --pdb <file> "
                "--output quality.json",
        }


def extract_plddt(pdb_text: str) -> list[ResidueConfidence]:
    """Extract per-residue CA B-factors, keyed by chain and residue.

    Only the first MODEL of a multi-model file is read, matching how
    docking and MD treat an NMR ensemble by default.
    """
    residues: dict[tuple[str, int, str], ResidueConfidence] = {}
    for line in pdb_text.splitlines():
        record = line[:6]
        if record == "ENDMDL":
            break  # only the first model
        if record not in ("ATOM  ", "HETATM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        try:
            chain = line[21].strip() or "_"
            resseq = int(line[22:26])
            icode = line[26].strip()
            resname = line[17:20].strip()
            bfactor = float(line[60:66])
        except (ValueError, IndexError):
            continue
        key = (chain, resseq, icode)
        # Keep the first altloc encountered, which is conventionally 'A'.
        if key not in residues:
            residues[key] = ResidueConfidence(chain, resseq, icode, resname,
                                              bfactor)
    return sorted(residues.values(), key=lambda r: (r.chain, r.resseq, r.icode))


def find_low_confidence_regions(residues: list[ResidueConfidence],
                                threshold: float = 50.0,
                                min_length: int = 1) -> list[dict[str, Any]]:
    """Find contiguous runs of residues below *threshold*, per chain.

    A run is broken by a chain change or a gap in residue numbering, so a
    region is always genuinely contiguous in sequence.
    """
    regions: list[dict[str, Any]] = []
    run: list[ResidueConfidence] = []

    def _flush() -> None:
        if len(run) >= min_length:
            regions.append({
                "chain": run[0].chain,
                "start_residue": run[0].resseq,
                "end_residue": run[-1].resseq,
                "length": len(run),
                "mean_plddt": round(statistics.fmean(r.plddt for r in run), 1),
            })
        run.clear()

    previous: ResidueConfidence | None = None
    for residue in residues:
        contiguous = (
            previous is not None
            and residue.chain == previous.chain
            and residue.resseq in (previous.resseq, previous.resseq + 1)
        )
        if residue.plddt < threshold:
            if run and not contiguous:
                _flush()
            run.append(residue)
        elif run:
            _flush()
        previous = residue
    if run:
        _flush()
    return regions


#: How the B-factor column was interpreted.
Scale = Literal["plddt_0_100", "plddt_0_1", "bfactor"]


#: Header records that identify an *experimental* structure. ``EXPDTA`` is
#: mandatory in a PDB-format entry from the wwPDB, and ``REMARK   2`` carries
#: the diffraction resolution.
_EXPERIMENTAL_MARKERS = (
    "EXPDTA",
    "REMARK   2 RESOLUTION",
    "REMARK   3   PROGRAM",
    "ANISOU",
)

#: Strings that identify a *predicted* model. The AlphaFold Database and the
#: ESMFold API both say so in ``TITLE``.
_PREDICTED_MARKERS = (
    "ALPHAFOLD",
    "ESMFOLD",
    "COLABFOLD",
    "OPENFOLD",
    "ROSETTAFOLD",
    "BOLTZ",
    "CHAI-1",
    "PREDICTED STRUCTURE",
)


def provenance(pdb_text: str) -> str | None:
    """Classify a coordinate file from its header records.

    Returns ``"experimental"``, ``"predicted"`` or ``None`` when the header
    says neither - which happens whenever a file has been stripped of
    everything but coordinates.

    Only the header is inspected, so this is cheap and unaffected by the
    coordinates themselves.
    """
    head = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM", "MODEL")):
            break
        head.append(line.upper())
    header = "\n".join(head)
    if not header:
        # ANISOU records appear among the coordinates, not before them, and
        # only an experimental refinement produces them.
        return "experimental" if "\nANISOU" in pdb_text else None
    if any(m in header for m in _PREDICTED_MARKERS):
        return "predicted"
    if any(m in header for m in _EXPERIMENTAL_MARKERS):
        return "experimental"
    if "ANISOU" in pdb_text:
        return "experimental"
    return None


def detect_scale(residues: list[ResidueConfidence],
                 pdb_text: str | None = None) -> Scale:
    """Decide what the B-factor column actually contains.

    Three cases occur in practice and they are not interchangeable:

    ``plddt_0_100``
        A predicted model from the AlphaFold Database or ColabFold. Higher
        is better, bounded at 100.
    ``plddt_0_1``
        A predicted model from the **ESMFold API**, which writes the same
        quantity as a fraction. Verified against the live endpoint: a
        120-residue fold came back with B-factors from 0.40 to 0.94, mean
        0.84 - that is an 84 pLDDT model, not a 0.84 one.
    ``bfactor``
        An experimental structure. The column holds atomic displacement
        parameters: unbounded and **lower is better**, so the pLDDT bands
        do not apply at all.

    **Magnitude alone cannot decide this**, which is why *pdb_text* is read
    when it is available. The two distributions overlap in the range that
    matters most: a well-refined crystal structure has a mean B-factor of
    roughly 20-50 A^2, and so does a poor predicted model's pLDDT. Measured
    on six real entries, only two (1M17 at 120.0 and 6HEZ at 155.8) have any
    B-factor above 100; 1HSG (2.0 A, mean 29.3), 4EY7 (mean 31.2), 6LU7
    (mean 40.9) and 2HYY (mean 42.2) all stay under it. Judged on numbers
    alone every one of those four is read as a failed prediction and
    reported as "not reliable enough for docking" - the better the
    crystallography, the lower the B-factors, and the more confidently the
    structure is rejected.

    The header settles it: a wwPDB entry carries ``EXPDTA`` and
    ``REMARK   2 RESOLUTION``, while the AlphaFold DB and ESMFold both
    announce themselves in ``TITLE``. Magnitude is used only as a fallback,
    for files stripped to bare coordinates.

    Args:
        residues: Per-residue values from :func:`extract_plddt`.
        pdb_text: The file's full text. Supply it whenever possible.
    """
    if not residues:
        return "bfactor"
    values = [r.plddt for r in residues]
    highest = max(values)

    # A fraction: every value in [0, 1]. No crystal structure looks like
    # this, so it is decisive on its own.
    if highest <= 1.0 and len(values) >= 5:
        return "plddt_0_1"

    if pdb_text is not None:
        source = provenance(pdb_text)
        if source == "experimental":
            return "bfactor"
        if source == "predicted":
            # Trust the header over the numbers, but a fraction is still a
            # fraction.
            return "plddt_0_1" if highest <= 1.0 else "plddt_0_100"

    # --- no usable header: fall back on magnitude ----------------------
    # Unbounded values can only be B-factors.
    if highest > 100.0:
        return "bfactor"
    # A predicted model worth docking into is mostly confident, so its mean
    # sits well above the 20-50 band where crystallographic B-factors live.
    # 70 is the AlphaFold "confident" threshold and the same number the
    # PROCEED verdict uses, so a file that cannot be identified is read as a
    # prediction only when it would actually pass.
    if statistics.fmean(values) >= 70.0:
        return "plddt_0_100"
    return "bfactor"


def assess_structure(pdb_path: str | Path) -> PlddtReport:
    """Assess a structure's per-residue confidence and recommend an action.

    The verdict is one of:

    ``PROCEED``
        Mean pLDDT >= 70 with no large disordered region; suitable for docking.
    ``REFINE``
        Usable in part. Specific remedial actions are listed.
    ``STOP``
        Mean pLDDT < 50, or over half the residues below 50.

    Raises:
        InvalidInputError: If the file contains no CA atoms.
    """
    path = require_file(pdb_path)
    text = path.read_text(encoding="utf-8", errors="replace")
    residues = extract_plddt(text)

    if not residues:
        raise InvalidInputError(
            f"{path} contains no CA atoms, so confidence cannot be assessed.",
            hint="Check the file is a protein coordinate file and not an "
                 "mmCIF, an error page, or a ligand-only extract.",
        )

    scale = detect_scale(residues, text)
    if scale == "plddt_0_1":
        # Rescale to the 0-100 convention everything downstream assumes,
        # so the bands and the verdict mean what they say.
        residues = [
            ResidueConfidence(r.chain, r.resseq, r.icode, r.resname,
                              r.plddt * 100.0)
            for r in residues
        ]
    values = [r.plddt for r in residues]
    total = len(values)
    mean = statistics.fmean(values)
    chains = sorted({r.chain for r in residues})

    bands = {
        "very_high": sum(1 for v in values if v >= 90),
        "confident": sum(1 for v in values if 70 <= v < 90),
        "low": sum(1 for v in values if 50 <= v < 70),
        "very_low": sum(1 for v in values if v < 50),
    }
    regions = find_low_confidence_regions(residues, threshold=50.0)

    predicted = scale in ("plddt_0_100", "plddt_0_1")
    caveat = {
        "plddt_0_100": (
            "B-factors are on a 0-100 scale with a high mean, consistent "
            "with pLDDT from a predicted model."
        ),
        "plddt_0_1": (
            "B-factors were all between 0 and 1, which is how the ESMFold "
            "API writes pLDDT. They have been rescaled by 100, so the "
            "figures below are on the usual 0-100 scale. Without that "
            "rescaling an 84 pLDDT model reads as 0.84 and would be "
            "reported as very low confidence."
        ),
        "bfactor": (
            "B-factors do not look like pLDDT (values above 100, or a low "
            "mean). This is most likely an experimental structure, where "
            "the B-factor column holds atomic displacement parameters and "
            "LOWER is better. The pLDDT verdict below does not apply; "
            "judge the structure by its resolution and R-free instead."
        ),
    }[scale]

    per_chain = {
        chain: {
            "n_residues": sum(1 for r in residues if r.chain == chain),
            "mean_plddt": round(
                statistics.fmean(r.plddt for r in residues if r.chain == chain), 2
            ),
        }
        for chain in chains
    }

    pct_very_low = bands["very_low"] / total * 100
    pct_usable = (bands["very_high"] + bands["confident"]) / total * 100
    actions: list[str] = []

    if not predicted:
        verdict = "NOT_APPLICABLE"
        reason = ("The B-factor column does not contain pLDDT, so a "
                  "confidence verdict would be meaningless.")
        actions = [
            "Treat this as an experimental structure: check resolution, "
            "R-free and the occupancy of the binding site residues.",
            "Remove waters and crystallisation additives with "
            "'prepare-receptor' before docking.",
        ]
    elif mean < 50 or pct_very_low > 50:
        verdict = "STOP"
        reason = (f"Mean pLDDT {mean:.1f} with {pct_very_low:.0f}% of residues "
                  "below 50. Not reliable enough for docking.")
        actions = [
            "Re-predict with an MSA-based method: 'generate-af2-script "
            "--method colabfold'.",
            "Search the PDB for a homologous experimental structure.",
            "If no structure can be obtained, consider ligand-based virtual "
            "screening instead of docking.",
        ]
    else:
        long_regions = [r for r in regions if r["length"] >= 5]
        needs_work = mean < 70 or bool(long_regions) or (
            bands["low"] / total > 0.30
        )
        if needs_work:
            verdict = "REFINE"
            reason = (f"Mean pLDDT {mean:.1f}; {pct_usable:.0f}% of residues at "
                      f"70 or above, {len(long_regions)} disordered "
                      "region(s) of 5+ residues.")
            if long_regions:
                actions.append(
                    "Truncate or exclude the disordered regions: "
                    + ", ".join(
                        f"chain {r['chain']} {r['start_residue']}-{r['end_residue']}"
                        for r in long_regions[:6]
                    )
                )
            actions.append(
                "Energy-minimise before docking to relieve strain: "
                "md-simulation gmxapi-setup --production-ns 0"
            )
            if mean < 80:
                actions.append(
                    "ESMFold is single-sequence; an MSA-based AlphaFold2 run "
                    "usually improves accuracy materially."
                )
        else:
            verdict = "PROCEED"
            reason = (f"Mean pLDDT {mean:.1f} with {pct_usable:.0f}% of "
                      "residues at 70 or above.")
            if bands["very_low"]:
                actions.append(
                    f"{bands['very_low']} residue(s) are below 50; keep the "
                    "docking box away from them."
                )

    return PlddtReport(
        path=str(path),
        n_residues=total,
        chains=chains,
        mean_plddt=mean,
        median_plddt=statistics.median(values),
        min_plddt=min(values),
        max_plddt=max(values),
        residues_very_high=bands["very_high"],
        residues_confident=bands["confident"],
        residues_low=bands["low"],
        residues_very_low=bands["very_low"],
        low_confidence_regions=regions,
        verdict=verdict,
        reason=reason,
        scale=scale,
        actions=actions,
        interpretation=(
            "pLDDT bands: >=90 very high, 70-90 confident, 50-70 low, "
            "<50 very low (AlphaFold DB convention)."
        ),
        b_factor_caveat=caveat,
        per_chain=per_chain,
    )
