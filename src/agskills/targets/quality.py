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


def detect_scale(residues: list[ResidueConfidence]) -> Scale:
    """Decide what the B-factor column actually contains.

    Three cases occur in practice and they are not interchangeable:

    ``plddt_0_100``
        A predicted model from the AlphaFold Database or ColabFold. Higher
        is better, bounded at 100.
    ``plddt_0_1``
        A predicted model from the **ESMFold API**, which writes the same
        quantity as a fraction. Verified against the live endpoint: a
        120-residue fold came back with B-factors from 0.40 to 0.94, mean
        0.84 - that is an 84 pLDDT model, not a 0.84 one. Treating the
        fraction as a 0-100 score reports a good model as "very low
        confidence, not recommended for docking", which is how the code
        this replaces behaved.
    ``bfactor``
        An experimental structure. The column holds atomic displacement
        parameters: unbounded, routinely above 100, and **lower is
        better**, so the pLDDT bands do not apply at all.

    Guessing wrong in either direction produces a confident and wrong
    verdict, so the detected scale is reported alongside it.
    """
    if not residues:
        return "bfactor"
    values = [r.plddt for r in residues]
    highest = max(values)

    # A fraction: every value in [0, 1] with a plausible spread. A crystal
    # structure never has every B-factor below 1.
    if highest <= 1.0 and len(values) >= 5:
        return "plddt_0_1"
    if highest > 100.0:
        return "bfactor"
    # Crystallographic B-factors cluster low; a usable predicted model does
    # not. 20 separates the two in practice.
    if statistics.fmean(values) >= 20.0:
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

    scale = detect_scale(residues)
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
