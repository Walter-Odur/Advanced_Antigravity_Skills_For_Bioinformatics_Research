"""Analysis of REINVENT 4 output CSVs.

Deliberately implemented on the standard library's ``csv`` module rather
than pandas, so analysing a run does not require the optional analysis
extra. A staged-learning CSV from a long run can be large, so rows are
streamed and only the top N are retained.

What this adds over simply sorting by score:

* **Duplicate collapsing.** An RL run re-samples the same molecule many
  times across steps. A "top 20" list taken straight from the CSV is often
  the same handful of structures repeatedly, which overstates how much was
  found. Molecules are collapsed on canonical SMILES and the number of
  times each was sampled is reported.
* **Scaffold diversity.** The single most useful diagnostic for an RL run
  is whether the agent collapsed onto one scaffold. Bemis-Murcko scaffolds
  are counted and the ratio reported.
* **Score progression.** Mean score by step shows whether learning actually
  happened, which a flat ranking cannot show.

Reference:
    Bemis GW, Murcko MA. *J Med Chem* 1996;39:2887-2893 (molecular frameworks).
"""

from __future__ import annotations

import csv
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from ..errors import InvalidInputError
from ..io_utils import require_file

__all__ = ["analyze_results", "SMILES_COLUMNS", "SCORE_COLUMNS"]

#: Column names REINVENT uses for structures, across run modes and versions.
SMILES_COLUMNS: tuple[str, ...] = (
    "SMILES", "smiles", "canonical_smiles", "Smiles", "molecule",
)
#: Column names that hold the aggregate score.
SCORE_COLUMNS: tuple[str, ...] = (
    "Score", "total_score", "score", "Total_Score", "aggregate_score",
)
#: Column names that hold the RL step index.
STEP_COLUMNS: tuple[str, ...] = ("step", "Step", "epoch", "Epoch")


def _find_column(fieldnames: list[str], candidates: tuple[str, ...]
                 ) -> str | None:
    lookup = {(f or "").strip().lower(): f for f in fieldnames}
    for candidate in candidates:
        match = lookup.get(candidate.lower())
        if match:
            return match
    return None


def _as_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(result) else result


def analyze_results(csv_path: str | Path, *, top_n: int = 20,
                    sort_by: str | None = None,
                    min_score: float | None = None,
                    compute_properties: bool = True) -> dict[str, Any]:
    """Rank and summarise a REINVENT output CSV.

    Args:
        csv_path: The CSV REINVENT wrote.
        top_n: How many unique molecules to return.
        sort_by: Column to rank on. Defaults to the detected score column.
        min_score: Discard molecules below this score before ranking.
        compute_properties: Attach RDKit descriptors to the top molecules.

    Raises:
        InvalidInputError: If the file has no header, no SMILES column, or
            the requested sort column does not exist.
    """
    path = require_file(csv_path)
    if top_n < 1:
        raise InvalidInputError("--top-n must be at least 1")

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise InvalidInputError(f"{path} has no header row")
        fieldnames = list(reader.fieldnames)

        smiles_column = _find_column(fieldnames, SMILES_COLUMNS)
        if smiles_column is None:
            raise InvalidInputError(
                f"{path} has no SMILES column.",
                hint=f"Columns present: {', '.join(fieldnames)}",
            )
        score_column = (sort_by if sort_by in fieldnames
                        else _find_column(fieldnames, SCORE_COLUMNS))
        if sort_by and sort_by not in fieldnames:
            # Be explicit rather than silently falling back.
            resolved = _find_column(fieldnames, (sort_by,))
            if resolved is None:
                raise InvalidInputError(
                    f"{path} has no column {sort_by!r}.",
                    hint=f"Columns present: {', '.join(fieldnames)}",
                )
            score_column = resolved
        step_column = _find_column(fieldnames, STEP_COLUMNS)

        # Component score columns are everything numeric that is not the
        # aggregate, the step, or bookkeeping.
        unique: dict[str, dict[str, Any]] = {}
        rows_read = 0
        invalid_rows = 0
        below_threshold = 0
        scores_by_step: dict[int, list[float]] = defaultdict(list)
        all_scores: list[float] = []

        for row in reader:
            rows_read += 1
            smiles = (row.get(smiles_column) or "").strip()
            if not smiles:
                invalid_rows += 1
                continue
            score = _as_float(row.get(score_column)) if score_column else None
            if score is not None:
                all_scores.append(score)
                if step_column:
                    step = _as_float(row.get(step_column))
                    if step is not None:
                        scores_by_step[int(step)].append(score)
                if min_score is not None and score < min_score:
                    below_threshold += 1
                    continue
            entry = unique.get(smiles)
            if entry is None:
                unique[smiles] = {
                    "smiles": smiles,
                    "score": score,
                    "times_sampled": 1,
                    "first_step": int(_as_float(row.get(step_column)) or 0)
                    if step_column else None,
                    "row": {k: v for k, v in row.items() if k != smiles_column},
                }
            else:
                entry["times_sampled"] += 1
                # Keep the best observed score for a repeatedly-sampled molecule.
                if score is not None and (entry["score"] is None
                                          or score > entry["score"]):
                    entry["score"] = score
                    entry["row"] = {k: v for k, v in row.items()
                                    if k != smiles_column}

    if not unique:
        raise InvalidInputError(
            f"{path} contained no usable rows "
            f"({rows_read} row(s) read, {invalid_rows} without a structure)."
        )

    molecules = list(unique.values())
    molecules.sort(key=lambda m: (m["score"] if m["score"] is not None else -1e9),
                   reverse=True)
    top = molecules[:top_n]

    # --- scaffold diversity --------------------------------------------
    scaffold_counts: Counter[str] = Counter()
    valid_structures = 0
    enriched: list[dict[str, Any]] = []
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold

        from ..chem.descriptors import compute_properties as descriptors
        from ..chem.smiles import rdkit_quiet

        with rdkit_quiet():
            for entry in molecules:
                mol = Chem.MolFromSmiles(entry["smiles"])
                if mol is None:
                    continue
                valid_structures += 1
                try:
                    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
                    scaffold_counts[scaffold or "(acyclic)"] += 1
                except Exception:
                    pass
            for entry in top:
                record = dict(entry)
                mol = Chem.MolFromSmiles(entry["smiles"])
                if mol is not None and compute_properties:
                    props = descriptors(mol, smiles=Chem.MolToSmiles(mol))
                    record["properties"] = props.as_dict()
                    try:
                        record["murcko_scaffold"] = \
                            MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
                    except Exception:
                        record["murcko_scaffold"] = None
                enriched.append(record)
        rdkit_available = True
    except ImportError:
        enriched = [dict(entry) for entry in top]
        rdkit_available = False

    n_unique = len(unique)
    n_scaffolds = len(scaffold_counts)
    scaffold_ratio = (n_scaffolds / valid_structures) if valid_structures else None

    progression: list[dict[str, Any]] = []
    if scores_by_step:
        for step in sorted(scores_by_step):
            values = scores_by_step[step]
            progression.append({
                "step": step,
                "n": len(values),
                "mean_score": round(statistics.fmean(values), 4),
                "max_score": round(max(values), 4),
            })

    learning_note = ""
    if len(progression) >= 4:
        first_quarter = progression[:max(1, len(progression) // 4)]
        last_quarter = progression[-max(1, len(progression) // 4):]
        start = statistics.fmean(p["mean_score"] for p in first_quarter)
        end = statistics.fmean(p["mean_score"] for p in last_quarter)
        delta = end - start
        if delta > 0.05:
            learning_note = (
                f"Mean score rose from {start:.3f} to {end:.3f} over the run, "
                "so the agent learned."
            )
        elif delta < -0.05:
            learning_note = (
                f"Mean score fell from {start:.3f} to {end:.3f}. Check the "
                "scoring definition and the learning rate."
            )
        else:
            learning_note = (
                f"Mean score was flat ({start:.3f} to {end:.3f}). The agent "
                "may already have converged, sigma may be too low, or the "
                "objective may be unreachable."
            )

    diversity_note = ""
    if scaffold_ratio is not None:
        if scaffold_ratio < 0.05:
            diversity_note = (
                f"Only {n_scaffolds} distinct Bemis-Murcko scaffold(s) across "
                f"{valid_structures} molecules ({scaffold_ratio:.1%}). The "
                "agent has collapsed onto very few frameworks; enable or "
                "strengthen the diversity filter."
            )
        elif scaffold_ratio < 0.20:
            diversity_note = (
                f"{n_scaffolds} scaffolds across {valid_structures} molecules "
                f"({scaffold_ratio:.1%}): moderate diversity."
            )
        else:
            diversity_note = (
                f"{n_scaffolds} scaffolds across {valid_structures} molecules "
                f"({scaffold_ratio:.1%}): good scaffold diversity."
            )

    return {
        "file": str(path),
        "columns": fieldnames,
        "smiles_column": smiles_column,
        "score_column": score_column,
        "step_column": step_column,
        "rows_read": rows_read,
        "rows_without_structure": invalid_rows,
        "rows_below_min_score": below_threshold,
        "unique_molecules": n_unique,
        "duplicate_samples": rows_read - invalid_rows - n_unique,
        "valid_structures": valid_structures,
        "distinct_scaffolds": n_scaffolds,
        "scaffold_diversity_ratio": (round(scaffold_ratio, 4)
                                     if scaffold_ratio is not None else None),
        "score_statistics": {
            "n": len(all_scores),
            "mean": round(statistics.fmean(all_scores), 4) if all_scores else None,
            "median": round(statistics.median(all_scores), 4) if all_scores else None,
            "max": round(max(all_scores), 4) if all_scores else None,
            "min": round(min(all_scores), 4) if all_scores else None,
            "stdev": (round(statistics.stdev(all_scores), 4)
                      if len(all_scores) > 1 else None),
        },
        "score_progression": progression,
        "most_common_scaffolds": [
            {"scaffold": s, "count": c}
            for s, c in scaffold_counts.most_common(10)
        ],
        "interpretation": {
            "learning": learning_note,
            "diversity": diversity_note,
            "note": (
                "times_sampled counts how often the agent produced each "
                "molecule; a high count on few structures means the search "
                "converged rather than explored. A REINVENT score is an "
                "aggregate of the configured objectives, not a predicted "
                "activity."
            ),
        },
        "rdkit_available": rdkit_available,
        "top_molecules": enriched,
    }
