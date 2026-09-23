"""De novo protein structure prediction for Drug Discovery Wizard.

Predicts protein structures from amino acid sequences when no experimental
structure exists in PDB or AlphaFold Database.

Two-tier approach:
  1. ESMFold API (fast, free, single-sequence) — for quick predictions
  2. ColabFold/AF2 HPC script generation — for high-quality MSA-based predictions

Usage:
    # Quick prediction via ESMFold API
    python predict_structure.py predict --sequence MKTLLILAVV... --output structure.pdb

    # From FASTA file
    python predict_structure.py predict --fasta protein.fasta --output structure.pdb

    # Generate AF2/ColabFold HPC script for full MSA-based prediction
    python predict_structure.py generate-af2-script --fasta protein.fasta \
        --output submit_af2.sh --gpu a100 --db-path /data/af2_dbs

    # Assess prediction quality (pLDDT from B-factor column)
    python predict_structure.py assess --pdb structure.pdb --output quality.json
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ESMFOLD_API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ESMFOLD_MAX_LENGTH = 400  # ESMFold API limit for single sequences
REQUEST_TIMEOUT = 120  # seconds (folding can take a while)

# Standard amino acid codes
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


# ---------------------------------------------------------------------------
# Sequence Utilities
# ---------------------------------------------------------------------------

def read_fasta(fasta_path: str) -> tuple[str, str]:
    """Read a FASTA file and return (header, sequence)."""
    header = ""
    sequence_lines: list[str] = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header and sequence_lines:
                    break  # Only read first sequence
                header = line[1:].strip()
            elif line:
                sequence_lines.append(line.upper())

    sequence = "".join(sequence_lines)
    return header, sequence


def validate_sequence(sequence: str) -> tuple[bool, str]:
    """Validate an amino acid sequence.

    Returns:
        (is_valid, message)
    """
    sequence = sequence.upper().replace(" ", "").replace("\n", "")

    if not sequence:
        return False, "Empty sequence"

    if len(sequence) < 10:
        return False, f"Sequence too short ({len(sequence)} aa). Minimum: 10"

    if len(sequence) > 2500:
        return False, (
            f"Sequence too long ({len(sequence)} aa) for API prediction. "
            "Use 'generate-af2-script' for HPC-based prediction."
        )

    invalid_chars = set(sequence) - VALID_AA
    if invalid_chars:
        return False, f"Invalid amino acid characters: {invalid_chars}"

    return True, f"Valid sequence: {len(sequence)} residues"


# ---------------------------------------------------------------------------
# ESMFold Prediction
# ---------------------------------------------------------------------------

def predict_esmfold(sequence: str, output_path: str) -> dict[str, Any]:
    """Predict structure using ESMFold API.

    Args:
        sequence: Amino acid sequence (single letter codes).
        output_path: Path to save the PDB file.

    Returns:
        Dict with prediction metadata.
    """
    sequence = sequence.upper().replace(" ", "").replace("\n", "")

    # Validate
    valid, msg = validate_sequence(sequence)
    if not valid:
        raise ValueError(msg)

    print(f"[INFO] Predicting structure for {len(sequence)} residues via ESMFold...")

    # Call ESMFold API
    req = urllib.request.Request(
        ESMFOLD_API_URL,
        data=sequence.encode("utf-8"),
        headers={
            "Content-Type": "text/plain",
            "User-Agent": "DrugDiscoveryWizard/1.0",
        },
    )

    try:
        with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as resp:
            pdb_text = resp.read().decode("utf-8")
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"ESMFold API returned HTTP {e.code}: {body}") from e
    except urllib.error.URLError as e:
        raise RuntimeError(
            f"Cannot reach ESMFold API: {e.reason}. "
            "Use 'generate-af2-script' for offline prediction."
        ) from e

    # Save PDB
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(pdb_text)

    # Extract pLDDT from B-factor column
    plddt_values = _extract_plddt(pdb_text)
    mean_plddt = sum(plddt_values) / len(plddt_values) if plddt_values else 0.0

    # Identify low-confidence regions (stretches of residues with pLDDT < 50)
    low_conf_regions = _find_low_confidence_regions(plddt_values)

    # Determine refinement decision
    refinement = _assess_refinement(mean_plddt, plddt_values, low_conf_regions)

    result = {
        "method": "ESMFold",
        "sequence_length": len(sequence),
        "output_pdb": output_path,
        "mean_plddt": round(mean_plddt, 2),
        "median_plddt": round(sorted(plddt_values)[len(plddt_values) // 2], 2) if plddt_values else 0.0,
        "min_plddt": round(min(plddt_values), 2) if plddt_values else 0.0,
        "max_plddt": round(max(plddt_values), 2) if plddt_values else 0.0,
        "residues_above_90": sum(1 for p in plddt_values if p >= 90),
        "residues_above_70": sum(1 for p in plddt_values if p >= 70),
        "residues_50_to_70": sum(1 for p in plddt_values if 50 <= p < 70),
        "residues_below_50": sum(1 for p in plddt_values if p < 50),
        "low_confidence_regions": low_conf_regions,
        "confidence_assessment": _assess_confidence(mean_plddt),
        "refinement": refinement,
        "pdb_size_bytes": len(pdb_text),
    }

    # Auto-generate quality report markdown alongside PDB
    report_path = Path(output_path).with_name(
        Path(output_path).stem + "_quality_report.md"
    )
    report_md = _generate_quality_report(result, sequence)
    report_path.write_text(report_md)
    result["quality_report"] = str(report_path)

    # Print prominent quality summary to console
    _print_quality_summary(result)

    return result


def _extract_plddt(pdb_text: str) -> list[float]:
    """Extract per-residue pLDDT from B-factor column of PDB ATOM records."""
    plddt_by_residue: dict[int, float] = {}

    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                res_num = int(line[22:26].strip())
                bfactor = float(line[60:66].strip())
                plddt_by_residue[res_num] = bfactor
            except (ValueError, IndexError):
                continue

    return list(plddt_by_residue.values())


def _assess_confidence(mean_plddt: float) -> str:
    """Assess overall confidence of the prediction."""
    if mean_plddt >= 90:
        return "Very high confidence — suitable for docking"
    elif mean_plddt >= 70:
        return "High confidence — usable for docking, check binding site region"
    elif mean_plddt >= 50:
        return "Low confidence — use with caution, consider experimental structure"
    else:
        return "Very low confidence — not recommended for docking"


def _assess_refinement(
    mean_plddt: float,
    plddt_values: list[float],
    low_conf_regions: list[dict],
) -> dict:
    """Determine whether the structure needs refinement before docking.

    Returns a dict with:
        verdict: "PROCEED" | "REFINE" | "STOP"
        reason: Human-readable explanation
        actions: List of recommended refinement actions
    """
    total = len(plddt_values)
    if total == 0:
        return {"verdict": "STOP", "reason": "No residues found", "actions": []}

    below_50 = sum(1 for p in plddt_values if p < 50)
    pct_below_50 = below_50 / total * 100
    above_70 = sum(1 for p in plddt_values if p >= 70)
    pct_above_70 = above_70 / total * 100

    actions: list[str] = []

    # --- STOP: Structure is too poor to use ---
    if mean_plddt < 50 or pct_below_50 > 50:
        return {
            "verdict": "STOP",
            "reason": (
                f"Structure quality is too low (mean pLDDT={mean_plddt:.1f}, "
                f"{pct_below_50:.0f}% residues below 50). "
                "Do NOT proceed to docking."
            ),
            "actions": [
                "Re-predict with AlphaFold2 (MSA-based): "
                "generate-af2-script --method colabfold",
                "Search PDB for homologous experimental structures",
                "If no structure available, consider ligand-based virtual screening instead",
            ],
        }

    # --- REFINE: Structure needs work before docking ---
    needs_refinement = False

    # Check for significant disordered loops
    if low_conf_regions:
        long_loops = [r for r in low_conf_regions if r["length"] >= 5]
        if long_loops:
            needs_refinement = True
            actions.append(
                f"Remove or truncate {len(long_loops)} disordered loop(s): "
                + ", ".join(f"res {r['start_residue']}-{r['end_residue']}" for r in long_loops)
            )

    # Check for moderate overall quality
    if 50 <= mean_plddt < 70:
        needs_refinement = True
        actions.append(
            "Run energy minimization before docking: "
            "gmxapi-setup --pdb predicted.pdb --production-ns 0 (EM only)"
        )
        actions.append(
            "Consider re-prediction with AlphaFold2 for better accuracy: "
            "generate-af2-script --method colabfold"
        )

    # Check if many residues are in the 50-70 zone
    pct_50_70 = sum(1 for p in plddt_values if 50 <= p < 70) / total * 100
    if pct_50_70 > 30:
        needs_refinement = True
        actions.append(
            f"{pct_50_70:.0f}% of residues are in the 50-70 pLDDT range. "
            "Run a short MD equilibration (1-5 ns) to relax the structure"
        )

    # ESMFold-specific: single-sequence models miss co-evolutionary info
    if mean_plddt < 80:
        actions.append(
            "ESMFold uses single-sequence prediction (no MSA). "
            "AlphaFold2 with MSA may significantly improve quality"
        )

    if needs_refinement:
        # Always recommend energy minimization for refinement cases
        if not any("energy minimization" in a.lower() for a in actions):
            actions.insert(0,
                "Run energy minimization to fix steric clashes: "
                "gmxapi-setup --pdb predicted.pdb"
            )
        return {
            "verdict": "REFINE",
            "reason": (
                f"Structure has usable regions but needs refinement "
                f"(mean pLDDT={mean_plddt:.1f}, "
                f"{len(low_conf_regions)} disordered region(s))"
            ),
            "actions": actions,
        }

    # --- PROCEED: Structure is good enough ---
    # Still might have minor recommendations
    if pct_below_50 > 0:
        actions.append(
            f"Minor: {below_50} residues ({pct_below_50:.0f}%) have pLDDT < 50. "
            f"Avoid placing binding site in these regions"
        )
    if mean_plddt < 90:
        actions.append(
            "Optional: Run energy minimization for cleaner geometry before docking"
        )

    return {
        "verdict": "PROCEED",
        "reason": (
            f"Structure quality is sufficient for docking "
            f"(mean pLDDT={mean_plddt:.1f}, "
            f"{pct_above_70:.0f}% above 70)"
        ),
        "actions": actions,
    }

def _find_low_confidence_regions(plddt_values: list[float], threshold: float = 50.0) -> list[dict]:
    """Identify contiguous stretches of low-confidence residues."""
    regions: list[dict] = []
    in_region = False
    start = 0

    for i, p in enumerate(plddt_values):
        if p < threshold:
            if not in_region:
                start = i + 1  # 1-indexed residue number
                in_region = True
        else:
            if in_region:
                regions.append({
                    "start_residue": start,
                    "end_residue": i,  # 1-indexed, last low-conf residue
                    "length": i - start + 1,
                    "mean_plddt": round(sum(plddt_values[start-1:i]) / (i - start + 1), 1),
                })
                in_region = False

    if in_region:
        regions.append({
            "start_residue": start,
            "end_residue": len(plddt_values),
            "length": len(plddt_values) - start + 1,
            "mean_plddt": round(sum(plddt_values[start-1:]) / (len(plddt_values) - start + 1), 1),
        })

    return regions


def _generate_quality_report(result: dict, sequence: str) -> str:
    """Generate a detailed quality report as markdown."""
    total = result["sequence_length"]
    above_90 = result["residues_above_90"]
    above_70 = result["residues_above_70"]
    r50_70 = result["residues_50_to_70"]
    below_50 = result["residues_below_50"]
    pct_above_70 = (above_70 / total * 100) if total else 0
    pct_below_50 = (below_50 / total * 100) if total else 0

    # Build pLDDT visual bar (10 segments)
    bar_chars = []
    segment_size = max(1, total // 40)
    for i in range(0, total, segment_size):
        chunk = [result.get("mean_plddt", 50)]  # fallback
        # simplified: just use overall confidence for the bar
    # Use a text-based confidence bar
    bar = _plddt_confidence_bar(result["mean_plddt"])

    lines = [
        "# Structure Prediction Quality Report",
        "",
        f"**Method**: {result['method']}",
        f"**Sequence length**: {total} residues",
        f"**PDB file**: `{result['output_pdb']}`",
        "",
        "---",
        "",
        "## Overall Confidence",
        "",
        f"```",
        f"Mean pLDDT:   {result['mean_plddt']:6.1f}  {bar}",
        f"Median pLDDT: {result['median_plddt']:6.1f}",
        f"Min pLDDT:    {result['min_plddt']:6.1f}",
        f"Max pLDDT:    {result['max_plddt']:6.1f}",
        f"```",
        "",
        f"**Assessment**: {result['confidence_assessment']}",
        "",
        "## pLDDT Distribution",
        "",
        "| Confidence Level | pLDDT Range | Residues | % | Interpretation |",
        "|---|---|---:|---:|---|",
        f"| Very High | ≥90 | {above_90} | {above_90/total*100:.0f}% | Excellent — highly reliable |",
        f"| High | 70-89 | {above_70 - above_90} | {(above_70 - above_90)/total*100:.0f}% | Good — suitable for docking |",
        f"| Low | 50-69 | {r50_70} | {r50_70/total*100:.0f}% | Caution — possible disorder |",
        f"| Very Low | <50 | {below_50} | {below_50/total*100:.0f}% | Unreliable — likely disordered |",
        "",
    ]

    # Warnings
    if result["mean_plddt"] < 50:
        lines.extend([
            "> ⛔ **CRITICAL**: Mean pLDDT is below 50. This structure is unreliable.",
            "> **Do NOT use this for docking.** Consider:",
            "> - Using AlphaFold2 with MSA: `generate-af2-script --method colabfold`",
            "> - Searching for a homologous experimental structure in PDB",
            "",
        ])
    elif result["mean_plddt"] < 70:
        lines.extend([
            "> ⚠️ **WARNING**: Mean pLDDT is below 70. Use this structure with caution.",
            "> - Only dock against regions with pLDDT > 70",
            "> - Consider AlphaFold2 for higher accuracy: `generate-af2-script`",
            "",
        ])
    else:
        lines.extend([
            "> ✅ **GOOD**: Mean pLDDT is above 70. This structure is suitable for",
            "> docking in the high-confidence regions.",
            "",
        ])

    # Low-confidence regions
    low_regions = result.get("low_confidence_regions", [])
    if low_regions:
        lines.extend([
            "## Low-Confidence Regions (pLDDT < 50)",
            "",
            "These regions are likely disordered or incorrectly folded.",
            "**Do NOT place your docking binding site in these regions.**",
            "",
            "| Region | Residues | Length | Mean pLDDT |",
            "|---|---|---:|---:|",
        ])
        for r in low_regions:
            lines.append(
                f"| {r['start_residue']}-{r['end_residue']} | "
                f"res {r['start_residue']}–{r['end_residue']} | "
                f"{r['length']} | {r['mean_plddt']:.1f} |"
            )
        lines.append("")

    # Refinement Decision
    refinement = result.get("refinement", {})
    verdict = refinement.get("verdict", "UNKNOWN")
    reason = refinement.get("reason", "")
    actions = refinement.get("actions", [])

    verdict_label = {
        "PROCEED": "[PROCEED] -- Ready for docking",
        "REFINE": "[REFINE] -- Refinement needed before docking",
        "STOP": "[STOP] -- Do NOT proceed to docking",
    }.get(verdict, verdict)

    lines.extend([
        "## Refinement Decision",
        "",
        f"### {verdict_label}",
        "",
        f"{reason}",
        "",
    ])

    if actions:
        lines.append("**Recommended actions:**")
        lines.append("")
        for i, action in enumerate(actions, 1):
            lines.append(f"{i}. {action}")
        lines.append("")

    if verdict == "PROCEED":
        lines.extend([
            "> Structure is ready for the next step in the pipeline.",
            "> Proceed to **Step 1b** (receptor preparation):",
            "> `python3 scripts/drug_discovery_api.py prepare-receptor "
            "--input-file predicted.pdb --output receptor.json`",
            "",
        ])
    elif verdict == "REFINE":
        lines.extend([
            "> Complete the recommended refinement actions above before proceeding.",
            "> After refinement, re-assess quality:",
            "> `python3 scripts/drug_discovery_api.py assess-structure "
            "--pdb refined.pdb --output quality.json`",
            "",
        ])
    elif verdict == "STOP":
        lines.extend([
            "> This structure is not suitable for structure-based drug discovery.",
            "> Follow the recommended actions above to obtain a better structure.",
            "",
        ])

    # Recommendations
    lines.extend([
        "## Recommendations for Downstream Use",
        "",
        "### For Molecular Docking (Step 4)",
        f"- Use regions with pLDDT > 70 ({pct_above_70:.0f}% of structure)",
        "- Define binding site only in high-confidence regions",
        "- Avoid flexible loops (pLDDT < 50) as docking targets",
        "",
        "### For Molecular Dynamics (Step 6)",
        "- Low-confidence regions may need extended equilibration",
        "- Consider restraining poorly predicted loops during NVT/NPT",
        "",
        "### If Quality Is Insufficient",
        "- Run full AlphaFold2 with MSA on HPC:",
        "  `python3 scripts/drug_discovery_api.py generate-af2-script --fasta seq.fasta --method colabfold --output submit_af2.sh`",
        "- Search for experimental structures of homologs in PDB",
        "",
        "---",
        "",
        "*pLDDT (predicted Local Distance Difference Test) measures per-residue*",
        "*prediction confidence on a 0-100 scale. Values above 70 are generally*",
        "*considered reliable for structural biology applications.*",
    ])

    return "\n".join(lines) + "\n"


def _plddt_confidence_bar(mean_plddt: float) -> str:
    """Generate a text-based confidence bar (ASCII-safe for Windows)."""
    filled = int(mean_plddt / 5)  # 20 chars for 0-100
    empty = 20 - filled
    if mean_plddt >= 90:
        label = "### EXCELLENT"
    elif mean_plddt >= 70:
        label = "##- GOOD"
    elif mean_plddt >= 50:
        label = "#-- CAUTION"
    else:
        label = "--- POOR"
    return f"[{'#' * filled}{'-' * empty}] {label}"


def _print_quality_summary(result: dict) -> None:
    """Print an unmissable quality summary to console."""
    mean = result["mean_plddt"]
    total = result["sequence_length"]
    below_50 = result["residues_below_50"]
    above_70 = result["residues_above_70"]

    print("")
    print("=" * 60)
    print("  STRUCTURE PREDICTION QUALITY REPORT")
    print("=" * 60)
    print(f"  PDB file:    {result['output_pdb']}")
    print(f"  Quality report: {result.get('quality_report', 'N/A')}")
    print(f"  Residues:    {total}")
    print(f"  Mean pLDDT:  {mean:.1f} / 100")
    print(f"  {_plddt_confidence_bar(mean)}")
    print("")
    print(f"  High confidence (>70):  {above_70:>4d} residues ({above_70/total*100:.0f}%)")
    print(f"  Low confidence (<50):   {below_50:>4d} residues ({below_50/total*100:.0f}%)")

    low_regions = result.get("low_confidence_regions", [])
    if low_regions:
        print("")
        print(f"  [WARNING] {len(low_regions)} disordered region(s) detected:")
        for r in low_regions:
            print(f"     - Residues {r['start_residue']}-{r['end_residue']} "
                  f"(pLDDT={r['mean_plddt']:.1f})")
        print("     Do NOT place binding site in these regions!")

    # Refinement Decision
    refinement = result.get("refinement", {})
    verdict = refinement.get("verdict", "")
    if verdict:
        print("")
        print("-" * 60)
        verdict_display = {
            "PROCEED": "  VERDICT: [PROCEED] Ready for docking",
            "REFINE": "  VERDICT: [REFINE] Refinement needed before docking",
            "STOP": "  VERDICT: [STOP] Do NOT proceed to docking",
        }.get(verdict, f"  VERDICT: {verdict}")
        print(verdict_display)
        print(f"  {refinement.get('reason', '')}")

        actions = refinement.get("actions", [])
        if actions:
            print("")
            print("  Recommended actions:")
            for i, action in enumerate(actions, 1):
                print(f"    {i}. {action}")

        if verdict == "PROCEED":
            print("")
            print("  Next step: prepare-receptor --input-file predicted.pdb --output receptor.json")
        elif verdict == "REFINE":
            print("")
            print("  After refinement, re-assess: assess-structure --pdb refined.pdb --output quality.json")

    print("")
    print("=" * 60)
    print("")

def assess_structure(pdb_path: str) -> dict[str, Any]:
    """Assess prediction quality from a PDB file's B-factor column."""
    with open(pdb_path) as f:
        pdb_text = f.read()

    plddt_values = _extract_plddt(pdb_text)

    if not plddt_values:
        return {
            "file": pdb_path,
            "error": "No CA atoms found — cannot assess pLDDT",
        }

    mean_plddt = sum(plddt_values) / len(plddt_values)
    low_conf_regions = _find_low_confidence_regions(plddt_values)

    result = {
        "file": pdb_path,
        "output_pdb": pdb_path,
        "method": "Assessed from PDB",
        "total_residues": len(plddt_values),
        "sequence_length": len(plddt_values),
        "mean_plddt": round(mean_plddt, 2),
        "median_plddt": round(sorted(plddt_values)[len(plddt_values) // 2], 2),
        "min_plddt": round(min(plddt_values), 2),
        "max_plddt": round(max(plddt_values), 2),
        "residues_above_90": sum(1 for p in plddt_values if p >= 90),
        "residues_above_70": sum(1 for p in plddt_values if p >= 70),
        "residues_50_to_70": sum(1 for p in plddt_values if 50 <= p < 70),
        "residues_below_50": sum(1 for p in plddt_values if p < 50),
        "low_confidence_regions": low_conf_regions,
        "confidence_assessment": _assess_confidence(mean_plddt),
        "refinement": _assess_refinement(mean_plddt, plddt_values, low_conf_regions),
        "docking_recommendation": (
            "Use high-confidence regions (pLDDT > 70) for binding site definition. "
            "Avoid flexible loops (pLDDT < 50) as docking targets."
        ),
    }

    # Print quality summary to console
    _print_quality_summary(result)

    return result


# ---------------------------------------------------------------------------
# AF2 / ColabFold HPC Script Generation
# ---------------------------------------------------------------------------

def generate_af2_hpc_script(
    fasta_path: str,
    output_path: str,
    method: str = "colabfold",
    gpu: str = "a100",
    ngpu: int = 1,
    mem: str = "64G",
    time: str = "12:00:00",
    db_path: str = "/data/colabfold_dbs",
    conda_env: str = "colabfold",
) -> str:
    """Generate a SLURM HPC script for AF2/ColabFold prediction."""

    if method == "colabfold":
        predict_cmd = (
            f"colabfold_batch {fasta_path} output/ "
            f"--num-recycle 3 --amber --use-gpu-relax"
        )
        module_load = "# module load colabfold  # Uncomment if using modules"
    elif method == "alphafold2":
        predict_cmd = (
            f"python3 run_alphafold.py "
            f"--fasta_paths={fasta_path} "
            f"--output_dir=output/ "
            f"--data_dir={db_path} "
            f"--model_preset=monomer "
            f"--max_template_date=2024-01-01 "
            f"--use_gpu_relax=true"
        )
        module_load = "# module load alphafold/2.3  # Uncomment if using modules"
    else:
        raise ValueError(f"Unknown method '{method}'. Choose 'colabfold' or 'alphafold2'.")

    script = f"""#!/bin/bash
#SBATCH --job-name=af2_predict
#SBATCH --partition=gpu
#SBATCH --gres=gpu:{gpu}:{ngpu}
#SBATCH --cpus-per-task=8
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --output=af2_%j.out
#SBATCH --error=af2_%j.err

# ============================================================
# AlphaFold2 / ColabFold Structure Prediction
# Generated by Drug Discovery Wizard
# Method: {method}
# ============================================================

set -euo pipefail

{module_load}

# Activate environment
conda activate {conda_env}

echo "=== Structure Prediction ==="
echo "Method: {method}"
echo "Input: {fasta_path}"
echo "GPU: {gpu} x {ngpu}"
echo "Start: $(date)"

mkdir -p output/

# Run prediction
{predict_cmd}

echo "=== Prediction Complete ==="
echo "End: $(date)"
echo "Output files in: output/"

# Assess quality
echo ""
echo "=== Quality Assessment ==="
for pdb in output/*.pdb; do
    if [ -f "$pdb" ]; then
        echo "Structure: $pdb"
        # Extract mean pLDDT from B-factor column
        mean_plddt=$(grep "^ATOM" "$pdb" | awk '{{sum+=$11; n++}} END {{if(n>0) printf "%.1f", sum/n; else print "N/A"}}')
        echo "  Mean pLDDT: $mean_plddt"
    fi
done

echo ""
echo "Next steps:"
echo "  1. Check pLDDT scores (>70 for reliable regions)"
echo "  2. Prepare receptor: python3 scripts/drug_discovery_api.py prepare-receptor --input-file output/best.pdb --output receptor.json"
echo "  3. Define binding site from high-confidence regions"
"""

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="\n") as f:
        f.write(script)

    return script


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="De novo protein structure prediction (ESMFold / AF2)"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # --- predict ---
    p_pred = sub.add_parser("predict", help="Predict structure via ESMFold API")
    seq_group = p_pred.add_mutually_exclusive_group(required=True)
    seq_group.add_argument("--sequence", help="Amino acid sequence")
    seq_group.add_argument("--fasta", help="FASTA file path")
    p_pred.add_argument("--output", required=True, help="Output PDB path")
    p_pred.add_argument("--json-output", help="Optional JSON metadata output")

    # --- generate-af2-script ---
    p_af2 = sub.add_parser("generate-af2-script",
                            help="Generate AF2/ColabFold HPC prediction script")
    p_af2.add_argument("--fasta", required=True, help="Input FASTA file")
    p_af2.add_argument("--method", default="colabfold",
                        choices=["colabfold", "alphafold2"])
    p_af2.add_argument("--gpu", default="a100")
    p_af2.add_argument("--ngpu", type=int, default=1)
    p_af2.add_argument("--mem", default="64G")
    p_af2.add_argument("--time", default="12:00:00")
    p_af2.add_argument("--db-path", default="/data/colabfold_dbs")
    p_af2.add_argument("--conda-env", default="colabfold")
    p_af2.add_argument("--output", required=True, help="Output script path")

    # --- assess ---
    p_assess = sub.add_parser("assess", help="Assess prediction quality from PDB")
    p_assess.add_argument("--pdb", required=True, help="PDB file to assess")
    p_assess.add_argument("--output", required=True, help="Output JSON report")

    args = parser.parse_args()

    if args.command == "predict":
        if args.fasta:
            header, sequence = read_fasta(args.fasta)
            print(f"[INFO] Read sequence from FASTA: {header}")
        else:
            sequence = args.sequence

        result = predict_esmfold(sequence, args.output)

        if args.json_output:
            with open(args.json_output, "w") as f:
                json.dump(result, f, indent=2)

        # Also write JSON alongside PDB
        json_path = Path(args.output).with_suffix(".json")
        with open(json_path, "w") as f:
            json.dump(result, f, indent=2)
        print(f"[OK] Metadata saved to {json_path}")

    elif args.command == "generate-af2-script":
        generate_af2_hpc_script(
            fasta_path=args.fasta,
            output_path=args.output,
            method=args.method,
            gpu=args.gpu,
            ngpu=args.ngpu,
            mem=args.mem,
            time=args.time,
            db_path=args.db_path,
            conda_env=args.conda_env,
        )
        print(f"[OK] AF2 HPC script saved to {args.output}")
        print(f"     Submit with: sbatch {args.output}")

    elif args.command == "assess":
        result = assess_structure(args.pdb)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            json.dump(result, f, indent=2)
        print(f"[OK] Quality report saved to {args.output}")
        print(f"     Mean pLDDT: {result.get('mean_plddt', 'N/A')}")
        print(f"     Assessment: {result.get('confidence_assessment', 'N/A')}")


if __name__ == "__main__":
    main()
