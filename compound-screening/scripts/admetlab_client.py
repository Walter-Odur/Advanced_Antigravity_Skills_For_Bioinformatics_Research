"""ADMETlab 3.0 REST API Client.

Queries the ADMETlab 3.0 web service for comprehensive ADMET predictions
(119 endpoints). Falls back to local RDKit-based filtering when the API
is unreachable.

Usage:
    python admetlab_client.py predict --smiles "CCO" "c1ccccc1" --output results.json
    python admetlab_client.py filter  --smiles "CCO" --output filtered.json
    python admetlab_client.py summarize --input results.json --output summary.md
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ADMETLAB_API_URL = "https://admetlab3.scbdd.com/server/api/aio"
RATE_LIMIT_DELAY = 0.25  # 4 req/sec max (conservative vs 5/sec limit)
MAX_BATCH = 100  # SMILES per batch request
REQUEST_TIMEOUT = 60  # seconds

# ---------------------------------------------------------------------------
# Pass/fail thresholds — TWO-TIER STRATEGY
#
# Tier 1 (RELAXED): Post-docking. Keeps more hits to seed REINVENT.
#   Wider physicochemical ranges, skips CYP/metabolic checks, allows
#   minor structural alerts. Goal: maximize diversity of training seeds.
#
# Tier 2 (STRICT): Post-REINVENT. Rigorous filtering for compounds
#   proceeding to MD and downstream analysis.
# ---------------------------------------------------------------------------

THRESHOLDS_STRICT: dict[str, dict[str, Any]] = {
    # --- Physicochemical ---
    "MW": {"max": 500, "desc": "Molecular weight (Lipinski: <=500)"},
    "LogP": {"max": 5.0, "desc": "Octanol-water partition (Lipinski: <=5)"},
    "nHD": {"max": 5, "desc": "H-bond donors (Lipinski: <=5)"},
    "nHA": {"max": 10, "desc": "H-bond acceptors (Lipinski: <=10)"},
    "TPSA": {"max": 140, "desc": "Topological polar surface area (<=140 for oral drugs)"},
    "LogS": {"min": -6.0, "desc": "Aqueous solubility (>-6 = acceptable)"},
    # --- Absorption ---
    "Caco-2": {"min": -5.15, "desc": "Caco-2 permeability (>-5.15 = good)"},
    "HIA": {"category": "positive", "desc": "Human intestinal absorption (positive = absorbed)"},
    "Pgp-inh": {"category": "negative", "desc": "P-gp inhibitor (negative = not inhibitor)"},
    "Pgp-sub": {"category": "negative", "desc": "P-gp substrate (negative = not substrate)"},
    "F20": {"category": "positive", "desc": "Oral bioavailability >=20% (positive = yes)"},
    "F30": {"category": "positive", "desc": "Oral bioavailability >=30% (positive = yes)"},
    # --- Distribution ---
    "BBB": {"category": "positive", "desc": "Blood-brain barrier penetration"},
    "PPB": {"max": 95.0, "desc": "Plasma protein binding (<=95% = acceptable)"},
    "VDss": {"min": 0.04, "max": 20.0, "desc": "Volume of distribution (0.04-20 L/kg)"},
    # --- Metabolism ---
    "CYP1A2-inh": {"category": "negative", "desc": "CYP1A2 inhibitor (negative = safe)"},
    "CYP2C9-inh": {"category": "negative", "desc": "CYP2C9 inhibitor (negative = safe)"},
    "CYP2C19-inh": {"category": "negative", "desc": "CYP2C19 inhibitor (negative = safe)"},
    "CYP2D6-inh": {"category": "negative", "desc": "CYP2D6 inhibitor (negative = safe)"},
    "CYP3A4-inh": {"category": "negative", "desc": "CYP3A4 inhibitor (negative = safe)"},
    "CYP1A2-sub": {"category": "negative", "desc": "CYP1A2 substrate"},
    "CYP2C9-sub": {"category": "negative", "desc": "CYP2C9 substrate"},
    "CYP2C19-sub": {"category": "negative", "desc": "CYP2C19 substrate"},
    "CYP2D6-sub": {"category": "negative", "desc": "CYP2D6 substrate"},
    "CYP3A4-sub": {"category": "negative", "desc": "CYP3A4 substrate"},
    # --- Excretion ---
    "CL": {"min": 1.0, "desc": "Clearance (>1 mL/min/kg = acceptable)"},
    "T12": {"min": 1.0, "desc": "Half-life (>1 hour = acceptable)"},
    # --- Toxicity ---
    "hERG": {"category": "negative", "desc": "hERG inhibition (negative = safe)"},
    "AMES": {"category": "negative", "desc": "Ames mutagenicity (negative = safe)"},
    "DILI": {"category": "negative", "desc": "Drug-induced liver injury (negative = safe)"},
    "Carcinogenicity": {"category": "negative", "desc": "Carcinogenicity (negative = safe)"},
    "SkinSen": {"category": "negative", "desc": "Skin sensitization (negative = safe)"},
    "EI": {"category": "negative", "desc": "Eye irritation (negative = safe)"},
    "RespiratoryTox": {"category": "negative", "desc": "Respiratory toxicity (negative = safe)"},
    # --- Medicinal Chemistry ---
    "Lipinski": {"category": "accepted", "desc": "Lipinski Rule of Five"},
    "Pfizer": {"category": "accepted", "desc": "Pfizer 3/75 Rule"},
    "GoldenTriangle": {"category": "accepted", "desc": "Golden Triangle (MW vs LogD)"},
    "PAINS": {"max": 0, "desc": "PAINS alerts (0 = clean)"},
    "Brenk": {"max": 0, "desc": "Brenk structural alerts (0 = clean)"},
}

THRESHOLDS_RELAXED: dict[str, dict[str, Any]] = {
    # --- Physicochemical (widened for seed diversity) ---
    "MW": {"max": 700, "desc": "Molecular weight (relaxed: <=700 for seed diversity)"},
    "LogP": {"max": 7.0, "desc": "Octanol-water partition (relaxed: <=7)"},
    "nHD": {"max": 7, "desc": "H-bond donors (relaxed: <=7)"},
    "nHA": {"max": 15, "desc": "H-bond acceptors (relaxed: <=15)"},
    "TPSA": {"max": 200, "desc": "Topological polar surface area (relaxed: <=200)"},
    "LogS": {"min": -8.0, "desc": "Aqueous solubility (relaxed: >-8)"},
    # --- Absorption (keep key absorption, loosen permeability) ---
    "HIA": {"category": "positive", "desc": "Human intestinal absorption (positive = absorbed)"},
    "F20": {"category": "positive", "desc": "Oral bioavailability >=20% (positive = yes)"},
    # NOTE: Caco-2, Pgp-inh, Pgp-sub, F30 — SKIPPED in relaxed mode
    # --- Distribution (skip BBB — not relevant for all targets) ---
    "PPB": {"max": 99.0, "desc": "Plasma protein binding (relaxed: <=99%)"},
    # NOTE: BBB, VDss — SKIPPED in relaxed mode
    # --- Metabolism — SKIPPED entirely in relaxed mode ---
    # CYP interactions are tunable by REINVENT. No point filtering seeds.
    # --- Excretion — SKIPPED in relaxed mode ---
    # Clearance and half-life are optimizable.
    # --- Toxicity (keep only critical safety) ---
    "AMES": {"category": "negative", "desc": "Ames mutagenicity (negative = safe)"},
    "Carcinogenicity": {"category": "negative", "desc": "Carcinogenicity (negative = safe)"},
    # NOTE: hERG, DILI, SkinSen, EI, RespiratoryTox — SKIPPED
    # hERG and DILI are tunable by REINVENT scoring functions
    # --- Medicinal Chemistry (loosened) ---
    "PAINS": {"max": 1, "desc": "PAINS alerts (relaxed: <=1 for seed diversity)"},
    "Brenk": {"max": 2, "desc": "Brenk structural alerts (relaxed: <=2)"},
    # NOTE: Lipinski, Pfizer, GoldenTriangle — SKIPPED (already widened individual props)
}

# Backward compatibility alias
THRESHOLDS = THRESHOLDS_STRICT


def get_thresholds(strictness: str = "strict") -> dict[str, dict[str, Any]]:
    """Get threshold profile by strictness level.

    Args:
        strictness: "strict" (default, for post-REINVENT) or
                    "relaxed" (for post-docking REINVENT seeding)

    Returns:
        Dict of property thresholds.
    """
    profiles = {
        "strict": THRESHOLDS_STRICT,
        "relaxed": THRESHOLDS_RELAXED,
    }
    if strictness not in profiles:
        raise ValueError(
            f"Unknown strictness '{strictness}'. Choose: {list(profiles.keys())}"
        )
    return profiles[strictness]


# ---------------------------------------------------------------------------
# API Client
# ---------------------------------------------------------------------------

def _post_json(url: str, payload: dict, timeout: int = REQUEST_TIMEOUT) -> dict:
    """POST JSON to ADMETlab and return parsed response."""
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url,
        data=data,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "DrugDiscoveryWizard/1.0",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", errors="replace")
        raise RuntimeError(
            f"ADMETlab API returned HTTP {e.code}: {body}"
        ) from e
    except urllib.error.URLError as e:
        raise RuntimeError(f"Cannot reach ADMETlab API: {e.reason}") from e


def predict_admet(smiles_list: list[str]) -> list[dict]:
    """Call ADMETlab 3.0 API for ADMET predictions.

    Args:
        smiles_list: List of SMILES strings (max 100 per batch).

    Returns:
        List of dicts, one per molecule, with all ADMET properties.
    """
    all_results: list[dict] = []

    for batch_start in range(0, len(smiles_list), MAX_BATCH):
        batch = smiles_list[batch_start: batch_start + MAX_BATCH]
        payload = {"smiles": "\n".join(batch)}

        resp = _post_json(ADMETLAB_API_URL, payload)

        # Parse response — ADMETlab returns various formats
        if isinstance(resp, dict):
            if "data" in resp:
                results = resp["data"]
            elif "results" in resp:
                results = resp["results"]
            else:
                results = [resp]
        elif isinstance(resp, list):
            results = resp
        else:
            results = [{"raw_response": resp}]

        all_results.extend(results if isinstance(results, list) else [results])

        # Rate limit between batches
        if batch_start + MAX_BATCH < len(smiles_list):
            time.sleep(RATE_LIMIT_DELAY)

    return all_results


def apply_thresholds(
    predictions: list[dict],
    smiles_list: list[str],
    names: list[str] | None = None,
    strictness: str = "strict",
) -> dict:
    """Apply pass/fail thresholds to ADMETlab predictions.

    Args:
        predictions: List of prediction dicts from ADMETlab
        smiles_list: Corresponding SMILES strings
        names: Optional compound names
        strictness: "strict" (post-REINVENT) or "relaxed" (post-docking seeds)

    Returns structured results with per-property pass/fail status.
    """
    active_thresholds = get_thresholds(strictness)
    names = names or [f"compound_{i}" for i in range(len(smiles_list))]
    compounds = []

    for i, pred in enumerate(predictions):
        smi = smiles_list[i] if i < len(smiles_list) else "unknown"
        name = names[i] if i < len(names) else f"compound_{i}"

        property_results: dict[str, dict] = {}
        flags: list[str] = []
        passes = 0
        total_checked = 0

        for prop_key, thresh in active_thresholds.items():
            val = pred.get(prop_key)
            if val is None:
                continue

            total_checked += 1
            passed = True
            reason = ""

            if "category" in thresh:
                expected = thresh["category"]
                actual = str(val).lower().strip()
                if expected == "negative":
                    passed = actual in ("negative", "non-inhibitor",
                                        "non-substrate", "non-toxic",
                                        "no", "0", "false", "safe",
                                        "rejected", "-")
                elif expected == "positive":
                    passed = actual in ("positive", "yes", "1", "true",
                                        "absorbed", "permeable", "accepted", "+")
                elif expected == "accepted":
                    passed = actual in ("accepted", "yes", "pass", "1", "true", "+")
                reason = f"expected={expected}, got={val}"
            else:
                try:
                    num_val = float(val)
                except (ValueError, TypeError):
                    continue
                if "min" in thresh and num_val < thresh["min"]:
                    passed = False
                    reason = f"{num_val} < {thresh['min']}"
                if "max" in thresh and num_val > thresh["max"]:
                    passed = False
                    reason = f"{num_val} > {thresh['max']}"

            if passed:
                passes += 1
            else:
                flags.append(f"{prop_key}: {reason}")

            property_results[prop_key] = {
                "value": val,
                "pass": passed,
                "description": thresh.get("desc", ""),
            }

        compounds.append({
            "name": name,
            "smiles": smi,
            "properties": property_results,
            "flags": flags,
            "checked": total_checked,
            "passed": passes,
            "failed": total_checked - passes,
            "overall_pass": len(flags) == 0,
            "raw_prediction": pred,
        })

    tier_label = {
        "strict": "Tier 2 (strict) — post-REINVENT, for downstream analysis",
        "relaxed": "Tier 1 (relaxed) — post-docking, for REINVENT seeding",
    }.get(strictness, strictness)

    return {
        "source": "ADMETlab 3.0",
        "strictness": strictness,
        "tier": tier_label,
        "properties_checked": len(active_thresholds),
        "total_compounds": len(compounds),
        "total_passed": sum(1 for c in compounds if c["overall_pass"]),
        "total_failed": sum(1 for c in compounds if not c["overall_pass"]),
        "compounds": compounds,
    }


def generate_summary_markdown(results: dict) -> str:
    """Generate a ranked markdown summary table from filtered results."""
    lines = [
        "# ADMET Prediction Summary (ADMETlab 3.0)",
        "",
        f"**Total**: {results['total_compounds']} compounds | "
        f"**Passed**: {results['total_passed']} | "
        f"**Failed**: {results['total_failed']}",
        "",
        "| Rank | Name | SMILES | Checks Passed | Flags | Overall |",
        "|---|---|---|---|---|---|",
    ]

    sorted_compounds = sorted(
        results["compounds"],
        key=lambda c: (-c["passed"], c["failed"]),
    )

    for rank, comp in enumerate(sorted_compounds, 1):
        status = "PASS" if comp["overall_pass"] else "FAIL"
        flag_str = "; ".join(comp["flags"][:3])
        if len(comp["flags"]) > 3:
            flag_str += f" (+{len(comp['flags']) - 3} more)"
        smi_short = comp["smiles"][:40] + "..." if len(comp["smiles"]) > 40 else comp["smiles"]
        lines.append(
            f"| {rank} | {comp['name']} | `{smi_short}` | "
            f"{comp['passed']}/{comp['checked']} | {flag_str or '-'} | {status} |"
        )

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# RDKit Fallback
# ---------------------------------------------------------------------------

def _rdkit_fallback(smiles_list: list[str], names: list[str] | None = None) -> dict:
    """Fallback to RDKit-based ADMET when API is unreachable."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, FilterCatalog
        from rdkit.Contrib.SA_Score import sascorer
    except ImportError:
        return {
            "error": "Neither ADMETlab API nor RDKit available",
            "compounds": [],
        }

    names = names or [f"compound_{i}" for i in range(len(smiles_list))]

    # Build PAINS catalog
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)

    compounds = []
    for i, smi in enumerate(smiles_list):
        name = names[i] if i < len(names) else f"compound_{i}"
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            compounds.append({"name": name, "smiles": smi, "error": "Invalid SMILES"})
            continue

        mw = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
        entry = catalog.GetFirstMatch(mol)
        pains_alerts = []
        if entry:
            pains_alerts = [entry.GetDescription()]
        sa = sascorer.calculateScore(mol)

        compounds.append({
            "name": name,
            "smiles": smi,
            "mw": round(mw, 2),
            "logp": round(logp, 2),
            "hbd": hbd,
            "hba": hba,
            "lipinski_violations": violations,
            "lipinski_pass": violations == 0,
            "pains_pass": len(pains_alerts) == 0,
            "pains_alerts": pains_alerts,
            "sa_score": round(sa, 2),
            "sa_pass": sa <= 6.0,
            "overall_pass": violations == 0 and len(pains_alerts) == 0 and sa <= 6.0,
        })

    return {
        "source": "RDKit (offline fallback)",
        "total_compounds": len(compounds),
        "total_passed": sum(1 for c in compounds if c.get("overall_pass", False)),
        "total_failed": sum(1 for c in compounds if not c.get("overall_pass", True)),
        "compounds": compounds,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="ADMETlab 3.0 Client — comprehensive ADMET prediction"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # --- predict ---
    p_pred = sub.add_parser("predict", help="Full 119-property ADMET prediction")
    p_pred.add_argument("--smiles", nargs="+", required=True, help="SMILES strings")
    p_pred.add_argument("--names", nargs="+", help="Compound names")
    p_pred.add_argument("--output", required=True, help="Output JSON path")

    # --- filter ---
    p_filt = sub.add_parser("filter", help="Predict + apply pass/fail thresholds")
    p_filt.add_argument("--smiles", nargs="+", required=True, help="SMILES strings")
    p_filt.add_argument("--names", nargs="+", help="Compound names")
    p_filt.add_argument("--strictness", default="strict",
                         choices=["strict", "relaxed"],
                         help="'relaxed' = post-docking seeds for REINVENT, "
                              "'strict' = post-REINVENT for downstream (default)")
    p_filt.add_argument("--output", required=True, help="Output JSON path")

    # --- summarize ---
    p_sum = sub.add_parser("summarize", help="Generate markdown summary from results")
    p_sum.add_argument("--input", required=True, help="Input JSON from filter")
    p_sum.add_argument("--output", required=True, help="Output markdown path")

    args = parser.parse_args()

    if args.command == "predict":
        try:
            results = predict_admet(args.smiles)
            output = {
                "source": "ADMETlab 3.0",
                "total_compounds": len(results),
                "predictions": results,
            }
        except RuntimeError as e:
            print(f"[WARNING] API unavailable: {e}", file=sys.stderr)
            print("[INFO] Falling back to RDKit offline ADMET", file=sys.stderr)
            output = _rdkit_fallback(args.smiles, args.names)

        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            json.dump(output, f, indent=2, default=str)
        print(f"[OK] Predictions saved to {args.output}")

    elif args.command == "filter":
        try:
            predictions = predict_admet(args.smiles)
            output = apply_thresholds(
                predictions, args.smiles, args.names,
                strictness=args.strictness,
            )
        except RuntimeError as e:
            print(f"[WARNING] API unavailable: {e}", file=sys.stderr)
            print("[INFO] Falling back to RDKit offline ADMET", file=sys.stderr)
            output = _rdkit_fallback(args.smiles, args.names)

        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            json.dump(output, f, indent=2, default=str)
        print(f"[OK] Filtered results saved to {args.output}")

    elif args.command == "summarize":
        with open(args.input) as f:
            data = json.load(f)
        md = generate_summary_markdown(data)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            f.write(md)
        print(f"[OK] Summary saved to {args.output}")


if __name__ == "__main__":
    main()
