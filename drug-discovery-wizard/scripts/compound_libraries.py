"""Compound library clients for Drug Discovery Wizard.

Provides programmatic access to:
  1. COCONUT (COlleCtion of Open NatUral producTs) — 400k+ natural products
     Includes: ANPDB, CMNPD, NANPDB, AfroDB, SANCDB, ConMedNP, and more
  2. ZINC20 — purchasable drug-like compounds (via tranche downloads)

Usage via CLI:
    # Search COCONUT for natural products
    python3 scripts/drug_discovery_api.py query-coconut \
        --query "manzamine" --limit 20 --output results.json

    # Search COCONUT by SMILES (similarity)
    python3 scripts/drug_discovery_api.py query-coconut \
        --smiles "CC(=O)Oc1ccccc1C(=O)O" --limit 10 --output results.json

    # Generate ZINC20 download script for drug-like tranches
    python3 scripts/drug_discovery_api.py query-zinc \
        --subset drug-like --mw-range 250 500 --logp-range -1 5 \
        --output zinc_download.sh
"""

from __future__ import annotations

import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

COCONUT_BASE = "https://coconut.naturalproducts.net"
COCONUT_SEARCH_URL = f"{COCONUT_BASE}/api/search"
REQUEST_TIMEOUT = 30

# ZINC20 tranche download base
ZINC20_BASE = "https://zinc20.docking.org"
ZINC20_TRANCHES_URL = f"{ZINC20_BASE}/tranches/download"

# Standard drug-like property ranges for ZINC tranches
ZINC_SUBSETS = {
    "drug-like": {"mw": (250, 500), "logp": (-1, 5), "reactive": "anodyne"},
    "lead-like": {"mw": (250, 350), "logp": (-1, 3.5), "reactive": "anodyne"},
    "fragment-like": {"mw": (0, 250), "logp": (-1, 3.5), "reactive": "anodyne"},
    "all-purchasable": {"mw": (0, 600), "logp": (-4, 6), "reactive": "anodyne"},
}


# ---------------------------------------------------------------------------
# COCONUT Client
# ---------------------------------------------------------------------------

def search_coconut(
    query: str | None = None,
    smiles: str | None = None,
    limit: int = 25,
    page: int = 1,
) -> dict[str, Any]:
    """Search COCONUT natural products database.

    Args:
        query: Text search (compound name, organism, etc.)
        smiles: SMILES string for structure search
        limit: Max results per page (max 50)
        page: Page number for pagination

    Returns:
        Dict with search results including SMILES, names, identifiers.
    """
    if not query and not smiles:
        raise ValueError("Provide either --query or --smiles")

    search_term = query or smiles
    limit = min(limit, 50)

    # Build request
    payload = json.dumps({
        "query": search_term,
        "limit": limit,
        "page": page,
    }).encode("utf-8")

    req = urllib.request.Request(
        COCONUT_SEARCH_URL,
        data=payload,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "DrugDiscoveryWizard/1.0",
        },
    )

    print(f"[INFO] Searching COCONUT for '{search_term}' (limit={limit}, page={page})...")

    try:
        with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as resp:
            raw = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"COCONUT API error HTTP {e.code}: {body}") from e
    except urllib.error.URLError as e:
        raise RuntimeError(f"Cannot reach COCONUT: {e.reason}") from e

    # Parse response
    page_data = raw.get("data", {})
    compounds = page_data.get("data", [])

    results = []
    for c in compounds:
        results.append({
            "id": c.get("identifier", ""),
            "name": c.get("name", ""),
            "smiles": c.get("canonical_smiles", ""),
            "iupac_name": c.get("iupac_name", ""),
            "annotation_level": c.get("annotation_level", 0),
            "organism_count": c.get("organism_count", 0),
            "citation_count": c.get("citation_count", 0),
            "active": c.get("active", False),
        })

    output = {
        "source": "COCONUT",
        "query": search_term,
        "total_results": page_data.get("total", len(results)),
        "page": page_data.get("current_page", page),
        "last_page": page_data.get("last_page", 1),
        "compounds": results,
    }

    print(f"[OK] Found {len(results)} compounds (page {output['page']}/{output['last_page']}, "
          f"total: {output['total_results']})")

    return output


def search_coconut_multi_page(
    query: str,
    max_results: int = 100,
) -> dict[str, Any]:
    """Search COCONUT with automatic pagination.

    Args:
        query: Search term
        max_results: Maximum total results to fetch

    Returns:
        Combined results from all pages.
    """
    all_compounds: list[dict] = []
    page = 1
    per_page = min(max_results, 50)

    while len(all_compounds) < max_results:
        result = search_coconut(query=query, limit=per_page, page=page)
        all_compounds.extend(result["compounds"])

        if page >= result["last_page"]:
            break
        page += 1
        time.sleep(0.5)  # Rate limiting

    all_compounds = all_compounds[:max_results]

    return {
        "source": "COCONUT",
        "query": query,
        "total_fetched": len(all_compounds),
        "compounds": all_compounds,
    }


# ---------------------------------------------------------------------------
# ZINC20 Tranche Download Script Generator
# ---------------------------------------------------------------------------

def generate_zinc_download_script(
    output_path: str,
    subset: str = "drug-like",
    mw_range: tuple[int, int] | None = None,
    logp_range: tuple[float, float] | None = None,
    format: str = "smi",
    max_compounds: int = 10000,
) -> str:
    """Generate a download script for ZINC20 drug-like tranches.

    ZINC20 organizes compounds into tranches by MW and LogP. This generates
    a shell script to download the appropriate tranches.

    Args:
        output_path: Path to save the download script
        subset: Predefined subset (drug-like, lead-like, fragment-like)
        mw_range: (min_mw, max_mw) override
        logp_range: (min_logp, max_logp) override
        format: Download format (smi, sdf, mol2)
        max_compounds: Approximate max compounds to download

    Returns:
        The generated script content.
    """
    # Get property ranges
    if subset in ZINC_SUBSETS:
        props = ZINC_SUBSETS[subset]
        mw = mw_range or props["mw"]
        logp = logp_range or props["logp"]
    else:
        mw = mw_range or (250, 500)
        logp = logp_range or (-1, 5)

    # ZINC20 tranche codes:
    # MW tranches: A=0-200, B=200-250, C=250-300, D=300-325, E=325-350,
    #              F=350-375, G=375-400, H=400-425, I=425-450, J=450-500, K=500+
    # LogP tranches: A=-inf to -4, B=-4 to -2, C=-2 to 0, D=0 to 1,
    #                E=1 to 2, F=2 to 2.5, G=2.5 to 3, H=3 to 3.5,
    #                I=3.5 to 4, J=4 to 4.5, K=4.5 to 5, L=5+
    mw_codes = _mw_to_tranche_codes(mw[0], mw[1])
    logp_codes = _logp_to_tranche_codes(logp[0], logp[1])

    # Build download URLs
    urls = []
    for mw_code in mw_codes:
        for logp_code in logp_codes:
            tranche_id = f"{mw_code}{logp_code}"
            url = f"{ZINC20_BASE}/tranches/{tranche_id}/download.{format}"
            urls.append((tranche_id, url))

    script = f"""#!/bin/bash
# ==============================================================
# ZINC20 Compound Library Download Script
# Generated by Drug Discovery Wizard
# Subset: {subset}
# MW range: {mw[0]}-{mw[1]} Da
# LogP range: {logp[0]} to {logp[1]}
# Format: {format}
# Tranches: {len(urls)}
# ==============================================================

set -euo pipefail

OUTPUT_DIR="zinc20_{subset}"
mkdir -p "$OUTPUT_DIR"

echo "Downloading ZINC20 {subset} compounds..."
echo "  MW: {mw[0]}-{mw[1]}, LogP: {logp[0]} to {logp[1]}"
echo "  Tranches: {len(urls)}"
echo ""

TOTAL=0
"""

    for tranche_id, url in urls:
        script += f"""
# Tranche {tranche_id}
echo "Downloading tranche {tranche_id}..."
curl -sS -o "$OUTPUT_DIR/{tranche_id}.{format}" \\
    "{url}" 2>/dev/null || echo "  [SKIP] Tranche {tranche_id} not available"
if [ -f "$OUTPUT_DIR/{tranche_id}.{format}" ]; then
    COUNT=$(wc -l < "$OUTPUT_DIR/{tranche_id}.{format}" 2>/dev/null || echo 0)
    TOTAL=$((TOTAL + COUNT))
    echo "  Got $COUNT compounds from {tranche_id}"
fi
"""

    script += f"""
echo ""
echo "=== Download Complete ==="
echo "Total compounds: $TOTAL"
echo "Output directory: $OUTPUT_DIR/"
echo ""
echo "Next steps:"
echo "  1. Merge: cat $OUTPUT_DIR/*.{format} > zinc20_{subset}_all.{format}"
echo "  2. Filter: python3 scripts/drug_discovery_api.py admet-predict --smiles-file zinc20_{subset}_all.smi --output admet.json"
echo "  3. Dock: python3 scripts/drug_discovery_api.py dock --receptor receptor.pdbqt --smiles-file zinc20_{subset}_all.smi --output docking.json"
"""

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="\n") as f:
        f.write(script)

    print(f"[OK] ZINC20 download script saved to {output_path}")
    print(f"     Subset: {subset}, {len(urls)} tranches")
    print(f"     MW: {mw[0]}-{mw[1]}, LogP: {logp[0]} to {logp[1]}")
    print(f"     Run with: bash {output_path}")

    return script


def _mw_to_tranche_codes(mw_min: int, mw_max: int) -> list[str]:
    """Convert MW range to ZINC20 tranche letter codes."""
    boundaries = [
        (0, 200, "A"), (200, 250, "B"), (250, 300, "C"),
        (300, 325, "D"), (325, 350, "E"), (350, 375, "F"),
        (375, 400, "G"), (400, 425, "H"), (425, 450, "I"),
        (450, 500, "J"), (500, 9999, "K"),
    ]
    codes = []
    for low, high, code in boundaries:
        if mw_min < high and mw_max > low:
            codes.append(code)
    return codes or ["C", "D", "E", "F", "G", "H", "I", "J"]


def _logp_to_tranche_codes(logp_min: float, logp_max: float) -> list[str]:
    """Convert LogP range to ZINC20 tranche letter codes."""
    boundaries = [
        (-99, -4, "A"), (-4, -2, "B"), (-2, 0, "C"),
        (0, 1, "D"), (1, 2, "E"), (2, 2.5, "F"),
        (2.5, 3, "G"), (3, 3.5, "H"), (3.5, 4, "I"),
        (4, 4.5, "J"), (4.5, 5, "K"), (5, 99, "L"),
    ]
    codes = []
    for low, high, code in boundaries:
        if logp_min < high and logp_max > low:
            codes.append(code)
    return codes or ["C", "D", "E", "F", "G", "H", "I"]


# ---------------------------------------------------------------------------
# CLI (standalone)
# ---------------------------------------------------------------------------

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Compound library search (COCONUT / ZINC20)")
    sub = parser.add_subparsers(dest="command", required=True)

    # --- coconut ---
    p_coco = sub.add_parser("coconut", help="Search COCONUT natural products")
    q_group = p_coco.add_mutually_exclusive_group(required=True)
    q_group.add_argument("--query", help="Text search (name, organism, etc.)")
    q_group.add_argument("--smiles", help="SMILES for structure search")
    p_coco.add_argument("--limit", type=int, default=25)
    p_coco.add_argument("--page", type=int, default=1)
    p_coco.add_argument("--output", required=True)

    # --- zinc ---
    p_zinc = sub.add_parser("zinc", help="Generate ZINC20 download script")
    p_zinc.add_argument("--subset", default="drug-like",
                         choices=["drug-like", "lead-like", "fragment-like", "all-purchasable"])
    p_zinc.add_argument("--mw-range", type=int, nargs=2, metavar=("MIN", "MAX"))
    p_zinc.add_argument("--logp-range", type=float, nargs=2, metavar=("MIN", "MAX"))
    p_zinc.add_argument("--format", default="smi", choices=["smi", "sdf", "mol2"])
    p_zinc.add_argument("--output", required=True)

    args = parser.parse_args()

    if args.command == "coconut":
        result = search_coconut(
            query=args.query,
            smiles=args.smiles,
            limit=args.limit,
            page=args.page,
        )
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            json.dump(result, f, indent=2)
        print(f"[OK] Results saved to {args.output}")

    elif args.command == "zinc":
        generate_zinc_download_script(
            output_path=args.output,
            subset=args.subset,
            mw_range=tuple(args.mw_range) if args.mw_range else None,
            logp_range=tuple(args.logp_range) if args.logp_range else None,
            format=args.format,
        )


if __name__ == "__main__":
    main()
