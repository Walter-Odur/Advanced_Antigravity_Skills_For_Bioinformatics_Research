"""Compound Screening API.\n\nCompound library search (ChEMBL, COCONUT, ZINC20), molecular docking\n(AutoDock Vina), and ADMET filtering (two-tier: relaxed/strict).

Usage:
    python3 scripts/compound_screening_api.py <subcommand> --output <file> [options]
"""

import argparse
import json
import os
import re
import sys
import time
import urllib.parse
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any

# Locate http_client relative to this script.
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

import http_client

# ---------------------------------------------------------------------------
# API Clients
# ---------------------------------------------------------------------------

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
PDB_BASE = "https://data.rcsb.org/rest/v1/core"
ALPHAFOLD_BASE = "https://alphafold.ebi.ac.uk/api"
RCSB_FILES_BASE = "https://files.rcsb.org"
AF_FILES_BASE = "https://alphafold.ebi.ac.uk"

_CHEMBL_CLIENT = http_client.HttpClient(CHEMBL_BASE, qps=5.0)
_PDB_CLIENT = http_client.HttpClient(PDB_BASE, qps=5.0)
_ALPHAFOLD_CLIENT = http_client.HttpClient(ALPHAFOLD_BASE, qps=5.0)


def _build_url(base: str, path: str, params: dict | None = None) -> str:
    """Build a full URL from base + path + optional query parameters."""
    url = f"{base}/{path.lstrip('/')}"
    if params:
        qs = urllib.parse.urlencode(params)
        url = f"{url}?{qs}"
    return url


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_output(data: Any, output_file: str) -> None:
    """Write data to a JSON file."""
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=str)
    print(f"Success! Data written to: {output_file}")


def write_text_output(text: str, output_file: str) -> None:
    """Write plain text to a file."""
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    with open(output_file, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    print(f"Success! Script written to: {output_file}")


# ---------------------------------------------------------------------------


def cmd_query_chembl(args):
    """Query ChEMBL for bioactive compounds against a target."""
    target_id = args.target_id
    pchembl_min = args.pchembl_min
    assay_type = args.assay_type
    limit = args.limit

    print(f"Querying ChEMBL for {target_id} "
          f"(pChEMBL >= {pchembl_min}, assay_type={assay_type}, "
          f"limit={limit})...")

    params = {
        "target_chembl_id": target_id,
        "pchembl_value__gte": pchembl_min,
        "assay_type": assay_type,
        "limit": limit,
        "format": "json",
    }

    try:
        url = _build_url(CHEMBL_BASE, "activity.json", params)
        data = _CHEMBL_CLIENT.fetch_json(url)
    except Exception as e:
        print(f"ERROR: ChEMBL query failed: {e}", file=sys.stderr)
        sys.exit(1)

    # Extract and clean results
    activities = data.get("activities", data.get("results", []))
    if isinstance(data, dict) and "activities" not in data:
        # API may return differently depending on version
        activities = data.get("activities", [])

    results = []
    for act in activities:
        smiles = act.get("canonical_smiles")
        pchembl = act.get("pchembl_value")
        if smiles and pchembl:
            results.append({
                "molecule_chembl_id": act.get("molecule_chembl_id"),
                "canonical_smiles": smiles,
                "pchembl_value": float(pchembl),
                "standard_type": act.get("standard_type"),
                "standard_value": act.get("standard_value"),
                "standard_units": act.get("standard_units"),
                "assay_chembl_id": act.get("assay_chembl_id"),
            })

    # Sort by pChEMBL descending
    results.sort(key=lambda x: x["pchembl_value"], reverse=True)

    output = {
        "query": {
            "target_chembl_id": target_id,
            "pchembl_min": pchembl_min,
            "assay_type": assay_type,
        },
        "total_hits": len(results),
        "compounds": results,
    }

    print(f"  Retrieved {len(results)} compounds")
    if results:
        best = results[0]
        print(f"  Top hit: {best['molecule_chembl_id']} "
              f"(pChEMBL={best['pchembl_value']})")

    write_output(output, args.output)

# ---------------------------------------------------------------------------


def cmd_search_chembl_target(args):
    """Search ChEMBL for targets by name."""
    query = args.query
    limit = args.limit

    print(f"Searching ChEMBL targets for: {query} (limit={limit})...")

    try:
        url = _build_url(CHEMBL_BASE, "target/search.json",
                         {"q": query, "limit": limit, "format": "json"})
        data = _CHEMBL_CLIENT.fetch_json(url)
    except Exception as e:
        print(f"ERROR: ChEMBL target search failed: {e}", file=sys.stderr)
        sys.exit(1)

    targets = data.get("targets", [])
    results = []
    for t in targets:
        results.append({
            "target_chembl_id": t.get("target_chembl_id"),
            "pref_name": t.get("pref_name"),
            "target_type": t.get("target_type"),
            "organism": t.get("organism"),
        })

    output = {
        "query": query,
        "total_hits": len(results),
        "targets": results,
    }

    print(f"  Found {len(results)} targets")
    write_output(output, args.output)

# ---------------------------------------------------------------------------


def cmd_admet_filter(args):
    """Run ADMET filtering cascade on compounds."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, FilterCatalog, rdMolDescriptors
    except ImportError:
        print("ERROR: RDKit not installed. Run: uv pip install rdkit",
              file=sys.stderr)
        sys.exit(1)

    # Load compounds
    if args.csv:
        try:
            import csv
            compounds = []
            with open(args.csv, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    smiles = row.get(args.smiles_col, "")
                    name = row.get(args.name_col, f"compound_{len(compounds)+1}")
                    if smiles:
                        compounds.append((smiles, name))
        except Exception as e:
            print(f"ERROR: Failed to read CSV: {e}", file=sys.stderr)
            sys.exit(1)
    elif args.smiles:
        names = args.names or [f"compound_{i+1}"
                                for i in range(len(args.smiles))]
        compounds = list(zip(args.smiles, names))
    else:
        print("ERROR: Provide --smiles or --csv", file=sys.stderr)
        sys.exit(1)

    # PAINS filter setup
    pains_params = FilterCatalog.FilterCatalogParams()
    pains_params.AddCatalog(
        FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A
    )
    pains_params.AddCatalog(
        FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B
    )
    pains_params.AddCatalog(
        FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C
    )
    pains_catalog = FilterCatalog.FilterCatalog(pains_params)

    # SA score (Ertl & Schuffenhauer)
    try:
        from rdkit.Chem import RDConfig
        sa_module_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if sa_module_path not in sys.path:
            sys.path.insert(0, sa_module_path)
        import sascorer
        sa_available = True
    except (ImportError, OSError):
        sa_available = False

    sa_threshold = args.sa_threshold
    strictness = getattr(args, "strictness", "strict")

    # Tier-specific thresholds for the RDKit cascade
    if strictness == "relaxed":
        mw_max, logp_max, hbd_max, hba_max = 700, 7.0, 7, 15
        max_violations = 3
        max_pains = 1
        if sa_threshold == 6.0:  # only override if user didn't set custom
            sa_threshold = 8.0
        tier_label = "Tier 1 (relaxed) -- post-docking, for REINVENT seeding"
    else:
        mw_max, logp_max, hbd_max, hba_max = 500, 5.0, 5, 10
        max_violations = 1
        max_pains = 0
        tier_label = "Tier 2 (strict) -- post-REINVENT, for downstream analysis"

    print(f"ADMET filtering {len(compounds)} compounds "
          f"(SA threshold={sa_threshold}, strictness={strictness})...")

    results = []
    pass_count = 0

    for smiles, name in compounds:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            results.append({
                "name": name, "smiles": smiles, "error": "Invalid SMILES",
                "overall_pass": False,
            })
            continue

        # Lipinski
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        violations = sum([mw > mw_max, logp > logp_max,
                          hbd > hbd_max, hba > hba_max])
        lipinski_pass = violations <= max_violations

        # PAINS
        pains_matches = pains_catalog.GetMatches(mol)
        pains_alerts = [m.GetDescription() for m in pains_matches]
        pains_pass = len(pains_alerts) <= max_pains

        # SA Score
        sa_score = None
        sa_pass = True
        if sa_available:
            sa_score = round(sascorer.calculateScore(mol), 2)
            sa_pass = sa_score <= sa_threshold

        overall = lipinski_pass and pains_pass and sa_pass
        if overall:
            pass_count += 1

        results.append({
            "name": name,
            "smiles": Chem.MolToSmiles(mol),
            "mw": round(mw, 1),
            "logp": round(logp, 2),
            "hbd": hbd,
            "hba": hba,
            "lipinski_violations": violations,
            "lipinski_pass": lipinski_pass,
            "pains_pass": pains_pass,
            "pains_alerts": pains_alerts,
            "sa_score": sa_score,
            "sa_pass": sa_pass,
            "overall_pass": overall,
        })

    output = {
        "strictness": strictness,
        "tier": tier_label,
        "total_compounds": len(compounds),
        "passed": pass_count,
        "failed": len(compounds) - pass_count,
        "sa_threshold": sa_threshold,
        "compounds": results,
    }

    print(f"  {pass_count}/{len(compounds)} compounds passed all filters")
    print(f"  Tier: {tier_label}")
    write_output(output, args.output)

# ---------------------------------------------------------------------------


def cmd_dock(args):
    """Dock compounds using AutoDock Vina."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        print("ERROR: RDKit not installed. Run: uv pip install rdkit",
              file=sys.stderr)
        sys.exit(1)

    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
    except ImportError:
        print("ERROR: Meeko not installed. Run: uv pip install meeko",
              file=sys.stderr)
        sys.exit(1)

    # Parse binding site
    if args.site:
        with open(args.site) as f:
            site = json.load(f)
        center = (site["center_x"], site["center_y"], site["center_z"])
        box_size = (site["size_x"], site["size_y"], site["size_z"])
    elif args.center:
        center = tuple(args.center)
        box_size = tuple(args.box_size)
    else:
        print("ERROR: Provide --site or --center", file=sys.stderr)
        sys.exit(1)

    receptor = args.receptor
    if not os.path.exists(receptor):
        print(f"ERROR: Receptor file not found: {receptor}", file=sys.stderr)
        sys.exit(1)

    # Collect compounds
    if args.smiles:
        names = args.names or [f"compound_{i+1}"
                                for i in range(len(args.smiles))]
        compounds = list(zip(args.smiles, names))
    elif args.csv:
        import csv
        compounds = []
        with open(args.csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                smi = row.get(args.smiles_col, "")
                name = row.get(args.name_col, f"compound_{len(compounds)+1}")
                if smi:
                    compounds.append((smi, name))
    else:
        print("ERROR: Provide --smiles or --csv", file=sys.stderr)
        sys.exit(1)

    output_dir = os.path.dirname(args.output) or "."
    results = []

    print(f"Docking {len(compounds)} compounds against {receptor}")
    print(f"  Center: ({center[0]}, {center[1]}, {center[2]})")
    print(f"  Box: {box_size[0]} x {box_size[1]} x {box_size[2]} A")

    import shutil
    import subprocess

    vina_cmd = shutil.which("vina") or shutil.which("vina.exe")
    use_cli = True

    # Try Python API first
    try:
        from vina import Vina
        use_cli = False
        print("  Using Vina Python API")
    except ImportError:
        if vina_cmd:
            print(f"  Using Vina CLI: {vina_cmd}")
        else:
            print("ERROR: Vina not found. Install via conda or download "
                  "vina.exe from github.com/ccsb-scripps/AutoDock-Vina/releases",
                  file=sys.stderr)
            sys.exit(1)

    for i, (smiles, name) in enumerate(compounds, 1):
        print(f"  [{i}/{len(compounds)}] {name}...", end=" ", flush=True)
        start = time.time()

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            results.append({"name": name, "smiles": smiles,
                            "error": "Invalid SMILES", "score": None})
            print("FAILED (invalid SMILES)")
            continue

        # Prepare ligand PDBQT
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) == -1:
            results.append({"name": name, "smiles": smiles,
                            "error": "3D embedding failed", "score": None})
            print("FAILED (3D)")
            continue

        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)

        try:
            preparator = MoleculePreparation()
            mol_setups = preparator.prepare(mol)
            pdbqt_str, is_ok, err = PDBQTWriterLegacy.write_string(
                mol_setups[0]
            )
            if not is_ok:
                results.append({"name": name, "smiles": smiles,
                                "error": f"PDBQT error: {err}", "score": None})
                print(f"FAILED (PDBQT)")
                continue
        except Exception as e:
            results.append({"name": name, "smiles": smiles,
                            "error": str(e), "score": None})
            print(f"FAILED (Meeko)")
            continue

        # Dock
        out_pdbqt = os.path.join(output_dir, f"{name}_docked.pdbqt")
        score = None
        scores = []

        if not use_cli:
            try:
                v = Vina(sf_name="vina")
                v.set_receptor(receptor)
                v.set_ligand_from_string(pdbqt_str)
                v.compute_vina_maps(center=list(center),
                                     box_size=list(box_size))
                v.dock(exhaustiveness=args.exhaustiveness,
                       n_poses=args.n_poses)
                energies = v.energies()
                scores = [e[0] for e in energies] if energies else []
                score = scores[0] if scores else None
                v.write_poses(out_pdbqt, n_poses=args.n_poses, overwrite=True)
            except Exception as e:
                results.append({"name": name, "smiles": smiles,
                                "error": str(e), "score": None})
                print(f"FAILED (Vina API)")
                continue
        else:
            lig_file = os.path.join(output_dir, f"{name}_ligand.pdbqt")
            with open(lig_file, "w") as f:
                f.write(pdbqt_str)
            cmd = [
                vina_cmd,
                "--receptor", receptor, "--ligand", lig_file,
                "--center_x", str(center[0]),
                "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x", str(box_size[0]),
                "--size_y", str(box_size[1]),
                "--size_z", str(box_size[2]),
                "--exhaustiveness", str(args.exhaustiveness),
                "--num_modes", str(args.n_poses),
                "--out", out_pdbqt,
            ]
            try:
                proc = subprocess.run(cmd, capture_output=True, text=True,
                                       timeout=300)
                if proc.returncode != 0:
                    results.append({"name": name, "smiles": smiles,
                                    "error": proc.stderr[:200], "score": None})
                    print(f"FAILED (Vina CLI)")
                    continue
                for line in proc.stdout.splitlines():
                    m = re.match(
                        r"\s+(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", line
                    )
                    if m:
                        scores.append(float(m.group(2)))
                score = scores[0] if scores else None
            except Exception as e:
                results.append({"name": name, "smiles": smiles,
                                "error": str(e), "score": None})
                print(f"FAILED (CLI)")
                continue

        elapsed = time.time() - start
        results.append({
            "name": name,
            "smiles": smiles,
            "score_kcal_mol": score,
            "num_modes": len(scores),
            "all_scores": scores,
            "output_pdbqt": out_pdbqt,
        })
        print(f"{score:.1f} kcal/mol ({elapsed:.1f}s)")

    # Sort by score
    scored = [r for r in results if r.get("score_kcal_mol") is not None]
    failed = [r for r in results if r.get("score_kcal_mol") is None]
    scored.sort(key=lambda r: r["score_kcal_mol"])

    output = {
        "receptor": receptor,
        "binding_site": {"center": list(center), "box_size": list(box_size)},
        "exhaustiveness": args.exhaustiveness,
        "total_compounds": len(compounds),
        "successful": len(scored),
        "failed": len(failed),
        "ranked_results": scored + failed,
    }

    print(f"\n  {len(scored)}/{len(compounds)} docked successfully")
    if scored:
        print(f"  Best: {scored[0]['name']} ({scored[0]['score_kcal_mol']:.1f} "
              f"kcal/mol)")

    write_output(output, args.output)

# ---------------------------------------------------------------------------


def cmd_admet_predict(args):
    """ADMETlab 3.0 prediction — delegates to admetlab_client.py."""
    from admetlab_client import predict_admet, apply_thresholds, _rdkit_fallback
    try:
        predictions = predict_admet(args.smiles)
        output = apply_thresholds(predictions, args.smiles, args.names)
    except RuntimeError as e:
        print(f"[WARNING] ADMETlab API unavailable: {e}", file=sys.stderr)
        print("[INFO] Falling back to RDKit offline ADMET", file=sys.stderr)
        output = _rdkit_fallback(args.smiles, args.names)
    write_output(output, args.output)
    print(f"[OK] ADMET predictions saved to {args.output}")


def cmd_query_coconut(args):
    """Search COCONUT natural products — delegates to compound_libraries.py."""
    from compound_libraries import search_coconut
    result = search_coconut(
        query=args.query,
        smiles=args.smiles,
        limit=args.limit,
        page=args.page,
    )
    write_output(result, args.output)
    print(f"[OK] {len(result['compounds'])} compounds saved to {args.output}")


def cmd_query_zinc(args):
    """Generate ZINC20 download script — delegates to compound_libraries.py."""
    from compound_libraries import generate_zinc_download_script
    generate_zinc_download_script(
        output_path=args.output,
        subset=args.subset,
        mw_range=tuple(args.mw_range) if args.mw_range else None,
        logp_range=tuple(args.logp_range) if args.logp_range else None,
        format=args.format,
    )


# ---------------------------------------------------------------------------
# CLI definition
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compound Screening — library search, docking, ADMET filtering"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # --- query-chembl ---
    p_chembl = sub.add_parser("query-chembl", help="Query ChEMBL for bioactive compounds")
    p_chembl.add_argument("--target-id", required=True, help="ChEMBL target ID")
    p_chembl.add_argument("--pchembl-min", type=float, default=5.0)
    p_chembl.add_argument("--limit", type=int, default=25)
    p_chembl.add_argument("--assay-type", default="B",
                           help="Assay type (B=binding, F=functional)")
    p_chembl.add_argument("--output", required=True)
    p_chembl.set_defaults(func=cmd_query_chembl)

    # --- search-chembl-target ---
    p_search = sub.add_parser("search-chembl-target", help="Search ChEMBL targets by name")
    p_search.add_argument("--query", required=True, help="Search term")
    p_search.add_argument("--limit", type=int, default=10)
    p_search.add_argument("--output", required=True)
    p_search.set_defaults(func=cmd_search_chembl_target)

    # --- admet-filter ---
    p_admet = sub.add_parser("admet-filter", help="Run ADMET filtering cascade on compounds")
    inp = p_admet.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles", nargs="+", help="SMILES strings")
    inp.add_argument("--csv", help="CSV file with compounds")
    p_admet.add_argument("--names", nargs="+", help="Compound names (with --smiles)")
    p_admet.add_argument("--smiles-col", default="SMILES", help="SMILES column in CSV")
    p_admet.add_argument("--name-col", default="Name", help="Name column in CSV")
    p_admet.add_argument("--sa-threshold", type=float, default=6.0,
                          help="SA score threshold (default: 6.0)")
    p_admet.add_argument("--strictness", default="strict",
                          choices=["strict", "relaxed"],
                          help="'relaxed' = post-docking seeds for REINVENT, "
                               "'strict' = post-REINVENT for downstream (default)")
    p_admet.add_argument("--output", required=True, help="Output JSON file path")
    p_admet.set_defaults(func=cmd_admet_filter)

    # --- dock ---
    p_dock = sub.add_parser("dock", help="Dock compounds using AutoDock Vina")
    p_dock.add_argument("--receptor", required=True, help="Receptor PDBQT file")
    dinp = p_dock.add_mutually_exclusive_group(required=True)
    dinp.add_argument("--smiles", nargs="+", help="SMILES strings")
    dinp.add_argument("--csv", help="CSV file with compounds")
    p_dock.add_argument("--names", nargs="+")
    p_dock.add_argument("--smiles-col", default="SMILES")
    p_dock.add_argument("--name-col", default="Name")
    p_dock.add_argument("--site", help="Binding site JSON file")
    p_dock.add_argument("--center", nargs=3, type=float, metavar=("X", "Y", "Z"))
    p_dock.add_argument("--size", nargs=3, type=float,
                         default=[20, 20, 20], metavar=("X", "Y", "Z"))
    p_dock.add_argument("--exhaustiveness", type=int, default=32)
    p_dock.add_argument("--n-poses", type=int, default=9)
    p_dock.add_argument("--output", required=True)
    p_dock.set_defaults(func=cmd_dock)

    # --- admet-predict ---
    p_admet_api = sub.add_parser("admet-predict",
        help="Full ADMET prediction via ADMETlab 3.0 API (119 endpoints)")
    p_admet_api.add_argument("--smiles", nargs="+", required=True,
                              help="SMILES strings to predict")
    p_admet_api.add_argument("--names", nargs="+", help="Compound names")
    p_admet_api.add_argument("--output", required=True,
                              help="Output JSON file path")
    p_admet_api.set_defaults(func=cmd_admet_predict)

    # --- query-coconut ---
    p_coco = sub.add_parser("query-coconut",
        help="Search COCONUT 2.0 for natural products")
    coco_inp = p_coco.add_mutually_exclusive_group(required=True)
    coco_inp.add_argument("--query", help="Search term")
    coco_inp.add_argument("--smiles", help="SMILES for similarity search")
    p_coco.add_argument("--limit", type=int, default=25)
    p_coco.add_argument("--page", type=int, default=1)
    p_coco.add_argument("--output", required=True)
    p_coco.set_defaults(func=cmd_query_coconut)

    # --- query-zinc ---
    p_zinc = sub.add_parser("query-zinc",
        help="Generate ZINC20 tranche download script")
    p_zinc.add_argument("--subset", default="drug-like",
                         choices=["drug-like", "lead-like", "fragment-like",
                                  "all-purchasable"])
    p_zinc.add_argument("--mw-range", nargs=2, type=float,
                         metavar=("MIN", "MAX"))
    p_zinc.add_argument("--logp-range", nargs=2, type=float,
                         metavar=("MIN", "MAX"))
    p_zinc.add_argument("--format", default="smi",
                         choices=["smi", "sdf", "mol2"],
                         help="Download format (default: smi)")
    p_zinc.add_argument("--output", required=True)
    p_zinc.set_defaults(func=cmd_query_zinc)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
