"""Drug Discovery Wizard CLI.

Unified CLI tool for structure-based drug discovery. Provides subcommands for
receptor preparation, docking, ADMET filtering, database queries, and HPC
script generation. All subcommands write output to a file specified by
--output.

Usage:
    python3 scripts/drug_discovery_api.py <subcommand> --output <file> [options]
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
# Subcommand: prepare-receptor
# ---------------------------------------------------------------------------

def cmd_prepare_receptor(args):
    """Fetch and prepare a receptor for docking."""
    try:
        from Bio.PDB import PDBParser, PDBIO, Select
    except ImportError:
        print("ERROR: BioPython not installed. Run: uv pip install biopython",
              file=sys.stderr)
        sys.exit(1)

    result = {"status": "success", "receptor_id": None, "source": None,
              "chains": [], "output_files": []}

    pdb_id = args.pdb_id
    uniprot_id = args.uniprot_id
    input_file = args.input_file

    if pdb_id:
        # Fetch from RCSB PDB
        print(f"Fetching PDB structure: {pdb_id}")
        pdb_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        try:
            pdb_text = http_client.HttpClient(
                RCSB_FILES_BASE, qps=3.0
            ).fetch_text(f"/download/{pdb_id.upper()}.pdb")
            input_file = os.path.join(
                os.path.dirname(args.output), f"{pdb_id.upper()}.pdb"
            )
            os.makedirs(os.path.dirname(input_file) or ".", exist_ok=True)
            with open(input_file, "w") as f:
                f.write(pdb_text)
            result["source"] = "RCSB PDB"
            result["receptor_id"] = pdb_id.upper()
            print(f"  Downloaded: {input_file}")
        except Exception as e:
            print(f"ERROR: Failed to fetch PDB {pdb_id}: {e}", file=sys.stderr)
            sys.exit(1)

    elif uniprot_id:
        # Fetch from AlphaFold DB
        print(f"Fetching AlphaFold model: {uniprot_id}")
        af_url = (f"https://alphafold.ebi.ac.uk/files/"
                  f"AF-{uniprot_id.upper()}-F1-model_v4.pdb")
        try:
            pdb_text = http_client.HttpClient(
                AF_FILES_BASE, qps=3.0
            ).fetch_text(f"/files/AF-{uniprot_id.upper()}-F1-model_v4.pdb")
            input_file = os.path.join(
                os.path.dirname(args.output),
                f"AF-{uniprot_id.upper()}-F1.pdb"
            )
            os.makedirs(os.path.dirname(input_file) or ".", exist_ok=True)
            with open(input_file, "w") as f:
                f.write(pdb_text)
            result["source"] = "AlphaFold DB"
            result["receptor_id"] = uniprot_id.upper()
            result["warning"] = (
                "AlphaFold models lack experimental ligand data. "
                "Binding site coordinates must be provided manually or "
                "inferred from homologous PDB structures."
            )
            print(f"  Downloaded: {input_file}")
            print(f"  WARNING: {result['warning']}")
        except Exception as e:
            print(f"ERROR: Failed to fetch AlphaFold model: {e}",
                  file=sys.stderr)
            sys.exit(1)

    elif input_file:
        result["source"] = "local file"
        result["receptor_id"] = Path(input_file).stem
    else:
        print("ERROR: Provide --pdb-id, --uniprot-id, or --input-file",
              file=sys.stderr)
        sys.exit(1)

    # Parse structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", input_file)
    model = structure[0]

    # List chains
    chains = [chain.id for chain in model.get_chains()]
    result["chains"] = chains
    print(f"  Chains found: {', '.join(chains)}")

    # Select chain if specified
    target_chain = args.chain

    class ChainAndCleanSelect(Select):
        def accept_chain(self, chain):
            if target_chain:
                return chain.id == target_chain
            return True
        def accept_residue(self, residue):
            hetflag = residue.id[0]
            if hetflag == "W":
                return False  # Remove waters
            return True

    # Save cleaned PDB
    clean_pdb = os.path.join(
        os.path.dirname(args.output),
        f"{result['receptor_id']}_clean.pdb"
    )
    io = PDBIO()
    io.set_structure(structure)
    io.save(clean_pdb, ChainAndCleanSelect())
    result["output_files"].append(clean_pdb)
    print(f"  Cleaned PDB saved: {clean_pdb}")

    # Try PDBQT conversion with Open Babel
    pdbqt_file = clean_pdb.replace(".pdb", ".pdbqt")
    try:
        from openbabel import openbabel
        obconv = openbabel.OBConversion()
        obconv.SetInAndOutFormats("pdb", "pdbqt")
        mol = openbabel.OBMol()
        obconv.ReadFile(mol, clean_pdb)
        mol.AddHydrogens(True)  # polar only
        obconv.WriteFile(mol, pdbqt_file)
        result["output_files"].append(pdbqt_file)
        print(f"  PDBQT saved: {pdbqt_file}")
    except ImportError:
        try:
            from meeko import PDBQTWriterLegacy
            result["note"] = (
                "Open Babel not available. PDBQT conversion requires "
                "obabel or ADFR Suite's prepare_receptor script."
            )
            print(f"  NOTE: {result['note']}")
        except ImportError:
            result["note"] = (
                "Neither Open Babel nor Meeko available for PDBQT conversion. "
                "Install: conda install -c conda-forge openbabel"
            )
            print(f"  NOTE: {result['note']}")

    write_output(result, args.output)


# ---------------------------------------------------------------------------
# Subcommand: query-chembl
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
# Subcommand: search-chembl-target
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
# Subcommand: admet-filter
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
# Subcommand: dock
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
# Subcommand: generate-reinvent-config
# ---------------------------------------------------------------------------

def cmd_generate_reinvent_config(args):
    """Generate REINVENT 4 TOML configuration."""
    from datetime import datetime

    seeds = args.seeds
    if args.seeds_csv:
        import csv
        with open(args.seeds_csv) as f:
            reader = csv.DictReader(f)
            seeds = [row[args.smiles_col] for row in reader
                     if row.get(args.smiles_col)]

    if not seeds:
        print("ERROR: No seed compounds provided", file=sys.stderr)
        sys.exit(1)

    print(f"Generating REINVENT 4 config with {len(seeds)} seed compounds")

    seeds_toml = ", ".join(f'"{s}"' for s in seeds)

    # Load binding site if provided
    docking_block = ""
    if args.receptor and args.site:
        with open(args.site) as f:
            site = json.load(f)
        docking_block = f"""
[[scoring.component]]
[scoring.component.DockStream]
[[scoring.component.DockStream.endpoint]]
name = "Vina docking"
weight = 0.4
[scoring.component.DockStream.endpoint.params]
receptor_path = "{args.receptor}"
center_x = {site['center_x']}
center_y = {site['center_y']}
center_z = {site['center_z']}
size_x = {site['size_x']}
size_y = {site['size_y']}
size_z = {site['size_z']}
exhaustiveness = 8
[scoring.component.DockStream.endpoint.transform]
type = "reverse_sigmoid"
high = -8.0
low = -12.0
k = 0.3
"""

    config = f"""# REINVENT 4 Configuration - Reinforcement Learning
# Generated by Drug Discovery Wizard on {datetime.now().isoformat()}

[run_type]
type = "reinforcement_learning"

[parameters]
summary_csv_prefix = "{args.output_dir}/rl_summary"
use_checkpoint = false
num_steps = {args.num_steps}
batch_size = {args.batch_size}

[parameters.reinforcement_learning]
prior = "{args.prior}"
agent = "{args.prior}"
learning_rate = 0.0001
reset = 0
reset_score_cutoff = 0.5
margin_threshold = 50

[parameters.reinforcement_learning.diversity_filter]
type = "IdenticalMurckoScaffold"
bucket_size = 25
minscore = 0.4
minsimilarity = 0.4

[[stage]]
max_steps = {args.num_steps}
chkpt_file = "{args.output_dir}/stage1.chkpt"

[stage.inception]
smiles = [{seeds_toml}]
deduplicate = true
memory_size = 100

[scoring]
type = "geometric_mean"

[[scoring.component]]
[scoring.component.custom_alerts]
[[scoring.component.custom_alerts.endpoint]]
name = "Unwanted SMARTS"
weight = 1.0
params.smarts = [
    "[#6](=[#16])(-[#7])-[#16]",
]
[scoring.component.custom_alerts.endpoint.transform]
type = "step"
low = 0.0
high = 0.5

[[scoring.component]]
[scoring.component.QED]
[[scoring.component.QED.endpoint]]
name = "QED score"
weight = 0.3
[scoring.component.QED.endpoint.transform]
type = "sigmoid"
high = {args.target_qed}
k = 0.5

[[scoring.component]]
[scoring.component.MolecularWeight]
[[scoring.component.MolecularWeight.endpoint]]
name = "MW filter"
weight = 0.1
[scoring.component.MolecularWeight.endpoint.transform]
type = "double_sigmoid"
low = 200.0
high = 500.0
coef_div = 500.0
coef_si = 20.0
coef_se = 20.0

[[scoring.component]]
[scoring.component.custom_property]
[[scoring.component.custom_property.endpoint]]
name = "LogP"
weight = 0.2
[scoring.component.custom_property.endpoint.params]
property = "logp"
[scoring.component.custom_property.endpoint.transform]
type = "double_sigmoid"
low = 0.5
high = {args.target_logp}
coef_div = 500.0
coef_si = 10.0
coef_se = 10.0
{docking_block}
[logging]
type = "local"
logging_path = "{args.output_dir}/logs"
job_name = "drug_discovery_wizard_rl"
"""

    write_text_output(config, args.output)


# ---------------------------------------------------------------------------
# Subcommand: generate-hpc-script
# ---------------------------------------------------------------------------

def cmd_generate_hpc_script(args):
    """Generate SLURM submission script for HPC jobs."""
    from datetime import datetime

    job_type = args.job_type

    # Build SLURM header
    header_lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_type}",
        f"#SBATCH --partition={'gpu' if args.ngpu > 0 else 'cpu'}",
        f"#SBATCH --time={args.time}",
        f"#SBATCH --mem={args.mem}",
        f"#SBATCH --cpus-per-task={args.ncpus}",
        f"#SBATCH --output={job_type}_%j.out",
        f"#SBATCH --error={job_type}_%j.err",
    ]
    if args.ngpu > 0:
        header_lines.append(
            f"#SBATCH --gres=gpu:{args.gpu}:{args.ngpu}"
        )
    if args.array:
        header_lines.append(f"#SBATCH --array={args.array}")

    header = "\n".join(header_lines)
    header += f"\n\n# Generated by Drug Discovery Wizard on {datetime.now().isoformat()}\n"

    if job_type == "reinvent":
        if not args.config:
            print("ERROR: --config required for reinvent jobs", file=sys.stderr)
            sys.exit(1)
        body = f"""
echo "Job started at: $(date)"
nvidia-smi

module purge
module load anaconda3/2024.02
conda activate {args.conda_env}

which reinvent || {{ echo "ERROR: reinvent not found"; exit 1; }}

reinvent {args.config}

echo "Job completed at: $(date)"
"""
    elif job_type == "gromacs-md":
        nsteps = args.production_ns * 500000
        body = f"""
echo "Job started at: $(date)"
nvidia-smi

module purge
module load gromacs/2024.3-gpu

WORKDIR=$SLURM_SUBMIT_DIR
cd $WORKDIR

# Energy Minimisation
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em -nb gpu -v

# NVT Equilibration (100 ps)
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -nb gpu -pme gpu -v

# NPT Equilibration (100 ps)
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -nb gpu -pme gpu -v

# Production MD ({args.production_ns} ns = {nsteps} steps @ 2 fs)
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o production.tpr
gmx mdrun -deffnm production -nb gpu -pme gpu -bonded gpu -update gpu -ntomp $SLURM_CPUS_PER_TASK -v

# Post-processing
echo "1 1" | gmx rms -s production.tpr -f production.xtc -o rmsd.xvg
echo "1" | gmx rmsf -s production.tpr -f production.xtc -o rmsf.xvg
echo "1" | gmx gyrate -s production.tpr -f production.xtc -o gyrate.xvg

echo "Job completed at: $(date)"
"""
    elif job_type == "vina-screen":
        if not args.receptor or not args.ligand_dir:
            print("ERROR: --receptor and --ligand-dir required",
                  file=sys.stderr)
            sys.exit(1)
        body = f"""
echo "Array task $SLURM_ARRAY_TASK_ID"

module purge
module load autodock-vina/1.2.5

RECEPTOR="{args.receptor}"
LIGAND_DIR="{args.ligand_dir}"
SITE_CONFIG="{args.site_config}"
OUTPUT_DIR="results"
mkdir -p $OUTPUT_DIR

LIGAND=$(ls $LIGAND_DIR/*.pdbqt | sed -n "${{SLURM_ARRAY_TASK_ID}}p")
BASENAME=$(basename $LIGAND .pdbqt)

if [ -z "$LIGAND" ]; then exit 0; fi

CENTER_X=$(python3 -c "import json; print(json.load(open('$SITE_CONFIG'))['center_x'])")
CENTER_Y=$(python3 -c "import json; print(json.load(open('$SITE_CONFIG'))['center_y'])")
CENTER_Z=$(python3 -c "import json; print(json.load(open('$SITE_CONFIG'))['center_z'])")
SIZE_X=$(python3 -c "import json; print(json.load(open('$SITE_CONFIG'))['size_x'])")
SIZE_Y=$(python3 -c "import json; print(json.load(open('$SITE_CONFIG'))['size_y'])")
SIZE_Z=$(python3 -c "import json; print(json.load(open('$SITE_CONFIG'))['size_z'])")

vina --receptor $RECEPTOR --ligand $LIGAND \\
     --center_x $CENTER_X --center_y $CENTER_Y --center_z $CENTER_Z \\
     --size_x $SIZE_X --size_y $SIZE_Y --size_z $SIZE_Z \\
     --exhaustiveness 32 --num_modes 9 \\
     --out $OUTPUT_DIR/${{BASENAME}}_out.pdbqt \\
     --log $OUTPUT_DIR/${{BASENAME}}_log.txt

echo "Done: $BASENAME"
"""
    else:
        print(f"ERROR: Unknown job type: {job_type}", file=sys.stderr)
        sys.exit(1)

    script = header + "\n" + body
    write_text_output(script, args.output)


# ---------------------------------------------------------------------------
# Subcommand: generate-md-setup
# ---------------------------------------------------------------------------

def cmd_generate_md_setup(args):
    """Generate GROMACS MD system setup script."""
    from datetime import datetime

    nsteps = args.production_ns * 500000

    script = f"""#!/bin/bash
# GROMACS MD System Setup Script
# Generated by Drug Discovery Wizard on {datetime.now().isoformat()}
set -euo pipefail

PDB_FILE="{args.pdb}"
FF="{args.ff}"
WATER="{args.water}"
BOX_DIST={args.box_dist}

echo "=== GROMACS MD Setup ==="
echo "Input: $PDB_FILE | FF: $FF | Water: $WATER"

# Generate topology
gmx pdb2gmx -f $PDB_FILE -o protein.gro -p topol.top -ignh -ff $FF -water $WATER

# Simulation box
gmx editconf -f protein.gro -o box.gro -c -d $BOX_DIST -bt dodecahedron

# Solvate
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top

# Add ions
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

echo "=== Setup Complete ==="
echo "Files: protein.gro, box.gro, solv.gro, solv_ions.gro, topol.top"

# Generate MDP files
cat > em.mdp << 'EOF'
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
EOF

cat > nvt.mdp << 'EOF'
define      = -DPOSRES
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout-compressed = 5000
nstenergy   = 5000
nstlog      = 5000
continuation = no
constraint_algorithm = lincs
constraints = h-bonds
cutoff-scheme = Verlet
nstlist     = 10
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1   0.1
ref_t       = 300   300
pcoupl      = no
pbc         = xyz
gen_vel     = yes
gen_temp    = 300
gen_seed    = -1
EOF

cat > npt.mdp << 'EOF'
define      = -DPOSRES
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout-compressed = 5000
nstenergy   = 5000
nstlog      = 5000
continuation = yes
constraint_algorithm = lincs
constraints = h-bonds
cutoff-scheme = Verlet
nstlist     = 10
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1   0.1
ref_t       = 300   300
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
pbc         = xyz
gen_vel     = no
EOF

cat > md.mdp << EOFMD
integrator  = md
nsteps      = {nsteps}
dt          = 0.002
nstxout-compressed = 5000
nstenergy   = 5000
nstlog      = 5000
continuation = yes
constraint_algorithm = lincs
constraints = h-bonds
cutoff-scheme = Verlet
nstlist     = 10
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1   0.1
ref_t       = 300   300
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
pbc         = xyz
gen_vel     = no
EOFMD

cat > ions.mdp << 'EOF'
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000
nstlist     = 1
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = cutoff
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
EOF

echo "MDP files: em.mdp, nvt.mdp, npt.mdp, md.mdp, ions.mdp"
echo "Ready for SLURM submission."
"""

    write_text_output(script, args.output)


# ---------------------------------------------------------------------------
# CLI definition
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Drug Discovery Wizard CLI"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- prepare-receptor ---
    p_prep = subparsers.add_parser(
        "prepare-receptor",
        help="Fetch and prepare a protein receptor for docking",
    )
    src = p_prep.add_mutually_exclusive_group(required=True)
    src.add_argument("--pdb-id", help="PDB ID to fetch (e.g., 1M17)")
    src.add_argument("--uniprot-id",
                      help="UniProt ID for AlphaFold model (e.g., P00533)")
    src.add_argument("--input-file", help="Local PDB file path")
    p_prep.add_argument("--chain", default=None,
                         help="Chain ID to extract (default: all)")
    p_prep.add_argument("--output", required=True,
                         help="Output JSON file path")
    p_prep.set_defaults(func=cmd_prepare_receptor)

    # --- query-chembl ---
    p_chembl = subparsers.add_parser(
        "query-chembl",
        help="Query ChEMBL for bioactive compounds against a target",
    )
    p_chembl.add_argument("--target-id", required=True,
                           help="ChEMBL target ID (e.g., CHEMBL203)")
    p_chembl.add_argument("--pchembl-min", type=float, required=True,
                           help="Minimum pChEMBL value")
    p_chembl.add_argument("--assay-type", default="B",
                           choices=["B", "F", "A"],
                           help="Assay type: B=binding, F=functional")
    p_chembl.add_argument("--limit", type=int, required=True,
                           help="Maximum results to return")
    p_chembl.add_argument("--output", required=True,
                           help="Output JSON file path")
    p_chembl.set_defaults(func=cmd_query_chembl)

    # --- search-chembl-target ---
    p_search = subparsers.add_parser(
        "search-chembl-target",
        help="Search ChEMBL for targets by name",
    )
    p_search.add_argument("--query", required=True,
                           help="Target name to search")
    p_search.add_argument("--limit", type=int, required=True,
                           help="Maximum results to return")
    p_search.add_argument("--output", required=True,
                           help="Output JSON file path")
    p_search.set_defaults(func=cmd_search_chembl_target)

    # --- admet-filter ---
    p_admet = subparsers.add_parser(
        "admet-filter",
        help="Run ADMET filtering cascade on compounds",
    )
    inp = p_admet.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles", nargs="+", help="SMILES strings")
    inp.add_argument("--csv", help="CSV file with compounds")
    p_admet.add_argument("--names", nargs="+",
                          help="Compound names (with --smiles)")
    p_admet.add_argument("--smiles-col", default="SMILES",
                          help="SMILES column in CSV")
    p_admet.add_argument("--name-col", default="Name",
                          help="Name column in CSV")
    p_admet.add_argument("--sa-threshold", type=float, default=6.0,
                          help="SA score threshold (default: 6.0)")
    p_admet.add_argument("--strictness", default="strict",
                          choices=["strict", "relaxed"],
                          help="'relaxed' = post-docking seeds for REINVENT, "
                               "'strict' = post-REINVENT for downstream (default)")
    p_admet.add_argument("--output", required=True,
                          help="Output JSON file path")
    p_admet.set_defaults(func=cmd_admet_filter)

    # --- dock ---
    p_dock = subparsers.add_parser(
        "dock",
        help="Dock compounds using AutoDock Vina",
    )
    p_dock.add_argument("--receptor", required=True,
                         help="Receptor PDBQT file")
    dinp = p_dock.add_mutually_exclusive_group(required=True)
    dinp.add_argument("--smiles", nargs="+", help="SMILES to dock")
    dinp.add_argument("--csv", help="CSV file with compounds")
    p_dock.add_argument("--names", nargs="+",
                         help="Compound names (with --smiles)")
    p_dock.add_argument("--smiles-col", default="SMILES")
    p_dock.add_argument("--name-col", default="Name")
    site_grp = p_dock.add_mutually_exclusive_group(required=True)
    site_grp.add_argument("--site", help="Binding site JSON config")
    site_grp.add_argument("--center", nargs=3, type=float,
                           metavar=("X", "Y", "Z"))
    p_dock.add_argument("--box-size", nargs=3, type=float,
                         default=[20, 20, 20], metavar=("X", "Y", "Z"))
    p_dock.add_argument("--exhaustiveness", type=int, default=32)
    p_dock.add_argument("--n-poses", type=int, default=9)
    p_dock.add_argument("--output", required=True,
                         help="Output JSON file path")
    p_dock.set_defaults(func=cmd_dock)

    # --- generate-reinvent-config ---
    p_reinv = subparsers.add_parser(
        "generate-reinvent-config",
        help="Generate REINVENT 4 TOML configuration",
    )
    sinp = p_reinv.add_mutually_exclusive_group(required=True)
    sinp.add_argument("--seeds", nargs="+", help="Seed SMILES")
    sinp.add_argument("--seeds-csv", help="CSV with seed compounds")
    p_reinv.add_argument("--smiles-col", default="SMILES")
    p_reinv.add_argument("--target-logp", type=float, default=3.0)
    p_reinv.add_argument("--target-qed", type=float, default=0.6)
    p_reinv.add_argument("--receptor", default=None)
    p_reinv.add_argument("--site", default=None)
    p_reinv.add_argument("--prior", default="reinvent.prior")
    p_reinv.add_argument("--num-steps", type=int, default=500)
    p_reinv.add_argument("--batch-size", type=int, default=128)
    p_reinv.add_argument("--output-dir", default="reinvent_output")
    p_reinv.add_argument("--output", required=True,
                          help="Output TOML config file path")
    p_reinv.set_defaults(func=cmd_generate_reinvent_config)

    # --- generate-hpc-script ---
    p_hpc = subparsers.add_parser(
        "generate-hpc-script",
        help="Generate SLURM submission script",
    )
    p_hpc.add_argument("--job-type", required=True,
                        choices=["reinvent", "gromacs-md", "vina-screen"])
    p_hpc.add_argument("--config", help="REINVENT config file")
    p_hpc.add_argument("--conda-env", default="reinvent4")
    p_hpc.add_argument("--receptor", help="Receptor PDBQT for vina-screen")
    p_hpc.add_argument("--ligand-dir", help="Ligand dir for vina-screen")
    p_hpc.add_argument("--site-config", default="binding_site.json")
    p_hpc.add_argument("--production-ns", type=int, default=100)
    p_hpc.add_argument("--gpu", default="a100")
    p_hpc.add_argument("--ngpu", type=int, default=1)
    p_hpc.add_argument("--ncpus", type=int, default=8)
    p_hpc.add_argument("--mem", default="64G")
    p_hpc.add_argument("--time", default="24:00:00")
    p_hpc.add_argument("--array", default=None)
    p_hpc.add_argument("--output", required=True,
                        help="Output script file path")
    p_hpc.set_defaults(func=cmd_generate_hpc_script)

    # --- generate-md-setup ---
    p_md = subparsers.add_parser(
        "generate-md-setup",
        help="Generate GROMACS MD system setup script",
    )
    p_md.add_argument("--pdb", required=True, help="Input PDB file")
    p_md.add_argument("--ff", default="charmm36", help="Force field")
    p_md.add_argument("--water", default="tip3p", help="Water model")
    p_md.add_argument("--box-dist", type=float, default=1.2,
                       help="Box distance in nm")
    p_md.add_argument("--production-ns", type=int, default=100)
    p_md.add_argument("--output", required=True,
                       help="Output script file path")
    p_md.set_defaults(func=cmd_generate_md_setup)

    # --- admet-predict (ADMETlab 3.0) ---
    p_admet_api = subparsers.add_parser(
        "admet-predict",
        help="Full ADMET prediction via ADMETlab 3.0 API (119 endpoints)",
    )
    p_admet_api.add_argument("--smiles", nargs="+", required=True,
                              help="SMILES strings to predict")
    p_admet_api.add_argument("--names", nargs="+", help="Compound names")
    p_admet_api.add_argument("--output", required=True,
                              help="Output JSON file path")
    p_admet_api.set_defaults(func=cmd_admet_predict)

    # --- gmxapi-setup ---
    p_gmxapi = subparsers.add_parser(
        "gmxapi-setup",
        help="Generate gmxapi Python workflow + MDP files for GROMACS MD",
    )
    p_gmxapi.add_argument("--pdb", required=True, help="Input PDB file")
    p_gmxapi.add_argument("--ff", default="charmm36", help="Force field")
    p_gmxapi.add_argument("--water", default="tip3p", help="Water model")
    p_gmxapi.add_argument("--production-ns", type=float, default=100.0,
                           help="Production run length in ns")
    p_gmxapi.add_argument("--box-dist", type=float, default=1.2,
                           help="Box distance in nm")
    p_gmxapi.add_argument("--temperature", type=float, default=300.0,
                           help="Temperature in K")
    p_gmxapi.add_argument("--output", required=True,
                           help="Output Python workflow script path")
    p_gmxapi.set_defaults(func=cmd_gmxapi_setup)

    # --- gmxapi-validate ---
    p_gval = subparsers.add_parser(
        "gmxapi-validate",
        help="Validate GROMACS MDP parameter files",
    )
    p_gval.add_argument("--mdp", nargs="+", required=True, help="MDP files")
    p_gval.add_argument("--output", required=True, help="Output JSON report")
    p_gval.set_defaults(func=cmd_gmxapi_validate)

    # --- predict-structure (ESMFold) ---
    p_predict = subparsers.add_parser(
        "predict-structure",
        help="De novo protein structure prediction via ESMFold API",
    )
    seq_group = p_predict.add_mutually_exclusive_group(required=True)
    seq_group.add_argument("--sequence", help="Amino acid sequence")
    seq_group.add_argument("--fasta", help="FASTA file path")
    p_predict.add_argument("--output", required=True,
                            help="Output PDB file path")
    p_predict.set_defaults(func=cmd_predict_structure)

    # --- generate-af2-script ---
    p_af2 = subparsers.add_parser(
        "generate-af2-script",
        help="Generate AF2/ColabFold HPC prediction script",
    )
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
    p_af2.set_defaults(func=cmd_generate_af2_script)

    # --- assess-structure ---
    p_assess = subparsers.add_parser(
        "assess-structure",
        help="Assess predicted structure quality via pLDDT",
    )
    p_assess.add_argument("--pdb", required=True, help="PDB file to assess")
    p_assess.add_argument("--output", required=True, help="Output JSON report")
    p_assess.set_defaults(func=cmd_assess_structure)

    # --- query-coconut ---
    p_coco = subparsers.add_parser(
        "query-coconut",
        help="Search COCONUT natural products database (ANPDB, CMNPD, etc.)",
    )
    coco_q = p_coco.add_mutually_exclusive_group(required=True)
    coco_q.add_argument("--query", help="Text search (name, organism, etc.)")
    coco_q.add_argument("--smiles", help="SMILES for structure search")
    p_coco.add_argument("--limit", type=int, default=25)
    p_coco.add_argument("--page", type=int, default=1)
    p_coco.add_argument("--output", required=True, help="Output JSON")
    p_coco.set_defaults(func=cmd_query_coconut)

    # --- query-zinc ---
    p_zinc = subparsers.add_parser(
        "query-zinc",
        help="Generate ZINC20 drug-like compound download script",
    )
    p_zinc.add_argument("--subset", default="drug-like",
                         choices=["drug-like", "lead-like", "fragment-like",
                                  "all-purchasable"])
    p_zinc.add_argument("--mw-range", type=int, nargs=2,
                         metavar=("MIN", "MAX"), help="MW range in Da")
    p_zinc.add_argument("--logp-range", type=float, nargs=2,
                         metavar=("MIN", "MAX"), help="LogP range")
    p_zinc.add_argument("--format", default="smi",
                         choices=["smi", "sdf", "mol2"])
    p_zinc.add_argument("--output", required=True, help="Output script path")
    p_zinc.set_defaults(func=cmd_query_zinc)

    args = parser.parse_args()
    args.func(args)


# ---------------------------------------------------------------------------
# New command handlers (delegate to standalone scripts)
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


def cmd_gmxapi_setup(args):
    """gmxapi workflow generation — delegates to gmxapi_setup.py."""
    from gmxapi_setup import generate_gmxapi_workflow, generate_mdp
    script = generate_gmxapi_workflow(
        pdb=args.pdb,
        ff=args.ff,
        water=args.water,
        production_ns=args.production_ns,
        box_dist=args.box_dist,
        temperature=args.temperature,
    )
    out_dir = Path(args.output).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    for mdp_type in ["ions", "em", "nvt", "npt", "production"]:
        mdp_content = generate_mdp(mdp_type, args.production_ns)
        mdp_name = "md.mdp" if mdp_type == "production" else f"{mdp_type}.mdp"
        (out_dir / mdp_name).write_text(mdp_content)
    Path(args.output).write_text(script)
    print(f"[OK] gmxapi workflow + MDP files saved to {out_dir}")


def cmd_gmxapi_validate(args):
    """MDP validation — delegates to gmxapi_setup.py."""
    from gmxapi_setup import validate_mdp
    results = []
    for mdp_file in args.mdp:
        results.append(validate_mdp(mdp_file))
    report = {
        "total_files": len(results),
        "all_valid": all(r["valid"] for r in results),
        "files": results,
    }
    write_output(report, args.output)
    for r in results:
        status = "VALID" if r["valid"] else "INVALID"
        print(f"  {r['file']}: {status}")


def cmd_predict_structure(args):
    """De novo structure prediction via ESMFold — delegates to predict_structure.py."""
    from predict_structure import predict_esmfold, read_fasta
    if args.fasta:
        header, sequence = read_fasta(args.fasta)
        print(f"[INFO] Read sequence from FASTA: {header}")
    else:
        sequence = args.sequence
    result = predict_esmfold(sequence, args.output)
    json_out = Path(args.output).with_suffix(".json")
    write_output(result, str(json_out))
    print(f"[OK] Prediction metadata saved to {json_out}")


def cmd_generate_af2_script(args):
    """Generate AF2/ColabFold HPC script — delegates to predict_structure.py."""
    from predict_structure import generate_af2_hpc_script
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


def cmd_assess_structure(args):
    """Assess predicted structure quality — delegates to predict_structure.py."""
    from predict_structure import assess_structure
    result = assess_structure(args.pdb)
    write_output(result, args.output)
    print(f"[OK] Quality report: mean pLDDT = {result.get('mean_plddt', 'N/A')}")


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


if __name__ == "__main__":
    main()
