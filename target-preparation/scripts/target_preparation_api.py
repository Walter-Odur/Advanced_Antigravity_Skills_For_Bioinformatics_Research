"""Target Preparation API.\n\nStructure prediction (ESMFold, AF2), receptor fetching (PDB, AlphaFold DB),\nand receptor preparation (cleaning, PDBQT conversion).

Usage:
    python3 scripts/target_preparation_api.py <subcommand> --output <file> [options]
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


# ---------------------------------------------------------------------------
# CLI definition
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Target Preparation — structure prediction & receptor setup"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # --- prepare-receptor ---
    p_prep = sub.add_parser(
        "prepare-receptor",
        help="Fetch and prepare a receptor for docking",
    )
    inp = p_prep.add_mutually_exclusive_group()
    inp.add_argument("--pdb-id", help="PDB ID to fetch")
    inp.add_argument("--uniprot-id", help="UniProt ID for AlphaFold fetch")
    inp.add_argument("--input-file", help="Local PDB file path")
    p_prep.add_argument("--chain", default=None,
                         help="Chain ID to extract (default: all)")
    p_prep.add_argument("--output", required=True,
                         help="Output JSON file path")
    p_prep.set_defaults(func=cmd_prepare_receptor)

    # --- predict-structure ---
    p_predict = sub.add_parser(
        "predict-structure",
        help="Predict protein structure from sequence (ESMFold)",
    )
    p_predict.add_argument("--sequence", help="Amino acid sequence")
    p_predict.add_argument("--fasta", help="FASTA file path")
    p_predict.add_argument("--output", required=True,
                            help="Output PDB file path")
    p_predict.set_defaults(func=cmd_predict_structure)

    # --- generate-af2-script ---
    p_af2 = sub.add_parser(
        "generate-af2-script",
        help="Generate AlphaFold2/ColabFold HPC submission script",
    )
    p_af2.add_argument("--fasta", required=True, help="FASTA file")
    p_af2.add_argument("--method", default="colabfold",
                        choices=["colabfold", "alphafold2"])
    p_af2.add_argument("--gpu", default="a100")
    p_af2.add_argument("--ngpu", type=int, default=1)
    p_af2.add_argument("--mem", default="64G")
    p_af2.add_argument("--time", default="12:00:00")
    p_af2.add_argument("--db-path", default="/data/colabfold_dbs")
    p_af2.add_argument("--conda-env", default="colabfold")
    p_af2.add_argument("--output", required=True)
    p_af2.set_defaults(func=cmd_generate_af2_script)

    # --- assess-structure ---
    p_assess = sub.add_parser(
        "assess-structure",
        help="Assess predicted structure quality (pLDDT)",
    )
    p_assess.add_argument("--pdb", required=True, help="PDB file to assess")
    p_assess.add_argument("--output", required=True,
                           help="Output JSON file path")
    p_assess.set_defaults(func=cmd_assess_structure)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
