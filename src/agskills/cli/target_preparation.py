"""CLI for the target-preparation skill."""

from __future__ import annotations

import argparse
from pathlib import Path

from ..errors import UsageError
from ..io_utils import write_json
from ..targets.clean import clean_structure
from ..targets.esmfold import predict_structure, read_fasta, validate_sequence
from ..targets.fetch import fetch_alphafold, fetch_pdb
from ..targets.pdbqt import available_converters, receptor_to_pdbqt
from ..targets.quality import assess_structure
from ..targets.site import (
    list_ligands,
    site_from_hetatm,
    site_from_residues,
    site_from_structure,
)
from ._common import (
    add_common_flags,
    add_output,
    build_root_parser,
    emit_json,
    emit_text,
    run_cli,
)

DESCRIPTION = """\
Prepare a protein target for structure-based drug discovery.

Typical sequence:
  1. prepare-receptor   fetch a structure, clean it, convert to PDBQT
  2. define-site        produce the docking box the screening skill needs
  3. assess-structure   check a predicted model before trusting it
"""


# ---------------------------------------------------------------------------
# Handlers
# ---------------------------------------------------------------------------


def cmd_prepare_receptor(args: argparse.Namespace) -> None:
    out_dir = Path(args.output).parent
    sources = [bool(args.pdb_id), bool(args.uniprot_id), bool(args.input_file)]
    if sum(sources) != 1:
        raise UsageError(
            "Provide exactly one of --pdb-id, --uniprot-id or --input-file."
        )

    result: dict = {"steps": []}
    if args.pdb_id:
        download = fetch_pdb(args.pdb_id, out_dir)
        result["source"] = download.as_dict()
        structure = download.path
    elif args.uniprot_id:
        download = fetch_alphafold(args.uniprot_id, out_dir)
        result["source"] = download.as_dict()
        structure = download.path
    else:
        structure = Path(args.input_file)
        result["source"] = {"path": str(structure), "source": "local file"}
    result["steps"].append(f"obtained structure {structure.name}")

    cleaned_path = out_dir / f"{structure.stem}_clean.pdb"
    clean = clean_structure(
        structure, cleaned_path, chain=args.chain,
        keep_ligands=args.keep_ligand, keep_metals=not args.remove_metals,
        keep_waters=args.keep_waters, remove_hydrogens=args.remove_hydrogens,
    )
    result["cleaning"] = clean.as_dict()
    result["steps"].append(
        f"cleaned to {cleaned_path.name}: kept {clean.atoms_kept} atoms in "
        f"chain(s) {','.join(clean.chains_kept)}"
    )

    quality = assess_structure(structure)
    result["quality"] = quality.as_dict()

    # Derive the site before converting, so the converter can report
    # whether any residue it had to approximate lies inside the box.
    ligands = list_ligands(structure)
    result["candidate_ligands"] = ligands
    site = None
    if ligands:
        try:
            site = site_from_hetatm(structure)
            site_payload = site.as_dict()
            site_file = out_dir / "binding_site.json"
            write_json(site_payload, site_file)
            result["binding_site"] = site_payload
            result["binding_site_file"] = str(site_file)
            result["steps"].append(
                f"binding site defined from {site.reference}"
            )
        except Exception as exc:  # non-fatal: the receptor is still usable
            result["binding_site_error"] = str(exc)

    conversion = receptor_to_pdbqt(cleaned_path, site=site)
    result["pdbqt"] = conversion.as_dict()
    result["pdbqt_converters_available"] = available_converters()
    if conversion.ok:
        result["steps"].append(
            f"converted to PDBQT with {conversion.converter}"
            + ("" if conversion.strict else " (permissive mode)")
        )
        result["receptor_for_docking"] = conversion.output_path
    else:
        result["steps"].append("PDBQT conversion unavailable")
        result["receptor_for_docking"] = None

    result["next_step"] = (
        "Dock against it: ag-compound-screening dock --receptor "
        f"{Path(conversion.output_path).name if conversion.output_path else '<receptor.pdbqt>'} "
        f"--site binding_site.json --smiles ..."
        if site and conversion.ok else
        "No bound ligand was found, so define the docking box from pocket "
        "residues: define-site --pdb <file> --residues A:123 A:124 ..."
        if not ligands else
        "Install a PDBQT converter (pip install meeko scipy gemmi) to "
        "produce a docking-ready receptor."
    )

    emit_json(result, args, summary=[
        f"Receptor: {result['source'].get('identifier', structure.stem)}",
        f"  cleaned PDB: {cleaned_path}",
        f"  atoms kept: {clean.atoms_kept}, waters removed: {clean.waters_removed}",
        f"  heterogens removed: {len(clean.heterogens_removed)}",
        f"  PDBQT: {conversion.output_path or 'not produced - ' + conversion.converter}",
        *(f"  pdbqt warning: {w}" for w in conversion.warnings),
        f"  confidence verdict: {quality.verdict}",
        f"  bound ligands: {len(ligands)}"
        + (f"; site from {site.reference}" if site else ""),
        f"  next: {result['next_step']}",
    ])


def cmd_define_site(args: argparse.Namespace) -> None:
    selectors = [args.from_ligand, bool(args.residues), args.whole_structure]
    if sum(bool(s) for s in selectors) != 1:
        raise UsageError(
            "Choose exactly one of --from-ligand, --residues or "
            "--whole-structure.",
            hint="--from-ligand is the most reliable when the structure has "
                 "a bound ligand.",
        )
    if args.from_ligand:
        site = site_from_hetatm(args.pdb, resname=args.ligand_resname,
                                chain=args.ligand_chain,
                                resseq=args.ligand_resseq,
                                padding=args.padding)
    elif args.residues:
        site = site_from_residues(args.pdb, args.residues, padding=args.padding)
    else:
        site = site_from_structure(args.pdb, padding=args.padding)

    payload = site.as_dict()
    payload["pdb"] = str(args.pdb)
    payload["usage"] = (
        "Pass this file to docking: ag-compound-screening dock "
        f"--receptor <receptor.pdbqt> --site {args.output} --smiles ..."
    )
    emit_json(payload, args, summary=[
        f"Binding site from {site.method}",
        f"  reference: {site.reference}",
        f"  centre: ({site.center_x:.2f}, {site.center_y:.2f}, {site.center_z:.2f})",
        f"  box: {site.size_x:.1f} x {site.size_y:.1f} x {site.size_z:.1f} A "
        f"= {site.volume:.0f} A^3",
        f"  confidence: {site.confidence}",
        *(f"  warning: {w}" for w in site.warnings),
    ])


def cmd_list_ligands(args: argparse.Namespace) -> None:
    ligands = list_ligands(args.pdb)
    emit_json({"pdb": str(args.pdb), "n_candidates": len(ligands),
               "ligands": ligands,
               "note": "Solvent, buffer additives and lone ions are excluded. "
                       "The largest heterogen is usually the biologically "
                       "relevant ligand."},
              args,
              summary=[f"Candidate ligands in {Path(args.pdb).name}: {len(ligands)}"]
                      + [f"  {entry['label']}  ({entry['n_atoms']} atoms)"
                         for entry in ligands[:10]])


def cmd_assess_structure(args: argparse.Namespace) -> None:
    report = assess_structure(args.pdb)
    emit_json(report.as_dict(), args, summary=[
        f"Structure: {Path(args.pdb).name}",
        f"  residues: {report.n_residues} in chain(s) {','.join(report.chains)}",
        f"  mean B-factor/pLDDT: {report.mean_plddt:.1f} "
        f"(median {report.median_plddt:.1f})",
        f"  bands: {report.residues_very_high} very high, "
        f"{report.residues_confident} confident, {report.residues_low} low, "
        f"{report.residues_very_low} very low",
        f"  VERDICT: {report.verdict} - {report.reason}",
        *(f"  action: {a}" for a in report.actions),
        f"  note: {report.b_factor_caveat}",
    ])


def cmd_predict_structure(args: argparse.Namespace) -> None:
    if bool(args.sequence) == bool(args.fasta):
        raise UsageError("Provide exactly one of --sequence or --fasta.")
    if args.fasta:
        info = read_fasta(args.fasta)
        sequence = info.sequence
    else:
        sequence = args.sequence
    result = predict_structure(sequence, args.output)
    # The PDB goes to --output; its metadata goes alongside it.
    metadata_path = Path(args.output).with_suffix(".json")
    write_json(result, metadata_path)
    quality = result["quality"]
    if not args.quiet:
        print(f"Predicted structure: {args.output}")
        print(f"  method: {result['method']}")
        print(f"  residues: {result['sequence_length']}")
        print(f"  mean pLDDT: {quality['mean_plddt']}")
        print(f"  VERDICT: {quality['verdict']} - {quality['reason']}")
        for action in quality["actions"]:
            print(f"  action: {action}")
        print(f"Written: {args.output}")
        print(f"Written: {metadata_path}")


def cmd_validate_sequence(args: argparse.Namespace) -> None:
    if bool(args.sequence) == bool(args.fasta):
        raise UsageError("Provide exactly one of --sequence or --fasta.")
    if args.fasta:
        info = read_fasta(args.fasta)
        validated = validate_sequence(info.sequence)
        warnings = info.warnings
        header = info.header
    else:
        validated = validate_sequence(args.sequence)
        warnings, header = [], ""
    emit_json({"valid": True, "length": validated.length, "header": header,
               "sequence": validated.sequence, "warnings": warnings},
              args,
              summary=[f"Sequence is valid: {validated.length} residues"])


def cmd_generate_af2_script(args: argparse.Namespace) -> None:
    from ..hpc import SlurmResources, generate_slurm_script
    resources = SlurmResources(
        job_name=args.job_name, gpu=args.gpu, ngpu=args.ngpu,
        ncpus=args.ncpus, mem=args.mem, time=args.time,
        conda_env=args.conda_env, partition=args.partition,
        account=args.account,
    )
    script = generate_slurm_script(
        args.method, resources,
        {"fasta": args.fasta, "db_path": args.db_path,
         "out_dir": args.out_dir},
    )
    emit_text(script, args, executable=True, summary=[
        f"{args.method} submission script for {args.fasta}",
        f"  resources: {args.ngpu}x{args.gpu}, {args.mem}, {args.time}",
        f"  submit with: sbatch {args.output}",
    ])


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser, sub = build_root_parser(
        "ag-target-preparation", DESCRIPTION,
        epilog="Every subcommand requires --output.",
    )

    p = sub.add_parser("prepare-receptor",
                       help="Fetch, clean and convert a receptor for docking.")
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--pdb-id", help="RCSB PDB entry id, e.g. 6HEZ.")
    source.add_argument("--uniprot-id",
                        help="UniProt accession for an AlphaFold model, "
                             "e.g. P9WJA5.")
    source.add_argument("--input-file", help="A local PDB file.")
    p.add_argument("--chain", help="Keep only this chain id.")
    p.add_argument("--keep-ligand", nargs="+", default=None, metavar="RESNAME",
                   help="Residue names to retain, e.g. a cofactor (FAD, NAD).")
    p.add_argument("--keep-waters", action="store_true",
                   help="Retain water molecules (off by default).")
    p.add_argument("--remove-metals", action="store_true",
                   help="Also remove metal ions, which are kept by default "
                        "because they are often catalytic.")
    p.add_argument("--remove-hydrogens", action="store_true",
                   help="Strip hydrogens so the docking tool adds its own.")
    add_output(p, help_text="Output JSON report. Structure files are written "
                            "next to it.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_prepare_receptor)

    p = sub.add_parser("define-site",
                       help="Produce the docking box JSON that 'dock' needs.")
    p.add_argument("--pdb", required=True, help="Structure to analyse.")
    p.add_argument("--from-ligand", action="store_true",
                   help="Centre the box on a co-crystallised ligand "
                        "(most reliable).")
    p.add_argument("--ligand-resname", help="Select a specific ligand by "
                                            "residue name.")
    p.add_argument("--ligand-chain", help="Select a ligand by chain.")
    p.add_argument("--ligand-resseq", type=int,
                   help="Select a ligand by residue number.")
    p.add_argument("--residues", nargs="+", metavar="CHAIN:RESSEQ",
                   help="Centre the box on these pocket residues, "
                        "e.g. A:230 A:231.")
    p.add_argument("--whole-structure", action="store_true",
                   help="Enclose the whole protein for blind docking "
                        "(low confidence).")
    p.add_argument("--padding", type=float, default=4.0,
                   help="Angstroms added around the reference (default: 4).")
    add_output(p, help_text="Output binding-site JSON.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_define_site)

    p = sub.add_parser("list-ligands",
                       help="List candidate co-crystallised ligands.")
    p.add_argument("--pdb", required=True)
    add_output(p, help_text="Output JSON list.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_list_ligands)

    p = sub.add_parser("assess-structure",
                       help="Assess per-residue confidence (pLDDT) and advise.")
    p.add_argument("--pdb", required=True, help="Structure to assess.")
    add_output(p, help_text="Output JSON quality report.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_assess_structure)

    p = sub.add_parser("predict-structure",
                       help="Predict a structure from sequence with ESMFold.")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--sequence", help="Amino-acid sequence.")
    group.add_argument("--fasta", help="FASTA file.")
    add_output(p, help_text="Output PDB path. Metadata is written alongside "
                            "it as .json.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_predict_structure)

    p = sub.add_parser("validate-sequence",
                       help="Check a sequence is foldable before submitting it.")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--sequence")
    group.add_argument("--fasta")
    add_output(p, help_text="Output JSON report.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_validate_sequence)

    p = sub.add_parser("generate-af2-script",
                       help="Generate a ColabFold or AlphaFold2 SLURM script.")
    p.add_argument("--fasta", required=True)
    p.add_argument("--method", default="colabfold",
                   choices=["colabfold", "alphafold2"])
    p.add_argument("--out-dir", default="predictions",
                   help="Where the job writes predictions.")
    p.add_argument("--db-path", default="/data/colabfold_dbs",
                   help="Sequence database directory on the cluster.")
    p.add_argument("--job-name", default="fold")
    p.add_argument("--gpu", default="a100")
    p.add_argument("--ngpu", type=int, default=1)
    p.add_argument("--ncpus", type=int, default=8)
    p.add_argument("--mem", default="64G")
    p.add_argument("--time", default="12:00:00")
    p.add_argument("--conda-env", default="colabfold")
    p.add_argument("--partition")
    p.add_argument("--account")
    add_output(p, help_text="Output .sh script.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_generate_af2_script)

    return parser


def main(argv: list[str] | None = None) -> int:
    return run_cli(build_parser(), argv)


if __name__ == "__main__":
    raise SystemExit(main())
