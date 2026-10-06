"""CLI for the md-simulation skill."""

from __future__ import annotations

import argparse
from pathlib import Path

from ..hpc import JOB_TYPES, SlurmResources, generate_slurm_script
from ..io_utils import write_json, write_text
from ..md.mdp import FORCE_FIELDS, WATER_MODELS, validate_mdp_file
from ..md.workflow import (
    generate_gmxapi_workflow,
    generate_setup_shell,
    workflow_bundle,
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
Generate validated GROMACS workflows and cluster submission scripts.

This skill produces files; it does not run simulations. Production MD is
GPU work that belongs in a batch job.

Non-bonded settings are derived from the force field you choose, because
they are not interchangeable: CHARMM36 requires a force-switched van der
Waals potential and no dispersion correction, while the AMBER family
requires a plain cutoff with one. Mixing them changes the physics.
"""


def cmd_gmxapi_setup(args: argparse.Namespace) -> None:
    bundle = workflow_bundle(
        pdb=args.pdb, force_field=args.ff, water=args.water,
        temperature=args.temperature, pressure=args.pressure,
        production_ns=args.production_ns,
        equilibration_ps=args.equilibration_ps, dt_ps=args.dt,
        box_distance=args.box_dist, box_type=args.box_type,
        ion_concentration=args.ion_concentration,
        has_ligand=args.has_ligand, seed=args.seed,
    )
    out_path = Path(args.output)
    out_dir = out_path.parent
    written: list[str] = []

    for filename, content in bundle["mdp_files"].items():
        written.append(str(write_text(content, out_dir / filename)))
    written.append(str(write_text(generate_gmxapi_workflow(bundle), out_path)))
    shell = out_dir / "setup_md.sh"
    written.append(str(write_text(generate_setup_shell(bundle), shell,
                                  executable=True)))

    report = {**bundle, "written_files": written,
               "workflow_script": str(out_path),
               "shell_script": str(shell)}
    report.pop("mdp_files")
    metadata = out_path.with_suffix(".meta.json")
    write_json(report, metadata)

    settings = bundle["settings"]
    lines = [
        f"MD setup for {args.pdb}",
        f"  force field: {settings['force_field']} "
        f"({settings['gromacs_force_field']}), water {settings['water_model']}",
        f"  temperature: {settings['temperature_K']:g} K, "
        f"pressure {settings['pressure_bar']:g} bar",
        f"  production: {settings['production_ns']:g} ns at "
        f"{settings['timestep_ps']:g} ps/step",
        f"  {bundle['force_field_note']}",
        f"  MDP files: {', '.join(sorted(bundle['stages'][i]['file'] for i in range(len(bundle['stages']))))}",
        f"  shell script: {shell}",
    ]
    for warning in bundle["warnings"]:
        lines.append(f"  warning: {warning}")
    lines.append("  next: bash setup_md.sh, then submit the production run")
    if not args.quiet:
        for line in lines:
            print(line)
        print(f"Written: {out_path}")
        print(f"Written: {metadata}")


def cmd_gmxapi_validate(args: argparse.Namespace) -> None:
    reports = [validate_mdp_file(path, force_field=args.ff)
               for path in args.mdp]
    all_valid = all(r["valid"] for r in reports)
    payload = {
        "all_valid": all_valid,
        "force_field_checked": args.ff,
        "n_files": len(reports),
        "n_invalid": sum(1 for r in reports if not r["valid"]),
        "reports": reports,
        "note": ("An MDP file does not record its intended force field, so "
                 "pass --ff to have the non-bonded settings cross-checked "
                 "against it."),
    }
    lines = [f"MDP validation: {len(reports)} file(s), "
             f"{'all valid' if all_valid else 'PROBLEMS FOUND'}"]
    for report in reports:
        lines.append(f"  {Path(report['file']).name}: "
                     f"{'valid' if report['valid'] else 'INVALID'} "
                     f"({report['n_parameters']} parameters)")
        for error in report["errors"]:
            lines.append(f"      ERROR: {error}")
        for warning in report["warnings"]:
            lines.append(f"      warning: {warning}")
        for note in report["notes"]:
            lines.append(f"      {note}")
    emit_json(payload, args, summary=lines)


def cmd_generate_md_setup(args: argparse.Namespace) -> None:
    bundle = workflow_bundle(
        pdb=args.pdb, force_field=args.ff, water=args.water,
        temperature=args.temperature, pressure=args.pressure,
        production_ns=args.production_ns,
        equilibration_ps=args.equilibration_ps, dt_ps=args.dt,
        box_distance=args.box_dist, box_type=args.box_type,
        ion_concentration=args.ion_concentration,
        has_ligand=args.has_ligand, seed=args.seed,
    )
    out_dir = Path(args.output).parent
    written = [str(write_text(content, out_dir / filename))
               for filename, content in bundle["mdp_files"].items()]
    emit_text(generate_setup_shell(bundle), args, executable=True, summary=[
        f"GROMACS setup script for {args.pdb}",
        f"  force field {bundle['settings']['force_field']}, "
        f"{bundle['settings']['temperature_K']:g} K",
        f"  MDP files written: {len(written)}",
        *(f"  warning: {w}" for w in bundle["warnings"]),
        f"  run with: bash {args.output}",
    ])


def cmd_generate_hpc_script(args: argparse.Namespace) -> None:
    resources = SlurmResources(
        job_name=args.job_name or args.job_type, gpu=args.gpu, ngpu=args.ngpu,
        ncpus=args.ncpus, mem=args.mem, time=args.time, array=args.array,
        conda_env=args.conda_env, partition=args.partition,
        account=args.account, email=args.email,
    )
    settings = {
        "config": args.config,
        "receptor": args.receptor,
        "ligand_dir": args.ligand_dir,
        "site_config": args.site_config,
        "production_ns": args.production_ns,
        "tpr": args.tpr,
        "deffnm": args.deffnm,
        "exhaustiveness": args.exhaustiveness,
        "fasta": args.fasta,
    }
    script = generate_slurm_script(args.job_type, resources,
                                   {k: v for k, v in settings.items()
                                    if v is not None})
    emit_text(script, args, executable=True, summary=[
        f"SLURM script: {args.job_type} - {JOB_TYPES[args.job_type]}",
        f"  {args.ngpu}x{args.gpu or 'gpu'}, {args.ncpus} CPU, {args.mem}, "
        f"{args.time}" + (f", array {args.array}" if args.array else ""),
        f"  submit with: sbatch {args.output}",
    ])


def _add_system_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--pdb", required=True, help="Cleaned input PDB.")
    parser.add_argument("--ff", default="charmm36", choices=sorted(FORCE_FIELDS),
                        help="Force field (default: charmm36). Determines the "
                             "non-bonded settings.")
    parser.add_argument("--water", default="tip3p", choices=sorted(WATER_MODELS),
                        help="Water model (default: tip3p).")
    parser.add_argument("--temperature", type=float, default=300.0,
                        help="Reference temperature in K (default: 300). "
                             "Body temperature is 310.")
    parser.add_argument("--pressure", type=float, default=1.0,
                        help="Reference pressure in bar (default: 1).")
    parser.add_argument("--production-ns", type=float, default=100.0,
                        help="Production length in ns (default: 100). "
                             "0 generates equilibration only.")
    parser.add_argument("--equilibration-ps", type=float, default=100.0,
                        help="Length of each of NVT and NPT, in ps.")
    parser.add_argument("--dt", type=float, default=0.002,
                        help="Timestep in ps (default: 0.002).")
    parser.add_argument("--box-dist", type=float, default=1.2,
                        help="Solute-to-box distance in nm (default: 1.2).")
    parser.add_argument("--box-type", default="dodecahedron",
                        choices=["cubic", "dodecahedron", "octahedron",
                                 "triclinic"],
                        help="A dodecahedron needs ~29%% fewer waters than a "
                             "cube for the same clearance.")
    parser.add_argument("--ion-concentration", type=float, default=0.15,
                        help="Salt concentration in M (default: 0.15, "
                             "physiological).")
    parser.add_argument("--has-ligand", action="store_true",
                        help="Couple a ligand with the protein's temperature "
                             "group.")
    parser.add_argument("--seed", type=int, default=-1,
                        help="Velocity seed; -1 is random, a positive value "
                             "makes the run reproducible.")


def build_parser() -> argparse.ArgumentParser:
    parser, sub = build_root_parser("ag-md-simulation", DESCRIPTION)

    p = sub.add_parser("gmxapi-setup",
                       help="Generate a gmxapi workflow, a shell script and "
                            "all MDP files.")
    _add_system_args(p)
    add_output(p, help_text="Output workflow .py path. MDP files and "
                            "setup_md.sh are written beside it.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_gmxapi_setup)

    p = sub.add_parser("generate-md-setup",
                       help="Generate the GROMACS setup shell script "
                            "and MDP files.")
    _add_system_args(p)
    add_output(p, help_text="Output .sh script.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_generate_md_setup)

    p = sub.add_parser("gmxapi-validate",
                       help="Validate MDP files, optionally against a "
                            "force field.")
    p.add_argument("--mdp", nargs="+", required=True, help="MDP files.")
    p.add_argument("--ff", choices=sorted(FORCE_FIELDS),
                   help="Cross-check the non-bonded settings against this "
                        "force field. Strongly recommended: the file itself "
                        "does not record which force field it is for.")
    add_output(p, help_text="Output JSON validation report.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_gmxapi_validate)

    p = sub.add_parser("generate-hpc-script",
                       help="Generate a SLURM submission script.")
    p.add_argument("--job-type", required=True, choices=sorted(JOB_TYPES),
                   help="; ".join(f"{k}: {v}" for k, v in JOB_TYPES.items()))
    p.add_argument("--job-name")
    p.add_argument("--config", help="REINVENT TOML, for --job-type reinvent.")
    p.add_argument("--tpr", help="GROMACS .tpr, for gromacs-md.")
    p.add_argument("--deffnm", help="GROMACS output basename.")
    p.add_argument("--receptor", help="Receptor PDBQT, for vina-screen.")
    p.add_argument("--ligand-dir", help="Ligand directory, for vina-screen.")
    p.add_argument("--site-config", help="Binding-site JSON, for vina-screen.")
    p.add_argument("--exhaustiveness", type=int, help="Vina exhaustiveness.")
    p.add_argument("--fasta", help="FASTA, for colabfold or alphafold2.")
    p.add_argument("--production-ns", type=int, default=100)
    p.add_argument("--gpu", default="a100")
    p.add_argument("--ngpu", type=int, default=1)
    p.add_argument("--ncpus", type=int, default=8)
    p.add_argument("--mem", default="64G")
    p.add_argument("--time", default="24:00:00")
    p.add_argument("--array", help="SLURM array spec, e.g. 1-100%%10. "
                                   "Required for vina-screen.")
    p.add_argument("--conda-env")
    p.add_argument("--partition")
    p.add_argument("--account")
    p.add_argument("--email")
    add_output(p, help_text="Output .sh script.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_generate_hpc_script)

    return parser


def main(argv: list[str] | None = None) -> int:
    return run_cli(build_parser(), argv)


if __name__ == "__main__":
    raise SystemExit(main())
