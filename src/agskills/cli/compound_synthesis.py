"""CLI for the compound-synthesis (REINVENT 4) skill."""

from __future__ import annotations

import argparse
from pathlib import Path

from ..errors import UsageError
from ..hpc import SlurmResources, generate_slurm_script
from ..io_utils import write_json
from ..reinvent.analyze import analyze_results
from ..reinvent.config import RUN_MODES, build_config
from ..reinvent.priors import GENERATORS, check_setup
from ..reinvent.runner import prepare_seeds, preflight, run_reinvent
from ..reinvent.scoring import AGGREGATION_TYPES, SCORING_PROFILES
from ._common import (
    add_common_flags,
    add_compound_inputs,
    add_output,
    build_root_parser,
    compounds_from_args,
    emit_json,
    emit_text,
    run_cli,
)

DESCRIPTION = """\
Generate and optimise molecules with REINVENT 4.

Typical sequence:
  1. check-setup        confirm the install, priors and device
  2. prepare-seeds      validate and de-duplicate input structures
  3. generate-config    write a validated TOML for the run you want
  4. run  /  generate-hpc-script
  5. analyze-results    rank the output and diagnose the run

REINVENT 4 is located via --reinvent-dir, then $REINVENT_DIR, then the
installed package. Prior models are located via $REINVENT_PRIOR_BASE or
<REINVENT4>/priors. No path is hardcoded.
"""


def cmd_check_setup(args: argparse.Namespace) -> None:
    report = check_setup(args.reinvent_dir)
    lines = [
        f"REINVENT 4 setup: {'READY' if report['ready'] else 'NOT READY'}",
        f"  checkout: {report['reinvent_dir'] or 'not found'}",
        f"  version: {report['reinvent_version'] or 'unknown'}",
        f"  importable: {report['reinvent_importable']}  "
        f"executable: {report['reinvent_executable'] or 'not on PATH'}",
        f"  priors: {len(report['priors_found'])} found, "
        f"{len(report['priors_missing'])} missing "
        f"(in {report['priors_dir'] or 'no priors directory'})",
        f"  device: {report['recommended_device']}",
    ]
    for problem in report["problems"]:
        lines.append(f"  PROBLEM: {problem}")
    for advice in report["advice"]:
        lines.append(f"  -> {advice}")
    emit_json(report, args, summary=lines)


def cmd_prepare_seeds(args: argparse.Namespace) -> None:
    records = compounds_from_args(args)
    result = prepare_seeds(records, args.output, generator=args.generator,
                           include_names=args.include_names)
    lines = [
        f"Seed preparation for {args.generator}",
        f"  input {result['input_count']}  written {result['written']}  "
        f"rejected {result['rejected']}",
    ]
    for entry in result["rejected_detail"][:8]:
        lines.append(f"    dropped {entry['name']}: {entry['reason']}")
    for warning in result["warnings"]:
        lines.append(f"  warning: {warning}")
    if args.report:
        write_json(result, args.report)
        lines.append(f"  report: {args.report}")
    if not args.quiet:
        for line in lines:
            print(line)
        print(f"Written: {result['output_file']}")


def cmd_generate_config(args: argparse.Namespace) -> None:
    result = build_config(
        mode=args.mode, generator=args.generator, prior=args.prior,
        device=args.device, scoring_profile=args.scoring_profile,
        components=args.component, aggregation=args.aggregation,
        smiles_file=args.smiles_file,
        validation_smiles_file=args.validation_smiles_file,
        model_file=args.model_file, num_steps=args.num_steps,
        max_steps=args.max_steps, min_steps=args.min_steps,
        max_score=args.max_score, num_epochs=args.num_epochs,
        num_smiles=args.num_smiles, batch_size=args.batch_size,
        sigma=args.sigma, learning_rate=args.learning_rate,
        csv_prefix=args.csv_prefix, output_csv=args.output_csv,
        output_model=args.output_model, checkpoint_file=args.checkpoint_file,
        diversity_filter=not args.no_diversity_filter,
        inception_smiles_file=args.inception_smiles_file,
        reinvent_dir=args.reinvent_dir, sample_strategy=args.sample_strategy,
    )
    lines = [
        f"REINVENT config: {result.run_type} with the {result.generator} "
        "generator",
        *(f"  {k}: {v}" for k, v in result.settings.items()),
        f"  expected outputs: {', '.join(result.expected_outputs) or 'none'}",
    ]
    for warning in result.warnings:
        lines.append(f"  warning: {warning}")
    for step in result.next_steps:
        lines.append(f"  next: {step}")
    emit_text(result.toml, args, summary=lines)
    metadata = Path(args.output).with_suffix(".meta.json")
    write_json(result.as_dict(), metadata)
    if not args.quiet:
        print(f"Written: {metadata}")


def cmd_list_profiles(args: argparse.Namespace) -> None:
    payload = {
        "profiles": {
            name: {
                "description": spec["description"],
                "aggregation": spec["aggregation"],
                "n_components": len(spec["components"]),
                "components": [
                    {
                        "component": c["component"],
                        "name": c["endpoint"]["name"],
                        "weight": c["endpoint"].get("weight", 1.0),
                        "transform": c["endpoint"].get("transform", {}).get("type"),
                    }
                    for c in spec["components"]
                ],
                **({"caveat": spec["caveat"]} if spec.get("caveat") else {}),
            }
            for name, spec in SCORING_PROFILES.items()
        },
        "aggregation_types": AGGREGATION_TYPES,
        "generators": GENERATORS,
        "run_modes": RUN_MODES,
    }
    lines = []
    for name, spec in SCORING_PROFILES.items():
        lines.append(f"{name} ({len(spec['components'])} components, "
                     f"{spec['aggregation']})")
        lines.append(f"  {spec['description']}")
        if spec.get("caveat"):
            lines.append(f"  CAVEAT: {spec['caveat']}")
    emit_json(payload, args, summary=lines)


def cmd_preflight(args: argparse.Namespace) -> None:
    report = preflight(args.config, device=args.device,
                       reinvent_dir=args.reinvent_dir)
    emit_json(report, args, summary=[
        f"Preflight for {args.config}: "
        f"{'READY' if report['ready'] else 'NOT READY'}",
        f"  run_type: {report['run_type']}",
        f"  device: {report['device_requested']} - {report['device_note']}",
        f"  files present: {len(report['files_present'])}",
        f"  reinvent executable: "
        f"{report['reinvent_executable'] or 'not on PATH'}",
    ])


def cmd_run(args: argparse.Namespace) -> None:
    result = run_reinvent(args.config, log_path=args.log, device=args.device,
                          reinvent_dir=args.reinvent_dir,
                          timeout=args.timeout, dry_run=args.dry_run)
    if result.get("dry_run"):
        emit_json(result, args, summary=[
            "Dry run: preflight passed, nothing was executed.",
            f"  would run: {' '.join(result['would_run'])}",
        ])
        return
    emit_json(result, args, summary=[
        f"REINVENT run: {result['status']}",
        f"  elapsed: {result['elapsed_seconds']}s",
        f"  outputs: {', '.join(Path(f).name for f in result['output_files']) or 'none'}",
        f"  {result.get('message', '')}",
    ])


def cmd_analyze_results(args: argparse.Namespace) -> None:
    result = analyze_results(args.csv, top_n=args.top_n, sort_by=args.sort_by,
                             min_score=args.min_score,
                             compute_properties=not args.no_properties)
    statistics = result["score_statistics"]
    lines = [
        f"REINVENT results: {Path(result['file']).name}",
        f"  rows {result['rows_read']}  unique molecules "
        f"{result['unique_molecules']}  duplicate samples "
        f"{result['duplicate_samples']}",
        f"  score: mean {statistics['mean']} median {statistics['median']} "
        f"max {statistics['max']}",
        f"  scaffolds: {result['distinct_scaffolds']} distinct "
        f"({result['scaffold_diversity_ratio']})",
    ]
    if result["interpretation"]["learning"]:
        lines.append(f"  learning: {result['interpretation']['learning']}")
    if result["interpretation"]["diversity"]:
        lines.append(f"  diversity: {result['interpretation']['diversity']}")
    lines.append(f"  top {min(args.top_n, len(result['top_molecules']))}:")
    for rank, molecule in enumerate(result["top_molecules"][:args.show], 1):
        properties = molecule.get("properties", {})
        lines.append(
            f"    {rank:>3d}. score {molecule['score']} "
            f"(x{molecule['times_sampled']}) "
            f"MW {properties.get('mw', '?')} QED {properties.get('qed', '?')} "
            f"{molecule['smiles'][:44]}"
        )
    emit_json(result, args, summary=lines)


def cmd_score_molecules(args: argparse.Namespace) -> None:
    """Score molecules with REINVENT's own scoring framework."""
    records = compounds_from_args(args)
    work_dir = Path(args.output).parent
    smiles_path = work_dir / "score_input.smi"
    from ..io_utils import write_text
    write_text("\n".join(r.smiles for r in records) + "\n", smiles_path)

    config = build_config(
        mode="scoring", generator="reinvent",
        scoring_profile=args.scoring_profile, components=args.component,
        aggregation=args.aggregation, smiles_file=str(smiles_path),
        output_csv=str(work_dir / "scoring.csv"), device=args.device,
        reinvent_dir=args.reinvent_dir,
    )
    config_path = work_dir / "scoring_config.toml"
    write_text(config.toml, config_path)

    result = run_reinvent(config_path, device=args.device,
                          reinvent_dir=args.reinvent_dir,
                          timeout=args.timeout, dry_run=args.dry_run)
    payload = {"config": str(config_path), "smiles_file": str(smiles_path),
               "n_molecules": len(records), "run": result}
    csv_path = work_dir / "scoring.csv"
    if result.get("status") == "ok" and csv_path.is_file():
        payload["analysis"] = analyze_results(csv_path, top_n=args.top_n)
    emit_json(payload, args, summary=[
        f"Scored {len(records)} molecules with the "
        f"{args.scoring_profile or 'custom'} profile",
        f"  status: {result.get('status', 'dry-run')}",
    ])


def cmd_generate_hpc_script(args: argparse.Namespace) -> None:
    resources = SlurmResources(
        job_name=args.job_name, gpu=args.gpu, ngpu=args.ngpu,
        ncpus=args.ncpus, mem=args.mem, time=args.time,
        conda_env=args.conda_env, partition=args.partition,
        account=args.account, email=args.email,
    )
    script = generate_slurm_script("reinvent", resources, {
        "config": args.config,
        "log": args.log or "reinvent.log",
    })
    emit_text(script, args, executable=True, summary=[
        f"REINVENT SLURM script for {args.config}",
        f"  {args.ngpu}x{args.gpu}, {args.ncpus} CPU, {args.mem}, {args.time}",
        f"  submit with: sbatch {args.output}",
    ])


def build_parser() -> argparse.ArgumentParser:
    parser, sub = build_root_parser("ag-compound-synthesis", DESCRIPTION)

    p = sub.add_parser("check-setup",
                       help="Verify the REINVENT install, priors and device.")
    p.add_argument("--reinvent-dir", help="REINVENT 4 checkout.")
    add_output(p, help_text="Output JSON report.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_check_setup)

    p = sub.add_parser("prepare-seeds",
                       help="Validate and de-duplicate seed structures.")
    add_compound_inputs(p)
    p.add_argument("--generator", default="reinvent", choices=sorted(GENERATORS),
                   help="Generator the seeds are for; enables format checks.")
    p.add_argument("--include-names", action="store_true",
                   help="Write a second, name column.")
    p.add_argument("--report", help="Optional JSON report path.")
    add_output(p, help_text="Output .smi file.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_prepare_seeds)

    p = sub.add_parser("generate-config",
                       help="Generate a validated REINVENT TOML config.")
    p.add_argument("--mode", required=True, choices=sorted(RUN_MODES))
    p.add_argument("--generator", default="reinvent", choices=sorted(GENERATORS))
    p.add_argument("--prior", help="Prior filename, or a Mol2Mol variant "
                                   "such as scaffold_generic.")
    p.add_argument("--device", help="Torch device, e.g. cpu or cuda:0.")
    p.add_argument("--scoring-profile", choices=sorted(SCORING_PROFILES),
                   help="Named multi-objective profile.")
    p.add_argument("--component", action="append", metavar="SPEC",
                   help="Extra scoring component, e.g. "
                        "'SlogP:low=1:high=4:weight=2'. Repeatable.")
    p.add_argument("--aggregation", choices=sorted(AGGREGATION_TYPES))
    p.add_argument("--smiles-file", help="Seed, training or scoring SMILES.")
    p.add_argument("--validation-smiles-file",
                   help="Held-out set for transfer learning.")
    p.add_argument("--model-file", help="Model file for sampling.")
    p.add_argument("--num-steps", type=int, help="RL steps.")
    p.add_argument("--max-steps", type=int)
    p.add_argument("--min-steps", type=int)
    p.add_argument("--max-score", type=float,
                   help="Stop the stage once this total score is exceeded.")
    p.add_argument("--num-epochs", type=int, help="Transfer-learning epochs.")
    p.add_argument("--num-smiles", type=int, help="Molecules to sample.")
    p.add_argument("--batch-size", type=int)
    p.add_argument("--sigma", type=int, help="DAP sigma (default 128).")
    p.add_argument("--learning-rate", type=float)
    p.add_argument("--csv-prefix", help="Prefix for staged-learning CSVs.")
    p.add_argument("--output-csv", help="Output CSV for sampling or scoring.")
    p.add_argument("--output-model", help="Output model for transfer learning.")
    p.add_argument("--checkpoint-file")
    p.add_argument("--no-diversity-filter", action="store_true",
                   help="Disable the diversity filter. Not advised: it is "
                        "what stops the agent collapsing onto one scaffold.")
    p.add_argument("--inception-smiles-file",
                   help="Known-good SMILES to guide early RL (Reinvent only).")
    p.add_argument("--sample-strategy", choices=["multinomial", "beamsearch"])
    p.add_argument("--reinvent-dir")
    add_output(p, help_text="Output .toml config. Metadata is written "
                            "alongside it as .meta.json.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_generate_config)

    p = sub.add_parser("list-profiles",
                       help="Show scoring profiles, generators and run modes.")
    add_output(p, help_text="Output JSON.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_list_profiles)

    p = sub.add_parser("preflight",
                       help="Check a config and its inputs before running.")
    p.add_argument("--config", required=True)
    p.add_argument("--device")
    p.add_argument("--reinvent-dir")
    add_output(p, help_text="Output JSON report.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_preflight)

    p = sub.add_parser("run", help="Run REINVENT against a config.")
    p.add_argument("--config", required=True)
    p.add_argument("--log", help="REINVENT log file path.")
    p.add_argument("--device")
    p.add_argument("--reinvent-dir")
    p.add_argument("--timeout", type=float,
                   help="Wall-clock limit in seconds.")
    p.add_argument("--dry-run", action="store_true",
                   help="Preflight only; do not execute.")
    add_output(p, help_text="Output JSON run report.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_run)

    p = sub.add_parser("analyze-results",
                       help="Rank output molecules and diagnose the run.")
    p.add_argument("--csv", required=True, help="REINVENT output CSV.")
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument("--sort-by", help="Column to rank on.")
    p.add_argument("--min-score", type=float,
                   help="Discard molecules below this score.")
    p.add_argument("--no-properties", action="store_true",
                   help="Skip RDKit descriptor computation.")
    p.add_argument("--show", type=int, default=10)
    add_output(p, help_text="Output JSON analysis.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_analyze_results)

    p = sub.add_parser("score-molecules",
                       help="Score existing molecules with REINVENT scoring.")
    add_compound_inputs(p)
    p.add_argument("--scoring-profile", default="drug-like",
                   choices=sorted(SCORING_PROFILES))
    p.add_argument("--component", action="append", metavar="SPEC")
    p.add_argument("--aggregation", choices=sorted(AGGREGATION_TYPES))
    p.add_argument("--device", default="cpu")
    p.add_argument("--reinvent-dir")
    p.add_argument("--top-n", type=int, default=50)
    p.add_argument("--timeout", type=float, default=3600.0)
    p.add_argument("--dry-run", action="store_true")
    add_output(p, help_text="Output JSON results.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_score_molecules)

    p = sub.add_parser("generate-hpc-script",
                       help="Generate a SLURM script for a REINVENT run.")
    p.add_argument("--config", required=True)
    p.add_argument("--log")
    p.add_argument("--job-name", default="reinvent")
    p.add_argument("--gpu", default="a100")
    p.add_argument("--ngpu", type=int, default=1)
    p.add_argument("--ncpus", type=int, default=8)
    p.add_argument("--mem", default="64G")
    p.add_argument("--time", default="24:00:00")
    p.add_argument("--conda-env", default="reinvent4")
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
