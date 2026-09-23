"""Compound Synthesis CLI -- REINVENT4 Wrapper.

Unified CLI for generative molecular design with REINVENT 4. Provides
subcommands for setup verification, config generation, seed preparation,
execution, result analysis, and HPC script generation.

Usage:
    python compound_synthesis.py <subcommand> --output <file> [options]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import asdict
from pathlib import Path

# Locate reinvent_runner relative to this script.
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

import reinvent_runner as runner

# Default REINVENT4 location: <workspace_root>/REINVENT4
# Derived from script location, overridable via REINVENT_DIR env var.
_WORKSPACE_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(_SCRIPT_DIR))))
DEFAULT_REINVENT_DIR = os.environ.get(
    "REINVENT_DIR",
    os.path.join(_WORKSPACE_ROOT, "REINVENT4"),
)


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_json(data, output_file: str) -> None:
    """Write data to JSON file."""
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    with open(output_file, "w") as f:
        json.dump(data, f, indent=2, default=str)
    print(f"Success! Data written to: {output_file}")


def write_text(text: str, output_file: str) -> None:
    """Write text to file."""
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    with open(output_file, "w") as f:
        f.write(text)
    print(f"Success! Config written to: {output_file}")


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------

def cmd_check_setup(args: argparse.Namespace) -> None:
    """Check REINVENT4 installation and dependencies."""
    reinvent_dir = args.reinvent_dir or DEFAULT_REINVENT_DIR
    print(f"Checking REINVENT4 installation at: {reinvent_dir}")

    report = runner.check_reinvent_installation(reinvent_dir)
    data = asdict(report)

    # Print summary
    print(f"\n  Python:    {report.python_version}")
    print(f"  PyTorch:   {'[OK] ' + report.torch_version if report.torch_available else '[MISSING] Not installed'}")
    print(f"  RDKit:     {'[OK] ' + report.rdkit_version if report.rdkit_available else '[MISSING] Not installed'}")
    print(f"  REINVENT:  {'[OK] v' + report.reinvent_version if report.reinvent_importable else '[MISSING] Not importable'}")
    print(f"  CLI:       {report.reinvent_cli or '[MISSING] Not found'}")
    print(f"  GPU:       {'[OK] ' + report.gpu_name + f' ({report.gpu_memory_mb} MB)' if report.gpu_available else '[WARN] CPU only'}")
    print(f"  Priors:    {len(report.priors_found)} found, {len(report.priors_missing)} missing")

    if report.errors:
        print(f"\n  Errors:")
        for e in report.errors:
            print(f"    [ERROR] {e}")

    print(f"\n  {'[OK] READY' if report.ready else '[FAIL] NOT READY -- fix errors above'}")

    write_json(data, args.output)


def cmd_download_priors(args: argparse.Namespace) -> None:
    """Download prior models from Zenodo."""
    reinvent_dir = args.reinvent_dir or DEFAULT_REINVENT_DIR
    models = args.models if args.models else None
    print(f"Downloading prior models to: {reinvent_dir}/priors/")

    results = runner.download_priors(reinvent_dir, models)

    if results["downloaded"]:
        print(f"\n  Downloaded {len(results['downloaded'])} model(s)")
    if results["skipped"]:
        print(f"  Skipped {len(results['skipped'])} model(s) (already exist)")
    if results["errors"]:
        print(f"  Errors: {len(results['errors'])}")
        for e in results["errors"]:
            print(f"    [ERROR] {e}")

    write_json(results, args.output)


def cmd_generate_config(args: argparse.Namespace) -> None:
    """Generate REINVENT4 TOML configuration file."""
    reinvent_dir = args.reinvent_dir or DEFAULT_REINVENT_DIR
    mode = args.mode

    # Smart device detection if user didn't specify --device
    if args.device:
        device = args.device
        print(f"Generating REINVENT4 config: mode={mode}, device={device} (user-specified)")
    else:
        num_steps = args.num_steps or args.num_epochs or 100
        batch_size = args.batch_size or 128
        dev_info = runner.detect_device(
            run_mode=mode,
            num_steps=num_steps,
            batch_size=batch_size,
        )
        device = dev_info["device"]

        if dev_info["gpu_available"]:
            print(f"  GPU detected: {dev_info['gpu_name']} ({dev_info['gpu_memory_mb']} MB VRAM)")
        else:
            print(f"  No GPU detected, using CPU")

        print(f"  Auto-selected device: {device}")

        if dev_info["recommend_hpc"]:
            print(f"  [HPC RECOMMENDED] {dev_info['hpc_reason']}")
            # Auto-generate HPC script alongside the config
            hpc_output = os.path.splitext(args.output)[0] + "_hpc.sh"
            hpc_script = runner.generate_hpc_script(
                config_file=args.output,
                gpu="a100",
                ngpu=1,
                mem="64G",
                time="24:00:00",
            )
            os.makedirs(os.path.dirname(hpc_output) or ".", exist_ok=True)
            with open(hpc_output, "w") as f:
                f.write(hpc_script)
            print(f"  HPC script auto-generated: {hpc_output}")

    if mode == "staged_learning":
        # Resolve prior path
        prior = args.prior or ".reinvent"
        prior_path = runner.resolve_prior(prior, reinvent_dir)

        # Get scoring components
        if args.scoring_profile:
            components = runner.get_scoring_profile(args.scoring_profile)
        elif args.component:
            components = [runner.parse_custom_component(c) for c in args.component]
        else:
            components = runner.get_scoring_profile("drug-like")
            print("  Using default scoring profile: drug-like")

        # Handle seeds - write to smiles_file if needed
        smiles_file = args.smiles_file
        if args.seeds and not smiles_file:
            seeds_path = os.path.splitext(args.output)[0] + "_seeds.smi"
            with open(seeds_path, "w") as f:
                for s in args.seeds:
                    f.write(s.strip() + "\n")
            smiles_file = seeds_path
            print(f"  Seeds written to: {seeds_path}")

        toml_str = runner.generate_staged_learning_toml(
            prior_file=prior_path,
            device=device,
            batch_size=args.batch_size or 128,
            num_steps=args.num_steps or 300,
            max_score=args.max_score or 0.9,
            min_steps=args.min_steps or 50,
            sigma=args.sigma or 128,
            scoring_components=components,
            smiles_file=smiles_file,
            summary_csv_prefix=args.csv_prefix or "staged_learning",
        )
        write_text(toml_str, args.output)

    elif mode == "transfer_learning":
        prior = args.prior or ".reinvent"
        prior_path = runner.resolve_prior(prior, reinvent_dir)

        if not args.smiles_file:
            print("ERROR: --smiles-file is required for transfer_learning mode")
            sys.exit(1)

        output_model = args.output_model or os.path.splitext(args.output)[0] + ".model"

        toml_str = runner.generate_transfer_learning_toml(
            input_model_file=prior_path,
            output_model_file=output_model,
            smiles_file=args.smiles_file,
            device=device,
            num_epochs=args.num_epochs or 50,
            batch_size=args.batch_size or 64,
            validation_smiles_file=args.validation_smiles_file,
        )
        write_text(toml_str, args.output)

    elif mode == "sampling":
        model_file = args.model_file
        if not model_file:
            prior = args.prior or ".reinvent"
            model_file = runner.resolve_prior(prior, reinvent_dir)

        toml_str = runner.generate_sampling_toml(
            model_file=model_file,
            device=device,
            num_smiles=args.num_smiles or 1000,
            output_file=args.sample_output or "samples.csv",
            smiles_file=args.smiles_file,
        )
        write_text(toml_str, args.output)

    elif mode == "scoring":
        if not args.smiles_file:
            print("ERROR: --smiles-file is required for scoring mode")
            sys.exit(1)

        if args.scoring_profile:
            components = runner.get_scoring_profile(args.scoring_profile)
        elif args.component:
            components = [runner.parse_custom_component(c) for c in args.component]
        else:
            components = runner.get_scoring_profile("drug-like")

        toml_str = runner.generate_scoring_toml(
            smiles_file=args.smiles_file,
            scoring_components=components,
            output_csv=args.scoring_csv or "scoring_results.csv",
        )
        write_text(toml_str, args.output)

    else:
        print(f"ERROR: Unsupported mode '{mode}'")
        sys.exit(1)


def cmd_prepare_seeds(args: argparse.Namespace) -> None:
    """Validate and prepare seed SMILES."""
    smiles = args.smiles or []

    # Also read from file if provided
    if args.smiles_file:
        with open(args.smiles_file) as f:
            for line in f:
                s = line.strip().split()[0] if line.strip() else ""
                if s:
                    smiles.append(s)

    if not smiles:
        print("ERROR: No SMILES provided. Use --smiles or --smiles-file.")
        sys.exit(1)

    print(f"Preparing {len(smiles)} seed SMILES...")

    report = runner.prepare_seeds(smiles, args.output)
    data = asdict(report)

    print(f"  Valid: {report.valid}")
    print(f"  Invalid: {report.invalid}")
    print(f"  Duplicates removed: {report.duplicates_removed}")
    if report.invalid_smiles:
        print(f"  Invalid SMILES: {report.invalid_smiles[:5]}")

    if args.report:
        write_json(data, args.report)
    else:
        print(f"Success! Seeds written to: {args.output}")


def cmd_run(args: argparse.Namespace) -> None:
    """Execute REINVENT4 locally."""
    reinvent_dir = args.reinvent_dir or DEFAULT_REINVENT_DIR
    config = args.config

    if not os.path.exists(config):
        print(f"ERROR: Config file not found: {config}")
        sys.exit(1)

    print(f"Running REINVENT4...")
    print(f"  Config: {config}")
    print(f"  REINVENT dir: {reinvent_dir}")
    print(f"  Log: {args.log or 'stderr'}")

    result = runner.run_reinvent(
        config_path=config,
        reinvent_dir=reinvent_dir,
        log_file=args.log,
        device=args.device,
    )

    if result["success"]:
        print(f"\n  [OK] REINVENT4 completed successfully")
    else:
        print(f"\n  [FAIL] REINVENT4 failed (return code: {result['return_code']})")
        if result["stderr_tail"]:
            print("  Last error lines:")
            for line in result["stderr_tail"][-5:]:
                print(f"    {line}")

    write_json(result, args.output)


def cmd_analyze_results(args: argparse.Namespace) -> None:
    """Analyze REINVENT4 CSV output."""
    csv_path = args.csv
    if not os.path.exists(csv_path):
        print(f"ERROR: CSV file not found: {csv_path}")
        sys.exit(1)

    print(f"Analyzing: {csv_path}")

    results = runner.analyze_results(
        csv_path=csv_path,
        top_n=args.top_n or 20,
        sort_by=args.sort_by or "total_score",
    )

    if "error" in results:
        print(f"  ERROR: {results['error']}")
    else:
        stats = results.get("stats", {})
        print(f"  Total compounds: {stats.get('total_compounds', 0)}")
        print(f"  Unique molecules: {stats.get('unique_molecules', 'N/A')}")
        print(f"  Mean score: {stats.get('mean_score', 0):.4f}")
        print(f"  Max score: {stats.get('max_score', 0):.4f}")
        print(f"  Top {len(results.get('top_compounds', []))} compounds extracted")

    write_json(results, args.output)


def cmd_generate_hpc_script(args: argparse.Namespace) -> None:
    """Generate SLURM HPC submission script."""
    config = args.config
    print(f"Generating SLURM script for: {config}")

    script = runner.generate_hpc_script(
        config_file=config,
        gpu=args.gpu or "a100",
        ngpu=args.ngpu or 1,
        mem=args.mem or "64G",
        time=args.time or "24:00:00",
        conda_env=args.conda_env or "reinvent4",
        job_name=args.job_name or "reinvent4",
    )

    write_text(script, args.output)


def cmd_score_molecules(args: argparse.Namespace) -> None:
    """Score existing molecules using REINVENT4 scoring framework."""
    reinvent_dir = args.reinvent_dir or DEFAULT_REINVENT_DIR

    if not args.smiles_file and not args.smiles:
        print("ERROR: Provide --smiles-file or --smiles")
        sys.exit(1)

    # If inline SMILES, write to temp file
    smiles_file = args.smiles_file
    if not smiles_file and args.smiles:
        smiles_file = os.path.splitext(args.output)[0] + "_input.smi"
        with open(smiles_file, "w") as f:
            for s in args.smiles:
                f.write(s.strip() + "\n")

    # Get scoring components
    if args.scoring_profile:
        components = runner.get_scoring_profile(args.scoring_profile)
    else:
        components = runner.get_scoring_profile("drug-like")

    # Generate scoring config
    scoring_csv = os.path.splitext(args.output)[0] + "_scores.csv"
    config_path = os.path.splitext(args.output)[0] + "_scoring.toml"

    toml_str = runner.generate_scoring_toml(
        smiles_file=smiles_file,
        scoring_components=components,
        output_csv=scoring_csv,
    )

    with open(config_path, "w") as f:
        f.write(toml_str)

    print(f"  Config: {config_path}")
    print(f"  Running REINVENT scoring...")

    result = runner.run_reinvent(
        config_path=config_path,
        reinvent_dir=reinvent_dir,
        device=args.device or "cpu",
    )

    if result["success"] and os.path.exists(scoring_csv):
        analysis = runner.analyze_results(scoring_csv, top_n=args.top_n or 50)
        write_json(analysis, args.output)
    else:
        result["config_file"] = config_path
        write_json(result, args.output)


# ---------------------------------------------------------------------------
# Main CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        prog="compound_synthesis",
        description="Compound Synthesis — REINVENT4 generative chemistry wrapper",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- check-setup ---
    p_setup = subparsers.add_parser("check-setup", help="Verify REINVENT4 installation")
    p_setup.add_argument("--reinvent-dir", default=None, help="Path to REINVENT4 repo")
    p_setup.add_argument("--output", required=True, help="Output JSON report")
    p_setup.set_defaults(func=cmd_check_setup)

    # --- download-priors ---
    p_dl = subparsers.add_parser("download-priors", help="Download prior models from Zenodo")
    p_dl.add_argument("--reinvent-dir", default=None, help="Path to REINVENT4 repo")
    p_dl.add_argument("--models", nargs="+", help="Model names to download (e.g. reinvent libinvent)")
    p_dl.add_argument("--output", required=True, help="Output JSON report")
    p_dl.set_defaults(func=cmd_download_priors)

    # --- generate-config ---
    p_cfg = subparsers.add_parser("generate-config", help="Generate REINVENT4 TOML config")
    p_cfg.add_argument("--mode", required=True,
                        choices=["staged_learning", "transfer_learning", "sampling", "scoring"],
                        help="REINVENT4 run mode")
    p_cfg.add_argument("--generator", default="reinvent",
                        choices=["reinvent", "libinvent", "linkinvent", "mol2mol", "pepinvent"],
                        help="Generator architecture")
    p_cfg.add_argument("--prior", default=None,
                        help="Prior model key (e.g. .reinvent) or file path")
    p_cfg.add_argument("--model-file", default=None, help="Model file for sampling")
    p_cfg.add_argument("--scoring-profile", default=None,
                        choices=["drug-like", "lead-like", "fragment-like",
                                 "kinase-inhibitor", "anti-tb"],
                        help="Built-in scoring profile")
    p_cfg.add_argument("--component", action="append",
                        help="Custom component spec: 'Name:weight=1.0,transform=sigmoid,...'")
    p_cfg.add_argument("--seeds", nargs="+", help="Seed SMILES for inception/TL")
    p_cfg.add_argument("--smiles-file", help="SMILES file (for TL, scoring, or seeds)")
    p_cfg.add_argument("--validation-smiles-file", help="Validation SMILES (TL only)")
    p_cfg.add_argument("--num-steps", type=int, help="RL steps (staged_learning)")
    p_cfg.add_argument("--num-epochs", type=int, help="Training epochs (transfer_learning)")
    p_cfg.add_argument("--num-smiles", type=int, help="SMILES to sample (sampling)")
    p_cfg.add_argument("--batch-size", type=int, help="Batch size")
    p_cfg.add_argument("--max-score", type=float, help="Early termination score")
    p_cfg.add_argument("--min-steps", type=int, help="Min steps before early stop")
    p_cfg.add_argument("--sigma", type=int, help="DAP sigma parameter")
    p_cfg.add_argument("--device", default=None, help="Torch device (cpu/cuda:0)")
    p_cfg.add_argument("--sample-output", default=None, help="Sampling output CSV filename")
    p_cfg.add_argument("--scoring-csv", default=None, help="Scoring output CSV filename")
    p_cfg.add_argument("--output-model", default=None, help="Output model file (TL)")
    p_cfg.add_argument("--csv-prefix", default=None, help="CSV prefix for staged learning")
    p_cfg.add_argument("--reinvent-dir", default=None, help="Path to REINVENT4 repo")
    p_cfg.add_argument("--output", required=True, help="Output TOML config file")
    p_cfg.set_defaults(func=cmd_generate_config)

    # --- prepare-seeds ---
    p_seed = subparsers.add_parser("prepare-seeds", help="Validate and prepare seed SMILES")
    p_seed.add_argument("--smiles", nargs="+", help="SMILES strings")
    p_seed.add_argument("--smiles-file", help="File with SMILES (one per line)")
    p_seed.add_argument("--report", help="Optional JSON report file")
    p_seed.add_argument("--output", required=True, help="Output .smi file")
    p_seed.set_defaults(func=cmd_prepare_seeds)

    # --- run ---
    p_run = subparsers.add_parser("run", help="Execute REINVENT4 locally")
    p_run.add_argument("--config", required=True, help="TOML config file")
    p_run.add_argument("--reinvent-dir", default=None, help="Path to REINVENT4 repo")
    p_run.add_argument("--log", default=None, help="Log file path")
    p_run.add_argument("--device", default=None, help="Override device (cpu/cuda:0)")
    p_run.add_argument("--output", required=True, help="Output JSON report")
    p_run.set_defaults(func=cmd_run)

    # --- analyze-results ---
    p_analyze = subparsers.add_parser("analyze-results", help="Parse REINVENT4 CSV output")
    p_analyze.add_argument("--csv", required=True, help="REINVENT4 output CSV file")
    p_analyze.add_argument("--top-n", type=int, default=20, help="Top N compounds")
    p_analyze.add_argument("--sort-by", default="total_score", help="Sort column")
    p_analyze.add_argument("--output", required=True, help="Output JSON analysis")
    p_analyze.set_defaults(func=cmd_analyze_results)

    # --- generate-hpc-script ---
    p_hpc = subparsers.add_parser("generate-hpc-script", help="Generate SLURM script")
    p_hpc.add_argument("--config", required=True, help="REINVENT4 config file")
    p_hpc.add_argument("--gpu", default="a100", help="GPU type")
    p_hpc.add_argument("--ngpu", type=int, default=1, help="Number of GPUs")
    p_hpc.add_argument("--mem", default="64G", help="Memory")
    p_hpc.add_argument("--time", default="24:00:00", help="Walltime")
    p_hpc.add_argument("--conda-env", default="reinvent4", help="Conda environment")
    p_hpc.add_argument("--job-name", default="reinvent4", help="SLURM job name")
    p_hpc.add_argument("--output", required=True, help="Output .sh script")
    p_hpc.set_defaults(func=cmd_generate_hpc_script)

    # --- score-molecules ---
    p_score = subparsers.add_parser("score-molecules", help="Score SMILES using REINVENT4")
    p_score.add_argument("--smiles-file", help="SMILES file")
    p_score.add_argument("--smiles", nargs="+", help="Inline SMILES")
    p_score.add_argument("--scoring-profile", default="drug-like",
                          choices=["drug-like", "lead-like", "fragment-like",
                                   "kinase-inhibitor", "anti-tb"])
    p_score.add_argument("--top-n", type=int, default=50)
    p_score.add_argument("--device", default="cpu")
    p_score.add_argument("--reinvent-dir", default=None)
    p_score.add_argument("--output", required=True, help="Output JSON results")
    p_score.set_defaults(func=cmd_score_molecules)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
