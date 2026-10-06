"""CLI for the drug-discovery-wizard orchestrator.

This is a true orchestrator. The implementation it replaces was a
1,533-line monolith that *re-implemented* every subcommand of the other
four skills instead of calling them, and the copies had already diverged
from their originals: the docking argument bug (``--size`` defined,
``args.box_size`` read) existed in the screening skill but not in the
wizard's own copy, so the two behaved differently.

Nothing here duplicates science. Each stage calls the same library
functions the individual skills call, so a fix reaches every caller at once.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Any

from ..chem.admet import TIERS, screen_compounds, tier_summary
from ..chem.smiles import CompoundRecord
from ..errors import AgSkillsError, UsageError
from ..io_utils import write_json, write_text
from ._common import (
    add_common_flags,
    add_output,
    build_root_parser,
    emit_json,
    run_cli,
)

_TIERS = tier_summary()

DESCRIPTION = f"""\
Orchestrate the full structure-based discovery pipeline.

Stages:
  1. target      fetch, clean and assess the receptor; define the site
  2. sourcing    retrieve known actives from ChEMBL and natural products
  3. triage      relaxed ADMET screening ({_TIERS['relaxed']} checks)
  4. docking     AutoDock Vina against the prepared receptor
  5. generative  REINVENT 4 configuration and submission script
  6. selection   strict ADMET screening ({_TIERS['strict']} checks)
  7. dynamics    GROMACS setup and submission script

Use 'plan' to see what a pipeline would do without running it, and
'run' to execute the stages that can run locally. Stages needing a GPU
cluster produce submission scripts rather than pretending to run.

For a single step, call that skill directly: ag-target-preparation,
ag-compound-screening, ag-compound-synthesis or ag-md-simulation.
"""

STAGES: dict[str, dict[str, Any]] = {
    "target": {
        "skill": "target-preparation",
        "local": True,
        "needs": "a PDB id, a UniProt accession, or a local PDB file",
        "produces": "cleaned PDB, receptor PDBQT, binding-site JSON, "
                    "quality report",
    },
    "sourcing": {
        "skill": "compound-screening",
        "local": True,
        "needs": "a ChEMBL target id or a search term",
        "produces": "known actives with potency, natural-product hits",
    },
    "triage": {
        "skill": "compound-screening",
        "local": True,
        "needs": "compounds from sourcing",
        "produces": f"relaxed ADMET report ({_TIERS['relaxed']} checks) and "
                    "a seed set",
    },
    "docking": {
        "skill": "compound-screening",
        "local": False,
        "needs": "receptor PDBQT, binding site, AutoDock Vina",
        "produces": "ranked binding affinities and pose files",
    },
    "generative": {
        "skill": "compound-synthesis",
        "local": False,
        "needs": "seed structures, REINVENT 4, a GPU",
        "produces": "REINVENT TOML and a SLURM script",
    },
    "selection": {
        "skill": "compound-screening",
        "local": True,
        "needs": "generated or docked compounds",
        "produces": f"strict ADMET report ({_TIERS['strict']} checks)",
    },
    "dynamics": {
        "skill": "md-simulation",
        "local": False,
        "needs": "a cleaned PDB, GROMACS, a GPU",
        "produces": "MDP files, workflow and submission scripts",
    },
}


def _stage_dir(root: Path, name: str) -> Path:
    path = root / name
    path.mkdir(parents=True, exist_ok=True)
    return path


def cmd_plan(args: argparse.Namespace) -> None:
    """Describe the pipeline for a given target without executing it."""
    from ..docking.vina import find_vina
    from ..reinvent.priors import check_setup
    from ..targets.pdbqt import available_converters

    vina_kind, vina_detail = find_vina()
    reinvent = check_setup(args.reinvent_dir)
    converters = available_converters()

    # Auto-install GROMACS if absent
    try:
        import shutil
        from .._deps import ensure_gromacs
        gromacs = shutil.which("gmx") or shutil.which("gmx_mpi")
        if not gromacs:
            ensure_gromacs()
            gromacs = shutil.which("gmx") or shutil.which("gmx_mpi")
    except Exception:
        gromacs = None

    # Auto-install REINVENT if absent
    if not reinvent["ready"]:
        try:
            from .._deps import ensure_reinvent
            if ensure_reinvent():
                reinvent = check_setup(args.reinvent_dir)
        except Exception:
            pass

    conv_note = (", ".join(converters)
                 if converters
                 else "none (a cleaned PDB is still produced)")
    readiness = {
        "target": {"ready": bool(converters) or True,
                   "note": f"PDBQT converters available: {conv_note}"},
        "sourcing": {"ready": True, "note": "Needs internet access."},
        "triage": {"ready": True, "note": "Offline RDKit battery."},
        "docking": {"ready": vina_kind != "none",
                    "note": f"Vina: {vina_kind} {vina_detail}"
                            if vina_kind != "none"
                            else "AutoDock Vina is not installed; this stage "
                                 "will produce a cluster script instead."},
        "generative": {"ready": reinvent["ready"],
                       "note": reinvent["summary"]},
        "selection": {"ready": True, "note": "Offline RDKit battery."},
        "dynamics": {"ready": bool(gromacs),
                     "note": f"GROMACS: {gromacs}" if gromacs
                             else "GROMACS is not installed; setup files and "
                                  "a cluster script are still produced."},
    }

    selected = args.stages or list(STAGES)
    unknown = [s for s in selected if s not in STAGES]
    if unknown:
        raise UsageError(
            f"Unknown stage(s): {', '.join(unknown)}",
            hint=f"Available: {', '.join(STAGES)}",
        )

    plan = []
    for index, name in enumerate(selected, 1):
        stage = STAGES[name]
        plan.append({
            "order": index,
            "stage": name,
            "skill": stage["skill"],
            "runs_locally": stage["local"],
            "needs": stage["needs"],
            "produces": stage["produces"],
            "ready": readiness[name]["ready"],
            "readiness_note": readiness[name]["note"],
        })

    payload = {
        "target": {
            "pdb_id": args.pdb_id, "uniprot_id": args.uniprot_id,
            "input_file": args.input_file, "chembl_target": args.chembl_target,
        },
        "workdir": str(args.workdir),
        "stages": plan,
        "environment": {
            "vina": f"{vina_kind}: {vina_detail}",
            "pdbqt_converters": converters,
            "gromacs": gromacs,
            "reinvent_ready": reinvent["ready"],
            "reinvent_problems": reinvent["problems"],
        },
        "two_tier_rationale": (
            "Relaxed screening is applied before generative chemistry so a "
            "diverse seed set survives; strict screening is applied "
            "afterwards, when the aim is a short, defensible candidate list. "
            "Applying strict thresholds to seeds discards chemistry the "
            "generator could have repaired."
        ),
        "note": "This skill (drug-discovery-wizard) produced this plan.",
    }
    lines = ["Pipeline plan:"]
    for entry in plan:
        mark = "ready" if entry["ready"] else "NOT READY"
        where = "local" if entry["runs_locally"] else "cluster"
        lines.append(f"  {entry['order']}. {entry['stage']:<11s} "
                     f"[{where:<7s}] [{mark}]  {entry['skill']}")
        lines.append(f"       produces: {entry['produces']}")
        if not entry["ready"]:
            lines.append(f"       {entry['readiness_note']}")
    emit_json(payload, args, summary=lines)


def cmd_run(args: argparse.Namespace) -> None:
    """Execute the pipeline stages that can run here."""
    root = Path(args.workdir)
    root.mkdir(parents=True, exist_ok=True)
    selected = args.stages or ["target", "sourcing", "triage"]
    unknown = [s for s in selected if s not in STAGES]
    if unknown:
        raise UsageError(f"Unknown stage(s): {', '.join(unknown)}",
                         hint=f"Available: {', '.join(STAGES)}")

    results: dict[str, Any] = {}
    artefacts: dict[str, str] = {}
    compounds: list[CompoundRecord] = []
    started = time.monotonic()
    lines = [f"Pipeline in {root}"]

    def record(stage: str, status: str, detail: Any) -> None:
        results[stage] = {"status": status, **(detail if isinstance(detail, dict)
                                               else {"detail": detail})}

    # ---- 1. target ---------------------------------------------------
    if "target" in selected:
        stage_dir = _stage_dir(root, "01_target")
        try:
            from ..targets.clean import clean_structure
            from ..targets.fetch import fetch_alphafold, fetch_pdb
            from ..targets.pdbqt import receptor_to_pdbqt
            from ..targets.quality import assess_structure
            from ..targets.site import list_ligands, site_from_hetatm

            if args.pdb_id:
                download = fetch_pdb(args.pdb_id, stage_dir)
                structure = download.path
                source = download.as_dict()
            elif args.uniprot_id:
                download = fetch_alphafold(args.uniprot_id, stage_dir)
                structure = download.path
                source = download.as_dict()
            elif args.input_file:
                structure = Path(args.input_file)
                source = {"path": str(structure), "source": "local file"}
            else:
                raise UsageError(
                    "The target stage needs a structure.",
                    hint="Pass --pdb-id, --uniprot-id or --input-file.",
                )

            cleaned = stage_dir / f"{structure.stem}_clean.pdb"
            clean = clean_structure(structure, cleaned, chain=args.chain)
            quality = assess_structure(structure)
            conversion = receptor_to_pdbqt(cleaned)
            ligands = list_ligands(structure)

            site_payload = None
            if ligands:
                site = site_from_hetatm(structure)
                site_payload = site.as_dict()
                site_file = stage_dir / "binding_site.json"
                write_json(site_payload, site_file)
                artefacts["binding_site"] = str(site_file)

            artefacts["cleaned_pdb"] = str(cleaned)
            if conversion.ok and conversion.output_path:
                artefacts["receptor_pdbqt"] = conversion.output_path

            record("target", "ok", {
                "source": source, "cleaning": clean.as_dict(),
                "quality": quality.as_dict(), "pdbqt": conversion.as_dict(),
                "candidate_ligands": ligands, "binding_site": site_payload,
            })
            lines += [
                f"  target: {structure.name} -> {clean.atoms_kept} atoms, "
                f"{len(clean.heterogens_removed)} heterogen type(s) removed",
                f"          verdict {quality.verdict}; "
                f"{len(ligands)} candidate ligand(s)",
            ]
            if site_payload:
                lines.append(
                    f"          site centre "
                    f"({site_payload['center_x']}, {site_payload['center_y']}, "
                    f"{site_payload['center_z']})"
                )
        except AgSkillsError as exc:
            record("target", "failed", {"error": str(exc)})
            lines.append(f"  target: FAILED - {exc.message}")
            if not args.keep_going:
                emit_json({"workdir": str(root), "stages": results,
                           "artefacts": artefacts}, args, summary=lines)
                raise

    # ---- 2. sourcing -------------------------------------------------
    if "sourcing" in selected:
        stage_dir = _stage_dir(root, "02_sourcing")
        sourced: list[CompoundRecord] = []
        detail: dict[str, Any] = {}
        try:
            if args.chembl_target:
                from ..libraries.chembl import query_activities
                actives = query_activities(args.chembl_target,
                                           pchembl_min=args.pchembl_min,
                                           limit=args.limit)
                write_json(actives, stage_dir / "chembl_actives.json")
                detail["chembl"] = {
                    "distinct_molecules": actives["distinct_molecules_found"],
                    "returned": actives["returned"],
                }
                sourced += [
                    CompoundRecord(name=c["molecule_chembl_id"],
                                   smiles=c["canonical_smiles"],
                                   source="ChEMBL")
                    for c in actives["compounds"]
                ]
                lines.append(
                    f"  sourcing: ChEMBL {args.chembl_target} -> "
                    f"{actives['returned']} actives"
                )
            if args.natural_product_query:
                from ..libraries.coconut import search_coconut
                natural = search_coconut(args.natural_product_query,
                                         limit=args.limit)
                write_json(natural, stage_dir / "coconut.json")
                detail["coconut"] = {"returned": natural["returned"],
                                      "with_smiles": natural["with_smiles"]}
                sourced += [
                    CompoundRecord(name=c["id"] or c["name"],
                                   smiles=c["smiles"], source="COCONUT")
                    for c in natural["compounds"] if c["smiles"]
                ]
                lines.append(
                    f"  sourcing: COCONUT {args.natural_product_query!r} -> "
                    f"{natural['with_smiles']} structures"
                )
            if not args.chembl_target and not args.natural_product_query:
                record("sourcing", "skipped",
                       {"reason": "neither --chembl-target nor "
                                  "--natural-product-query was given"})
                lines.append("  sourcing: skipped (no source specified)")
            else:
                compounds += sourced
                detail["total_sourced"] = len(sourced)
                record("sourcing", "ok", detail)
        except AgSkillsError as exc:
            record("sourcing", "failed", {"error": str(exc)})
            lines.append(f"  sourcing: FAILED - {exc.message}")
            if not args.keep_going:
                emit_json({"workdir": str(root), "stages": results,
                           "artefacts": artefacts}, args, summary=lines)
                raise

    # ---- 3. triage ---------------------------------------------------
    if "triage" in selected:
        stage_dir = _stage_dir(root, "03_triage")
        if not compounds:
            record("triage", "skipped",
                   {"reason": "no compounds were sourced"})
            lines.append("  triage: skipped (no compounds)")
        else:
            report = screen_compounds(compounds, strictness="relaxed")
            path = write_json(report, stage_dir / "relaxed_admet.json")
            artefacts["triage_report"] = str(path)
            seeds = [c for c in report["compounds"] if c["overall_pass"]]
            if seeds:
                seed_file = stage_dir / "seeds.smi"
                write_text("\n".join(c["smiles"] for c in seeds) + "\n",
                            seed_file)
                artefacts["seeds"] = str(seed_file)
            record("triage", "ok", {
                "tier": report["tier"],
                "checks_per_compound": report["checks_per_compound"],
                "total": report["total_compounds"], "passed": report["passed"],
                "failed": report["failed"], "errored": report["errored"],
            })
            lines.append(
                f"  triage: {report['tier']} "
                f"({report['checks_per_compound']} checks) -> "
                f"{report['passed']}/{report['total_compounds']} passed"
            )

    # ---- 4. docking --------------------------------------------------
    if "docking" in selected:
        stage_dir = _stage_dir(root, "04_docking")
        receptor = args.receptor or artefacts.get("receptor_pdbqt")
        site_file = args.site or artefacts.get("binding_site")
        if not receptor or not site_file:
            record("docking", "skipped", {
                "reason": "needs a receptor PDBQT and a binding site",
                "receptor": receptor, "site": site_file,
            })
            lines.append("  docking: skipped (no receptor PDBQT or no "
                          "binding site)")
        elif not compounds:
            record("docking", "skipped", {"reason": "no compounds"})
            lines.append("  docking: skipped (no compounds)")
        else:
            from ..docking.vina import dock_compounds, find_vina
            from ..io_utils import read_json
            from ..targets.site import BindingSite
            kind, _ = find_vina()
            if kind == "none":
                from ..hpc import SlurmResources, generate_slurm_script
                script = generate_slurm_script(
                    "vina-screen",
                    SlurmResources(job_name="vina", ngpu=0, gpu=None,
                                   ncpus=16, mem="32G", time="12:00:00",
                                   array="1-50"),
                    {"receptor": receptor, "ligand_dir": "ligands",
                     "site_config": site_file},
                )
                path = write_text(script, stage_dir / "submit_vina.sh",
                                   executable=True)
                artefacts["vina_script"] = str(path)
                record("docking", "deferred", {
                    "reason": "AutoDock Vina is not installed here",
                    "script": str(path),
                })
                lines.append(f"  docking: deferred to cluster -> {path}")
            else:
                data = read_json(site_file)
                site = BindingSite(
                    center_x=data["center_x"], center_y=data["center_y"],
                    center_z=data["center_z"], size_x=data["size_x"],
                    size_y=data["size_y"], size_z=data["size_z"],
                    method=data.get("method", "from file"),
                )
                subset = compounds[:args.dock_limit]
                result = dock_compounds(subset, receptor, site,
                                        output_dir=stage_dir / "poses",
                                        exhaustiveness=args.exhaustiveness)
                write_json(result, stage_dir / "docking.json")
                record("docking", "ok", {
                    "docked": result["docked"], "failed": result["failed"],
                    "best_affinity_kcal_mol": result["best_affinity_kcal_mol"],
                })
                lines.append(
                    f"  docking: {result['docked']}/{len(subset)} docked, "
                    f"best {result['best_affinity_kcal_mol']} kcal/mol"
                )

    # ---- 5. generative ----------------------------------------------
    if "generative" in selected:
        stage_dir = _stage_dir(root, "05_generative")
        from ..hpc import SlurmResources, generate_slurm_script
        from ..reinvent.config import build_config
        seed_file = artefacts.get("seeds")
        config = build_config(
            mode="staged_learning", generator="reinvent",
            scoring_profile=args.scoring_profile,
            inception_smiles_file=seed_file,
            num_steps=args.num_steps, batch_size=args.batch_size,
            device=args.device, reinvent_dir=args.reinvent_dir,
            csv_prefix="generated",
        )
        config_path = write_text(config.toml, stage_dir / "reinvent.toml")
        write_json(config.as_dict(), stage_dir / "reinvent.meta.json")
        script = generate_slurm_script(
            "reinvent",
            SlurmResources(job_name="reinvent", gpu="a100", ngpu=1,
                           ncpus=8, mem="64G", time="24:00:00",
                           conda_env="reinvent4"),
            {"config": config_path.name, "log": "reinvent.log"},
        )
        script_path = write_text(script, stage_dir / "submit_reinvent.sh",
                                  executable=True)
        artefacts["reinvent_config"] = str(config_path)
        artefacts["reinvent_script"] = str(script_path)
        record("generative", "prepared", {
            "config": str(config_path), "script": str(script_path),
            "scoring_profile": args.scoring_profile,
            "warnings": config.warnings,
        })
        lines.append(f"  generative: config and submission script prepared "
                      f"-> {script_path}")

    # ---- 6. selection -----------------------------------------------
    if "selection" in selected:
        stage_dir = _stage_dir(root, "06_selection")
        if not compounds:
            record("selection", "skipped", {"reason": "no compounds"})
            lines.append("  selection: skipped (no compounds)")
        else:
            report = screen_compounds(compounds, strictness="strict")
            path = write_json(report, stage_dir / "strict_admet.json")
            artefacts["selection_report"] = str(path)
            record("selection", "ok", {
                "tier": report["tier"],
                "checks_per_compound": report["checks_per_compound"],
                "total": report["total_compounds"], "passed": report["passed"],
                "failed": report["failed"], "errored": report["errored"],
            })
            lines.append(
                f"  selection: {report['tier']} "
                f"({report['checks_per_compound']} checks) -> "
                f"{report['passed']}/{report['total_compounds']} passed"
            )

    # ---- 7. dynamics -------------------------------------------------
    if "dynamics" in selected:
        stage_dir = _stage_dir(root, "07_dynamics")
        cleaned = args.input_file or artefacts.get("cleaned_pdb")
        if not cleaned:
            record("dynamics", "skipped",
                   {"reason": "no cleaned PDB from the target stage"})
            lines.append("  dynamics: skipped (no structure)")
        else:
            from ..hpc import SlurmResources, generate_slurm_script
            from ..md.workflow import (
                generate_gmxapi_workflow,
                generate_setup_shell,
                workflow_bundle,
            )
            bundle = workflow_bundle(
                pdb=Path(cleaned).name, force_field=args.force_field,
                water=args.water, temperature=args.temperature,
                production_ns=args.production_ns,
            )
            for filename, content in bundle["mdp_files"].items():
                write_text(content, stage_dir / filename)
            write_text(generate_gmxapi_workflow(bundle),
                        stage_dir / "run_md.py")
            write_text(generate_setup_shell(bundle), stage_dir / "setup_md.sh",
                        executable=True)
            script = generate_slurm_script(
                "gromacs-md",
                SlurmResources(job_name="md", gpu="a100", ngpu=1, ncpus=8,
                               mem="64G", time="48:00:00", conda_env="gromacs"),
                {"production_ns": args.production_ns},
            )
            script_path = write_text(script, stage_dir / "submit_md.sh",
                                      executable=True)
            artefacts["md_script"] = str(script_path)
            record("dynamics", "prepared", {
                "force_field": bundle["settings"]["force_field"],
                "temperature_K": bundle["settings"]["temperature_K"],
                "production_ns": bundle["settings"]["production_ns"],
                "script": str(script_path),
                "warnings": bundle["warnings"],
            })
            lines.append(
                f"  dynamics: {bundle['settings']['force_field']} at "
                f"{bundle['settings']['temperature_K']:g} K, "
                f"{bundle['settings']['production_ns']:g} ns -> {script_path}"
            )

    elapsed = time.monotonic() - started
    payload = {
        "skill": "drug-discovery-wizard",
        "workdir": str(root),
        "stages_requested": selected,
        "stages": results,
        "artefacts": artefacts,
        "compounds_in_flight": len(compounds),
        "elapsed_seconds": round(elapsed, 1),
        "note": ("This pipeline was orchestrated by the "
                 "drug-discovery-wizard skill. Stages marked 'deferred' or "
                 "'prepared' produced cluster scripts rather than running "
                 "GPU work locally."),
    }
    lines.append(f"  completed in {elapsed:.1f}s")
    emit_json(payload, args, summary=lines)


def cmd_describe(args: argparse.Namespace) -> None:
    payload = {
        "skill": "drug-discovery-wizard",
        "stages": STAGES,
        "tiers": {name: tier.describe() for name, tier in TIERS.items()},
        "sub_skills": {
            "target-preparation": "ag-target-preparation",
            "compound-screening": "ag-compound-screening",
            "compound-synthesis": "ag-compound-synthesis",
            "md-simulation": "ag-md-simulation",
        },
        "do_not_use_when": [
            "Only one step is needed - call that skill directly.",
            "The task is wet-lab protocol design.",
            "The task is clinical trial design.",
            "The task is retrosynthetic route planning - use AiZynthFinder, "
            "ASKCOS or IBM RXN.",
        ],
    }
    lines = ["drug-discovery-wizard stages:"]
    for name, stage in STAGES.items():
        lines.append(f"  {name:<11s} {stage['skill']:<20s} "
                     f"{'local' if stage['local'] else 'cluster'}")
        lines.append(f"              needs: {stage['needs']}")
        lines.append(f"              produces: {stage['produces']}")
    emit_json(payload, args, summary=lines)


def build_parser() -> argparse.ArgumentParser:
    parser, sub = build_root_parser("ag-wizard", DESCRIPTION)

    def add_target_args(p: argparse.ArgumentParser) -> None:
        source = p.add_mutually_exclusive_group()
        source.add_argument("--pdb-id", help="RCSB PDB entry id.")
        source.add_argument("--uniprot-id", help="UniProt accession.")
        source.add_argument("--input-file", help="Local PDB file.")
        p.add_argument("--chain", help="Keep only this chain.")
        p.add_argument("--chembl-target", help="ChEMBL target id for sourcing.")
        p.add_argument("--natural-product-query",
                       help="COCONUT search term for sourcing.")

    p = sub.add_parser("plan",
                       help="Describe the pipeline and check readiness.")
    add_target_args(p)
    p.add_argument("--stages", nargs="+", choices=sorted(STAGES),
                   help="Restrict the plan to these stages.")
    p.add_argument("--workdir", default="pipeline",
                   help="Where a run would write its stages.")
    p.add_argument("--reinvent-dir")
    add_output(p, help_text="Output JSON plan.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_plan)

    p = sub.add_parser("run", help="Execute the pipeline stages.")
    add_target_args(p)
    p.add_argument("--stages", nargs="+", choices=sorted(STAGES),
                   help="Stages to run (default: target sourcing triage).")
    p.add_argument("--workdir", default="pipeline",
                   help="Directory for stage outputs (default: pipeline).")
    p.add_argument("--keep-going", action="store_true",
                   help="Continue after a stage fails.")
    p.add_argument("--pchembl-min", type=float, default=6.0)
    p.add_argument("--limit", type=int, default=25,
                   help="Compounds to retrieve per source.")
    p.add_argument("--receptor", help="Receptor PDBQT, if not from the "
                                      "target stage.")
    p.add_argument("--site", help="Binding-site JSON, if not from the "
                                  "target stage.")
    p.add_argument("--dock-limit", type=int, default=10,
                   help="How many compounds to dock locally.")
    p.add_argument("--exhaustiveness", type=int, default=32)
    p.add_argument("--scoring-profile", default="drug-like",
                   help="REINVENT scoring profile for the generative stage.")
    p.add_argument("--num-steps", type=int, default=300)
    p.add_argument("--batch-size", type=int, default=128)
    p.add_argument("--device", default="cuda:0")
    p.add_argument("--reinvent-dir")
    p.add_argument("--force-field", default="charmm36")
    p.add_argument("--water", default="tip3p")
    p.add_argument("--temperature", type=float, default=310.0,
                   help="MD temperature in K (default: 310, physiological).")
    p.add_argument("--production-ns", type=float, default=100.0)
    add_output(p, help_text="Output JSON pipeline report.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_run)

    p = sub.add_parser("describe",
                       help="Document the stages, tiers and sub-skills.")
    add_output(p, help_text="Output JSON.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_describe)

    return parser


def main(argv: list[str] | None = None) -> int:
    return run_cli(build_parser(), argv)


if __name__ == "__main__":
    raise SystemExit(main())
