"""CLI for the compound-screening skill."""

from __future__ import annotations

import argparse
from pathlib import Path

from ..chem.admet import TIERS, screen_compounds, tier_summary
from ..errors import UsageError
from ..io_utils import read_json
from ..libraries.chembl import ASSAY_TYPES, query_activities, search_targets
from ..libraries.coconut import search_coconut
from ..libraries.zinc import (
    PURCHASABILITY_CODES,
    REACTIVITY_CODES,
    SUBSETS,
    generate_zinc_script,
    select_tranches,
    verify_tranches,
)
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

_TIER_COUNTS = tier_summary()

DESCRIPTION = f"""\
Search compound libraries, dock molecules, and screen them on ADMET.

Two-tier screening:
  --strictness relaxed   {_TIER_COUNTS['relaxed']} checks, post-docking triage
                         and seed selection for generative chemistry
  --strictness strict    {_TIER_COUNTS['strict']} checks, final candidates
                         (the default)

The check counts above are derived from the check registry, so they always
match what is actually run. Every check is reported individually in the
output, with its limit and its literature citation.
"""


# ---------------------------------------------------------------------------
# Handlers
# ---------------------------------------------------------------------------


def cmd_admet_filter(args: argparse.Namespace) -> None:
    records = compounds_from_args(args)
    result = screen_compounds(records, strictness=args.strictness)
    if args.passing_only:
        result["compounds"] = [c for c in result["compounds"]
                               if c.get("overall_pass")]
    if not args.include_check_detail:
        # The per-check array is verbose; keep the failure reasons.
        for compound in result["compounds"]:
            compound.pop("checks", None)
        result["tier_definition"].pop("checks", None)

    lines = [
        f"ADMET screening: {result['tier']} "
        f"({result['checks_per_compound']} checks per compound)",
        f"  total {result['total_compounds']}  "
        f"passed {result['passed']}  failed {result['failed']}  "
        f"errored {result['errored']}",
    ]
    for compound in result["compounds"][:args.show]:
        properties = compound.get("properties", {})
        if compound["status"] == "error":
            lines.append(f"  [ERROR] {compound['name']}: {compound['error']}")
            continue
        mark = "PASS" if compound["overall_pass"] else "FAIL"
        lines.append(
            f"  [{mark}] {compound['name']:<22s} "
            f"MW {properties.get('mw', '?'):>6} "
            f"logP {properties.get('logp', '?'):>5} "
            f"QED {properties.get('qed', '?'):>5} "
            f"SA {properties.get('sa_score', '?'):>4} "
            f"({compound['checks_passed']}/{result['checks_per_compound']})"
        )
        if not compound["overall_pass"]:
            for reason in compound["failed_checks"][:3]:
                lines.append(f"            - {reason}")
    emit_json(result, args, summary=lines)


def cmd_describe_tiers(args: argparse.Namespace) -> None:
    payload = {
        "tiers": {name: tier.describe() for name, tier in TIERS.items()},
        "check_counts": tier_summary(),
        "note": ("These definitions are read from the check registry, so the "
                 "counts reported here are exactly what 'admet-filter' "
                 "performs."),
    }
    lines = []
    for name, tier in TIERS.items():
        lines.append(f"{tier.label} ({name}): {tier.n_checks} checks")
        lines.append(f"  {tier.purpose}")
        for entry in tier.describe()["checks"]:
            policy = entry["policy"]
            rendered = ", ".join(f"{k}={v}" for k, v in policy.items()) or "-"
            lines.append(f"    {entry['check']:<20s} {entry['category']:<22s} "
                         f"{rendered}")
    emit_json(payload, args, summary=lines)


def cmd_search_chembl_target(args: argparse.Namespace) -> None:
    result = search_targets(args.query, limit=args.limit,
                            organism=args.organism)
    lines = [f"ChEMBL target search: {args.query!r} "
             f"-> {result['total_returned']} target(s)"]
    for target in result["targets"]:
        lines.append(
            f"  {target['target_chembl_id']:<16s} "
            f"{(target['organism'] or '-')[:30]:<30s} "
            f"{(target['pref_name'] or '-')[:40]}"
        )
        if target["uniprot_accessions"]:
            lines.append(f"      UniProt: "
                         f"{', '.join(target['uniprot_accessions'])}")
    emit_json(result, args, summary=lines)


def cmd_query_chembl(args: argparse.Namespace) -> None:
    result = query_activities(args.target_id, pchembl_min=args.pchembl_min,
                              assay_type=args.assay_type, limit=args.limit)
    lines = [
        f"ChEMBL activities for {args.target_id} (pChEMBL >= {args.pchembl_min})",
        f"  examined {result['activities_examined']} activity records over "
        f"{result['pages_fetched']} page(s)",
        f"  {result['distinct_molecules_found']} distinct molecules, "
        f"returning {result['returned']}",
    ]
    for compound in result["compounds"][:args.show]:
        lines.append(
            f"  {compound['molecule_chembl_id']:<16s} "
            f"pChEMBL {compound['pchembl_median']:>5.2f} "
            f"(n={compound['n_measurements']}) "
            f"{compound['canonical_smiles'][:46]}"
        )
    emit_json(result, args, summary=lines)


def cmd_query_coconut(args: argparse.Namespace) -> None:
    result = search_coconut(query=args.query, smiles=args.smiles,
                            limit=args.limit, page=args.page,
                            max_pages=args.max_pages)
    lines = [
        f"COCONUT {result['search_type']} search: {result['query']!r}",
        f"  returned {result['returned']} ({result['with_smiles']} with a "
        f"structure) of {result['total_in_database']} in the database",
    ]
    for compound in result["compounds"][:args.show]:
        lines.append(f"  {compound['id']:<18s} "
                     f"{(compound['name'] or '-')[:28]:<28s} "
                     f"{compound['smiles'][:40]}")
    emit_json(result, args, summary=lines)


def cmd_query_zinc(args: argparse.Namespace) -> None:
    selection = select_tranches(
        args.subset,
        mw_range=tuple(args.mw_range) if args.mw_range else None,
        logp_range=tuple(args.logp_range) if args.logp_range else None,
        reactivity=args.reactivity, purchasability=args.purchasability,
        fmt=args.format,
    )
    verification = None
    if args.verify:
        verification = verify_tranches(selection, limit=args.verify_limit)
    script = generate_zinc_script(selection, out_dir=args.download_dir)

    lines = [
        f"ZINC tranche selection: {selection.subset}",
        f"  MW {selection.mw_range[0]:g}-{selection.mw_range[1]:g} Da "
        f"(letters {''.join(selection.size_letters)})",
        f"  logP {selection.logp_range[0]:g}-{selection.logp_range[1]:g} "
        f"(letters {''.join(selection.logp_letters)})",
        f"  {len(selection.codes)} tranches, format {selection.fmt}",
        f"  example URL: {selection.urls[0][1]}",
    ]
    if verification:
        lines.append(f"  verified {verification['reachable']}/"
                     f"{verification['checked']} tranche URLs reachable")
        for entry in verification["unreachable"][:5]:
            lines.append(f"    unreachable: {entry['code']} - {entry['reason']}")
    lines.append(f"  run it with: bash {args.output}")

    emit_text(script, args, executable=True, summary=lines)
    # The selection metadata is useful downstream, so write it too.
    from ..io_utils import write_json
    metadata = Path(args.output).with_suffix(".json")
    write_json(selection.as_dict(), metadata)
    if not args.quiet:
        print(f"Written: {metadata}")


def cmd_dock(args: argparse.Namespace) -> None:
    from ..docking.vina import dock_compounds
    from ..targets.site import BindingSite

    records = compounds_from_args(args)
    if args.site:
        data = read_json(args.site)
        try:
            site = BindingSite(
                center_x=float(data["center_x"]), center_y=float(data["center_y"]),
                center_z=float(data["center_z"]), size_x=float(data["size_x"]),
                size_y=float(data["size_y"]), size_z=float(data["size_z"]),
                method=data.get("method", "from file"),
                reference=data.get("reference", str(args.site)),
                confidence=data.get("confidence", ""),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise UsageError(
                f"{args.site} is not a binding-site file: {exc}",
                hint="Produce one with 'ag-target-preparation define-site'.",
            ) from exc
    elif args.center:
        site = BindingSite(
            center_x=args.center[0], center_y=args.center[1],
            center_z=args.center[2], size_x=args.box_size[0],
            size_y=args.box_size[1], size_z=args.box_size[2],
            method="supplied on the command line",
            reference="--center/--box-size",
            confidence="unknown - supplied by the user",
        )
    else:
        raise UsageError(
            "A docking box is required.",
            hint="Pass --site (from 'ag-target-preparation define-site') or "
                 "--center X Y Z with --box-size X Y Z.",
        )

    result = dock_compounds(
        records, args.receptor, site,
        output_dir=args.output_dir or Path(args.output).parent,
        exhaustiveness=args.exhaustiveness, n_poses=args.n_poses,
        cpu=args.cpu, seed=args.seed, timeout=args.timeout,
    )
    lines = [
        f"Docking: {result['engine']}",
        f"  receptor {Path(result['receptor']).name}, "
        f"box {site.size_x:.0f}x{site.size_y:.0f}x{site.size_z:.0f} A",
        f"  {result['docked']}/{result['total_compounds']} docked, "
        f"{result['failed']} failed",
    ]
    for entry in result["ranked_results"][:args.show]:
        affinity = entry.get("best_affinity_kcal_mol")
        if affinity is None:
            lines.append(f"  [{entry['status']}] {entry['name']}: "
                         f"{entry.get('error', '')[:60]}")
        else:
            lines.append(f"  {affinity:>7.2f} kcal/mol  {entry['name']}")
    emit_json(result, args, summary=lines)


def cmd_admet_predict(args: argparse.Namespace) -> None:
    """Online ADMET enrichment, with an honest account of availability."""
    records = compounds_from_args(args)
    result = screen_compounds(records, strictness=args.strictness)
    result["online_prediction"] = {
        "attempted": False,
        "reason": (
            "No public ADMETlab 3.0 REST endpoint is currently reachable. "
            "The documented path /server/api/aio returned HTTP 404 when "
            "checked on 2026-10-06, as did seven other candidate paths, "
            "while the web interface itself responded normally. Rather than "
            "report a failure as a result, this command runs the offline "
            "RDKit battery, which is deterministic and reproducible."
        ),
        "how_to_add_it": (
            "If you have access to an ADMET prediction service, set "
            "AGSKILLS_ADMET_URL to its endpoint; this command will then POST "
            "{'smiles': [...]} to it and merge the response under "
            "'online_properties'."
        ),
        "web_interface": "https://admetlab3.scbdd.com/",
    }
    emit_json(result, args, summary=[
        "ADMET prediction (offline RDKit engine)",
        f"  {result['tier']}: {result['checks_per_compound']} checks",
        f"  total {result['total_compounds']}  passed {result['passed']}  "
        f"failed {result['failed']}  errored {result['errored']}",
        "  note: no public ADMETlab REST endpoint is reachable; see "
        "'online_prediction' in the output.",
    ])


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser, sub = build_root_parser("ag-compound-screening", DESCRIPTION)

    p = sub.add_parser("admet-filter",
                       help="Screen compounds on the two-tier ADMET battery.")
    add_compound_inputs(p)
    p.add_argument("--strictness", default="strict",
                   choices=sorted(TIERS),
                   help="Screening tier (default: strict).")
    p.add_argument("--passing-only", action="store_true",
                   help="Write only the compounds that passed.")
    p.add_argument("--include-check-detail", action="store_true",
                   help="Include every individual check in the output.")
    p.add_argument("--show", type=int, default=12,
                   help="How many compounds to summarise on stdout.")
    add_output(p, help_text="Output JSON report.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_admet_filter)

    p = sub.add_parser("describe-tiers",
                       help="Print exactly which checks each tier performs.")
    add_output(p, help_text="Output JSON tier definitions.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_describe_tiers)

    p = sub.add_parser("search-chembl-target",
                       help="Find a ChEMBL target id by name.")
    p.add_argument("--query", required=True, help="Target name or gene symbol.")
    p.add_argument("--organism", help="Filter by organism substring.")
    p.add_argument("--limit", type=int, default=10)
    p.add_argument("--show", type=int, default=10)
    add_output(p, help_text="Output JSON.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_search_chembl_target)

    p = sub.add_parser("query-chembl",
                       help="Retrieve known actives for a ChEMBL target.")
    p.add_argument("--target-id", required=True, help="e.g. CHEMBL203")
    p.add_argument("--pchembl-min", type=float, default=6.0,
                   help="Minimum pChEMBL; 6 is 1 uM (default: 6).")
    p.add_argument("--assay-type", default="B", choices=sorted(ASSAY_TYPES),
                   help="B=binding, F=functional (default: B).")
    p.add_argument("--limit", type=int, default=25)
    p.add_argument("--show", type=int, default=10)
    add_output(p, help_text="Output JSON.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_query_chembl)

    p = sub.add_parser("query-coconut",
                       help="Search COCONUT natural products.")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--query", help="Free-text term.")
    group.add_argument("--smiles", help="SMILES for a structure search.")
    p.add_argument("--limit", type=int, default=25)
    p.add_argument("--page", type=int, default=1)
    p.add_argument("--max-pages", type=int, default=1)
    p.add_argument("--show", type=int, default=10)
    add_output(p, help_text="Output JSON.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_query_coconut)

    p = sub.add_parser("query-zinc",
                       help="Generate a verified ZINC bulk download script.")
    p.add_argument("--subset", default="drug-like", choices=sorted(SUBSETS),
                   help="Named property window (default: drug-like).")
    p.add_argument("--mw-range", nargs=2, type=float, metavar=("MIN", "MAX"))
    p.add_argument("--logp-range", nargs=2, type=float, metavar=("MIN", "MAX"))
    p.add_argument("--reactivity", default="A", choices=sorted(REACTIVITY_CODES),
                   help="Reactivity filter (default: A, most conservative).")
    p.add_argument("--purchasability", default="A",
                   choices=sorted(PURCHASABILITY_CODES),
                   help="Availability (default: A, in-stock).")
    p.add_argument("--format", default="smi", choices=["smi", "sdf", "mol2"])
    p.add_argument("--download-dir",
                   help="Directory the generated script downloads into.")
    p.add_argument("--verify", action="store_true",
                   help="HEAD-check the tranche URLs before writing the script.")
    p.add_argument("--verify-limit", type=int, default=6,
                   help="How many tranches to check with --verify (0 = all).")
    add_output(p, help_text="Output .sh download script. Selection metadata "
                            "is written alongside it as .json.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_query_zinc)

    p = sub.add_parser("dock", help="Dock compounds with AutoDock Vina.")
    p.add_argument("--receptor", required=True, help="Receptor PDBQT file.")
    add_compound_inputs(p)
    box = p.add_mutually_exclusive_group(required=True)
    box.add_argument("--site", help="Binding-site JSON from 'define-site'.")
    box.add_argument("--center", nargs=3, type=float, metavar=("X", "Y", "Z"),
                     help="Box centre in Angstroms.")
    p.add_argument("--box-size", nargs=3, type=float, default=[22.0, 22.0, 22.0],
                   metavar=("X", "Y", "Z"),
                   help="Box dimensions for --center (default: 22 22 22).")
    p.add_argument("--exhaustiveness", type=int, default=32,
                   help="Vina search effort (default: 32).")
    p.add_argument("--n-poses", type=int, default=9)
    p.add_argument("--cpu", type=int, default=0, help="0 lets Vina choose.")
    p.add_argument("--seed", type=int, default=42,
                   help="Fixed seed, so a run is reproducible (default: 42).")
    p.add_argument("--timeout", type=float, default=600.0,
                   help="Per-compound limit in seconds.")
    p.add_argument("--output-dir", help="Where pose files are written.")
    p.add_argument("--show", type=int, default=12)
    add_output(p, help_text="Output JSON results.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_dock)

    p = sub.add_parser("admet-predict",
                       help="Full ADMET profile (offline engine; see output "
                            "for online availability).")
    add_compound_inputs(p)
    p.add_argument("--strictness", default="strict", choices=sorted(TIERS))
    p.add_argument("--passing-only", action="store_true")
    p.add_argument("--include-check-detail", action="store_true", default=True)
    p.add_argument("--show", type=int, default=12)
    add_output(p, help_text="Output JSON report.")
    add_common_flags(p)
    p.set_defaults(handler=cmd_admet_predict)

    return parser


def main(argv: list[str] | None = None) -> int:
    return run_cli(build_parser(), argv)


if __name__ == "__main__":
    raise SystemExit(main())
