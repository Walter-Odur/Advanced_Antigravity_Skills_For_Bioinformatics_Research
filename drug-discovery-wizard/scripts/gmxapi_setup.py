"""gmxapi Workflow Generator for Drug Discovery Wizard.

Generates Python scripts using the gmxapi package for programmatic GROMACS
molecular dynamics control. Falls back to shell script generation when
gmxapi is not available.

Reference: https://manual.gromacs.org/documentation/current/gmxapi/

Usage:
    python gmxapi_setup.py setup --pdb protein.pdb --ff charmm36 --water tip3p \
        --production-ns 100 --output md_workflow.py
    python gmxapi_setup.py mdp --type production --ns 100 --output md.mdp
    python gmxapi_setup.py validate --mdp em.mdp nvt.mdp npt.mdp md.mdp --output report.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from textwrap import dedent
from typing import Any

# ---------------------------------------------------------------------------
# MDP Templates (validated against GROMACS 2024+ defaults)
# ---------------------------------------------------------------------------

MDP_TEMPLATES: dict[str, dict[str, Any]] = {
    "ions": {
        "desc": "Ion placement (genion preparation)",
        "params": {
            "integrator": "steep",
            "emtol": 1000.0,
            "emstep": 0.01,
            "nsteps": 50000,
            "nstlist": 1,
            "cutoff-scheme": "Verlet",
            "ns_type": "grid",
            "coulombtype": "cutoff",
            "rcoulomb": 1.0,
            "rvdw": 1.0,
            "pbc": "xyz",
        },
    },
    "em": {
        "desc": "Energy minimization (steepest descent)",
        "params": {
            "integrator": "steep",
            "emtol": 1000.0,
            "emstep": 0.01,
            "nsteps": 50000,
            "nstlist": 10,
            "cutoff-scheme": "Verlet",
            "ns_type": "grid",
            "coulombtype": "PME",
            "rcoulomb": 1.0,
            "rvdw": 1.0,
            "pbc": "xyz",
        },
    },
    "nvt": {
        "desc": "NVT equilibration (100 ps, V-rescale thermostat)",
        "params": {
            "define": "-DPOSRES",
            "integrator": "md",
            "nsteps": 50000,
            "dt": 0.002,
            "nstxout-compressed": 5000,
            "nstenergy": 5000,
            "nstlog": 5000,
            "continuation": "no",
            "constraint_algorithm": "lincs",
            "constraints": "h-bonds",
            "lincs_iter": 1,
            "lincs_order": 4,
            "cutoff-scheme": "Verlet",
            "nstlist": 10,
            "rcoulomb": 1.0,
            "rvdw": 1.0,
            "coulombtype": "PME",
            "pme_order": 4,
            "fourierspacing": 0.16,
            "tcoupl": "V-rescale",
            "tc-grps": "Protein Non-Protein",
            "tau_t": "0.1 0.1",
            "ref_t": "300 300",
            "pcoupl": "no",
            "pbc": "xyz",
            "DispCorr": "EnerPres",
            "gen_vel": "yes",
            "gen_temp": 300,
            "gen_seed": -1,
        },
    },
    "npt": {
        "desc": "NPT equilibration (100 ps, Parrinello-Rahman barostat)",
        "params": {
            "define": "-DPOSRES",
            "integrator": "md",
            "nsteps": 50000,
            "dt": 0.002,
            "nstxout-compressed": 5000,
            "nstenergy": 5000,
            "nstlog": 5000,
            "continuation": "yes",
            "constraint_algorithm": "lincs",
            "constraints": "h-bonds",
            "lincs_iter": 1,
            "lincs_order": 4,
            "cutoff-scheme": "Verlet",
            "nstlist": 10,
            "rcoulomb": 1.0,
            "rvdw": 1.0,
            "coulombtype": "PME",
            "pme_order": 4,
            "fourierspacing": 0.16,
            "tcoupl": "V-rescale",
            "tc-grps": "Protein Non-Protein",
            "tau_t": "0.1 0.1",
            "ref_t": "300 300",
            "pcoupl": "Parrinello-Rahman",
            "pcoupltype": "isotropic",
            "tau_p": 2.0,
            "ref_p": 1.0,
            "compressibility": "4.5e-5",
            "pbc": "xyz",
            "DispCorr": "EnerPres",
            "gen_vel": "no",
        },
    },
    "production": {
        "desc": "Production MD run",
        "params": {
            "integrator": "md",
            "nsteps": None,  # Set dynamically from --production-ns
            "dt": 0.002,
            "nstxout-compressed": 5000,
            "nstenergy": 5000,
            "nstlog": 5000,
            "continuation": "yes",
            "constraint_algorithm": "lincs",
            "constraints": "h-bonds",
            "lincs_iter": 1,
            "lincs_order": 4,
            "cutoff-scheme": "Verlet",
            "nstlist": 10,
            "rcoulomb": 1.0,
            "rvdw": 1.0,
            "coulombtype": "PME",
            "pme_order": 4,
            "fourierspacing": 0.16,
            "tcoupl": "V-rescale",
            "tc-grps": "Protein Non-Protein",
            "tau_t": "0.1 0.1",
            "ref_t": "300 300",
            "pcoupl": "Parrinello-Rahman",
            "pcoupltype": "isotropic",
            "tau_p": 2.0,
            "ref_p": 1.0,
            "compressibility": "4.5e-5",
            "pbc": "xyz",
            "DispCorr": "EnerPres",
            "gen_vel": "no",
        },
    },
}

# Valid GROMACS parameters for validation
VALID_INTEGRATORS = {"md", "md-vv", "md-vv-avek", "steep", "cg", "l-bfgs", "nm", "tpi", "tpic", "sd"}
VALID_COULOMBTYPES = {"Cut-off", "cutoff", "Ewald", "PME", "P3M-AD", "Reaction-Field"}
VALID_TCOUPL = {"no", "berendsen", "nose-hoover", "andersen", "andersen-massive", "v-rescale", "V-rescale"}
VALID_PCOUPL = {"no", "berendsen", "Parrinello-Rahman", "MTTK", "C-rescale"}
VALID_CONSTRAINTS = {"none", "h-bonds", "all-bonds", "h-angles", "all-angles"}


# ---------------------------------------------------------------------------
# MDP Generation
# ---------------------------------------------------------------------------

def generate_mdp(mdp_type: str, production_ns: float = 100.0) -> str:
    """Generate a valid MDP parameter file.

    Args:
        mdp_type: One of ions, em, nvt, npt, production.
        production_ns: Production run length in nanoseconds.

    Returns:
        MDP file content as string.
    """
    if mdp_type not in MDP_TEMPLATES:
        raise ValueError(f"Unknown MDP type '{mdp_type}'. Choose from: {list(MDP_TEMPLATES.keys())}")

    template = MDP_TEMPLATES[mdp_type]
    params = dict(template["params"])

    # Calculate production nsteps
    if mdp_type == "production" and params.get("nsteps") is None:
        dt = params.get("dt", 0.002)
        params["nsteps"] = int(production_ns * 1_000_000 / (dt * 1000))

    lines = [f"; {template['desc']}", f"; Generated by Drug Discovery Wizard gmxapi setup", ""]
    for key, val in params.items():
        if val is not None:
            lines.append(f"{key:<25s} = {val}")

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# gmxapi Workflow Generation
# ---------------------------------------------------------------------------

def generate_gmxapi_workflow(
    pdb: str,
    ff: str = "charmm36",
    water: str = "tip3p",
    production_ns: float = 100.0,
    box_dist: float = 1.2,
    temperature: float = 300.0,
) -> str:
    """Generate a complete gmxapi Python workflow script.

    This script uses gmxapi.commandline_operation() to wrap GROMACS commands
    in a Python-native workflow with proper error handling.
    """
    nsteps_prod = int(production_ns * 1_000_000 / 2)  # dt=0.002 ps

    script = dedent(f'''\
        #!/usr/bin/env python3
        """GROMACS MD Workflow using gmxapi.

        Generated by Drug Discovery Wizard.
        Reference: https://manual.gromacs.org/documentation/current/gmxapi/

        Workflow:
            1. Topology generation (pdb2gmx)
            2. Box definition (editconf)
            3. Solvation (solvate)
            4. Ion placement (genion)
            5. Energy minimization
            6. NVT equilibration (100 ps)
            7. NPT equilibration (100 ps)
            8. Production MD ({production_ns} ns)

        Requirements:
            pip install gmxapi
            GROMACS must be installed and on PATH.
        """

        import os
        import sys
        import subprocess
        from pathlib import Path

        try:
            import gmxapi as gmx
            USE_GMXAPI = True
            print(f"[OK] gmxapi version: {{gmx.__version__}}")
        except ImportError:
            USE_GMXAPI = False
            print("[WARNING] gmxapi not available. Using subprocess fallback.")

        # ---------------------------------------------------------------
        # Configuration
        # ---------------------------------------------------------------
        PDB_FILE = "{pdb}"
        FF = "{ff}"
        WATER = "{water}"
        BOX_DIST = {box_dist}
        TEMPERATURE = {temperature}
        PRODUCTION_NS = {production_ns}
        PRODUCTION_NSTEPS = {nsteps_prod}


        def run_gmx(program: str, args: list[str], stdin_text: str | None = None):
            """Run a GROMACS command, using gmxapi if available."""
            if USE_GMXAPI:
                cmd = gmx.commandline_operation(
                    executable="gmx",
                    arguments=[program] + args,
                )
                cmd.run()
                if cmd.output.returncode.result() != 0:
                    raise RuntimeError(f"gmx {{program}} failed")
                return cmd
            else:
                full_cmd = ["gmx", program] + args
                proc = subprocess.run(
                    full_cmd,
                    input=stdin_text,
                    capture_output=True,
                    text=True,
                )
                if proc.returncode != 0:
                    print(proc.stderr, file=sys.stderr)
                    raise RuntimeError(f"gmx {{program}} failed (rc={{proc.returncode}})")
                return proc


        def main():
            print("=" * 60)
            print(f"GROMACS MD Workflow: {{PDB_FILE}}")
            print(f"Force field: {{FF}} | Water: {{WATER}} | Production: {{PRODUCTION_NS}} ns")
            print("=" * 60)

            # 1. Generate topology
            print("\\n[Step 1/8] Generating topology...")
            run_gmx("pdb2gmx", [
                "-f", PDB_FILE,
                "-o", "protein.gro",
                "-p", "topol.top",
                "-ignh",
                "-ff", FF,
                "-water", WATER,
            ])

            # 2. Define simulation box
            print("[Step 2/8] Defining simulation box...")
            run_gmx("editconf", [
                "-f", "protein.gro",
                "-o", "box.gro",
                "-c",
                "-d", str(BOX_DIST),
                "-bt", "dodecahedron",
            ])

            # 3. Solvate
            print("[Step 3/8] Adding solvent...")
            run_gmx("solvate", [
                "-cp", "box.gro",
                "-cs", "spc216.gro",
                "-o", "solv.gro",
                "-p", "topol.top",
            ])

            # 4. Add ions
            print("[Step 4/8] Adding ions...")
            run_gmx("grompp", [
                "-f", "ions.mdp",
                "-c", "solv.gro",
                "-p", "topol.top",
                "-o", "ions.tpr",
                "-maxwarn", "1",
            ])
            run_gmx("genion", [
                "-s", "ions.tpr",
                "-o", "solv_ions.gro",
                "-p", "topol.top",
                "-pname", "NA",
                "-nname", "CL",
                "-neutral",
            ], stdin_text="SOL\\n" if not USE_GMXAPI else None)

            # 5. Energy minimization
            print("[Step 5/8] Energy minimization...")
            run_gmx("grompp", [
                "-f", "em.mdp",
                "-c", "solv_ions.gro",
                "-p", "topol.top",
                "-o", "em.tpr",
            ])
            if USE_GMXAPI:
                em_input = gmx.read_tpr("em.tpr")
                md_em = gmx.mdrun(em_input)
                md_em.run()
            else:
                run_gmx("mdrun", ["-v", "-deffnm", "em"])

            # 6. NVT equilibration
            print("[Step 6/8] NVT equilibration (100 ps)...")
            run_gmx("grompp", [
                "-f", "nvt.mdp",
                "-c", "em.gro",
                "-r", "em.gro",
                "-p", "topol.top",
                "-o", "nvt.tpr",
            ])
            if USE_GMXAPI:
                nvt_input = gmx.read_tpr("nvt.tpr")
                md_nvt = gmx.mdrun(nvt_input)
                md_nvt.run()
            else:
                run_gmx("mdrun", ["-v", "-deffnm", "nvt"])

            # 7. NPT equilibration
            print("[Step 7/8] NPT equilibration (100 ps)...")
            run_gmx("grompp", [
                "-f", "npt.mdp",
                "-c", "nvt.gro",
                "-r", "nvt.gro",
                "-t", "nvt.cpt",
                "-p", "topol.top",
                "-o", "npt.tpr",
            ])
            if USE_GMXAPI:
                npt_input = gmx.read_tpr("npt.tpr")
                md_npt = gmx.mdrun(npt_input)
                md_npt.run()
            else:
                run_gmx("mdrun", ["-v", "-deffnm", "npt"])

            # 8. Production MD
            print(f"[Step 8/8] Production MD ({{PRODUCTION_NS}} ns)...")
            run_gmx("grompp", [
                "-f", "md.mdp",
                "-c", "npt.gro",
                "-t", "npt.cpt",
                "-p", "topol.top",
                "-o", "md.tpr",
            ])
            if USE_GMXAPI:
                md_input = gmx.read_tpr("md.tpr")
                md_prod = gmx.mdrun(md_input)
                md_prod.run()
            else:
                run_gmx("mdrun", ["-v", "-deffnm", "md"])

            print("\\n" + "=" * 60)
            print("WORKFLOW COMPLETE")
            print(f"Trajectory: md.xtc | Energy: md.edr | Log: md.log")
            print("=" * 60)


        if __name__ == "__main__":
            main()
    ''')

    return script


# ---------------------------------------------------------------------------
# MDP Validation
# ---------------------------------------------------------------------------

def validate_mdp(mdp_path: str) -> dict:
    """Validate an MDP file for common issues."""
    issues: list[str] = []
    warnings: list[str] = []
    params: dict[str, str] = {}

    with open(mdp_path) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith(";"):
                continue
            if "=" in line:
                key, _, val = line.partition("=")
                key = key.strip()
                val = val.strip().split(";")[0].strip()
                params[key] = val

    # Check integrator
    integrator = params.get("integrator", "").lower()
    if integrator and integrator not in VALID_INTEGRATORS:
        issues.append(f"Invalid integrator '{integrator}'. Valid: {VALID_INTEGRATORS}")

    # Check coulombtype
    ct = params.get("coulombtype", "")
    if ct and ct not in VALID_COULOMBTYPES:
        issues.append(f"Invalid coulombtype '{ct}'. Valid: {VALID_COULOMBTYPES}")

    # Check thermostat
    tc = params.get("tcoupl", "")
    if tc and tc not in VALID_TCOUPL:
        issues.append(f"Invalid tcoupl '{tc}'. Valid: {VALID_TCOUPL}")

    # Check barostat
    pc = params.get("pcoupl", "")
    if pc and pc not in VALID_PCOUPL:
        issues.append(f"Invalid pcoupl '{pc}'. Valid: {VALID_PCOUPL}")

    # Check constraints
    con = params.get("constraints", "")
    if con and con not in VALID_CONSTRAINTS:
        issues.append(f"Invalid constraints '{con}'. Valid: {VALID_CONSTRAINTS}")

    # Check nsteps
    nsteps = params.get("nsteps")
    if nsteps:
        try:
            n = int(nsteps)
            if n <= 0:
                issues.append(f"nsteps must be positive, got {n}")
            elif n > 500_000_000:
                warnings.append(f"Very large nsteps ({n}). Consider splitting into chunks.")
        except ValueError:
            issues.append(f"nsteps must be integer, got '{nsteps}'")

    # Check dt
    dt = params.get("dt")
    if dt:
        try:
            d = float(dt)
            if d <= 0 or d > 0.005:
                warnings.append(f"Unusual dt={d}. Standard is 0.001-0.002 ps.")
        except ValueError:
            issues.append(f"dt must be numeric, got '{dt}'")

    return {
        "file": mdp_path,
        "parameters": params,
        "valid": len(issues) == 0,
        "issues": issues,
        "warnings": warnings,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="gmxapi Workflow Generator — programmatic GROMACS MD"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # --- setup ---
    p_setup = sub.add_parser("setup", help="Generate gmxapi Python workflow script")
    p_setup.add_argument("--pdb", required=True, help="Input PDB file")
    p_setup.add_argument("--ff", default="charmm36", help="Force field (default: charmm36)")
    p_setup.add_argument("--water", default="tip3p", help="Water model (default: tip3p)")
    p_setup.add_argument("--production-ns", type=float, default=100.0, help="Production run (ns)")
    p_setup.add_argument("--box-dist", type=float, default=1.2, help="Box distance (nm)")
    p_setup.add_argument("--temperature", type=float, default=300.0, help="Temperature (K)")
    p_setup.add_argument("--output", required=True, help="Output Python script path")

    # --- mdp ---
    p_mdp = sub.add_parser("mdp", help="Generate MDP parameter file")
    p_mdp.add_argument("--type", required=True, choices=list(MDP_TEMPLATES.keys()), help="MDP type")
    p_mdp.add_argument("--ns", type=float, default=100.0, help="Production length (ns)")
    p_mdp.add_argument("--output", required=True, help="Output MDP path")

    # --- validate ---
    p_val = sub.add_parser("validate", help="Validate MDP file(s)")
    p_val.add_argument("--mdp", nargs="+", required=True, help="MDP files to validate")
    p_val.add_argument("--output", required=True, help="Output JSON report")

    args = parser.parse_args()

    if args.command == "setup":
        script = generate_gmxapi_workflow(
            pdb=args.pdb,
            ff=args.ff,
            water=args.water,
            production_ns=args.production_ns,
            box_dist=args.box_dist,
            temperature=args.temperature,
        )
        # Also generate MDP files alongside the script
        out_dir = Path(args.output).parent
        for mdp_type in ["ions", "em", "nvt", "npt", "production"]:
            mdp_content = generate_mdp(mdp_type, args.production_ns)
            mdp_name = "md.mdp" if mdp_type == "production" else f"{mdp_type}.mdp"
            mdp_path = out_dir / mdp_name
            mdp_path.write_text(mdp_content)
            print(f"[OK] Generated {mdp_path}")

        Path(args.output).write_text(script)
        print(f"[OK] gmxapi workflow script saved to {args.output}")

    elif args.command == "mdp":
        content = generate_mdp(args.type, args.ns)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text(content)
        print(f"[OK] MDP file saved to {args.output}")

    elif args.command == "validate":
        results = []
        for mdp_file in args.mdp:
            results.append(validate_mdp(mdp_file))
        report = {
            "total_files": len(results),
            "all_valid": all(r["valid"] for r in results),
            "files": results,
        }
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            json.dump(report, f, indent=2)
        print(f"[OK] Validation report saved to {args.output}")
        for r in results:
            status = "VALID" if r["valid"] else "INVALID"
            print(f"  {r['file']}: {status} ({len(r['issues'])} issues, {len(r['warnings'])} warnings)")


if __name__ == "__main__":
    main()
