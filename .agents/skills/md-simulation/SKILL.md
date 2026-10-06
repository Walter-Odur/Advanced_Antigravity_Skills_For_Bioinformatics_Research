---
name: md-simulation
description: >
  Use this skill when the user wants to set up molecular dynamics
  simulations. Generates GROMACS workflows and MDP files with
  force-field-aware non-bonded settings, validates existing MDP files
  including cross-checking them against a declared force field, and
  produces SLURM submission scripts for GROMACS MD, REINVENT, AutoDock Vina
  batch screening, ColabFold and AlphaFold2. Do not use for receptor
  preparation or compound screening - use the dedicated skills instead.
---

# MD Simulation

## Overview

Generates production-ready GROMACS workflows, validates simulation
parameters, and creates cluster submission scripts.

**This skill produces files; it does not run simulations.** Production MD
is GPU work that belongs in a batch job.

**Do NOT use when:**

- The user wants to prepare a receptor - use `target-preparation`.
- The user wants to dock compounds - use `compound-screening`.
- The user wants to generate molecules - use `compound-synthesis`.
- The user wants to *run* a simulation here. Advise submitting the
  generated script.

## Prerequisites

| Requirement | Needed for |
|---|---|
| Python >= 3.10 | everything |
| GROMACS | *running* the generated workflow, not generating it |
| A SLURM cluster | submitting the generated scripts |

Generation needs nothing but Python. The outputs are self-contained.

## Core Rules

1. **Never hand-write an MDP file or a GROMACS command.** Use
   `gmxapi-setup`, which derives the non-bonded settings from the force
   field and validates everything it writes.
2. **`--output` is required** for every subcommand.
3. **Pass `--ff` to `gmxapi-validate`.** An MDP file does not record which
   force field it was written for, so the most consequential class of error
   is invisible without being told.
4. **Production runs go to the cluster.** Use `generate-hpc-script`.
5. **State the force field and temperature you used** when reporting.
6. If this skill is used, mention it in the output.

## Subcommands

| Subcommand | Purpose |
|---|---|
| `gmxapi-setup` | Generate the workflow, the shell script and all MDP files |
| `generate-md-setup` | Generate the setup shell script and MDP files only |
| `gmxapi-validate` | Validate MDP files, optionally against a force field |
| `generate-hpc-script` | SLURM script for GROMACS, REINVENT, Vina, ColabFold or AlphaFold2 |

Invoke as either:

```bash
ag-md-simulation <subcommand> --output <file> [options]
python .agents/skills/md-simulation/scripts/md_simulation_api.py <subcommand> ...
```

## Quick Start

```bash
ag-md-simulation gmxapi-setup \
    --pdb protein_clean.pdb \
    --ff charmm36 --water tip3p \
    --temperature 310 --production-ns 100 \
    --output md_workflow/run_md.py
```

This writes into `md_workflow/`:

- `run_md.py` - a gmxapi-based Python driver
- `setup_md.sh` - the equivalent shell script, which is what most clusters
  actually run
- `ions.mdp`, `em.mdp`, `nvt.mdp`, `npt.mdp`, `md.mdp`
- `run_md.meta.json` - the resolved settings and any warnings

Then:

```bash
bash md_workflow/setup_md.sh          # topology, box, solvate, EM, NVT, NPT
ag-md-simulation generate-hpc-script \
    --job-type gromacs-md --ngpu 4 --mem 128G --time 48:00:00 \
    --output md_workflow/submit_md.sh
sbatch md_workflow/submit_md.sh       # production
```

## Force Fields Are Not Interchangeable

Each force field was parameterised with a particular treatment of the
Lennard-Jones tail, and the non-bonded settings are part of the model, not
a preference. This skill derives them from `--ff`:

| Force field | Non-bonded treatment |
|---|---|
| `charmm36` | force-switched LJ from 1.0 to 1.2 nm, `rcoulomb` 1.2, **no** dispersion correction |
| `amber99sb`, `amber99sb-ildn`, `amber14sb` | plain 1.0 nm cutoff with `DispCorr = EnerPres` |
| `oplsaa` | plain 1.0 nm cutoff with `DispCorr`, TIP4P water |

Mixing them is a real error with a real consequence: applying AMBER's
`rvdw = 1.0` plus `DispCorr = EnerPres` to CHARMM36 double-counts a
dispersion tail that CHARMM36's force switch already absorbs into its
parameters, shifting pressure and density. The skill warns if the water
model does not match the one the force field was fitted with.

## Validating MDP Files

```bash
ag-md-simulation gmxapi-validate \
    --mdp em.mdp nvt.mdp npt.mdp md.mdp \
    --ff charmm36 \
    --output validation.json
```

Checks performed, each corresponding to a real grompp abort or a silent
physics error:

- `ref_t` and `tau_t` must have one value per `tc-grps` entry. A mismatch
  is one of the most common grompp failures.
- Pressure coupling requires `tau_p`, `ref_p` and `compressibility`.
- Position restraints combined with pressure coupling require
  `refcoord_scaling`; without it grompp warns and the restraint reference
  is not scaled with the box.
- `DispCorr` must not accompany a force-switched van der Waals potential.
- A timestep above 0.002 ps requires bond constraints; above 0.005 ps it is
  unstable outright.
- `gen_vel = yes` with `continuation = yes` contradict each other.
- Deprecated couplings are flagged: the Berendsen thermostat does not
  produce a correct canonical ensemble.
- With `--ff`, the non-bonded settings are cross-checked against that
  force field's requirements.

Output reports `valid`, `errors`, `warnings` and `notes` per file. The
notes include the simulated time implied by `nsteps` and `dt`, which is a
useful sanity check on its own.

## Simulation Parameters

| Flag | Default | Notes |
|---|---|---|
| `--ff` | `charmm36` | Determines the non-bonded settings |
| `--water` | `tip3p` | `tip4p` for OPLS-AA |
| `--temperature` | `300` | Kelvin. 310 is physiological |
| `--pressure` | `1.0` | bar |
| `--production-ns` | `100` | 0 generates equilibration only |
| `--equilibration-ps` | `100` | Per NVT and NPT stage |
| `--dt` | `0.002` | ps. Above 0.0025 needs hydrogen mass repartitioning |
| `--box-dist` | `1.2` | nm from solute to box edge |
| `--box-type` | `dodecahedron` | ~29% fewer waters than a cube for the same clearance |
| `--ion-concentration` | `0.15` | M, physiological |
| `--has-ligand` | off | Couple the ligand with the protein's temperature group |
| `--seed` | `-1` | Random. Set a positive integer for a reproducible run |

`--temperature` reaches `ref_t` and `gen_temp` in the MDP files, which is
what grompp reads. Verify it in the generated `md.mdp` if the simulation
temperature matters to the result.

## The Protocol

The generated workflow follows the standard sequence, with the reasoning
recorded in each MDP file's header:

1. **Topology** - `pdb2gmx` with the chosen force field and water model.
2. **Box and solvation** - `editconf`, then `solvate`.
3. **Neutralisation** - `genion` to the requested salt concentration.
4. **Minimisation** - steepest descent, to remove clashes from solvation.
5. **NVT equilibration** - brings the system to temperature with the solute
   restrained, so solvent relaxes around a fixed solute. V-rescale
   thermostat.
6. **NPT equilibration** - equilibrates density, solute still restrained.
   **C-rescale** barostat: it is stochastically correct and stable far from
   equilibrium, whereas Parrinello-Rahman rings when started from a freshly
   solvated box.
7. **Production** - unrestrained, Parrinello-Rahman barostat.

Two temperature-coupling groups are used throughout, which is what avoids
the "hot solvent, cold solute" artefact of coupling everything together.

## Protein-Ligand Systems

```bash
ag-md-simulation gmxapi-setup \
    --pdb complex.pdb --has-ligand \
    --ff charmm36 --temperature 310 --production-ns 100 \
    --output md_workflow/run_md.py
```

`--has-ligand` couples the ligand with the protein rather than the solvent
bath, and the output states that the `Protein_LIG` and `Water_and_ions`
index groups must be created with `gmx make_ndx` before grompp.

A ligand also needs its own topology, which GROMACS cannot generate:
CGenFF for CHARMM36, or ACPYPE/GAFF for AMBER. The skill says so rather
than producing a workflow that fails at `grompp`.

## Cluster Submission

```bash
# GROMACS production
ag-md-simulation generate-hpc-script \
    --job-type gromacs-md --ngpu 4 --mem 128G --time 48:00:00 \
    --output submit_md.sh

# REINVENT
ag-md-simulation generate-hpc-script \
    --job-type reinvent --config reinvent.toml \
    --gpu a100 --ngpu 1 --mem 64G --time 24:00:00 \
    --conda-env reinvent4 --output submit_reinvent.sh

# Vina batch screening, as a job array
ag-md-simulation generate-hpc-script \
    --job-type vina-screen --receptor receptor.pdbqt \
    --ligand-dir ligands/ --site-config binding_site.json \
    --array 1-100 --ngpu 0 --ncpus 16 --output submit_screen.sh

# Structure prediction
ag-md-simulation generate-hpc-script \
    --job-type colabfold --fasta target.fasta --output submit_fold.sh
```

Generated scripts are deliberately defensive, because a silent failure at
hour three of a 24-hour allocation is expensive:

- `set -euo pipefail`, so a failed step aborts.
- Inputs and executables are checked **before** the long step starts.
- `conda activate` goes through the shell hook, which a non-interactive
  SLURM shell otherwise lacks.
- GROMACS jobs use `-maxh` so a walltime limit produces a final checkpoint
  rather than a kill mid-write, and `-cpi` resumes from it automatically.
- Array jobs name their output per task (`%A_%a`), partition the ligand set
  by a disjoint stride so every ligand is docked exactly once, and skip
  work already done so a requeued task resumes.
- `--time`, `--mem` and array specifications are validated at generation
  time, so a typo fails immediately rather than at submission.

Review the `#SBATCH` directives against your cluster's partition and GPU
names before submitting.

## Reference

- gmxapi: <https://manual.gromacs.org/documentation/current/gmxapi/>
- GROMACS MDP options:
  <https://manual.gromacs.org/documentation/current/user-guide/mdp-options.html>
- CHARMM36: Huang & MacKerell, *J Comput Chem* 2013;34:2135-2145
- C-rescale barostat: Bernetti & Bussi, *J Chem Phys* 2020;153:114107
- V-rescale thermostat: Bussi, Donadio & Parrinello,
  *J Chem Phys* 2007;126:014101

## Workflow

1. Obtain a cleaned PDB, from `target-preparation` or elsewhere.
2. `gmxapi-setup` with the force field, temperature and length you want.
3. `gmxapi-validate --ff <the same force field>`, and read the warnings.
4. `bash setup_md.sh` to build and equilibrate the system.
5. `generate-hpc-script --job-type gromacs-md`, then `sbatch`.
6. Report the force field, water model, temperature and production length
   alongside any result.
