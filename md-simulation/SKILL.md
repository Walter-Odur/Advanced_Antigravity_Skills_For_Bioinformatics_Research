---
name: md-simulation
description: >
  Use this skill when the user wants to set up molecular dynamics simulations.
  Handles GROMACS gmxapi Python workflow generation, MDP file creation and
  validation, shell script fallback, and SLURM HPC job script generation
  for GROMACS MD, REINVENT, and Vina batch screening. Do not use for receptor
  preparation or compound screening — use the dedicated skills instead.
---

# MD Simulation

## Prerequisites

1.  **Python 3.10+** available.
2.  **GROMACS** (optional): Not required for script *generation*, but required
    to *run* the generated workflows. The scripts produce ready-to-submit files.
3.  **HPC cluster with SLURM** (optional): Required only for
    `generate-hpc-script`. Generated scripts are self-contained.

## Overview

Generates production-ready GROMACS workflows, validates MDP parameters, and
creates SLURM submission scripts. This skill produces files — it does not
execute simulations directly.

**Do NOT use when:**

-   The user wants to prepare a protein receptor (use `target-preparation`).
-   The user wants to dock compounds (use `compound-screening`).
-   The user wants to generate novel molecules (use `compound-synthesis`).
-   The user wants to *run* a simulation locally — this skill only generates
    setup files. Advise the user to submit to an HPC cluster.

## Core Rules

-   **Use the Wrapper**: ALWAYS execute the provided utility scripts to generate
    MD workflows rather than writing GROMACS commands or MDP files by hand. The
    scripts validate all parameters and apply production-ready defaults.
-   Do not attempt to write MDP files yourself; always use the `gmxapi-setup`
    subcommand which generates validated MDP files automatically.
-   **HPC for production runs**: GROMACS MD is GPU-intensive. Always use
    `generate-hpc-script` rather than advising the user to run locally.
-   If this skill is used, ensure this is mentioned in the output.

## Utility Scripts

All MD simulation setup uses one script with subcommands:

```bash
uv run scripts/md_simulation_api.py <subcommand> --output <file> [options]
```

--------------------------------------------------------------------------------

### 1. Generate gmxapi Workflow

Generates a complete gmxapi Python workflow script and all required MDP files
(energy minimization, NVT equilibration, NPT equilibration, production MD).

```bash
uv run scripts/md_simulation_api.py gmxapi-setup \
    --pdb protein_clean.pdb --ff charmm36 --water tip3p \
    --production-ns 100 --temperature 300 \
    --output /path/to/md_workflow/run_md.py
```

This produces:
-   `run_md.py` — Complete Python workflow using `gmxapi`
-   `ions.mdp`, `em.mdp`, `nvt.mdp`, `npt.mdp`, `md.mdp` — Validated MDP files

> **Note**: The input PDB file can come from any source — a PDB download,
> AlphaFold prediction, or homology model. It does not need to come from
> another skill.

--------------------------------------------------------------------------------

### 2. Validate MDP Files

Validates existing MDP files for correctness (integrator, nsteps, cutoffs,
thermostat/barostat settings).

```bash
uv run scripts/md_simulation_api.py gmxapi-validate \
    --mdp em.mdp nvt.mdp npt.mdp md.mdp \
    --output /path/to/validation.json
```

--------------------------------------------------------------------------------

### 3. Generate MD Setup Shell Script

Generates a GROMACS MD setup shell script as a fallback when gmxapi is not
desired.

```bash
uv run scripts/md_simulation_api.py generate-md-setup \
    --pdb protein_clean.pdb --ff charmm36 --water tip3p \
    --production-ns 100 --output /path/to/setup_md.sh
```

--------------------------------------------------------------------------------

### 4. Generate HPC Submission Script

Generates a SLURM HPC submission script. Supports three job types:

```bash
# GROMACS MD
uv run scripts/md_simulation_api.py generate-hpc-script \
    --job-type gromacs-md --ngpu 4 --mem 128G --time 48:00:00 \
    --production-ns 100 --output /path/to/submit_md.sh

# REINVENT 4
uv run scripts/md_simulation_api.py generate-hpc-script \
    --job-type reinvent --config reinvent.toml \
    --gpu a100 --ngpu 1 --mem 64G --time 24:00:00 \
    --output /path/to/submit_reinvent.sh

# Vina batch screening
uv run scripts/md_simulation_api.py generate-hpc-script \
    --job-type vina-screen --receptor receptor.pdbqt \
    --ligand-dir ligands/ --array 1-100 \
    --output /path/to/submit_screen.sh
```

## Interpreting the Output

1.  **gmxapi Workflow**: The generated `run_md.py` is a self-contained Python
    script. Read it to verify the simulation parameters match the user's intent
    (force field, water model, temperature, production length).
2.  **MDP Validation**: The output JSON reports `valid` (true/false) for each
    MDP file, with specific warnings for invalid parameters.
3.  **HPC Scripts**: Review the generated script for correct `#SBATCH`
    directives (GPU count, memory, walltime). Advise the user to submit with
    `sbatch <script>.sh`.

## Common Options

-   `--output FILE`: **Required.** Output file path.
-   `--pdb FILE`: Input PDB file.
-   `--ff {amber99sb,charmm36,oplsaa}`: Force field (default: charmm36).
-   `--water {spc,tip3p,tip4p}`: Water model (default: tip3p).
-   `--production-ns N`: Production run length in nanoseconds.
-   `--temperature N`: Temperature in Kelvin (default: 300).
-   `--job-type {gromacs-md,reinvent,vina-screen}`: HPC job type.
-   `--ngpu N`: Number of GPUs.
-   `--mem SIZE`: Memory allocation (e.g., 128G).
-   `--time HH:MM:SS`: Walltime.

## Reference

-   **gmxapi**: https://manual.gromacs.org/documentation/current/gmxapi/
-   **GROMACS Manual**: https://manual.gromacs.org/

## Workflow

1.  Provide a clean PDB file (from any source).
2.  Generate a gmxapi workflow with `gmxapi-setup`.
3.  Validate the MDP files with `gmxapi-validate`.
4.  Generate an HPC submission script with `generate-hpc-script`.
5.  Advise the user to submit with `sbatch submit_md.sh`.
