---
name: md-simulation
description: >
  Set up molecular dynamics simulations. Generates GROMACS workflows and MDP
  files with force-field-aware non-bonded settings, validates existing MDP
  files including a cross-check against a declared force field, and produces
  SLURM scripts for GROMACS MD, REINVENT, AutoDock Vina batch screening,
  ColabFold and AlphaFold2. Use when the user wants an MD setup or a cluster
  submission script. Do not use for receptor prep or compound screening.
---

# MD Simulation

This skill **produces files; it does not run simulations.** Production MD is
GPU work that belongs in a batch job.

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions. This
    skill declares no dependencies — every subcommand runs on the standard
    library alone.
2.  **User Notification**: If `.licenses/md_simulation_LICENSE.txt` does not
    exist in the workspace root, prominently notify the user to review the
    terms at https://www.gromacs.org/ and
    http://mackerell.umaryland.edu/charmm_ff.shtml, then create that file
    recording the notification text and a timestamp.
3.  **GROMACS**: needed to *run* the generated workflow, not to generate it.

## Core Rules

-   [!IMPORTANT] **Never hand-write an MDP file or a GROMACS command.** Use
    `gmxapi-setup`, which derives the non-bonded settings from the force
    field and validates everything it writes.

-   **Output to File (Required)**: `--output` is required for every
    subcommand.

-   **Pass `--ff` to `gmxapi-validate`.** An MDP file does not record which
    force field it was written for, so the most consequential class of error
    is invisible without being told.

-   **Production runs go to the cluster.** Use `generate-hpc-script`.

-   **State the force field and temperature** when reporting any result.

-   **Notification**: If this skill is used, mention it in the output.

## Utility Script

```bash
uv run scripts/md_simulation_api.py <subcommand> --output <file> [options]
```

| Subcommand | Purpose |
|---|---|
| `gmxapi-setup` | Generate the workflow, shell script and all MDP files |
| `generate-md-setup` | Generate the setup shell script and MDP files only |
| `gmxapi-validate` | Validate MDP files, optionally against a force field |
| `generate-hpc-script` | SLURM script for five job types |

--------------------------------------------------------------------------------

### 1. Generate a Workflow

```bash
uv run scripts/md_simulation_api.py gmxapi-setup \
    --pdb protein_clean.pdb --ff charmm36 --water tip3p \
    --temperature 310 --production-ns 100 \
    --output md_workflow/run_md.py
```

Writes into `md_workflow/`: `run_md.py` (a gmxapi driver), `setup_md.sh`
(the equivalent shell script, which is what most clusters run),
`ions.mdp`, `em.mdp`, `nvt.mdp`, `npt.mdp`, `md.mdp`, and
`run_md.meta.json`.

Then:

```bash
bash md_workflow/setup_md.sh      # topology, box, solvate, EM, NVT, NPT
uv run scripts/md_simulation_api.py generate-hpc-script \
    --job-type gromacs-md --ngpu 4 --mem 128G --time 48:00:00 \
    --module gromacs/2024.2 \
    --output md_workflow/submit_md.sh
sbatch md_workflow/submit_md.sh   # production
```

Each call writes **one script plus one guide** — `submit_md.sh` and
`submit_md_SUBMIT.md`. Point the user at the guide rather than explaining
the cluster yourself: it lists the required tools with probe and install
commands, the files to stage, the `sbatch`/`squeue`/`sacct` commands, and
what each failure message means. Pass `--module NAME` (repeatable) so the
job loads the cluster's own build; without it the script says so, because a
module loaded on a login node does not carry into the job.

> [!IMPORTANT] **Non-bonded settings belong to the force field.** Each was
> parameterised with a particular treatment of the Lennard-Jones tail; it is
> part of the model, not a preference. CHARMM36 needs a force-switched
> potential from 1.0 to 1.2 nm with **no** dispersion correction; the AMBER
> family needs a plain 1.0 nm cutoff **with** one. Mixing them
> double-counts a dispersion tail and shifts pressure and density. The skill
> derives them from `--ff`.

--------------------------------------------------------------------------------

### 2. Validate MDP Files

```bash
uv run scripts/md_simulation_api.py gmxapi-validate \
    --mdp em.mdp nvt.mdp npt.mdp md.mdp --ff charmm36 \
    --output validation.json
```

Reports `valid`, `errors`, `warnings` and `notes` per file. The notes
include the simulated time implied by `nsteps` and `dt`, which is a useful
sanity check on its own. Checks include group-count consistency
(`ref_t` per `tc-grps`, the commonest grompp abort), barostat completeness,
`refcoord_scaling` with position restraints, timestep stability, and the
force-field cross-check. Full list:
[references/mdp_validation.md](references/mdp_validation.md).

--------------------------------------------------------------------------------

### 3. Generate a Cluster Script

```bash
# GROMACS production
... generate-hpc-script --job-type gromacs-md --ngpu 4 --mem 128G \
    --time 48:00:00 --output submit_md.sh

# REINVENT
... generate-hpc-script --job-type reinvent --config reinvent.toml \
    --conda-env reinvent4 --output submit_reinvent.sh

# Vina batch screening, as a job array
... generate-hpc-script --job-type vina-screen --receptor receptor.pdbqt \
    --ligand-dir ligands/ --site-config binding_site.json \
    --array 1-100 --ngpu 0 --ncpus 16 --output submit_screen.sh

# Structure prediction
... generate-hpc-script --job-type colabfold --fasta target.fasta \
    --output submit_fold.sh
```

Generated scripts are defensive, because a silent failure at hour three of a
24-hour allocation is expensive: `set -euo pipefail`; inputs and executables
checked **before** the long step; `conda activate` through the shell hook a
non-interactive SLURM shell otherwise lacks; GROMACS jobs use `-maxh` so a
walltime limit writes a final checkpoint and `-cpi` resumes from it; array
jobs name output per task and partition the ligand set by a disjoint stride,
skipping work already done. `--time`, `--mem` and array specifications are
validated at generation time.

Review the `#SBATCH` directives against your cluster's partitions and GPU
names before submitting.

## Common Options

-   `--output FILE`: **Required.** Output path.
-   `--pdb FILE`: cleaned input structure.
-   `--ff charmm36|amber99sb|amber99sb-ildn|amber14sb|oplsaa` (default:
    `charmm36`). Determines the non-bonded settings.
-   `--water tip3p|tip4p|tip4pew|spc|spce|none` (default: `tip3p`).
-   `--temperature K` (default: 300; 310 is physiological) — reaches `ref_t`
    and `gen_temp` in the files grompp reads.
-   `--pressure BAR`, `--production-ns N`, `--equilibration-ps N`, `--dt PS`.
-   `--box-dist NM`, `--box-type dodecahedron|cubic|octahedron|triclinic`.
-   `--ion-concentration M` (default: 0.15), `--has-ligand`, `--seed N`.
-   `--job-type`, `--ngpu`, `--ncpus`, `--mem`, `--time`, `--array`.
-   `--detail compact|standard|full`: `full` echoes every parsed MDP
    parameter back; the default reports only errors, warnings and notes.

## Common Mistakes

1.  **Omitting `--ff` from `gmxapi-validate`.** The file cannot declare its
    own force field, so the cross-check is skipped silently.
2.  **Assuming `--temperature` is cosmetic.** It reaches `ref_t`; verify it
    in `md.mdp` if the simulation temperature matters.
3.  **Forgetting the ligand topology.** GROMACS cannot generate one. Use
    CGenFF for CHARMM36 or ACPYPE/GAFF for AMBER, and include the `.itp` in
    `topol.top` before grompp. `--has-ligand` says so in its output.
4.  **Running production locally.** Submit it.

## Reference

-   [references/questions.md](references/questions.md) — **read this first when the user's request is vague, phrased as a question, or mentions a target, organism or compound class not covered above.** Indexed by question: force-field choice, the four stages, walltime and checkpointing, and what to check afterwards.
-   [references/mdp_validation.md](references/mdp_validation.md) — every
    check, the protocol stage by stage, and the force-field table.
-   [references/citation.bib](references/citation.bib) — cite these for any
    simulation.

## Workflow

1.  Obtain a cleaned PDB, from `target-preparation` or elsewhere.
2.  `gmxapi-setup` with the force field, temperature and length you want.
3.  `gmxapi-validate --ff <the same force field>`; read the warnings.
4.  `bash setup_md.sh` to build and equilibrate.
5.  `generate-hpc-script --job-type gromacs-md`, then `sbatch`.
6.  Report force field, water model, temperature and production length
    alongside any result.
