# TEST.md — md-simulation

A verification protocol for the `md-simulation` skill. Work through it to
satisfy yourself the skill is accurate and functional on your machine.

**Not shipped with the skill.** This lives under `verification/` so the
bundle under `.agents/skills/` stays lean.

Every expected value below was **measured** by running the command, not
copied from the documentation.

## Setup

```bash
cd <workshop root>
export S=.agents/skills/md_simulation/scripts/md_simulation_api.py
export W=/tmp/md_test && mkdir -p $W
```

**This skill needs nothing but Python.** Its PEP 723 header declares no
dependencies at all:

```bash
head -4 $S
# #!/usr/bin/env python3
# # /// script
# # requires-python = ">=3.10"
# # ///
```

Everything below runs on a laptop with no GROMACS, because the skill
*generates* simulation input rather than running it. GROMACS is needed only
to execute what it writes.

| Need | For |
|---|---|
| `uv` (or any Python ≥3.10) | all cases |
| GROMACS | T9 only, which is not run here |

--------------------------------------------------------------------------------

## T1 — a complete CHARMM36 setup at body temperature

```bash
uv run $S gmxapi-setup --pdb $W/6HEZ_clean.pdb --ff charmm36 \
    --temperature 310 --production-ns 100 --output $W/run_md.py
```

**Measured:**

```
MD setup for .../6HEZ_clean.pdb
  force field: charmm36 (charmm36-jul2022), water tip3p
  temperature: 310 K, pressure 1 bar
  production: 100 ns at 0.002 ps/step
  CHARMM36 requires a force-switched LJ potential between 1.0 and 1.2 nm
    and no dispersion correction.
  MDP files: em.mdp, ions.mdp, md.mdp, npt.mdp, nvt.mdp
  shell script: .../setup_md.sh
  next: bash setup_md.sh, then submit the production run
```

**Pass criteria**

1.  Five MDP files, a shell script and a `run_md.py` gmxapi workflow.
2.  The force field's specific requirement is **stated**, not just applied.
3.  `charmm36-jul2022` — the dated release, since `charmm36` alone is
    ambiguous between GROMACS ports.
4.  The step count matches the request exactly:

```bash
grep -E "^(nsteps|dt)" $W/md.mdp
# dt       = 0.002
# nsteps   = 50000000
```

    50,000,000 × 0.002 ps = 100,000 ps = **100 ns**. If `nsteps` does not
    equal `production_ns / dt × 1000`, the run is not the length requested
    and nobody notices until the trajectory is short.

--------------------------------------------------------------------------------

## T2 — non-bonded settings belong to the force field

The central correctness claim of this skill. CHARMM36 and AMBER require
*different, non-interchangeable* van der Waals treatments, and mixing them
changes the physics while running perfectly happily.

```bash
uv run $S gmxapi-setup --pdb $W/6HEZ_clean.pdb --ff amber99sb-ildn \
    --temperature 300 --production-ns 50 --output $W/amber/run_md.py
grep -iE "rvdw|vdw-modifier|DispCorr|rcoulomb" $W/md.mdp $W/amber/md.mdp
```

**Measured:**

| Setting | CHARMM36 | AMBER99SB-ILDN |
|---|---|---|
| `vdw-modifier` | `force-switch` | `potential-shift-verlet` |
| `rvdw_switch` | `1.0` | *(absent)* |
| `rvdw` | `1.2` | `1.0` |
| `rcoulomb` | `1.2` | `1.0` |
| `DispCorr` | *absent — 0 lines* | `EnerPres` |

**Pass criteria**

1.  CHARMM36 gets a **force-switched** potential from 1.0 to 1.2 nm and
    **no** dispersion correction. The CHARMM36 parameters were fitted with
    that switching function; adding a dispersion correction double-counts
    the long-range attraction.
2.  AMBER gets a plain 1.0 nm cutoff **with** `DispCorr = EnerPres`. The
    AMBER parameters were fitted assuming the analytic correction is
    applied.
3.  `rcoulomb` tracks `rvdw` in both. Mismatched cutoffs are a performance
    and accuracy trap.
4.  The announcement differs per force field — AMBER prints *"AMBER force
    fields use a plain 1.0 nm cutoff with an analytic dispersion
    correction."*

If both force fields produce the same non-bonded block, the skill is
applying one force field's settings to the other's parameters.

--------------------------------------------------------------------------------

## T3 — the four-stage protocol is staged correctly

```bash
for f in em nvt npt md; do
  echo "== $f"
  grep -iE "integrator|^define|tcoupl|pcoupl |refcoord|gen_vel|continuation" $W/$f.mdp
done
```

**Measured:**

| | `em` | `nvt` | `npt` | `md` |
|---|---|---|---|---|
| `integrator` | `steep` | `md` | `md` | `md` |
| `define` | — | `-DPOSRES` | `-DPOSRES` | — |
| `tcoupl` | — | `V-rescale` | `V-rescale` | `V-rescale` |
| `pcoupl` | — | `no` | `C-rescale` | `Parrinello-Rahman` |
| `refcoord_scaling` | — | — | `com` | — |
| `gen_vel` | — | `yes` | `no` | `no` |
| `continuation` | — | `no` | `yes` | `yes` |

**Pass criteria** — each row is a distinct way to get this wrong:

1.  **`gen_vel = yes` appears only in NVT**, with `continuation = no`.
    Regenerating velocities in NPT or production discards the
    equilibration that was just paid for.
2.  **`refcoord_scaling = com` is present in `npt.mdp`.** Position
    restraints combined with pressure coupling and *no* `refcoord_scaling`
    is a classic silent error: GROMACS scales the box but not the restraint
    reference coordinates, so restraints pull atoms toward positions that
    drift away from them. `grompp` warns; most people pass `-maxwarn` and
    never read it.
3.  **`C-rescale` for equilibration, `Parrinello-Rahman` for production.**
    C-rescale (Bernetti & Bussi 2020) is stochastic and samples the correct
    ensemble while the box is still far from equilibrium, where
    Parrinello-Rahman can oscillate. Berendsen, the old default, does not
    sample the correct ensemble at all.
4.  **Position restraints are released for production.** `md.mdp` must have
    no `define = -DPOSRES`, or the protein is frozen through the run that
    is supposed to measure its motion.
5.  `em` uses `steep` with `emtol = 1000.0` kJ/mol/nm.

--------------------------------------------------------------------------------

## T4 — the validator catches a mismatched force field

The case that makes T2 enforceable rather than aspirational. Validate the
AMBER production file *against CHARMM36*:

```bash
uv run $S gmxapi-validate --mdp $W/amber/md.mdp --ff charmm36 \
    --output $W/val_bad.json
```

**Measured:**

```
MDP validation: 1 file(s), PROBLEMS FOUND
  md.mdp: INVALID (32 parameters)
      ERROR: charmm36 requires vdw-modifier = force-switch, but the file has
             potential-shift-verlet. CHARMM36 requires a force-switched LJ
             potential between 1.0 and 1.2 nm and no dispersion correction.
      ERROR: charmm36 must not use a dispersion correction, but
             DispCorr = EnerPres is set. ...
      ERROR: charmm36 expects rvdw = 1.2 nm, but the file has 1.0 nm.
      ERROR: charmm36 expects rcoulomb = 1.2 nm, but the file has 1.0 nm.
      Simulated time: 50 ns (2.5e+07 steps x 0.002 ps).
```

**Pass criteria**

1.  **All four mismatches** are reported, not just the first.
2.  Each says what the file has, what the force field needs, and **why**.
3.  The matched file validates clean:
    `gmxapi-validate --mdp $W/md.mdp --ff charmm36` →
    `md.mdp: valid (32 parameters)`, with
    `Simulated time: 100 ns (5e+07 steps x 0.002 ps)` — the simulated time
    recomputed from the file rather than taken on trust.

### T4b — the restraint trap

```bash
grep -v -i refcoord $W/npt.mdp > $W/npt_norefcoord.mdp
uv run $S gmxapi-validate --mdp $W/npt_norefcoord.mdp --output $W/v2.json
```

**Measured:**

```
  npt_norefcoord.mdp: INVALID (33 parameters)
      ERROR: Position restraints (-DPOSRES) are combined with pressure
             coupling but refcoord_scaling is not set. grompp warns, and
             restraint reference coordinates are not scaled with the box.
             Set refcoord_scaling = com.
```

**Pass criterion:** detected from the *combination* of two settings, neither
of which is wrong alone. The intact `npt.mdp` validates clean (34
parameters) — one more parameter, the one that was removed.

### T4c — warnings are not errors

```bash
sed 's/^dt .*= 0.002/dt                   = 0.005/' $W/md.mdp > $W/md_bigdt.mdp
uv run $S gmxapi-validate --mdp $W/md_bigdt.mdp --output $W/v4.json

grep -v -i "^tcoupl" $W/md.mdp > $W/md_notc.mdp
uv run $S gmxapi-validate --mdp $W/md_notc.mdp --output $W/v5.json
```

**Measured — both `valid`, each with a warning:**

```
  md_bigdt.mdp: valid (32 parameters)
      warning: dt = 0.005 ps needs hydrogen mass repartitioning or virtual
               sites to be safe.
      Simulated time: 250 ns (5e+07 steps x 0.005 ps).

  md_notc.mdp: valid (31 parameters)
      warning: No thermostat is set, so this run samples the microcanonical
               ensemble and its temperature will drift.
```

**Pass criterion:** these are **warnings, not errors**. A 5 fs timestep is
legitimate with hydrogen mass repartitioning, and an NVE run is a valid
choice. Only settings that are inconsistent with the declared force field
are errors. A validator that rejects every unusual choice gets switched off.

Note also that the 5 fs file reports 250 ns: changing `dt` without changing
`nsteps` silently changed the run length by 2.5×, and the validator says so.

--------------------------------------------------------------------------------

## T5 — the setup script is cluster-safe

```bash
bash -n $W/setup_md.sh && echo "syntax OK"
python -c "b=open('$W/setup_md.sh','rb').read(); print('CR', b.count(b'\r'), 'LF', b.count(b'\n'))"
grep -nE "set -eu|command -v gmx|pdb2gmx|solvate|genion|grompp" $W/setup_md.sh
```

**Measured:** syntax OK, **CR 0, LF 79**.

```
10:set -euo pipefail
19:if ! command -v gmx >/dev/null 2>&1 && ! command -v gmx_mpi >/dev/null 2>&1; then
23:GMX=$(command -v gmx || command -v gmx_mpi)
37:echo "== 1/8 Topology (pdb2gmx) =="
38:"$GMX" pdb2gmx -f "$PDB" -o processed.gro -p topol.top \
46:"$GMX" solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
49:"$GMX" grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2
51:echo SOL | "$GMX" genion -s ions.tpr -o solvated_ions.gro -p topol.top \
55:"$GMX" grompp -f em.mdp -c solvated_ions.gro -p topol.top -o em.tpr
```

**Pass criteria**

1.  `set -euo pipefail` — an eight-step pipeline must not continue past a
    failed step with half-written files.
2.  **Zero CR bytes.** A CRLF script fails on Linux with
    `bash\r: bad interpreter` before running a line. Check bytes, not
    `grep -c $'\r'`: under Git Bash that can match every line of a clean
    file (it reported 49 of 49 on a file with 0 CR bytes).
3.  It accepts **either** `gmx` or `gmx_mpi`, which is what clusters
    actually install, and fails with a clear message if neither is present.
4.  Every invocation is `"$GMX"`, quoted — module paths contain spaces.
5.  Steps are numbered `1/8`, so a failure says where it stopped.
6.  `genion` is fed its group on stdin (`echo SOL |`) rather than waiting
    for an interactive answer — this is the step that hangs a batch job
    forever when done wrong.

`run_md.py` must also be importable Python:

```bash
python -m py_compile $W/run_md.py && echo "compiles"
```

**Measured:** compiles; CR 0, LF 132.

--------------------------------------------------------------------------------

## T6 — the production SLURM script survives the walltime

```bash
uv run $S generate-hpc-script --job-type gromacs-md --tpr md.tpr \
    --output $W/submit.sh
sed -n '33,52p' $W/submit.sh
```

**Measured:**

```
MAXH=23.52
RESUME=""
if [ -f "${DEFFNM}.cpt" ]; then
  echo "Found ${DEFFNM}.cpt - resuming"
  RESUME="-cpi ${DEFFNM}.cpt"
fi

"$GMX" mdrun \
  -s "$TPR" \
  -deffnm "$DEFFNM" \
  -nb gpu -pme gpu -bonded gpu \
  -ntomp 8 \
  -maxh "$MAXH" \
  $RESUME
```

**Pass criteria**

1.  **`-maxh` is below the walltime and scales with it.** GROMACS stops
    cleanly and writes a final checkpoint when it hits `-maxh`; if that
    value exceeds the walltime, SLURM kills the job mid-write and the run
    is unresumable. Measured across walltimes:

| `--time` | `MAXH` | Margin |
|---|---:|---|
| `00:30:00` | 0.48 | 1 min (the floor) |
| `04:00:00` | 3.92 | ~5 min |
| `24:00:00` | 23.52 | ~29 min |
| `48:00:00` | 47.5 | 30 min (the cap) |
| `2-00:00:00` | 47.5 | 30 min, `D-HH:MM:SS` parsed |

    **This was a real defect.** `MAXH` was hardcoded to `23.5` for every
    walltime, while the comment beside it claimed it matched. At
    `--time 04:00:00` GROMACS was told it could run 23.5 hours and was
    killed at 4 without a checkpoint; at `--time 48:00:00` it stopped at
    23.5 and wasted half the allocation. It is now derived from the
    requested walltime.

2.  **`-cpi` resume is conditional on the checkpoint existing.** A 100 ns
    run does not fit in 24 hours, so it has to be resubmittable; an
    unconditional `-cpi` fails on the first submission.
3.  `-nb gpu -pme gpu -bonded gpu` appears only when a GPU was requested.
4.  `bash -n` clean, 0 CR bytes, and the `#SBATCH` block carries
    `--job-name`, `--output`, `--error`, `--time`, `--mem`,
    `--cpus-per-task`, `--gres=gpu:a100:1`.
5.  A bad walltime is refused at generation time: `--time 24h` exits **2**
    with `Use HH:MM:SS (for example 24:00:00) or D-HH:MM:SS.`

--------------------------------------------------------------------------------

## T7 — a CPU-only job must not request a GPU

```bash
uv run $S generate-hpc-script --job-type vina-screen --array 1-100 \
    --output $W/vs.sh
head -9 $W/vs.sh
grep -c gres $W/vs.sh
```

**Measured:**

```
SLURM script: vina-screen - AutoDock Vina batch docking as a SLURM job array (CPU).
  CPU only, 8 CPU, 64G, 24:00:00, array 1-100
#!/bin/bash
#SBATCH --job-name=vina-screen
#SBATCH --output=vina-screen-%A_%a.out
#SBATCH --error=vina-screen-%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-100
```

`grep -c gres` → **0**.

**Pass criteria**

1.  **No `--gres=gpu` line.** AutoDock Vina is CPU-only. **This was a real
    defect**: the script requested `--gres=gpu:a100:1`, and with
    `--array=1-100` that reserved *one hundred A100s* for a workload that
    cannot use one — queueing behind GPU availability that would never
    help, and burning an allocation if it ever started.
2.  **The summary agrees with the file.** It now reads `CPU only` rather
    than `1xa100`; a console line contradicting the script beside it is its
    own defect.
3.  `gromacs-md` is unaffected and still gets `--gres=gpu:a100:1`.
4.  Array output uses `%A_%a` (job and task), not `%j` — with `%j` all 100
    tasks write to one filename.
5.  `vina-screen` without `--array` is refused: exit **2**, *"A vina-screen
    job needs --array to partition the ligand set"*, with
    `--array 1-100` as the example.

All five job types generate `bash -n`-clean scripts: `gromacs-md`,
`reinvent`, `vina-screen`, `colabfold`, `alphafold2`.

--------------------------------------------------------------------------------

## T7b - each script comes with instructions for running it

A script alone does not say what the cluster needs, how to install it, what
to stage, or what to type. Every generation writes a companion guide.

```bash
uv run $S generate-hpc-script --job-type gromacs-md --tpr md.tpr     --module gromacs/2024.2 --output $W/submit_md.sh
ls $W/submit_md*
```

**Measured:** `submit_md.sh` (1,895 B) and `submit_md_SUBMIT.md` (142
lines). All five job types produce the pair:

| Job type | Script | Guide |
|---|---:|---:|
| `gromacs-md` | 1,895 B | 142 lines |
| `reinvent` | 1,441 B | 189 lines |
| `vina-screen` | 2,671 B | 154 lines |
| `colabfold` | 1,468 B | 138 lines |
| `alphafold2` | 1,528 B | 139 lines |

**Pass criteria**

1.  **One script, one guide, per call.** Every `.sh` has a matching
    `_SUBMIT.md`; `bash -n` clean and 0 CR bytes for all five.
2.  The guide answers four questions under numbered headings: what the job
    asks SLURM for, which tools the cluster must have, which files to stage,
    and how to submit. Plus a failure table.
3.  **The install commands are real.** Verified live: `gromacs` (2026.3),
    `vina` (1.2.7) and `openbabel` (3.2.1) all exist on conda-forge, and
    every URL returns HTTP 200 - the Vina releases page, the REINVENT
    repository, the Zenodo priors DOI, localcolabfold, deepmind/alphafold,
    pytorch.org and the Miniforge installer. The two specific filenames are
    real too: `install_colabbatch_linux.sh` in localcolabfold and
    `run_alphafold.py` in deepmind/alphafold.
4.  **The guide reflects the job, not a template.** `vina-screen` reports
    `GPUs | none - this job is CPU-only`, never mentions `a100`, and adds
    array-specific advice (`squeue -j <jobid> -r`, the `MaxArray` limit,
    `--array 1-100%20` to throttle). `reinvent` instead names the conda
    environment to create and the Zenodo priors.
5.  **`--module` works and the script says when it is absent.** With
    `--module gromacs/2024.2 --module cuda/12.1` the script contains both
    `module load` lines; without any, it states *"No module is loaded, so
    the tools must already be on PATH inside the job"* - the commonest
    reason a job cannot find its tools, since loading a module on a login
    node does not carry into the job.
6.  **The script names its guide.** The GROMACS tool check prints
    `See submit_md_SUBMIT.md section 2 for how to install it.`, so a job log
    read an hour later says where the answer is.
7.  The guide's quoted resources match the directives: walltime, memory,
    CPUs, GPUs and array all read from the same `SlurmResources`.

### T7c - the wizard pairs them too

```bash
uv run $WIZ run --pdb-id 1M17 --chain A --chembl-target CHEMBL203     --stages target sourcing triage docking generative selection dynamics     --workdir $W/pipeline --output $W/pipeline.json
find $W/pipeline -name "submit_*"
```

**Measured**, 14.1 s with no Vina, REINVENT or GROMACS installed:

```
04_docking/submit_vina.sh        04_docking/submit_vina_SUBMIT.md
05_generative/submit_reinvent.sh 05_generative/submit_reinvent_SUBMIT.md
07_dynamics/submit_md.sh         07_dynamics/submit_md_SUBMIT.md
```

**Pass criterion:** all three cluster stages produce the pair, each `bash -n`
clean with 0 CR bytes. `07_dynamics/setup_md.sh` has no guide, correctly -
it runs locally before submission, not on the cluster.

--------------------------------------------------------------------------------

## T8 — generate-md-setup, the shell-only path

```bash
uv run $S generate-md-setup --pdb $W/6HEZ_clean.pdb --ff charmm36 \
    --temperature 310 --production-ns 100 --output $W/setup_only.sh
```

**Measured:** `MDP files written: 5`, `bash -n` clean, **CR 0, LF 79** —
byte-for-byte the same line count as T1's `setup_md.sh`.

**Pass criteria**

1.  It produces the same five MDP files as T1, so a user who does not want
    the gmxapi Python workflow is not given a lesser protocol.
2.  `bash -n` clean, 0 CR bytes.
3.  **The physics is identical.** Two code paths emitting different
    settings for the same request would be the worst failure this skill
    could have:

```bash
python -c "
import re
keys={'rvdw','rvdw_switch','vdw-modifier','rcoulomb','dispcorr','nsteps','dt','ref_t','pcoupl'}
def nb(p):
    return {m.group(1).lower(): m.group(2)
            for line in open(p)
            if (m:=re.match(r'\s*([\w-]+)\s*=\s*(.+?)\s*$', line))
            and m.group(1).lower() in keys}
print('identical:', nb('$W/md.mdp') == nb('$W/shellonly/md.mdp'))"
# identical: True
```

    Measured identical across `dt`, `nsteps`, `ref_t`, `pcoupl`,
    `rcoulomb`, `rvdw`, `rvdw_switch` and `vdw-modifier`, with `DispCorr`
    absent from both.

--------------------------------------------------------------------------------

## T9 — running it, on a machine with GROMACS

Not executed here; this machine has no GROMACS. Recorded so the protocol is
complete.

```bash
cd $W && bash setup_md.sh            # topology, box, solvate, ions, EM, NVT, NPT
sbatch submit.sh                     # production
```

**Pass criteria**

1.  `setup_md.sh` completes all eight steps, ending with `npt.gro` and a
    `topol.top`.
2.  `grompp` on `md.mdp` emits **no warning about `refcoord_scaling`** —
    T4b's check, confirmed by GROMACS itself.
3.  The NPT stage's average density is near 1000 kg/m³ for a solvated
    protein. A density far from that means the box or the ion count is
    wrong, and no MDP check can detect it.
4.  Production writes `md.xtc`, `md.edr`, `md.gro`; resubmitting the same
    script resumes from `md.cpt` rather than restarting.

--------------------------------------------------------------------------------

## Results

| Case | What it checks | Result |
|---|---|---|
| T1 | complete CHARMM36 setup | pass — 5 MDPs, nsteps gives exactly 100 ns |
| T2 | force-field-specific non-bonded | pass — force-switch/no-DispCorr vs cutoff/EnerPres |
| T3 | four-stage protocol | pass — gen_vel, refcoord_scaling, C-rescale all correct |
| T4 | validator on a mismatch | pass — all 4 errors, each with its reason |
| T4b | the restraint trap | pass — caught from the combination |
| T4c | warnings are not errors | pass — 5 fs and NVE warn, stay valid |
| T5 | setup script cluster-safe | pass — 0 CR, `set -euo pipefail`, stdin-fed genion |
| T6 | `-maxh` and resume | pass — **after a fix**; see below |
| T7 | CPU job requests no GPU | pass — **after a fix**; see below |
| T7b | script + instructions per call | pass - 5/5 pairs, install routes verified live |
| T7c | wizard pairs them too | pass - 3/3 cluster stages |
| T8 | the shell-only path | pass — same physics as T1 |
| T9 | an actual run | **not executed** — no GROMACS on this machine |

**Two real defects found here.**

-   **T6:** `MAXH` was hardcoded at `23.5` for every walltime, while the
    generated comment claimed it matched the walltime. Now derived — a 2%
    margin, floored at 1 minute and capped at 30.
-   **T7:** the CPU-only `vina-screen` job type requested
    `--gres=gpu:a100:1`, which under `--array 1-100` would reserve 100
    GPUs. CPU-only job types now have the GPU request dropped, and the
    console summary reports what the file actually asks for.

**Environment:** Windows 11 Pro, Python 3.13.13. GROMACS **not installed**,
which is what T1–T8 verify the skill does not need. The input structure is
the cleaned chain A of PDB 6HEZ produced by the *target-preparation* skill.

## Automated equivalent

```bash
pytest tests/test_md_mdp.py tests/test_md_workflow.py -v
pytest tests/test_hpc.py -v                              # T6, T7
pytest --doctest-modules src/agskills/hpc.py             # the -maxh arithmetic
pytest tests/test_md_run.py -m requires_gromacs          # T9, if installed
```

The manual protocol exists because it is readable, needs no test
environment, and shows what a correct result looks like. T2's table and
T3's staging table are worth reading before trusting any MDP file —
including one you wrote yourself.
