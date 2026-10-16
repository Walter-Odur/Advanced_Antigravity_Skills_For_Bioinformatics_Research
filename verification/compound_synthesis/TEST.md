# TEST.md — compound-synthesis

A verification protocol for the `compound-synthesis` skill. Work through it
to satisfy yourself the skill is accurate and functional on your machine.

**Not shipped with the skill.** This lives under `verification/` so the
bundle under `.agents/skills/` stays lean.

Every expected value below was **measured** by running the command, not
copied from the documentation.

## Setup

```bash
cd <workshop root>
export S=.agents/skills/compound_synthesis/scripts/compound_synthesis_api.py
export W=/tmp/syn_test && mkdir -p $W
```

This skill's launcher declares only `rdkit`, because everything except
`run` and `score-molecules` works without REINVENT installed. **That is the
point of most of the cases below**: config generation, seed preparation,
preflight and result analysis are all useful on a laptop, while the run
itself belongs on a GPU node.

| Need | For |
|---|---|
| `uv` | all cases |
| REINVENT 4 + priors + PyTorch | T10 only; T9 is the case for *not* having it |

--------------------------------------------------------------------------------

## T1 — check-setup tells the truth about a machine

```bash
uv run $S check-setup --output $W/setup.json
```

**Measured on a machine with no REINVENT:**

```
REINVENT 4 setup: NOT READY
  checkout: E:\ANTIGRAVITY_WORKSHOP\REINVENT4
  version: unknown
  importable: False  executable: not on PATH
  priors: 0 found, 10 missing (in E:\ANTIGRAVITY_WORKSHOP\REINVENT4\priors)
  device: cpu
  PROBLEM: The 'reinvent' package is not importable.
  PROBLEM: The 'reinvent' command is not on PATH.
  PROBLEM: 10 of 10 prior model files are missing.
  PROBLEM: PyTorch is not installed.
  -> Install REINVENT 4: git clone https://github.com/MolecularAI/REINVENT4 && pip install -e REINVENT4, then set REINVENT_DIR.
  -> Download the priors from https://doi.org/10.5281/zenodo.15641296 into <checkout>/priors, or set REINVENT_PRIOR_BASE to where they live.
  -> REINVENT needs PyTorch. Install the build matching your CUDA version from https://pytorch.org.
```

**Pass criteria**

1.  Each of the four independent failures is reported **separately**. A
    single "REINVENT not available" would leave a user who has the package
    but not the priors with nothing to act on.
2.  Three distinct remedies, each a command or a DOI — not "install the
    dependencies".
3.  `device: cpu` is reported honestly rather than defaulting to `cuda`.
4.  **Exit code 0.** This is a diagnostic; reporting that a machine is not
    ready is a successful run. A non-zero code here would break any
    pipeline that calls `check-setup` to decide what to do next.

On a ready machine the same command prints `READY`, the detected version,
and the priors it found.

--------------------------------------------------------------------------------

## T2 — the scoring profiles, and their caveats

```bash
uv run $S list-profiles --output $W/profiles.json
```

**Measured:** six profiles — `drug-like` (9 components), `lead-like` (8),
`fragment-like` (6), `kinase-inhibitor` (9), `anti-tb` (7), `cns` (8), all
aggregating by `geometric_mean`.

**Pass criteria**

1.  `geometric_mean`, not arithmetic. With a geometric mean a single
    component scoring zero zeroes the total, so a molecule cannot buy its
    way past an unmet requirement by excelling elsewhere. An arithmetic
    mean would let a beautifully drug-like molecule with a reactive warhead
    score well.
2.  **The profiles that overreach say so.** `anti-tb` carries:

    > A property profile is not an activity model. Whole-cell
    > *M. tuberculosis* potency depends on envelope permeability and efflux,
    > which no physicochemical descriptor captures. Add a target-based
    > component (docking through DockStream, or a QSAR model through
    > ChemProp or Qptuna) before treating output as anti-tubercular leads.

    and `cns` states it does not model active efflux at the blood-brain
    barrier. A profile named after a disease invites being mistaken for an
    activity model; the caveat is what stops that.
3.  Generators and modes are listed: 5 generators (`reinvent`, `libinvent`,
    `linkinvent`, `mol2mol`, `pepinvent`) and 4 run modes (`sampling`,
    `scoring`, `staged_learning`, `transfer_learning`).

--------------------------------------------------------------------------------

## T3 — the generated TOML is valid REINVENT 4

```bash
uv run $S generate-config --mode staged_learning --generator reinvent \
    --scoring-profile anti-tb --num-steps 300 --output $W/reinvent.toml
```

**Measured stdout:**

```
REINVENT config: staged_learning with the reinvent generator
  run_type: staged_learning
  generator_architecture: RNN
  prior_file: priors/reinvent.prior
  device: cuda:0
  max_steps: 300   batch_size: 64   sigma: 128
  diversity_filter: True
  scoring_components: 7
  aggregation: geometric_mean
  expected outputs: staged_learning_1.csv, staged_learning.chkpt, tb_logs/
  warning: The prior model file 'reinvent.prior' was not found locally ...
  warning: device is 'cuda:0' but PyTorch is not installed in this
           interpreter, so the setting could not be checked.
  warning: A property profile is not an activity model. ...
  next: This run will evaluate roughly 19,200 molecules (300 steps x batch 64).
```

**Pass criteria**

1.  **It parses, and matches REINVENT's schema.** The nesting is unusual
    and easy to get subtly wrong:

```bash
python -c "
import tomllib
d=tomllib.load(open('$W/reinvent.toml','rb'))
print('top:', sorted(d.keys()))
print('stage:', sorted(d['stage'][0].keys()))
print('components:', len(d['stage'][0]['scoring']['component']))"
```

**Measured:**

```
top: ['device', 'diversity_filter', 'json_out_config', 'learning_strategy',
      'parameters', 'run_type', 'stage', 'tb_logdir']
stage: ['chkpt_file', 'max_score', 'max_steps', 'min_steps', 'scoring',
        'termination']
components: 7
```

    Each component is a single-key table inside the `component` list, and
    each carries an `endpoint` list — `[[stage.scoring.component.<Name>.endpoint]]`.
    A flat `[scoring.components]` table would be rejected by REINVENT at
    load time, after the job has been queued.

2.  **The arithmetic in the advice is right.** 300 steps × batch 64 =
    19,200 molecule evaluations. That is the number that decides whether a
    run fits in a walltime, so it is stated rather than left to the user.

3.  Expected output filenames are stated up front, so `analyze-results`
    can be pointed at the right file without guessing.

--------------------------------------------------------------------------------

## T4 — transform directions, which silently invert an objective

The subtlest correctness question in the skill. A transform maps a raw
descriptor onto a 0–1 reward. Pick the wrong one and the agent optimises
*away* from what you asked for, while every log looks healthy.

```bash
python -c "
import tomllib
for c in tomllib.load(open('$W/reinvent.toml','rb'))['stage'][0]['scoring']['component']:
    n=list(c)[0]; ep=c[n]['endpoint'][0]
    print('%-18s weight=%-5s %s' % (n, ep.get('weight'),
          (ep.get('transform') or {}).get('type')))"
```

**Measured:**

```
MolecularWeight    weight=1     double_sigmoid
SlogP              weight=1     double_sigmoid
TPSA               weight=1     double_sigmoid
QED                weight=0.5   None
NumRotBond         weight=1     reverse_sigmoid
SAScore            weight=1.5   reverse_sigmoid
custom_alerts      weight=1     None
```

**Pass criteria**

| Component | Required transform | Why |
|---|---|---|
| MolecularWeight, SlogP, TPSA | `double_sigmoid` | A **window** — too low is as bad as too high |
| NumRotBond, SAScore | `reverse_sigmoid` | **Lower is better** |
| QED | none | Already 0–1 with higher better |
| custom_alerts | none | Binary |

**`SAScore` must be `reverse_sigmoid`.** A plain `sigmoid` would reward
*harder*-to-synthesise molecules, and the run would converge happily on
unmakeable chemistry. Its weight of 1.5 — the highest in the profile —
says synthesisability outranks any single property window.

```bash
python -c "
import tomllib, json
for c in tomllib.load(open('$W/reinvent.toml','rb'))['stage'][0]['scoring']['component']:
    if list(c)[0]=='SAScore': print(json.dumps(c['SAScore']['endpoint'][0]))"
# {"name": "Synthetic accessibility", "weight": 1.5,
#  "transform": {"type": "reverse_sigmoid", "high": 4.5, "low": 1, "k": 0.5}}
```

### T4b — the profiles are genuinely different

A profile that merely renames `drug-like` is worse than no profile.

```bash
uv run $S generate-config --mode staged_learning --generator reinvent \
    --scoring-profile drug-like --num-steps 300 --output $W/dl.toml
python -c "
import tomllib
for f,l in [('dl.toml','drug-like'),('reinvent.toml','anti-tb')]:
    for c in tomllib.load(open('$W/'+f,'rb'))['stage'][0]['scoring']['component']:
        if list(c)[0]=='SlogP':
            t=c['SlogP']['endpoint'][0]['transform']
            print('%-10s logP window %s to %s' % (l, t['low'], t['high']))"
```

**Measured:**

```
drug-like  logP window 1 to 4
anti-tb    logP window 2 to 6
```

**Pass criterion:** the `anti-tb` window is shifted upward. The profile's
own rationale cites bedaquiline at 555 Da and clogP above 7, and T2 of the
*compound-screening* protocol measures bedaquiline's logP at 7.13 — so the
window has to accommodate chemistry that `drug-like` would reject. If both
profiles print the same window, the rationale is decoration.

The unwanted-substructure SMARTS are REINVENT's documented default set —
macrocycles `[*;r{8-17}]`, peroxides `[#8][#8]`, carbocations, disulfides,
acetylenes and the rest, 14 patterns — carried as `custom_alerts` params
rather than silently dropped.

--------------------------------------------------------------------------------

## T5 — seeds are de-duplicated chemically, not textually

```bash
uv run $S prepare-seeds \
  --smiles "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC" \
           "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC" \
           "c1ccccc" "CCO" "OCC" \
  --output $W/seeds.smi
```

**Measured:**

```
Seed preparation for reinvent
  input 5  written 2  rejected 3
    dropped compound_2: duplicate of compound_1
    dropped compound_3: invalid SMILES
    dropped compound_5: duplicate of compound_4
```

**Pass criteria**

1.  **`CCO` and `OCC` are recognised as the same molecule.** They are
    different strings for ethanol. A textual de-duplication keeps both, and
    the agent then sees a seed set with a hidden duplicate — which biases
    transfer learning towards it. De-duplication must happen on canonical
    SMILES.
2.  The invalid structure is dropped with its reason, not silently skipped.
3.  Each rejection names the compound it duplicates, so a user can tell a
    genuine duplicate from a misparse.

--------------------------------------------------------------------------------

## T6 — preflight refuses before the GPU bill starts

```bash
uv run $S preflight --config $W/reinvent.toml --output $W/pre.json \
    >/dev/null 2>&1; echo $?
```

**Measured:**

```
error: The config references 2 file(s) that do not exist:
       parameters.prior_file = priors/reinvent.prior,
       parameters.agent_file = priors/reinvent.prior
hint: Run 'check-setup' to see whether the prior models are installed, and
      'prepare-seeds' to create a SMILES file.
rc=6
```

**Pass criteria**

1.  Exit **6** — a missing local file, distinct from a missing dependency
    (4) or bad usage (2).
2.  **Every** missing path is listed with the config key that references it,
    not just the first. Fixing one at a time across a queue wait is the
    expensive failure mode this prevents.
3.  `run` performs the same check and refuses identically, so a job cannot
    reach a scheduler only to die after loading. Measured: `run` on the
    same config also exits 6 with the same message.

--------------------------------------------------------------------------------

## T7 — analysing a real run: 96 rows are not 96 molecules

`tests/data/reinvent_staged_learning_1.csv` is authentic REINVENT 4
staged-learning output, in the real column format:

```
Agent,Prior,Target,Score,SMILES,SMILES_state,QED drug-likeness,
Molecular weight,Lipophilicity (SlogP),Polar surface area
```

```bash
uv run $S analyze-results --csv tests/data/reinvent_staged_learning_1.csv \
    --top-n 5 --output $W/top.json
```

**Measured:**

```
REINVENT results: reinvent_staged_learning_1.csv
  rows 96  unique molecules 13  duplicate samples 83
  score: mean 0.7761 median 0.8831 max 0.9487
  scaffolds: 6 distinct (0.4615)
  diversity: 6 scaffold(s) across 13 molecules (46.2%). Too few molecules to
             judge diversity by ratio - read the scaffold list rather than
             this percentage. ...
  top 5:
      1. score 0.94875 (x6) MW 295.34 QED 0.795 COc1cc2ncnc(Nc3ccc(C)cc3)c2cc1OC
      2. score 0.9485  (x4) MW 281.32 QED 0.794 COc1cc2ncnc(Nc3ccccc3)c2cc1OC
      3. score 0.946725 (x7) MW 344.22 QED 0.787 CN(C)c1cc2ncnc(Nc3cccc(Br)c3)c2cn1
      4. score 0.94215 (x7) MW 330.19 QED 0.769 CNc1cc2ncnc(Nc3cccc(Br)c3)c2cn1
      5. score 0.929075 (x7) MW 315.76 QED 0.787 COc1cc2ncnc(Nc3ccc(Cl)cc3)c2cc1OC
```

**Pass criteria**

1.  **96 rows, 13 unique molecules, 83 duplicate samples.** This is the
    headline number and the one most easily misreported. A report claiming
    "96 compounds generated" overstates the yield by 7×. REINVENT samples
    with replacement, so the row count is a sampling count, not a
    molecule count.
2.  The sample count per molecule is shown — `(x6)`, `(x7)`. A molecule
    sampled seven times is a converged mode, which is information about the
    run, not noise.
3.  Ranking is by score descending, and properties are recomputed from the
    SMILES rather than trusted from the CSV.
4.  The top hits are 4-anilinoquinazolines — a real EGFR chemotype. The
    agent having converged on recognisable chemistry is the sanity check
    that the scoring function is pointing somewhere real.

--------------------------------------------------------------------------------

## T8 — mode collapse must be called mode collapse

The diagnosis a generative-chemistry run most needs, and the one most
easily botched by using a ratio alone.

```bash
python -c "
import csv
rows=[('COc1cc2ncnc(Nc3ccc(C)cc3)c2cc1OC',0.94),('COc1cc2ncnc(Nc3ccccc3)c2cc1OC',0.93),
      ('COc1cc2ncnc(Nc3ccc(Cl)cc3)c2cc1OC',0.92),('COc1cc2ncnc(Nc3ccc(F)cc3)c2cc1OC',0.91),
      ('COc1cc2ncnc(Nc3ccc(Br)cc3)c2cc1OC',0.90)]
with open('$W/collapsed.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['Agent','Prior','Target','Score','SMILES','SMILES_state'])
    for i in range(80):
        s,sc=rows[i%5]; w.writerow([0,0,0,sc,s,1])"

uv run $S analyze-results --csv $W/collapsed.csv --top-n 2 --output $W/col.json
```

Eighty samples, five molecules, all five sharing one Bemis-Murcko
framework — textbook collapse.

**Measured:**

```
  rows 80  unique molecules 5  duplicate samples 75
  scaffolds: 1 distinct (0.2)
  diversity: 1 scaffold(s) across 5 molecules (20.0%). Every molecule shares
             a single Bemis-Murcko framework: the agent has collapsed onto
             one scaffold and is only decorating it. Enable the diversity
             filter, or raise its bucket size if it is already on.
```

**Pass criteria**

1.  The word used is **collapse**, with a remedy.
2.  **The absolute scaffold count decides, not the ratio.** 1 scaffold
    across 5 molecules is a ratio of 0.20 — numerically identical to 20
    scaffolds across 100 molecules, which is a healthy run. This case
    previously reported *"good scaffold diversity"* for a single-scaffold
    run, because only the ratio was consulted. A pure-ratio diagnosis
    cannot distinguish the two, so the count is checked first.
3.  A short run says so instead of guessing. With 13 molecules (T7) the
    diagnosis declines to grade a percentage and points at the scaffold
    list.
4.  **A diverse run is still graded good** — the fix must not simply make
    everything a warning. Build a run spanning many ring systems:

```bash
python -c "
import csv, random
random.seed(0)
cores=['c1ccc2ncnc(N)c2c1','c1ccc2[nH]cnc2c1','c1ccc2oc(N)nc2c1','c1ccc2sc(N)nc2c1',
 'c1ccc2nc(N)nc(N)c2c1','C1CCNCC1','c1ccncc1','c1ccc(-c2ccccc2)cc1','c1cc2ccccc2cn1',
 'C1CC2CCC1C2','c1ccc2[nH]ncc2c1','O=C1NCCN1','c1ncc2ccccc2n1','C1COCCN1','c1ccsc1']
mols=[c.replace('N','N'+s,1) if 'N' in c else c+s
      for c in cores for s in ('C','CC')]
with open('$W/healthy.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['Agent','Prior','Target','Score','SMILES','SMILES_state'])
    for m in mols:
        for _ in range(2):
            w.writerow([0,0,0,round(0.6+0.3*random.random(),4),m,1])"

uv run $S analyze-results --csv $W/healthy.csv --top-n 2 --output $W/h.json
```

**Measured:**

```
  rows 60  unique molecules 30  duplicate samples 30
  scaffolds: 15 distinct (0.5769)
  diversity: 15 scaffold(s) across 26 molecules (57.7%): good scaffold diversity.
```

    26 of the 30 strings parse (the crude substitution above produces four
    invalid ones), 15 distinct frameworks, graded good — as it should be.

--------------------------------------------------------------------------------

## T9 — the run paths fail usefully without REINVENT

```bash
uv run $S score-molecules --smiles "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC" \
    --scoring-profile drug-like --output $W/scored.json >/dev/null 2>&1; echo $?
```

**Measured:**

```
error: Required dependency not available: the 'reinvent' command
hint: Install it with: clone and install REINVENT 4, then re-run
      'check-setup': git clone https://github.com/MolecularAI/REINVENT4 &&
      pip install -e REINVENT4 && export REINVENT_DIR=$PWD/REINVENT4.
      Prior models: https://doi.org/10.5281/zenodo.15641296
rc=4
```

**Pass criteria**

1.  Exit **4**, and only the two subcommands that genuinely need REINVENT
    (`run`, `score-molecules`) return it. Everything in T1–T8 works on a
    machine without it, which is what makes the skill usable before the
    cluster is set up.
2.  The hint is a complete recipe: clone, install, set the environment
    variable, fetch the priors.

Bad arguments are refused at generation time, not at run time:

```bash
uv run $S generate-config --mode staged_learning --generator libinvent \
    --scoring-profile drug-like --num-steps 100 --output $W/bad.toml
```

**Measured:**

```
error: The libinvent generator needs seed structures.
hint: Pass --smiles-file with one scaffold per line with attachment points,
      written as '*'. Use 'prepare-seeds' to validate and de-duplicate them
      first.
rc=2
```

**Pass criterion:** the message names the **attachment-point syntax**.
"libinvent requires seeds" would leave the commonest mistake — scaffolds
written without `*` — undiagnosed.

--------------------------------------------------------------------------------

## T10 — a real run, on a machine that has REINVENT

Not executed here; this machine has no REINVENT, no priors and no PyTorch
(T1). Recorded so the protocol is complete.

```bash
uv run $S check-setup --output $W/setup.json            # must say READY
uv run $S preflight --config $W/reinvent.toml --output $W/pre.json
uv run $S run --config $W/reinvent.toml --output $W/run_report.json
uv run $S analyze-results --csv staged_learning_1.csv --top-n 25 \
    --output $W/top.json
```

**Pass criteria**

1.  `preflight` passes before `run` is attempted.
2.  `run` writes the three artefacts the config named in T3:
    `staged_learning_1.csv`, `staged_learning.chkpt`, `tb_logs/`.
3.  Row count in the CSV is near 19,200 (300 × 64), the figure T3 predicted.
4.  `analyze-results` on real output reproduces T7's distinctions — rows
    against unique molecules, and a diversity verdict consistent with
    whether the diversity filter was on.

--------------------------------------------------------------------------------

## T11 — the SLURM script runs on a cluster

```bash
uv run $S generate-hpc-script --config $W/reinvent.toml --output $W/submit.sh
```

**Measured:** `1xa100, 8 CPU, 64G, 24:00:00`.

**Pass criteria**

1.  `bash -n $W/submit.sh` exits 0.
2.  No CR bytes —
    `python -c "print(open('$W/submit.sh','rb').read().count(b'\r'))"`
    prints `0`. Measured: 0 CR, 49 LF. A CRLF script dies on a Linux
    cluster with `bash\r: bad interpreter` before running a line. Check
    bytes rather than `grep -c $'\r'`, which under Git Bash can match every
    line of a perfectly clean file.
3.  The `#SBATCH` block is complete: `--job-name`, `--output`, `--error`,
    `--time`, `--mem`, `--cpus-per-task`, `--gres=gpu:a100:1`.
4.  `conda activate` goes through `eval "$(conda shell.bash hook)"`, which
    a non-interactive SLURM shell otherwise lacks — a bare
    `conda activate` fails there with "your shell has not been properly
    configured".
5.  It checks `reinvent` is on PATH *inside the job*, before the long step,
    so a mis-set environment fails in seconds rather than after the queue
    wait.

--------------------------------------------------------------------------------

## Results

| Case | What it checks | Result |
|---|---|---|
| T1 | check-setup diagnosis | pass — 4 problems, 3 remedies, rc 0 |
| T2 | profiles and caveats | pass — 6 profiles, geometric mean, anti-tb caveat present |
| T3 | TOML validity and schema | pass — parses, correct nesting, 19,200 arithmetic right |
| T4 | transform directions | pass — SAScore reverse_sigmoid, windows double_sigmoid |
| T4b | profiles genuinely differ | pass — logP 1–4 vs 2–6 |
| T5 | chemical de-duplication | pass — `CCO` and `OCC` collapsed |
| T6 | preflight before the GPU bill | pass — rc 6, both missing paths named |
| T7 | real run analysis | pass — 96 rows read as 13 molecules |
| T8 | mode collapse | pass — **after a fix**; see below |
| T9 | missing REINVENT | pass — rc 4 with a full recipe; rc 2 for bad args |
| T10 | a real run | **not executed** — no REINVENT on this machine |
| T11 | SLURM script | pass — `bash -n` clean, 0 CR bytes, hook present |

**T8 found a real defect.** A single-scaffold run was reported as *"good
scaffold diversity"*, because the diagnosis compared only
`n_scaffolds / n_molecules` against 0.20 and 1/5 lands exactly on the
boundary. `analyze.py` now checks the absolute scaffold count first and
only applies ratio bands once there are at least 20 molecules.

**Environment:** Windows 11 Pro, Python 3.13.13, RDKit 2026.03.6. REINVENT
4, its priors and PyTorch **not installed** — which is what T1, T6, T9 and
T11 verify the skill handles, and why T10 was not run.

## Automated equivalent

```bash
pytest tests/test_reinvent_config.py tests/test_reinvent_analyze.py \
       tests/test_reinvent_setup.py -v
pytest tests/test_hpc.py -v                            # T11
pytest tests/test_reinvent_run.py -m requires_reinvent  # T10, if installed
```

The manual protocol exists because it is readable, needs no test
environment, and shows what a correct result looks like — T4's transform
table and T7's row-versus-molecule distinction in particular are worth
reading before trusting any generative output.
