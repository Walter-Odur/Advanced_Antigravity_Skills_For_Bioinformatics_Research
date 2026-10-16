# TEST.md — drug-discovery-wizard

A verification protocol for the `drug-discovery-wizard` skill. Work through
it to satisfy yourself the skill is accurate and functional on your machine.

**Not shipped with the skill.** This lives under `verification/` so the
bundle under `.agents/skills/` stays lean.

Every expected value below was **measured** by running the command, not
copied from the documentation.

This skill orchestrates the other four, so its correctness is mostly about
three things: it must **delegate rather than reimplement**, it must be
**honest about what it cannot run**, and the artefacts it leaves behind must
be **usable somewhere else**.

## Setup

```bash
cd <workshop root>
export S=.agents/skills/drug_discovery_wizard/scripts/drug_discovery_api.py
export W=/tmp/wz_test && mkdir -p $W
```

| Need | For |
|---|---|
| `uv` | all cases |
| Internet | T3 onwards (RCSB, ChEMBL) |
| Vina / REINVENT / GROMACS | none of them — T4 is the case for *not* having them |

--------------------------------------------------------------------------------

## T1 — plan reports readiness per stage, and changes nothing

```bash
uv run $S plan --pdb-id 6HEZ --chain A --chembl-target CHEMBL3804751 \
    --output $W/plan.json
```

**Measured:**

```
Pipeline plan:
  1. target      [local  ] [ready]  target-preparation
       produces: cleaned PDB, receptor PDBQT, binding-site JSON, quality report
  2. sourcing    [local  ] [ready]  compound-screening
       produces: known actives with potency, natural-product hits
  3. triage      [local  ] [ready]  compound-screening
       produces: relaxed ADMET report (14 checks) and a seed set
  4. docking     [cluster] [NOT READY]  compound-screening
       AutoDock Vina is not installed; this stage will produce a cluster
       script instead.
  5. generative  [cluster] [NOT READY]  compound-synthesis
       4 problem(s) must be resolved before REINVENT can run. Config
       generation and seed preparation work regardless.
  6. selection   [local  ] [ready]  compound-screening
       produces: strict ADMET report (24 checks)
  7. dynamics    [cluster] [NOT READY]  md-simulation
       GROMACS is not installed; setup files and a cluster script are still
       produced.
```

**Pass criteria**

1.  Readiness is **per stage**, not one verdict for the pipeline. Four
    stages are ready here and three are not, and a user who only wants
    triage should not be told the pipeline is unavailable.
2.  Each unready stage says **what is missing** and **what it will do
    anyway**. "Config generation and seed preparation work regardless" is
    the difference between a blocked user and a productive one.
3.  **`plan` is read-only.** This is the most important property, and the
    easiest to lose:

```bash
ls $W | wc -l          # 1  -- only plan.json
time uv run $S plan --pdb-id 6HEZ --chain A --output $W/plan2.json
# real 0m1.969s
```

    **Measured: one file written, 1.97 s.** No structure is fetched, no
    package is installed, nothing is cloned. An earlier version of this
    codebase ran `sudo apt-get install` and `git clone`d REINVENT (pulling
    PyTorch) from inside `plan` — a read-only query that could rewrite the
    machine. If `plan` ever takes minutes or touches the network, that
    regression is back.

--------------------------------------------------------------------------------

## T2 — describe documents the stage graph

```bash
uv run $S describe --output $W/desc.json
```

**Measured:** all seven stages, each with its sub-skill, its locality, what
it **needs** and what it **produces**:

```
  target      target-preparation   local
              needs: a PDB id, a UniProt accession, or a local PDB file
  sourcing    compound-screening   local
              needs: a ChEMBL target id or a search term
  triage      compound-screening   local
              needs: compounds from sourcing
  docking     compound-screening   cluster
              needs: receptor PDBQT, binding site, AutoDock Vina
  generative  compound-synthesis   cluster
              needs: seed structures, REINVENT 4, a GPU
  selection   compound-screening   local
              needs: generated or docked compounds
  dynamics    md-simulation        cluster
```

**Pass criteria**

1.  The `needs` of each stage names the **output of an earlier stage**
    ("compounds from sourcing", "receptor PDBQT"), so the dependency graph
    is readable without running anything.
2.  The check counts match the screening skill's registry — 14 relaxed,
    24 strict. If the wizard states its own numbers they will drift from
    `describe-tiers`.

--------------------------------------------------------------------------------

## T3 — the local pipeline, end to end

```bash
uv run $S run --pdb-id 6HEZ --chain A --chembl-target CHEMBL3804751 \
    --stages target sourcing triage --workdir $W/pipeline \
    --output $W/pipeline.json
```

**Measured:**

```
Pipeline in .../pipeline
  target: 6HEZ.pdb -> 3231 atoms, 4 heterogen type(s) removed
          verdict NOT_APPLICABLE; 4 candidate ligand(s)
          site centre (13.869, -20.617, 37.131)
  sourcing: ChEMBL CHEMBL3804751 -> 25 actives
  triage: Tier 1 (relaxed) (14 checks) -> 25/25 passed
  completed in 16.7s
```

```
01_target/6HEZ.pdb              02_sourcing/chembl_actives.json
01_target/6HEZ_clean.pdb        03_triage/relaxed_admet.json
01_target/6HEZ_clean.pdbqt      03_triage/seeds.smi
01_target/binding_site.json
```

**Pass criteria**

1.  **Numbered stage directories.** `01_target`, `02_sourcing`,
    `03_triage` — the execution order is legible from the filesystem months
    later.
2.  **25/25 real actives survive triage.** These are measured DprE1
    inhibitors from ChEMBL; a relaxed tier that rejected them would be
    useless as a seed-selection step. Compare the *compound-screening*
    protocol's T5, where the same property holds for EGFR.
3.  `4 heterogen type(s) removed` — "type(s)", not a molecule count. The
    standalone skill prints `heterogens removed: 4`, which is the same
    number read less precisely; the wizard's phrasing is the clearer one.
4.  16.7 s for a fetch, a clean, a PDBQT conversion, a ChEMBL query and a
    14-check screen of 25 compounds.

### T3b — the wizard delegates; it does not reimplement

The single most valuable check in this file. Compare the wizard's binding
site against the one the *target-preparation* skill produces alone:

```bash
python -c "
import json
b=json.load(open('$W/pipeline/01_target/binding_site.json'))
s=json.load(open('/tmp/tp_test/site.json'))   # from target-preparation T3
print('wizard    :', tuple(round(b['center_'+k],4) for k in 'xyz'))
print('standalone:', tuple(round(s['center_'+k],4) for k in 'xyz'))
print('identical centre:', all(abs(b['center_'+k]-s['center_'+k])<1e-12 for k in 'xyz'),
      '| box:', all(abs(b['size_'+k]-s['size_'+k])<1e-12 for k in 'xyz'),
      '| ref:', b['reference']==s['reference'])"
```

**Measured:**

```
wizard    : (13.869, -20.617, 37.131)
standalone: (13.869, -20.617, 37.131)
identical centre: True | box: True | ref: True
```

Same centre to 1e-12, same box dimensions, same reference ligand
(`0SK:A:501`). **Pass criterion: exact equality, not approximate.** Any
difference means the wizard has its own copy of the site logic, and the two
will diverge the first time one is fixed.

--------------------------------------------------------------------------------

## T4 — cluster stages produce scripts instead of pretending

Run every stage on a machine with **no Vina, no REINVENT, no GROMACS**.

```bash
uv run $S run --pdb-id 6HEZ --chain A --chembl-target CHEMBL3804751 \
    --stages target sourcing triage docking generative selection dynamics \
    --workdir $W/full --output $W/full.json
```

**Measured:**

```
  target: 6HEZ.pdb -> 3231 atoms, 4 heterogen type(s) removed
  sourcing: ChEMBL CHEMBL3804751 -> 25 actives
  triage: Tier 1 (relaxed) (14 checks) -> 25/25 passed
  docking: deferred to cluster -> .../04_docking/submit_vina.sh
  generative: config and submission script prepared -> .../05_generative/submit_reinvent.sh
  selection: Tier 2 (strict) (24 checks) -> 5/25 passed
  dynamics: charmm36 at 310 K, 100 ns -> .../07_dynamics/submit_md.sh
  completed in 21.3s
```

**Pass criteria**

1.  **Nothing is faked.** `docking: deferred to cluster` — not a score, not
    a placeholder number. A pipeline that invented affinities would be
    worse than one that stopped.
2.  The whole seven-stage pipeline completes in **21.3 s** without any of
    the three heavy tools, producing 20 artefacts.
3.  **The funnel narrows for a reason:** 25 actives → **25/25** pass
    relaxed → **5/25** pass strict. A two-tier screen where both tiers give
    the same count is not two tiers. See the *compound-screening* protocol
    T3 for why marketed drugs also fail the strict tier.
4.  Every generated script is valid and portable:

```bash
for f in $(find $W/full -name "*.sh"); do bash -n $f && echo "OK ${f##*/full/}"; done
python -c "
import glob
for f in sorted(glob.glob('$W/full/**/*.sh', recursive=True)):
    print(open(f,'rb').read().count(b'\r'), f.split('full/')[-1])"
```

**Measured:** four scripts, all `bash -n` clean, **all 0 CR bytes** —
`04_docking/submit_vina.sh`, `05_generative/submit_reinvent.sh`,
`07_dynamics/setup_md.sh`, `07_dynamics/submit_md.sh`.

5.  The dynamics stage emits the full MDP set (`em`, `ions`, `md`, `npt`,
    `nvt`), a `run_md.py` and a `setup_md.sh` — the same artefacts the
    *md-simulation* skill produces alone, so its whole protocol applies
    here too.
6.  `submit_vina.sh` carries **no `--gres=gpu`** line. Vina is CPU-only;
    this was a real defect in the shared HPC module, recorded in the
    *md-simulation* protocol T7.

--------------------------------------------------------------------------------

## T5 — data actually flows between stages

An orchestrator whose stages do not feed each other is a batch file.

```bash
grep -E "center_|size_|SLURM_ARRAY" $W/full/04_docking/submit_vina.sh | head -6
grep -E "smiles_file" $W/full/05_generative/reinvent.toml
wc -l < $W/full/03_triage/seeds.smi
```

**Measured:**

```
print(s['center_x'],s['center_y'],s['center_z'],s['size_x'],s['size_y'],s['size_z'])
echo "Box centre ($CX $CY $CZ) size ($SX $SY $SZ)"
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
    --center_x "$CX" --center_y "$CY" --center_z "$CZ" \
    --size_x "$SX" --size_y "$SY" --size_z "$SZ" \

smiles_file = "../03_triage/seeds.smi"
25
```

**Pass criteria**

1.  **The docking script reads the box from `binding_site.json` at job
    time** rather than having the six numbers baked in as literals. If the
    site is corrected, the script does not need regenerating.
2.  **The generative config points at the triage stage's seeds**, all 25 of
    them. Stage 5 consumes stage 3's output.
3.  `--array=1-50` partitions the ligand set, with each task taking every
    Nth ligand.

### T5b — the artefacts are relocatable

The scripts target a Linux cluster but were generated here. The pipeline
directory therefore has to survive being copied.

```bash
python -c "
import os, re
cfg='$W/full/05_generative'
rel=re.search(r'\"(.+?)\"',
      [l for l in open(cfg+'/reinvent.toml') if 'smiles_file' in l][0]).group(1)
print('reference :', rel)
print('resolves  :', os.path.isfile(os.path.join(cfg, rel)))"
```

**Measured:**

```
reference : ../03_triage/seeds.smi
resolves  : True
```

**Pass criteria**

1.  The reference is **relative, with forward slashes**, and resolves from
    the config's own directory.
2.  `submit_reinvent.sh` refers to its config as `CONFIG="reinvent.toml"`,
    relative, so the stage directory is self-contained.

**This was a real defect.** The config embedded an absolute
`C:\Users\...\03_triage\seeds.smi`. Copied to a cluster it would not
resolve, and REINVENT would fail at **load** — after the queue wait, not at
submission. Seed paths are now written relative to the config.

--------------------------------------------------------------------------------

## T6 — refusals and empty stages

```bash
uv run $S run --pdb-id 6HEZ --stages target nonsense \
    --workdir $W/p2 --output $W/p2.json >/dev/null 2>&1; echo $?
# 2

uv run $S run --pdb-id 6HEZ --chain A --stages triage \
    --workdir $W/p3 --output $W/p3.json; echo $?
```

**Measured:**

```
  triage: skipped (no compounds)
  completed in 0.0s
rc=0
```

**Pass criteria**

1.  An unknown stage name exits **2** and the usage message lists the valid
    stages.
2.  A stage with nothing to work on is **skipped with its reason**, not
    crashed and not silently reported as having passed. `triage` without
    `sourcing` has no compounds; `skipped (no compounds)` with exit 0 is
    the right answer, because the run did what was asked.
3.  No partial `relaxed_admet.json` is left behind implying a screen ran.

--------------------------------------------------------------------------------

## T7 — output size

The pipeline report is read back by whatever called the wizard.

```bash
python -c "
import json; d=json.load(open('$W/full.json'))
print('keys:', sorted(d.keys()))
print('artefacts:', len(d['artefacts']), '| in flight:', d['compounds_in_flight'])"
ls -l $W/full.json
```

**Measured:**

| Report | Size | Artefacts | Stages |
|---|---:|---:|---:|
| `pipeline.json` (T3) | 15,110 B | 5 | 3 |
| `full.json` (T4) | 17,705 B | 10 | 7 |

**Pass criteria**

1.  **The report barely grows with the work done** — 2,595 bytes more for
    four extra stages, two of which produced a dozen files. That holds
    only because the report keeps **verdicts, counts and paths** and not
    the contents of what it points at. The whole selection stage is seven
    numbers:

```bash
python -c "
import json; print(json.dumps(json.load(open('$W/full.json'))['stages']['selection']))"
# {"status": "ok", "tier": "Tier 2 (strict)", "checks_per_compound": 24,
#  "total": 25, "passed": 5, "failed": 20, "errored": 0}
```

    The per-compound ADMET detail lives in
    `06_selection/strict_admet.json`. Inlining it here would duplicate
    every compound record and make reading either file cost twice as much.

2.  `compounds_in_flight` (measured: 25) and `artefacts` (5 for the
    three-stage run, 10 for the seven-stage one) summarise the run in two
    integers.
3.  A `note` field states which stages were deferred:

    > This pipeline was orchestrated by the drug-discovery-wizard skill.
    > Stages marked 'deferred' or 'prepared' produced cluster scripts
    > rather than running GPU work locally.

--------------------------------------------------------------------------------

## Results

| Case | What it checks | Result |
|---|---|---|
| T1 | per-stage readiness, plan is read-only | pass — 4 ready / 3 not, 1 file, 1.97 s |
| T2 | stage graph documented | pass — 7 stages with needs and produces |
| T3 | local pipeline end to end | pass — 25/25 triaged in 16.7 s |
| T3b | delegates, does not reimplement | pass — site identical to 1e-12 |
| T4 | cluster stages defer honestly | pass — 20 artefacts, nothing faked, 21.3 s |
| T5 | data flows between stages | pass — box read at job time, 25 seeds consumed |
| T5b | artefacts relocatable | pass — **after a fix**; see below |
| T6 | refusals and empty stages | pass — rc 2 unknown stage, skip with reason |
| T7 | report size | pass — paths and counts, not inlined contents |

**T5b found a real defect.** The generated REINVENT config referenced the
seed file by absolute Windows path, so the pipeline directory could not be
copied to the cluster its own submission script targets. Seed paths are now
written relative to the config, with forward slashes.

**Environment:** Windows 11 Pro, Python 3.13.13, RDKit 2026.03.6,
Meeko 0.8.0. AutoDock Vina, REINVENT 4 and GROMACS **all absent** — which
is what T1 and T4 verify the wizard handles. Live services: RCSB PDB,
ChEMBL.

## Automated equivalent

```bash
pytest tests/test_wizard.py -v
pytest tests/test_conventions.py -v        # packaging and token budgets
pytest tests/test_wizard.py -m network     # T3, T4
```

The manual protocol exists because it is readable, needs no test
environment, and shows what a correct result looks like. T1's read-only
check and T3b's exact-equality check are the two that would catch the most
damaging regressions.
