---
name: drug-discovery-wizard
description: >
  Orchestrator for the complete structure-based drug discovery pipeline.
  Coordinates four sub-skills - target-preparation, compound-screening,
  compound-synthesis and md-simulation - across seven stages from receptor
  preparation to molecular dynamics, with two-tier ADMET screening either
  side of generative chemistry. Use this skill when the user wants the
  end-to-end pipeline. For a single step, use that sub-skill directly.
---

# Drug Discovery Wizard

## Overview

The orchestrator for the full structure-based pipeline. It coordinates four
sub-skills, each of which also works standalone.

This skill **delegates**; it does not reimplement. Every stage calls the
same library functions the individual skills call, so a correction reaches
all callers at once.

**Do NOT use when:**

- The user wants one pipeline step. Call that skill directly:
  `target-preparation`, `compound-screening`, `compound-synthesis` or
  `md-simulation`.
- The user wants wet-lab protocol design.
- The user wants clinical trial design - use the
  `clinical-trials-database` skill.
- The user wants retrosynthetic route planning - suggest AiZynthFinder,
  ASKCOS or IBM RXN.

## Sub-Skills

| Skill | Console command | Responsibility |
|---|---|---|
| `target-preparation` | `ag-target-preparation` | Structures, cleaning, confidence, binding sites |
| `compound-screening` | `ag-compound-screening` | Libraries, docking, ADMET |
| `compound-synthesis` | `ag-compound-synthesis` | REINVENT 4 generative chemistry |
| `md-simulation` | `ag-md-simulation` | GROMACS setup, cluster scripts |

## Core Rules

1. **Always call the sub-skills' scripts.** Never write raw API calls,
   GROMACS commands or docking invocations by hand.
2. **`--output` is required** for every subcommand everywhere.
3. **Two-tier ADMET, in the right order.** Relaxed before generative
   chemistry, strict after. Applying strict thresholds to seeds discards
   chemistry the generator could have repaired.
4. **Heavy compute goes to the cluster.** REINVENT, GROMACS and batch Vina
   produce submission scripts. A stage that cannot run locally reports
   `prepared` or `deferred` rather than pretending.
5. **Report which stages ran, which were skipped, and why.** Every stage
   records a status.
6. If this skill is used, mention it in the output.

## Subcommands

| Subcommand | Purpose |
|---|---|
| `plan` | Describe the pipeline and check what is ready, without running |
| `run` | Execute the selected stages |
| `describe` | Document the stages, tiers and sub-skills |

```bash
ag-wizard <subcommand> --output <file> [options]
python .agents/skills/drug-discovery-wizard/scripts/drug_discovery_api.py <subcommand> ...
```

## The Pipeline

```
1. TARGET          target-preparation      local
   Fetch from PDB or AlphaFold, or predict with ESMFold. Clean. Assess
   confidence. Convert to PDBQT. Define the docking box from the bound
   ligand.
   -> cleaned PDB, receptor PDBQT, binding_site.json, quality report

2. SOURCING        compound-screening      local
   Known actives from ChEMBL, aggregated per molecule. Natural products
   from COCONUT. Purchasable compounds from ZINC.
   -> compounds with potency and provenance

3. TRIAGE          compound-screening      local
   Relaxed ADMET: 14 checks. Wide windows, tolerant of flagged
   liabilities, because the goal is a diverse seed set.
   -> triage report, seeds.smi

4. DOCKING         compound-screening      needs Vina
   AutoDock Vina against the prepared receptor and box. Seeded, so the
   run is reproducible.
   -> ranked affinities, pose files

5. GENERATIVE      compound-synthesis      cluster
   REINVENT 4 config and submission script, seeded from stage 3.
   -> reinvent.toml, submit_reinvent.sh

6. SELECTION       compound-screening      local
   Strict ADMET: 24 checks. Published thresholds at face value, no
   structural alerts.
   -> final candidate report

7. DYNAMICS        md-simulation           cluster
   GROMACS workflow, MDP files and submission script, with non-bonded
   settings derived from the force field.
   -> md_workflow/, submit_md.sh
```

## Why Two Tiers

The pipeline makes two different decisions and they want different
thresholds.

**Before generative chemistry** the question is "what is worth learning
from?" A liability the generator can optimise away - a borderline logP, a
single structural alert - is not a reason to discard a scaffold. Narrowing
too early leaves nothing to seed with. The relaxed tier therefore runs 14
checks with wide windows, and treats published risk *flags* as advisory
rather than disqualifying.

**After generative chemistry** the question is "what do we take forward?"
Here the strict tier's 24 checks apply published thresholds directly, and
structural alerts matter.

The counts are derived from the check registry, not asserted in prose. See
[reference/admet_thresholds.md](reference/admet_thresholds.md), which is
generated from the same source, and
`ag-compound-screening describe-tiers`.

## Quick Start

```bash
# What would run, and what is ready?
ag-wizard plan --pdb-id 6HEZ --chembl-target CHEMBL3804751 \
    --output plan.json

# Run the local stages
ag-wizard run \
    --pdb-id 6HEZ --chain A \
    --chembl-target CHEMBL3804751 --pchembl-min 7 \
    --natural-product-query "antimycobacterial" \
    --stages target sourcing triage \
    --workdir pipeline --output pipeline.json

# Prepare the cluster work
ag-wizard run \
    --input-file pipeline/01_target/6HEZ_clean.pdb \
    --stages generative dynamics \
    --temperature 310 --production-ns 100 \
    --workdir pipeline --output cluster.json
```

Start with `plan`. It reports, per stage, whether this machine can run it,
and names what is missing: whether a PDBQT converter is installed, whether
Vina is available, whether REINVENT and its priors are present, whether
GROMACS is on PATH. Nothing is attempted before you know that.

## Stage Outputs

`run` writes numbered directories under `--workdir`:

```
pipeline/
  01_target/      structure, cleaned PDB, PDBQT, binding_site.json
  02_sourcing/    chembl_actives.json, coconut.json
  03_triage/      relaxed_admet.json, seeds.smi
  04_docking/     docking.json, poses/, or submit_vina.sh
  05_generative/  reinvent.toml, reinvent.meta.json, submit_reinvent.sh
  06_selection/   strict_admet.json
  07_dynamics/    *.mdp, run_md.py, setup_md.sh, submit_md.sh
```

The report records `stages` (status per stage), `artefacts` (paths other
stages and you will need), and `compounds_in_flight`.

Stage statuses:

| Status | Meaning |
|---|---|
| `ok` | Ran to completion here |
| `prepared` | Cluster scripts written; submit them |
| `deferred` | Cannot run here, so a cluster script was written instead |
| `skipped` | Inputs absent; the reason is recorded |
| `failed` | Errored; the message is recorded |

By default a failing stage stops the pipeline. Pass `--keep-going` to
continue, which is useful when a later stage does not depend on the failure.

## A Worked Example: Anti-TB, DprE1

DprE1 is a validated *M. tuberculosis* target, and 6HEZ is a good test of
the pipeline because it exposes two traps.

```bash
ag-target-preparation prepare-receptor \
    --pdb-id 6HEZ --chain A --output work/receptor.json
```

**Trap one: the cofactor.** 6HEZ contains both an FAD cofactor (53 atoms)
and the inhibitor 0SK (28 atoms). Ranking heterogens by size centres the
docking box on the flavin site rather than the inhibitor pocket, so every
affinity afterwards would be for the wrong cavity. The skill prefers the
drug-like ligand, labels FAD as a cofactor, and lists the alternatives.
Check `selected_ligand` before docking.

**Trap two: incomplete side chains.** 16 residues in chain A have missing
side-chain atoms, which is routine in crystal structures. The receptor is
still built, but the report names the approximated residues and flags any
that lie **inside** the docking box - one does. A defect inside the box
makes affinities for that site less reliable; one outside it is harmless.

Then:

```bash
ag-compound-screening query-chembl \
    --target-id CHEMBL3804751 --pchembl-min 7 --output actives.json

ag-compound-screening admet-filter \
    --csv actives.csv --smiles-col canonical_smiles \
    --strictness relaxed --output triage.json

ag-compound-synthesis generate-config \
    --mode staged_learning --generator reinvent \
    --scoring-profile anti-tb --inception-smiles-file seeds.smi \
    --num-steps 300 --output reinvent.toml
```

The `anti-tb` scoring profile permits higher lipophilicity than the general
drug-like one, because the mycobacterial envelope is lipid-rich and known
actives sit there - bedaquiline is 555 Da at clogP above 7. It is a
**property profile, not an activity model**: whole-cell potency depends on
envelope permeability and efflux, which no descriptor captures. The output
says so. For target-directed generation, add a docking or QSAR component.

## Interpreting the Pipeline

A ranked candidate list rests on three independent signals, and none alone
is sufficient:

- **Binding affinity** from docking. Triage only. Differences under about
  1 kcal/mol are not meaningful, and a Vina score is a scoring function,
  not a binding free energy.
- **ADMET profile** from the strict tier. Computed, not measured. Report
  which checks failed rather than only a pass or fail, since a single
  structural alert can sink an otherwise sound molecule - penicillin G
  fails the strict tier on its beta-lactam alone.
- **Dynamic stability** from MD. Whether the pose survives unrestrained
  dynamics, which a static score cannot tell you.

Report the provenance of each: the structure used, its confidence or
resolution, the docking box and its origin, the tier applied, and the
force field, temperature and length of any simulation.

## Reference

- [reference/admet_thresholds.md](reference/admet_thresholds.md) -
  generated from the check registry

## Workflow

1. `plan`, and resolve anything reported as not ready.
2. `run --stages target sourcing triage`.
3. Confirm the binding site is the pocket you intend.
4. Dock locally, or submit the generated array job.
5. `run --stages generative`, then submit.
6. `run --stages selection` on what comes back.
7. `run --stages dynamics`, then submit.
8. Report ranked candidates with the provenance of every signal.
