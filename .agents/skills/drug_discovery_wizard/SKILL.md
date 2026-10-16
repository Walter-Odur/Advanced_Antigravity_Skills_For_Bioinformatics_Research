---
name: drug-discovery-wizard
description: >
  Orchestrate the complete structure-based drug discovery pipeline across
  seven stages, from receptor preparation to molecular dynamics, with
  two-tier ADMET screening either side of generative chemistry. Coordinates
  the target-preparation, compound-screening, compound-synthesis and
  md-simulation skills. Use when the user wants the end-to-end pipeline. For
  a single step, use that sub-skill directly.
---

# Drug Discovery Wizard

This skill **delegates**; it does not reimplement. Every stage calls the
same code the individual skills call, so a correction reaches all callers.

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions.
2.  **The four sub-skills**: `target-preparation`, `compound-screening`,
    `compound-synthesis`, `md-simulation`. Read the relevant one before
    interpreting that stage's output.
3.  **User Notification**: If `.licenses/drug_discovery_wizard_LICENSE.txt`
    does not exist in the workspace root, prominently notify the user to
    review the terms listed in
    [SKILL_LICENSES.md](../SKILL_LICENSES.md), then create that file
    recording the notification text and a timestamp.

## Core Rules

-   [!IMPORTANT] **Always call the sub-skills' scripts.** Never write raw API
    calls, GROMACS commands or docking invocations by hand.

-   **Output to File (Required)**: `--output` is required for every
    subcommand, here and in every sub-skill.

-   **Two-tier ADMET, in order**: `relaxed` before generative chemistry,
    `strict` after. Applying strict thresholds to seeds discards chemistry
    the generator could have repaired.

-   **Heavy compute goes to the cluster.** A stage that cannot run locally
    reports `prepared` or `deferred` and emits a submission script rather
    than pretending.

-   **`plan` installs nothing.** It reports readiness only.

-   **Report which stages ran, which were skipped, and why.**

-   **Notification**: If this skill is used, mention it in the output.

## Utility Script

```bash
uv run scripts/drug_discovery_api.py <subcommand> --output <file> [options]
```

| Subcommand | Purpose |
|---|---|
| `plan` | Describe the pipeline and check readiness, without running |
| `run` | Execute the selected stages |
| `describe` | Document the stages, tiers and sub-skills |

--------------------------------------------------------------------------------

### 1. Plan First

```bash
uv run scripts/drug_discovery_api.py plan \
    --pdb-id 6HEZ --chembl-target CHEMBL3804751 --output plan.json
```

Reports, per stage, whether this machine can run it and what is missing:
PDBQT converters, AutoDock Vina, REINVENT and its priors, GROMACS. Nothing
is attempted before you know that.

--------------------------------------------------------------------------------

### 2. The Seven Stages

| # | Stage | Sub-skill | Where | Produces |
|---|---|---|---|---|
| 1 | `target` | target-preparation | local | cleaned PDB, receptor PDBQT, `binding_site.json`, quality report |
| 2 | `sourcing` | compound-screening | local | actives with potency, natural products |
| 3 | `triage` | compound-screening | local | relaxed ADMET report (14 checks), `seeds.smi` |
| 4 | `docking` | compound-screening | needs Vina | ranked affinities, pose files |
| 5 | `generative` | compound-synthesis | cluster | `reinvent.toml`, `submit_reinvent.sh` |
| 6 | `selection` | compound-screening | local | strict ADMET report (24 checks) |
| 7 | `dynamics` | md-simulation | cluster | MDP files, workflow, `submit_md.sh` |

--------------------------------------------------------------------------------

### 3. Run

```bash
# Local stages
uv run scripts/drug_discovery_api.py run \
    --pdb-id 6HEZ --chain A \
    --chembl-target CHEMBL3804751 --pchembl-min 7 \
    --natural-product-query "antimycobacterial" \
    --stages target sourcing triage \
    --workdir pipeline --output pipeline.json

# Cluster preparation
uv run scripts/drug_discovery_api.py run \
    --input-file pipeline/01_target/6HEZ_clean.pdb \
    --stages generative dynamics --temperature 310 --production-ns 100 \
    --workdir pipeline --output cluster.json
```

Writes numbered stage directories under `--workdir`: `01_target/`,
`02_sourcing/`, `03_triage/`, `04_docking/`, `05_generative/`,
`06_selection/`, `07_dynamics/`.

| Stage status | Meaning |
|---|---|
| `ok` | Ran to completion here |
| `prepared` | Cluster scripts written; submit them |
| `deferred` | Cannot run here, so a script was written instead |
| `skipped` | Inputs absent; the reason is recorded |
| `failed` | Errored; the message is recorded |

A failing stage stops the pipeline unless `--keep-going` is passed.

--------------------------------------------------------------------------------

### 4. Why Two Tiers

Before generative chemistry the question is *"what is worth learning
from?"* A liability the generator can optimise away is no reason to discard
a scaffold, and narrowing too early leaves nothing to seed with. After, the
question is *"what do we take forward?"* — and there the published
thresholds apply directly.

Measured on the eight real DprE1 actives ChEMBL returns at pChEMBL >= 6:
8 of 8 pass `relaxed`, 1 of 8 passes `strict`.

See
[references/admet_thresholds.md](references/admet_thresholds.md), generated
from the check registry.

## Common Options

-   `--output FILE`: **Required.** Output path.
-   `--pdb-id` | `--uniprot-id` | `--input-file`: the target structure.
-   `--chain ID`, `--chembl-target ID`, `--natural-product-query TEXT`.
-   `--stages S [S ...]`: which stages to plan or run.
-   `--workdir DIR`: stage output root (default: `pipeline`).
-   `--keep-going`: continue after a stage fails.
-   `--detail compact|standard|full`: how much to write.

## Common Mistakes

1.  **Using the wizard for one step.** Call that sub-skill directly.
2.  **Docking into the cofactor site.** Check
    `stages.target.binding_site.selected_ligand.is_cofactor` before stage 4.
3.  **Reading a stage report at `--detail full` into context.** A
    200-compound screen is ~280,000 tokens that way. Keep the default.
4.  **Treating `prepared` as done.** Those stages wrote a script; submit it.

## Interpreting the Pipeline

A ranked candidate list rests on three independent signals, none sufficient
alone:

-   **Binding affinity** — triage only; sub-kcal/mol differences are not
    meaningful and a Vina score is not a binding free energy.
-   **ADMET profile** — computed, not measured. Report which checks failed;
    penicillin G fails `strict` on its beta-lactam alone.
-   **Dynamic stability** — whether the pose survives unrestrained
    dynamics, which no static score can tell you.

Report the provenance of each: the structure and its confidence or
resolution, the docking box and its origin, the tier applied, and the force
field, temperature and length of any simulation.

## Do Not Use When

-   The user wants one pipeline step — use that sub-skill.
-   The task is wet-lab protocol design.
-   The task is clinical trial design — use the `clinical-trials-database`
    skill.
-   The task is retrosynthetic route planning — suggest AiZynthFinder,
    ASKCOS or IBM RXN.

## Reference

-   [references/questions.md](references/questions.md) — **read this first when the user's request is vague, phrased as a question, or mentions a target, organism or compound class not covered above.** Indexed by question, and it links to each sub-skill's own index.
-   [references/admet_thresholds.md](references/admet_thresholds.md)
-   [references/citation.bib](references/citation.bib) — the union of the
    four sub-skills' citations.

## Workflow

1.  `plan`, and resolve anything reported not ready.
2.  `run --stages target sourcing triage`.
3.  Confirm the binding site is the intended pocket.
4.  Dock locally, or submit the generated array job.
5.  `run --stages generative`, then submit.
6.  `run --stages selection` on what comes back.
7.  `run --stages dynamics`, then submit.
8.  Report ranked candidates with the provenance of every signal.
