---
name: compound-synthesis
description: >
  Generate and optimise novel molecules with REINVENT 4. Covers de novo
  design, scaffold hopping, R-group replacement, linker and peptide design,
  transfer learning onto known actives, and multi-objective reinforcement
  learning. Writes validated TOML configs, prepares seeds, runs locally or
  on a cluster, and analyses output with scaffold-diversity diagnostics. Use
  when the user wants new molecules designed or an existing lead optimised.
---

# Compound Synthesis (REINVENT 4)

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions. Only
    RDKit is declared in the script header, because that is all the local
    subcommands need.
2.  **User Notification**: If `.licenses/compound_synthesis_LICENSE.txt` does
    not exist in the workspace root, prominently notify the user to review
    the terms at https://github.com/MolecularAI/REINVENT4, then create that
    file recording the notification text and a timestamp.
3.  **REINVENT 4** (`run` and `score-molecules` only):

    ```bash
    git clone https://github.com/MolecularAI/REINVENT4
    pip install -e REINVENT4
    export REINVENT_DIR=$PWD/REINVENT4
    ```

    Prior models: https://doi.org/10.5281/zenodo.15641296. This is a git
    checkout plus a multi-gigabyte PyTorch install, so it is **never**
    installed implicitly. Run `check-setup` to see what is missing.

## Core Rules

-   [!IMPORTANT] **Never hand-write a REINVENT TOML.** Use
    `generate-config`: it validates the combination, resolves the prior,
    refuses a generator whose seeds are missing, and warns when a request
    will be impractically slow.

-   **Output to File (Required)**: `--output` is required for every
    subcommand.

-   **`check-setup` before any run.** Config generation, seed preparation and
    result analysis all work without REINVENT installed; only `run` needs it.

-   **`prepare-seeds` before using seed structures.** It validates,
    canonicalises and de-duplicates, reporting every rejection with a reason.

-   **Use the cluster for real work.** Reinforcement learning is GPU work.
    Above roughly 100 steps, generate a submission script.

-   **Keep the diversity filter on.** It is what stops the agent collapsing
    onto a single scaffold — the most common failure of an unconstrained run.

-   **Notification**: If this skill is used, mention it in the output.

## Utility Script

```bash
uv run scripts/compound_synthesis_api.py <subcommand> --output <file> [options]
```

| Subcommand | Purpose | Needs REINVENT? |
|---|---|:---:|
| `check-setup` | Verify install, priors, device | no |
| `list-profiles` | Scoring profiles, generators, run modes | no |
| `prepare-seeds` | Validate and de-duplicate seed structures | no |
| `generate-config` | Write a validated TOML | no |
| `preflight` | Check a config and its inputs before running | no |
| `run` | Execute REINVENT | yes |
| `score-molecules` | Score existing molecules | yes |
| `analyze-results` | Rank output and diagnose the run | no |
| `generate-hpc-script` | SLURM submission script | no |

--------------------------------------------------------------------------------

### 1. Check the Installation

```bash
uv run scripts/compound_synthesis_api.py check-setup --output setup.json
```

Reports what is present, what is missing, and how to obtain it — without
failing. REINVENT is located via `--reinvent-dir`, then `$REINVENT_DIR`, then
the installed package. Priors via `$REINVENT_PRIOR_BASE` or
`<REINVENT4>/priors`. Nothing is hardcoded.

--------------------------------------------------------------------------------

### 2. Generate a Config

```bash
uv run scripts/compound_synthesis_api.py generate-config \
    --mode staged_learning --generator reinvent \
    --scoring-profile drug-like --num-steps 300 --batch-size 128 \
    --output run.toml
```

Also writes `run.meta.json` with the resolved settings, expected outputs,
warnings and the next commands to run.

| Mode | Purpose | Requires |
|---|---|---|
| `staged_learning` | Reinforcement learning — the main mode | a scoring definition |
| `transfer_learning` | Focus a prior onto known actives | `--smiles-file` |
| `sampling` | Draw molecules without optimising | a model |
| `scoring` | Score an existing set | `--smiles-file` + scoring |

| Generator | Seeds | Use for |
|---|:---:|---|
| `reinvent` (RNN) | no | Unconstrained de novo design |
| `mol2mol` | yes | Optimising or hopping from a known compound |
| `libinvent` | yes | R-group decoration of a fixed scaffold |
| `linkinvent` | yes | Linker design between two fragments |
| `pepinvent` | yes | Peptide design |

A generator that needs seeds is **refused** without them. Mol2Mol variants
(`scaffold_generic`, `similarity`, `mmp`, …) select with `--prior`.

--------------------------------------------------------------------------------

### 3. Scoring Profiles

```bash
uv run scripts/compound_synthesis_api.py list-profiles --output profiles.json
```

`drug-like`, `lead-like`, `fragment-like`, `kinase-inhibitor`, `anti-tb`,
`cns`. Add or override components:

```bash
    --component "NumAromaticRings:low=2:high=4:weight=2" \
    --component "SAScore:transform=reverse_sigmoid:low=1:high=4"
```

Every numeric component carries a transform mapping the property into
[0, 1]; without one the aggregate is meaningless and a hard threshold gives
the optimiser no gradient. Aggregation is the geometric mean, so one
unsatisfied objective drags the total down.

> [!IMPORTANT] **A property profile is not an activity model.** It shapes
> *where in property space* the generator searches; it says nothing about
> target binding. For `anti-tb`, whole-cell *M. tuberculosis* potency depends
> on envelope permeability and efflux, which no descriptor captures. Add a
> docking (DockStream) or QSAR (ChemProp, Qptuna) component for
> target-directed generation.

See [references/scoring_components.md](references/scoring_components.md).

--------------------------------------------------------------------------------

### 4. Prepare Seeds

```bash
uv run scripts/compound_synthesis_api.py prepare-seeds \
    --csv chembl_actives.csv --smiles-col canonical_smiles \
    --report seeds_report.json --output actives.smi
```

Canonicalises before de-duplicating, so differently-written forms of one
molecule collapse to a single entry. Every input is accounted for: kept,
invalid, or duplicate, with the reason.

--------------------------------------------------------------------------------

### 5. Run, or Submit

```bash
uv run scripts/compound_synthesis_api.py preflight --config run.toml \
    --output preflight.json

uv run scripts/compound_synthesis_api.py generate-hpc-script \
    --config run.toml --gpu a100 --ngpu 1 --mem 64G --time 24:00:00 \
    --output submit.sh
# sbatch submit.sh
```

One script plus one guide (`submit_SUBMIT.md`): the tools the cluster needs with probe and install commands, files to stage, the `sbatch`/`squeue`/`sacct` commands, and what each failure means. Point the user at it. `--module NAME` is repeatable and makes the job load the cluster's own build.

`run` performs preflight first: the TOML parses, every referenced file
exists, and the device is available. A missing prior is reported in seconds
rather than after minutes of startup.

--------------------------------------------------------------------------------

### 6. Analyse Results

```bash
uv run scripts/compound_synthesis_api.py analyze-results \
    --csv staged_learning_1.csv --top-n 25 --output top.json
```

Three diagnostics matter more than the ranking:

-   **`unique_molecules` against `rows_read`** — an RL run re-samples the
    same molecule many times, so a raw "top 20" is often a handful of
    structures repeated. `times_sampled` shows how often each appeared.
-   **`scaffold_diversity_ratio`** — distinct Bemis-Murcko scaffolds over
    valid structures. Below ~5% the agent has collapsed; strengthen the
    diversity filter.
-   **`score_progression`** — mean score by step. A flat trend means early
    convergence, sigma too low, or an unreachable objective, none of which a
    ranking alone reveals.

## Common Options

-   `--output FILE`: **Required.** Output path.
-   `--mode`, `--generator`, `--prior`, `--scoring-profile`: see above.
-   `--num-steps N`, `--batch-size N`, `--num-epochs N`, `--num-smiles N`.
-   `--device cpu|cuda:0`, `--reinvent-dir DIR`.
-   `--top-n N`: molecules to return from `analyze-results`.
-   `--dry-run`: preflight only; execute nothing.
-   `--detail compact|standard|full`: how much to write.

## Common Mistakes

1.  **Treating a scoring profile as an activity model.** It is not. Add a
    target-based component.
2.  **Reading a raw output CSV as a discovery list.** Use `analyze-results`;
    96 rows collapsed to 13 unique molecules in a typical run.
3.  **Disabling the diversity filter.** Mode collapse follows.
4.  **Running RL on CPU.** `generate-config` warns; submit to a GPU instead.

## Reference

-   [references/questions.md](references/questions.md) — **read this first when the user's request is vague, phrased as a question, or mentions a target, organism or compound class not covered above.** Indexed by question: generator choice, transforms, mode collapse, and what to do when no profile fits the target.
-   [references/scoring_components.md](references/scoring_components.md) —
    components, parameters and transforms, taken from the REINVENT 4
    checkout's own `configs/SCORING.md`.
-   [examples/](examples/) — five worked configs, each valid TOML.
-   [references/citation.bib](references/citation.bib) — cite these for any
    generated molecule.

## Workflow

1.  `check-setup`.
2.  `prepare-seeds`, if the generator or mode needs structures.
3.  `generate-config`; read its warnings and `.meta.json`.
4.  `preflight`.
5.  `generate-hpc-script` and submit, or `run` for a small job.
6.  `analyze-results`; read the diversity and learning diagnostics.
7.  Hand the top molecules to `compound-screening admet-filter
    --strictness strict`.
