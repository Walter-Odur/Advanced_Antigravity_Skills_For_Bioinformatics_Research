---
name: compound-synthesis
description: >
  Use this skill when the user wants to generate novel molecules with
  REINVENT 4 generative chemistry. Covers de novo design, scaffold hopping,
  R-group replacement, linker design, peptide design, molecule
  optimisation, transfer learning onto known actives, and reinforcement
  learning with multi-objective scoring. Generates validated TOML configs,
  prepares seeds, runs locally or on a cluster, and analyses the output
  with scaffold-diversity and learning diagnostics.
---

# Compound Synthesis (REINVENT 4)

## Overview

Generates novel drug-like molecules with REINVENT 4, AstraZeneca's
generative chemistry framework. This skill writes validated configs,
prepares inputs, submits or runs jobs, and analyses results.

**Do NOT use when:**

- The user wants to filter or dock existing compounds - use
  `compound-screening`.
- The user wants retrosynthetic route planning - suggest AiZynthFinder,
  ASKCOS or IBM RXN.
- The user wants molecular dynamics - use `md-simulation`.

## Prerequisites

| Requirement | Needed for |
|---|---|
| Python >= 3.10, `rdkit` | config generation, seed preparation, analysis |
| A REINVENT 4 checkout, installed | `run`, `score-molecules` |
| PyTorch, and a GPU in practice | running anything |
| Prior model files | running anything |

**Everything except `run` works without REINVENT installed.** Config
generation, seed preparation and result analysis are local.

### Locating REINVENT

Nothing is hardcoded. Resolution order:

1. `--reinvent-dir`
2. `$REINVENT_DIR`
3. the installed `reinvent` package's location
4. a `REINVENT4` directory beside the working tree

Prior models are found via `$REINVENT_PRIOR_BASE`, or `<REINVENT4>/priors`.
Download them from [Zenodo](https://doi.org/10.5281/zenodo.15641296).

**Run `check-setup` first.** It reports what is present, what is missing,
and exactly how to obtain it, without failing.

## Core Rules

1. **Never hand-write a REINVENT TOML.** Use `generate-config`: it
   validates the combination, resolves the prior, refuses a generator
   whose seeds are missing, and warns when a request will be impractically
   slow.
2. **`--output` is required** for every subcommand.
3. **`check-setup` before any run.**
4. **`prepare-seeds` before using seed structures.** It validates,
   canonicalises and de-duplicates, and reports every rejection with a
   reason.
5. **Use the cluster for real work.** Reinforcement learning is GPU work.
   Above about 100 steps, generate a submission script.
6. **Keep the diversity filter on.** It is what stops the agent collapsing
   onto a single scaffold, the most common failure of an unconstrained run.
7. If this skill is used, mention it in the output.

## Subcommands

| Subcommand | Purpose | Needs REINVENT? |
|---|---|:---:|
| `check-setup` | Verify install, priors, device | no |
| `list-profiles` | Scoring profiles, generators, run modes | no |
| `prepare-seeds` | Validate and de-duplicate seed structures | no |
| `generate-config` | Write a validated TOML | no |
| `preflight` | Check a config and its inputs before running | no |
| `run` | Execute REINVENT | yes |
| `score-molecules` | Score existing molecules with REINVENT scoring | yes |
| `analyze-results` | Rank output and diagnose the run | no |
| `generate-hpc-script` | SLURM submission script | no |

Invoke as either:

```bash
ag-compound-synthesis <subcommand> --output <file> [options]
python .agents/skills/compound-synthesis/scripts/compound_synthesis_api.py <subcommand> ...
```

## Quick Start

```bash
# 1. What is installed?
ag-compound-synthesis check-setup --output setup.json

# 2. A de novo reinforcement-learning config
ag-compound-synthesis generate-config \
    --mode staged_learning --generator reinvent \
    --scoring-profile drug-like \
    --num-steps 300 --batch-size 128 \
    --output run.toml

# 3. Check before committing a GPU allocation
ag-compound-synthesis preflight --config run.toml --output preflight.json

# 4. Submit
ag-compound-synthesis generate-hpc-script \
    --config run.toml --gpu a100 --ngpu 1 --mem 64G --time 24:00:00 \
    --output submit.sh
# sbatch submit.sh

# 5. Analyse
ag-compound-synthesis analyze-results \
    --csv staged_learning_1.csv --top-n 25 --output top.json
```

`generate-config` also writes `<name>.meta.json` beside the TOML,
recording the resolved settings, the expected output files, any warnings,
and the next commands to run.

## Generators

| Generator | Architecture | Seeds | Use for |
|---|---|:---:|---|
| `reinvent` | RNN | no | Unconstrained de novo design |
| `mol2mol` | Transformer | yes | Optimising or hopping from a known compound |
| `libinvent` | Transformer | yes | R-group decoration of a fixed scaffold |
| `linkinvent` | Transformer | yes | Linker design between two fragments |
| `pepinvent` | Transformer | yes | Peptide design |

Seed formats differ, and a mismatch is reported: LibInvent wants scaffolds
with `*` attachment points, LinkInvent wants two warheads per line
separated by `|`. A generator that needs seeds is **refused** without them
rather than silently running unconditioned.

Mol2Mol has variants selectable by short name with `--prior`:
`similarity`, `medium_similarity`, `high_similarity`, `mmp`, `scaffold`,
`scaffold_generic`. Use `scaffold_generic` for scaffold hopping.

## Run Modes

| Mode | Purpose | Requires |
|---|---|---|
| `staged_learning` | Reinforcement learning, the main optimisation mode | a scoring definition |
| `transfer_learning` | Focus a prior onto known actives | `--smiles-file` |
| `sampling` | Draw molecules without optimising | a model |
| `scoring` | Score an existing set | `--smiles-file` and a scoring definition |

## Scoring Profiles

```bash
ag-compound-synthesis list-profiles --output profiles.json
```

| Profile | Targets |
|---|---|
| `drug-like` | Oral small-molecule space: Lipinski-compatible, good QED, synthesisable |
| `lead-like` | Smaller and less lipophilic, leaving room to optimise |
| `fragment-like` | Rule-of-three fragment space |
| `kinase-inhibitor` | Larger, flatter, aromatic systems for the hinge |
| `anti-tb` | Higher lipophilicity, as the mycobacterial envelope demands |
| `cns` | Small, low polar surface area, few donors |

Every numeric component carries a **transform** mapping the raw property
into [0, 1]. This is not decoration: without it the aggregate is
meaningless, and a hard threshold gives the optimiser no gradient to
follow. `double_sigmoid` for a preferred range, `reverse_sigmoid` where
lower is better, `sigmoid` where higher is better. Aggregation is the
geometric mean, so one unsatisfied objective drags the total down rather
than being averaged away.

Add or override components:

```bash
ag-compound-synthesis generate-config \
    --mode staged_learning --generator reinvent \
    --scoring-profile kinase-inhibitor \
    --component "NumAromaticRings:low=2:high=4:weight=2" \
    --component "SAScore:transform=reverse_sigmoid:low=1:high=4" \
    --output run.toml
```

### A property profile is not an activity model

This matters and the skill says so in its output. A physicochemical
profile shapes *where in property space* the generator searches; it says
nothing about whether a molecule binds the target. The `anti-tb` profile
is explicit about this: whole-cell *M. tuberculosis* potency depends on
envelope permeability and efflux, which no descriptor captures. For
target-directed generation, add a structure- or model-based component -
docking through DockStream, or a QSAR model through ChemProp or Qptuna.

## Common Workflows

### De novo design

```bash
ag-compound-synthesis generate-config \
    --mode staged_learning --generator reinvent \
    --scoring-profile drug-like --num-steps 300 --batch-size 128 \
    --output denovo.toml
```

### Transfer learning onto known actives

```bash
# Seeds first: validated, canonicalised, de-duplicated
ag-compound-synthesis prepare-seeds \
    --csv chembl_actives.csv --smiles-col canonical_smiles \
    --report seeds_report.json --output actives.smi

ag-compound-synthesis generate-config \
    --mode transfer_learning --generator reinvent \
    --smiles-file actives.smi --validation-smiles-file holdout.smi \
    --num-epochs 50 --output tl.toml
```

Always pass `--validation-smiles-file`. Without a held-out set there is no
way to detect over-fitting to the training molecules, and the command warns
accordingly. Transfer learning typically needs several hundred molecules to
shift a prior meaningfully.

The focused model it produces becomes the `--prior` of a subsequent
`staged_learning` run, so optimisation happens inside the relevant
chemical space.

### Scaffold hopping

```bash
ag-compound-synthesis generate-config \
    --mode staged_learning --generator mol2mol \
    --prior scaffold_generic --smiles-file lead.smi \
    --scoring-profile drug-like --num-steps 200 --output hop.toml
```

### Seeding from screening results

Relaxed-tier survivors from `compound-screening` make good inception
seeds, which guide the early phase of a reinforcement-learning run:

```bash
ag-compound-synthesis generate-config \
    --mode staged_learning --generator reinvent \
    --scoring-profile drug-like \
    --inception-smiles-file seeds.smi \
    --num-steps 300 --output guided.toml
```

## Interpreting `analyze-results`

```bash
ag-compound-synthesis analyze-results \
    --csv staged_learning_1.csv --top-n 25 --output top.json
```

Three diagnostics matter more than the ranking itself.

**`unique_molecules` against `rows_read`.** A reinforcement-learning run
re-samples the same molecule many times. A "top 20" taken straight from the
CSV is often a handful of structures repeated, which overstates how much
was discovered. Molecules are collapsed on canonical SMILES, and
`times_sampled` shows how often each was produced.

**`scaffold_diversity_ratio`.** Distinct Bemis-Murcko scaffolds over valid
structures. Below about 5 per cent the agent has collapsed onto very few
frameworks; strengthen the diversity filter. `interpretation.diversity`
states the reading.

**`score_progression` and `interpretation.learning`.** Mean score by step.
A flat trend means the agent converged early, sigma is too low, or the
objective is unreachable - none of which a ranking alone would reveal.

The top molecules carry full descriptors and their Murcko scaffold, so they
can be handed straight to strict-tier screening.

## Running

```bash
# Locally, for a small job
ag-compound-synthesis run --config run.toml \
    --log reinvent.log --output report.json

# Check only
ag-compound-synthesis run --config run.toml --dry-run --output check.json
```

`run` performs preflight checks first: the TOML parses, every file the
config references exists, and the requested device is available. A missing
prior is reported in seconds rather than after several minutes of startup.

Reinforcement learning on CPU is impractically slow, and
`generate-config` warns when the request implies it. The molecule budget
is stated up front: 300 steps at batch 128 evaluates roughly 38,400
molecules.

## Reference

- [reference/scoring_components.md](reference/scoring_components.md) -
  components, parameters and transforms, taken from the REINVENT 4
  checkout's own `configs/SCORING.md`
- [examples/](examples/) - worked configs
- REINVENT 4: <https://github.com/MolecularAI/REINVENT4>
- Prior models: <https://doi.org/10.5281/zenodo.15641296>

## Workflow

1. `check-setup`.
2. `prepare-seeds`, if the generator or mode needs structures.
3. `generate-config`, and read its warnings and `.meta.json`.
4. `preflight`.
5. `generate-hpc-script` and submit, or `run` for a small job.
6. `analyze-results`, and read the diversity and learning diagnostics.
7. Hand the top molecules to `compound-screening admet-filter
   --strictness strict`.
