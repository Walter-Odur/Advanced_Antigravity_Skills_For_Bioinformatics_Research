---
name: compound-synthesis
description: >
  Use this skill when the user wants to generate novel molecules using
  REINVENT 4 generative chemistry. Covers de novo design, scaffold hopping,
  R-group replacement, linker design, molecule optimization, transfer
  learning, and RL-based multi-objective optimization. Wraps the REINVENT4
  codebase for config generation, execution, and result analysis.
---

# Compound Synthesis — REINVENT 4 Skill

## Overview

Generate novel drug-like molecules with REINVENT 4, an AI-driven generative
chemistry framework from AstraZeneca. This skill wraps the full REINVENT4
pipeline: config generation, seed preparation, local/HPC execution, and
result analysis across all five run modes.

**REINVENT4 location**: `E:\ANTIGRAVITY_WORKSHOP\REINVENT4` (cloned repo)

## Core Rules

1.  **Use the Wrapper**: ALWAYS execute `compound_synthesis.py` rather than
    writing REINVENT TOML configs by hand. The script validates all
    parameters, applies sensible defaults, and generates correct TOML.

2.  **Output Flag Required**: The `--output` flag is **required** for every
    subcommand. Results are written to the specified file.

3.  **Check Setup First**: Before any REINVENT run, call `check-setup` to
    verify the installation, prior models, and device availability.

4.  **HPC for Heavy Jobs**: For runs with >100 RL steps, docking-in-the-loop
    scoring, or transfer learning with >200 epochs, generate an HPC script
    with `generate-hpc-script` rather than running locally.

5.  **Prior Models**: REINVENT4 requires pre-trained prior models from
    [Zenodo](https://doi.org/10.5281/zenodo.15641296). Use `download-priors`
    to fetch them, or set `REINVENT_PRIOR_BASE` to a custom directory.

## Dependencies

- Python >= 3.11
- PyTorch >= 2.12.0
- RDKit >= 2025.09.1
- REINVENT4 installed (`pip install -e .` in the REINVENT4 directory)

## Subcommand Reference

| Subcommand | Purpose | Local? |
|---|---|:---:|
| `check-setup` | Verify install, priors, GPU | ✅ |
| `download-priors` | Fetch models from Zenodo | ✅ |
| `generate-config` | Generate TOML config for any run mode | ✅ |
| `prepare-seeds` | Validate/deduplicate seed SMILES | ✅ |
| `run` | Execute REINVENT locally | ✅ |
| `analyze-results` | Parse CSV output, rank compounds | ✅ |
| `generate-hpc-script` | SLURM submission script | ✅ |
| `score-molecules` | Score SMILES with REINVENT scoring | ✅ |

## Quick Start

```bash
# 1. Check installation
python compound_synthesis.py check-setup \
    --reinvent-dir E:\ANTIGRAVITY_WORKSHOP\REINVENT4 \
    --output setup.json

# 2. Generate de novo RL config with drug-like scoring
python compound_synthesis.py generate-config \
    --mode staged_learning \
    --generator reinvent \
    --prior .reinvent \
    --scoring-profile drug-like \
    --seeds "Cc1ccnc(NC(=O)c2ccsc2NC(=O)c2ccn3cnnc3c2)n1" \
    --num-steps 300 --batch-size 128 \
    --device cpu \
    --output run_config.toml

# 3. Run locally (small job)
python compound_synthesis.py run \
    --config run_config.toml \
    --reinvent-dir E:\ANTIGRAVITY_WORKSHOP\REINVENT4 \
    --log reinvent.log \
    --output run_report.json

# 4. Analyze results
python compound_synthesis.py analyze-results \
    --csv staged_learning_1.csv \
    --top-n 20 --sort-by total_score \
    --output top_compounds.json
```

## Common Workflows

### Workflow 1: De Novo RL (drug-likeness only)

```bash
python compound_synthesis.py generate-config \
    --mode staged_learning --generator reinvent --prior .reinvent \
    --scoring-profile drug-like \
    --num-steps 300 --batch-size 128 --device cpu \
    --output denovo.toml

python compound_synthesis.py run --config denovo.toml \
    --reinvent-dir E:\ANTIGRAVITY_WORKSHOP\REINVENT4 \
    --output denovo_report.json
```

### Workflow 2: Transfer Learning (focus on target actives)

```bash
python compound_synthesis.py prepare-seeds \
    --smiles "CCO" "c1ccccc1" "CC(=O)Oc1ccccc1C(=O)O" \
    --output actives.smi --report seeds_report.json

python compound_synthesis.py generate-config \
    --mode transfer_learning --generator reinvent --prior .reinvent \
    --smiles-file actives.smi \
    --num-epochs 50 --batch-size 64 --device cpu \
    --output tl_config.toml

python compound_synthesis.py run --config tl_config.toml \
    --reinvent-dir E:\ANTIGRAVITY_WORKSHOP\REINVENT4 \
    --output tl_report.json
```

### Workflow 3: Scaffold Hopping (Mol2Mol)

```bash
python compound_synthesis.py generate-config \
    --mode staged_learning --generator mol2mol \
    --prior .m2m_scaffold_generic \
    --smiles-file lead_compound.smi \
    --scoring-profile drug-like \
    --num-steps 200 --batch-size 64 \
    --output scaffold_hop.toml
```

### Workflow 4: Score Existing Molecules

```bash
python compound_synthesis.py score-molecules \
    --smiles-file candidates.smi \
    --scoring-profile drug-like \
    --reinvent-dir E:\ANTIGRAVITY_WORKSHOP\REINVENT4 \
    --output scores.json
```

### Workflow 5: HPC Submission

```bash
python compound_synthesis.py generate-hpc-script \
    --config run_config.toml \
    --gpu a100 --ngpu 1 --mem 64G --time 24:00:00 \
    --conda-env reinvent4 \
    --output submit.sh
```

## Scoring Profiles

Built-in scoring presets for `--scoring-profile`:

| Profile | MW | LogP | QED | HBA | HBD | SA | Alerts |
|---|---|---|---|---|---|---|---|
| `drug-like` | 200–500 | <5 | ✅ | ≤10 | ≤5 | — | PAINS |
| `lead-like` | 200–350 | 0–3 | ✅ | ≤6 | ≤3 | — | PAINS |
| `fragment-like` | 100–250 | -1–3 | — | ≤3 | ≤3 | — | — |
| `kinase-inhibitor` | 300–550 | 1–5 | ✅ | 3–8 | — | <5 | PAINS |
| `anti-tb` | 250–600 | 1–5 | ✅ | — | — | <4 | PAINS |
| `custom` | User-defined components via `--component` flags |

## Run Modes

| Mode | TOML `run_type` | Generator | Needs Seeds? |
|---|---|---|:---:|
| De novo generation | `staged_learning` | Reinvent | No |
| R-group replacement | `staged_learning` | LibInvent | Yes (scaffolds) |
| Linker design | `staged_learning` | LinkInvent | Yes (warheads) |
| Molecule optimization | `staged_learning` | Mol2Mol | Yes (molecules) |
| Peptide design | `staged_learning` | Pepinvent | Yes (peptides) |
| Transfer learning | `transfer_learning` | Any | Yes (SMILES) |
| Sampling | `sampling` | Any | Depends |
| Scoring only | `scoring` | — | Yes (SMILES) |

## Generator Selection Guide

- **Reinvent** (RNN): Best for unconstrained de novo design. No seed needed.
- **Mol2Mol** (Transformer): Best for optimizing a known compound. Needs seed SMILES.
- **LibInvent** (Transformer): Best for finding R-groups on a fixed scaffold.
- **LinkInvent** (Transformer): Best for finding linkers between two fragments.
- **Pepinvent** (Transformer): Best for peptide design.

## File Reference

For full parameter documentation, see:
- `reference/scoring_components.md` — All 24+ scoring components with transforms
- `reference/parameters.md` — All run mode parameters with defaults
- `reference/prior_models.md` — Available models and download instructions
