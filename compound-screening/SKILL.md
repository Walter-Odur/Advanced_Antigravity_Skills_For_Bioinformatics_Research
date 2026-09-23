---
name: compound-screening
description: >
  Use this skill when the user wants to search compound libraries, dock
  molecules, or filter compounds by ADMET properties. Handles ChEMBL queries,
  COCONUT natural product search, ZINC20 tranche downloads, AutoDock Vina
  docking, and two-tier ADMET filtering (relaxed for REINVENT seeds, strict
  for downstream analysis). Do not use for generative chemistry or de novo
  design — use the compound-synthesis skill instead.
---

# Compound Screening

## Prerequisites

1.  **Python 3.10+** with `requests` and `rdkit` available.
2.  **AutoDock Vina** (optional): Required only for the `dock` subcommand. If
    unavailable, the script will report an error with installation instructions.
3.  **Internet**: Required for `query-chembl`, `query-coconut`, and
    `admet-predict` (live API calls). Offline subcommands: `admet-filter`,
    `query-zinc`, `dock`.

## Overview

Compound sourcing, molecular docking, and ADMET filtering. Supports a two-tier
filtering strategy: relaxed thresholds to maximize hit diversity and strict
thresholds for rigorous downstream selection.

**Do NOT use when:**

-   The user wants to generate novel molecules (use `compound-synthesis`).
-   The user wants to set up MD simulations (use `md-simulation`).
-   The user wants to prepare a receptor (use `target-preparation`).

## Core Rules

-   **Use the Wrapper**: ALWAYS execute the provided utility script for all
    compound queries, docking, and ADMET filtering. NEVER use `curl` or custom
    Python requests to query APIs directly. The script enforces rate limits and
    handles retries.
-   **Output to File (Required)**: The `--output` flag is **required** for every
    subcommand. All JSON results are written to the specified file. After running
    the command, read the output file to extract the data.
-   **Two-Tier Filtering**: Use `--strictness relaxed` for broad filtering
    (14 checks) and `--strictness strict` (default) for rigorous filtering
    (34 checks). See `reference/admet_thresholds.md` for full threshold details.
-   If this skill is used, ensure this is mentioned in the output.

## Utility Script

All compound screening uses one script with subcommands:

```bash
uv run scripts/compound_screening_api.py <subcommand> --output <file> [options]
```

--------------------------------------------------------------------------------

### 1. Search ChEMBL Targets

```bash
uv run scripts/compound_screening_api.py search-chembl-target \
    --query "DprE1" --limit 10 --output /path/to/targets.json
```

--------------------------------------------------------------------------------

### 2. Query ChEMBL Compounds

```bash
uv run scripts/compound_screening_api.py query-chembl \
    --target-id CHEMBL203 --pchembl-min 7 --limit 20 \
    --output /path/to/chembl_hits.json
```

--------------------------------------------------------------------------------

### 3. Query COCONUT Natural Products

Searches COCONUT 2.0 (400k+ natural products, includes ANPDB, CMNPD, NANPDB,
AfroDB, SANCDB, ConMedNP).

```bash
uv run scripts/compound_screening_api.py query-coconut \
    --query "manzamine" --limit 20 --output /path/to/coconut.json
```

--------------------------------------------------------------------------------

### 4. Generate ZINC20 Download Script

Generates a shell script with `wget` commands for bulk ZINC20 tranche downloads.

```bash
uv run scripts/compound_screening_api.py query-zinc \
    --subset drug-like --output /path/to/download_zinc.sh

# Custom range
uv run scripts/compound_screening_api.py query-zinc \
    --subset drug-like --mw-range 300 450 --logp-range 0 4 \
    --format sdf --output /path/to/download_zinc.sh
```

--------------------------------------------------------------------------------

### 5. Dock Compounds

Docks compounds against a prepared receptor using AutoDock Vina.

```bash
uv run scripts/compound_screening_api.py dock \
    --receptor receptor.pdbqt \
    --smiles "CCO" "CC(=O)Oc1ccccc1C(=O)O" \
    --names "ethanol" "aspirin" \
    --site binding_site.json \
    --exhaustiveness 32 \
    --output /path/to/docking.json
```

--------------------------------------------------------------------------------

### 6. ADMET Filter (Offline, RDKit-Based)

RDKit-based ADMET filtering with Lipinski rules, PAINS alerts, and SA scoring.
Supports two strictness tiers:

```bash
# Relaxed (broad filtering — 14 checks)
uv run scripts/compound_screening_api.py admet-filter \
    --smiles "CC(=O)Oc1ccccc1C(=O)O" --names "aspirin" \
    --strictness relaxed --output /path/to/seeds.json

# Strict (rigorous filtering — 34 checks, default)
uv run scripts/compound_screening_api.py admet-filter \
    --smiles "CC(=O)Oc1ccccc1C(=O)O" --names "aspirin" \
    --strictness strict --output /path/to/candidates.json
```

> **Important**: The output includes `strictness` and `tier` metadata. Always
> report the tier used to the user.

--------------------------------------------------------------------------------

### 7. Full ADMET Prediction (Online, ADMETlab 3.0)

Full ADMET prediction via ADMETlab 3.0 API (119 endpoints). Requires internet.

```bash
uv run scripts/compound_screening_api.py admet-predict \
    --smiles "CC(=O)Oc1ccccc1C(=O)O" --names "aspirin" \
    --output /path/to/admet_full.json
```

## Interpreting the Output

1.  **ADMET Filter**: Each compound gets `overall_pass` (true/false), individual
    checks (`lipinski_pass`, `pains_pass`, `sa_pass`), and computed properties
    (`mw`, `logp`, `hbd`, `hba`, `tpsa`, `rotatable_bonds`, `sa_score`).
    The `tier` field indicates which filtering level was applied.
2.  **Docking**: Results include binding affinity (kcal/mol), pose coordinates,
    and RMSD values. Lower (more negative) affinity = stronger binding.
3.  **Library queries**: Results include SMILES, compound names, and source
    metadata. Always check `total_results` for result completeness.

## Common Options

-   `--output FILE`: **Required.** Output file path for JSON results.
-   `--smiles SMILES [SMILES ...]`: One or more SMILES strings.
-   `--names NAME [NAME ...]`: Compound names (must match `--smiles` count).
-   `--strictness {relaxed,strict}`: ADMET filtering tier (default: strict).
-   `--limit N`: Max results for database queries.

## Reference

-   **ADMET Thresholds**: See
    [reference/admet_thresholds.md](reference/admet_thresholds.md) for the full
    list of thresholds and relaxed vs strict comparison.

## Workflow

1.  Search compound libraries with `search-chembl-target`, `query-chembl`,
    `query-coconut`, or `query-zinc`.
2.  Read the output JSON to extract SMILES and compound names.
3.  Dock candidates with `dock` (small set) or generate an HPC script for
    batch screening.
4.  Filter hits with `admet-filter` using the appropriate strictness tier.
5.  Use `admet-predict` for full ADMET profiling of top candidates.
