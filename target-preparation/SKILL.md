---
name: target-preparation
description: >
  Use this skill when the user wants to prepare a protein target for
  structure-based drug discovery. Handles structure prediction (ESMFold, AF2),
  receptor fetching from PDB or AlphaFold Database, quality assessment (pLDDT),
  and receptor preparation for docking (cleaning, chain extraction, PDBQT
  conversion). Do not use if the user only has a gene name — ask for a PDB ID,
  UniProt ID, or amino acid sequence first.
---

# Target Preparation

## Prerequisites

1.  **Python 3.10+** with `requests` and `biopython` available.
2.  **Open Babel or Meeko** (optional): Required only for PDBQT conversion. If
    unavailable, the script will produce cleaned PDB files and advise on
    installing Open Babel.

## Overview

Prepares protein targets for structure-based drug discovery. Takes a PDB ID,
UniProt ID, or amino acid sequence and produces a docking-ready receptor file.

**Do NOT use when:**

-   The user only has a gene or protein name (no PDB/UniProt ID) — ask them
    to look up the ID on [RCSB PDB](https://www.rcsb.org) or
    [UniProt](https://www.uniprot.org).
-   The user wants to run molecular dynamics (use the `md-simulation` skill).
-   The user wants to screen compounds (use the `compound-screening` skill).

## Core Rules

-   **Use the Wrapper**: ALWAYS execute the provided utility scripts to prepare
    receptors rather than writing your own scripts. The scripts handle PDB
    fetching, chain extraction, water/ligand removal, and PDBQT conversion.
-   Do not attempt to assess structure quality yourself; always rely on the
    output provided by the `assess-structure` subcommand.
-   If this skill is used, ensure this is mentioned in the output.

## Utility Scripts

All target preparation uses one script with subcommands:

```bash
uv run scripts/target_preparation_api.py <subcommand> --output <file> [options]
```

--------------------------------------------------------------------------------

### 1. Prepare Receptor

Fetches a structure from PDB or AlphaFold Database, cleans it (removes water,
ligands, non-standard residues), extracts a single chain, and converts to PDBQT.

```bash
# From PDB ID
uv run scripts/target_preparation_api.py prepare-receptor \
    --pdb-id 6HEZ --chain A --output /path/to/receptor.json

# From AlphaFold Database (UniProt ID)
uv run scripts/target_preparation_api.py prepare-receptor \
    --uniprot-id P00533 --output /path/to/receptor.json

# From local file
uv run scripts/target_preparation_api.py prepare-receptor \
    --input-file protein.pdb --chain A --output /path/to/receptor.json
```

Always specify `--output` with an absolute path or a path relative to the user's
project root, never a path relative to the skill directory.

--------------------------------------------------------------------------------

### 2. Predict Structure

Predicts protein structure from amino acid sequence using the ESMFold API.

```bash
uv run scripts/target_preparation_api.py predict-structure \
    --sequence "MKTLLILAVVAAALA..." --output /path/to/predicted.pdb

# From FASTA file
uv run scripts/target_preparation_api.py predict-structure \
    --fasta target.fasta --output /path/to/predicted.pdb
```

> **Note**: ESMFold predictions are fast (~30s) but less accurate than
> AlphaFold2 for multi-domain proteins. Always assess quality before proceeding.

--------------------------------------------------------------------------------

### 3. Generate AlphaFold2 HPC Script

Generates a SLURM submission script for AlphaFold2 or ColabFold on an HPC
cluster.

```bash
uv run scripts/target_preparation_api.py generate-af2-script \
    --fasta target.fasta --method colabfold \
    --gpu a100 --mem 64G --time 12:00:00 \
    --output /path/to/submit_af2.sh
```

--------------------------------------------------------------------------------

### 4. Assess Structure Quality

Reads pLDDT confidence metrics from a PDB file (B-factor column) and prints a
heuristic quality assessment.

```bash
uv run scripts/target_preparation_api.py assess-structure \
    --pdb predicted.pdb --output /path/to/quality.json
```

## Interpreting the Output

The `assess-structure` output contains quality metrics. Read it carefully and
synthesize the results for the user:

1.  **pLDDT Score**: Mean per-residue confidence.
    -   **> 70**: PROCEED with docking. Good confidence.
    -   **50–70**: REFINE. Consider using only high-confidence regions for
        binding site definition.
    -   **< 50**: STOP. Structure is likely unreliable for docking.
2.  **Low-confidence regions**: Report any contiguous stretches with pLDDT < 50.
    Advise the user to exclude these from binding site searches.
3.  **Recommendation**: Always state whether the structure is suitable for
    downstream docking and suggest refinement if needed.

## Workflow

1.  Check RCSB PDB or AlphaFold Database for an existing structure.
2.  If not found: predict with `predict-structure` (ESMFold) or
    `generate-af2-script` (HPC).
3.  Assess quality with `assess-structure` (pLDDT > 70 to proceed).
4.  Prepare receptor with `prepare-receptor` (clean, extract chain, PDBQT).
