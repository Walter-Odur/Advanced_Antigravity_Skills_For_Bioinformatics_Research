---
name: drug-discovery-wizard
description: >
  Orchestrator skill for the complete structure-based drug discovery pipeline.
  Coordinates four sub-skills: target-preparation, compound-screening,
  compound-synthesis (REINVENT 4), and md-simulation. Use this skill when the
  user wants to run the full end-to-end pipeline.
---

# Drug Discovery Wizard (Orchestrator)

## Overview

This is the **orchestrator** for the full structure-based drug discovery
pipeline. It coordinates four independent sub-skills, each of which can also
be used standalone.

**Do NOT use when:**

-   The user only wants one pipeline step (e.g., just docking, just ADMET, just
    MD setup) — use the individual sub-skill directly.
-   The user wants wet lab protocol design, clinical trial design (use the
    `clinical-trials-database` skill), or retrosynthetic route planning
    (suggest ASKCOS or IBM RXN).

## Sub-Skills

| Skill | Purpose | When to Use |
|---|---|---|
| `target-preparation` | Structure prediction, receptor prep | User has a PDB ID, UniProt ID, or sequence |
| `compound-screening` | Library search, docking, ADMET | User wants to find/filter/dock compounds |
| `compound-synthesis` | REINVENT 4 generative chemistry | User wants de novo design or optimization |
| `md-simulation` | GROMACS MD setup, HPC scripts | User wants to simulate protein-ligand dynamics |

## Core Rules

-   **Always use sub-skill utility scripts** — never write raw API calls or
    GROMACS commands by hand.
-   **`--output` is required** for every subcommand across all sub-skills.
-   **Two-tier ADMET**: Use `--strictness relaxed` post-docking (broad), then
    `--strictness strict` post-REINVENT (rigorous).
-   **HPC for heavy compute**: REINVENT 4, GROMACS MD, and batch Vina screening
    should always generate HPC submission scripts.
-   If this skill is used, ensure this is mentioned in the output.

## Full Pipeline (Two-Tier Strategy)

```
1. TARGET PREPARATION (target-preparation skill)
   ├── Check PDB / AlphaFold DB for existing structure
   ├── If missing: predict-structure (ESMFold) or generate-af2-script (HPC)
   ├── assess-structure (pLDDT > 70 to proceed)
   └── prepare-receptor (clean, extract chain, PDBQT)

2. COMPOUND SOURCING (compound-screening skill)
   ├── search-chembl-target / query-chembl (known actives)
   ├── query-coconut (natural products: ANPDB, CMNPD, NANPDB, AfroDB, SANCDB)
   └── query-zinc (purchasable drug-like compounds)

3. MOLECULAR DOCKING (compound-screening skill)
   ├── dock (small set, local)
   └── generate-hpc-script --job-type vina-screen (large set, HPC)

4. TIER 1 ADMET — RELAXED (compound-screening skill)
   └── admet-filter --strictness relaxed  →  broad hits (14 checks)

5. GENERATIVE CHEMISTRY (compound-synthesis skill)
   ├── generate-reinvent-config
   └── generate-hpc-script --job-type reinvent

6. TIER 2 ADMET — STRICT (compound-screening skill)
   └── admet-filter --strictness strict  →  final candidates (34 checks)

7. MD SIMULATION (md-simulation skill)
   ├── gmxapi-setup (workflow + MDP files)
   └── generate-hpc-script --job-type gromacs-md

8. RANKED CANDIDATES
   └── Binding affinity + ADMET profile + MD stability
```

## Quick Start Example (Anti-TB, DprE1)

```bash
# Step 1: Prepare receptor
uv run .agents/skills/target-preparation/scripts/target_preparation_api.py \
    prepare-receptor --pdb-id 6HEZ --chain A --output receptor.json

# Step 2: Search libraries
uv run .agents/skills/compound-screening/scripts/compound_screening_api.py \
    search-chembl-target --query "DprE1" --limit 10 --output targets.json

uv run .agents/skills/compound-screening/scripts/compound_screening_api.py \
    query-coconut --query "antimycobacterial" --limit 50 --output coconut.json

# Step 3: Relaxed ADMET (Tier 1)
uv run .agents/skills/compound-screening/scripts/compound_screening_api.py \
    admet-filter --smiles "CCO" --strictness relaxed --output seeds.json

# Step 4: REINVENT config (compound-synthesis skill)

# Step 5: Strict ADMET (Tier 2)
uv run .agents/skills/compound-screening/scripts/compound_screening_api.py \
    admet-filter --smiles "generated" --strictness strict --output candidates.json

# Step 6: MD setup
uv run .agents/skills/md-simulation/scripts/md_simulation_api.py \
    gmxapi-setup --pdb protein_clean.pdb --ff charmm36 --production-ns 100 \
    --output md_workflow/run_md.py
```
