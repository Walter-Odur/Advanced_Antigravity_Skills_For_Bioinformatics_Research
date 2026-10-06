---
name: target-preparation
description: >
  Use this skill when the user wants to prepare a protein target for
  structure-based drug discovery. Fetches structures from the RCSB PDB or
  the AlphaFold Database, predicts them from sequence with ESMFold, cleans
  receptors for docking, converts to PDBQT, assesses model confidence from
  pLDDT, and defines the docking search box. Do not use if the user only
  has a gene or protein name - ask for a PDB ID, a UniProt accession, or an
  amino-acid sequence first.
---

# Target Preparation

## Overview

Turns an identifier or a sequence into a docking-ready receptor plus the
search box that docking needs. Produces files and JSON reports; it never
mutates its inputs.

**Do NOT use when:**

- The user has only a gene or protein name. Ask them to look up the
  identifier at [RCSB PDB](https://www.rcsb.org) or
  [UniProt](https://www.uniprot.org), or use
  `compound-screening search-chembl-target`, which reports UniProt
  accessions alongside each hit.
- The user wants to screen or dock compounds - use `compound-screening`.
- The user wants molecular dynamics - use `md-simulation`.

## Prerequisites

| Requirement | Needed for |
|---|---|
| Python >= 3.10, `rdkit`, `biopython` | everything |
| Internet access | `prepare-receptor` with an ID, `predict-structure` |
| `meeko` + `scipy` + `gemmi`, or Open Babel, or the ADFR Suite | PDBQT conversion |

Install everything with `pip install -e ".[all]"` from the workshop root.
Without a PDBQT converter the skill still produces a cleaned PDB and says
so plainly rather than reporting success.

## Core Rules

1. **Always run the script; never hand-write the steps.** It handles PDB
   fetching, version resolution, chain extraction, solvent and heterogen
   removal, alternate-location selection, and conversion.
2. **`--output` is required** for every subcommand. Read the JSON it
   writes; the console summary is a convenience, not the result.
3. **Define the binding site before docking.** `prepare-receptor` writes
   `binding_site.json` automatically when the structure has a bound ligand.
   If it does not, use `define-site --residues`.
4. **Do not assess structure quality yourself.** Use `assess-structure` and
   report its verdict. It distinguishes a predicted model's pLDDT from a
   crystal structure's B-factors, which are not the same quantity and run
   in opposite directions.
5. If this skill is used, mention it in the output.

## Subcommands

| Subcommand | Purpose |
|---|---|
| `prepare-receptor` | Fetch, clean, convert to PDBQT, assess, and define the site |
| `define-site` | Produce the docking box JSON that `dock` consumes |
| `list-ligands` | List candidate ligands, labelling cofactors |
| `assess-structure` | Per-residue confidence and a PROCEED/REFINE/STOP verdict |
| `predict-structure` | Fold a sequence with ESMFold |
| `validate-sequence` | Check a sequence is foldable before submitting it |
| `generate-af2-script` | SLURM script for ColabFold or AlphaFold2 |

Invoke as either:

```bash
ag-target-preparation <subcommand> --output <file> [options]
python .agents/skills/target-preparation/scripts/target_preparation_api.py <subcommand> ...
```

## Quick Start

```bash
# An experimental structure, chain A, with everything derived in one step
ag-target-preparation prepare-receptor \
    --pdb-id 6HEZ --chain A --output work/receptor.json

# A predicted model, by UniProt accession
ag-target-preparation prepare-receptor \
    --uniprot-id P9WJA5 --output work/receptor.json

# A local file
ag-target-preparation prepare-receptor \
    --input-file protein.pdb --chain A --output work/receptor.json
```

`prepare-receptor` writes, alongside the JSON report:

- `<ID>.pdb` - the structure as downloaded
- `<ID>_clean.pdb` - solvent, additives and heterogens removed
- `<ID>_clean.pdbqt` - the docking receptor, when a converter is available
- `binding_site.json` - the search box, when a ligand was present

## Defining the Binding Site

Docking needs a box. Three ways to get one, in descending reliability:

```bash
# 1. From a co-crystallised ligand. The gold standard.
ag-target-preparation define-site \
    --pdb 6HEZ.pdb --from-ligand --output binding_site.json

# 2. From pocket residues, when the site is known from mutagenesis
#    or homology. Each selector is CHAIN:RESSEQ.
ag-target-preparation define-site \
    --pdb model.pdb --residues A:766 A:769 A:820 \
    --output binding_site.json

# 3. The whole protein, for blind docking. Honest but expensive and
#    much less accurate; the output is labelled low-confidence.
ag-target-preparation define-site \
    --pdb model.pdb --whole-structure --output binding_site.json
```

**Cofactors are not ligands.** `--from-ligand` prefers a drug-like ligand
over a known cofactor, then ranks by size. This matters: in 6HEZ
(*M. tuberculosis* DprE1) the FAD cofactor has 53 atoms and the inhibitor
0SK has 28, so ranking by size alone centres the box on the flavin site
rather than the inhibitor pocket. Check the `selected_ligand` and
`alternative_ligands` fields, and override with `--ligand-resname` if the
automatic choice is not the site you want.

Inspect what is available first with:

```bash
ag-target-preparation list-ligands --pdb 6HEZ.pdb --output ligands.json
```

## Interpreting `assess-structure`

The output carries a `verdict`, a `reason`, and specific `actions`.

**For a predicted model**, the B-factor column holds pLDDT on a 0-100
scale, where higher is better (AlphaFold DB convention):

| Mean pLDDT | Verdict | Meaning |
|---|---|---|
| >= 70, no long disordered region | `PROCEED` | Suitable for docking |
| 50-70, or a disordered loop of 5+ residues | `REFINE` | Usable in part; follow the listed actions |
| < 50, or over half the residues below 50 | `STOP` | Not reliable for docking |

Bands: >= 90 very high, 70-90 confident, 50-70 low, < 50 very low.

**For an experimental structure** the verdict is `NOT_APPLICABLE`. A
crystallographic B-factor is an atomic displacement parameter: it is
unbounded, **lower is better**, and the pLDDT bands do not apply. Judge
such a structure by its resolution and R-free instead. The output states
which case it detected in `b_factor_caveat`.

Always report low-confidence regions, and keep the docking box away from
them.

## Structure Prediction

```bash
# Quick, single-sequence prediction
ag-target-preparation predict-structure \
    --sequence "MKTLLILAVVAAALA..." --output predicted.pdb

# From a FASTA file
ag-target-preparation predict-structure \
    --fasta target.fasta --output predicted.pdb

# Check first, to avoid a wasted submission
ag-target-preparation validate-sequence \
    --fasta target.fasta --output check.json
```

The PDB goes to `--output`; its metadata, including the pLDDT assessment,
is written beside it as `.json`.

ESMFold predicts from a single sequence with **no multiple sequence
alignment**. It is fast but typically less accurate than AlphaFold2 on
multi-domain proteins and on targets with shallow evolutionary coverage.
The public endpoint accepts at most 400 residues. For anything longer, or
when accuracy matters:

```bash
ag-target-preparation generate-af2-script \
    --fasta target.fasta --method colabfold \
    --gpu a100 --mem 64G --time 12:00:00 \
    --output submit_af2.sh
# then: sbatch submit_af2.sh
```

## Cleaning Options

Removal is a scientific judgement, so it is explicit and itemised in the
report:

| Flag | Effect |
|---|---|
| `--chain A` | Keep one chain. Omit to keep all. |
| `--keep-ligand FAD NAD` | Retain named heterogens, e.g. a cofactor you intend to model |
| `--keep-waters` | Retain water. Off by default: structural waters must be chosen deliberately |
| `--remove-metals` | Also drop metal ions, which are kept by default because they are often catalytic |
| `--remove-hydrogens` | Strip hydrogens so the docking tool adds its own consistently |

Buffer counter-ions, cryoprotectants and crystallisation additives are
always removed. Only the primary alternate location is kept, so the
receptor has one unambiguous conformation.

## PDBQT Conversion Caveat

Meeko matches every residue against a chemical template. Residues with
missing side-chain atoms fail that match, which is routine in crystal
structures - 6HEZ has 16 such residues in chain A alone. The skill retries
in permissive mode and reports which residues were approximated, and
critically **whether any of them lie inside the docking box**. A
compromised residue outside the box is harmless; one inside it makes
affinities for that site less reliable. Check `pdbqt.warnings`.

## Workflow

1. Search for an existing structure. Prefer an experimental holo structure:
   it gives both the receptor and the box.
2. If none exists, `predict-structure` (fast) or `generate-af2-script`
   (accurate).
3. `assess-structure`, and act on the verdict.
4. `prepare-receptor`, which also writes `binding_site.json`.
5. Confirm the site is the pocket you intend, then hand the receptor and
   the site to `compound-screening dock`.
