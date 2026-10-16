---
name: target-preparation
description: >
  Prepare a protein target for structure-based drug discovery. Fetches
  structures from the RCSB PDB or AlphaFold DB, predicts them from sequence
  with ESMFold, cleans receptors, converts to PDBQT, assesses model
  confidence from pLDDT, and defines the docking search box. Use when the
  user has a PDB ID, a UniProt accession, or an amino-acid sequence. Do not
  use with only a gene or protein name.
---

# Target Preparation

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions to ensure
    `uv` is installed and on PATH. Script dependencies resolve automatically
    from each script's PEP 723 header; there is nothing to install by hand.
2.  **User Notification**: If `.licenses/target_preparation_LICENSE.txt` does
    not exist in the workspace root, prominently notify the user to review
    the terms at https://www.rcsb.org/pages/usage-policy and
    https://alphafold.ebi.ac.uk/, then create that file recording the
    notification text and a timestamp.

## Core Rules

-   [!IMPORTANT] **Use the Utility Script**: ALWAYS use
    `scripts/target_preparation_api.py`. NEVER fetch structures with `curl`
    or custom requests: the script enforces a shared rate limit, retries
    transient failures, and gives up within a bounded deadline.

-   **Output to File (Required)**: `--output` is required for every
    subcommand. Read the JSON it writes with `jq` or a short Python
    expression — do not print whole files into context.

-   **Define the site before docking**: `prepare-receptor` writes
    `binding_site.json` automatically when the structure has a bound ligand.
    Without one, use `define-site --residues`.

-   **Do not judge structure quality yourself**: use `assess-structure` and
    report its verdict. It distinguishes a predicted model's pLDDT from a
    crystal structure's B-factors, which run in opposite directions.

-   **Notification**: If this skill is used, mention it in the output.

## Utility Script

```bash
uv run scripts/target_preparation_api.py <subcommand> --output <file> [options]
```

| Subcommand | Purpose |
|---|---|
| `prepare-receptor` | Fetch, clean, convert to PDBQT, assess, define the site |
| `define-site` | Produce the docking box JSON that `dock` consumes |
| `list-ligands` | List candidate ligands, labelling cofactors |
| `assess-structure` | Per-residue confidence and a PROCEED/REFINE/STOP verdict |
| `predict-structure` | Fold a sequence with ESMFold |
| `validate-sequence` | Check a sequence is foldable before submitting it |
| `generate-af2-script` | SLURM script for ColabFold or AlphaFold2 |

--------------------------------------------------------------------------------

### 1. Prepare a Receptor

```bash
uv run scripts/target_preparation_api.py prepare-receptor \
    --pdb-id 6HEZ --chain A --output work/receptor.json
```

Accepts one of `--pdb-id`, `--uniprot-id` (AlphaFold) or `--input-file`.
Writes beside the JSON: `<ID>.pdb`, `<ID>_clean.pdb`,
`<ID>_clean.pdbqt` (when a converter is available), and
`binding_site.json` (when a ligand was present).

--------------------------------------------------------------------------------

### 2. Define the Binding Site

Docking needs a box. Three routes, most reliable first:

```bash
# From a co-crystallised ligand -- the gold standard
uv run scripts/target_preparation_api.py define-site \
    --pdb 6HEZ.pdb --from-ligand --output binding_site.json

# From pocket residues, as CHAIN:RESSEQ
uv run scripts/target_preparation_api.py define-site \
    --pdb model.pdb --residues A:766 A:769 A:820 --output binding_site.json

# Whole protein, for blind docking (labelled low-confidence)
uv run scripts/target_preparation_api.py define-site \
    --pdb model.pdb --whole-structure --output binding_site.json
```

> [!IMPORTANT] **Cofactors are not ligands.** `--from-ligand` prefers a
> drug-like ligand over a known cofactor before ranking by size. In 6HEZ
> (*M. tuberculosis* DprE1) the FAD cofactor has 53 atoms and the inhibitor
> 28, so size-only ranking centres the box on the flavin site instead of the
> inhibitor pocket. Check `selected_ligand`; override with
> `--ligand-resname`.

--------------------------------------------------------------------------------

### 3. Assess Structure Quality

```bash
uv run scripts/target_preparation_api.py assess-structure \
    --pdb predicted.pdb --output quality.json
```

Report the `verdict` and its `actions`:

| Verdict | Meaning |
|---|---|
| `PROCEED` | Mean pLDDT >= 70, no long disordered region |
| `REFINE` | Usable in part; follow the listed actions |
| `STOP` | Mean pLDDT < 50, or over half the residues below 50 |
| `NOT_APPLICABLE` | An experimental structure — judge it by resolution and R-free, not pLDDT |

The `scale` field says which of three B-factor conventions was detected.
See [references/confidence.md](references/confidence.md).

--------------------------------------------------------------------------------

### 4. Predict a Structure

```bash
uv run scripts/target_preparation_api.py validate-sequence \
    --fasta target.fasta --output check.json

uv run scripts/target_preparation_api.py predict-structure \
    --fasta target.fasta --output predicted.pdb
```

ESMFold is single-sequence (no MSA), fast, and capped at 400 residues. For
longer chains or better accuracy, generate a cluster script:

```bash
uv run scripts/target_preparation_api.py generate-af2-script \
    --fasta target.fasta --method colabfold --output submit_af2.sh
```

One script plus one guide (`submit_af2_SUBMIT.md`): the tools the cluster needs with probe and install commands, files to stage, the `sbatch`/`squeue`/`sacct` commands, and what each failure means. Point the user at it. `--module NAME` is repeatable and makes the job load the cluster's own build.

## Common Options

-   `--output FILE`: **Required.** Output path.
-   `--chain ID`: Keep one chain. Omit to keep all.
-   `--keep-ligand RESNAME ...`: Retain named heterogens, e.g. a cofactor.
-   `--keep-waters`, `--remove-metals`, `--remove-hydrogens`: Cleaning
    overrides. Metals are kept by default because they are often catalytic;
    buffer ions and cryoprotectants are always removed.
-   `--detail compact|standard|full`: How much to write. `compact` is the
    default and carries the verdict plus the numbers behind it; `full` adds
    every region and descriptor and costs several times the context.

## Common Mistakes

1.  **Docking into the cofactor site.** Confirm `selected_ligand.is_cofactor`
    is `false` before docking, or pass `--ligand-resname`.
2.  **Applying a pLDDT verdict to a crystal structure.** If `scale` is
    `bfactor`, the confidence bands do not apply; lower B-factors are better
    and the scale is unbounded.
3.  **Ignoring `pdbqt.warnings`.** Crystal structures routinely have
    incomplete side chains. The report names any approximated residue that
    lies **inside** the docking box; one there makes affinities for that
    site unreliable, whereas one elsewhere is harmless.

## Reference

-   [references/questions.md](references/questions.md) — **read this first when the user's request is vague, phrased as a question, or mentions a target, organism or compound class not covered above.** Indexed by question, with worked examples across human, viral, bacterial, fungal and protozoan targets.
-   [references/best_practices.md](references/best_practices.md) — structural biology practice for SBDD: resolution and R-free, which structure to choose, waters and metals, oligomeric state, protonation and flips, box size, and the limits of predicted models.
-   [references/confidence.md](references/confidence.md) — pLDDT bands, the
    three B-factor scales, and how the verdict is reached.
-   [references/citation.bib](references/citation.bib) — cite these for any
    structure or prediction used.

## Workflow

1.  Prefer an experimental holo structure: it supplies both the receptor and
    the box.
2.  If none exists, `validate-sequence` then `predict-structure`, or
    `generate-af2-script` for a cluster run.
3.  `assess-structure` and act on the verdict.
4.  `prepare-receptor`, which also writes `binding_site.json`.
5.  Confirm the site is the intended pocket, then hand the receptor and site
    to the `compound-screening` skill's `dock` subcommand.
