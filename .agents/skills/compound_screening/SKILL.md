---
name: compound-screening
description: >
  Search compound libraries, dock molecules, and screen them on ADMET
  properties. Queries ChEMBL for targets and bioactivities, COCONUT for
  natural products, and generates verified ZINC bulk downloads; docks with
  AutoDock Vina; screens with a two-tier ADMET battery. Use when the user
  asks about compounds, actives, docking, binding affinity, drug-likeness or
  ADMET filtering. Do not use for de novo design — use compound-synthesis.
---

# Compound Screening

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions.
    Dependencies resolve from the script's PEP 723 header automatically.
2.  **User Notification**: If `.licenses/compound_screening_LICENSE.txt` does
    not exist in the workspace root, prominently notify the user to review
    the terms at
    https://chembl.gitbook.io/chembl-interface-documentation/about,
    https://coconut.naturalproducts.net/about and
    https://zinc.docking.org/, then create that file recording the
    notification text and a timestamp.
3.  **AutoDock Vina** (`dock` only): not pip-installable on Windows. Install
    with `conda install -c conda-forge vina`. Everything else is offline or
    needs only internet.

## Core Rules

-   [!IMPORTANT] **Use the Utility Script**: ALWAYS use
    `scripts/compound_screening_api.py`. NEVER query these APIs with `curl`
    or custom requests: the script enforces a shared cross-process rate
    limit, retries transient failures with backoff, honours `Retry-After`,
    and gives up within a bounded deadline instead of hanging.

-   **Output to File (Required)**: `--output` is required for every
    subcommand. Read the JSON with `jq` or a short Python expression. A
    1,000-compound screen is a large file — never print one into context.

-   **Use the right tier, and report which**: `relaxed` before generative
    chemistry, `strict` for final candidates. The output carries
    `strictness`, `tier` and `checks_per_compound`.

-   **A docking box is required**: obtain `binding_site.json` from the
    `target-preparation` skill's `define-site`. Never invent coordinates.

-   **Report advisories, not just failures**: a compound can pass while
    carrying a flagged liability.

-   **Notification**: If this skill is used, mention it in the output.

## Utility Script

```bash
uv run scripts/compound_screening_api.py <subcommand> --output <file> [options]
```

| Subcommand | Purpose | Offline? |
|---|---|:---:|
| `admet-filter` | Two-tier ADMET screening | yes |
| `describe-tiers` | Print exactly which checks each tier performs | yes |
| `admet-predict` | Full profile; see the note on online services | yes |
| `search-chembl-target` | Find a ChEMBL target ID by name | no |
| `query-chembl` | Known actives for a target, with potency | no |
| `query-coconut` | Natural-product search | no |
| `query-zinc` | Generate a verified bulk download script | yes |
| `dock` | AutoDock Vina docking | needs Vina |

--------------------------------------------------------------------------------

### 1. ADMET Screening (Two Tiers)

The pipeline makes two different decisions, and they want different
thresholds.

```bash
# Tier 1 -- post-docking triage, and seed selection for generative chemistry
uv run scripts/compound_screening_api.py admet-filter \
    --csv docked_hits.csv --smiles-col SMILES \
    --strictness relaxed --output triage.json

# Tier 2 -- final candidates (the default)
uv run scripts/compound_screening_api.py admet-filter \
    --smiles "CC(=O)Oc1ccccc1C(=O)O" --names aspirin \
    --strictness strict --output candidates.json

# What exactly does each tier check?
uv run scripts/compound_screening_api.py describe-tiers --output tiers.json
```

`relaxed` runs 14 checks with wide windows: the goal is a diverse seed set,
and a liability the generator can optimise away is no reason to discard a
scaffold. `strict` runs 24 checks at published thresholds with no structural
alerts. Counts are derived from the check registry, so they cannot drift
from what is performed.

**Gating versus advisory.** A gating check decides the verdict; an advisory
check is reported but does not reject. Published *risk flags* are advisory:
the Pfizer 3/75 criterion catches logP > 3 with TPSA < 75, which describes
much of lead-like space — ibuprofen and the whole 4-anilinoquinazoline
kinase series trip it. Read each compound's `advisories`.

Reading the output:

| Field | Meaning |
|---|---|
| `overall_pass` | Passed every gating check |
| `failed_checks` | The gating failures — why it did not pass |
| `advisories` | Flags raised without blocking |
| `properties` | Headline descriptors (all of them at `--detail standard`) |

`passed + failed + errored` always equals `total_compounds`. An unparseable
structure is `errored`, never dropped. A check that could not run is
`skip`, never counted as a pass.

> Marketed drugs routinely fail `strict`, and that is correct: aspirin falls
> below its size floor and penicillin G fails only on the beta-lactam alert.
> Report **which** check failed so the user can judge relevance.

See [references/admet_thresholds.md](references/admet_thresholds.md).

--------------------------------------------------------------------------------

### 2. Find Compounds

```bash
# Targets. The output includes UniProt accessions, which feed
# target-preparation's prepare-receptor --uniprot-id.
uv run scripts/compound_screening_api.py search-chembl-target \
    --query "DprE1" --organism "tuberculosis" --limit 10 --output targets.json

# Known actives. pChEMBL 6 is 1 uM, 7 is 100 nM, 9 is 1 nM.
uv run scripts/compound_screening_api.py query-chembl \
    --target-id CHEMBL3804751 --pchembl-min 7 --limit 25 --output actives.json

# Natural products, including AfroDB, ANPDB, NANPDB, SANCDB and CMNPD.
uv run scripts/compound_screening_api.py query-coconut \
    --query "antimycobacterial" --limit 50 --output natural.json
```

`query-chembl` **aggregates per molecule** and follows paging. ChEMBL holds
several measurements of the same compound against one target, so an
unaggregated "top 20" can be one molecule twenty times. Results rank on
`pchembl_median`, which resists a single optimistic measurement;
`n_measurements` shows how well supported each value is.

--------------------------------------------------------------------------------

### 3. Bulk Libraries from ZINC

```bash
uv run scripts/compound_screening_api.py query-zinc \
    --subset drug-like --verify --output download_zinc.sh
bash download_zinc.sh
```

Subsets: `drug-like`, `lead-like`, `fragment-like`, `all-purchasable`.
Override with `--mw-range` and `--logp-range`.

> [!IMPORTANT] A ZINC tranche that does not exist is served as **HTML with
> HTTP 200**. The generated script inspects content and rejects HTML, skips
> files already present so an interrupted download resumes, and exits
> non-zero if nothing arrived. `--verify` HEAD-checks the URLs first.

--------------------------------------------------------------------------------

### 4. Docking

```bash
uv run scripts/compound_screening_api.py dock \
    --receptor work/6HEZ_clean.pdbqt \
    --site work/binding_site.json \
    --csv seeds.csv --smiles-col SMILES \
    --exhaustiveness 32 --seed 42 --output docking.json
```

`--site` is the reliable route; `--center X Y Z` with `--box-size X Y Z`
works when coordinates are known independently. Runs are **seeded** so they
reproduce: conformer generation and the Vina search are both stochastic.

For more than a few hundred compounds, generate a job array with the
`md-simulation` skill's `generate-hpc-script --job-type vina-screen`.

**Interpreting affinities.** More negative is stronger, in kcal/mol. Treat
the ranking as triage: differences under about 1 kcal/mol are not
meaningful, and a Vina score is a scoring function, not a binding free
energy. Confirm top poses visually and by MD.

--------------------------------------------------------------------------------

### 5. Online ADMET Prediction

`admet-predict` runs the offline engine and says so. No public ADMETlab 3.0
REST endpoint is currently reachable: `/server/api/aio` returns HTTP 404, as
do seven other candidate paths, while the web interface responds normally.
Rather than present a failure as a result, the command returns the
deterministic offline battery and records the situation under
`online_prediction`. Set `AGSKILLS_ADMET_URL` if you have access to a
service.

## Common Options

-   `--output FILE`: **Required.** Output path.
-   `--smiles S [S ...]` | `--csv FILE` | `--smiles-file FILE`: compound
    input. `--names` must match `--smiles` one-to-one or it is an error.
-   `--strictness relaxed|strict`: screening tier (default: `strict`).
-   `--limit N`: maximum results for a database query.
-   `--passing-only`: write only the compounds that passed.
-   `--show N`: how many compounds to summarise on stdout (default: 12).
-   `--detail compact|standard|full`: how much to write. `compact` is the
    default; `full` adds every individual check and costs roughly ten times
    the context.

## Common Mistakes

1.  **Printing a whole screening file into context.** At `--detail full` a
    200-compound screen is ~280,000 tokens. Read fields with `jq`, or keep
    the default `compact`.
2.  **Treating a sub-kcal/mol docking difference as real.** It is not.
3.  **Applying `strict` to generative seeds.** It discards chemistry the
    generator could have repaired. Use `relaxed` first.

## Reference

-   [references/questions.md](references/questions.md) — **read this first when the user's request is vague, phrased as a question, or mentions a target, organism or compound class not covered above.** Indexed by question: library choice, pChEMBL thresholds, which tier to believe, what a docking score means.
-   [references/admet_thresholds.md](references/admet_thresholds.md) —
    generated from the check registry, so it cannot drift.
-   [references/citation.bib](references/citation.bib) — cite these for any
    screening or docking result.

## Workflow

1.  `search-chembl-target`, then `query-chembl` for known actives.
2.  `query-coconut` or `query-zinc` to broaden the set.
3.  Prepare the receptor and box with the `target-preparation` skill.
4.  `dock` a small set locally, or generate a job array for a large one.
5.  `admet-filter --strictness relaxed` to triage and pick generative seeds.
6.  Hand the seeds to `compound-synthesis`.
7.  `admet-filter --strictness strict` on what comes back.
8.  Hand survivors to `md-simulation`.
