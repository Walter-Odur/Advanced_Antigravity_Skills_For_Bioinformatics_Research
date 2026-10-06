---
name: compound-screening
description: >
  Use this skill when the user wants to search compound libraries, dock
  molecules against a prepared receptor, or screen compounds on ADMET
  properties. Handles ChEMBL target and bioactivity queries, COCONUT
  natural-product search, ZINC bulk download scripts, AutoDock Vina
  docking, and two-tier ADMET screening (relaxed for triage and generative
  seeding, strict for final candidates). Do not use for generative
  chemistry or de novo design - use the compound-synthesis skill instead.
---

# Compound Screening

## Overview

Compound sourcing, molecular docking, and ADMET screening. The screening
engine is offline and deterministic, so a result is reproducible and
citable.

**Do NOT use when:**

- The user wants to generate novel molecules - use `compound-synthesis`.
- The user wants to set up molecular dynamics - use `md-simulation`.
- The user wants to prepare a receptor or define a binding site - use
  `target-preparation`.

## Prerequisites

| Requirement | Needed for |
|---|---|
| Python >= 3.10, `rdkit` | ADMET screening, all local work |
| Internet access | `search-chembl-target`, `query-chembl`, `query-coconut` |
| AutoDock Vina, plus `meeko`/`scipy`/`gemmi` | `dock` |

`admet-filter`, `describe-tiers` and `query-zinc` are fully offline.

## Core Rules

1. **Always run the script.** Never query these APIs with `curl` or ad-hoc
   requests: the client enforces a shared, cross-process rate limit,
   retries transient failures with backoff, honours `Retry-After`, and
   gives up within a bounded deadline instead of hanging.
2. **`--output` is required** for every subcommand. Read the JSON; the
   console summary is a convenience.
3. **Use the right tier, and report which one you used.** Relaxed before
   generative chemistry, strict for final candidates. The output carries
   `strictness`, `tier` and `checks_per_compound`.
4. **A docking box is required.** Get `binding_site.json` from
   `target-preparation define-site`. Never invent coordinates.
5. **Report advisories, not just failures.** A compound can pass while
   carrying a flagged liability.
6. If this skill is used, mention it in the output.

## Subcommands

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

Invoke as either:

```bash
ag-compound-screening <subcommand> --output <file> [options]
python .agents/skills/compound-screening/scripts/compound_screening_api.py <subcommand> ...
```

## Two-Tier ADMET Screening

The pipeline has two decision points, and they want different thresholds.

**Tier 1, `--strictness relaxed` (14 checks).** Applied after docking, to
choose seeds for generative chemistry. Wide windows, tolerant of alerts.
The goal is a diverse seed set: a liability the generator can optimise away
is not a reason to discard a scaffold.

**Tier 2, `--strictness strict` (24 checks, the default).** Applied to
final candidates before molecular dynamics and reporting. Published
thresholds at face value, no structural alerts.

```bash
# Triage, and seed selection
ag-compound-screening admet-filter \
    --csv docked_hits.csv --smiles-col SMILES \
    --strictness relaxed --output triage.json

# Final candidates
ag-compound-screening admet-filter \
    --smiles "CC(=O)Oc1ccccc1C(=O)O" --names aspirin \
    --strictness strict --output candidates.json

# What exactly does each tier check?
ag-compound-screening describe-tiers --output tiers.json
```

The check counts above are **derived from the registry**, not asserted in
prose, and `describe-tiers` prints the live definition. See
[reference/admet_thresholds.md](reference/admet_thresholds.md), which is
generated from the same source.

### Gating versus advisory checks

A **gating** check decides the verdict. An **advisory** check is computed
and reported but does not. This distinction is scientific, not cosmetic:
the Pfizer 3/75 criterion flags logP above 3 together with TPSA below 75,
which describes much of lead-like space - ibuprofen and the entire
4-anilinoquinazoline kinase-inhibitor series trip it. It was published as a
signal to look harder at toxicity, not as an exclusion filter. Read each
compound's `advisories` field and judge whether the flag matters for the
series in hand.

### Reading the output

| Field | Meaning |
|---|---|
| `overall_pass` | Passed every gating check |
| `failed_checks` | The gating failures: why it did not pass |
| `advisories` | Flags raised without blocking |
| `checks_passed` / `_failed` / `_skipped` | Per-compound tally; sums to the tier's check count |
| `properties` | Every computed descriptor |
| `checks` | Each check with its value, limit, role and citation |

`passed + failed + errored` always equals `total_compounds`. A compound
whose structure could not be parsed is `errored`, never silently dropped.
A check that could not be computed is `skip` and is excluded from the
verdict - it is never treated as a pass.

Marketed drugs routinely fail the strict tier, and that is correct
behaviour rather than a defect: aspirin and caffeine fall below its size
floor, and penicillin G fails only on the beta-lactam structural alert.
Report *which* check failed so the user can judge relevance.

## Finding Compounds

```bash
# 1. Find the target. The output includes UniProt accessions, which feed
#    'target-preparation prepare-receptor --uniprot-id'.
ag-compound-screening search-chembl-target \
    --query "DprE1" --organism "tuberculosis" --limit 10 \
    --output targets.json

# 2. Retrieve known actives. pChEMBL 6 is 1 uM, 7 is 100 nM, 9 is 1 nM.
ag-compound-screening query-chembl \
    --target-id CHEMBL3804751 --pchembl-min 7 --limit 25 \
    --output actives.json

# 3. Natural products, including the African and marine collections
#    (AfroDB, ANPDB, NANPDB, SANCDB, CMNPD).
ag-compound-screening query-coconut \
    --query "antimycobacterial" --limit 50 --output natural.json
```

`query-chembl` **aggregates per molecule**. ChEMBL often holds several
measurements of the same compound against the same target, so a ranking
taken straight from the activity list can be one molecule repeated.
Results are ranked on `pchembl_median`, which is robust to a single
optimistic measurement, and `n_measurements` shows how well supported each
value is. Paging is followed, so a request for more than one page returns a
complete set rather than a silently truncated one.

## Bulk Libraries from ZINC

```bash
ag-compound-screening query-zinc \
    --subset drug-like --output download_zinc.sh
bash download_zinc.sh
```

Subsets: `drug-like`, `lead-like`, `fragment-like`, `all-purchasable`.
Override with `--mw-range` and `--logp-range`, and select reactivity and
availability with `--reactivity` and `--purchasability`.

Pass `--verify` to HEAD-check the tranche URLs before writing the script.
This is worth doing: a tranche that does not exist is served as an **HTML
page with HTTP 200**, so a script that only checks whether a file arrived
will happily save web pages as `.smi` files. The generated script inspects
the content of every download and rejects HTML, skips files already
present so an interrupted bulk download resumes, and exits non-zero if
nothing was retrieved.

## Docking

```bash
# Prepare the receptor and the box first
ag-target-preparation prepare-receptor \
    --pdb-id 6HEZ --chain A --output work/receptor.json

ag-compound-screening dock \
    --receptor work/6HEZ_clean.pdbqt \
    --site work/binding_site.json \
    --csv seeds.csv --smiles-col SMILES \
    --exhaustiveness 32 --seed 42 \
    --output docking.json
```

`--site` is the reliable route. `--center X Y Z` with `--box-size X Y Z` is
available when the coordinates are known independently.

Docking is **seeded** (`--seed`, default 42) so a run is reproducible.
Conformer generation and the Vina search are both stochastic; without a
fixed seed the same input gives different affinities on every run.

For more than a few hundred compounds, generate a job array rather than
docking locally:

```bash
ag-md-simulation generate-hpc-script \
    --job-type vina-screen --receptor receptor.pdbqt \
    --ligand-dir ligands/ --site-config binding_site.json \
    --array 1-100 --ncpus 16 --output submit_screen.sh
```

### Interpreting affinities

Affinity is in kcal/mol and **more negative is stronger**. Treat the
ranking as triage:

- Differences under about 1 kcal/mol are not meaningful.
- A Vina score is a scoring function, not a binding free energy.
- Confirm the top poses visually and by molecular dynamics.
- Check whether the receptor had approximated residues inside the box
  (`pdbqt.warnings` from `prepare-receptor`); if so, affinities for that
  site are less reliable.

Every compound appears in `ranked_results`, including failures, each with a
`status` (`ok`, `preparation_failed`, `docking_failed`, `timeout`,
`no_pose`) and an explanation.

## A Note on Online ADMET Prediction

`admet-predict` runs the offline engine and says so. No public ADMETlab 3.0
REST endpoint is currently reachable: the path the previous implementation
used, `/server/api/aio`, returns HTTP 404, as do seven other candidate
paths, while the web interface itself responds normally. Rather than
present a failure as a result, this command returns the deterministic
offline battery and records the situation under `online_prediction`.

If you have access to a prediction service, set `AGSKILLS_ADMET_URL` and
the command will post to it and merge the response.

For reproducible work the offline engine is preferable anyway: it is
deterministic, versioned with the package, and every threshold carries a
citation.

## Workflow

1. `search-chembl-target`, then `query-chembl` for known actives.
2. `query-coconut` or `query-zinc` to broaden the set.
3. Prepare the receptor and the box with `target-preparation`.
4. `dock` a small set locally, or generate a job array for a large one.
5. `admet-filter --strictness relaxed` to triage and to pick generative
   seeds.
6. Hand the seeds to `compound-synthesis`.
7. `admet-filter --strictness strict` on what comes back.
8. Hand the survivors to `md-simulation`.
