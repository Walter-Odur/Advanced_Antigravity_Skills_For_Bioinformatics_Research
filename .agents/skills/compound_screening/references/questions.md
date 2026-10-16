# Questions this skill answers

Find your question, run the command. The commands are chemistry-general:
nothing here is specific to a target class, a disease area or an organism.

## Start here

| You want | Run |
|---|---|
| Known actives for a target | `search-chembl-target` then `query-chembl` |
| Natural products | `query-coconut` |
| A big purchasable library | `query-zinc` |
| To know if a compound is drug-like | `admet-filter` |
| To know what the filters actually check | `describe-tiers` |
| Binding affinities | `dock` (needs a receptor from `target-preparation`) |
| A full ADMET profile | `admet-predict` |

--------------------------------------------------------------------------------

## Finding compounds

**Q1. I have a target name. How do I get compounds known to hit it?**

```bash
uv run scripts/compound_screening_api.py search-chembl-target \
    --query "DprE1" --output target.json
```

This returns candidates **with organism and UniProt accession**, which
matters more than it sounds: a name usually matches several species.
Searching `DprE1` returns both *M. smegmatis* (CHEMBL3797019) and
*M. tuberculosis* (CHEMBL3804751). Pick the organism you mean, then:

```bash
uv run scripts/compound_screening_api.py query-chembl \
    --target-id CHEMBL3804751 --pchembl-min 7 --output actives.json
```

**Q2. What is pChEMBL and what threshold should I use?**

pChEMBL is −log₁₀ of a molar potency (IC₅₀, Kᵢ, Kd or EC₅₀), so 7 means
100 nM and 9 means 1 nM. Higher is more potent. Rough guidance:

| `--pchembl-min` | Keeps | Use for |
|---|---|---|
| 5 | ≤ 10 µM | hit lists, fragments, weak-but-real chemistry |
| 6 | ≤ 1 µM | a broad active set |
| 7 | ≤ 100 nM | a solid set for seeding generative design |
| 8 | ≤ 10 nM | leads and tool compounds only |

Start at 6 or 7. If you get too few compounds, lower it rather than
abandoning the target.

**Q3. Why does it return fewer molecules than there are activity records?**

Because one molecule measured several times must be **one row**. Measured on
CHEMBL3804751: 34 activity records collapse to 32 distinct molecules, and
the aggregate is recorded rather than an arbitrary pick:

```json
{"molecule_chembl_id": "CHEMBL3330226", "pchembl_median": 8.0,
 "pchembl_max": 8.0, "pchembl_min": 7.92, "n_measurements": 3,
 "n_assays": 3, "standard_types": ["IC50"], "units": ["nM"]}
```

Without this, a compound assayed three times appears three times and is
triaged three times, inflating every hit rate you compute afterwards.

**Q4. The compounds I got back look wrong for my target.**

Then the target ID is wrong, not the query. Sanity-check the chemotypes
against what you know of the literature: a DprE1 query should return
benzothiazinones and nitroaromatics; an EGFR query should return
4-anilinoquinazolines. A set of unrelated scaffolds means you have the wrong
CHEMBL ID, or the target is a protein family rather than a single protein.

**Q5. I want natural products.**

```bash
uv run scripts/compound_screening_api.py query-coconut \
    --query "flavonoid" --output np.json
```

COCONUT aggregates open natural-product collections [16]. Natural products
systematically violate Lipinski-style rules and that is expected, not a
defect — screen them at `--strictness relaxed` and read Q11.

**Q6. I want a large purchasable library for virtual screening.**

```bash
uv run scripts/compound_screening_api.py query-zinc \
    --logp-range -1 3 --mw-range 250 350 --output download.sh
```

This writes a download script rather than fetching gigabytes. ZINC encodes
its tranches as four letters — molecular weight, logP, reactivity,
purchasability [14, 15] — and the script expands your ranges into the right
set. The generated URLs are real: `https://files.docking.org/2D/CB/CBAA.smi`
returns a live file with `smiles zinc_id` columns.

For ultra-large screening, read [7] before committing compute; the practical
guides matter more than the library size.

--------------------------------------------------------------------------------

## Screening and filtering

**Q7. Is my compound drug-like?**

```bash
uv run scripts/compound_screening_api.py admet-filter \
    --smiles "CC(=O)OC1=CC=CC=C1C(=O)O" --names aspirin \
    --strictness relaxed --output screen.json
```

Every check is reported individually with its threshold and its citation.

**Q8. Which strictness should I use?**

| Tier | Checks | Use when |
|---|---|---|
| `relaxed` | 14 (13 gating + 1 advisory) | triaging a screen, choosing seeds for generative design, fragments, natural products |
| `strict` | 24 (22 gating + 2 advisory) | selecting final candidates for MD and reporting |

`describe-tiers` prints exactly what each performs, read from the registry
rather than from documentation, so the counts cannot drift.

**Q9. Aspirin and caffeine FAIL the strict tier. Is the filter broken?**

No — and this is the most misread output in the skill. Measured: aspirin,
caffeine and imatinib all fail `strict`, while all three pass `relaxed`.

Aspirin fails on `Heavy atoms: 13 < 20`. It is a 13-heavy-atom fragment,
below the lead-like window, which is a statement about **where it sits in
chemical space**, not about whether it is a drug. The strict tier is a
design-space filter for lead optimisation [8, 9]; read `passed 0` as "none
of these are lead-like", never as "none of these are drugs".

If you are screening marketed drugs, repurposing candidates, fragments or
natural products, use `relaxed`.

**Q10. Why does the Lipinski check say 500 Da when I asked for a 700 Da
limit?**

Because a published rule always reports **its authors' thresholds**. A tier
changes which rules must pass and how many violations are tolerated; it
never redefines a violation. Measured on rifampicin:

```
relaxed | mw       | 100 Da to 700 Da  | 822.953 Da > 700 Da
relaxed | lipinski | <= 2 violation(s) | MW 822.953 Da > 500 Da; HBD 6 > 5 ...
strict  | mw       | 150 Da to 500 Da  | 822.953 Da > 500 Da
strict  | lipinski | <= 1 violation(s) | MW 822.953 Da > 500 Da; HBD 6 > 5 ...
```

The tier's own window moves; the Lipinski text is identical. A tool that
rewrote Lipinski's 500 Da to its own 700 Da would report rifampicin as
rule-of-five compliant, which no medicinal chemist would accept.

**Q11. My compounds are big and greasy — natural products, macrocycles,
PROTACs. Are they hopeless?**

No. There is a well-characterised oral drug space beyond the rule of 5,
occupied by macrocycles, many natural products and most of the
protein–protein-interaction inhibitors [10, 11]. Use `relaxed`, read the
individual violations rather than the overall verdict, and judge against
comparable chemistry rather than against oral small molecules.

The same caution applies in reverse for antibacterials: compounds that
cross the Gram-negative envelope are far more polar and more charged than
typical oral drugs, so oral-drug property filters select *against* them
[12, 13].

**Q12. All five of my known actives were flagged. Should I drop them?**

Check whether the flag is **advisory**. Measured on five real EGFR
inhibitors from ChEMBL: all five pass at 13 of 14 checks, and all five are
flagged for Pfizer 3/75 (logP > 3 **and** TPSA < 75) [5].

```
passed: 5 | compounds_with_advisories: 5
egfr1 13/14 | ['Pfizer 3/75 toxicity risk: logP 4.15 > 3 and TPSA 56.3 < 75']
```

Pfizer 3/75 is a published *risk flag*, not an exclusion criterion. Gating
on it would reject 100% of a validated active series. It is reported so you
can weigh it; the Golden Triangle is advisory for the same reason.

**Q13. What is PAINS, and why did my compound fail it?**

Pan-assay interference compounds: substructures that produce apparent
activity by mechanisms other than specific binding — redox cycling, metal
chelation, covalent reaction, aggregation, fluorescence [3]. They are a
reason to be suspicious of an assay result, not proof a compound is
worthless [4].

Measured on benzylidene rhodanine, the most-cited PAINS scaffold, three
independent catalogues flag it:

```
 - PAINS alerts: ene_rhod_A(235)
 - Brenk alerts: Michael_acceptor_1; Thiocarbonyl_group
 - NIH alerts: rhodanine
```

`relaxed` tolerates ≤1 PAINS and ≤3 Brenk deliberately — a generative model
can optimise an alert away, so excluding alerts at the seed stage discards
usable chemistry. `strict` tolerates none.

**Q14. What is QED? What is an SA score? Which numbers should I care about?**

| Metric | Range | Reading |
|---|---|---|
| **QED** | 0–1, higher better | Quantitative estimate of drug-likeness; a weighted desirability across eight properties [1] |
| **SA score** | 1–10, lower better | Synthetic accessibility; ~1–4 routine, >6 hard [2] |
| **TPSA** | Å² | Polar surface area; ≲140 for oral absorption, ≲90 for CNS [6] |
| **logP** | — | Lipophilicity. See Q15 |
| **Ligand efficiency** | — | Potency per heavy atom; the right way to compare molecules of different size [17] |

A single number should never decide anything. QED around 0.7 is good for an
oral small molecule and meaningless for a fragment.

**Q15. Your logP does not match PubChem's.**

It is not meant to. This skill reports Wildman–Crippen `MolLogP`; PubChem
reports XLogP3. They are **different estimators**, and measured across ten
marketed drugs they differ by up to 1.3 units — atorvastatin 6.31 against
5.0, imatinib 4.59 against 3.5.

A threshold is only meaningful against the estimator it was authored for.
Expecting the two to agree is the mistake, not the gap. The same applies to
HBA counts: three of the ten disagree with PubChem, and all three carry a
carboxylic acid, whose hydroxyl oxygen PubChem counts as an acceptor and
RDKit's refined definition does not. Both conventions are reported, and the
Lipinski rule uses Lipinski's own.

**Q16. One bad SMILES in my file — will the whole run die?**

No. Invalid structures are counted as `errored`, separately from `failed`,
and the run continues with exit code 0:

```
  total 4  passed 1  failed 1  errored 2
  [ERROR] broken: invalid SMILES: RDKit could not parse the structure
```

`errored` is kept distinct because an unparseable string is not a compound
that failed screening, and collapsing the two corrupts every hit rate.

--------------------------------------------------------------------------------

## Docking

**Q17. How do I dock something?**

You need a receptor and a box from the `target-preparation` skill:

```bash
uv run scripts/compound_screening_api.py dock \
    --receptor work/6HEZ_clean.pdbqt --site work/binding_site.json \
    --smiles "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC" --output dock.json
```

**Q18. It says AutoDock Vina is not available.**

Exit code 4, with an install route. Vina is not pip-installable on Windows
because the Python module needs Boost headers:

```bash
conda install -c conda-forge vina
# or a binary from https://github.com/ccsb-scripps/AutoDock-Vina/releases
```

Alternatively generate a cluster script with
`ag-md-simulation generate-hpc-script --job-type vina-screen --array 1-100`.

**Q19. My score is −9.4 kcal/mol. Is that good? Is it the binding energy?**

It is a reasonable score for a drug-like ligand in a real pocket, and it is
**not** a binding free energy. Practical reading:

- Typical drug-like range: roughly −6 to −12.
- Differences under about **1 kcal/mol do not rank compounds** [18, 19].
- Scoring functions predict poses better than affinities, and enrichment
  better than either [18].
- A positive score, or one beyond about −15, usually means the box is
  misplaced or the receptor is not protonated — not an extraordinary hit.

**Q20. All my compounds score about the same.**

Usually the box is too large, so everything finds a mediocre surface
minimum, or the pocket is shallow. Tighten the box
([target-preparation](../../target_preparation/references/best_practices.md)
section 6) and confirm it is on the intended ligand, not a cofactor.

**Q21. How do I know my docking setup is any good?**

Validate it before trusting it:

1.  **Re-dock the cognate ligand** from the crystal structure and measure
    RMSD to the crystal pose. Under 2 Å is the usual pass mark.
2.  **Screen known actives against decoys** and measure enrichment.
    DUD-E is the common source; LIT-PCBA exists because DUD-E-style decoys
    can be separated on properties alone rather than on binding [20, 21].
3.  Expect worse performance against a predicted model than an experimental
    one [22].

--------------------------------------------------------------------------------

## ADMET

**Q22. What is the difference between `admet-filter` and `admet-predict`?**

`admet-filter` applies published physicochemical rules and structural alerts
— fast, transparent, every threshold cited.

`admet-predict` is the same battery plus a fuller descriptor profile, and it
says plainly what it could not do:

```
ADMET prediction (offline RDKit engine)
  Tier 2 (strict): 24 checks
  note: no public ADMETlab REST endpoint is reachable; see
        'online_prediction' in the output.
```

There is no free, stable REST API for the published ADMET model servers, so
nothing is fetched and nothing is faked. Read `online_prediction` in the
JSON for what it tried.

Neither command is a measured assay. Both are property-based estimates;
absorption, metabolism and toxicity need experiment.

**Q23. Can it predict toxicity / hERG / CYP inhibition?**

Only insofar as published physicochemical surrogates go — Pfizer 3/75 for
in-vivo toxicity risk [5], structural alerts for reactivity [3, 4]. Treat
them as flags that direct experiments, not as predictions.

--------------------------------------------------------------------------------

## Output and cost

**Q24. The JSON is huge.**

`--detail compact` is the default and is roughly 13× smaller than `full`.
Measured on 300 compounds:

| Detail | JSON | Per compound | stdout |
|---|---|---|---|
| `compact` | 187,156 B | 624 B | 15 lines |
| `standard` | 329,537 B | 1,098 B | 15 lines |
| `full` | 2,475,036 B | 8,250 B | 15 lines |

Verdicts are identical at every level; only the per-compound payload
changes. stdout stays at 15 lines whether you screen 50 compounds or 300.

**Q25. How do I get just the passing compounds?**

```bash
uv run scripts/compound_screening_api.py admet-filter \
    --smiles-file library.smi --strictness strict --passing-only \
    --output hits.json
```

Or read them out:

```bash
python -c "
import json; d=json.load(open('screen.json'))
print([c['name'] for c in d['compounds'] if c['status']=='pass'])"
```

--------------------------------------------------------------------------------

## When something goes wrong

| Symptom | Cause | Fix |
|---|---|---|
| Exit 4 on `dock` | Vina missing | Q18 |
| Exit 3 | Invalid input data | Check SMILES or identifiers |
| Exit 5 | ChEMBL/COCONUT unreachable | Retry; the client already backs off |
| 0 compounds from ChEMBL | pChEMBL threshold too high, or wrong target | Q2, Q4 |
| Everything fails `strict` | Wrong tier for your chemistry | Q9, Q11 |
| `invalid SMILES` | Unclosed ring, bad element symbol | Q16 |
| Mode 1 missing from a Vina log | Old parser | Fixed; modes report bare `0` RMSD |

--------------------------------------------------------------------------------

## References

Full entries in [citation.bib](citation.bib). Thresholds and their sources
are also tabulated in [admet_thresholds.md](admet_thresholds.md).

1.  Bickerton GR et al. *Nat Chem* 2012;4:90–98 — QED.
2.  Ertl P, Schuffenhauer A. *J Cheminform* 2009;1:8 — synthetic
    accessibility.
3.  Baell JB, Holloway GA. *J Med Chem* 2010;53:2719–2740 — PAINS.
4.  Brenk R et al. *ChemMedChem* 2008;3:435–444 — unwanted substructures.
5.  Hughes JD et al. *Bioorg Med Chem Lett* 2008;18:4872–4875 — Pfizer 3/75.
6.  Veber DF et al. *J Med Chem* 2002;45:2615–2623 — TPSA and rotatable
    bonds.
7.  Bender BJ et al. *Nat Protoc* 2021;16:4799–4832 — large-scale docking.
8.  Leeson PD, Springthorpe B. *Nat Rev Drug Discov* 2007;6:881–890 —
    drug-like concepts in decision-making.
9.  Hann MM, Keserü GM. *Nat Rev Drug Discov* 2012;11:355–365 — the sweet
    spot.
10. Doak BC, Over B, Giordanetto F, Kihlberg J. *Chem Biol*
    2014;21:1115–1142 — oral druggable space beyond the rule of 5.
11. Doak BC et al. *J Med Chem* 2016;59:2312–2327 — how bRo5 drugs bind.
12. O'Shea R, Moser HE. *J Med Chem* 2008;51:2871–2878 — physicochemical
    properties of antibacterials.
13. Richter MF et al. *Nature* 2017;545:299–304 — compound accumulation
    rules for Gram-negative bacteria.
14. Sterling T, Irwin JJ. *J Chem Inf Model* 2015;55:2324–2337 — ZINC 15.
15. Irwin JJ et al. *J Chem Inf Model* 2020;60:6065–6073 — ZINC20.
16. Sorokina M et al. *J Cheminform* 2021;13:2 — COCONUT.
17. Hopkins AL, Groom CR, Alex A. *Drug Discov Today* 2004;9:430–431 —
    ligand efficiency.
18. Warren GL et al. *J Med Chem* 2006;49:5912–5931 — docking assessed.
19. Eberhardt J et al. *J Chem Inf Model* 2021;61:3891–3898 — Vina 1.2.0.
20. Mysinger MM et al. *J Med Chem* 2012;55:6582–6594 — DUD-E.
21. Tran-Nguyen VK, Jacquemard C, Rognan D. *J Chem Inf Model*
    2020;60:4263–4273 — LIT-PCBA.
22. Scardino V, Di Filippo JI, Cavasotto CN. *iScience* 2023;26:105920 —
    predicted structures in virtual screening.
23. Zdrazil B et al. *Nucleic Acids Res* 2024;52:D1180–D1192 — ChEMBL.
