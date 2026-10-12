# ADMET Screening Reference

**This file is generated from the check registry in**
`src/agskills/chem/admet.py`. Do not edit it by hand: regenerate it
with `python tools/generate_admet_reference.py`. It exists so the
documented thresholds cannot drift away from the ones the code
applies, which is what happened to the previous '14 checks' and
'34 checks' claims - the implementation behind them ran three.

## Tiers

| Tier | Checks | Gating | Advisory | Purpose |
|---|---:|---:|---:|---|
| `relaxed` (Tier 1 (relaxed)) | 14 | 13 | 1 | Post-docking triage. Selects a diverse seed set for generative chemistry; liabilities that a generator can optimise away are tolerated here. |
| `strict` (Tier 2 (strict)) | 24 | 22 | 2 | Final candidate selection before molecular dynamics and reporting. Published thresholds applied at face value, no structural alerts. |

A **gating** check decides the verdict. An **advisory** check is
computed and reported but does not. Advisory checks are published
risk *flags* rather than exclusion criteria: the Pfizer 3/75
criterion flags logP above 3 together with TPSA below 75, which
describes much of lead-like space, so gating on it would reject
most real medicinal chemistry.

## Tier `relaxed`: Tier 1 (relaxed)

Post-docking triage. Selects a diverse seed set for generative chemistry; liabilities that a generator can optimise away are tolerated here.

14 checks per compound (13 gating, 1 advisory).

| Check | Role | Category | Limit at this tier | Source |
|---|---|---|---|---|
| `mw` | gating | physicochemical | 100 to 700 | project policy |
| `logp` | gating | physicochemical | -3 to 7 | Wildman & Crippen, J Chem Inf Comput Sci 1999 |
| `tpsa` | gating | physicochemical | <= 200 | Ertl et al. J Med Chem 2000;43:3714-3717 |
| `hbd` | gating | physicochemical | <= 7 | project policy |
| `hba` | gating | physicochemical | <= 15 | project policy |
| `rotatable_bonds` | gating | physicochemical | <= 15 | project policy |
| `heavy_atoms` | gating | physicochemical | 8 to 90 | project policy |
| `formal_charge` | gating | physicochemical | -2 to 2 | project policy |
| `lipinski` | gating | drug-likeness rule | <= 2 violation(s) | project policy |
| `pains` | gating | structural alert | <= 1 | Baell & Holloway, J Med Chem 2010;53:2719-2740 |
| `brenk` | gating | structural alert | <= 3 | Brenk et al. ChemMedChem 2008;3:435-444 |
| `sa_score` | gating | synthesis | <= 8 | Ertl & Schuffenhauer, J Cheminform 2009;1:8 |
| `qed` | gating | composite | >= 0.1 | Bickerton et al. Nat Chem 2012;4:90-98 |
| `pfizer_3_75` | advisory | drug-likeness rule | rule default | project policy |

## Tier `strict`: Tier 2 (strict)

Final candidate selection before molecular dynamics and reporting. Published thresholds applied at face value, no structural alerts.

24 checks per compound (22 gating, 2 advisory).

| Check | Role | Category | Limit at this tier | Source |
|---|---|---|---|---|
| `mw` | gating | physicochemical | 150 to 500 | project policy |
| `logp` | gating | physicochemical | -0.4 to 5 | Wildman & Crippen, J Chem Inf Comput Sci 1999 |
| `tpsa` | gating | physicochemical | 20 to 140 | Ertl et al. J Med Chem 2000;43:3714-3717 |
| `hbd` | gating | physicochemical | <= 5 | project policy |
| `hba` | gating | physicochemical | <= 10 | project policy |
| `rotatable_bonds` | gating | physicochemical | <= 10 | project policy |
| `heavy_atoms` | gating | physicochemical | 20 to 70 | project policy |
| `molar_refractivity` | gating | physicochemical | 40 to 130 | project policy |
| `fraction_csp3` | gating | physicochemical | >= 0.1 | project policy |
| `formal_charge` | gating | physicochemical | -1 to 1 | project policy |
| `stereocenters` | gating | synthesis | <= 4 | project policy |
| `lipinski` | gating | drug-likeness rule | <= 1 violation(s) | project policy |
| `veber` | gating | drug-likeness rule | rule default | project policy |
| `egan` | gating | drug-likeness rule | rule default | project policy |
| `ghose` | gating | drug-likeness rule | rule default | project policy |
| `muegge` | gating | drug-likeness rule | rule default | project policy |
| `pfizer_3_75` | advisory | drug-likeness rule | rule default | project policy |
| `gsk_4_400` | gating | drug-likeness rule | rule default | project policy |
| `golden_triangle` | advisory | drug-likeness rule | rule default | project policy |
| `pains` | gating | structural alert | <= 0 | Baell & Holloway, J Med Chem 2010;53:2719-2740 |
| `brenk` | gating | structural alert | <= 0 | Brenk et al. ChemMedChem 2008;3:435-444 |
| `nih` | gating | structural alert | <= 0 | Doveston et al. / NIH filter set |
| `sa_score` | gating | synthesis | <= 6 | Ertl & Schuffenhauer, J Cheminform 2009;1:8 |
| `qed` | gating | composite | >= 0.3 | Bickerton et al. Nat Chem 2012;4:90-98 |

## Conventions

**Molecular weight** is the average mass (RDKit `MolWt`), which is
what Lipinski's rule means and what PubChem reports. The
monoisotopic mass is reported separately as `exact_mw` and is
never substituted: for aspirin those are 180.16 and 180.04 Da.

**Donors and acceptors** are reported under both conventions,
because they disagree and different rules were authored against
different ones. `hbd`/`hba` use RDKit's refined definitions, which
reproduce PubChem. `hbd_lipinski`/`hba_lipinski` use Lipinski's own
(OH + NH groups; all N and O atoms). The rule-of-five
implementation uses the Lipinski counts, as the 1997 paper
specifies. For metformin these differ substantially: 3/1 refined
against 5/5 Lipinski.

**logP** is the Wildman-Crippen estimate. Lipinski's paper states
the calculated-logP cutoff as 5; that is the form applied here.

**A rule always reports its authors' thresholds.** A tier changes
which rules must pass and how many violations are tolerated; it
never redefines what a violation is. A 650 Da compound reports a
molecular-weight violation at every tier.

## Accounting

Every compound ends in exactly one of `passed`, `failed` or
`errored`, and the three always sum to the number of inputs. A
check that could not be computed is reported as `skip` and
excluded from the verdict - it is never treated as a pass.
