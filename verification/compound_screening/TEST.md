# TEST.md — compound-screening

A verification protocol for the `compound-screening` skill. Work through it
to satisfy yourself the skill is accurate and functional on your machine.

**Not shipped with the skill.** This lives under `verification/` so the
bundle under `.agents/skills/` stays lean.

Every expected value below was **measured** by running the command, not
copied from the documentation. Where the skill disagrees with an external
reference, the disagreement is recorded and explained rather than hidden.

## Setup

```bash
cd <workshop root>
export S=.agents/skills/compound_screening/scripts/compound_screening_api.py
export W=/tmp/cs_test && mkdir -p $W
```

`uv run $S ...` builds its own environment from the launcher's PEP 723
header. `python $S ...` needs an interpreter that already has `rdkit` and
`polite_http`.

| Need | For |
|---|---|
| `uv` | all cases |
| Internet | T8, T9, T10 |
| AutoDock Vina | T11 only; T12 is the case for *not* having it |

--------------------------------------------------------------------------------

## T1 — the tiers are what the skill says they are

```bash
uv run $S describe-tiers --output $W/tiers.json
```

**Expected:** `Tier 1 (relaxed): 14 checks`, `Tier 2 (strict): 24 checks`,
each check listed with its threshold.

**Pass criteria**

1.  The printed count equals the number of checks listed. The counts come
    from the registry, not from a constant in the documentation, so a check
    added without updating the docs cannot go unnoticed.
2.  The gating/advisory split is explicit:

```bash
python -c "
import json; d=json.load(open('$W/tiers.json'))
for n,t in d['tiers'].items():
    print(n, t['n_gating_checks'], t['n_advisory_checks'],
          [c['check'] for c in t['checks'] if c['role']!='gating'])"
```

**Measured:**

```
relaxed 13 1 ['pfizer_3_75']
strict  22 2 ['pfizer_3_75', 'golden_triangle']
```

T5 shows why that split matters.

--------------------------------------------------------------------------------

## T2 — descriptors against an independent reference

The most important accuracy case in the skill.
`tests/data/pubchem_reference.json` holds PubChem's own computed values for
ten marketed drugs, fetched from PubChem rather than produced by this code.

```bash
python -c "
import json, sys; sys.path.insert(0,'src')
from rdkit import Chem
from agskills.chem.descriptors import compute_properties
ref=json.load(open('tests/data/pubchem_reference.json'))
for k,v in ref.items():
    m=Chem.MolFromSmiles(v['SMILES']); p=compute_properties(m, smiles=v['SMILES'])
    print('%-12s%9.2f%9.2f%7.2f%7.1f%7.1f%7.1f%5d%5d%5d%5d%5d%5d' % (
        k, p.mw, float(v['MolecularWeight']), p.logp, v['XLogP'],
        p.tpsa, float(v['TPSA']), p.hbd, v['HBondDonorCount'],
        p.hba, v['HBondAcceptorCount'], p.rotatable_bonds, v['RotatableBondCount']))"
```

**Measured** (skill value then PubChem's, for each descriptor):

```
drug               MW      ref   logP  XLogP   TPSA    ref  HBD  ref  HBA  ref  rot  ref
aspirin        180.16   180.16   1.31    1.2   63.6   63.6    1    1    3    4    2    3
caffeine       194.19   194.19  -1.03   -0.1   61.8   58.4    0    0    3    3    0    0
imatinib       493.62   493.60   4.59    3.5   86.3   86.3    2    2    7    7    7    7
bedaquiline    555.52   555.50   7.13    7.2   45.6   45.6    1    1    4    4    8    8
isoniazid      137.14   137.14  -0.31   -0.7   68.0   68.0    2    2    3    3    1    1
rifampicin     822.95   822.90   4.34    4.9  220.1  220.0    6    6   15   15    4    5
atorvastatin   558.65   558.60   6.31    5.0  111.8  112.0    4    4    4    6   12   12
ibuprofen      206.28   206.28   3.07    3.5   37.3   37.3    1    1    1    2    4    4
metformin      129.17   129.16  -1.24   -1.3   91.5   91.5    3    3    1    1    0    2
ethambutol     204.31   204.31  -0.29   -0.1   64.5   64.5    4    4    4    4    9    9
```

**Pass criteria**

| Descriptor | Agreement | Verdict |
|---|---|---|
| **MW** | 10/10 within 0.05 Da | Must match. A mismatch means the monoisotopic mass was substituted for the average mass |
| **HBD** | 10/10 exact | Must match |
| **TPSA** | 8/10 exact, 2 within 0.2 | Must match except caffeine |
| **HBA** | 7/10 exact | Convention difference, below |
| **Rotatable bonds** | 7/10 exact | Convention difference |
| **logP** | differs by up to 1.3 | **Different method**, below |

Four of these are to be understood, not treated as noise:

-   **logP is not comparable.** The skill reports Wildman–Crippen
    `MolLogP`; PubChem reports XLogP3. They are different estimators, and
    atorvastatin (6.31 vs 5.0) and imatinib (4.59 vs 3.5) show how far
    apart they land. A threshold is only meaningful against the estimator
    it was authored for. Expecting these two columns to agree is the
    mistake, not the 1.3-unit gap.
-   **HBA disagrees on carboxylic acids.** The three mismatches — aspirin
    (3 vs 4), ibuprofen (1 vs 2), atorvastatin (4 vs 6) — all carry a
    `-COOH`, whose hydroxyl oxygen PubChem counts as an acceptor and
    RDKit's `CalcNumHBA` does not.
-   **Caffeine's TPSA** is 61.8 here against PubChem's 58.4, a known
    implementation difference over its aromatic nitrogens.
-   **Rotatable bonds** uses RDKit's strict definition, which excludes
    amide C–N bonds. Metformin (0 vs 2) is the clearest case.

Both HBD/HBA conventions must be present and must differ where they should:

```bash
python -c "
import sys; sys.path.insert(0,'src')
from rdkit import Chem
from agskills.chem.descriptors import compute_properties
s='CN(C)C(=N)N=C(N)N'
p=compute_properties(Chem.MolFromSmiles(s), smiles=s)
print('refined ', p.hbd, p.hba)
print('lipinski', p.hbd_lipinski, p.hba_lipinski)"
```

**Measured** for metformin:

```
refined  3 1
lipinski 5 5
```

If these ever print the same numbers, one convention has been aliased to
the other. `rules.py` applies the Lipinski counts, as the 1997 paper
specifies.

--------------------------------------------------------------------------------

## T3 — the strict tier is a design filter, not a drug detector

```bash
uv run $S admet-filter \
  --smiles "CC(=O)OC1=CC=CC=C1C(=O)O" \
           "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" \
           "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5" \
  --names aspirin caffeine imatinib --strictness strict --output $W/drugs.json
```

**Measured:** `total 3  passed 0  failed 3`.

Aspirin fails on `Heavy atoms: 13 < 20`. Caffeine fails on
`Lipophilicity (clogP): -1.0293 < -0.4`. Imatinib — a marketed kinase
inhibitor — fails on `Molar refractivity: 146.894 > 130` and GSK 4/400.

**This is the correct result and the most misreadable output in the skill.**
The strict tier is the lead-like window used for final candidate selection.
Aspirin is below it because it is a 13-heavy-atom fragment, not because
aspirin is not a drug. Read `passed 0` as "none of these sit in the
lead-like design space", never as "none of these are drugs".

The same three at relaxed:

```bash
uv run $S admet-filter --smiles <the same three> \
    --names aspirin caffeine imatinib --strictness relaxed --output $W/rel.json
```

**Measured:** `total 3  passed 3  failed 0`, each at 14/14.

**Pass criterion:** 0/3 strict and 3/3 relaxed. If both tiers agree on this
set, the tiers are not actually distinct.

--------------------------------------------------------------------------------

## T4 — a published rule always reports its authors' thresholds

The convention most implementations get wrong. A tier changes *which* rules
must pass and *how many* violations are tolerated. It never rewrites what a
published rule considers a violation.

```bash
R="C[C@H]1/C=C/C=C(\C(=O)NC2=C(C(=C3C(=C2O)C(=C(C4=C3C(=O)[C@](O4)(O/C=C/[C@@H]([C@H]([C@H]([C@@H]([C@@H]([C@@H]([C@H]1O)C)O)C)OC(=O)C)C)OC)C)C)O)O)/C=N/N5CCN(CC5)C)/C"
for t in relaxed strict; do
  uv run $S admet-filter --smiles "$R" --names rifampicin \
      --strictness $t --detail full --output $W/rif_$t.json
done
python -c "
import json
for t in ['relaxed','strict']:
    c=json.load(open('$W/rif_'+t+'.json'))['compounds'][0]
    for ch in c['checks']:
        if ch['check'] in ('mw','lipinski'):
            print(t, '|', ch['check'], '|', ch['policy'], '|', ch['detail'])"
```

**Measured:**

```
relaxed | mw       | 100 Da to 700 Da  | 822.953 Da > 700 Da
relaxed | lipinski | <= 2 violation(s) | MW 822.953 Da > 500 Da; HBD 6 > 5; HBA 16 > 10
strict  | mw       | 150 Da to 500 Da  | 822.953 Da > 500 Da
strict  | lipinski | <= 1 violation(s) | MW 822.953 Da > 500 Da; HBD 6 > 5; HBA 16 > 10
```

**Pass criteria**

1.  The tier's own `mw` window moves: 700 Da at relaxed, 500 Da at strict.
2.  The Lipinski `detail` string is **byte-identical** at both tiers. Only
    the tolerance changes, from ≤2 violations to ≤1.
3.  Lipinski reports `HBA 16` — the Lipinski count, not the refined count
    of 15 — as T2 requires.

A tier that rewrote Lipinski's 500 Da to its own 700 Da would report
rifampicin as rule-of-five compliant, which no medicinal chemist would
accept.

--------------------------------------------------------------------------------

## T5 — advisory checks do not reject compounds

The case that made the distinction necessary. These are five real
4-anilinoquinazoline EGFR inhibitors, as ChEMBL returns them for CHEMBL203
at pChEMBL ≥ 8.

```bash
uv run $S admet-filter \
  --smiles "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC" \
           "CN(C)c1cc2ncnc(Nc3cccc(Br)c3)c2cn1" \
           "Nc1ccc2c(Nc3cccc(Br)c3)ncnc2c1" \
           "CNc1cc2ncnc(Nc3cccc(Br)c3)c2cn1" \
           "Cc1ccc(Nc2ncnc3cc(OC)c(OC)cc23)cc1" \
  --names egfr1 egfr2 egfr3 egfr4 egfr5 \
  --strictness relaxed --output $W/egfr.json
```

**Measured:** `total 5  passed 5  failed 0`, every compound at **13/14**:

```
  [PASS] egfr5  MW 295.34 logP   3.7 QED 0.795 SA 1.86 (13/14)
  [PASS] egfr2  MW 344.22 logP   3.6 QED 0.787 SA 2.38 (13/14)
  [PASS] egfr4  MW 330.19 logP  3.57 QED 0.769 SA 2.34 (13/14)
  [PASS] egfr1  MW 360.21 logP  4.15 QED 0.757 SA 1.99 (13/14)
  [PASS] egfr3  MW 315.17 logP  3.72 QED 0.708 SA 2.05 (13/14)
```

```bash
python -c "
import json; d=json.load(open('$W/egfr.json'))
print(d['passed'], d['compounds_with_advisories'])
for c in d['compounds']: print(c['name'], c['advisories'])"
```

**Measured:** `5 5`, and each compound carries

```
Pfizer 3/75 toxicity risk: logP 3.70 > 3 and TPSA 56.3 A^2 < 75
  (elevated in-vivo toxicity risk)
```

**Pass criteria**

1.  `passed: 5`. These are validated actives against a real target; a
    screen that rejects them is not usable for triage.
2.  `compounds_with_advisories: 5` — every one is flagged.
3.  Each scores **13/14, not 14/14**. The flag is reported, not suppressed.

Pfizer 3/75 (logP > 3 **and** TPSA < 75) is a published *risk flag*. Were
it gating, this screen would reject 100% of a validated active series —
which is exactly what it did before the gating/advisory split existed.

Ordering is QED descending (0.795 → 0.708), so the most drug-like compound
is the first line of output.

--------------------------------------------------------------------------------

## T6 — structural alerts, and three catalogues agreeing

Benzylidene rhodanine is the most-cited PAINS scaffold in the literature.

```bash
uv run $S admet-filter --smiles "O=C1NC(=S)SC1=Cc1ccccc1" --names rhodanine \
    --strictness strict --output $W/rhod.json
python -c "
import json
for f in json.load(open('$W/rhod.json'))['compounds'][0]['failed_checks']:
    print(' -', f)"
```

**Measured — 6 of 24 checks fail:**

```
 - Heavy atoms: 14 < 20
 - Fraction of sp3 carbons: 0 < 0.1
 - Ghose drug-like range: heavy atoms 14 outside 20-70
 - PAINS alerts: ene_rhod_A(235)
 - Brenk alerts: Michael_acceptor_1; Thiocarbonyl_group
 - NIH alerts: rhodanine
```

**Pass criteria**

1.  **All three catalogues fire independently** — PAINS names the filter
    `ene_rhod_A(235)`, Brenk names the Michael acceptor and the
    thiocarbonyl, NIH names the scaffold. They are separate published
    catalogues; their agreement here is what shows each is really being
    consulted rather than one aliased to another.
2.  At **relaxed** the same compound **passes 14/14**, with the alerts
    still recorded in the check details. Relaxed tolerates ≤1 PAINS and
    ≤3 Brenk deliberately: a generative model can optimise an alert away,
    so excluding alerts at the seed stage discards usable chemistry.
    Strict tolerates none.
3.  A clean compound reports nothing. For
    `COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC` the PAINS and Brenk details are
    empty — a filter that fires on everything is as useless as one that
    never fires.

--------------------------------------------------------------------------------

## T7 — a bad structure does not take the batch down

```bash
uv run $S admet-filter \
  --smiles "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC" "c1ccccc" "[Xx]" "CCO" \
  --names good broken notanelement ethanol \
  --strictness relaxed --output $W/mixed.json; echo "rc=$?"
```

**Measured:**

```
  total 4  passed 1  failed 1  errored 2
  [PASS]  good     MW 360.21 logP  4.15 QED 0.757 SA 1.99 (13/14)
  [FAIL]  ethanol  MW  46.07 logP   0.0 QED 0.407 SA 1.98 (12/14)
  [ERROR] broken: invalid SMILES: RDKit could not parse the structure
  [ERROR] notanelement: invalid SMILES: RDKit could not parse the structure
rc=0
```

**Pass criteria**

1.  `errored` is counted **separately** from `failed`. An unparseable
    structure is not a compound that failed screening, and collapsing the
    two corrupts every hit rate computed from the output.
2.  The two valid compounds are still screened.
3.  Exit code **0** — the run succeeded; individual compounds did not. A
    screen of 10,000 library compounds must not abort on row 7.

--------------------------------------------------------------------------------

## T8 — ChEMBL target search, live

```bash
uv run $S search-chembl-target --query "DprE1" --output $W/tgt.json
```

**Measured:**

```
ChEMBL target search: 'DprE1' -> 2 target(s)
  CHEMBL3797019  Mycobacterium smegmatis (strai  Decaprenylphosphoryl-beta-D-ribose oxida
      UniProt: A0R607
  CHEMBL3804751  Mycobacterium tuberculosis (st  Decaprenylphosphoryl-beta-D-ribose oxida
      UniProt: P9WJF1
```

**Pass criteria**

1.  Both organisms are returned, not just the top hit. They differ only by
    species, and picking the wrong one silently yields actives against the
    wrong organism.
2.  The organism and the UniProt accession are both shown, which is what
    makes the choice possible without a second query.

--------------------------------------------------------------------------------

## T9 — ChEMBL actives, and measurement aggregation

```bash
uv run $S query-chembl --target-id CHEMBL3804751 --pchembl-min 7 \
    --output $W/actives.json
```

**Measured:**

```
ChEMBL activities for CHEMBL3804751 (pChEMBL >= 7.0)
  examined 34 activity records over 1 page(s)
  32 distinct molecules, returning 25
  CHEMBL3330226  pChEMBL  8.00 (n=3) O=c1nc(N2CCN(CC3CCCCC3)CC2)sc2c([N+](=O)[O-])c
  CHEMBL5180878  pChEMBL  8.00 (n=1) O=C(NC(=O)c1ccsc1NC(=O)c1nc2ccccc2s1)OCC1CC1
  ...
```

**Pass criteria**

1.  **34 records collapse to 32 molecules.** One molecule measured in
    several assays must be one row. Otherwise a compound assayed three
    times appears three times and is triaged three times.
2.  The aggregate is recorded rather than one arbitrary value:

```bash
python -c "
import json; d=json.load(open('$W/actives.json'))
c=[x for x in d['compounds'] if x['n_measurements']>1][0]
print({k:c[k] for k in ('molecule_chembl_id','pchembl_median','pchembl_max',
       'pchembl_min','n_measurements','n_assays','standard_types','units')})"
```

**Measured:**

```
{'molecule_chembl_id': 'CHEMBL3330226', 'pchembl_median': 8.0,
 'pchembl_max': 8.0, 'pchembl_min': 7.92, 'n_measurements': 3,
 'n_assays': 3, 'standard_types': ['IC50'], 'units': ['nM']}
```

The median comes with its spread and its provenance, so a consensus of
three IC50 measurements is distinguishable from a single one.

3.  **A chemistry check on the result.** The compounds returned should be
    recognisable DprE1 chemotypes. Benzothiazinones and nitroaromatics
    dominate the measured output, which is what the DprE1 literature
    describes. A set of unrelated scaffolds would mean the target id, not
    the query, is wrong.

--------------------------------------------------------------------------------

## T10 — ZINC tranche encoding, checked against the live server

The tranche letters are an encoding. Getting one wrong downloads a
different slice of chemical space and reports no error.

```bash
uv run $S query-zinc --logp-range -1 3 --mw-range 250 350 --output $W/zinc.sh
```

**Measured:**

```
ZINC tranche selection: drug-like
  MW 250-350 Da (letters CDE)
  logP -1-3 (letters BCDEF)
  15 tranches, format smi
  example URL: https://files.docking.org/2D/CB/CBAA.smi
```

**Pass criteria**

1.  3 MW letters × 5 logP letters = **15** tranches. If the count is not
    the product, a range boundary is being dropped.
2.  Each code is four letters — MW, logP, reactivity, purchasability — and
    the defaults are recorded in words, not as bare letters:
    `"reactivity": "A (anodyne - the most conservative filter, reactive
    groups excluded)"`.
3.  `bash -n $W/zinc.sh` exits 0, and the file has no CR bytes:
    `python -c "print(open('$W/zinc.sh','rb').read().count(b'\r'))"` prints
    `0`. Measured: 0. Check bytes rather than `grep -c $'\r'`, which is
    unreliable under Git Bash — it can match every line in a clean file.
4.  **The generated URL is real.** This is what shows the encoding is
    correct rather than merely self-consistent:

```bash
curl -sI https://files.docking.org/2D/CB/CBAA.smi | head -1
curl -s  https://files.docking.org/2D/CB/CBAA.smi | head -2
```

**Measured:**

```
HTTP/1.1 200 OK
smiles zinc_id
CN(C)CC1(O)CCN(C(=O)CN2CCCCC2=O)CC1 23582000
```

A real 200, the real `smiles zinc_id` header, real ZINC ids.

--------------------------------------------------------------------------------

## T11 — docking, with Vina present

Skip to T12 if Vina is not installed. That is the commoner case, and it is
also a test.

```bash
uv run $S dock --receptor $W/6HEZ_clean.pdbqt --site $W/binding_site.json \
    --smiles "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC" --output $W/dock.json
```

**Pass criteria**

1.  The box in `dock.json` is the box from `binding_site.json`, unchanged.
2.  `affinity_kcal_mol` is negative and between about −5 and −12 for a
    drug-like ligand in a real pocket. A positive score, or one beyond
    −15, means the box is misplaced or the receptor is not protonated.
3.  The output states that a Vina score is **not** a binding free energy
    and that differences under about 1 kcal/mol do not rank compounds.

**The log parser is verified without Vina**, against authentic captured
output in `tests/data/vina_output.txt`:

```bash
python -c "
import sys; sys.path.insert(0,'src')
from agskills.docking.vina import parse_vina_log
r=parse_vina_log(open('tests/data/vina_output.txt').read())
print('modes:', len(r)); print(r[0])"
```

**Measured:**

```
modes: 9
PoseScore(mode=1, affinity_kcal_mol=-9.428, rmsd_lb=0.0, rmsd_ub=0.0)
```

**Pass criterion:** 9 modes, and **mode 1 is present**. Real Vina v1.2.5
prints mode 1's RMSD columns as a bare `0`, not `0.000`:

```
   1       -9.428          0          0
   2       -9.106      1.892      2.771
```

A parser whose pattern requires a decimal point silently drops mode 1 — the
best pose, and the only one most pipelines read. This was a real defect,
found only by parsing authentic output instead of a handwritten sample.

--------------------------------------------------------------------------------

## T12 — docking without Vina fails usefully

```bash
uv run $S dock --receptor $W/6HEZ_clean.pdbqt --site $W/binding_site.json \
    --smiles "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC" --output $W/dock.json \
    >/dev/null 2>&1; echo $?
```

**Measured:**

```
error: Required dependency not available: AutoDock Vina
hint: Install it with: conda install -c conda-forge vina, or download a
      binary from https://github.com/ccsb-scripps/AutoDock-Vina/releases
rc=4
```

**Pass criteria**

1.  Exit code **4**, the documented code for a missing dependency, which is
    distinguishable from bad input (3) and a service failure (5).
2.  The hint names a command that works on this platform. Vina is not
    pip-installable on Windows, because the Python module needs Boost
    headers, so a hint reading `pip install vina` would be wrong here.
3.  No traceback, and no empty `dock.json` left behind.

| Code | Meaning |
|---|---|
| 0 | Success |
| 2 | Bad usage |
| 3 | Invalid input data |
| 4 | A required dependency or binary is missing |
| 5 | A remote service failed |
| 6 | A required local file does not exist |

--------------------------------------------------------------------------------

## T13 — output size stays bounded

The JSON is read back by whoever called the skill, so its size is a running
cost paid on every screen.

```bash
python -c "
a=['COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC','CN(C)c1cc2ncnc(Nc3cccc(Br)c3)c2cn1',
   'Nc1ccc2c(Nc3cccc(Br)c3)ncnc2c1','CNc1cc2ncnc(Nc3cccc(Br)c3)c2cn1',
   'Cc1ccc(Nc2ncnc3cc(OC)c(OC)cc23)cc1']
open('$W/panel300.smi','w').write(''.join('%s c%d\n' % (a[i%5],i) for i in range(300)))"

for d in compact standard full; do
  uv run $S admet-filter --smiles-file $W/panel300.smi --strictness strict \
      --detail $d --output $W/p300_$d.json > $W/so_$d.txt 2>&1
  echo "$d: $(stat -c%s $W/p300_$d.json) bytes, stdout $(wc -l < $W/so_$d.txt) lines"
done
```

**Measured, 300 compounds:**

| Detail | JSON | Per compound | ~Tokens/compound | stdout |
|---|---:|---:|---:|---:|
| `compact` (default) | 187,156 B | 624 B | 156 | 15 lines |
| `standard` | 329,537 B | 1,098 B | 275 | 15 lines |
| `full` | 2,475,036 B | 8,250 B | 2,063 | 15 lines |

**Pass criteria**

1.  `compact` is roughly **13× smaller** than `full` — a 92% saving on
    detail that is rarely read.
2.  **stdout is 15 lines at every detail level and every input size.** A
    50-compound screen also prints 15. The console summary does not grow
    with the input, so a 10,000-compound screen cannot flood the caller's
    context.
3.  **The verdicts are identical across all three levels.** Cheaper output
    must never mean a different answer:

```bash
python -c "
import json
v=[tuple((c['name'],c['status']) for c in
   json.load(open('$W/p300_'+d+'.json'))['compounds'])
   for d in ('compact','standard','full')]
print('verdicts identical:', v[0]==v[1]==v[2])"
# verdicts identical: True
```

4.  `compact` carries `tier_definition_command` in place of the full tier
    definition, and lists only the first few failed checks with a count.
    The rest is one command away rather than inlined 300 times.

--------------------------------------------------------------------------------

## Results

| Case | What it checks | Result |
|---|---|---|
| T1 | tier counts, gating/advisory split | pass — 14 (13+1), 24 (22+2) |
| T2 | descriptors vs PubChem, 10 drugs | pass — MW 10/10, HBD 10/10; other differences documented |
| T3 | strict is a design filter | pass — 0/3 strict, 3/3 relaxed |
| T4 | published thresholds never rewritten | pass — Lipinski detail identical across tiers |
| T5 | advisory does not reject | pass — 5/5 real actives pass, all 5 flagged |
| T6 | PAINS / Brenk / NIH | pass — three catalogues flag rhodanine independently |
| T7 | invalid SMILES isolated | pass — 2 errored, 2 screened, rc 0 |
| T8 | ChEMBL target search | pass — both organisms, with accessions |
| T9 | activity aggregation | pass — 34 records to 32 molecules, median and spread |
| T10 | ZINC encoding vs live server | pass — 15 tranches, URL returns 200 and real ids |
| T11 | Vina log parser | pass — 9 modes, mode 1's bare `0` parsed |
| T12 | missing Vina | pass — rc 4, platform-correct hint |
| T13 | output size and stdout bound | pass — 13× saving, 15 lines, identical verdicts |

**Environment:** Windows 11 Pro, Python 3.13.13, RDKit 2026.03.6,
Meeko 0.8.0, SciPy 1.18.1. AutoDock Vina **not installed**, so T11's
docking run was not executed; its parser case was. Live services queried:
ChEMBL, files.docking.org.

## Automated equivalent

```bash
pytest tests/test_chem_descriptors.py tests/test_chem_rules.py \
       tests/test_chem_alerts.py tests/test_chem_admet.py -v
pytest tests/test_docking.py -v                       # parser, no Vina needed
pytest tests/test_libraries_and_fetch.py -m network    # T8, T9, T10
pytest tests/test_docking.py -m requires_vina          # T11, if Vina is present
```

The manual protocol exists because it is readable, needs no test
environment, and shows what a correct result looks like — including the
places where the skill *should* disagree with an outside reference.
