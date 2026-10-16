# TEST.md — target-preparation

A verification protocol for the `target-preparation` skill. Work through it
to satisfy yourself the skill is accurate and functional on your machine.

**Not shipped with the skill.** This lives under `verification/` so the
bundle under `.agents/skills/` stays lean.

Every expected value below was **measured** by running the command against
the real data source, not copied from the documentation.

## Setup

```bash
cd <workshop root>
export S=.agents/skills/target_preparation/scripts/target_preparation_api.py
export W=/tmp/tp_test && mkdir -p $W
```

Run each case with `uv run $S ...`, which is what the SKILL.md documents:
`uv` reads the PEP 723 header at the top of the launcher and builds the
environment itself, so nothing needs installing first.

`python $S ...` works identically **only** in an interpreter that already has
the dependencies. Note that this is a real distinction, not a formality — on
the machine these results came from, the `python` first on `PATH` is a conda
install without `polite_http`, and the launcher fails on import there. If you
want the plain-`python` form, use the interpreter you installed the package
into (`.venv/Scripts/python.exe` on Windows, `.venv/bin/python` elsewhere).

| Need | For |
|---|---|
| `uv` | all cases |
| Internet | T1, T10 |
| Meeko or Open Babel | the PDBQT assertion in T1 |

--------------------------------------------------------------------------------

## T1 — prepare-receptor from the PDB

The core case: fetch, clean, convert, assess and derive the docking box in
one step.

```bash
uv run $S prepare-receptor --pdb-id 6HEZ --chain A --output $W/receptor.json
```

PDB **6HEZ** is *M. tuberculosis* DprE1, chosen because it exercises three
things a simpler structure would not: two chains, an FAD cofactor *and* a
drug-like inhibitor, and 16 residues with incomplete side chains.

**Expected stdout:**

```
Receptor: 6HEZ
  cleaned PDB: .../6HEZ_clean.pdb
  atoms kept: 3231, waters removed: 41
  heterogens removed: 4
  PDBQT: .../6HEZ_clean.pdbqt
  pdbqt warning: 16 residue(s) did not match a chemical template, most
    often because side-chain atoms are missing from the crystal structure.
    Meeko padded them from its templates: A:7, A:18, A:37, A:41, A:254,
    A:259, A:263, A:266, A:267, A:294, A:299, A:317 (+4 more)
  pdbqt warning: IMPORTANT: 1 of those residue(s) lie inside the docking
    box (A:317). Their geometry is approximated, so affinities for this
    site are less reliable.
  confidence verdict: NOT_APPLICABLE
  bound ligands: 4; site from 0SK:A:501
```

**Pass criteria**

| # | Check | Why it matters |
|---|---|---|
| 1 | `atoms kept: 3231` | Chain A only; chain B excluded |
| 2 | `waters removed: 41` | Solvent is actually removed, not merely documented as removed |
| 3 | `heterogens removed: 4` | Four species gone: the inhibitor `0SK`, the cofactor `FAD`, 41 `HOH`, and `IMD` - imidazole, classified as a buffer component rather than a ligand |
| 4 | `site from 0SK:A:501` | **The box is on the inhibitor, not the cofactor** |
| 5 | `A:317` named as inside the box | Approximated geometry *in the pocket* makes affinities unreliable; elsewhere it is harmless |
| 6 | `NOT_APPLICABLE` | A crystal structure must not receive a pLDDT verdict |
| 7 | Four files written | `6HEZ.pdb`, `6HEZ_clean.pdb`, `6HEZ_clean.pdbqt`, `binding_site.json` |

**Verify the receptor is a genuine PDBQT** — AutoDock atom types in the
last column, partial charges before them:

```bash
head -2 $W/6HEZ_clean.pdbqt
# ATOM      1  C   THR A   8     -12.379  -3.407  49.404  1.00  0.00    +0.245 C
# ATOM      2  O   THR A   8     -13.195  -3.213  50.306  1.00  0.00    -0.273 OA
```

Measured: 309,339 bytes, 3,819 `ATOM` records. That exceeds the 3,231 atoms
kept, because Meeko adds polar hydrogens - a receptor PDBQT with exactly as
many atoms as the cleaned PDB has not been protonated.

> **If the warnings say Open Babel was used instead of Meeko**, Meeko is
> not importable in the interpreter running the script. That is a valid
> fallback and the receptor is still usable, but Meeko assigns better
> partial charges. Check with `uv run python -c "import meeko"`. Note that
> a `mk_prepare_receptor` on `PATH` belonging to a *different* environment
> cannot be used — the skill now prefers its own interpreter's copy for
> exactly this reason.

--------------------------------------------------------------------------------

## T2 — list-ligands distinguishes cofactors from drug candidates

```bash
uv run $S list-ligands --pdb $W/6HEZ.pdb --output $W/ligands.json
```

**Expected stdout:**

```
Candidate ligands in 6HEZ.pdb: 4
  0SK:A:501  (28 atoms)
  0SK:B:501  (28 atoms)
  FAD:A:502  (53 atoms)
  FAD:B:502  (53 atoms)
```

**Pass criteria**

1.  Exactly 4 candidates — waters and buffer components excluded.
2.  **0SK is listed before FAD despite being smaller.** Ordering is by
    druggability, not size.
3.  In the JSON, `FAD` carries `"is_cofactor": true`; `0SK` carries
    `"is_cofactor": false` and `"role": "drug-like ligand"`.

```bash
python -c "
import json; d=json.load(open('$W/ligands.json'))
print([(l['resname'], l['is_cofactor']) for l in d['ligands']])"
# [('0SK', False), ('0SK', False), ('FAD', True), ('FAD', True)]
```

--------------------------------------------------------------------------------

## T3 — define-site from a bound ligand

```bash
uv run $S define-site --pdb $W/6HEZ.pdb --from-ligand --output $W/site.json
```

**Expected stdout:**

```
Binding site from co-crystallised ligand centroid
  reference: 0SK:A:501
  centre: (13.87, -20.62, 37.13)
  box: 18.1 x 18.0 x 14.9 A = 4870 A^3
  confidence: high - the box is centred on an experimentally observed
              bound ligand
  warning: 3 other heterogen(s) are present: 0SK:B:501 (28 atoms,
           drug-like ligand), FAD:A:502 (53 atoms, cofactor or prosthetic
           group), FAD:B:502 (53 atoms, cofactor or prosthetic group).
           If one of those is your target ligand, re-run with
           --ligand-resname.
```

**Pass criteria**

1.  Centre is `(13.87, -20.62, 37.13)` to two decimals — the inhibitor
    site.
2.  `n_reference_atoms` is 28: one 0SK copy, not both averaged together.
    Averaging across chains would centre the box in the space between them.
3.  The warning names the alternatives, so a user who wants the FAD site
    knows how to ask.

**Independent check that the box encloses the ligand** — arithmetic from
the coordinates, not the skill's own claim:

```bash
python -c "
import json,sys; sys.path.insert(0,'src')
from agskills.targets.site import parse_atoms
s=json.load(open('$W/site.json'))
a=[x for x in parse_atoms(open('$W/6HEZ.pdb').read())
   if x.is_hetatm and x.resname=='0SK' and x.chain=='A' and x.element!='H']
print('atoms:', len(a), 'all inside box:',
  all(abs(x.x-s['center_x'])<=s['size_x']/2 and
      abs(x.y-s['center_y'])<=s['size_y']/2 and
      abs(x.z-s['center_z'])<=s['size_z']/2 for x in a))"
# atoms: 28 all inside box: True
```

### T3b — the cofactor site is a genuinely different pocket

```bash
uv run $S define-site --pdb $W/6HEZ.pdb --from-ligand \
    --ligand-resname FAD --ligand-chain A --output $W/site_fad.json
```

```bash
python -c "
import json
a=json.load(open('$W/site.json')); b=json.load(open('$W/site_fad.json'))
d=sum((a[f'center_{k}']-b[f'center_{k}'])**2 for k in 'xyz')**.5
print(f'0SK box vs FAD box: {d:.1f} A  -- pass if > 5')"
# 0SK box vs FAD box: 13.0 A  -- pass if > 5
```

**Measured: 13.0 Å** — the FAD box centres at `(6.71, -9.95, 39.04)` and is
8835 Å³ against the inhibitor box's 4870 Å³, a cofactor-shaped cavity rather
than a drug-sized one. Confirms these are two distinct pockets, so getting
the choice wrong is not a cosmetic difference — it docks into the wrong
pocket.

--------------------------------------------------------------------------------

## T4 — define-site from pocket residues

Uses PDB **1M17**, EGFR kinase with erlotinib bound. Thr766, Met769 and
Leu820 are the gatekeeper and hinge residues of the ATP pocket.

```bash
uv run $S define-site --pdb tests/data/1M17.pdb \
    --residues A:766 A:769 A:820 --output $W/site_res.json
```

**Expected stdout:**

```
Binding site from pocket residue centroid
  reference: A:766, A:769, A:820
  centre: (23.26, -0.57, 57.99)
  box: 19.0 x 19.4 x 13.8 A = 5094 A^3
  confidence: medium - depends on the residue selection being correct
```

**Pass criterion — the structural-biology check.** A box built from those
pocket residues must land near the box built from erlotinib itself:

```bash
uv run $S define-site --pdb tests/data/1M17.pdb --from-ligand \
    --output $W/site_aq4.json
python -c "
import json
a=json.load(open('$W/site_res.json')); b=json.load(open('$W/site_aq4.json'))
d=sum((a[f'center_{k}']-b[f'center_{k}'])**2 for k in 'xyz')**.5
print(f'separation: {d:.1f} A  -- pass if < 10')"
# separation: 5.4 A  -- pass if < 10
```

**Measured: 5.4 Å.** The erlotinib box itself is centred at
`(22.01, 0.25, 52.79)` from `AQ4:A:999`, 29 atoms. A separation above 10 Å
would mean the residue numbering or the centroid arithmetic was wrong.

--------------------------------------------------------------------------------

## T5 — assess-structure on a predicted model

```bash
uv run $S assess-structure --pdb tests/data/AF-P9WJA5-F1.pdb \
    --output $W/af_q.json
```

**Expected stdout:**

```
Structure: AF-P9WJA5-F1.pdb
  residues: 125 in chain(s) A
  mean B-factor/pLDDT: 78.7 (median 82.9)
  bands: 3 very high, 100 confident, 18 low, 4 very low
  VERDICT: PROCEED - Mean pLDDT 78.7 with 82% of residues at 70 or above.
  action: 4 residue(s) are below 50; keep the docking box away from them.
  note: B-factors are on a 0-100 scale with a high mean, consistent with
        pLDDT from a predicted model.
```

**Pass criteria**

1.  125 residues — the real length of this AlphaFold entry.
2.  The bands sum to 125 (3 + 100 + 18 + 4): they partition the residues,
    so nothing is double-counted or lost.
3.  `scale` is `plddt_0_100` in the JSON.
4.  Verdict `PROCEED`, with the four low-confidence residues called out.

--------------------------------------------------------------------------------

## T6 — assess-structure refuses to misjudge a crystal structure

The single most important accuracy case.

```bash
uv run $S assess-structure --pdb tests/data/1M17.pdb --output $W/x_q.json
```

**Expected stdout:**

```
  residues: 312 in chain(s) A
  mean B-factor/pLDDT: 60.6 (median 50.6)
  bands: 56 very high, 30 confident, 73 low, 153 very low
  VERDICT: NOT_APPLICABLE - The B-factor column does not contain pLDDT,
           so a confidence verdict would be meaningless.
  action: Treat this as an experimental structure: check resolution,
          R-free and the occupancy of the binding site residues.
  action: Remove waters and crystallisation additives with 'prepare-receptor'
          before docking.
  note: ... LOWER is better. The pLDDT verdict below does not apply;
        judge the structure by its resolution and R-free instead.
```

**Pass criteria**

1.  Verdict is `NOT_APPLICABLE` and `scale` is `bfactor`.
2.  The note states **LOWER is better**.

A naive implementation reads mean 60.6 against a pLDDT scale and reports
"low confidence, use with caution" for a perfectly serviceable 2.6 Å
crystal structure. `NOT_APPLICABLE` here is the skill declining to answer a
question the data cannot support.

### T6b — a well-ordered structure must not be called unreliable

The most damaging bug found in this skill, and the reason T6 cannot be
judged on 1M17 alone.

```bash
for id in 1HSG 6LU7 2HYY 4EY7 5TZ1 1LDG 3QGT; do
  uv run $S assess-structure --pdb $W/$id.pdb --output $W/q.json \
    | grep -E "VERDICT"
done
```

**Measured — all seven `NOT_APPLICABLE`.** Previously four of them reported

```
VERDICT: STOP - Mean pLDDT 29.3 with 96% of residues below 50.
                Not reliable enough for docking.
```

for 1HSG, a 2.0 Å structure of HIV-1 protease.

**Pass criteria**

1.  Every experimental structure returns `NOT_APPLICABLE`, whatever its
    B-factor distribution.
2.  The predicted model still returns a pLDDT verdict —
    `tests/data/AF-P9WJA5-F1.pdb` gives `PROCEED`.

**Why magnitude alone cannot decide this.** Scale detection used to classify
any file whose B-factors stayed under 100 with a mean above 20 as a
predicted model. Measured across nine real entries:

| Entry | mean | max | Under the old rule |
|---|---:|---:|---|
| 1HSG | 29.3 | 59.7 | **misread as pLDDT → STOP** |
| 4EY7 | 31.2 | 93.6 | **misread as pLDDT → STOP** |
| 5TZ1 | 39.4 | 256.5 | correct (max > 100) |
| 6LU7 | 40.9 | 83.0 | **misread as pLDDT → STOP** |
| 2HYY | 42.2 | 98.3 | **misread as pLDDT → STOP** |
| 3QGT | 42.3 | 119.5 | correct (max > 100) |
| 1M17 | 60.6 | 120.0 | correct (max > 100) |
| 6HEZ | 62.4 | 155.8 | correct (max > 100) |
| 1LDG | 9.6 | 53.3 | correct (mean < 20) |

A well-refined structure has low B-factors, so **the better the
crystallography, the more confidently it was rejected.** The two fixtures
this skill shipped with, 1M17 and 6HEZ, both happen to have a value above
100 — which is why the original T6 passed and the bug survived.

The fix reads provenance instead. A wwPDB entry carries `EXPDTA` and
`REMARK 2 RESOLUTION`; the AlphaFold DB and ESMFold announce themselves in
`TITLE`:

```bash
python -c "
import sys; sys.path.insert(0,'src')
from agskills.targets.quality import provenance
print(provenance(open('$W/1HSG.pdb').read()))
print(provenance(open('tests/data/AF-P9WJA5-F1.pdb').read()))"
# experimental
# predicted
```

3.  **Cleaning preserves provenance**, or the fix would not survive the
    pipeline's own output. `PDBIO` writes coordinates only, so `HEADER`,
    `TITLE`, `EXPDTA`, `CRYST1` and `REMARK 2 RESOLUTION` are carried
    across:

```bash
head -5 $W/1HSG_clean.pdb
# HEADER    HYDROLASE (ACID PROTEINASE)             31-MAR-95   1HSG
# TITLE     CRYSTAL STRUCTURE AT 1.9 ANGSTROMS RESOLUTION OF HUMAN
# ...
# EXPDTA    X-RAY DIFFRACTION
# REMARK   2 RESOLUTION.    2.00 ANGSTROMS.
```

    `CRYST1` is kept for a second reason: GROMACS `pdb2gmx` and `editconf`
    read the unit cell when one is present.

4.  **A headerless file declines rather than rejecting.** Stripped to bare
    coordinates, a mean of 20–50 is genuinely ambiguous, so the fallback
    reads it as B-factors and reports `NOT_APPLICABLE`. Reporting `STOP`
    would discard good structures, which is the worse of the two errors. A
    headerless file is read as pLDDT only above a mean of 70, where it
    would pass anyway.

### T6c — the third scale: ESMFold's fractional pLDDT

The ESMFold API writes pLDDT between 0 and 1. Read as a 0–100 score, an 84
pLDDT model becomes 0.84 and is graded unusable.

```bash
python -c "
import sys; sys.path.insert(0,'src')
from agskills.targets.quality import detect_scale, extract_plddt
pdb=''.join(
  f'ATOM  {i:5d}  CA  ALA A{i:4d}       0.000   0.000   0.000  1.00{v:6.2f}'
  '           C\n'
  for i,v in enumerate([.40,.55,.78,.84,.88,.91,.94,.86,.80,.75],1))
print(detect_scale(extract_plddt(pdb)))"
# plddt_0_1
```

Running a real sequence through T10 and checking `quality.scale` gives the
same answer against the live endpoint.

--------------------------------------------------------------------------------

## T6d — a site ligand in a removed chain is flagged

```bash
uv run $S prepare-receptor --pdb-id 1HSG --chain A --output $W/mono.json
```

**Measured:**

```
  atoms kept: 757, waters removed: 44
  bound ligands: 1; site from MK1:B:902
  site warning: The site was defined from MK1:B:902, but chain B is not in
    the cleaned receptor (kept: A). The box is in the right place, yet part
    of the pocket has been removed with that chain. If the site lies at a
    subunit interface - as it does in an obligate dimer such as HIV-1
    protease - re-run without --chain, or keep every chain that
    contributes to the pocket.
```

**Pass criteria**

1.  The warning fires. HIV-1 protease is an obligate homodimer: the
    catalytic Asp25 pair and the flaps that close over the ligand come from
    both monomers, and the inhibitor is recorded as chain B. `--chain A`
    leaves the box correctly placed over half a site.
2.  Measured: 115 of the 757 remaining receptor atoms lie inside that box,
    so it is not empty — just incomplete, which is harder to notice.
3.  **It does not fire when it should not.** `1HSG` without `--chain`, and
    `1M17 --chain A` (whose ligand `AQ4:A:999` is in the chain kept), both
    print no warning.

--------------------------------------------------------------------------------

## T7 — cleaning options are honoured

```bash
# Keep the cofactor deliberately
uv run $S prepare-receptor --input-file $W/6HEZ.pdb --chain A \
    --keep-ligand FAD --output $W/with_fad.json
```

**Pass criteria**

1.  `cleaning.heterogens_kept` contains `FAD`, with a reason naming
    `--keep-ligand`.
2.  `FAD` no longer appears in `heterogens_removed`.
3.  Without the flag (T1), `FAD` appears in `heterogens_removed`. Retention
    is a choice that is recorded either way.

```bash
python -c "
import json; d=json.load(open('$W/with_fad.json'))['cleaning']
print('kept   :', [(h['resname'],h['reason'][:28]) for h in d['heterogens_kept']])
print('removed:', [h['resname'] for h in d['heterogens_removed']])"
```

--------------------------------------------------------------------------------

## T8 — validate-sequence accepts protein and refuses the rest

```bash
uv run $S validate-sequence --sequence \
"MKIEEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDIIFWAHDRFGGYAQSGLLAEITPDKAFQDKLYPFTWDAVRYNGKLIAYPIAVEALSLIYNKDL" \
    --output $W/seq.json
# Sequence is valid: 122 residues
```

Then each of these must be **refused**, with a message that says what to do:

| Input | Expected message | Exit |
|---|---|---:|
| `ATGGCGCATTGCAGCTAGCTAGCTAGCATCG` | `This looks like a DNA sequence, not a protein: it uses only the letters ACGT.` | 3 |
| 500 × `A` | `has 500 residues but the public ESMFold endpoint accepts at most 400` + points at `generate-af2-script` | 3 |
| `MKIEEGKLVIWINGXDKGYNG` | names `X (any residue)` as an ambiguity code | 3 |
| `MKT` | `at least 10` | 3 |

The DNA case is worth dwelling on: A, C, G and T are alanine, cysteine,
glycine and threonine, and N is asparagine, so a nucleotide sequence is
*valid protein* by character. It is caught on composition instead — and a
real protein must not trip that check:

```bash
uv run $S validate-sequence --sequence "ACGTNACGTN" --output $W/s.json
# Sequence is valid: 10 residues   <- short, so no composition verdict
```

--------------------------------------------------------------------------------

## T9 — failure paths return documented exit codes

Exit codes are part of the interface. Run each **unpiped**, because a pipe
reports the pipe's status rather than the command's.

```bash
uv run $S validate-sequence --sequence "ATGGCGCATTGCAGCTAGCTAGCTAGCATCG" \
    --output $W/x.json >/dev/null 2>&1; echo $?     # 3
uv run $S prepare-receptor --pdb-id NOTANID --output $W/x.json \
    >/dev/null 2>&1; echo $?                        # 3
uv run $S prepare-receptor --uniprot-id EGFR --output $W/x.json \
    >/dev/null 2>&1; echo $?                        # 3
uv run $S define-site --pdb tests/data/1M17.pdb --from-ligand \
    --whole-structure --output $W/x.json >/dev/null 2>&1; echo $?   # 2
uv run $S assess-structure --pdb nope.pdb --output $W/x.json \
    >/dev/null 2>&1; echo $?                        # 6
```

**Measured: 3, 3, 3, 2, 6.**

| Code | Meaning |
|---|---|
| 0 | Success |
| 2 | Bad usage — conflicting or missing flags |
| 3 | Invalid input data |
| 4 | A required dependency or binary is missing |
| 5 | A remote service failed |
| 6 | A required local file does not exist |

**Pass criteria**

1.  The codes match the table.
2.  No traceback is printed — the exception types are part of the interface,
    so a deliberate failure is one line plus a hint.
3.  No output file is created on failure.

`--uniprot-id EGFR` is worth its own case: a gene symbol is the commonest
user error, and the message must point at uniprot.org rather than simply
failing to fetch.

--------------------------------------------------------------------------------

## T10 — AlphaFold and ESMFold, live

```bash
uv run $S prepare-receptor --uniprot-id P9WJA5 --output $W/af/receptor.json
```

**Pass criteria**

1.  A file named `AF-P9WJA5-F1-v<N>.pdb` is written, where `<N>` is the
    version the API reported — **not** a hardcoded one. At the time of
    writing that is v6; v4 and v5 both return HTTP 404.
2.  `source.metadata.model_version` is present.
3.  `source.metadata.caveat` warns that a predicted model carries no
    ligands and therefore no crystallographic binding site.
4.  `quality.verdict` is a pLDDT verdict, since this *is* a predicted model.

**Measured:** `AF-P9WJA5-F1-v6.pdb`, `model_version: 6`, 1023 atoms kept, 0
waters, 0 heterogens, `verdict: PROCEED`, `scale: plddt_0_100`, mean 78.68.
The caveat reads *"A predicted model carries no ligands, cofactors, metals or
waters, and no crystallographic binding site."* and the next step offered is
`define-site --residues`, not `--from-ligand` — which is the right advice for
a structure that has no ligand to centre on.

```bash
uv run $S predict-structure --sequence \
"MKIEEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDIIFWAHDRFGGYAQSGLLAEITPDKAFQDKLYPFTWDAVRYNGKLIAYPIAVEALSLIYNKDL" \
    --output $W/predicted.pdb
```

5.  `predicted.pdb` is a real PDB (begins `HEADER` or has `ATOM` records),
    not an HTML error page served with HTTP 200.
6.  `predicted.json` beside it reports `quality.scale` as `plddt_0_1` —
    the fractional case from T6c, now against the live service.

**Measured:** 122 residues, `mean pLDDT: 84.16`, `PROCEED`. The raw file
confirms why the scale check matters:

```bash
grep '^ATOM' $W/predicted.pdb | head -2
# ATOM      1  N   MET A   1     -10.307  -5.714 -13.303  1.00  0.36           N
# ATOM      2  CA  MET A   1      -9.505  -4.532 -13.001  1.00  0.40           C
```

The endpoint serves **0.36 and 0.40**, not 36 and 40. Taken at face value on
a 0–100 scale this 84-pLDDT model reports a mean of 0.84 and is graded
unusable. `detect_scale` catches it and rescales, which is the difference
between using this model and discarding it.

--------------------------------------------------------------------------------

## T11 — generate-af2-script

```bash
uv run $S generate-af2-script --fasta target.fasta --method colabfold \
    --output $W/submit_af2.sh
```

**Pass criteria**

1.  `bash -n $W/submit_af2.sh` exits 0 — it is valid shell.
2.  It has **no CR bytes**. A CRLF script fails on a Linux cluster with
    `bash\r: bad interpreter` before it runs a single line. Check it at the
    byte level:

    ```bash
    python -c "print(open('$W/submit_af2.sh','rb').read().count(b'\r'))"
    # 0
    ```

    Do **not** use `grep -c $'\r'` for this. Under Git Bash on Windows the
    shell may leave grep an empty pattern that matches every line — it
    reported 49 of 49 lines on a file with zero CR bytes. Measured here:
    **0**. `write_text` normalises CRLF on every generated artefact, so this
    holds whatever platform produced the file.
3.  It checks its input exists and that `colabfold_batch` is on PATH
    *before* the long step.
4.  `conda activate` goes through `conda shell.bash hook`, which a
    non-interactive SLURM shell otherwise lacks.
5.  An invalid walltime is refused at generation time:
    `--time 24h` exits 2 with a message naming `HH:MM:SS`.

--------------------------------------------------------------------------------

## T12 — output size stays bounded

The JSON is read back by whatever called the skill, so its size is a cost.

```bash
uv run $S prepare-receptor --pdb-id 6HEZ --chain A \
    --output $W/compact.json
uv run $S prepare-receptor --pdb-id 6HEZ --chain A --detail full \
    --output $W/full.json
ls -l $W/compact.json $W/full.json
```

**Pass criteria**

1.  Measured: `compact.json` 5,713 bytes, `full.json` 11,481 bytes - a 2.0x
    difference, or roughly 1,400 tokens saved per call.
2.  `compact.json` carries `quality.verdict` and `quality.reason` but not
    the per-region breakdown — `assess-structure` provides that on demand.
3.  Both report the same verdict. Cheaper output must not mean a different
    answer.

--------------------------------------------------------------------------------

## Results

Recorded from a full pass. Reproduce by working through the cases above.

| Case | What it checks | Result |
|---|---|---|
| T1 | prepare-receptor on 6HEZ | pass — 3231 atoms, 41 waters, site on 0SK |
| T2 | cofactor vs drug-like ranking | pass — 0SK before FAD |
| T3 | site from ligand, box encloses it | pass — 28/28 atoms inside |
| T3b | cofactor site is distinct | pass — 13.0 Å apart |
| T4 | site from pocket residues | pass — 5.4 Å from the erlotinib box |
| T5 | predicted model assessed | pass — PROCEED, bands sum to 125 |
| T6 | crystal structure not misjudged | pass — NOT_APPLICABLE |
| T6b | well-ordered structures not rejected | pass — **after a fix**; 9/9 `NOT_APPLICABLE` |
| T6c | fractional pLDDT detected | pass — `plddt_0_1` |
| T6d | interface site flagged | pass — **after a fix**; fires on 1HSG `--chain A` only |
| T7 | cleaning options honoured | pass |
| T8 | sequence validation | pass — DNA, length and ambiguity all refused |
| T9 | exit codes | pass — 3, 3, 3, 2, 6 |
| T10 | AlphaFold version resolution, ESMFold | pass — v6 resolved, not hardcoded |
| T11 | generated script valid | pass — `bash -n` clean, no CRLF |
| T12 | output size | pass — ~6 KB compact |

**Environment:** Windows 11 Pro, Python 3.13.13, RDKit 2026.03.6,
Meeko 0.8.0, BioPython 1.88, SciPy 1.18.1, gemmi 0.7.5, Open Babel 3.2.1.
Data sources queried live: RCSB PDB, AlphaFold DB (v6), ESM Atlas.

## Automated equivalent

These cases are also asserted programmatically, so a regression is caught
without a manual pass:

```bash
pytest tests/test_targets_site.py tests/test_targets_quality.py \
       tests/test_targets_clean.py tests/test_targets_pdbqt.py -v
pytest tests/test_libraries_and_fetch.py -v          # identifier validation
pytest tests/test_libraries_and_fetch.py -m network  # T10, live
```

148 tests cover this skill. The manual protocol exists because it is
readable, needs no test environment, and shows *what good looks like*.
