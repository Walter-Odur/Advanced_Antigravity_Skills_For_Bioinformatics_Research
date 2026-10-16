# Questions this skill answers

Find your question, run the command. If your question is not here, the
nearest one usually is — the commands are general, not target-specific.

These skills are **organism-agnostic**. Every worked example below was run
against a real structure, and the set spans human, viral, bacterial, fungal
and protozoan targets deliberately: nothing in the code knows or cares which
organism a structure came from.

## Start here

| You have | Run |
|---|---|
| A PDB ID (e.g. `1M17`) | `prepare-receptor --pdb-id 1M17` |
| A UniProt accession (e.g. `P00533`) | `prepare-receptor --uniprot-id P00533` |
| A PDB/mmCIF file on disk | `prepare-receptor --input-file my.pdb` |
| An amino-acid sequence | `validate-sequence` then `predict-structure` |
| Only a gene or protein name | **Nothing yet** — see Q2 |
| A structure and no idea if it is usable | `assess-structure` |

--------------------------------------------------------------------------------

## Getting started

**Q1. I have never used this. What is the shortest path from a PDB ID to
something I can dock into?**

```bash
uv run scripts/target_preparation_api.py prepare-receptor \
    --pdb-id 1M17 --output work/receptor.json
```

One command fetches, cleans, converts to PDBQT, assesses quality and writes
`binding_site.json`. Read `work/receptor.json`; the console summary is a
convenience, not the result.

**Q2. I only know my target's name — "EGFR", "main protease", "DprE1".**

The skill will refuse a name, on purpose: gene symbols are ambiguous across
organisms and map to many structures of differing quality. Resolve it first:

- Search <https://www.rcsb.org> for an experimental structure. Prefer one
  with a bound ligand resembling your chemistry (Q10).
- Or search <https://www.uniprot.org> for the accession, then use
  `--uniprot-id` to pull the AlphaFold model.
- `ag-compound-screening search-chembl-target --query "<name>"` returns
  ChEMBL targets **with their organism and UniProt accession**, which is
  usually the fastest way to turn a name into an identifier.

**Q3. Which should I use, a crystal structure or an AlphaFold model?**

An experimental holo structure, whenever one exists. It gives you the
receptor *and* the binding site, in a conformation a real ligand selected.
A predicted model gives neither: no ligand, no waters, no metals, and
side-chain placement that was never validated against density [12, 15].

Use a predicted model when there is no experimental structure, and read Q18
before docking into one.

**Q4. How do I tell whether a structure is good enough?**

```bash
uv run scripts/target_preparation_api.py assess-structure \
    --pdb my.pdb --output quality.json
```

It reports which scale the B-factor column is on and gives a verdict. For a
**predicted model** it grades pLDDT (`PROCEED` / `REFINE` / `STOP`). For an
**experimental structure** it returns `NOT_APPLICABLE` and tells you to judge
by resolution and R-free instead — see [best_practices.md](best_practices.md),
which covers what numbers to look for.

**Q5. What do I do with the files it writes?**

| File | Goes to |
|---|---|
| `<ID>_clean.pdbqt` | `ag-compound-screening dock --receptor` |
| `binding_site.json` | `ag-compound-screening dock --site` |
| `<ID>_clean.pdb` | `ag-md-simulation gmxapi-setup --pdb` |
| `receptor.json` | Read it; it records every decision taken |

--------------------------------------------------------------------------------

## Choosing and checking the binding site

**Q6. Where will it dock? How do I know the box is in the right place?**

If the structure has a bound ligand, the box is centred on it and that is the
most reliable option available [19]. Check the reference it chose:

```bash
python -c "
import json; s=json.load(open('work/binding_site.json'))
print(s['reference'], s['confidence'])"
```

**Q7. My protein has several heterogens. Did it pick the right one?**

This is the commonest silent error in the whole pipeline. Cofactors are
large, drug-like and bind tightly, so a naive "biggest ligand" rule picks
the cofactor and docks into the wrong pocket.

```bash
uv run scripts/target_preparation_api.py list-ligands \
    --pdb my.pdb --output ligands.json
```

Candidates are ranked by **druggability, not size**, and cofactors are
flagged. Measured on real entries:

| Structure | Ranked first | Cofactors flagged | Note |
|---|---|---|---|
| 6HEZ | `0SK` (28 atoms) | `FAD` (53 atoms) | the larger heterogen is the cofactor |
| 5TZ1 | `VT1` (37) | `HEM` (43) | heme is not a docking site |
| 1LDG | `OXM` (6) | `NAI` (NADH, 44) | a 6-atom fragment outranks NADH |
| 3QGT | `CP6` (17) | `NDP`, `UMP` | two cofactors, one drug |

To override: `define-site --from-ligand --ligand-resname FAD --ligand-chain A`.

**Q8. There is no ligand at all. Now what?**

```bash
uv run scripts/target_preparation_api.py define-site \
    --pdb apo.pdb --residues A:766 A:769 A:820 --output site.json
```

Give residues you have grounds to believe line the pocket — catalytic
residues, a known hotspot, or residues equivalent to those contacting a
ligand in a homologue. Then sanity-check it: a box built from pocket
residues should land close to one built from a ligand in a related holo
structure. Measured on EGFR, the pocket-residue box sits **5.4 Å** from the
erlotinib box, which is within a box width.

Verified apo example: `2J3Z` (*Clostridium botulinum* C2 toxin) reports
`bound ligands: 0`, so `--residues` is the only route.

If you have no idea where the site is, this skill cannot help — run a pocket
detection tool (fpocket, SiteMap, PrankWeb) first, or superpose a holo
homologue.

**Q9. How big should the box be?**

The default is the ligand's extent plus padding, which is the right habit: a
box much larger than the pocket costs search accuracy, because the sampling
problem grows with volume [19]. As a rule of thumb keep it within roughly
2.9× the ligand's radius of gyration [19]. `--whole-structure` exists for
blind docking but expect worse results.

**Q10. The ligand in my chosen structure is nothing like my compounds.**

Then prefer a different structure. Binding sites are plastic, and a pocket
shaped around a small fragment may not accommodate your series at all; the
reverse holds too. Choose the holo structure whose ligand most resembles
your chemistry, and treat any single rigid receptor as one hypothesis among
several [20, 21].

--------------------------------------------------------------------------------

## Cleaning decisions

**Q11. Should I keep one chain or all of them?**

Keep every chain that contributes to the pocket. The trap is an obligate
oligomer whose site sits at a subunit interface: the ligand is recorded
under one chain, but the pocket is built from several.

HIV-1 protease is the textbook case. In `1HSG` the inhibitor is chain B, and
the catalytic Asp25 and the flaps come from both monomers:

```bash
# Wrong for this target
prepare-receptor --pdb-id 1HSG --chain A
#   atoms kept: 757
#   site warning: The site was defined from MK1:B:902, but chain B is not
#   in the cleaned receptor (kept: A). ... If the site lies at a subunit
#   interface - as it does in an obligate dimer such as HIV-1 protease -
#   re-run without --chain.
```

The skill warns whenever the site ligand's chain was removed. Heed it.
Use `--chain` to drop a genuinely separate copy of a monomeric protein, not
to reduce an oligomer.

**Q12. Why were my waters deleted, and should I have kept some?**

Waters are removed by default because most are bulk solvent and a docking
program models solvation implicitly. But a few are structural: bridging
waters that mediate ligand binding should be treated as part of the receptor
[16, 17]. There is no reliable automatic rule, so keeping them is a
deliberate act: `--keep-waters` retains all of them, after which you delete
the ones you do not want by hand.

**Q13. What happened to the metals?**

Kept by default, because they are frequently catalytic — the zinc of a
metalloprotease, the magnesium of a kinase, the iron of a P450's heme.
Buffer counter-ions and cryoprotectants are always removed. Use
`--remove-metals` only if you know the ion is adventitious.

**Q14. There are warnings about residues that "did not match a chemical
template". Does that matter?**

It matters only if one is in your pocket, and the report tells you which:

```
pdbqt warning: 16 residue(s) did not match a chemical template, most often
  because side-chain atoms are missing from the crystal structure.
pdbqt warning: IMPORTANT: 1 of those residue(s) lie inside the docking box
  (A:317). Their geometry is approximated, so affinities for this site are
  less reliable.
```

Incomplete side chains are routine in crystal structures — the density was
not there, so the atoms were not modelled [1]. Padded from a template they
are a guess. A guess at the edge of the protein is harmless; a guess lining
the pocket undermines every affinity you compute for it. Rebuild those side
chains, or pick a structure with better occupancy there.

**Q15. Should I add hydrogens, and at what pH?**

The PDBQT converter adds the polar hydrogens AutoDock needs. For anything
more careful — and especially for MD — assign protonation states explicitly
with PROPKA or an equivalent [8, 9], and check the His tautomers and the
Asn/Gln side-chain amide orientations, which are routinely modelled the
wrong way round because the density cannot distinguish N from O [10, 11].
A mis-protonated catalytic histidine will quietly ruin both docking and
simulation.

--------------------------------------------------------------------------------

## Predicted structures

**Q16. I have a sequence, not a structure.**

```bash
uv run scripts/target_preparation_api.py validate-sequence \
    --sequence MKIEEGKLVIW... --output seq.json
uv run scripts/target_preparation_api.py predict-structure \
    --sequence MKIEEGKLVIW... --output predicted.pdb
```

Validate first: it catches a DNA sequence pasted by mistake, ambiguity
codes, and anything over the 400-residue ESMFold limit, before you wait on
the service.

**Q17. My sequence is too long, or I want a better model.**

```bash
uv run scripts/target_preparation_api.py generate-af2-script \
    --fasta target.fasta --method colabfold --output submit_af2.sh
```

ESMFold is single-sequence and fast; AlphaFold2 with an MSA is more accurate
[12, 13]. The script is ready to `sbatch`.

**Q18. Can I dock into an AlphaFold model?**

Sometimes, with care, and the literature is explicit about the limits.
Predicted models give near-experimental backbones for well-covered families
but side-chain placement in the pocket is the weak point, and virtual
screening enrichment against predicted structures is systematically worse
than against experimental ones [14, 15]. Treat a predicted model as a
hypothesis [16].

Practical rules:

- Check pLDDT **in the pocket**, not the global mean. `assess-structure`
  lists low-confidence regions; a 90-pLDDT protein with a 40-pLDDT binding
  loop is not usable for docking there.
- There is no ligand, so you must define the site yourself (Q8).
- Expect to need ensemble or flexible-side-chain docking rather than a
  single rigid run.

**Q19. The verdict says `STOP`. Is the structure useless?**

For docking at that site, effectively yes — but check it is really a
prediction. `STOP` is a pLDDT verdict; it should never appear for a crystal
structure. If it does on an experimental file, the B-factor column has been
misread, and you should report it. The skill reads the file's `EXPDTA` and
`REMARK 2 RESOLUTION` records to decide, so this only arises for files
stripped of their headers.

**Q20. What is pLDDT and why does my crystal structure have "low" values?**

Because they are not pLDDT. A predicted model stores per-residue confidence
in the B-factor column, where **higher is better**, capped at 100. A crystal
structure stores atomic displacement parameters in the same column, where
**lower is better** and the scale is unbounded. A beautifully ordered 2.0 Å
structure has B-factors around 20–30, which read as pLDDT would mean
"unusable". `assess-structure` reports which scale it detected; see
[confidence.md](confidence.md).

--------------------------------------------------------------------------------

## Worked examples across organisms

Every row was run through `prepare-receptor`; the ligand named is the one the
skill selected for the box.

| PDB | Organism | Target | Site ligand | Why it is here |
|---|---|---|---|---|
| `1M17` | *Homo sapiens* | EGFR kinase | AQ4 (erlotinib) | clean single-chain baseline |
| `2HYY` | *Homo sapiens* | ABL1 kinase | STI (imatinib) | multi-chain, 4 heterogens |
| `4EY7` | *Homo sapiens* | acetylcholinesterase | E20 (donepezil) | 9 heterogens, 552 waters |
| `1HSG` | HIV-1 | protease | MK1 (indinavir) | obligate dimer; interface site (Q11) |
| `6LU7` | SARS-CoV-2 | main protease | PJE (N3) | peptidomimetic recorded as its own chain |
| `6HEZ` | *M. tuberculosis* | DprE1 | 0SK | FAD cofactor trap (Q7) |
| `5TZ1` | *C. albicans* | CYP51 | VT1 (oteseconazole) | heme cofactor |
| `1LDG` | *P. falciparum* | lactate dehydrogenase | OXM (oxamate) | NADH cofactor; fragment-sized site ligand |
| `3QGT` | *P. falciparum* | DHFR-TS | CP6 (pyrimethamine) | bifunctional, two cofactors |
| `2J3Z` | *C. botulinum* | C2 toxin I | none | apo: `--residues` required (Q8) |
| `P9WJA5` | *M. tuberculosis* | HigB1 | none | AlphaFold model, pLDDT path |

--------------------------------------------------------------------------------

## When something goes wrong

| Symptom | Cause | Fix |
|---|---|---|
| `'EGFR' is not a UniProt accession` | A gene name was passed | Q2 |
| `'NOTANID' is not a PDB entry id` | Four characters, first a digit | Check at rcsb.org |
| `This looks like a DNA sequence` | Nucleotides pasted as protein | Translate first |
| `has 500 residues but ... at most 400` | ESMFold limit | Q17 |
| `File not found` (exit 6) | Earlier step not run | Run it |
| `no CA atoms` | mmCIF, an HTML error page, or a ligand-only file | Supply a protein PDB |
| PDBQT not produced | No converter importable | `pip install meeko scipy gemmi` |
| Open Babel used instead of Meeko | Meeko not importable in *this* interpreter | Q21 |

**Q21. Why does it say Open Babel when I have Meeko installed?**

Because Meeko is not importable in the interpreter running the script. A
`mk_prepare_receptor` on `PATH` belonging to a *different* environment
cannot be used: it would import that environment's RDKit, not yours. Check
with `uv run python -c "import meeko"`. Open Babel is a valid fallback, but
Meeko assigns better partial charges.

**Q22. Exit codes?**

`0` success · `2` bad usage · `3` invalid input · `4` missing dependency ·
`5` remote service failed · `6` missing local file.

**Q23. The output is enormous and fills my context.**

`--detail compact` is already the default. Read fields out of the file
rather than printing it:

```bash
python -c "
import json; d=json.load(open('work/receptor.json'))
print(d['quality']['verdict'], d['binding_site']['reference'])"
```

`--detail full` adds every region and descriptor and costs several times as
much; ask for it only when you need it.

--------------------------------------------------------------------------------

## References

Numbers cited above. Full entries in [citation.bib](citation.bib).

1.  Wlodawer A, Minor W, Dauter Z, Jaskolski M. *FEBS J* 2008;275:1–21 —
    how to read a published structure without over-trusting it.
2.  Brünger AT. *Nature* 1992;355:472–475 — R-free.
3.  Read RJ et al. *Structure* 2011;19:1395–1412 — wwPDB validation.
4.  Davis IW et al. *Nucleic Acids Res* 2007;35:W375–W383 — MolProbity.
5.  Berman HM et al. *Nucleic Acids Res* 2000;28:235–242 — the PDB.
6.  Madhavi Sastry G et al. *J Comput Aided Mol Des* 2013;27:221–234 —
    receptor and ligand preparation, and its effect on screening.
7.  Forli S et al. *Nat Protoc* 2016;11:905–919 — the AutoDock protocol.
8.  Olsson MHM et al. *J Chem Theory Comput* 2011;7:525–537 — PROPKA3.
9.  Søndergaard CR et al. *J Chem Theory Comput* 2011;7:2284–2295 — pKa.
10. Word JM et al. *J Mol Biol* 1999;285:1735–1747 — Asn/Gln flips.
11. Hooft RWW, Sander C, Vriend G. *Proteins* 1996;26:363–376 — H-bond
    network optimisation.
12. Jumper J et al. *Nature* 2021;596:583–589 — AlphaFold.
13. Lin Z et al. *Science* 2023;379:1123–1130 — ESMFold.
14. Karelina M, Noh JJ, Dror RO. *eLife* 2023;12:RP89386 — how accurately
    AlphaFold models predict binding modes.
15. Scardino V, Di Filippo JI, Cavasotto CN. *iScience* 2023;26:105920 —
    AlphaFold models in docking-based virtual screening.
16. Terwilliger TC et al. *Nat Methods* 2024;21:110–116 — predictions are
    hypotheses.
17. Roberts BC, Mancera RL. *J Chem Inf Model* 2008;48:397–408 — docking
    with water molecules.
18. Amadasi A et al. *J Med Chem* 2008;51:1063–1067 — classifying relevant
    waters.
19. Feinstein WP, Brylinski M. *J Cheminform* 2015;7:18 — optimal docking
    box size.
20. Warren GL et al. *J Med Chem* 2006;49:5912–5931 — docking programs and
    scoring functions assessed.
21. Chen YC. *Trends Pharmacol Sci* 2015;36:78–95 — beware of docking.
22. Varadi M et al. *Nucleic Acids Res* 2022;50:D439–D444 — AlphaFold DB.
