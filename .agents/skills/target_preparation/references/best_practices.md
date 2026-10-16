# Structural biology practice for structure-based drug design

Read this before trusting a receptor you prepared. The skill automates the
mechanical steps; the judgements below stay with you. Nothing here is
organism-specific.

Numbered references are listed at the end and in [citation.bib](citation.bib).

--------------------------------------------------------------------------------

## 1. Judging an experimental structure

`assess-structure` deliberately returns `NOT_APPLICABLE` for experimental
structures, because the quantities that matter are not in the B-factor
column. Check these instead, from the entry page or the header:

| Quantity | What to want | Why |
|---|---|---|
| **Resolution** | ≤ 2.5 Å for SBDD; ≤ 2.0 Å to trust side-chain detail | Below ~2.8 Å, side-chain rotamers and waters are substantially modelled rather than observed [1] |
| **R-free** | Within ~0.05 of R-work, and roughly ≤ 0.28 | A large gap indicates over-fitting; R-free is the cross-validated statistic [2] |
| **Occupancy at the site** | 1.0 for pocket residues and the ligand | Partial occupancy means the atoms are only sometimes there |
| **Real-space fit of the ligand** | Good density, sensible geometry | Ligands are among the worst-modelled parts of deposited structures [3, 4] |
| **Alternate conformations** | Know which one you kept | `altloc` A is kept by convention, not by merit |
| **Missing residues** | None lining the pocket | A gap in a binding loop invalidates affinities for that site |

The wwPDB validation report for an entry gives percentile ranks for most of
this [3]. For local geometry, MolProbity's clashscore and rotamer outliers
are the standard check [4].

**A high-resolution structure has low B-factors.** That is the point of the
warning in [questions.md](questions.md) Q20: a mean B-factor near 25 Å² is a
sign of quality, and reading it as a confidence score inverts the meaning.

## 2. Which structure to choose

When several structures exist for one target:

1.  **Holo over apo.** A ligand-bound structure supplies the site and a
    pocket conformation that a real ligand selected.
2.  **Ligand chemistry close to yours.** Pockets are plastic. A site shaped
    around a fragment often cannot accommodate a larger series, and
    cross-docking between non-cognate receptor–ligand pairs loses accuracy
    measurably [5, 6].
3.  **Resolution and R-free**, per section 1.
4.  **Complete pocket.** Prefer the entry with no gaps and full occupancy
    where you intend to dock, even at slightly worse resolution.
5.  **The right oligomeric state.** See section 4.

A single rigid receptor encodes one conformational hypothesis. Where several
reasonable structures exist, docking against an ensemble is the honest
approach and usually the more predictive one [6, 7].

## 3. What to remove, and what not to

| Component | Default | Reasoning |
|---|---|---|
| Bulk water | Removed | Solvation is treated implicitly by docking scoring functions |
| Bridging / structural water | **Your call** | A water mediating a conserved ligand contact behaves as part of the receptor; retaining it can improve pose prediction markedly [8, 9] |
| Catalytic metals | Kept | Zn in a metalloprotease, Mg in a kinase, Fe in a P450 heme are part of the site |
| Buffer ions, cryoprotectants | Removed | Crystallisation artefacts: sulfate, glycerol, PEG, imidazole, acetate |
| Cofactors | Removed unless asked | Keep with `--keep-ligand` when the cofactor forms a wall of your pocket — NAD(P)H, FAD, heme and SAM commonly do |
| The co-crystallised ligand | Removed | It would occupy the site you are docking into; its coordinates are still used to place the box |

Deciding which waters are "relevant" has no reliable automatic answer;
published classifiers exist but should be applied knowingly [9]. The
defensible positions are to keep none or to keep a specific, justified few.

## 4. Oligomeric state and interface sites

Determine the **biological assembly**, not just the asymmetric unit. The
crystallographic asymmetric unit may contain more or fewer copies than the
functional molecule, and the PDB entry states the assembly separately.

A site at a subunit interface cannot be modelled from one subunit:

- HIV-1 protease is an obligate homodimer; the catalytic Asp25–Asp25' pair
  and the flaps that close over the ligand come from both chains.
- Many kinase and transporter sites, and most allosteric sites at dimer
  interfaces, behave the same way.

The skill warns when the site ligand's chain was dropped by `--chain`
([questions.md](questions.md) Q11), but it cannot know your target's
biology. Check the assembly.

## 5. Protonation, tautomers and flips

These change hydrogen-bonding and therefore every pose and score. X-ray
density at usual resolutions does not see hydrogens, so all of it is
assigned computationally and much of it is assigned wrongly by default:

- **Side-chain pKa shifts.** A buried aspartate or an active-site histidine
  routinely has a pKa far from its solution value. Use PROPKA or an
  equivalent rather than assuming standard states [10, 11].
- **Histidine tautomers.** Nδ1-H, Nε2-H and the charged form are three
  different hydrogen-bond donors/acceptors.
- **Asn/Gln amide flips and His ring flips.** The density cannot distinguish
  N from O (or C from N), so these are frequently modelled 180° out.
  Optimising the hydrogen-bond network corrects them [12, 13].
- **Ligand tautomer and charge state.** Decide it deliberately; a carboxylic
  acid docked neutral behaves nothing like its anion.

Receptor preparation choices of this kind measurably change virtual
screening enrichment [14], so they are not cosmetic.

## 6. Defining the search box

- Centre on a co-crystallised ligand where one exists — the most defensible
  choice available.
- **Keep the box small.** Enlarging it degrades results: the search space
  grows with volume, and spurious surface minima become reachable. Roughly
  2.9× the ligand's radius of gyration is a reasonable published guide [15].
- Blind docking over a whole protein is a last resort, not a default.
- Verify the box encloses the intended pocket, and that approximated side
  chains are not inside it ([questions.md](questions.md) Q14).

## 7. What a docking score is, and is not

- It is **not** a binding free energy, and correlation with measured
  affinity across unrelated chemotypes is weak [5, 16].
- Differences under about 1 kcal/mol do not rank compounds.
- Scoring functions are better at **pose prediction** than at affinity
  ranking, and better at enrichment than at either [5].
- Validate prospectively where you can: re-dock the cognate ligand and
  check the RMSD to the crystal pose; screen a set with known actives and
  decoys and measure enrichment [17, 18].
- Beware decoy bias when you do: DUD-E-style decoy sets can be gamed by
  property differences, and LIT-PCBA was built to address that [17, 18].

## 8. Predicted models

- Backbone accuracy for well-covered families is often near-experimental;
  **pocket side chains are the weak point** [19, 20].
- Check pLDDT locally in the pocket, not globally. Also read PAE/PTM for
  domain-domain placement if you have it.
- Predicted models carry no ligands, metals, waters or alternate
  conformations, and sit in a single, often apo-like, conformation.
- Virtual screening enrichment against predicted structures is
  systematically lower than against experimental ones [20, 21].
- Treat a prediction as a valuable hypothesis that accelerates but does not
  replace experiment [22].

## 9. Record what you did

Every choice above changes the result, so the cleaned file and the JSON
report both record provenance: the source database and version, the chains
kept, each heterogen removed with its reason, and the detected B-factor
scale. The cleaned PDB retains `HEADER`, `TITLE`, `EXPDTA`, `CRYST1` and
`REMARK 2 RESOLUTION` from the original so the file remains
self-describing — and so downstream tools, including this skill's own scale
detection and GROMACS' handling of the unit cell, can read it.

Cite the structure and the method. [citation.bib](citation.bib) has the
entries; the PDB asks that entries be cited by ID and primary reference [23].

--------------------------------------------------------------------------------

## References

1.  Wlodawer A, Minor W, Dauter Z, Jaskolski M. Protein crystallography for
    non-crystallographers. *FEBS J* 2008;275:1–21.
2.  Brünger AT. Free R value. *Nature* 1992;355:472–475.
3.  Read RJ et al. A new generation of crystallographic validation tools for
    the PDB. *Structure* 2011;19:1395–1412.
4.  Davis IW et al. MolProbity. *Nucleic Acids Res* 2007;35:W375–W383.
5.  Warren GL et al. A critical assessment of docking programs and scoring
    functions. *J Med Chem* 2006;49:5912–5931.
6.  Scior T et al. Recognizing pitfalls in virtual screening.
    *J Chem Inf Model* 2012;52:867–881.
7.  Bender BJ et al. A practical guide to large-scale docking.
    *Nat Protoc* 2021;16:4799–4832.
8.  Roberts BC, Mancera RL. Ligand–protein docking with water molecules.
    *J Chem Inf Model* 2008;48:397–408.
9.  Amadasi A et al. Robust classification of relevant water molecules in
    putative protein binding sites. *J Med Chem* 2008;51:1063–1067.
10. Olsson MHM et al. PROPKA3. *J Chem Theory Comput* 2011;7:525–537.
11. Søndergaard CR et al. Improved treatment of ligands and coupling
    effects. *J Chem Theory Comput* 2011;7:2284–2295.
12. Word JM et al. Asparagine and glutamine: using hydrogen atom contacts in
    the choice of side-chain amide orientation. *J Mol Biol*
    1999;285:1735–1747.
13. Hooft RWW, Sander C, Vriend G. Positioning hydrogen atoms by optimizing
    hydrogen-bond networks. *Proteins* 1996;26:363–376.
14. Madhavi Sastry G et al. Protein and ligand preparation: parameters,
    protocols, and influence on virtual screening enrichments.
    *J Comput Aided Mol Des* 2013;27:221–234.
15. Feinstein WP, Brylinski M. Calculating an optimal box size for ligand
    docking and virtual screening. *J Cheminform* 2015;7:18.
16. Chen YC. Beware of docking! *Trends Pharmacol Sci* 2015;36:78–95.
17. Mysinger MM et al. Directory of Useful Decoys, Enhanced (DUD-E).
    *J Med Chem* 2012;55:6582–6594.
18. Tran-Nguyen VK, Jacquemard C, Rognan D. LIT-PCBA.
    *J Chem Inf Model* 2020;60:4263–4273.
19. Jumper J et al. Highly accurate protein structure prediction with
    AlphaFold. *Nature* 2021;596:583–589.
20. Karelina M, Noh JJ, Dror RO. How accurately can one predict drug binding
    modes using AlphaFold models? *eLife* 2023;12:RP89386.
21. Scardino V, Di Filippo JI, Cavasotto CN. How good are AlphaFold models
    for docking-based virtual screening? *iScience* 2023;26:105920.
22. Terwilliger TC et al. AlphaFold predictions are valuable hypotheses and
    accelerate but do not replace experimental structure determination.
    *Nat Methods* 2024;21:110–116.
23. Berman HM et al. The Protein Data Bank. *Nucleic Acids Res*
    2000;28:235–242.
