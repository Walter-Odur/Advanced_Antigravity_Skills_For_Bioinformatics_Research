# Questions this skill answers

This skill runs the whole pipeline. If your question is about one step, the
sub-skill's own question index is more detailed:

- [target-preparation](../../target_preparation/references/questions.md) —
  structures, binding sites, model confidence
- [compound-screening](../../compound_screening/references/questions.md) —
  libraries, docking, ADMET
- [compound-synthesis](../../compound_synthesis/references/questions.md) —
  generative design
- [md-simulation](../../md_simulation/references/questions.md) — MD setup

## Start here

| You want | Run |
|---|---|
| To see what would happen, changing nothing | `plan` |
| To understand the stages | `describe` |
| To actually run it | `run --stages ...` |

--------------------------------------------------------------------------------

## Getting started

**Q1. I know nothing about any of this. Where do I begin?**

```bash
uv run scripts/drug_discovery_api.py plan \
    --pdb-id 1M17 --chembl-target CHEMBL203 --output plan.json
```

`plan` is **read-only**. It fetches nothing, installs nothing and takes
about two seconds. It tells you, stage by stage, what the pipeline would do
and whether this machine can do it:

```
  1. target      [local  ] [ready]     target-preparation
  2. sourcing    [local  ] [ready]     compound-screening
  3. triage      [local  ] [ready]     compound-screening
  4. docking     [cluster] [NOT READY] AutoDock Vina is not installed; this
                                       stage will produce a cluster script
  5. generative  [cluster] [NOT READY] 4 problem(s) ... Config generation
                                       and seed preparation work regardless
  6. selection   [local  ] [ready]     compound-screening
  7. dynamics    [cluster] [NOT READY] GROMACS is not installed; setup files
                                       and a cluster script are still produced
```

Read it, then run the stages you want.

**Q2. What are the seven stages for?**

| Stage | Does | Needs |
|---|---|---|
| `target` | fetch, clean, PDBQT, assess, define the site | a PDB ID, UniProt accession or file |
| `sourcing` | retrieve known actives and natural products | a ChEMBL target ID |
| `triage` | relaxed ADMET (14 checks), writes a seed set | compounds from `sourcing` |
| `docking` | rank by binding affinity | receptor, site, Vina |
| `generative` | REINVENT config and submission script | seeds, REINVENT, a GPU |
| `selection` | strict ADMET (24 checks) | compounds from anywhere upstream |
| `dynamics` | MDP files, workflow, submission script | a cleaned PDB |

**Q3. What can I run on an ordinary laptop?**

`target`, `sourcing`, `triage` and `selection` run locally. The three
cluster stages still produce their configuration and submission scripts, so
a full seven-stage pass completes in about 21 seconds with no Vina, no
REINVENT and no GROMACS installed. Nothing is faked: `docking` reports
`deferred to cluster` and a script path, never an invented affinity.

**Q4. Just run the useful part.**

```bash
uv run scripts/drug_discovery_api.py run \
    --pdb-id 1M17 --chain A --chembl-target CHEMBL203 \
    --stages target sourcing triage \
    --workdir pipeline --output pipeline.json
```

Measured: about 17 seconds for a fetch, clean, PDBQT conversion, ChEMBL
query and a 14-check screen.

**Q5. Where does the output go?**

Numbered stage directories, so the order is legible from the filesystem
months later:

```
01_target/      6HEZ.pdb, 6HEZ_clean.pdb, 6HEZ_clean.pdbqt, binding_site.json
02_sourcing/    chembl_actives.json
03_triage/      relaxed_admet.json, seeds.smi
04_docking/     submit_vina.sh
05_generative/  reinvent.toml, reinvent.meta.json, submit_reinvent.sh
06_selection/   strict_admet.json
07_dynamics/    em.mdp ions.mdp md.mdp npt.mdp nvt.mdp run_md.py setup_md.sh submit_md.sh
```

`pipeline.json` records verdicts, counts and paths — not the contents of
those files. A seven-stage report is only ~2.6 KB larger than a three-stage
one, because detail stays where it was written.

--------------------------------------------------------------------------------

## Using it on your own target

**Q6. My target is not a kinase, and not from the organism in your
examples.**

That makes no difference. Nothing in the pipeline knows what organism a
structure came from. The only inputs are an identifier and, optionally, a
ChEMBL target; both are organism-neutral.

Verified across human, viral, bacterial, fungal and protozoan targets — see
the table in
[target-preparation questions](../../target_preparation/references/questions.md#worked-examples-across-organisms).

What *does* change with your target:

- **The scoring profile** for the generative stage. The default suits oral
  small molecules; a Gram-negative bacterial target or a
  protein–protein-interaction site needs different windows. See
  [compound-synthesis Q7](../../compound_synthesis/references/questions.md).
- **The ADMET tier** you believe. The strict tier is a lead-like window;
  fragments, natural products and beyond-rule-of-5 chemistry should be
  judged at `relaxed`.
- **The force field** for dynamics, if your system has metals, cofactors or
  non-standard residues.

**Q7. I have a sequence, not a structure.**

Use the `target-preparation` skill directly to fold it first, then feed the
resulting PDB in with `--input-file`. Read
[target-preparation Q18](../../target_preparation/references/questions.md)
before docking into a predicted model — the limits are real and documented
[3, 4, 5].

**Q8. I have no target structure at all, only compounds.**

Then skip `target` and `docking`. `triage` and `selection` work on compounds
alone, and `generative` needs only seeds. Ligand-based work is a legitimate
route; this pipeline is simply structure-based where a structure exists.

**Q9. Can I stop partway and resume?**

Run the stages you want, then run the rest later against the same
`--workdir`. Each stage reads its predecessors' artefacts from disk. A stage
with nothing to work on is **skipped with its reason** rather than crashed:

```
  triage: skipped (no compounds)
```

--------------------------------------------------------------------------------

## Reading the results

**Q10. The pipeline says 25/25 passed triage but 5/25 passed selection. Did
something break?**

No — that is the funnel doing its job. `triage` is the relaxed tier (14
checks), `selection` the strict tier (24). A two-tier screen where both
tiers give the same count is not two tiers. See
[compound-screening Q8 and Q9](../../compound_screening/references/questions.md).

**Q11. How do I know the wizard is doing the same thing as the individual
skills?**

Because it delegates rather than reimplementing. Measured: the binding site
the wizard produces for 6HEZ is identical to the one
`target-preparation define-site` produces alone — same centre to 1e-12, same
box, same reference ligand. If they ever differ, the wizard has grown its
own copy of the logic and the two will diverge.

**Q12. Can I take the output to a cluster?**

Yes. The generated scripts are LF-only and valid shell, and the configs
reference their inputs **relatively** — the REINVENT config points at
`../03_triage/seeds.smi`, not at an absolute path from your machine, so the
whole pipeline directory can be copied.

**Q13. How long does a full run take, and what does it cost?**

Measured on a laptop with none of the heavy tools installed: `plan` about
2 seconds, three local stages about 17 seconds, all seven about 21 seconds.
The cost is in the cluster stages you then submit, and the generated scripts
state their own resource requests.

The reports are small by design — `--detail compact` is the default
everywhere, and a seven-stage `pipeline.json` is about 18 KB because it
holds verdicts, counts and paths rather than the contents of the files it
points at.

**Q14. Can I run just one stage, or use a different tool for one step?**

Yes to both. `--stages` takes any subset, and every stage writes ordinary
files — PDB, PDBQT, JSON, SMILES, TOML, MDP — so you can substitute your own
tool at any point and feed its output to the next stage. Nothing is locked
to this pipeline's own formats.

**Q15. Does it need the internet?**

Only `target` (when fetching by ID) and `sourcing` (ChEMBL, COCONUT). Give
`target` a local file with `--input-file` and supply compounds yourself, and
the whole pipeline runs offline.

**Q16. Which stage should I distrust most?**

`docking`, then `generative`. Both produce numbers that look authoritative
and are the least reliable things in the pipeline — a docking score is
triage, not affinity [7, 10], and a generative score measures the objective
you wrote, not activity [6]. `target` and `triage` are mechanical by
comparison: they fetch, clean and apply published thresholds.

**Q17. What do I do with the compounds at the end?**

They are hypotheses. The honest next steps:

1.  Look at them. A medicinal chemist's eye catches things no filter does.
2.  Check synthesisability properly — the SA score is a heuristic, and
    generated molecules are frequently not makeable [6].
3.  Validate the docking setup before trusting any ranking: re-dock the
    cognate ligand, and measure enrichment against actives and decoys
    [7, 8, 9].
4.  Run MD on your best few complexes to see whether the pose survives.
5.  Buy or synthesise and assay. Nothing above this line is evidence of
    activity.

--------------------------------------------------------------------------------

## What this pipeline is not

It is a **triage and hypothesis-generation** pipeline. Specifically:

- A docking score is not a binding free energy, and differences under about
  1 kcal/mol do not rank compounds [7, 10].
- An ADMET filter is a set of published physicochemical rules, not an assay
  [11, 12].
- A property-based generative objective produces molecules with the right
  properties, which is not the same as active molecules [6].
- A predicted structure is a hypothesis about a fold, and side chains in the
  pocket are its weakest part [3, 4, 5].
- One MD trajectory is one sample; conclusions need replicas [13].

Each stage reports its own caveats, and they are worth reading rather than
skipping to the numbers.

--------------------------------------------------------------------------------

## When something goes wrong

| Symptom | Cause | Fix |
|---|---|---|
| Exit 2 | Unknown stage name | The usage message lists the valid ones |
| `skipped (no compounds)` | Stage ran without its input | Run `sourcing` first, or pass compounds |
| A stage says NOT READY | The tool is absent | It still produces a cluster script |
| `plan` is slow or touches the network | It should not | Report it; `plan` is read-only by design |
| Scripts fail on the cluster | Paths or line endings | Q12 |

--------------------------------------------------------------------------------

## References

Full entries in [citation.bib](citation.bib), which collects the citations
from all four sub-skills.

1.  Anderson AC. The process of structure-based drug design. *Chem Biol*
    2003;10:787–797.
2.  Berman HM et al. The Protein Data Bank. *Nucleic Acids Res*
    2000;28:235–242.
3.  Jumper J et al. *Nature* 2021;596:583–589 — AlphaFold.
4.  Karelina M, Noh JJ, Dror RO. *eLife* 2023;12:RP89386 — binding modes
    from AlphaFold models.
5.  Terwilliger TC et al. *Nat Methods* 2024;21:110–116 — predictions as
    hypotheses.
6.  Gao W, Coley CW. *J Chem Inf Model* 2020;60:5714–5723 —
    synthesizability of generated molecules.
7.  Warren GL et al. *J Med Chem* 2006;49:5912–5931 — docking assessed.
8.  Mysinger MM et al. *J Med Chem* 2012;55:6582–6594 — DUD-E.
9.  Tran-Nguyen VK, Jacquemard C, Rognan D. *J Chem Inf Model*
    2020;60:4263–4273 — LIT-PCBA.
10. Chen YC. *Trends Pharmacol Sci* 2015;36:78–95 — beware of docking.
11. Lipinski CA et al. *Adv Drug Deliv Rev* 1997;23:3–25 — rule of five.
12. Leeson PD, Springthorpe B. *Nat Rev Drug Discov* 2007;6:881–890 —
    drug-like concepts in decision-making.
13. Knapp B, Ospina L, Deane CM. *J Chem Theory Comput* 2018;14:6127–6138 —
    replicas in molecular simulation.
