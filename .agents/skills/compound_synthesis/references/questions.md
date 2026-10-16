# Questions this skill answers

Find your question, run the command. The generators and scoring components
are chemistry-general; none of them is tied to a target class or organism.

## Start here

| You want | Run |
|---|---|
| To know whether this machine can run REINVENT | `check-setup` |
| To see the scoring profiles and generators | `list-profiles` |
| Brand-new molecules from scratch | `generate-config --mode staged_learning --generator reinvent` |
| To decorate a scaffold you already have | `--generator libinvent` |
| To join two fragments | `--generator linkinvent` |
| Molecules similar to a known active | `--generator mol2mol` |
| Peptides | `--generator pepinvent` |
| To bias a model toward your own series | `--mode transfer_learning` |
| To score molecules you already have | `--mode scoring` or `score-molecules` |
| To check a config before queueing it | `preflight` |
| To understand what a finished run produced | `analyze-results` |

--------------------------------------------------------------------------------

## Getting started

**Q1. What does this skill actually do? Does it invent drugs?**

It runs REINVENT 4 [1], a generative model that proposes molecules scoring
well against an objective you define. It proposes *structures*; it does not
establish that they bind, that they are active, or that anyone can make
them. Everything it produces is a hypothesis for screening and synthesis.

**Q2. Do I need a GPU?**

For `run` and `score-molecules`, yes in practice. Everything else —
`check-setup`, `list-profiles`, `generate-config`, `prepare-seeds`,
`preflight`, `analyze-results`, `generate-hpc-script` — works on a laptop
with nothing installed but Python. Design the run locally, execute it on a
cluster.

**Q3. Nothing works. What is missing?**

```bash
uv run scripts/compound_synthesis_api.py check-setup --output setup.json
```

It reports each failure separately with a remedy, and exits 0 because
reporting that a machine is not ready is a successful diagnosis:

```
REINVENT 4 setup: NOT READY
  PROBLEM: The 'reinvent' package is not importable.
  PROBLEM: 10 of 10 prior model files are missing.
  PROBLEM: PyTorch is not installed.
  -> git clone https://github.com/MolecularAI/REINVENT4 && pip install -e REINVENT4
  -> Download the priors from https://doi.org/10.5281/zenodo.15641296
  -> Install the PyTorch build matching your CUDA version
```

**Q4. Which generator do I want?**

| Generator | Input | Produces |
|---|---|---|
| `reinvent` | nothing | complete molecules, de novo |
| `libinvent` | a scaffold with `*` attachment points | R-group variations on it |
| `linkinvent` | two fragments with `*` | a linker joining them |
| `mol2mol` | one or more molecules | close analogues |
| `pepinvent` | peptide input | peptide variants |

If you have a lead and want to improve it, `mol2mol` or `libinvent`. If you
have nothing but a target, `reinvent`.

**Q5. What is "staged learning" versus "sampling"?**

| Mode | Does |
|---|---|
| `sampling` | draws molecules from the prior; no optimisation |
| `scoring` | scores molecules you supply, no generation |
| `staged_learning` | reinforcement learning against your objective — the usual choice |
| `transfer_learning` | retrains the prior on your compounds first |

`staged_learning` is reinforcement learning as introduced for this problem
in [2] and carried forward in REINVENT [1, 3].

--------------------------------------------------------------------------------

## Setting up a run

**Q6. Which scoring profile should I pick?**

```bash
uv run scripts/compound_synthesis_api.py list-profiles --output profiles.json
```

| Profile | Components | For |
|---|---|---|
| `drug-like` | 9 | oral small-molecule space |
| `lead-like` | 8 | leaving room for later optimisation |
| `fragment-like` | 6 | rule-of-three fragment space [4] |
| `kinase-inhibitor` | 9 | ATP-competitive kinase space |
| `anti-tb` | 7 | higher lipophilicity, as mycobacterial actives occupy |
| `cns` | 8 | CNS-penetrant space, after the CNS MPO consensus [5] |

These are **property profiles, not activity models**, and the ones that
could be mistaken for activity models say so. The `anti-tb` caveat reads:

> A property profile is not an activity model. Whole-cell potency depends on
> envelope permeability and efflux, which no physicochemical descriptor
> captures. Add a target-based component (docking through DockStream, or a
> QSAR model through ChemProp or Qptuna) before treating output as leads.

The same holds for every profile. Property optimisation alone converges on
molecules that look right and do nothing.

**Q7. My target isn't in the list. Is the skill useless to me?**

No — the profiles are property windows, not organism presets. Pick the one
matching the *physicochemical space* you want:

- A human intracellular enzyme → `drug-like` or `lead-like`.
- A kinase of any organism → `kinase-inhibitor`.
- A viral protease → `drug-like`; many antivirals are peptidomimetic and
  larger, so check the MW window suits your series.
- A Gram-negative bacterial target → none of these fits well. Compounds
  that cross the outer membrane are more polar and more charged than oral
  drugs [11, 12], so an oral-drug profile selects against them. Widen the
  TPSA and MW windows and add a target-based component.
- A fungal or protozoan intracellular target → `drug-like` is a reasonable
  starting point.

Then add a component that encodes **activity**, which is the part no profile
supplies.

**Q8. How do I know the profiles are genuinely different?**

Compare the windows. Measured:

```
drug-like  logP window 1 to 4
anti-tb    logP window 2 to 6
```

The `anti-tb` window is shifted up because known actives there sit at higher
lipophilicity — bedaquiline is 555 Da with a measured clogP of 7.13. A
profile that merely renamed `drug-like` would be worse than no profile.

**Q9. Generate me a config.**

```bash
uv run scripts/compound_synthesis_api.py generate-config \
    --mode staged_learning --generator reinvent \
    --scoring-profile drug-like --num-steps 300 --output reinvent.toml
```

It tells you the cost up front: 300 steps × batch 64 = **19,200 molecule
evaluations**. That is the number that decides whether a run fits a
walltime.

**Q10. How many steps do I need?**

Start at 100–300 and look at the score trajectory before spending more. A
run that has plateaued is not improved by more steps; it is more likely to
collapse (Q17). Sample efficiency in this setting is an active research
question, not a solved one [6].

**Q11. What are "transforms" and why should I care?**

They map a raw descriptor onto a 0–1 reward, and the **direction** matters
more than anything else in the config:

| Component | Transform | Because |
|---|---|---|
| MolecularWeight, SlogP, TPSA | `double_sigmoid` | a **window** — too low is as bad as too high |
| NumRotBond, SAScore | `reverse_sigmoid` | **lower is better** |
| QED | none | already 0–1, higher better |
| custom_alerts | none | binary |

`SAScore` with a plain `sigmoid` instead of `reverse_sigmoid` would reward
*harder*-to-make molecules, and every log would look healthy while the run
converged on unmakeable chemistry. The generated configs use weight 1.5 on
SAScore — the highest in the profile — because synthesisability outranks any
single property window [7].

**Q12. Aggregation: geometric or arithmetic mean?**

Geometric, in every profile here. With a geometric mean one component
scoring zero zeroes the total, so a molecule cannot buy its way past an
unmet requirement by excelling elsewhere. An arithmetic mean lets a
beautifully drug-like molecule with a reactive warhead score well.

**Q13. How do I prepare seeds?**

```bash
uv run scripts/compound_synthesis_api.py prepare-seeds \
    --smiles-file actives.smi --output seeds.smi
```

De-duplication is **chemical, not textual** — `CCO` and `OCC` are the same
molecule and collapse to one. A textual de-duplication keeps both and
quietly biases the run toward that structure. Invalid structures are dropped
with their reason.

For `libinvent` and `linkinvent`, seeds need attachment points written as
`*`. The error message says so if you forget.

**Q14. Check my config before I queue it.**

```bash
uv run scripts/compound_synthesis_api.py preflight \
    --config reinvent.toml --output pre.json
```

It lists **every** missing path with the config key referencing it, and
exits 6. `run` performs the same check and refuses identically, so a job
cannot reach a scheduler only to die after loading — which costs a queue
wait to discover.

**Q15. Submit it to a cluster.**

```bash
uv run scripts/compound_synthesis_api.py generate-hpc-script \
    --config reinvent.toml --output submit.sh
sbatch submit.sh
```

The script activates conda through `eval "$(conda shell.bash hook)"`,
because a non-interactive SLURM shell cannot `conda activate` otherwise, and
checks `reinvent` is on PATH before the long step.

--------------------------------------------------------------------------------

## Reading the results

**Q16. My run finished. What do I look at?**

```bash
uv run scripts/compound_synthesis_api.py analyze-results \
    --csv staged_learning_1.csv --top-n 25 --output top.json
```

**Q17. It generated 96 molecules — or did it?**

It did not. Measured on real REINVENT output:

```
rows 96  unique molecules 13  duplicate samples 83
```

REINVENT samples **with replacement**, so the row count is a sampling count.
Reporting "96 compounds generated" overstates the yield by sevenfold. The
per-molecule sample count is shown as `(x6)`, `(x7)` — a molecule sampled
seven times is a converged mode, which is information about the run.

**Q18. What is mode collapse and how do I spot it?**

The agent finds one scaffold that scores well and only decorates it. You get
many molecules and one idea. It is a recognised failure mode of generative
molecular design [8, 9].

The diagnosis reports Bemis–Murcko frameworks [10]:

```
scaffolds: 1 distinct (0.2)
diversity: 1 scaffold(s) across 5 molecules (20.0%). Every molecule shares a
  single Bemis-Murcko framework: the agent has collapsed onto one scaffold
  and is only decorating it. Enable the diversity filter, or raise its
  bucket size if it is already on.
```

Note the ratio alone cannot diagnose this: 1 scaffold over 5 molecules and
20 scaffolds over 100 molecules are both 0.20, and only the first is
collapse. The absolute count is checked first, and short runs are told so
rather than graded on a percentage.

Remedies: enable or strengthen the diversity filter, lower `sigma`, add a
diversity or similarity-penalty component, or stop the run earlier.

**Q19. The scores are high but the molecules look absurd.**

This is the expected outcome of optimising a property objective without an
activity term, and the central caveat of the field: generative models
readily produce molecules that satisfy the score and are either
unsynthesisable or biologically meaningless [7, 8]. Benchmarks exist
precisely because high scores are easy [9].

Checks worth running:

1.  Screen the output through `ag-compound-screening admet-filter
    --strictness strict` — an independent battery the generator was not
    optimising against.
2.  Look at the SA score distribution, and read [7] on why generated
    molecules are often not makeable.
3.  Show them to a synthetic chemist before promising anything.
4.  Add a structure-based or QSAR component so activity, not just
    properties, is being optimised.

**Q20. My scores never improved.**

`analyze-results` reports the trajectory and says which case it is: the mean
rose, fell, or stayed flat. Flat usually means the objective is unreachable,
`sigma` is too low, or the agent converged before the run started. A falling
mean points at the scoring definition or the learning rate.

**Q21. Nothing passes the strict ADMET tier afterwards.**

Normal, and not a failure. The relaxed tier is the seed tier and the strict
tier is the final-candidate tier; measured on a real pipeline, 25 actives
give 25/25 at relaxed and 5/25 at strict. If *nothing* passes, check your
scoring profile's windows against the strict tier's windows — you may be
optimising into a region the strict tier excludes by design.

--------------------------------------------------------------------------------

## When something goes wrong

| Symptom | Cause | Fix |
|---|---|---|
| Exit 4 | REINVENT or PyTorch missing | Q3 |
| Exit 6 on `preflight`/`run` | Prior files or seeds missing | Q14, then `check-setup` |
| `The libinvent generator needs seed structures` | No `--smiles-file` | Q13; attachment points are `*` |
| Exit 2 | Bad flag combination | The message names what is wrong |
| Config loads but REINVENT errors on paths | Absolute paths from another machine | Paths are now written relative to the config; regenerate it |
| Run dies at load after a queue wait | Config never preflighted | Q14 |
| `0 valid molecules` in output | Prior and generator mismatched | Check `--generator` against the prior file |

--------------------------------------------------------------------------------

## References

Full entries in [citation.bib](citation.bib). Scoring components are
documented in [scoring_components.md](scoring_components.md), and runnable
configs are in [../examples/](../examples/).

1.  Loeffler HH et al. *J Cheminform* 2024;16:20 — REINVENT 4.
2.  Olivecrona M et al. *J Cheminform* 2017;9:48 — molecular de novo design
    through deep reinforcement learning.
3.  Blaschke T et al. *J Chem Inf Model* 2020;60:5918–5922 — REINVENT 2.0.
4.  Congreve M et al. *Drug Discov Today* 2003;8:876–877 — rule of three.
5.  Wager TT et al. *ACS Chem Neurosci* 2010;1:435–449 — CNS MPO.
6.  Thomas M et al. *J Cheminform* 2022;14:68 — sample efficiency in de novo
    generation.
7.  Gao W, Coley CW. *J Chem Inf Model* 2020;60:5714–5723 — the
    synthesizability of molecules proposed by generative models.
8.  Renz P et al. *Drug Discov Today Technol* 2019;32–33:55–63 — failure
    modes in molecule generation and optimisation.
9.  Brown N et al. *J Chem Inf Model* 2019;59:1096–1108 — GuacaMol.
10. Bemis GW, Murcko MA. *J Med Chem* 1996;39:2887–2893 — molecular
    frameworks.
11. O'Shea R, Moser HE. *J Med Chem* 2008;51:2871–2878 — antibacterial
    physicochemical properties.
12. Richter MF et al. *Nature* 2017;545:299–304 — accumulation rules for
    Gram-negative bacteria.
13. Ertl P, Schuffenhauer A. *J Cheminform* 2009;1:8 — synthetic
    accessibility.
