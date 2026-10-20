# Antigravity Skills for Biomedical Research

**Agentic skills for structure-based drug discovery** — tell
[Antigravity](https://github.com/google-deepmind/science-skills) what you need
in plain language, and it handles the computational pipeline from protein
structure to docked, filtered, and optimised compounds.

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://docs.astral.sh/ruff/)

---

## How It Works

These are **agent skills**, not command-line tools you type yourself. You prompt
Antigravity in natural language — describe what you need — and the agent reads
the skill instructions, picks the right one, runs the right commands, interprets
the output, and reports back.

You never need to remember a flag, a subcommand or a file format. The agent
does.

### Example prompts

> **"Prepare PDB entry 1ITF for docking"**
>
> Antigravity activates the **target-preparation** skill — fetches the
> structure from RCSB PDB, cleans it, converts to PDBQT, assesses quality,
> defines the binding site from the co-crystallised ligand, and reports the
> result with every decision recorded.

> **"Find known active compounds against DprE1 from M. tuberculosis"**
>
> Antigravity activates **compound-screening** — searches ChEMBL for the
> target, retrieves bioactivities filtered by potency, and returns the hits.

> **"Screen these compounds for drug-likeness and toxicity alerts"**
>
> Antigravity runs the **two-tier ADMET battery** — Lipinski, Veber, Egan,
> Ghose and Muegge rules at the appropriate strictness, plus PAINS, Brenk and
> aggregator alerts — and reports which compounds passed, which failed, and why.

> **"Dock the top 10 hits against the receptor we just prepared"**
>
> Antigravity calls AutoDock Vina with the receptor and binding site from the
> earlier step, ranks the results, and flags any score differences under 1
> kcal/mol as not meaningful.

> **"Design new molecules around the best scaffold using REINVENT"**
>
> Antigravity activates **compound-synthesis** — configures a REINVENT 4
> generative run with multi-objective scoring, prepares the seeds from the
> active compounds, and either runs locally or writes a cluster submission
> script.

> **"Set up an MD simulation for the top hit at body temperature"**
>
> Antigravity activates **md-simulation** — generates a four-stage GROMACS
> workflow (minimisation → NVT → NPT → production) with force-field-aware
> settings, and writes a SLURM submission script.

> **"Run the full drug discovery pipeline for PDB 6HEZ"**
>
> Antigravity activates the **drug-discovery-wizard** — orchestrates all four
> skills across seven stages, with ADMET screening on both sides of generative
> chemistry, and reports progress at each stage.

---

## The Five Skills

```
                         ┌────────────────────────┐
                         │  drug-discovery-wizard  │
                         │    (orchestrator)       │
                         └───────────┬────────────┘
                                     │ delegates to
             ┌───────────┬───────────┼───────────┬───────────┐
             ▼           ▼           ▼           ▼           ▼
     ┌──────────┐ ┌────────────┐ ┌────────┐ ┌────────┐ ┌────────┐
     │  target  │ │ compound   │ │compound│ │  md    │ │ ADMET  │
     │  prep    │→│ screening  │→│synthesis│→│  sim   │ │ tiers  │
     └──────────┘ └────────────┘ └────────┘ └────────┘ └────────┘
       PDB/AF2     ChEMBL/ZINC   REINVENT 4  GROMACS    RDKit
       ESMFold     COCONUT/Vina              MDP gen    2-tier
```

| Skill | When Antigravity uses it |
|---|---|
| **Target Preparation** | You mention a PDB ID, UniProt accession, or amino-acid sequence and want a receptor ready for docking |
| **Compound Screening** | You ask about compounds, actives, docking, binding affinity, drug-likeness or ADMET filtering |
| **Compound Synthesis** | You want new molecules designed or an existing lead optimised |
| **MD Simulation** | You want a molecular dynamics setup or a cluster submission script |
| **Drug Discovery Wizard** | You want the end-to-end pipeline — the wizard delegates to the four skills above |

Each skill's full instructions live in `.agents/skills/<skill>/SKILL.md`. The
agent reads these automatically when it selects a skill. You do not need to read
them yourself unless you want to understand what is happening under the hood.

### Design principles

- **Every decision is recorded.** The report names the source, the chains
  retained, each heterogen removed with its reason, the B-factor scale
  detected and the interpretation applied. A receptor can be reconstructed or
  audited from the report alone.
- **Ambiguity is refused, not guessed.** A gene name is rejected because it
  maps to many structures across many organisms. A structure whose provenance
  cannot be established is reported as indeterminate rather than graded on a
  scale that may not apply.
- **Nothing is organism-specific.** No code path inspects the source organism.
  The worked examples span human, viral, bacterial, fungal and protozoan
  targets to make that concrete, and a test keeps it that way.
- **Output size is bounded.** Compact output by default so the agent's context
  is not flooded. Verdicts are identical at every detail level; only the
  per-compound payload changes.

---

## Setup

### Install

Requires **Python 3.10+**. Clone this repository into your Antigravity
workspace:

```bash
git clone https://github.com/Walter-Odur/Antigravity-Skills-For-Biomedical-Research.git
cd Antigravity-Skills-For-Biomedical-Research

# Full install (core + docking + analysis + tests)
pip install -e ".[all]"

# Or core only (rdkit + biopython)
pip install -e .
```

### Verify

```bash
pytest
```

The test suite needs no network access and no external binaries. It should
pass in under a minute.

### What is optional

The core runs on `rdkit` and `biopython` alone. Everything else degrades
gracefully — the agent gets a clear message about what is missing and suggests
alternatives:

| Tool | Needed for | Without it |
|---|---|---|
| `meeko`, `scipy`, `gemmi`, or Open Babel | PDBQT conversion | A cleaned PDB is still produced; the agent reports a converter is missing |
| AutoDock Vina | Docking | The agent suggests installing it or using a cluster script |
| GROMACS | *Running* an MD workflow | File generation is unaffected — it produces files, not simulations |
| REINVENT 4 + PyTorch + priors | Generative chemistry | Config generation, seed preparation and analysis still work |
| Internet access | Database queries, ESMFold | Offline capabilities are unaffected |

<details>
<summary><strong>Installing optional tools</strong></summary>

```bash
# PDBQT conversion
pip install meeko scipy gemmi

# AutoDock Vina (not pip-installable on Windows — needs Boost headers)
conda install -c conda-forge vina
# or: https://github.com/ccsb-scripps/AutoDock-Vina/releases

# GROMACS
conda install -c conda-forge gromacs

# REINVENT 4
git clone https://github.com/MolecularAI/REINVENT4
pip install -e REINVENT4
export REINVENT_DIR=$PWD/REINVENT4
# Prior models: https://doi.org/10.5281/zenodo.15641296

# Diagnostic figures
pip install "agskills[plots]"
```

</details>

---

## What Each Skill Can Do

### Target Preparation

Takes a protein target from a database identifier to a docking-ready receptor
with a defined binding site.

**Just ask Antigravity:**
- *"Prepare PDB 6HEZ chain A for docking"*
- *"Fetch the AlphaFold model for UniProt P9WJA5 and assess its quality"*
- *"I have this sequence — can you predict its structure and tell me if it's good enough for docking?"*
- *"What ligands are in PDB 4EY7? Which one should I use for the binding site?"*

**What happens behind the scenes:**

| Capability | Detail |
|---|---|
| Fetch structures | RCSB PDB by ID, AlphaFold DB by UniProt accession (resolves current model version), local files |
| Clean receptors | Removes water, cryoprotectants, buffer ions; keeps catalytic metals; records every removal with its reason |
| Preserve provenance | Retains `HEADER`, `TITLE`, `EXPDTA`, `CRYST1`, `REMARK 2 RESOLUTION` so the file stays self-describing |
| Convert to PDBQT | Meeko preferred (better charges), Open Babel fallback; warns about approximated side chains *inside the docking box* |
| Assess quality | Detects the B-factor scale (pLDDT 0–100, pLDDT 0–1, or experimental B-factors) from provenance, not magnitude; issues PROCEED / REFINE / STOP verdicts for predicted models |
| Define binding site | From a co-crystallised ligand (gold standard), pocket residues, or whole-structure blind docking |
| Cofactor awareness | Ranks heterogens by druggability, not size — FAD, heme and NADH are flagged, not selected |
| Interface warning | Reports when a chain contributing to the pocket is removed (the obligate-oligomer trap) |
| Structure prediction | ESMFold for ≤400 residues; generates SLURM scripts for ColabFold / AlphaFold2 for longer sequences |

**The B-factor problem.** One column of a coordinate file carries different
quantities depending on its source, and they run in opposite directions. A
well-refined crystal has B-factors of 20–50; so does a poor predicted model's
pLDDT. The skill reads provenance from the header rather than guessing from
magnitude, because a magnitude rule rejects good crystal structures and
approves bad predictions. The ESMFold API writes pLDDT between 0 and 1; the
skill detects this and rescales, so a model with mean 0.84 is correctly read as
84 pLDDT rather than being called unusable.

### Compound Screening

Searches compound libraries, docks molecules, and screens them on ADMET
properties.

**Just ask Antigravity:**
- *"Find actives against CHEMBL3804751 with pChEMBL above 7"*
- *"Search COCONUT for natural products with a benzimidazole core"*
- *"Screen these SMILES for drug-likeness"*
- *"Dock these 20 compounds against the receptor and rank them"*

**What happens behind the scenes:**

| Capability | Detail |
|---|---|
| ChEMBL queries | Target search by name or accession; bioactivity retrieval filtered by pChEMBL, assay type and organism |
| COCONUT | Natural product search by substructure, molecular formula or text |
| ZINC | Verified bulk download URLs for purchasable compound subsets |
| Docking | AutoDock Vina with automatic box placement from the binding site; re-docking validation |
| ADMET screening | Two-tier battery: *relaxed* (lead-like) and *strict* (drug-like). Published rules always report their authors' thresholds — a tier changes which rules must pass, never redefines a violation |
| Structural alerts | PAINS, Brenk and aggregator filters — flagging is distinct from rejecting |
| Molecular descriptors | MW, LogP, TPSA, HBD, HBA, rotatable bonds, QED, synthetic accessibility, Fsp3, aromatic rings |

### Compound Synthesis

Generates and optimises novel molecules with REINVENT 4.

**Just ask Antigravity:**
- *"Design new molecules around the best scaffold from the docking results"*
- *"Do scaffold hopping on this lead compound"*
- *"Set up a REINVENT run optimising for docking score and synthetic accessibility"*

**What happens behind the scenes:**

- **De novo design** — unrestricted generation from a prior
- **Scaffold hopping** — constrained to a core with varied decorations
- **R-group replacement** — systematic substitution at defined positions
- **Linker design** — connecting two fragments with generated linkers
- **Transfer learning** — fine-tuning the prior on known actives
- **Multi-objective RL** — simultaneous optimisation of docking score, QED,
  synthetic accessibility and custom scoring components
- **Scaffold-diversity diagnostics** — Murcko decomposition and Tanimoto
  clustering on REINVENT output

### MD Simulation

Generates GROMACS workflows and validated MDP files. This skill **produces
files; it does not run simulations** — production MD is GPU work that belongs in
a batch job.

**Just ask Antigravity:**
- *"Set up an MD simulation for this protein–ligand complex using CHARMM36"*
- *"Generate a SLURM submission script for a 100 ns production run"*
- *"Validate this MDP file against the AMBER force field"*

**What happens behind the scenes:**

- **Force-field-aware MDP generation** — CHARMM36 requires a force-switched
  Lennard-Jones potential and no dispersion correction; AMBER requires a plain
  cutoff with one. They are not interchangeable, and the skill enforces this.
- **Four-stage workflow** — energy minimisation, NVT equilibration, NPT
  equilibration, production MD, each with correct settings
- **MDP validation** — cross-checks an existing MDP against a declared force
  field
- **SLURM script generation** — for GROMACS MD, REINVENT, AutoDock Vina batch
  screening, ColabFold and AlphaFold2

### Drug Discovery Wizard

Orchestrates the four skills above across seven pipeline stages. It
**delegates; it does not reimplement** — every stage calls the same code the
individual skills call, so a correction reaches all callers.

**Just ask Antigravity:**
- *"Run the full drug discovery pipeline for PDB 6HEZ"*
- *"What can this machine run? Show me the readiness plan"*
- *"Resume the pipeline from the docking stage"*

**The seven stages:** target preparation → compound sourcing → relaxed ADMET
triage → docking → strict ADMET filtering → generative chemistry → MD setup

---

## Diagnostic Figures

Every subcommand that writes a report can draw one alongside it:

| Subcommand | What the figure shows |
|---|---|
| `prepare-receptor` | Per-residue trace, cleaning counts, ligand ranking with cofactors coloured apart, and the docking box drawn over the structure with the selected ligand in red |
| `assess-structure` | Per-residue confidence with the pLDDT bands shaded and low-confidence regions marked, beside the band counts |
| `admet-filter` | Outcome funnel, which checks compounds failed, how close compounds came to passing, and descriptor spread against the tier bounds |
| `dock` | Ranked affinities with the 1 kcal/mol band marked, and pose spread within one compound |
| `ag-wizard run` | Stage-by-stage progress with status |

Two of these catch most problems. **Which checks compounds failed** makes a
single check killing an entire set obvious at a glance. **The docking box
panel** shows the box sitting on the ligand with receptor around it, which
is the question most worth answering before committing a screen.

Each figure is built from the report the subcommand already wrote, so it
cannot disagree with the numbers. An advisory check is drawn in amber and
labelled as not rejecting, because presenting it as a rejection would
invert the meaning of the tier.

matplotlib is optional. Without it a run proceeds normally and reports that
the figure was skipped.

---

## Configuration

All optional. Nothing is hardcoded to a particular machine.

| Variable | Purpose |
|---|---|
| `AGSKILLS_USER_AGENT` | User-Agent sent to public APIs. Set to something that identifies you |
| `AGSKILLS_CACHE_DIR` | Where rate-limiter state is kept. Defaults to a per-user temp directory |
| `REINVENT_DIR` | The REINVENT 4 checkout |
| `REINVENT_PRIOR_BASE` | Where the `.prior` model files live |

---

## Repository Layout

```
pyproject.toml              package metadata, dependencies, pytest config
README.md                   this file
LICENSE                     MIT
src/agskills/               the shared implementation
  errors.py                 typed exceptions, each with an exit code
  io_utils.py               atomic UTF-8 / LF writes, filename safety
  http.py                   rate-limited HTTP client with bounded deadline
  hpc.py                    SLURM script generation
  chem/                     descriptors, published rules, alerts, ADMET tiers
    admet.py                two-tier ADMET battery
    alerts.py               PAINS, Brenk, aggregator filters
    descriptors.py          molecular property calculations
    rules.py                Lipinski, Veber, Egan, Ghose, Muegge
    smiles.py               SMILES handling and canonicalisation
  targets/                  fetch, clean, confidence, PDBQT, binding sites
    fetch.py                PDB / AlphaFold DB retrieval
    clean.py                receptor cleaning
    quality.py              B-factor scale detection and pLDDT assessment
    pdbqt.py                PDBQT conversion via Meeko or Open Babel
    site.py                 binding-site definition and cofactor ranking
    esmfold.py              ESMFold structure prediction
  libraries/                compound database clients
    chembl.py               ChEMBL target and bioactivity queries
    coconut.py              COCONUT natural product search
    zinc.py                 ZINC bulk download URL generation
  docking/                  molecular docking
    vina.py                 AutoDock Vina interface
  md/                       molecular dynamics setup
    mdp.py                  MDP generation and validation
    workflow.py             GROMACS four-stage workflow
  reinvent/                 generative chemistry
    config.py               REINVENT 4 TOML config generation
    priors.py               prior model management
    scoring.py              scoring component setup
    runner.py               REINVENT execution
    analyze.py              output analysis and diversity diagnostics
  cli/                      one entry point per skill
.agents/skills/             the skill bundle (what Antigravity reads)
  plugin.json               bundle manifest (v2.1.0)
  SKILL_LICENSES.md         third-party data source terms
  <skill>/SKILL.md          each skill's instructions for the agent
  <skill>/scripts/          PEP 723 launcher scripts
  <skill>/references/       citation.bib, questions.md, best_practices.md
tests/                      21 test modules, real-data fixtures
documentation/              detailed per-skill documentation (DOCX)
verification/<skill>/       hand-written verification protocols
tools/                      generators for derived documentation
```

---

## Testing

```bash
pytest                       # hermetic — no network, no external binaries
pytest -v                    # with test names
pytest -m network            # opt in to live API calls
pytest -m requires_vina      # opt in to real docking
pytest --cov=agskills        # coverage
```

The default run deselects anything needing the network, an external binary, or
more than a few seconds. Capability-marked tests skip cleanly when the tool is
absent.

Tests run against **real data**, not synthetic stand-ins: two experimental
structures from the PDB, a real AlphaFold model, independently computed
reference properties for ten marketed drugs, a REINVENT output CSV in the real
column format, and authentic AutoDock Vina output.

### Question indexes

Each skill ships a question index mapping the questions a newcomer actually asks
onto the command that answers them — 106 questions across the five skills, each
answer measured rather than asserted.

| Index | Questions | Covers |
|---|---:|---|
| [target_preparation](.agents/skills/target_preparation/references/questions.md) | 23 | structures, binding sites, cleaning choices, predicted models |
| [compound_screening](.agents/skills/compound_screening/references/questions.md) | 25 | libraries, pChEMBL, tiers, alerts, docking scores |
| [compound_synthesis](.agents/skills/compound_synthesis/references/questions.md) | 21 | generators, transforms, mode collapse, reading results |
| [md_simulation](.agents/skills/md_simulation/references/questions.md) | 20 | force fields, the four stages, walltime, post-run checks |
| [drug_discovery_wizard](.agents/skills/drug_discovery_wizard/references/questions.md) | 17 | stage selection, delegation, what to distrust |

### Verification protocols

Each skill has a hand-written verification protocol under `verification/`:

| Protocol | Cases |
|---|---|
| [target_preparation](verification/target_preparation/TEST.md) | 14 — fetch, clean, PDBQT, three B-factor scales, binding sites, exit codes |
| [compound_screening](verification/compound_screening/TEST.md) | 13 — descriptors vs PubChem, tier semantics, alerts, ChEMBL, ZINC, output size |
| [compound_synthesis](verification/compound_synthesis/TEST.md) | 12 — TOML schema, transform directions, de-duplication, mode collapse |
| [md_simulation](verification/md_simulation/TEST.md) | 11 — force-field non-bonded settings, staging, MDP validation, `-maxh` |
| [drug_discovery_wizard](verification/drug_discovery_wizard/TEST.md) | 9 — per-stage readiness, delegation, deferral, relocatable artefacts |

Every expected value was measured by running the command, including places
where the skill *should* disagree with an outside reference.

---

## Scientific Conventions

Worth knowing before interpreting any output:

- **Molecular weight** is the average mass, as Lipinski's rule and PubChem both
  mean. The monoisotopic mass is reported separately as `exact_mw`.
- **A published rule always reports its authors' thresholds.** A screening tier
  changes which rules must pass; it never redefines a violation.
- **Gating against advisory.** Published risk *flags* (e.g. the Pfizer 3/75
  criterion) are reported but do not reject a compound.
- **Non-bonded settings belong to the force field.** CHARMM36 requires a
  force-switched LJ potential and no dispersion correction; AMBER requires a
  plain cutoff with one.
- **A docking score is triage.** Differences under ~1 kcal/mol are not
  meaningful. Correlation with measured affinity across unrelated chemotypes is
  weak. Scoring functions predict poses better than affinities, and enrichment
  better than either.
- **pLDDT is not a B-factor.** For a predicted model, higher is better on a
  0–100 scale. For a crystal structure, the column holds atomic displacement
  parameters — lower is better, unbounded. The skills detect which they are
  looking at.

---

## Output and Context Cost

Every subcommand writes structured JSON to a file. The agent reads the fields
it needs rather than printing the whole thing into the conversation.

| Control | Effect |
|---|---|
| `--detail compact` (default) | Verdicts + headline numbers. ~160 tokens per compound |
| `--detail full` | Every check, parameter and descriptor. ~2,090 tokens per compound |

Verdicts, counts and ordering are identical at every detail level.

---

## Citations

The methods these skills wrap, and the thresholds they apply, come from:

- Lipinski CA et al. *Adv Drug Deliv Rev* 1997;23:3–25
- Veber DF et al. *J Med Chem* 2002;45:2615–2623
- Egan WJ et al. *J Med Chem* 2000;43:3867–3877
- Ghose AK et al. *J Comb Chem* 1999;1:55–68
- Muegge I et al. *J Med Chem* 2001;44:1841–1846
- Hughes JD et al. *Bioorg Med Chem Lett* 2008;18:4872–4875
- Baell JB, Holloway GA. *J Med Chem* 2010;53:2719–2740 (PAINS)
- Brenk R et al. *ChemMedChem* 2008;3:435–444
- Ertl P, Schuffenhauer A. *J Cheminform* 2009;1:8 (synthetic accessibility)
- Bickerton GR et al. *Nat Chem* 2012;4:90–98 (QED)
- Anderson AC. *Chem Biol* 2003;10:787–797 (SBDD process)
- Berman HM et al. *Nucleic Acids Res* 2000;28:235–242 (Protein Data Bank)
- Trott O, Olson AJ. *J Comput Chem* 2010;31:455–461 (AutoDock Vina)
- Eberhardt J et al. *J Chem Inf Model* 2021;61:3891–3898 (Vina 1.2)
- Huang J, MacKerell AD. *J Comput Chem* 2013;34:2135–2145 (CHARMM36)
- Bernetti M, Bussi G. *J Chem Phys* 2020;153:114107 (C-rescale barostat)
- Lin Z et al. *Science* 2023;379:1123–1130 (ESMFold)
- Jumper J et al. *Nature* 2021;596:583–589 (AlphaFold)
- Varadi M et al. *Nucleic Acids Res* 2024;52:D368–D375 (AlphaFold DB)
- Zdrazil B et al. *Nucleic Acids Res* 2024;52:D1180–D1192 (ChEMBL)
- Sorokina M et al. *J Cheminform* 2021;13:2 (COCONUT)
- Loeffler HH et al. *J Cheminform* 2024;16:20 (REINVENT 4)
- Warren GL et al. *J Med Chem* 2006;49:5912–5931 (docking assessment)
- Feinstein WP, Brylinski M. *J Cheminform* 2015;7:18 (box sizing)
- Madhavi Sastry G et al. *J Comput-Aided Mol Des* 2013;27:221–234 (preparation influence)
- Olsson MHM et al. *J Chem Theory Comput* 2011;7:525–537 (PROPKA3)

Individual thresholds carry their citation in the JSON output, so a result can
be traced to its source. Each skill also ships `references/citation.bib`.

In total the bundle carries **171 citation entries**, every one with a DOI
resolved through Crossref and checked against the registered title, journal,
volume and pages.

---

## Third-Party Licenses

The skills themselves are [MIT licensed](LICENSE). The tools and data sources
they reach carry their own terms — see
[SKILL_LICENSES.md](.agents/skills/SKILL_LICENSES.md) for the full table.

---

## Author

**Walter Odur**

Built during the [Google DeepMind Antigravity Workshop](https://github.com/google-deepmind/science-skills), 2026.
