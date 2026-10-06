# Antigravity Drug Discovery Skills

Five agent skills for structure-based drug discovery, sharing one tested
implementation.

| Skill | Console command | Does |
|---|---|---|
| `target-preparation` | `ag-target-preparation` | Fetch, clean and assess protein structures; define docking sites |
| `compound-screening` | `ag-compound-screening` | Search ChEMBL, COCONUT and ZINC; dock with Vina; two-tier ADMET screening |
| `compound-synthesis` | `ag-compound-synthesis` | REINVENT 4 generative chemistry |
| `md-simulation` | `ag-md-simulation` | GROMACS workflows and cluster submission scripts |
| `drug-discovery-wizard` | `ag-wizard` | Orchestrate the four above end to end |

Each skill's `SKILL.md` lives under `.agents/skills/<skill>/` and is the
reference for using it. The science lives in one package, `agskills`, under
`src/`.

## Install

Requires Python 3.10 or newer.

```bash
# From the workshop root
pip install -e ".[all]"
```

That installs the package plus every optional extra (docking, analysis and
test dependencies). For the core only:

```bash
pip install -e .
```

If you use `uv` rather than `pip`:

```bash
uv pip install -e ".[all]"
```

### Verify the installation

```bash
pytest
ag-compound-screening describe-tiers --output /tmp/tiers.json
```

The test suite should report all tests passing in well under a minute. It
needs no network access and no external binaries.

## What is optional

The core runs on `rdkit` and `biopython` alone. Everything else degrades
with a clear message rather than failing obscurely:

| Tool | Needed for | Without it |
|---|---|---|
| `meeko`, `scipy`, `gemmi`, or Open Babel | PDBQT conversion | A cleaned PDB is still produced, and the report says a converter is missing |
| AutoDock Vina | `dock` | Report says how to install it; use a cluster script instead |
| GROMACS | *Running* an MD workflow | Generation is unaffected; it produces files |
| REINVENT 4 + PyTorch + priors | `compound-synthesis run` | Config generation, seed preparation and result analysis all still work |
| Internet access | Database queries, ESMFold | Offline subcommands are unaffected |

Check what is available:

```bash
ag-compound-synthesis check-setup --output setup.json
ag-wizard plan --pdb-id 6HEZ --output plan.json
```

`plan` reports, per pipeline stage, whether this machine can run it and
what is missing.

### Installing the optional tools

```bash
# PDBQT conversion (cross-platform, via pip)
pip install meeko scipy gemmi

# AutoDock Vina. Not pip-installable on Windows: the vina module needs
# Boost headers. Use conda, or download a binary release.
conda install -c conda-forge vina
# or https://github.com/ccsb-scripps/AutoDock-Vina/releases

# GROMACS
conda install -c conda-forge gromacs

# REINVENT 4
git clone https://github.com/MolecularAI/REINVENT4
pip install -e REINVENT4
export REINVENT_DIR=$PWD/REINVENT4
# Prior models: https://doi.org/10.5281/zenodo.15641296
```

## Running the skills

Two equivalent ways. The console command, once installed:

```bash
ag-target-preparation prepare-receptor --pdb-id 6HEZ --chain A \
    --output work/receptor.json
```

Or the skill script directly, which works from a checkout even before
installing:

```bash
python .agents/skills/target-preparation/scripts/target_preparation_api.py \
    prepare-receptor --pdb-id 6HEZ --chain A --output work/receptor.json
```

`--output` is required for every subcommand of every skill. Read the file
it writes; the console summary is a convenience, not the result.

## A short worked example

```bash
# 1. Receptor, binding site and a confidence report, in one step
ag-target-preparation prepare-receptor \
    --pdb-id 6HEZ --chain A --output work/receptor.json

# 2. Known actives for the target
ag-compound-screening query-chembl \
    --target-id CHEMBL3804751 --pchembl-min 7 --output work/actives.json

# 3. Triage them, keeping a diverse set
ag-compound-screening admet-filter \
    --smiles "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC" --names hit1 \
    --strictness relaxed --output work/triage.json

# 4. A generative run, configured and ready to submit
ag-compound-synthesis generate-config \
    --mode staged_learning --generator reinvent \
    --scoring-profile anti-tb --num-steps 300 \
    --output work/reinvent.toml

# 5. An MD setup at body temperature
ag-md-simulation gmxapi-setup \
    --pdb work/6HEZ_clean.pdb --ff charmm36 --temperature 310 \
    --production-ns 100 --output work/md/run_md.py
```

Or let the orchestrator sequence it:

```bash
ag-wizard run --pdb-id 6HEZ --chain A \
    --chembl-target CHEMBL3804751 \
    --stages target sourcing triage \
    --workdir pipeline --output pipeline.json
```

## Configuration

All optional. Nothing is hardcoded to a particular machine.

| Variable | Purpose |
|---|---|
| `AGSKILLS_USER_AGENT` | The User-Agent sent to public APIs. Set it to something that identifies you and includes a contact address |
| `AGSKILLS_CACHE_DIR` | Where rate-limiter state is kept. Defaults to a per-user temporary directory |
| `REINVENT_DIR` | The REINVENT 4 checkout |
| `REINVENT_PRIOR_BASE` | Where the `.prior` model files live |

## Layout

```
pyproject.toml              package metadata, dependencies, pytest config
README.md                   this file
src/agskills/               the shared implementation
  errors.py                 typed exceptions, each with an exit code
  io_utils.py               atomic UTF-8 / LF writes, filename safety
  http.py                   rate-limited HTTP client with a bounded deadline
  hpc.py                    SLURM script generation
  chem/                     descriptors, published rules, alerts, ADMET tiers
  targets/                  fetch, clean, confidence, PDBQT, binding sites
  libraries/                ChEMBL, COCONUT, ZINC
  docking/                  AutoDock Vina
  md/                       MDP generation and validation, GROMACS workflows
  reinvent/                 REINVENT 4 setup, configs, scoring, analysis
  cli/                      one entry point per skill
.agents/skills/<skill>/     SKILL.md, reference docs, a thin launcher
tests/                      the test suite and its real data fixtures
tools/                      generators for documentation derived from code
```

## Testing

```bash
pytest                       # the default: hermetic, no network
pytest -v                    # with test names
pytest -m network            # opt in to live API calls
pytest -m requires_vina      # opt in to real docking, if Vina is installed
pytest --cov=agskills        # coverage
```

The default run deselects anything needing the network, an external binary,
or more than a few seconds. Capability-marked tests skip cleanly when the
tool is absent, so an opt-in run is still informative on a machine that
lacks it.

Tests run against real data rather than synthetic stand-ins: two
experimental structures from the PDB, a real AlphaFold model, independently
computed reference properties for ten marketed drugs, a REINVENT output
CSV in the real column format, and authentic AutoDock Vina output. See
`tests/conftest.py` for what each fixture is and why.

## Regenerating derived documentation

Threshold tables are generated from the code so they cannot drift:

```bash
python tools/generate_admet_reference.py
```

A test fails if a `SKILL.md` states a check count the registry does not
perform, so forgetting to regenerate is caught rather than shipped.

## Scientific conventions

Worth knowing before interpreting any output:

- **Molecular weight** is the average mass, as Lipinski's rule and PubChem
  both mean. The monoisotopic mass is reported separately as `exact_mw`.
- **A published rule always reports its authors' thresholds.** A screening
  tier changes which rules must pass; it never redefines a violation. A
  650 Da compound reports a molecular-weight violation at every tier.
- **Gating against advisory.** Published risk *flags*, such as the Pfizer
  3/75 criterion, are reported but do not reject a compound. Gating on them
  would discard most lead-like chemistry.
- **Non-bonded settings belong to the force field.** CHARMM36 requires a
  force-switched Lennard-Jones potential and no dispersion correction;
  AMBER requires a plain cutoff with one. They are not interchangeable.
- **A docking score is triage.** Differences under about 1 kcal/mol are not
  meaningful, and a Vina score is not a binding free energy.
- **pLDDT is not a B-factor.** For a predicted model higher is better on a
  0-100 scale; for a crystal structure the column holds atomic displacement
  parameters, where lower is better and the scale is unbounded. The skills
  detect which they are looking at.

## Citations

The methods these skills wrap, and the thresholds they apply, are from:

- Lipinski CA et al. *Adv Drug Deliv Rev* 1997;23:3-25
- Veber DF et al. *J Med Chem* 2002;45:2615-2623
- Egan WJ et al. *J Med Chem* 2000;43:3867-3877
- Ghose AK et al. *J Comb Chem* 1999;1:55-68
- Muegge I et al. *J Med Chem* 2001;44:1841-1846
- Hughes JD et al. *Bioorg Med Chem Lett* 2008;18:4872-4875
- Baell JB, Holloway GA. *J Med Chem* 2010;53:2719-2740 (PAINS)
- Brenk R et al. *ChemMedChem* 2008;3:435-444
- Ertl P, Schuffenhauer A. *J Cheminform* 2009;1:8 (synthetic accessibility)
- Bickerton GR et al. *Nat Chem* 2012;4:90-98 (QED)
- Eberhardt J et al. *J Chem Inf Model* 2021;61:3891-3898 (AutoDock Vina)
- Huang J, MacKerell AD. *J Comput Chem* 2013;34:2135-2145 (CHARMM36)
- Bernetti M, Bussi G. *J Chem Phys* 2020;153:114107 (C-rescale)
- Lin Z et al. *Science* 2023;379:1123-1130 (ESMFold)
- Varadi M et al. *Nucleic Acids Res* 2022;50:D439-D444 (AlphaFold DB)
- Zdrazil B et al. *Nucleic Acids Res* 2024;52:D1180-D1192 (ChEMBL)
- Sorokina M et al. *J Cheminform* 2021;13:2 (COCONUT)
- Loeffler HH et al. *J Cheminform* 2024;16:20 (REINVENT 4)

Individual thresholds carry their citation in the JSON output, so a result
can be traced to its source.
