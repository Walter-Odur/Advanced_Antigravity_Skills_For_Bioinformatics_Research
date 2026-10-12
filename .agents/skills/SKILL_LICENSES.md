# Skill Licenses and Terms of Use

You are responsible for ensuring that your use of these skills complies with
the terms applicable to the data sources and software they reach. Relevant
references are listed below to assist compliance; site addresses change from
time to time, so this list may become incomplete or outdated.

Each skill that touches a third-party source notifies you on first use and
records that notification under `.licenses/<skill-name>_LICENSE.txt` in the
workspace root.

Skill                  | Data source / software | Reference URL(s)
---------------------- | ---------------------- | ----------------
target-preparation     | RCSB PDB               | https://www.rcsb.org/pages/usage-policy
target-preparation     | AlphaFold DB           | https://alphafold.ebi.ac.uk/
target-preparation     | ESMFold / ESM Atlas    | https://esmatlas.com/about
target-preparation     | Meeko                  | https://github.com/forlilab/Meeko
target-preparation     | Open Babel             | https://openbabel.org/
compound-screening     | ChEMBL                 | https://chembl.gitbook.io/chembl-interface-documentation/about
compound-screening     | COCONUT                | https://coconut.naturalproducts.net/about
compound-screening     | ZINC                   | https://zinc.docking.org/
compound-screening     | AutoDock Vina          | https://github.com/ccsb-scripps/AutoDock-Vina
compound-screening     | RDKit                  | https://www.rdkit.org/docs/Overview.html#license
compound-synthesis     | REINVENT 4             | https://github.com/MolecularAI/REINVENT4
compound-synthesis     | REINVENT priors        | https://doi.org/10.5281/zenodo.15641296
md-simulation          | GROMACS                | https://www.gromacs.org/
md-simulation          | CHARMM36 force field   | http://mackerell.umaryland.edu/charmm_ff.shtml
md-simulation          | AMBER force fields     | https://ambermd.org/
drug-discovery-wizard  | all of the above       | see the rows above

## Software licenses

The skills themselves are MIT licensed. The tools and data they reach carry
their own terms, which are not granted by this bundle:

-   **RDKit** — BSD 3-Clause.
-   **AutoDock Vina** — Apache 2.0.
-   **GROMACS** — LGPL 2.1.
-   **REINVENT 4** — Apache 2.0.
-   **Meeko** — LGPL 2.1.
-   **Open Babel** — GPL 2.0.
-   **BioPython** — Biopython License Agreement.
-   **polite-http** — see https://pypi.org/project/polite-http/.

Academic use of a database or method generally carries a citation expectation
rather than a licence restriction. Each skill ships
`references/citation.bib` with the primary references for the methods it
wraps; cite them in any publication that relies on the skill's output.
