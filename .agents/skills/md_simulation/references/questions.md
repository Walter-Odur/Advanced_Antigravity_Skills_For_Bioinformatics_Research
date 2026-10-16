# Questions this skill answers

Find your question, run the command. Nothing here is specific to a protein
family or an organism — the settings follow from the force field and the
system, not the source of the sequence.

**This skill produces files; it does not run simulations.** Everything below
works on a laptop with no GROMACS installed.

## Start here

| You want | Run |
|---|---|
| A complete MD setup from a PDB | `gmxapi-setup` |
| The same, as shell only | `generate-md-setup` |
| To check MDP files you already have | `gmxapi-validate` |
| A cluster submission script | `generate-hpc-script` |

--------------------------------------------------------------------------------

## Getting started

**Q1. I have a cleaned PDB. Give me an MD setup.**

```bash
uv run scripts/md_simulation_api.py gmxapi-setup \
    --pdb work/receptor_clean.pdb --ff charmm36 \
    --temperature 310 --production-ns 100 --output work/md/run_md.py
```

You get five MDP files (`em`, `ions`, `nvt`, `npt`, `md`), a
`setup_md.sh` that runs the eight preparation steps, and a `run_md.py`
gmxapi workflow. Then `bash setup_md.sh`, then submit the production run.

**Q2. Do I need GROMACS installed to use this?**

Only to execute what it writes. Generation and validation need nothing but
Python — the launcher declares no dependencies at all.

**Q3. Which force field should I use?**

| Force field | Notes |
|---|---|
| `charmm36` | CHARMM36m is the common current choice for proteins [2, 3] |
| `amber99sb-ildn` | Long-established; ff14SB/ff19SB are its successors [4, 5] |

Both are well validated for folded proteins. The important thing is not
which you pick but that you do not mix its settings with another's (Q5).
If you are extending published work, use whatever that work used.

**Q4. What temperature and how long?**

310 K for physiological conditions, 300 K for comparison with most
published simulations. Length depends on the question: side-chain
rearrangement in nanoseconds, loop motion in tens to hundreds,
domain motion and binding events far beyond that.

**One long trajectory is weaker evidence than several short ones.** Run
replicas with different initial velocities; conclusions from a single
trajectory are frequently artefacts of it [9, 10].

--------------------------------------------------------------------------------

## The settings, and why they are not yours to mix

**Q5. Why does the skill overwrite my cutoffs?**

Because non-bonded settings belong to the force field. They were part of how
its parameters were fitted, and swapping them changes the physics while the
simulation runs perfectly happily.

| Setting | CHARMM36 | AMBER |
|---|---|---|
| `vdw-modifier` | `force-switch` | `potential-shift-verlet` |
| `rvdw_switch` | `1.0` | — |
| `rvdw` | `1.2` | `1.0` |
| `rcoulomb` | `1.2` | `1.0` |
| `DispCorr` | **absent** | `EnerPres` |

CHARMM36 was parameterised with a force-switched Lennard-Jones potential
between 1.0 and 1.2 nm [2, 6]; adding a dispersion correction on top
double-counts the long-range attraction. AMBER was parameterised with a
plain cutoff **and** the analytic correction, so omitting it is equally
wrong. Mixing them is one of the commonest silent errors in MD setup [1].

**Q6. I have MDP files from a paper / a colleague / an old project. Are they
right?**

```bash
uv run scripts/md_simulation_api.py gmxapi-validate \
    --mdp md.mdp --ff charmm36 --output check.json
```

Passing `--ff` cross-checks the file against that force field's
requirements. Measured, validating an AMBER file against CHARMM36:

```
md.mdp: INVALID
  ERROR: charmm36 requires vdw-modifier = force-switch, but the file has
         potential-shift-verlet. ...
  ERROR: charmm36 must not use a dispersion correction, but
         DispCorr = EnerPres is set. ...
  ERROR: charmm36 expects rvdw = 1.2 nm, but the file has 1.0 nm.
  ERROR: charmm36 expects rcoulomb = 1.2 nm, but the file has 1.0 nm.
```

All four, each with its reason. It also recomputes the simulated time from
`nsteps × dt`, which catches a file whose length is not what its author
thought.

**Q7. What is the difference between an ERROR and a warning?**

An **ERROR** is inconsistent with the declared force field. A **warning** is
unusual but defensible, and the run is still valid:

```
warning: dt = 0.005 ps needs hydrogen mass repartitioning or virtual sites
         to be safe.
warning: No thermostat is set, so this run samples the microcanonical
         ensemble and its temperature will drift.
```

A 5 fs timestep is legitimate with hydrogen mass repartitioning [7, 8]; an
NVE run is a deliberate choice. A validator that rejected every unusual
setting would simply be switched off.

**Q8. Why four stages? Can I skip to production?**

No. Each stage does something the next one depends on:

| Stage | Key settings | Purpose |
|---|---|---|
| `em` | `integrator = steep`, `emtol = 1000` | remove clashes from solvation |
| `nvt` | `-DPOSRES`, `pcoupl = no`, `gen_vel = yes` | equilibrate temperature with the protein restrained |
| `npt` | `-DPOSRES`, `pcoupl = C-rescale`, `refcoord_scaling = com` | equilibrate density and box |
| `md` | no restraints, `pcoupl = Parrinello-Rahman` | production |

Three details in that table are easy to get wrong and expensive to notice:

- **`gen_vel = yes` only in NVT**, with `continuation = no` there and `yes`
  afterwards. Regenerating velocities later discards the equilibration you
  just paid for.
- **`refcoord_scaling = com` whenever position restraints meet pressure
  coupling.** Without it GROMACS scales the box but not the restraint
  reference coordinates, so restraints pull atoms toward positions that
  drift. `grompp` warns; most people pass `-maxwarn` and never read it. The
  validator catches this from the *combination* of two settings, neither
  wrong alone.
- **C-rescale for equilibration, Parrinello–Rahman for production.**
  C-rescale is stochastic and samples the correct ensemble while the box is
  far from equilibrium [11]. Berendsen, the old default, does not sample the
  correct ensemble at all.

**Q9. Why V-rescale and not Berendsen for temperature?**

Velocity rescaling with a stochastic term samples the canonical ensemble
correctly [12]; Berendsen temperature coupling suppresses fluctuations and
can produce the "flying ice cube" artefact. Use `V-rescale` (or
Nosé–Hoover for production if you need strictly canonical dynamics).

**Q10. What is LINCS and why is my timestep 2 fs?**

LINCS constrains bond lengths [13], which removes the fastest vibrations and
lets you take a 2 fs step instead of 0.5 fs. `constraints = h-bonds` is the
usual setting. Going beyond 2 fs needs hydrogen mass repartitioning or
virtual sites [7, 8] — see Q7.

**Q11. Should I simulate with my ligand bound?**

Usually yes, if the question is about the complex. That needs ligand
parameters, which this skill does not generate: use CGenFF for CHARMM, or
GAFF2/`antechamber` for AMBER, and check the resulting parameters rather
than trusting them. Topology generation for a non-standard residue is the
single most error-prone part of setting up a complex.

**Q12. My protein has a metal / a cofactor / a disulfide / a non-standard
residue.**

`pdb2gmx` will stop or guess. Metals and cofactors need parameters;
disulfides are usually detected but confirm them; non-standard residues need
an explicit entry in the force field. Deal with these before running
`setup_md.sh`, not after.

--------------------------------------------------------------------------------

## Running on a cluster

**Q13. Give me a submission script.**

```bash
uv run scripts/md_simulation_api.py generate-hpc-script \
    --job-type gromacs-md --tpr md.tpr --time 24:00:00 --output submit.sh
```

**Q14. My run is longer than the walltime. Will I lose it?**

No, if you use the generated script. It sets `-maxh` just below the
walltime so GROMACS stops cleanly and writes a final checkpoint instead of
being killed mid-write, and resumes with `-cpi` when a checkpoint exists:

| `--time` | `MAXH` |
|---|---|
| `00:30:00` | 0.48 |
| `04:00:00` | 3.92 |
| `24:00:00` | 23.52 |
| `48:00:00` | 47.5 |

Resubmit the same script to continue. A 100 ns run rarely fits in one
allocation, so this is the normal path, not an edge case.

**Q15. Can I use it for jobs other than GROMACS?**

Yes — `--job-type` accepts `gromacs-md`, `reinvent`, `vina-screen`,
`colabfold` and `alphafold2`. A `vina-screen` job is a CPU job array and is
generated **without** a GPU request, because AutoDock Vina cannot use one;
it needs `--array` to partition the ligand set.

**Q16. The script fails on the cluster with `bash\r: bad interpreter`.**

That is a CRLF line ending, and the generated scripts have none —
everything goes through a writer that normalises them. If you see it, the
file was edited or transferred with conversion on. Check with:

```bash
python -c "print(open('submit.sh','rb').read().count(b'\r'))"
```

Expect `0`. Do not rely on `grep -c $'\r'`, which under Git Bash on Windows
can match every line of a clean file.

**Q17. `conda activate` fails inside the job.**

The generated script already handles this with
`eval "$(conda shell.bash hook)"`. A non-interactive SLURM shell has not
sourced the conda init block, so a bare `conda activate` fails there with
"your shell has not been properly configured".

--------------------------------------------------------------------------------

## After the simulation

**Q18. How do I know the simulation is any good?**

This skill does not analyse trajectories, but the checks are standard and
worth doing before you interpret anything:

| Check | Expect |
|---|---|
| Potential energy | plateaued after EM, stable in production |
| Temperature | at your setpoint, no drift |
| Density (NPT) | near 1000 kg/m³ for a solvated protein |
| Pressure | noisy; judge the running average, not instantaneous values |
| Backbone RMSD | plateaus; a monotonic rise means it is still equilibrating |
| RMSF | high only where you expect flexibility — loops, termini |
| Replica agreement | the same conclusion from each [9, 10] |

Discard the equilibration period before analysing, and say in any write-up
how much you discarded and how you chose it [9].

**Q19. How long should I discard?**

Judge it from the trajectory rather than by convention: discard until the
properties you care about stop trending. Quantifying sampling quality and
the uncertainty that follows from it has published best practice [10] — use
it rather than reporting a single number with no error estimate.

**Q20. My RMSD keeps climbing. Is the simulation broken?**

Possibly not — but it is not equilibrated, so do not analyse it as if it
were. Common causes: a genuine conformational change, a structure that was
missing loops that are now flopping, a bad protonation state, or an
insufficient equilibration.

--------------------------------------------------------------------------------

## When something goes wrong

| Symptom | Cause | Fix |
|---|---|---|
| Exit 2 on a walltime | Not `HH:MM:SS` or `D-HH:MM:SS` | `--time 24:00:00` |
| Exit 2 on `vina-screen` | No `--array` | `--array 1-100` |
| `grompp` warns about `refcoord_scaling` | Restraints with pressure coupling | Q8; the validator catches this |
| `pdb2gmx` cannot find a residue | Non-standard residue or ligand | Q11, Q12 |
| LINCS warnings, then a crash | Clashes, bad geometry, or too large a timestep | Re-run EM; check the input structure |
| Box too small warnings | Protein within the cutoff of its own image | Increase the box with `editconf -d` |
| `setup_md.sh` hangs | An interactive prompt | `genion` is fed its group on stdin in the generated script |

--------------------------------------------------------------------------------

## References

Full entries in [citation.bib](citation.bib). Validation rules and their
sources are tabulated in [mdp_validation.md](mdp_validation.md).

1.  Abraham MJ et al. *SoftwareX* 2015;1:19–25 — GROMACS.
2.  Huang J, MacKerell AD. *J Comput Chem* 2013;34:2135–2145 — CHARMM36.
3.  Huang J et al. *Nat Methods* 2017;14:71–73 — CHARMM36m.
4.  Maier JA et al. *J Chem Theory Comput* 2015;11:3696–3713 — ff14SB.
5.  Tian C et al. *J Chem Theory Comput* 2020;16:528–552 — ff19SB.
6.  Steinbach PJ, Brooks BR. *J Comput Chem* 1994;15:667–683 — spherical
    cutoff methods and force switching.
7.  Feenstra KA, Hess B, Berendsen HJC. *J Comput Chem* 1999;20:786–798 —
    hydrogen-rich systems and longer timesteps.
8.  Hopkins CW et al. *J Chem Theory Comput* 2015;11:1864–1874 — hydrogen
    mass repartitioning.
9.  Braun E et al. *Living J Comput Mol Sci* 2018;1:5957 — best practices
    for foundations in molecular simulations.
10. Grossfield A et al. *Living J Comput Mol Sci* 2018;1:5067 — best
    practices for uncertainty and sampling quality.
11. Bernetti M, Bussi G. *J Chem Phys* 2020;153:114107 — C-rescale
    barostat.
12. Bussi G, Donadio D, Parrinello M. *J Chem Phys* 2007;126:014101 —
    canonical sampling through velocity rescaling.
13. Hess B et al. *J Comput Chem* 1997;18:1463–1472 — LINCS.
14. Darden T, York D, Pedersen L. *J Chem Phys* 1993;98:10089–10092 — PME.
15. Essmann U et al. *J Chem Phys* 1995;103:8577–8593 — smooth PME.
16. Parrinello M, Rahman A. *J Appl Phys* 1981;52:7182–7190 — the
    Parrinello–Rahman barostat.
17. Knapp B, Ospina L, Deane CM. *J Chem Theory Comput* 2018;14:6127–6138 —
    the importance of replicas.
18. Lemkul JA. *Living J Comput Mol Sci* 2019;1:5068 — GROMACS tutorials.
