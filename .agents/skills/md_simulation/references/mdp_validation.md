# MDP Validation and Protocol Reference

Deep detail for `gmxapi-setup` and `gmxapi-validate`. The SKILL.md carries
what you need to run them; read this when a warning needs defending or a
parameter needs changing.

## Force-field non-bonded settings

Each force field was parameterised with a particular treatment of the
Lennard-Jones tail. These are part of the model, not preferences.

| Force field | GROMACS name | Non-bonded treatment |
|---|---|---|
| `charmm36` | `charmm36-jul2022` | `vdw-modifier = force-switch`, `rvdw_switch` 1.0, `rvdw` 1.2, `rcoulomb` 1.2, **no** `DispCorr` |
| `amber99sb` | `amber99sb` | plain cutoff, `rvdw` 1.0, `rcoulomb` 1.0, `DispCorr = EnerPres` |
| `amber99sb-ildn` | `amber99sb-ildn` | as above |
| `amber14sb` | `amber14sb` | as above |
| `oplsaa` | `oplsaa` | as above; parameterised with TIP4P water |

### Why mixing them is wrong

CHARMM36's force switch already accounts for the truncated dispersion tail.
Adding `DispCorr` on top double-counts it, shifting pressure and density.
The simulation still runs and still looks plausible, which is what makes the
error expensive.

Cross-checking an AMBER-style file against `charmm36` produces four
specific errors:

```
charmm36 requires vdw-modifier = force-switch, but the file has no vdw-modifier.
charmm36 must not use a dispersion correction, but DispCorr = EnerPres is set.
charmm36 expects rvdw = 1.2 nm, but the file has 1.0 nm.
charmm36 expects rcoulomb = 1.2 nm, but the file has 1.0 nm.
```

The same file validates clean against `amber99sb-ildn` — it is valid AMBER,
just mislabelled.

## Every check performed

Each corresponds to a real grompp abort or a silent physics error.

| Check | Why |
|---|---|
| `integrator` is a GROMACS integrator | Typos abort grompp |
| `coulombtype` valid; PME for dynamics | A plain Coulomb cutoff in a charged biomolecular system introduces serious artefacts |
| `vdw-modifier` valid | — |
| `ref_t` count equals `tc-grps` count | The commonest grompp abort |
| `tau_t` count equals `tc-grps` count | Same |
| `ref_t` present when `tcoupl` is on | — |
| `ref_t` within 250–400 K (warning) | Outside the biomolecular range |
| `tau_p`, `ref_p`, `compressibility` present when `pcoupl` is on | grompp aborts without them |
| `refcoord_scaling` set when `-DPOSRES` meets pressure coupling | Otherwise grompp warns and the restraint reference is not scaled with the box |
| `DispCorr` absent with `vdw-modifier = force-switch` | Double-counts the dispersion tail |
| `rvdw_switch` present and below `rvdw` for force-switch | — |
| `rcoulomb` and `rvdw` within 0.2 nm (warning) | A wider gap makes the Verlet buffer inefficient |
| `dt` <= 0.005 ps | Above that, unstable for an all-atom system |
| `dt` > 0.0025 ps requires constraints | Energy is not conserved without them |
| `dt` > 0.002 ps warns | Needs hydrogen mass repartitioning or virtual sites |
| `gen_vel = yes` with `continuation = yes` | They contradict: generating velocities discards the ones being continued from |
| Berendsen thermostat (warning) | Does not sample a correct canonical ensemble |
| `nsteps` non-negative; 0 warns for dynamics | — |

Notes report the simulated time implied by `nsteps` and `dt`.

## The protocol, stage by stage

| Stage | File | Purpose |
|---|---|---|
| Ion placement | `ions.mdp` | Steepest descent so `genion` has a sane configuration. A plain cutoff is adequate here |
| Minimisation | `em.mdp` | Removes clashes from solvation and from the input model |
| NVT | `nvt.mdp` | Brings the system to temperature with the solute restrained, so solvent relaxes around a fixed solute. V-rescale thermostat |
| NPT | `npt.mdp` | Equilibrates density, solute still restrained. **C-rescale** barostat |
| Production | `md.mdp` | Unrestrained. Parrinello-Rahman barostat |

### Why C-rescale for equilibration, Parrinello-Rahman for production

C-rescale (stochastic cell rescaling) is correct far from equilibrium.
Parrinello-Rahman oscillates when started from a freshly solvated box, which
is exactly the NPT stage's starting condition. Production begins from an
equilibrated box, where Parrinello-Rahman is appropriate.

### Temperature coupling groups

Two groups throughout — the solute and everything else. Coupling a small
solute together with a large solvent bath produces the "hot solvent, cold
solute" artefact.

With `--has-ligand` the groups become `Protein_LIG` and `Water_and_ions`, so
the ligand is coupled with the protein rather than the bath. Create them
with `gmx make_ndx` before grompp; the output says so.

## Timestep and constraints

`dt = 0.002` ps is safe with `constraints = h-bonds` and LINCS. Above
0.0025 ps you need hydrogen mass repartitioning or virtual sites. The
conversion is:

```
nsteps = nanoseconds * 1000 / dt_ps
100 ns at 2 fs  ->  50,000,000 steps
```

## Box type

A dodecahedron needs roughly 29% fewer water molecules than a cube for the
same solute-to-boundary clearance, so it is the default. `--box-dist` below
1.0 nm warns: periodic images begin to interact across a 1.2 nm cutoff.

## Ligand topologies

GROMACS cannot generate one. Use CGenFF for CHARMM36 or ACPYPE/GAFF for
AMBER, and include the resulting `.itp` in `topol.top` before grompp. The
skill reports this rather than emitting a workflow that fails at grompp.

## Output size

`--detail full` echoes every parsed MDP parameter back — roughly 30 entries
per file, which is simply the input repeated. The default reports only
errors, warnings and notes, which is the result.
