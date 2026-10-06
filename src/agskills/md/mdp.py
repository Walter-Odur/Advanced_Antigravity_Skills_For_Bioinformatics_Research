"""GROMACS MDP generation and validation, with force-field-aware physics.

Two substantive corrections over the code this replaces.

**The requested temperature was discarded.** ``--temperature`` was accepted
by the parser and interpolated into a comment in the generated Python
workflow, but ``generate_mdp()`` took only ``(mdp_type, production_ns)``.
Running with ``--temperature 310`` produced MDP files containing
``ref_t = 300 300`` and ``gen_temp = 300``. grompp reads the MDP, so the
simulation ran at 300 K while the user believed it ran at 310 K.

**The non-bonded settings contradicted the default force field.** The
default force field was ``charmm36``, but every template used
``rvdw = 1.0`` with ``DispCorr = EnerPres``. That is the AMBER convention.
CHARMM36 was parameterised with a force-switched Lennard-Jones potential
and must use ``vdw-modifier = force-switch``, ``rvdw_switch = 1.0``,
``rvdw = 1.2``, ``rcoulomb = 1.2`` and **no** dispersion correction; adding
one double-counts long-range dispersion that the force field already
absorbs into its parameters. Non-bonded settings are now derived from the
chosen force field.

Smaller fixes: ``refcoord_scaling = com`` is set whenever position
restraints are combined with pressure coupling (grompp warns otherwise and
the restraint reference is not scaled with the box), and equilibration uses
the C-rescale barostat, which is stochastically correct and stable far from
equilibrium, rather than Parrinello-Rahman, which oscillates when started
from an unequilibrated box.

References:
    Huang J, MacKerell AD. *J Comput Chem* 2013;34:2135-2145 (CHARMM36).
    Bernetti M, Bussi G. *J Chem Phys* 2020;153:114107 (C-rescale).
    Bussi G, Donadio D, Parrinello M. *J Chem Phys* 2007;126:014101
        (velocity rescaling thermostat).
    GROMACS manual, "Molecular dynamics parameters (.mdp options)".
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from ..errors import UsageError

__all__ = [
    "FORCE_FIELDS",
    "WATER_MODELS",
    "MdpStage",
    "mdp_stages",
    "generate_mdp",
    "steps_for_ns",
    "validate_mdp_text",
    "validate_mdp_file",
]


#: Force-field-specific non-bonded settings. These are not interchangeable:
#: each force field was fitted with a particular treatment of the
#: Lennard-Jones tail, and using another one changes the physics.
FORCE_FIELDS: dict[str, dict[str, Any]] = {
    "charmm36": {
        "gromacs_name": "charmm36-jul2022",
        "recommended_water": "tip3p",
        "nonbonded": {
            "cutoff-scheme": "Verlet",
            "coulombtype": "PME",
            "rcoulomb": 1.2,
            "vdwtype": "cutoff",
            "vdw-modifier": "force-switch",
            "rvdw_switch": 1.0,
            "rvdw": 1.2,
            # No DispCorr: CHARMM36's force switch already accounts for the
            # truncated dispersion tail.
        },
        "note": ("CHARMM36 requires a force-switched LJ potential between "
                 "1.0 and 1.2 nm and no dispersion correction."),
        "citation": "Huang J, MacKerell AD. J Comput Chem 2013;34:2135-2145",
    },
    "amber99sb-ildn": {
        "gromacs_name": "amber99sb-ildn",
        "recommended_water": "tip3p",
        "nonbonded": {
            "cutoff-scheme": "Verlet",
            "coulombtype": "PME",
            "rcoulomb": 1.0,
            "vdwtype": "cutoff",
            "vdw-modifier": "potential-shift-verlet",
            "rvdw": 1.0,
            "DispCorr": "EnerPres",
        },
        "note": ("AMBER force fields use a plain 1.0 nm cutoff with an "
                 "analytic dispersion correction."),
        "citation": "Lindorff-Larsen K et al. Proteins 2010;78:1950-1958",
    },
    "amber99sb": {
        "gromacs_name": "amber99sb",
        "recommended_water": "tip3p",
        "nonbonded": {
            "cutoff-scheme": "Verlet",
            "coulombtype": "PME",
            "rcoulomb": 1.0,
            "vdwtype": "cutoff",
            "vdw-modifier": "potential-shift-verlet",
            "rvdw": 1.0,
            "DispCorr": "EnerPres",
        },
        "note": "AMBER99SB; prefer amber99sb-ildn for side-chain torsions.",
        "citation": "Hornak V et al. Proteins 2006;65:712-725",
    },
    "amber14sb": {
        "gromacs_name": "amber14sb",
        "recommended_water": "tip3p",
        "nonbonded": {
            "cutoff-scheme": "Verlet",
            "coulombtype": "PME",
            "rcoulomb": 1.0,
            "vdwtype": "cutoff",
            "vdw-modifier": "potential-shift-verlet",
            "rvdw": 1.0,
            "DispCorr": "EnerPres",
        },
        "note": "AMBER14SB protein force field.",
        "citation": "Maier JA et al. J Chem Theory Comput 2015;11:3696-3713",
    },
    "oplsaa": {
        "gromacs_name": "oplsaa",
        "recommended_water": "tip4p",
        "nonbonded": {
            "cutoff-scheme": "Verlet",
            "coulombtype": "PME",
            "rcoulomb": 1.0,
            "vdwtype": "cutoff",
            "vdw-modifier": "potential-shift-verlet",
            "rvdw": 1.0,
            "DispCorr": "EnerPres",
        },
        "note": "OPLS-AA/L; TIP4P is the water model it was fitted with.",
        "citation": "Kaminski GA et al. J Phys Chem B 2001;105:6474-6487",
    },
}

#: Water models GROMACS ships, with the force fields they suit.
WATER_MODELS: dict[str, str] = {
    "tip3p": "3-site; standard for CHARMM and AMBER.",
    "tip4p": "4-site; the model OPLS-AA was parameterised with.",
    "tip4pew": "4-site Ewald-corrected variant.",
    "spc": "Simple point charge; common with GROMOS.",
    "spce": "Extended SPC; better bulk properties than SPC.",
    "none": "No solvent (vacuum or implicit).",
}


@dataclass
class MdpStage:
    """One stage of the standard equilibration and production protocol."""

    name: str
    filename: str
    description: str
    params: dict[str, Any] = field(default_factory=dict)
    purpose: str = ""

    def render(self) -> str:
        """Render the stage as MDP text."""
        lines = [
            f"; {self.description}",
            f"; Purpose: {self.purpose}" if self.purpose else "; ",
            "; Generated by agskills. Review before running production work.",
            "",
        ]
        width = max((len(k) for k in self.params), default=20)
        for key, value in self.params.items():
            if value is None:
                continue
            rendered = "yes" if value is True else "no" if value is False else value
            lines.append(f"{key:<{width}} = {rendered}")
        return "\n".join(lines) + "\n"


def steps_for_ns(nanoseconds: float, dt_ps: float) -> int:
    """Convert a simulation length in nanoseconds to a step count.

    >>> steps_for_ns(100, 0.002)
    50000000
    >>> steps_for_ns(1, 0.002)
    500000
    >>> steps_for_ns(0, 0.002)
    0
    """
    if dt_ps <= 0:
        raise UsageError("The timestep must be positive")
    if nanoseconds < 0:
        raise UsageError("The simulation length cannot be negative")
    return int(round(nanoseconds * 1000.0 / dt_ps))


def _resolve_force_field(force_field: str) -> tuple[str, dict[str, Any]]:
    key = (force_field or "").strip().lower()
    aliases = {"charmm": "charmm36", "amber": "amber99sb-ildn",
               "opls": "oplsaa", "amber99sb-ildn": "amber99sb-ildn"}
    key = aliases.get(key, key)
    if key not in FORCE_FIELDS:
        raise UsageError(
            f"Unknown force field {force_field!r}.",
            hint=f"Choose one of: {', '.join(sorted(FORCE_FIELDS))}",
        )
    return key, FORCE_FIELDS[key]


def mdp_stages(*, force_field: str = "charmm36", temperature: float = 300.0,
               pressure: float = 1.0, production_ns: float = 100.0,
               equilibration_ps: float = 100.0, dt_ps: float = 0.002,
               has_ligand: bool = False, seed: int = -1,
               save_interval_ps: float = 10.0) -> list[MdpStage]:
    """Build the full MDP stage set for a protein (or protein-ligand) system.

    Args:
        force_field: A key from :data:`FORCE_FIELDS`. Determines the
            non-bonded treatment.
        temperature: Reference temperature in kelvin. **Applied** to
            ``ref_t`` and ``gen_temp``.
        pressure: Reference pressure in bar.
        production_ns: Production length in nanoseconds. 0 generates
            minimisation and equilibration only.
        equilibration_ps: Length of each of the NVT and NPT stages.
        dt_ps: Integration timestep in picoseconds. 0.002 is safe with
            h-bond constraints; 0.004 needs heavy-hydrogen repartitioning.
        has_ligand: Couple the ligand to the protein's temperature group
            rather than to the solvent group.
        seed: Velocity-generation seed. -1 asks GROMACS for a random seed;
            set a positive integer for a reproducible run.
        save_interval_ps: How often to write coordinates and energies.

    Raises:
        UsageError: On an unknown force field or non-physical inputs.
    """
    ff_key, ff = _resolve_force_field(force_field)
    if temperature <= 0:
        raise UsageError("The temperature must be above 0 K")
    if temperature > 1000:
        raise UsageError(
            f"{temperature} K is not a biomolecular simulation temperature.",
            hint="Body temperature is 310 K; most protocols use 300 K.",
        )
    if pressure <= 0:
        raise UsageError("The reference pressure must be positive")
    if dt_ps <= 0 or dt_ps > 0.005:
        raise UsageError(
            f"A timestep of {dt_ps} ps is outside the usable range.",
            hint="Use 0.002 ps with h-bond constraints. Above 0.0025 ps you "
                 "need virtual sites or hydrogen mass repartitioning.",
        )

    nonbonded = dict(ff["nonbonded"])
    # Two temperature-coupling groups: the solute and everything else.
    # Coupling a small solute separately from a large solvent bath is what
    # avoids the "hot solvent, cold solute" artefact.
    if has_ligand:
        tc_groups = "Protein_LIG Water_and_ions"
        group_note = ("The ligand is coupled with the protein. Create the "
                      "Protein_LIG and Water_and_ions groups with "
                      "'gmx make_ndx' before grompp.")
    else:
        tc_groups = "Protein Non-Protein"
        group_note = "Default groups; no index file needed."
    n_groups = 2
    ref_t = " ".join([f"{temperature:g}"] * n_groups)
    tau_t = " ".join(["0.1"] * n_groups)

    output_every = max(1, steps_for_ns(save_interval_ps / 1000.0, dt_ps))
    equil_steps = max(1, steps_for_ns(equilibration_ps / 1000.0, dt_ps))

    common_md = {
        "integrator": "md",
        "dt": dt_ps,
        "nstxout-compressed": output_every,
        "nstenergy": output_every,
        "nstlog": output_every,
        "constraint_algorithm": "lincs",
        "constraints": "h-bonds",
        "lincs_iter": 1,
        "lincs_order": 4,
        "nstlist": 20,
        "pbc": "xyz",
        **nonbonded,
        "pme_order": 4,
        "fourierspacing": 0.12,
        "tcoupl": "V-rescale",
        "tc-grps": tc_groups,
        "tau_t": tau_t,
        "ref_t": ref_t,
    }

    stages: list[MdpStage] = [
        MdpStage(
            name="ions", filename="ions.mdp",
            description="Steepest-descent minimisation for ion placement",
            purpose="Prepares the topology so 'gmx genion' can neutralise "
                    "the system.",
            params={
                "integrator": "steep", "emtol": 1000.0, "emstep": 0.01,
                "nsteps": 50_000, "nstlist": 10,
                "cutoff-scheme": "Verlet",
                # A plain cutoff is adequate here: this run only needs a
                # sane starting configuration for genion.
                "coulombtype": "cutoff", "rcoulomb": 1.0, "rvdw": 1.0,
                "pbc": "xyz",
            },
        ),
        MdpStage(
            name="em", filename="em.mdp",
            description="Steepest-descent energy minimisation",
            purpose="Removes steric clashes from solvation and from the "
                    "input model before any dynamics.",
            params={
                "integrator": "steep", "emtol": 1000.0, "emstep": 0.01,
                "nsteps": 50_000, "nstlist": 10, "pbc": "xyz",
                **nonbonded,
            },
        ),
        MdpStage(
            name="nvt", filename="nvt.mdp",
            description=f"NVT equilibration, {equilibration_ps:g} ps at "
                        f"{temperature:g} K",
            purpose="Brings the system to temperature with the solute "
                    "restrained, so solvent relaxes around a fixed solute.",
            params={
                "define": "-DPOSRES",
                **common_md,
                "nsteps": equil_steps,
                "continuation": "no",
                "pcoupl": "no",
                "gen_vel": "yes",
                "gen_temp": f"{temperature:g}",
                "gen_seed": seed,
            },
        ),
        MdpStage(
            name="npt", filename="npt.mdp",
            description=f"NPT equilibration, {equilibration_ps:g} ps at "
                        f"{temperature:g} K and {pressure:g} bar",
            purpose="Equilibrates density at constant pressure, still with "
                    "the solute restrained.",
            params={
                "define": "-DPOSRES",
                **common_md,
                "nsteps": equil_steps,
                "continuation": "yes",
                # C-rescale is stochastic and correct far from equilibrium;
                # Parrinello-Rahman can ring badly when started from an
                # unequilibrated box, which is exactly this stage.
                "pcoupl": "C-rescale",
                "pcoupltype": "isotropic",
                "tau_p": 1.0,
                "ref_p": f"{pressure:g}",
                "compressibility": "4.5e-5",
                # Required whenever position restraints coexist with
                # pressure coupling, or grompp warns and the restraint
                # reference is not scaled with the box.
                "refcoord_scaling": "com",
                "gen_vel": "no",
            },
        ),
    ]

    if production_ns > 0:
        stages.append(MdpStage(
            name="production", filename="md.mdp",
            description=f"Production MD, {production_ns:g} ns at "
                        f"{temperature:g} K and {pressure:g} bar",
            purpose="The unrestrained trajectory used for analysis.",
            params={
                **common_md,
                "nsteps": steps_for_ns(production_ns, dt_ps),
                "continuation": "yes",
                "pcoupl": "Parrinello-Rahman",
                "pcoupltype": "isotropic",
                "tau_p": 2.0,
                "ref_p": f"{pressure:g}",
                "compressibility": "4.5e-5",
                "gen_vel": "no",
            },
        ))

    # Record the provenance so the physics choices are auditable from the
    # files themselves.
    for stage in stages:
        stage.description += (
            f" | force field {ff_key} ({ff['gromacs_name']}) | {group_note}"
        )
    return stages


def generate_mdp(stage: str, **kwargs) -> str:
    """Render a single named MDP stage.

    Raises:
        UsageError: If *stage* is not a known stage name.
    """
    wanted = (stage or "").strip().lower()
    aliases = {"md": "production", "prod": "production"}
    wanted = aliases.get(wanted, wanted)
    stages = {s.name: s for s in mdp_stages(**kwargs)}
    if wanted not in stages:
        raise UsageError(
            f"Unknown MDP stage {stage!r}.",
            hint=f"Choose one of: {', '.join(stages)}",
        )
    return stages[wanted].render()


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

VALID_INTEGRATORS = frozenset({
    "md", "md-vv", "md-vv-avek", "sd", "bd", "steep", "cg", "l-bfgs",
    "nm", "tpi", "tpic", "mimic",
})
VALID_COULOMBTYPES = frozenset({
    "cut-off", "cutoff", "ewald", "pme", "p3m-ad", "reaction-field",
    "user", "pme-switch", "pme-user", "pme-user-switchnb",
})
VALID_TCOUPL = frozenset({
    "no", "berendsen", "nose-hoover", "andersen", "andersen-massive",
    "v-rescale",
})
VALID_PCOUPL = frozenset({
    "no", "berendsen", "parrinello-rahman", "mttk", "c-rescale",
})
VALID_CONSTRAINTS = frozenset({
    "none", "h-bonds", "all-bonds", "h-angles", "all-angles",
})
VALID_VDW_MODIFIERS = frozenset({
    "potential-shift-verlet", "potential-shift", "none", "force-switch",
    "potential-switch",
})

#: Thermostats and barostats that do not sample the correct ensemble.
DEPRECATED_COUPLING = {
    "berendsen": ("The Berendsen thermostat does not produce a correct "
                  "canonical ensemble; use v-rescale. (GROMACS deprecates "
                  "it for tcoupl.)"),
}


def parse_mdp(text: str) -> dict[str, str]:
    """Parse MDP text into a parameter dict, normalising key spelling.

    GROMACS treats ``-`` and ``_`` in option names as equivalent and is
    case-insensitive, so keys are normalised to lower case with dashes.
    """
    params: dict[str, str] = {}
    for raw in text.splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, _, value = line.partition("=")
        normalised = key.strip().lower().replace("_", "-")
        params[normalised] = value.strip()
    return params


def _check_force_field_consistency(params: dict[str, str], force_field: str,
                                   errors: list[str], warnings: list[str]) -> None:
    """Cross-check non-bonded settings against a declared force field.

    An MDP file does not record which force field it is meant for, so this
    inconsistency is invisible to a file-only validator. It is also the most
    consequential error in the code this package replaces: the default force
    field was CHARMM36 while every template carried AMBER's ``rvdw = 1.0``
    plus ``DispCorr = EnerPres``.
    """
    key, spec = _resolve_force_field(force_field)
    expected = spec["nonbonded"]
    integrator = params.get("integrator", "").strip().lower()
    if integrator in ("steep", "cg", "l-bfgs"):
        return  # minimisation is insensitive to these details

    def get(option: str) -> str | None:
        return params.get(option.lower().replace("_", "-"))

    wants_switch = expected.get("vdw-modifier") == "force-switch"
    has_dispcorr = (get("DispCorr") or "no").strip().lower() not in ("no", "")
    modifier = (get("vdw-modifier") or "").strip().lower()

    if wants_switch:
        if modifier != "force-switch":
            errors.append(
                f"{key} requires vdw-modifier = force-switch, but the file "
                f"has {modifier or 'no vdw-modifier'}. {spec['note']}"
            )
        if has_dispcorr:
            errors.append(
                f"{key} must not use a dispersion correction, but DispCorr "
                f"= {get('DispCorr')} is set. {spec['note']}"
            )
    elif expected.get("DispCorr") and not has_dispcorr:
        warnings.append(
            f"{key} is normally run with DispCorr = EnerPres; the file sets "
            "none, which shifts pressure and density."
        )

    for option in ("rvdw", "rcoulomb"):
        want = expected.get(option)
        actual = get(option)
        if want is None or actual is None:
            continue
        try:
            if abs(float(actual) - float(want)) > 1e-6:
                errors.append(
                    f"{key} expects {option} = {want:g} nm, but the file has "
                    f"{actual} nm."
                )
        except ValueError:
            errors.append(f"{option} is not a number: {actual!r}")


def validate_mdp_text(text: str, *, name: str = "<mdp>",
                      force_field: str | None = None) -> dict[str, Any]:
    """Validate MDP content for internal consistency and physical sanity.

    Args:
        text: MDP file content.
        name: Label used in the report.
        force_field: If given, additionally cross-check the non-bonded
            settings against that force field's requirements. An MDP file
            does not state its intended force field, so this mismatch
            cannot be detected without being told.

    Checks that the previous validator did not perform, each of which
    corresponds to a real grompp failure or a silent physics error:

    * ``ref_t`` and ``tau_t`` must have one value per ``tc-grps`` entry.
      A mismatch is one of the most common grompp aborts.
    * Pressure coupling requires ``tau_p``, ``ref_p`` and
      ``compressibility``.
    * Position restraints plus pressure coupling require
      ``refcoord_scaling``.
    * ``DispCorr`` must not be combined with a force-switched van der Waals
      potential, which is the CHARMM36 setting.
    * A timestep above 0.002 ps requires constraints.
    """
    params = parse_mdp(text)
    errors: list[str] = []
    warnings: list[str] = []
    notes: list[str] = []

    def lower(key: str) -> str:
        return params.get(key, "").strip().lower()

    integrator = lower("integrator")
    if not integrator:
        warnings.append("integrator is not set; GROMACS defaults to 'md'.")
    elif integrator not in VALID_INTEGRATORS:
        errors.append(f"integrator '{params['integrator']}' is not a GROMACS "
                      f"integrator.")

    is_minimisation = integrator in ("steep", "cg", "l-bfgs")

    coulombtype = lower("coulombtype")
    if coulombtype and coulombtype not in VALID_COULOMBTYPES:
        errors.append(f"coulombtype '{params['coulombtype']}' is not valid.")
    if not is_minimisation and coulombtype in ("cutoff", "cut-off"):
        warnings.append(
            "A plain Coulomb cutoff in a dynamics run introduces serious "
            "artefacts in a charged biomolecular system; use PME."
        )

    vdw_modifier = lower("vdw-modifier")
    if vdw_modifier and vdw_modifier not in VALID_VDW_MODIFIERS:
        errors.append(f"vdw-modifier '{params['vdw-modifier']}' is not valid.")

    # --- numeric fields --------------------------------------------------
    def number(key: str) -> float | None:
        if key not in params:
            return None
        try:
            return float(params[key])
        except ValueError:
            errors.append(f"{key} is not a number: {params[key]!r}")
            return None

    nsteps = number("nsteps")
    if nsteps is not None:
        if nsteps < 0:
            errors.append(f"nsteps cannot be negative (got {nsteps:g}).")
        elif nsteps == 0 and not is_minimisation:
            warnings.append("nsteps is 0, so this run performs no dynamics.")

    dt = number("dt")
    if dt is not None and not is_minimisation:
        if dt <= 0:
            errors.append(f"dt must be positive (got {dt:g}).")
        elif dt > 0.005:
            errors.append(
                f"dt = {dt:g} ps is unstable for an all-atom biomolecular "
                "system."
            )
        elif dt > 0.0025 and lower("constraints") in ("", "none"):
            errors.append(
                f"dt = {dt:g} ps without bond constraints will not conserve "
                "energy. Set constraints = h-bonds."
            )
        elif dt > 0.002:
            warnings.append(
                f"dt = {dt:g} ps needs hydrogen mass repartitioning or "
                "virtual sites to be safe."
            )
        if nsteps is not None and dt > 0:
            notes.append(
                f"Simulated time: {nsteps * dt / 1000.0:g} ns "
                f"({nsteps:g} steps x {dt:g} ps)."
            )

    # --- temperature coupling -------------------------------------------
    tcoupl = lower("tcoupl")
    if tcoupl and tcoupl not in VALID_TCOUPL:
        errors.append(f"tcoupl '{params['tcoupl']}' is not valid.")
    if tcoupl in DEPRECATED_COUPLING:
        warnings.append(DEPRECATED_COUPLING[tcoupl])

    if tcoupl and tcoupl != "no":
        groups = params.get("tc-grps", "").split()
        ref_t = params.get("ref-t", "").split()
        tau_t = params.get("tau-t", "").split()
        if not groups:
            errors.append("tcoupl is on but tc-grps is not set.")
        if groups and ref_t and len(groups) != len(ref_t):
            errors.append(
                f"tc-grps has {len(groups)} group(s) but ref_t has "
                f"{len(ref_t)} value(s); grompp requires one per group."
            )
        if groups and tau_t and len(groups) != len(tau_t):
            errors.append(
                f"tc-grps has {len(groups)} group(s) but tau_t has "
                f"{len(tau_t)} value(s); grompp requires one per group."
            )
        if not ref_t:
            errors.append("tcoupl is on but ref_t is not set.")
        for value in ref_t:
            try:
                kelvin = float(value)
            except ValueError:
                errors.append(f"ref_t value {value!r} is not a number.")
                continue
            if kelvin <= 0:
                errors.append(f"ref_t = {value} K is not physical.")
            elif not 250 <= kelvin <= 400:
                warnings.append(
                    f"ref_t = {value} K is outside the usual biomolecular "
                    "range of 250-400 K."
                )
    elif not is_minimisation and not tcoupl:
        warnings.append(
            "No thermostat is set, so this run samples the microcanonical "
            "ensemble and its temperature will drift."
        )

    # --- pressure coupling ----------------------------------------------
    pcoupl = lower("pcoupl")
    if pcoupl and pcoupl not in VALID_PCOUPL:
        errors.append(f"pcoupl '{params['pcoupl']}' is not valid.")
    if pcoupl and pcoupl != "no":
        for required in ("tau-p", "ref-p", "compressibility"):
            if required not in params:
                errors.append(
                    f"pcoupl = {params['pcoupl']} requires "
                    f"{required.replace('-', '_')} to be set."
                )
        if "-dposres" in params.get("define", "").lower():
            if "refcoord-scaling" not in params:
                errors.append(
                    "Position restraints (-DPOSRES) are combined with "
                    "pressure coupling but refcoord_scaling is not set. "
                    "grompp warns, and restraint reference coordinates are "
                    "not scaled with the box. Set refcoord_scaling = com."
                )

    # --- constraints -----------------------------------------------------
    constraints = lower("constraints")
    if constraints and constraints not in VALID_CONSTRAINTS:
        errors.append(f"constraints '{params['constraints']}' is not valid.")

    # --- force-field consistency ----------------------------------------
    dispcorr = lower("dispcorr")
    if dispcorr and dispcorr != "no" and vdw_modifier == "force-switch":
        errors.append(
            "DispCorr is set together with vdw-modifier = force-switch. "
            "A force-switched potential (the CHARMM36 convention) already "
            "accounts for the dispersion tail, so this double-counts it."
        )
    rvdw = number("rvdw")
    rvdw_switch = number("rvdw-switch")
    if vdw_modifier == "force-switch":
        if rvdw_switch is None:
            errors.append("vdw-modifier = force-switch requires rvdw_switch.")
        elif rvdw is not None and rvdw_switch >= rvdw:
            errors.append(
                f"rvdw_switch ({rvdw_switch:g}) must be less than rvdw "
                f"({rvdw:g})."
            )
    rcoulomb = number("rcoulomb")
    if rcoulomb is not None and rvdw is not None and coulombtype == "pme":
        if abs(rcoulomb - rvdw) > 0.2001:
            warnings.append(
                f"rcoulomb ({rcoulomb:g}) and rvdw ({rvdw:g}) differ by more "
                "than 0.2 nm, which makes the Verlet buffer inefficient."
            )

    # --- velocity generation --------------------------------------------
    if lower("gen-vel") == "yes" and lower("continuation") == "yes":
        errors.append(
            "gen_vel = yes with continuation = yes discards the velocities "
            "being continued from. Use one or the other."
        )

    if force_field:
        _check_force_field_consistency(params, force_field, errors, warnings)

    return {
        "file": name,
        "valid": not errors,
        "force_field_checked": force_field,
        "n_parameters": len(params),
        "errors": errors,
        "warnings": warnings,
        "notes": notes,
        "parameters": params,
    }


def validate_mdp_file(path: str, *, force_field: str | None = None
                      ) -> dict[str, Any]:
    """Validate an MDP file on disk, optionally against a force field."""
    from ..io_utils import require_file
    file_path = require_file(path)
    return validate_mdp_text(
        file_path.read_text(encoding="utf-8", errors="replace"),
        name=str(file_path),
        force_field=force_field,
    )
