"""GROMACS MDP generation and validation.

Two substantive defects are regression-tested here.

**The requested temperature was discarded.** ``--temperature`` was parsed
and interpolated into a comment in the generated Python, but
``generate_mdp()`` took only ``(mdp_type, production_ns)``. A run requested
at 310 K produced ``ref_t = 300 300``, and grompp reads the MDP.

**The non-bonded settings contradicted the force field.** The default force
field was CHARMM36 while every template carried AMBER's ``rvdw = 1.0`` plus
``DispCorr = EnerPres``. CHARMM36 was parameterised with a force-switched
Lennard-Jones potential between 1.0 and 1.2 nm and no dispersion
correction; adding one double-counts a tail the parameters already absorb.
"""

from __future__ import annotations

import pytest

from agskills.errors import UsageError
from agskills.md.mdp import (
    FORCE_FIELDS,
    generate_mdp,
    mdp_stages,
    parse_mdp,
    steps_for_ns,
    validate_mdp_text,
)


def params(stage: str, **kwargs) -> dict[str, str]:
    return parse_mdp(generate_mdp(stage, **kwargs))


# ---------------------------------------------------------------------------
# Temperature propagation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("temperature", [280.0, 300.0, 310.0, 323.15])
def test_requested_temperature_reaches_every_coupled_stage(temperature):
    """ref_t must equal the requested temperature, for every group.

    This is the regression test for the discarded-temperature defect.
    """
    for stage in ("nvt", "npt", "production"):
        values = params(stage, temperature=temperature)
        groups = values["tc-grps"].split()
        reference = values["ref-t"].split()
        assert len(reference) == len(groups)
        for value in reference:
            assert float(value) == pytest.approx(temperature), (
                f"{stage}: ref_t {value} but {temperature} was requested"
            )


def test_velocity_generation_temperature_also_follows(tmp_path):
    """gen_temp seeds the initial velocities and must match too."""
    values = params("nvt", temperature=310.0)
    assert float(values["gen-temp"]) == pytest.approx(310.0)
    assert values["gen-vel"] == "yes"


def test_pressure_propagates_to_the_barostat():
    values = params("npt", pressure=1.5)
    assert float(values["ref-p"]) == pytest.approx(1.5)


def test_non_physical_temperature_is_refused():
    with pytest.raises(UsageError):
        mdp_stages(temperature=0)
    with pytest.raises(UsageError):
        mdp_stages(temperature=-10)
    with pytest.raises(UsageError) as excinfo:
        mdp_stages(temperature=5000)
    assert "310 K" in str(excinfo.value)


# ---------------------------------------------------------------------------
# Force-field-specific non-bonded settings
# ---------------------------------------------------------------------------


def test_charmm36_uses_force_switch_and_no_dispersion_correction():
    """CHARMM36's published non-bonded protocol, not AMBER's.

    This is the regression test for the force-field mismatch.
    """
    for stage in ("em", "nvt", "npt", "production"):
        values = params(stage, force_field="charmm36")
        assert values["vdw-modifier"] == "force-switch", stage
        assert float(values["rvdw-switch"]) == pytest.approx(1.0), stage
        assert float(values["rvdw"]) == pytest.approx(1.2), stage
        assert float(values["rcoulomb"]) == pytest.approx(1.2), stage
        assert "dispcorr" not in values, (
            f"{stage}: CHARMM36 must not set DispCorr"
        )


def test_amber_uses_a_plain_cutoff_with_a_dispersion_correction():
    for family in ("amber99sb-ildn", "amber14sb"):
        values = params("production", force_field=family)
        assert float(values["rvdw"]) == pytest.approx(1.0), family
        assert float(values["rcoulomb"]) == pytest.approx(1.0), family
        assert values["dispcorr"] == "EnerPres", family
        assert values["vdw-modifier"] == "potential-shift-verlet", family


def test_the_two_families_are_genuinely_different():
    """If both produced the same file, the distinction would be cosmetic."""
    charmm = params("production", force_field="charmm36")
    amber = params("production", force_field="amber99sb-ildn")
    assert charmm["rvdw"] != amber["rvdw"]
    assert ("dispcorr" in amber) and ("dispcorr" not in charmm)


def test_every_declared_force_field_generates_a_valid_stage_set():
    for name in FORCE_FIELDS:
        stages = mdp_stages(force_field=name, production_ns=1.0)
        assert len(stages) == 5
        for stage in stages:
            report = validate_mdp_text(stage.render(), name=stage.filename,
                                       force_field=name)
            assert report["valid"], (name, stage.filename, report["errors"])


def test_unknown_force_field_lists_the_options():
    with pytest.raises(UsageError) as excinfo:
        mdp_stages(force_field="made_up")
    message = str(excinfo.value)
    assert "charmm36" in message and "oplsaa" in message


def test_force_field_aliases_resolve():
    assert params("em", force_field="charmm")["rvdw"] == \
        params("em", force_field="charmm36")["rvdw"]


# ---------------------------------------------------------------------------
# Step counts
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("ns,dt,expected", [
    (100.0, 0.002, 50_000_000),
    (1.0, 0.002, 500_000),
    (0.1, 0.002, 50_000),
    (100.0, 0.001, 100_000_000),
    (10.0, 0.004, 2_500_000),
    (0.0, 0.002, 0),
])
def test_step_count_arithmetic(ns, dt, expected):
    """nsteps = nanoseconds x 1000 / timestep, with no rounding drift."""
    assert steps_for_ns(ns, dt) == expected


def test_production_nsteps_matches_the_requested_length():
    values = params("production", production_ns=100.0, dt_ps=0.002)
    assert int(values["nsteps"]) == 50_000_000
    assert float(values["dt"]) == pytest.approx(0.002)
    # And the round trip is self-consistent.
    simulated_ns = int(values["nsteps"]) * float(values["dt"]) / 1000.0
    assert simulated_ns == pytest.approx(100.0)


def test_equilibration_length_is_honoured():
    values = params("nvt", equilibration_ps=250.0, dt_ps=0.002)
    assert int(values["nsteps"]) == 125_000


def test_zero_production_omits_the_production_stage():
    """A setup-only request should not generate a production run input."""
    names = [s.name for s in mdp_stages(production_ns=0)]
    assert "production" not in names
    assert names == ["ions", "em", "nvt", "npt"]


def test_invalid_timestep_is_refused():
    with pytest.raises(UsageError):
        mdp_stages(dt_ps=0)
    with pytest.raises(UsageError) as excinfo:
        mdp_stages(dt_ps=0.01)
    assert "virtual sites" in str(excinfo.value) or "repartition" in \
        str(excinfo.value)


def test_negative_length_is_refused():
    with pytest.raises(UsageError):
        steps_for_ns(-1, 0.002)


# ---------------------------------------------------------------------------
# Position restraints and pressure coupling
# ---------------------------------------------------------------------------


def test_restrained_stages_set_refcoord_scaling():
    """POSRES with pressure coupling requires refcoord_scaling.

    Without it grompp warns and the restraint reference coordinates are not
    scaled with the box. The original npt.mdp omitted it.
    """
    npt = params("npt")
    assert "-DPOSRES" in npt["define"]
    assert npt["pcoupl"] != "no"
    assert npt["refcoord-scaling"] == "com"


def test_nvt_restrains_without_pressure_coupling():
    nvt = params("nvt")
    assert "-DPOSRES" in nvt["define"]
    assert nvt["pcoupl"] == "no"


def test_production_is_unrestrained():
    production = params("production")
    assert "define" not in production or "-DPOSRES" not in \
        production.get("define", "")
    assert "refcoord-scaling" not in production


def test_equilibration_uses_c_rescale_and_production_parrinello_rahman():
    """C-rescale is stable far from equilibrium; Parrinello-Rahman is not.

    Starting Parrinello-Rahman from a freshly solvated box makes the volume
    ring, which is why the equilibration stage uses the stochastic
    barostat and only production switches over.
    """
    assert params("npt")["pcoupl"] == "C-rescale"
    assert params("production")["pcoupl"] == "Parrinello-Rahman"


def test_continuation_flags_are_consistent():
    """gen_vel and continuation must not contradict each other."""
    nvt = params("nvt")
    assert nvt["gen-vel"] == "yes" and nvt["continuation"] == "no"
    for stage in ("npt", "production"):
        values = params(stage)
        assert values["gen-vel"] == "no"
        assert values["continuation"] == "yes"


def test_minimisation_stages_use_a_minimiser():
    for stage in ("ions", "em"):
        assert params(stage)["integrator"] == "steep"
        assert "dt" not in params(stage)


def test_ligand_systems_use_a_combined_temperature_group():
    """A ligand belongs with the solute, not the solvent bath."""
    with_ligand = params("production", has_ligand=True)
    assert "LIG" in with_ligand["tc-grps"]
    without = params("production", has_ligand=False)
    assert without["tc-grps"] == "Protein Non-Protein"
    # Group and value counts still agree.
    assert len(with_ligand["tc-grps"].split()) == \
        len(with_ligand["ref-t"].split())


def test_seed_is_recorded_for_reproducibility():
    assert params("nvt", seed=12345)["gen-seed"] == "12345"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_validator_accepts_our_own_output():
    for stage in mdp_stages(production_ns=10.0):
        report = validate_mdp_text(stage.render(), name=stage.filename)
        assert report["valid"], (stage.filename, report["errors"])
        assert report["n_parameters"] > 5


def test_validator_catches_a_group_count_mismatch():
    """One ref_t per tc-grps entry: the most common grompp abort."""
    bad = (
        "integrator = md\nnsteps = 1000\ndt = 0.002\n"
        "tcoupl = v-rescale\ntc-grps = Protein Non-Protein\n"
        "tau_t = 0.1 0.1\nref_t = 300\n"
    )
    report = validate_mdp_text(bad)
    assert not report["valid"]
    assert any("ref_t" in e and "group" in e for e in report["errors"])


def test_validator_catches_a_tau_t_mismatch():
    bad = (
        "integrator = md\nnsteps = 1000\ndt = 0.002\n"
        "tcoupl = v-rescale\ntc-grps = Protein Non-Protein\n"
        "tau_t = 0.1\nref_t = 300 300\n"
    )
    report = validate_mdp_text(bad)
    assert not report["valid"]
    assert any("tau_t" in e for e in report["errors"])


def test_validator_requires_barostat_parameters():
    bad = (
        "integrator = md\nnsteps = 1000\ndt = 0.002\n"
        "tcoupl = v-rescale\ntc-grps = System\ntau_t = 0.1\nref_t = 300\n"
        "pcoupl = Parrinello-Rahman\n"
    )
    report = validate_mdp_text(bad)
    assert not report["valid"]
    for required in ("tau_p", "ref_p", "compressibility"):
        assert any(required in e for e in report["errors"]), required


def test_validator_catches_missing_refcoord_scaling():
    bad = (
        "define = -DPOSRES\nintegrator = md\nnsteps = 1000\ndt = 0.002\n"
        "tcoupl = v-rescale\ntc-grps = System\ntau_t = 0.1\nref_t = 300\n"
        "pcoupl = C-rescale\ntau_p = 1.0\nref_p = 1.0\n"
        "compressibility = 4.5e-5\n"
    )
    report = validate_mdp_text(bad)
    assert not report["valid"]
    assert any("refcoord_scaling" in e for e in report["errors"])


def test_validator_catches_an_unstable_timestep():
    bad = "integrator = md\nnsteps = 1000\ndt = 0.01\nconstraints = h-bonds\n"
    report = validate_mdp_text(bad)
    assert not report["valid"]
    assert any("unstable" in e for e in report["errors"])


def test_validator_catches_a_long_timestep_without_constraints():
    bad = "integrator = md\nnsteps = 1000\ndt = 0.004\nconstraints = none\n"
    report = validate_mdp_text(bad)
    assert not report["valid"]
    assert any("constraints" in e for e in report["errors"])


def test_validator_catches_contradictory_velocity_settings():
    bad = (
        "integrator = md\nnsteps = 1000\ndt = 0.002\n"
        "gen_vel = yes\ncontinuation = yes\n"
    )
    report = validate_mdp_text(bad)
    assert not report["valid"]
    assert any("gen_vel" in e for e in report["errors"])


def test_validator_catches_dispcorr_with_force_switch():
    """The CHARMM/AMBER contradiction, detectable from the file alone."""
    bad = (
        "integrator = md\nnsteps = 1000\ndt = 0.002\n"
        "vdw-modifier = force-switch\nrvdw_switch = 1.0\nrvdw = 1.2\n"
        "DispCorr = EnerPres\n"
    )
    report = validate_mdp_text(bad)
    assert not report["valid"]
    assert any("double-counts" in e for e in report["errors"])


def test_validator_cross_checks_against_a_declared_force_field():
    """AMBER settings under a CHARMM36 label must be caught.

    An MDP file does not record its intended force field, so this
    inconsistency is invisible without being told which one was meant. It
    is exactly the defect the previous generator produced.
    """
    amber_style = (
        "integrator = md\nnsteps = 1000\ndt = 0.002\n"
        "cutoff-scheme = Verlet\ncoulombtype = PME\n"
        "rcoulomb = 1.0\nrvdw = 1.0\nDispCorr = EnerPres\n"
        "tcoupl = v-rescale\ntc-grps = System\ntau_t = 0.1\nref_t = 300\n"
    )
    # Valid on its own terms.
    assert validate_mdp_text(amber_style)["valid"]
    # Valid as AMBER.
    assert validate_mdp_text(amber_style, force_field="amber99sb-ildn")["valid"]
    # Invalid as CHARMM36, with four specific complaints.
    report = validate_mdp_text(amber_style, force_field="charmm36")
    assert not report["valid"]
    assert report["force_field_checked"] == "charmm36"
    errors = " ".join(report["errors"])
    assert "force-switch" in errors
    assert "DispCorr" in errors
    assert "rvdw" in errors
    assert "rcoulomb" in errors


def test_validator_warns_about_a_deprecated_thermostat():
    text = (
        "integrator = md\nnsteps = 1000\ndt = 0.002\n"
        "tcoupl = berendsen\ntc-grps = System\ntau_t = 0.1\nref_t = 300\n"
    )
    report = validate_mdp_text(text)
    assert any("canonical" in w for w in report["warnings"])


def test_validator_reports_simulated_time():
    text = "integrator = md\nnsteps = 500000\ndt = 0.002\n"
    report = validate_mdp_text(text)
    assert any("1 ns" in note for note in report["notes"])


def test_validator_rejects_unknown_enum_values():
    report = validate_mdp_text("integrator = bogus\n")
    assert not report["valid"]
    report = validate_mdp_text("integrator = md\ncoulombtype = magic\n")
    assert not report["valid"]


def test_validator_warns_about_a_plain_cutoff_in_dynamics():
    text = (
        "integrator = md\nnsteps = 1000\ndt = 0.002\ncoulombtype = cutoff\n"
        "tcoupl = v-rescale\ntc-grps = System\ntau_t = 0.1\nref_t = 300\n"
    )
    report = validate_mdp_text(text)
    assert any("PME" in w for w in report["warnings"])


def test_validator_tolerates_cutoff_in_minimisation():
    """The ions stage legitimately uses a plain cutoff."""
    report = validate_mdp_text(generate_mdp("ions"), force_field="charmm36")
    assert report["valid"]
    assert not any("PME" in w for w in report["warnings"])


def test_parse_mdp_normalises_key_spelling():
    """GROMACS treats dashes and underscores as equivalent, case-insensitively."""
    values = parse_mdp("REF_T = 300\nvdw_modifier = force-switch\n")
    assert values["ref-t"] == "300"
    assert values["vdw-modifier"] == "force-switch"


def test_parse_mdp_strips_comments():
    values = parse_mdp("nsteps = 1000 ; this many steps\n; a whole comment\n")
    assert values["nsteps"] == "1000"
    assert len(values) == 1


def test_unknown_stage_name_lists_the_options():
    with pytest.raises(UsageError) as excinfo:
        generate_mdp("equilibrate")
    assert "nvt" in str(excinfo.value)


def test_md_is_an_alias_for_production():
    assert parse_mdp(generate_mdp("md")) == parse_mdp(
        generate_mdp("production"))


def test_rendered_mdp_documents_its_purpose():
    text = generate_mdp("npt")
    assert text.startswith(";")
    assert "Purpose:" in text
    assert "force field" in text
