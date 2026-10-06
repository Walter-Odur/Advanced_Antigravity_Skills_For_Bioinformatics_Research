"""The two-tier ADMET engine: registry integrity, accounting, and tiers.

Three properties matter most here, and each corresponds to a defect in the
code this replaces:

* **The documented check count must be the performed check count.** The old
  documentation claimed 14 and 34 checks; the code ran three.
* **Every input must be accounted for.** The old counters used
  ``.get("overall_pass", False)`` for passes and
  ``not .get("overall_pass", True)`` for failures, so an invalid structure
  was counted in neither: three inputs reported ``passed=2, failed=0``.
* **A check that could not run must not count as a pass.**
"""

from __future__ import annotations

import json

import pytest

from agskills.chem.admet import (
    CHECKS,
    TIERS,
    CompoundScreen,
    screen_compound,
    screen_compounds,
    tier_summary,
)
from agskills.chem.smiles import CompoundRecord


# ---------------------------------------------------------------------------
# Registry integrity
# ---------------------------------------------------------------------------


def test_tier_check_counts_are_derived_from_the_registry():
    """The reported count must equal the number of checks actually run."""
    summary = tier_summary()
    assert set(summary) == {"relaxed", "strict"}
    for name, tier in TIERS.items():
        assert summary[name] == len(tier.checks)
        assert tier.n_checks == len(tier.checks)
        # And the count is what a screen reports.
        report = screen_compounds([CompoundRecord("x", "CCO")],
                                  strictness=name)
        assert report["checks_per_compound"] == len(tier.checks)
        performed = report["compounds"][0]
        total = (performed["checks_passed"] + performed["checks_failed"]
                 + performed["checks_skipped"])
        assert total == len(tier.checks)


def test_every_tier_check_exists_in_the_registry():
    """A tier cannot name a check that is not implemented."""
    for name, tier in TIERS.items():
        for check_id in tier.checks:
            assert check_id in CHECKS, f"{name} names unknown check {check_id}"


def test_every_tier_policy_key_is_a_check_the_tier_runs():
    """A policy entry for a check the tier does not run is dead configuration."""
    for name, tier in TIERS.items():
        for check_id in tier.policy:
            assert check_id in tier.checks, (
                f"{name} has a policy for {check_id}, which it does not run"
            )


def test_strict_is_a_superset_of_relaxed_in_spirit():
    """Strict must run more checks than relaxed, and share its core."""
    relaxed, strict = TIERS["relaxed"], TIERS["strict"]
    assert strict.n_checks > relaxed.n_checks
    core = {"mw", "logp", "hbd", "hba", "lipinski", "pains", "sa_score", "qed"}
    assert core <= set(relaxed.checks)
    assert core <= set(strict.checks)


def test_alert_catalogues_are_screened_for_every_alert_check():
    """An alert check must have its catalogue screened, or it would skip."""
    for name, tier in TIERS.items():
        for check_id in tier.checks:
            if CHECKS[check_id].category == "structural alert":
                assert check_id in tier.alert_catalogs, (
                    f"{name} runs the {check_id} check but does not screen "
                    "that catalogue"
                )


def test_strict_policy_is_tighter_than_relaxed_on_shared_windows():
    """Where both tiers bound a property, strict must not be looser."""
    relaxed, strict = TIERS["relaxed"], TIERS["strict"]
    for check_id in set(relaxed.policy) & set(strict.policy):
        r, s = relaxed.policy[check_id], strict.policy[check_id]
        if "max" in r and "max" in s:
            assert s["max"] <= r["max"], check_id
        if "min" in r and "min" in s:
            assert s["min"] >= r["min"], check_id
        if "max_violations" in r and "max_violations" in s:
            assert s["max_violations"] <= r["max_violations"], check_id


# ---------------------------------------------------------------------------
# Accounting
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("strictness", ["relaxed", "strict"])
def test_accounting_balances_with_invalid_input(strictness):
    """passed + failed + errored must equal the number of inputs.

    This is the regression test for the counter defect: the same three
    inputs previously produced total=3, passed=2, failed=0.
    """
    records = [
        CompoundRecord("ethanol", "CCO"),
        CompoundRecord("bogus", "NOT_A_MOLECULE"),
        CompoundRecord("benzene", "c1ccccc1"),
        CompoundRecord("also_bogus", "C(C(C"),
    ]
    report = screen_compounds(records, strictness=strictness)
    assert report["total_compounds"] == 4
    assert (report["passed"] + report["failed"]
            + report["errored"]) == report["total_compounds"]
    assert report["errored"] == 2
    assert len(report["compounds"]) == 4


def test_every_input_appears_in_the_output_exactly_once():
    """No compound may be silently dropped or duplicated."""
    names = [f"c{i}" for i in range(12)]
    records = [CompoundRecord(n, "CCO") for n in names[:6]]
    records += [CompoundRecord(n, "bad!!") for n in names[6:]]
    report = screen_compounds(records, strictness="relaxed")
    returned = [c["name"] for c in report["compounds"]]
    assert sorted(returned) == sorted(names)


def test_invalid_structure_is_an_error_not_a_pass():
    """An unparseable structure must never be reported as passing."""
    report = screen_compounds([CompoundRecord("bad", "@@@")],
                              strictness="relaxed")
    compound = report["compounds"][0]
    assert compound["status"] == "error"
    assert compound["overall_pass"] is False
    assert "invalid SMILES" in compound["error"]
    assert report["passed"] == 0


def test_a_skipped_check_does_not_count_as_a_pass(monkeypatch):
    """When a descriptor cannot be computed, the check skips and is reported.

    Silently treating an uncomputable check as satisfied is how a compound
    passes a filter it was never actually tested against.
    """
    monkeypatch.setattr("agskills.chem.descriptors._sascorer", lambda: None)
    report = screen_compounds([CompoundRecord("aspirin",
                                              "CC(=O)Oc1ccccc1C(=O)O")],
                              strictness="relaxed")
    compound = report["compounds"][0]
    sa_checks = [c for c in compound["checks"] if c["check"] == "sa_score"]
    assert len(sa_checks) == 1
    assert sa_checks[0]["status"] == "skip"
    assert compound["checks_skipped"] >= 1
    # A skip is not a pass, and not a failure either.
    assert sa_checks[0]["status"] not in ("pass", "fail")


# ---------------------------------------------------------------------------
# Tier behaviour
# ---------------------------------------------------------------------------


def test_relaxed_admits_small_drugs_that_strict_rejects(reference_properties):
    """The two tiers must genuinely differ, not just in name.

    Small marketed drugs such as aspirin and caffeine fall below the strict
    tier's size floor but are perfectly reasonable generative seeds, which
    is the whole purpose of a relaxed first pass.
    """
    records = [
        CompoundRecord(name, reference_properties[name]["SMILES"])
        for name in ("aspirin", "caffeine", "isoniazid")
    ]
    relaxed = screen_compounds(records, strictness="relaxed")
    strict = screen_compounds(records, strictness="strict")
    assert relaxed["passed"] > strict["passed"], (
        f"relaxed passed {relaxed['passed']}, strict passed "
        f"{strict['passed']}"
    )
    assert relaxed["passed"] == 3
    assert strict["passed"] == 0


def test_a_compound_passing_strict_also_passes_relaxed():
    """Strict must be the harder test for any given molecule.

    A molecule that satisfied the rigorous tier but failed the permissive
    one would mean the tiers were not ordered.
    """
    candidates = [
        "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC",
        "Cc1ccc(Nc2ncnc3cc(OC)c(OC)cc23)cc1",
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "c1ccc(Nc2ncnc3cc(OC)c(OC)cc23)cc1",
        "CN1CCN(c2ccc(Nc3ncnc4cc(OC)c(OC)cc34)cc2)CC1",
    ]
    records = [CompoundRecord(f"c{i}", s) for i, s in enumerate(candidates)]
    relaxed = {c["name"]: c["overall_pass"]
               for c in screen_compounds(records, "relaxed")["compounds"]}
    strict = {c["name"]: c["overall_pass"]
              for c in screen_compounds(records, "strict")["compounds"]}
    for name, passed_strict in strict.items():
        if passed_strict:
            assert relaxed[name], (
                f"{name} passed strict but failed relaxed: the tiers are not "
                "ordered"
            )


def test_structural_alert_can_be_the_sole_reason_for_failure():
    """Penicillin G fails strict only on the beta-lactam alert.

    This documents a real and important limitation: medicinal-chemistry
    alert sets flag the beta-lactam ring, yet penicillin is a cornerstone
    antibiotic. A screen must therefore report *which* check failed, so the
    user can judge whether the flag matters for their series.
    """
    penicillin = CompoundRecord(
        "penicillin_g",
        "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O")
    report = screen_compounds([penicillin], strictness="strict")
    compound = report["compounds"][0]
    assert compound["overall_pass"] is False
    failures = compound["failed_checks"]
    assert len(failures) == 1
    assert "alert" in failures[0].lower()
    # Everything else passed, so the molecule is otherwise well-behaved.
    assert compound["checks_passed"] == report["checks_per_compound"] - 1


def test_unknown_strictness_raises_with_the_valid_options():
    with pytest.raises(KeyError) as excinfo:
        screen_compounds([CompoundRecord("x", "CCO")], strictness="medium")
    message = str(excinfo.value)
    assert "medium" in message
    assert "relaxed" in message and "strict" in message


# ---------------------------------------------------------------------------
# Output contract
# ---------------------------------------------------------------------------


def test_report_is_json_serialisable_and_self_describing():
    """A report must record how it was produced, not just its verdicts."""
    report = screen_compounds([CompoundRecord("aspirin",
                                              "CC(=O)Oc1ccccc1C(=O)O")],
                              strictness="strict")
    json.dumps(report)
    for key in ("engine", "strictness", "tier", "tier_purpose",
                "checks_per_compound", "tier_definition", "total_compounds",
                "passed", "failed", "errored", "compounds"):
        assert key in report, key
    definition = report["tier_definition"]
    assert definition["n_checks"] == report["checks_per_compound"]
    assert len(definition["checks"]) == report["checks_per_compound"]


def test_each_failed_check_states_a_limit_and_a_measurement():
    """A failure must be actionable."""
    report = screen_compounds([CompoundRecord("ethanol", "CCO")],
                              strictness="strict")
    compound = report["compounds"][0]
    failed = [c for c in compound["checks"] if c["status"] == "fail"]
    assert failed
    for check in failed:
        assert check.get("limit"), check
        assert check.get("detail"), check


def test_numeric_checks_carry_citations_where_the_source_is_published():
    """A user must be able to trace a threshold to its paper."""
    report = screen_compounds([CompoundRecord("aspirin",
                                              "CC(=O)Oc1ccccc1C(=O)O")],
                              strictness="strict")
    checks = report["compounds"][0]["checks"]
    cited = [c for c in checks if c.get("citation")]
    # Every rule and alert check is published; the plain property windows
    # are project policy and need not be.
    rule_checks = [c for c in checks
                   if c["category"] in ("drug-likeness rule",
                                        "structural alert")]
    assert rule_checks
    for check in rule_checks:
        assert check.get("citation"), check["check"]
    assert len(cited) >= len(rule_checks)


def test_results_are_ranked_with_passing_compounds_first():
    """Ordering must put the useful compounds at the top."""
    records = [
        CompoundRecord("bad", "@@@"),
        CompoundRecord("too_small", "C"),
        CompoundRecord("good", "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC"),
    ]
    report = screen_compounds(records, strictness="relaxed")
    statuses = [c["status"] for c in report["compounds"]]
    # Passes before failures before errors.
    assert statuses.index("pass") < statuses.index("error")
    assert statuses[-1] == "error"


# ---------------------------------------------------------------------------
# Gating versus advisory
# ---------------------------------------------------------------------------


def test_a_risk_flag_does_not_block_a_lead_like_compound():
    """Pfizer 3/75 must flag without rejecting.

    The criterion catches logP above 3 together with TPSA below 75, which
    describes a large share of legitimate lead-like chemistry: ibuprofen
    and the whole 4-anilinoquinazoline kinase series trip it. It was
    published as a signal to look harder at toxicity, not as a filter, so
    gating on it would empty a triage pass of exactly the chemistry it is
    meant to retain.
    """
    series = CompoundRecord("quinazoline",
                            "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC")
    report = screen_compounds([series], strictness="relaxed")
    compound = report["compounds"][0]

    assert compound["overall_pass"] is True
    assert compound["failed_checks"] == []
    # The liability is still visible, just not disqualifying.
    assert any("Pfizer" in advisory for advisory in compound["advisories"])

    pfizer = next(c for c in compound["checks"] if c["check"] == "pfizer_3_75")
    assert pfizer["status"] == "fail"
    assert pfizer["role"] == "advisory"


def test_advisory_checks_are_declared_in_the_tier_definition():
    """A reader of the output must be able to tell which checks gate."""
    report = screen_compounds([CompoundRecord("x", "CCO")],
                              strictness="relaxed")
    definition = report["tier_definition"]
    assert definition["n_checks"] == report["checks_per_compound"]
    assert (definition["n_gating_checks"]
            + definition["n_advisory_checks"]) == definition["n_checks"]
    roles = {c["check"]: c["role"] for c in definition["checks"]}
    assert roles["pfizer_3_75"] == "advisory"
    assert roles["lipinski"] == "gating"
    assert report["advisory_note"]


def test_advisory_checks_are_a_subset_of_the_checks_run():
    """A tier cannot mark a check advisory that it does not run."""
    for name, tier in TIERS.items():
        for check_id in tier.advisory:
            assert check_id in tier.checks, (
                f"{name} marks {check_id} advisory but does not run it"
            )


def test_gating_failures_and_advisories_are_reported_separately():
    """The two lists must not be conflated.

    A compound that both fails a real limit and trips a flag must show the
    limit under ``failed_checks`` and the flag under ``advisories``, so the
    user can see which one actually caused the rejection.
    """
    # Very lipophilic and far too large: fails real windows, and also
    # trips the 3/75 flag.
    greasy = CompoundRecord(
        "greasy",
        "c1ccc(-c2ccc(-c3ccc(-c4ccc(-c5ccccc5)cc4)cc3)cc2)cc1")
    report = screen_compounds([greasy], strictness="relaxed")
    compound = report["compounds"][0]
    assert compound["overall_pass"] is False
    assert compound["failed_checks"], "a gating failure should be listed"
    assert all("Pfizer" not in reason
               for reason in compound["failed_checks"]), (
        "the advisory flag must not appear among the blocking reasons"
    )


def test_relaxed_retains_the_kinase_inhibitor_series():
    """A real generative seed set must survive the permissive tier.

    These are the actual 4-anilinoquinazolines ChEMBL returns for EGFR at
    pChEMBL 8 or better. If triage rejected them there would be nothing
    sensible left to seed a generative run with.
    """
    series = [
        "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC",
        "CN(C)c1cc2ncnc(Nc3cccc(Br)c3)c2cn1",
        "Nc1ccc2c(Nc3cccc(Br)c3)ncnc2c1",
        "CNc1cc2ncnc(Nc3cccc(Br)c3)c2cn1",
        "Cc1ccc(Nc2ncnc3cc(OC)c(OC)cc23)cc1",
    ]
    records = [CompoundRecord(f"egfr_{i}", s) for i, s in enumerate(series)]
    report = screen_compounds(records, strictness="relaxed")
    assert report["passed"] == len(series), (
        "relaxed triage rejected real EGFR actives: "
        f"{[c['failed_checks'] for c in report['compounds']]}"
    )


def test_screen_compound_returns_a_typed_result():
    screen = screen_compound(CompoundRecord("aspirin",
                                            "CC(=O)Oc1ccccc1C(=O)O"),
                             TIERS["relaxed"])
    assert isinstance(screen, CompoundScreen)
    assert screen.status == "pass"
    assert screen.properties is not None
    assert screen.n_passed == TIERS["relaxed"].n_checks


def test_canonical_smiles_is_returned_not_the_input_spelling():
    """Output SMILES must be canonical, so downstream de-duplication works."""
    report = screen_compounds([CompoundRecord("benzene", "C1=CC=CC=C1")],
                              strictness="relaxed")
    assert report["compounds"][0]["smiles"] == "c1ccccc1"


def test_screening_is_deterministic():
    """Two identical screens must agree exactly."""
    records = [CompoundRecord("a", "CC(=O)Oc1ccccc1C(=O)O"),
               CompoundRecord("b", "CN1C=NC2=C1C(=O)N(C)C(=O)N2C")]
    first = screen_compounds(records, "strict")
    second = screen_compounds(records, "strict")
    assert json.dumps(first, sort_keys=True) == json.dumps(second,
                                                            sort_keys=True)


def test_empty_input_produces_an_empty_but_valid_report():
    report = screen_compounds([], strictness="strict")
    assert report["total_compounds"] == 0
    assert report["passed"] == report["failed"] == report["errored"] == 0
    assert report["compounds"] == []
    assert report["checks_per_compound"] == TIERS["strict"].n_checks


def test_large_batch_completes_and_balances():
    """A realistic batch size must not degrade or lose compounds."""
    records = [CompoundRecord(f"c{i}", "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC")
               for i in range(200)]
    report = screen_compounds(records, strictness="strict")
    assert report["total_compounds"] == 200
    assert report["passed"] + report["failed"] + report["errored"] == 200
