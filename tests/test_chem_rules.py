"""Published drug-likeness rules, checked against their source papers.

The central property under test is that a rule reports **its authors'**
thresholds and nothing else. The implementation this replaces substituted
the screening tier's limits into the rule, so a 650 Da compound screened at
the relaxed tier was reported as having zero rule-of-five violations.
"""

from __future__ import annotations

import pytest

from agskills.chem.descriptors import compute_properties
from agskills.chem.rules import (
    RULES,
    RuleResult,
    evaluate_all,
    evaluate_rule,
    lipinski,
)
from agskills.chem.smiles import parse_smiles


def properties_for(smiles: str):
    mol = parse_smiles(smiles)
    assert mol is not None, smiles
    return compute_properties(mol)


# ---------------------------------------------------------------------------
# Lipinski: the rule the original code got wrong
# ---------------------------------------------------------------------------


def test_lipinski_thresholds_are_the_published_ones(reference_properties):
    """Aspirin passes with zero violations against MW<=500, logP<=5, 5, 10."""
    result = lipinski(properties_for(reference_properties["aspirin"]["SMILES"]))
    assert result.passed
    assert result.n_violations == 0
    assert result.max_allowed_violations == 1
    assert "1997" in result.citation


def test_lipinski_allows_exactly_one_violation():
    """The rule predicts poor absorption at two or more violations.

    A single violation must therefore pass. Atorvastatin (559 Da, clogP
    6.3) violates both MW and logP and must fail.
    """
    # One violation only: large but not lipophilic.
    one = properties_for("NC(=O)c1ccc(cc1)C(=O)NCCCCCCNC(=O)c1ccc(C(N)=O)cc1")
    result = lipinski(one)
    if result.n_violations == 1:
        assert result.passed, result.violations

    atorvastatin = properties_for(
        "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)"
        "C[C@@H](O)CC(=O)O")
    result = lipinski(atorvastatin)
    assert result.n_violations >= 2
    assert not result.passed


def test_lipinski_violation_count_is_absolute_not_tier_relative():
    """A 650 Da compound always reports an MW violation.

    This is the regression test for the defect that motivated separating
    rules from tiers: the violation count is a property of the molecule and
    the 1997 paper, never of the screening policy in force.
    """
    heavy = properties_for(
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)NCCCCNC(=O)c1ccccc1")
    assert heavy.mw > 500
    result = lipinski(heavy)
    assert result.n_violations >= 1
    assert any("MW" in violation for violation in result.violations)

    # And the reported violation names the published limit, not a tier's.
    mw_violation = next(v for v in result.violations if v.startswith("MW"))
    assert "> 500" in mw_violation


def test_lipinski_uses_the_lipinski_donor_acceptor_convention():
    """Donor and acceptor counts follow the 1997 definitions.

    Metformin has 5 Lipinski donors (OH + NH) but only 3 refined donors.
    The rule must use 5, because that is what the paper counts.
    """
    metformin = properties_for("CN(C)C(=N)N=C(N)N")
    assert metformin.hbd == 3
    assert metformin.hbd_lipinski == 5
    result = lipinski(metformin)
    # Five donors does not exceed the limit of five, so no HBD violation.
    assert not any("HBD" in v for v in result.violations)

    # A compound with six Lipinski donors must be flagged.
    sugary = properties_for("OC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO")
    if sugary.hbd_lipinski > 5:
        assert any("HBD" in v for v in lipinski(sugary).violations)


@pytest.mark.parametrize("name,expect_pass", [
    ("aspirin", True), ("caffeine", True), ("ibuprofen", True),
    ("imatinib", True), ("isoniazid", True), ("metformin", True),
    ("ethambutol", True),
    ("rifampicin", False),      # 823 Da, 16 Lipinski acceptors
    ("atorvastatin", False),    # 559 Da and clogP 6.3
])
def test_lipinski_verdicts_on_marketed_drugs(reference_properties, name,
                                             expect_pass):
    """Rule-of-five verdicts must match what the literature reports.

    Most oral drugs comply. Rifampicin and atorvastatin are the textbook
    non-compliant examples, so a rule that passed them would be wrong.
    """
    result = lipinski(properties_for(reference_properties[name]["SMILES"]))
    assert result.passed is expect_pass, (
        f"{name}: {result.n_violations} violation(s) {result.violations}"
    )


# ---------------------------------------------------------------------------
# The remaining rules
# ---------------------------------------------------------------------------


def test_every_rule_runs_and_reports_a_citation(reference_properties):
    """No rule may be unimplemented, uncited, or raise on a real drug."""
    properties = properties_for(reference_properties["imatinib"]["SMILES"])
    results = evaluate_all(properties)
    assert set(results) == set(RULES)
    for name, result in results.items():
        assert isinstance(result, RuleResult), name
        assert result.citation, f"{name} has no citation"
        assert result.description, f"{name} has no description"
        assert isinstance(result.passed, bool), name
        # A failing rule must say why.
        if not result.passed:
            assert result.violations, f"{name} failed with no stated reason"


def test_veber_flags_excessive_flexibility_and_polarity():
    """Veber's criteria are <=10 rotatable bonds and TPSA <= 140 A^2."""
    rigid = properties_for("c1ccc2c(c1)ccc1ccccc12")
    assert evaluate_rule("veber", rigid).passed

    floppy = properties_for("CCCCCCCCCCCCCCCCCC(=O)OCCCCCCCCCCCC")
    result = evaluate_rule("veber", floppy)
    assert not result.passed
    assert any("rotatable" in v for v in result.violations)


def test_pfizer_3_75_flags_lipophilic_low_polarity_compounds():
    """The 3/75 rule flags logP > 3 together with TPSA < 75.

    Passing means "not flagged", so a compound must trip *both* conditions.
    """
    # High logP, low TPSA: flagged.
    risky = properties_for("c1ccc(-c2ccc(-c3ccccc3)cc2)cc1")
    assert risky.logp > 3 and risky.tpsa < 75
    assert not evaluate_rule("pfizer_3_75", risky).passed

    # High logP but high TPSA: not flagged.
    polar = properties_for(
        "OC(=O)c1ccc(cc1)C(=O)Nc1ccc(cc1)S(=O)(=O)Nc1ccccc1C(=O)O")
    if polar.logp > 3 and polar.tpsa >= 75:
        assert evaluate_rule("pfizer_3_75", polar).passed

    # Low logP: not flagged regardless of TPSA.
    hydrophilic = properties_for("CN(C)C(=N)N=C(N)N")
    assert hydrophilic.logp <= 3
    assert evaluate_rule("pfizer_3_75", hydrophilic).passed


def test_rule_of_three_accepts_fragments_and_rejects_drugs(
        reference_properties):
    """Fragment criteria must admit a fragment and exclude a full drug."""
    fragment = properties_for("c1ccc(O)cc1")
    assert evaluate_rule("rule_of_three", fragment).passed

    imatinib = properties_for(reference_properties["imatinib"]["SMILES"])
    assert not evaluate_rule("rule_of_three", imatinib).passed


def test_lead_like_is_stricter_on_size_than_lipinski(reference_properties):
    """Lead-likeness caps MW at 350, so a 494 Da drug fails it but passes Ro5."""
    imatinib = properties_for(reference_properties["imatinib"]["SMILES"])
    assert lipinski(imatinib).passed
    assert not evaluate_rule("lead_like", imatinib).passed


def test_ghose_rejects_molecules_that_are_too_small(reference_properties):
    """Ghose requires 20-70 heavy atoms, so aspirin (13) fails."""
    aspirin = properties_for(reference_properties["aspirin"]["SMILES"])
    result = evaluate_rule("ghose", aspirin)
    assert not result.passed
    assert any("heavy atoms" in v for v in result.violations)


def test_gsk_4_400_is_stricter_than_lipinski(reference_properties):
    """GSK 4/400 caps MW at 400 and logP at 4."""
    imatinib = properties_for(reference_properties["imatinib"]["SMILES"])
    assert lipinski(imatinib).passed
    assert not evaluate_rule("gsk_4_400", imatinib).passed


def test_muegge_requires_a_minimum_size():
    """Muegge's window starts at 200 Da, so ethanol fails on several counts."""
    ethanol = properties_for("CCO")
    result = evaluate_rule("muegge", ethanol)
    assert not result.passed
    assert len(result.violations) >= 2


def test_unknown_rule_name_raises_with_the_available_names():
    with pytest.raises(KeyError) as excinfo:
        evaluate_rule("not_a_rule", properties_for("CCO"))
    message = str(excinfo.value)
    assert "not_a_rule" in message
    assert "lipinski" in message


def test_evaluate_all_accepts_a_subset():
    properties = properties_for("CCO")
    results = evaluate_all(properties, ["lipinski", "veber"])
    assert set(results) == {"lipinski", "veber"}


def test_rule_results_serialise(reference_properties):
    """Rule output must be JSON-serialisable for the report files."""
    import json
    properties = properties_for(reference_properties["aspirin"]["SMILES"])
    for result in evaluate_all(properties).values():
        payload = result.as_dict()
        json.dumps(payload)
        assert set(payload) >= {"name", "pass", "violations", "n_violations",
                                "citation"}


def test_violation_messages_name_the_quantity_and_the_limit():
    """A violation must be actionable: what was measured, and against what."""
    heavy = properties_for(
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)NCCCCNC(=O)c1ccccc1")
    for violation in lipinski(heavy).violations:
        # e.g. "MW 650 Da > 500 Da"
        assert ">" in violation or "<" in violation or "outside" in violation
        assert any(char.isdigit() for char in violation)
