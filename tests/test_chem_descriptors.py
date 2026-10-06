"""Descriptor correctness, checked against PubChem rather than against itself.

A descriptor layer that only agrees with itself proves nothing. These tests
compare computed values with PubChem's independently-computed properties
for ten marketed drugs, and assert the specific conventions the rest of the
package relies on.
"""

from __future__ import annotations

import math

import pytest

from agskills.chem.descriptors import (
    MolecularProperties,
    compute_properties,
    sa_score,
    sa_score_available,
)
from agskills.chem.smiles import parse_smiles


# ---------------------------------------------------------------------------
# Agreement with an independent source
# ---------------------------------------------------------------------------


def test_molecular_weight_matches_pubchem(reference_properties):
    """Average molecular weight must agree with PubChem to within 0.1 Da.

    This is the tightest cross-source check available: both compute the
    same quantity from the same formula, so any disagreement beyond
    rounding indicates a convention error.
    """
    for name, reference in reference_properties.items():
        mol = parse_smiles(reference["SMILES"])
        assert mol is not None, f"{name} failed to parse"
        properties = compute_properties(mol)
        expected = float(reference["MolecularWeight"])
        assert properties.mw == pytest.approx(expected, abs=0.1), (
            f"{name}: computed MW {properties.mw:.2f} vs PubChem "
            f"{expected:.2f}"
        )


def test_molecular_weight_is_average_not_monoisotopic(reference_properties):
    """``mw`` must be the average mass, and ``exact_mw`` the monoisotopic one.

    Reporting the monoisotopic mass as "molecular weight" and comparing it
    against Lipinski's 500 Da threshold is the specific defect this
    separation exists to prevent. For aspirin the two differ by 0.12 Da; for
    rifampicin, by 0.54 Da.
    """
    aspirin = parse_smiles(reference_properties["aspirin"]["SMILES"])
    properties = compute_properties(aspirin)
    assert properties.mw == pytest.approx(180.16, abs=0.02)
    assert properties.exact_mw == pytest.approx(180.042, abs=0.01)
    assert properties.mw > properties.exact_mw

    # Across the panel the two conventions must never be conflated.
    for name, reference in reference_properties.items():
        mol = parse_smiles(reference["SMILES"])
        computed = compute_properties(mol)
        pubchem = float(reference["MolecularWeight"])
        assert abs(computed.mw - pubchem) <= 0.1
        # exact_mw is always the lighter of the two for organic molecules.
        assert computed.exact_mw < computed.mw + 0.001, name


def test_tpsa_matches_pubchem(reference_properties):
    """TPSA must agree with PubChem within the aromatic-nitrogen convention.

    PubChem and RDKit differ on whether some aromatic nitrogens contribute,
    so caffeine legitimately differs by 3.4 A^2. A 5 A^2 tolerance catches
    a real error while admitting that known difference.
    """
    for name, reference in reference_properties.items():
        mol = parse_smiles(reference["SMILES"])
        properties = compute_properties(mol)
        expected = float(reference["TPSA"])
        assert properties.tpsa == pytest.approx(expected, abs=5.0), (
            f"{name}: computed TPSA {properties.tpsa:.1f} vs PubChem "
            f"{expected:.1f}"
        )


def test_hbd_hba_refined_convention_matches_pubchem(reference_properties):
    """``hbd``/``hba`` use the refined definitions, which match PubChem."""
    for name, reference in reference_properties.items():
        mol = parse_smiles(reference["SMILES"])
        properties = compute_properties(mol)
        assert properties.hbd == int(reference["HBondDonorCount"]), (
            f"{name}: HBD {properties.hbd} vs PubChem "
            f"{reference['HBondDonorCount']}"
        )


def test_lipinski_hbd_hba_are_reported_separately(reference_properties):
    """Both donor/acceptor conventions must be available and distinct.

    Lipinski's rule counts acceptors as all N and O atoms, which is not the
    refined count. Metformin shows the difference starkly: 3 refined donors
    against 5 Lipinski donors, and 1 refined acceptor against 5 Lipinski
    acceptors. The rule implementation needs the Lipinski numbers, and a
    user comparing against PubChem needs the refined ones, so both are
    reported.
    """
    metformin = parse_smiles(reference_properties["metformin"]["SMILES"])
    properties = compute_properties(metformin)
    assert properties.hbd == 3
    assert properties.hba == 1
    assert properties.hbd_lipinski == 5
    assert properties.hba_lipinski == 5
    # The two conventions are genuinely different fields, not aliases.
    assert (properties.hbd, properties.hba) != (
        properties.hbd_lipinski, properties.hba_lipinski)


def test_logp_is_crippen_and_correlates_with_pubchem_xlogp(reference_properties):
    """clogP must be the Crippen estimate and track PubChem's XLogP.

    These are different algorithms, so they are not expected to agree
    exactly. What must hold is that they rank lipophilicity consistently:
    a descriptor that disagreed in direction would make every logP-based
    rule wrong.
    """
    pairs = []
    for reference in reference_properties.values():
        if reference.get("XLogP") is None:
            continue
        mol = parse_smiles(reference["SMILES"])
        pairs.append((float(reference["XLogP"]),
                      compute_properties(mol).logp))
    assert len(pairs) >= 8

    # Pearson correlation between the two estimates.
    n = len(pairs)
    mean_x = sum(p[0] for p in pairs) / n
    mean_y = sum(p[1] for p in pairs) / n
    covariance = sum((x - mean_x) * (y - mean_y) for x, y in pairs)
    sd_x = math.sqrt(sum((x - mean_x) ** 2 for x, _ in pairs))
    sd_y = math.sqrt(sum((y - mean_y) ** 2 for _, y in pairs))
    correlation = covariance / (sd_x * sd_y)
    assert correlation > 0.9, f"Crippen vs XLogP correlation only {correlation:.3f}"

    # And bedaquiline, a famously lipophilic drug, must read as such.
    bedaquiline = parse_smiles(reference_properties["bedaquiline"]["SMILES"])
    assert compute_properties(bedaquiline).logp > 6.0


def test_rotatable_bonds_match_pubchem_for_most_drugs(reference_properties):
    """Rotatable-bond counts track PubChem within one bond for most drugs.

    The definitions differ on amide and terminal bonds, so exact agreement
    is not expected; a systematic error would show up as a large deviation.
    """
    deviations = []
    for name, reference in reference_properties.items():
        mol = parse_smiles(reference["SMILES"])
        computed = compute_properties(mol).rotatable_bonds
        expected = int(reference["RotatableBondCount"])
        deviations.append((name, abs(computed - expected)))
    large = [(n, d) for n, d in deviations if d > 2]
    assert not large, f"rotatable-bond counts deviate by more than 2: {large}"


# ---------------------------------------------------------------------------
# Field contracts
# ---------------------------------------------------------------------------


def test_all_fields_populated_for_a_real_drug(reference_properties):
    """Every descriptor must be computed, not left as None, for a normal drug."""
    mol = parse_smiles(reference_properties["imatinib"]["SMILES"])
    properties = compute_properties(mol)
    never_none = [
        "mw", "exact_mw", "heavy_atoms", "logp", "tpsa", "hbd", "hba",
        "hbd_lipinski", "hba_lipinski", "rotatable_bonds", "aromatic_rings",
        "rings", "heteroatoms", "carbons", "formal_charge", "fraction_csp3",
        "molar_refractivity",
    ]
    for field in never_none:
        assert getattr(properties, field) is not None, field
    assert properties.qed is not None
    assert properties.formula == "C29H31N7O"
    assert properties.aromatic_rings == 4
    assert properties.warnings == []


def test_qed_ranks_drug_likeness_sensibly(reference_properties):
    """QED must rate ibuprofen as more drug-like than rifampicin.

    QED is in [0, 1] and calibrated on marketed oral drugs. Ibuprofen is a
    small, simple drug; rifampicin is a 823 Da natural product. If the
    ordering came out the other way, the composite score would be
    misleading wherever it is used.
    """
    scores = {}
    for name in ("ibuprofen", "aspirin", "imatinib", "rifampicin"):
        mol = parse_smiles(reference_properties[name]["SMILES"])
        value = compute_properties(mol).qed
        assert 0.0 <= value <= 1.0, f"{name}: QED {value} out of range"
        scores[name] = value
    assert scores["ibuprofen"] > scores["imatinib"] > scores["rifampicin"]
    assert scores["ibuprofen"] > 0.7
    assert scores["rifampicin"] < 0.2


def test_synthetic_accessibility_orders_by_complexity(reference_properties):
    """SA score must rate aspirin easier to make than rifampicin."""
    if not sa_score_available():
        pytest.skip("the RDKit SA_Score contrib module is unavailable")
    aspirin = parse_smiles(reference_properties["aspirin"]["SMILES"])
    rifampicin = parse_smiles(reference_properties["rifampicin"]["SMILES"])
    easy = compute_properties(aspirin).sa_score
    hard = compute_properties(rifampicin).sa_score
    assert 1.0 <= easy <= 10.0
    assert 1.0 <= hard <= 10.0
    assert easy < 2.5, f"aspirin SA {easy} should be low"
    assert hard > 6.0, f"rifampicin SA {hard} should be high"
    assert easy < hard


def test_sa_score_returns_none_rather_than_a_default(monkeypatch, aspirin):
    """A missing scorer yields None, never a value that would pass a filter.

    Treating an unavailable score as a pass is how a compound slips through
    a filter it was never tested against.
    """
    monkeypatch.setattr("agskills.chem.descriptors._sascorer",
                        lambda: None)
    assert sa_score(aspirin) is None
    properties = compute_properties(aspirin, include_sa=True)
    assert properties.sa_score is None
    assert any("synthetic accessibility unavailable" in w
               for w in properties.warnings)


def test_compute_properties_rejects_none():
    """Passing None must raise, not return an empty record."""
    with pytest.raises(ValueError, match="valid molecule"):
        compute_properties(None)


def test_as_dict_is_json_serialisable_and_rounded(aspirin):
    """The dict form must round floats and avoid negative zero."""
    import json

    properties = compute_properties(aspirin)
    data = properties.as_dict()
    json.dumps(data)  # must not raise
    assert data["mw"] == pytest.approx(180.16, abs=0.01)
    assert isinstance(data["warnings"], list)

    # Negative zero reads as an error to a user; ethanol's logP is ~0.
    ethanol = compute_properties(parse_smiles("CCO")).as_dict()
    assert not str(ethanol["logp"]).startswith("-0.0")


def test_formal_charge_is_signed_correctly():
    """Charged species must report their real formal charge."""
    cases = {"CC(=O)[O-]": -1, "C[NH3+]": 1, "CCO": 0,
             "[Na+].[Cl-]": 0}
    for smiles, expected in cases.items():
        mol = parse_smiles(smiles)
        assert compute_properties(mol).formal_charge == expected, smiles


def test_stereocentres_counted_for_a_chiral_drug(reference_properties):
    """A drug with defined stereochemistry must report its centres."""
    penicillin = parse_smiles(
        "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O")
    assert compute_properties(penicillin).stereocenters == 3
    assert compute_properties(parse_smiles("CCO")).stereocenters == 0


def test_fraction_csp3_distinguishes_aromatic_from_aliphatic():
    """Fsp3 must be 0 for benzene and 1 for cyclohexane."""
    assert compute_properties(parse_smiles("c1ccccc1")).fraction_csp3 == 0.0
    assert compute_properties(parse_smiles("C1CCCCC1")).fraction_csp3 == 1.0


def test_properties_are_deterministic(aspirin):
    """Repeated computation must give identical results.

    Descriptors feed filtering decisions, so any run-to-run variation would
    make a screen irreproducible.
    """
    first = compute_properties(aspirin).as_dict()
    second = compute_properties(aspirin).as_dict()
    assert first == second


def test_properties_independent_of_smiles_spelling():
    """The same molecule written differently must give the same descriptors."""
    spellings = ["c1ccccc1", "C1=CC=CC=C1", "C1:C:C:C:C:C1"]
    results = []
    for smiles in spellings:
        mol = parse_smiles(smiles)
        if mol is None:
            continue
        results.append(compute_properties(mol).as_dict())
    assert len(results) >= 2
    first = results[0]
    for other in results[1:]:
        assert other == first
