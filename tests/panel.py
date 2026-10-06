"""Shared reference molecules, importable by both conftest and test modules.

Kept separate from ``conftest.py`` so test modules can import the constants
for ``@pytest.mark.parametrize`` without importing the fixture module.
"""

from __future__ import annotations

#: Marketed drugs, spanning the property ranges the screening tiers care
#: about: small and simple (aspirin, ethanol), mid-sized and compliant
#: (imatinib), and deliberately non-compliant (rifampicin at 823 Da,
#: bedaquiline above clogP 7).
DRUG_SMILES: dict[str, str] = {
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "imatinib": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
    "bedaquiline": "CCN(CC)CC[C@@H](O)[C@@H](c1ccc2ccccc2n1)c1cc2ccccc2cc1OC",
    "rifampicin": (
        "C[C@H]1/C=C/C=C(\\C)C(=O)Nc2c(O)c3c(c(O)c2/C=N/N2CCN(C)CC2)"
        "C(=O)[C@](C)(O3)O/C=C/[C@@H](OC)[C@H](C)[C@@H](OC(C)=O)"
        "[C@H](C)[C@H](O)[C@H](C)C1"
    ),
    "penicillin_g": "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",
    "ethanol": "CCO",
    "benzene": "c1ccccc1",
}

#: Structures that must fail to parse, each for a different reason, so a
#: single over-permissive change cannot pass them all.
INVALID_SMILES: list[str] = [
    "",                 # empty
    "   ",              # whitespace only
    "NOT_A_MOLECULE",   # not SMILES at all
    "C(C(C",            # unbalanced parentheses
    "c1ccccc",          # unclosed aromatic ring
    "[Xx]",             # not an element symbol
    "C1CC",             # unclosed ring bond
]

#: Real 4-anilinoquinazoline EGFR inhibitors, as returned by ChEMBL for
#: CHEMBL203 at pChEMBL 8 or better. Used wherever a test needs genuine
#: medicinal chemistry rather than toy molecules.
EGFR_ACTIVES: list[str] = [
    "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC",
    "CN(C)c1cc2ncnc(Nc3cccc(Br)c3)c2cn1",
    "Nc1ccc2c(Nc3cccc(Br)c3)ncnc2c1",
    "CNc1cc2ncnc(Nc3cccc(Br)c3)c2cn1",
    "Cc1ccc(Nc2ncnc3cc(OC)c(OC)cc23)cc1",
]
