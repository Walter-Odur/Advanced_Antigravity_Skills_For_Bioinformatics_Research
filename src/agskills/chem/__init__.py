"""Cheminformatics core: one implementation of every molecular property.

The code this replaces computed the same quantity in several places with
different conventions - molecular weight was ``Descriptors.MolWt`` in the
docking filter but ``Descriptors.ExactMolWt`` in the ADMET fallback, so
aspirin was reported as 180.16 Da on one path and 180.04 Da on the other
and compared against Lipinski's 500 Da cutoff either way. Everything now
flows through :mod:`agskills.chem.descriptors`.
"""

from __future__ import annotations

from .descriptors import MolecularProperties, compute_properties
from .smiles import (
    CompoundRecord,
    load_compounds,
    parse_smiles,
    standardize_smiles,
)

__all__ = [
    "MolecularProperties",
    "compute_properties",
    "CompoundRecord",
    "load_compounds",
    "parse_smiles",
    "standardize_smiles",
]
