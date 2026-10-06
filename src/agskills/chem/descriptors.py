"""Molecular descriptors - the single source of truth for the whole package.

Conventions are fixed here and documented, because the quantities below are
reported to users and compared against published thresholds.

Molecular weight
    ``Descriptors.MolWt`` - the **average** molecular weight, which is what
    "molecular weight" means in Lipinski's rule of five and what PubChem's
    ``MolecularWeight`` field reports. ``ExactMolWt`` (the monoisotopic
    mass) is reported separately as ``exact_mw`` and is never substituted
    for ``mw``: for aspirin those are 180.16 and 180.04 Da respectively.

Hydrogen bond donors and acceptors
    Two conventions exist and both are reported, because they disagree and
    different rules were authored against different ones:

    * ``hbd`` / ``hba`` use RDKit's refined definitions
      (``CalcNumHBD`` / ``CalcNumHBA``). These reproduce PubChem's
      ``HBondDonorCount`` / ``HBondAcceptorCount``.
    * ``hbd_lipinski`` / ``hba_lipinski`` use Lipinski's own definitions  - 
      donors are the sum of OH and NH groups, acceptors the sum of all N
      and O atoms (``CalcNumLipinskiHBD`` / ``CalcNumLipinskiHBA``).

    The rule-of-five implementation in :mod:`agskills.chem.rules` uses the
    Lipinski counts, as the 1997 paper specifies. For metformin these
    differ substantially (3/1 refined versus 5/5 Lipinski), so the choice
    is not cosmetic.

logP
    Wildman-Crippen ``MolLogP``. Lipinski's paper states the cutoff for
    calculated logP as 5 (its MLOGP variant uses 4.15); the clogP form with
    a cutoff of 5 is the one in common use and the one applied here.

References:
    Lipinski CA et al. *Adv Drug Deliv Rev* 1997;23:3-25.
    Wildman SA, Crippen GM. *J Chem Inf Comput Sci* 1999;39:868-873.
    Ertl P, Rohde B, Selzer P. *J Med Chem* 2000;43:3714-3717 (TPSA).
    Ertl P, Schuffenhauer A. *J Cheminform* 2009;1:8 (synthetic accessibility).
    Bickerton GR et al. *Nat Chem* 2012;4:90-98 (QED).
"""

from __future__ import annotations

import functools
import os
import sys
from dataclasses import asdict, dataclass, field
from typing import Any

from ..errors import MissingDependencyError

__all__ = [
    "MolecularProperties",
    "compute_properties",
    "sa_score",
    "sa_score_available",
]


@dataclass
class MolecularProperties:
    """Computed descriptors for one molecule.

    ``None`` for a field means "could not be computed", never "zero".
    """

    smiles: str
    mw: float
    exact_mw: float
    heavy_atoms: int
    logp: float
    tpsa: float
    hbd: int
    hba: int
    hbd_lipinski: int
    hba_lipinski: int
    rotatable_bonds: int
    aromatic_rings: int
    rings: int
    heteroatoms: int
    carbons: int
    formal_charge: int
    fraction_csp3: float
    molar_refractivity: float
    qed: float | None = None
    sa_score: float | None = None
    formula: str = ""
    stereocenters: int = 0
    warnings: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        """Return a JSON-serialisable dict with floats rounded for reporting."""
        data = asdict(self)
        for key, places in (
            ("mw", 2), ("exact_mw", 4), ("logp", 2), ("tpsa", 2),
            ("fraction_csp3", 3), ("molar_refractivity", 2),
            ("qed", 3), ("sa_score", 2),
        ):
            value = data.get(key)
            if isinstance(value, float):
                data[key] = round(value, places)
                # Avoid "-0.0" appearing in reports.
                if data[key] == 0:
                    data[key] = abs(data[key])
        return data


@functools.lru_cache(maxsize=1)
def _sascorer():
    """Load RDKit's contributed synthetic-accessibility scorer.

    RDKit ships ``sascorer.py`` under ``RDConfig.RDContribDir`` rather than
    as an importable module, so the directory has to be put on ``sys.path``
    first. ``rdkit.Contrib.SA_Score`` happens to work on some builds and
    not others, so it is tried second rather than relied upon.
    """
    try:
        from rdkit.Chem import RDConfig
    except ImportError as exc:  # pragma: no cover
        raise MissingDependencyError("rdkit", install="pip install rdkit") from exc

    contrib = os.path.join(RDConfig.RDContribDir, "SA_Score")
    if os.path.isdir(contrib) and contrib not in sys.path:
        sys.path.append(contrib)
    try:
        import sascorer  # type: ignore[import-not-found]
        return sascorer
    except ImportError:
        pass
    try:
        from rdkit.Contrib.SA_Score import sascorer as contrib_sascorer  # type: ignore
        return contrib_sascorer
    except ImportError:
        return None


def sa_score_available() -> bool:
    """Whether the synthetic-accessibility scorer could be loaded."""
    return _sascorer() is not None


def sa_score(mol) -> float | None:
    """Ertl-Schuffenhauer synthetic accessibility, 1 (easy) to 10 (hard).

    Returns ``None`` if the scorer is unavailable, so a caller can say so
    rather than treating a missing score as a pass.
    """
    scorer = _sascorer()
    if scorer is None:
        return None
    try:
        return float(scorer.calculateScore(mol))
    except Exception:
        return None


def compute_properties(mol, *, smiles: str | None = None,
                       include_qed: bool = True,
                       include_sa: bool = True) -> MolecularProperties:
    """Compute every descriptor the package uses for one RDKit molecule.

    Args:
        mol: A sanitised RDKit ``Mol``.
        smiles: Canonical SMILES to record. Computed from *mol* if omitted.
        include_qed: Compute the quantitative estimate of drug-likeness.
        include_sa: Compute the synthetic-accessibility score.

    Raises:
        ValueError: If *mol* is ``None``. Callers must filter invalid
            structures first so that they can be reported individually.
    """
    if mol is None:
        raise ValueError("compute_properties() requires a valid molecule")

    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

    warnings: list[str] = []

    qed_value: float | None = None
    if include_qed:
        try:
            from rdkit.Chem import QED
            qed_value = float(QED.qed(mol))
        except Exception as exc:
            warnings.append(f"QED unavailable: {exc}")

    sa_value: float | None = None
    if include_sa:
        sa_value = sa_score(mol)
        if sa_value is None:
            warnings.append(
                "synthetic accessibility unavailable (RDKit SA_Score contrib "
                "module not found)"
            )

    try:
        stereocenters = len(Chem.FindMolChiralCenters(
            mol, includeUnassigned=True, useLegacyImplementation=False))
    except Exception:
        stereocenters = 0

    return MolecularProperties(
        smiles=smiles or Chem.MolToSmiles(mol),
        mw=float(Descriptors.MolWt(mol)),
        exact_mw=float(Descriptors.ExactMolWt(mol)),
        heavy_atoms=int(mol.GetNumHeavyAtoms()),
        logp=float(Crippen.MolLogP(mol)),
        tpsa=float(rdMolDescriptors.CalcTPSA(mol)),
        hbd=int(rdMolDescriptors.CalcNumHBD(mol)),
        hba=int(rdMolDescriptors.CalcNumHBA(mol)),
        hbd_lipinski=int(rdMolDescriptors.CalcNumLipinskiHBD(mol)),
        hba_lipinski=int(rdMolDescriptors.CalcNumLipinskiHBA(mol)),
        rotatable_bonds=int(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        aromatic_rings=int(rdMolDescriptors.CalcNumAromaticRings(mol)),
        rings=int(rdMolDescriptors.CalcNumRings(mol)),
        heteroatoms=int(rdMolDescriptors.CalcNumHeteroatoms(mol)),
        carbons=sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6),
        formal_charge=int(Chem.GetFormalCharge(mol)),
        fraction_csp3=float(rdMolDescriptors.CalcFractionCSP3(mol)),
        molar_refractivity=float(Crippen.MolMR(mol)),
        qed=qed_value,
        sa_score=sa_value,
        formula=rdMolDescriptors.CalcMolFormula(mol),
        stereocenters=stereocenters,
        warnings=warnings,
    )
