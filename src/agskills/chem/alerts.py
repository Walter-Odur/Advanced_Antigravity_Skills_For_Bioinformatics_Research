"""Structural-alert catalogues (PAINS, Brenk, NIH, ZINC).

The replaced code used two different PAINS configurations on two code
paths - ``PAINS_A|PAINS_B|PAINS_C`` with ``GetMatches`` in the docking
filter, but the combined ``PAINS`` catalogue with ``GetFirstMatch`` in the
ADMET fallback. The second form can only ever report one alert, so the same
molecule could be described as having one alert or four depending on which
subcommand the user happened to run.

Catalogues are built once and cached, because constructing a RDKit
``FilterCatalog`` is expensive relative to screening a single molecule and
screening runs iterate over thousands of compounds.

References:
    Baell JB, Holloway GA. *J Med Chem* 2010;53:2719-2740 (PAINS).
    Brenk R et al. *ChemMedChem* 2008;3:435-444 (Brenk).
"""

from __future__ import annotations

import functools
from dataclasses import dataclass, field
from typing import Any

__all__ = [
    "AlertHit",
    "AlertReport",
    "CATALOG_NAMES",
    "screen_alerts",
    "available_catalogs",
]

#: Catalogue identifiers accepted by :func:`screen_alerts`, mapped to the
#: RDKit ``FilterCatalogs`` members they enable.
CATALOG_NAMES: dict[str, tuple[str, ...]] = {
    # PAINS subsets A, B and C together, rather than the combined catalogue,
    # so that the subset a hit came from is recoverable.
    "pains": ("PAINS_A", "PAINS_B", "PAINS_C"),
    "brenk": ("BRENK",),
    "nih": ("NIH",),
    "zinc": ("ZINC",),
    "chembl": ("CHEMBL",),
}


@dataclass(frozen=True)
class AlertHit:
    """One structural alert matched by a molecule."""

    catalog: str
    description: str

    def as_dict(self) -> dict[str, str]:
        return {"catalog": self.catalog, "description": self.description}


@dataclass
class AlertReport:
    """All structural alerts found for one molecule."""

    hits: list[AlertHit] = field(default_factory=list)
    catalogs_screened: list[str] = field(default_factory=list)

    @property
    def n_alerts(self) -> int:
        return len(self.hits)

    def count(self, catalog: str) -> int:
        """Number of alerts from one catalogue."""
        return sum(1 for h in self.hits if h.catalog == catalog)

    def descriptions(self, catalog: str | None = None) -> list[str]:
        return [h.description for h in self.hits
                if catalog is None or h.catalog == catalog]

    def as_dict(self) -> dict[str, Any]:
        return {
            "n_alerts": self.n_alerts,
            "catalogs_screened": self.catalogs_screened,
            "counts": {c: self.count(c) for c in self.catalogs_screened},
            "hits": [h.as_dict() for h in self.hits],
        }


def available_catalogs() -> list[str]:
    """Catalogue names this RDKit build actually provides."""
    try:
        from rdkit.Chem import FilterCatalog
    except ImportError:  # pragma: no cover
        return []
    members = FilterCatalog.FilterCatalogParams.FilterCatalogs
    return [
        name for name, required in CATALOG_NAMES.items()
        if all(hasattr(members, r) for r in required)
    ]


@functools.lru_cache(maxsize=16)
def _catalog(name: str):
    """Build and cache the RDKit FilterCatalog for one catalogue name."""
    from rdkit.Chem import FilterCatalog

    try:
        required = CATALOG_NAMES[name]
    except KeyError:
        raise KeyError(
            f"Unknown alert catalogue {name!r}. "
            f"Available: {', '.join(sorted(CATALOG_NAMES))}"
        ) from None

    members = FilterCatalog.FilterCatalogParams.FilterCatalogs
    params = FilterCatalog.FilterCatalogParams()
    added = 0
    for member in required:
        if hasattr(members, member):
            params.AddCatalog(getattr(members, member))
            added += 1
    if added == 0:  # pragma: no cover - depends on the RDKit build
        raise KeyError(f"RDKit build provides none of {required} for {name!r}")
    return FilterCatalog.FilterCatalog(params)


def screen_alerts(mol, catalogs: tuple[str, ...] | list[str] = ("pains",)) -> AlertReport:
    """Screen *mol* against the named alert catalogues.

    Uses ``GetMatches`` so that **every** alert is reported, not just the
    first one.

    Args:
        mol: A sanitised RDKit molecule.
        catalogs: Catalogue names from :data:`CATALOG_NAMES`.

    Raises:
        ValueError: If *mol* is ``None``.
        KeyError: If a catalogue name is not recognised.
    """
    if mol is None:
        raise ValueError("screen_alerts() requires a valid molecule")

    report = AlertReport(catalogs_screened=list(catalogs))
    for name in catalogs:
        catalog = _catalog(name)
        for match in catalog.GetMatches(mol):
            report.hits.append(
                AlertHit(catalog=name, description=match.GetDescription())
            )
    return report
