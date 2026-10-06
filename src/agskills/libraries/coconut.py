"""COCONUT 2.0 natural-product client.

COCONUT aggregates open natural-product collections, including the African
and marine sets relevant to anti-infective work (AfroDB, ANPDB, NANPDB,
SANCDB, CMNPD, ConMedNP).

The endpoint and request shape were verified live on 2026-10-06: the search
API accepts a JSON **POST** to ``/api/search`` and rejects GET with HTTP 405,
so the method matters. The previous implementation had this right, but it
bypassed the package's HTTP client entirely and used raw ``urllib`` with a
hand-rolled ``time.sleep(0.5)``, which meant no rate-limit coordination
between concurrent runs, no retry on a transient 5xx, and no robots
handling - while the skill documentation stated that the wrapper "enforces
rate limits and handles retries".

Reference:
    Sorokina M, Merseburger P, Rajan K, Yirik MA, Steinbeck C.
    *J Cheminform* 2021;13:2.
"""

from __future__ import annotations

import time
from typing import Any

from ..errors import InvalidInputError, RemoteServiceError
from ..http import HttpClient

__all__ = ["COCONUT_BASE", "search_coconut", "SOURCE_COLLECTIONS"]

COCONUT_BASE = "https://coconut.naturalproducts.net"
_SEARCH_PATH = "/api/search"

#: Page size ceiling the API honours.
MAX_PAGE_SIZE = 50

#: Collections aggregated by COCONUT that are commonly of interest here.
SOURCE_COLLECTIONS: dict[str, str] = {
    "AfroDB": "African medicinal plant compounds",
    "ANPDB": "African Natural Products Database",
    "NANPDB": "Northern African Natural Products Database",
    "SANCDB": "South African Natural Compounds Database",
    "CMNPD": "Comprehensive Marine Natural Products Database",
    "ConMedNP": "Congo medicinal plants",
}


def _normalise_compound(raw: dict[str, Any]) -> dict[str, Any]:
    """Map one COCONUT record onto a stable schema.

    COCONUT has changed field names between releases, so each value is
    looked up under the names observed across versions rather than one.
    """
    def first(*keys: str, default: Any = None) -> Any:
        for key in keys:
            value = raw.get(key)
            if value not in (None, "", []):
                return value
        return default

    return {
        "id": first("identifier", "coconut_id", "id", default=""),
        "name": first("name", "preferred_name", default=""),
        "smiles": first("canonical_smiles", "smiles", "SMILES", default=""),
        "iupac_name": first("iupac_name", "iupac", default=""),
        "inchikey": first("standard_inchi_key", "inchikey", default=""),
        "molecular_formula": first("molecular_formula", default=""),
        "molecular_weight": first("molecular_weight", "exact_molecular_weight"),
        "annotation_level": first("annotation_level", default=0),
        "organism_count": first("organism_count", default=0),
        "citation_count": first("citation_count", default=0),
        "np_likeness": first("np_likeness", "np_likeness_score"),
    }


def search_coconut(query: str | None = None, *, smiles: str | None = None,
                   limit: int = 25, page: int = 1,
                   max_pages: int = 1,
                   client: HttpClient | None = None) -> dict[str, Any]:
    """Search COCONUT by free text or by structure.

    Args:
        query: Free-text term - a compound name, an organism, or an activity
            keyword such as "antimycobacterial".
        smiles: A SMILES string for a structure search.
        limit: Results per page, capped at :data:`MAX_PAGE_SIZE`.
        page: First page to fetch (1-based).
        max_pages: How many consecutive pages to fetch. The previous
            multi-page helper could loop until it had made as many requests
            as the reported page count even when pages came back empty; this
            stops on the first empty page.

    Raises:
        InvalidInputError: If neither or both of *query* and *smiles* are given.
        RemoteServiceError: If the API is unreachable or returns an
            unexpected shape.
    """
    if bool(query) == bool(smiles):
        raise InvalidInputError(
            "Provide exactly one of query or smiles.",
            hint="Use --query for a text search or --smiles for a structure "
                 "search.",
        )
    if limit < 1:
        raise InvalidInputError("--limit must be at least 1")
    if page < 1:
        raise InvalidInputError("--page must be at least 1")

    term = (query or smiles or "").strip()
    page_size = min(limit, MAX_PAGE_SIZE)
    http = client or HttpClient(COCONUT_BASE, qps=2.0)

    compounds: list[dict[str, Any]] = []
    total_reported: int | None = None
    last_page: int | None = None
    pages_fetched = 0

    current = page
    while pages_fetched < max_pages and len(compounds) < limit:
        payload = http.fetch_json(
            _SEARCH_PATH,
            method="POST",
            json_body={"query": term, "limit": page_size, "page": current},
        )
        if not isinstance(payload, dict):
            raise RemoteServiceError(
                f"COCONUT returned {type(payload).__name__}, expected an object"
            )
        # The payload nests the page under "data", and the rows under
        # "data.data". Tolerate a flattened shape too.
        container = payload.get("data")
        if isinstance(container, dict):
            rows = container.get("data") or []
            total_reported = container.get("total", total_reported)
            last_page = container.get("last_page", last_page)
        elif isinstance(container, list):
            rows = container
        else:
            rows = payload.get("results") or []

        if not rows:
            break
        compounds.extend(_normalise_compound(r) for r in rows if isinstance(r, dict))
        pages_fetched += 1
        if last_page is not None and current >= last_page:
            break
        current += 1
        if pages_fetched < max_pages:
            time.sleep(0.1)

    compounds = compounds[:limit]
    with_structure = [c for c in compounds if c["smiles"]]

    return {
        "source": "COCONUT 2.0",
        "search_type": "structure" if smiles else "text",
        "query": term,
        "pages_fetched": pages_fetched,
        "first_page": page,
        "total_in_database": total_reported,
        "returned": len(compounds),
        "with_smiles": len(with_structure),
        "compounds": compounds,
        "note": (
            "annotation_level (0-5) reflects how well characterised an entry "
            "is; organism_count and citation_count indicate how often the "
            "compound has been reported. Prefer well-annotated, "
            "multiply-reported entries. Records without a SMILES cannot be "
            "docked or screened."
        ),
        "citation": "Sorokina M et al. J Cheminform 2021;13:2",
    }
