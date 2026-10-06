"""ChEMBL client for target search and bioactivity retrieval.

Fixes over the code this replaces:

* **Pagination.** The old query passed ``limit`` straight to the API and
  took whatever came back. ChEMBL caps a page at 1000 records and returns
  the rest behind ``page_meta.next``, so any request for more than a page
  silently returned a truncated set that was then presented as the complete
  answer. :func:`query_activities` follows ``next`` until the requested
  count is satisfied.
* **Unit handling.** Activities were kept whenever ``pchembl_value`` was
  truthy and cast with a bare ``float()``, which raises on a malformed
  record and so aborts the whole query. Values are parsed defensively and
  bad records are reported, not fatal.
* **Duplicate activities.** ChEMBL frequently holds several measurements
  for the same compound against the same target. The old code returned
  them all, so a "top 20 compounds" list could be the same molecule twenty
  times. Results are now aggregated per molecule.

Reference:
    Zdrazil B et al. *Nucleic Acids Res* 2024;52:D1180-D1192.
"""

from __future__ import annotations

import statistics
from typing import Any

from ..errors import InvalidInputError, RemoteServiceError
from ..http import HttpClient

__all__ = [
    "CHEMBL_BASE",
    "search_targets",
    "query_activities",
    "ASSAY_TYPES",
]

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"

#: ChEMBL assay type codes.
ASSAY_TYPES: dict[str, str] = {
    "B": "binding - direct target affinity (Ki, Kd, IC50)",
    "F": "functional - a functional response in an assay",
    "A": "ADMET",
    "T": "toxicity",
    "P": "physicochemical",
    "U": "unclassified",
}

#: ChEMBL's hard page size.
_MAX_PAGE = 1000


def _client(client: HttpClient | None) -> HttpClient:
    return client or HttpClient(CHEMBL_BASE, qps=4.0)


def search_targets(query: str, *, limit: int = 10,
                   organism: str | None = None,
                   client: HttpClient | None = None) -> dict[str, Any]:
    """Search ChEMBL targets by free text.

    Args:
        query: Target name, gene symbol or keyword.
        limit: Maximum targets to return.
        organism: Optional case-insensitive substring filter on organism.

    Raises:
        InvalidInputError: If *query* is blank or *limit* is not positive.
    """
    if not (query or "").strip():
        raise InvalidInputError("--query cannot be empty")
    if limit < 1:
        raise InvalidInputError("--limit must be at least 1")

    http = _client(client)
    payload = http.fetch_json(
        f"target/search.json?q={_quote(query)}&limit={min(limit * 3, _MAX_PAGE)}"
    )
    targets = payload.get("targets", [])

    results = []
    for entry in targets:
        organism_name = entry.get("organism") or ""
        if organism and organism.lower() not in organism_name.lower():
            continue
        components = entry.get("target_components") or []
        accessions = [
            c.get("accession") for c in components if c.get("accession")
        ]
        results.append({
            "target_chembl_id": entry.get("target_chembl_id"),
            "pref_name": entry.get("pref_name"),
            "target_type": entry.get("target_type"),
            "organism": organism_name,
            "uniprot_accessions": accessions,
            "species_group_flag": entry.get("species_group_flag"),
            # The search endpoint's relevance score, useful for ranking.
            "score": entry.get("score"),
        })
        if len(results) >= limit:
            break

    return {
        "source": "ChEMBL",
        "query": query,
        "organism_filter": organism,
        "total_returned": len(results),
        "targets": results,
        "note": (
            "Pick a target_chembl_id from this list and pass it to "
            "'query-chembl --target-id'. The uniprot_accessions field feeds "
            "'prepare-receptor --uniprot-id'."
            if results else
            "No targets matched. Try a gene symbol (EGFR) or the full protein "
            "name, and check the organism filter."
        ),
        "citation": "Zdrazil B et al. Nucleic Acids Res 2024;52:D1180-D1192",
    }


def _quote(value: str) -> str:
    import urllib.parse
    return urllib.parse.quote(value.strip())


def _as_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    # ChEMBL occasionally carries nulls encoded as strings.
    if result != result:  # NaN
        return None
    return result


def query_activities(target_id: str, *, pchembl_min: float = 5.0,
                     assay_type: str = "B", limit: int = 25,
                     client: HttpClient | None = None) -> dict[str, Any]:
    """Retrieve bioactivities for a ChEMBL target, aggregated per molecule.

    Args:
        target_id: A ChEMBL target identifier such as ``CHEMBL203``.
        pchembl_min: Minimum pChEMBL value. 5 is 10 uM, 6 is 1 uM,
            7 is 100 nM; 6 or above is a reasonable hit threshold.
        assay_type: A code from :data:`ASSAY_TYPES`.
        limit: Maximum distinct molecules to return.

    Raises:
        InvalidInputError: On a malformed target id or assay type.
        RemoteServiceError: If ChEMBL cannot be reached.
    """
    target = (target_id or "").strip().upper()
    if not target.startswith("CHEMBL") or not target[6:].isdigit():
        raise InvalidInputError(
            f"{target_id!r} is not a ChEMBL target id.",
            hint="Ids look like CHEMBL203. Find one with "
                 "'search-chembl-target --query <name>'.",
        )
    assay = (assay_type or "B").strip().upper()
    if assay not in ASSAY_TYPES:
        raise InvalidInputError(
            f"Unknown assay type {assay_type!r}.",
            hint="Options: " + ", ".join(f"{k} = {v}"
                                         for k, v in ASSAY_TYPES.items()),
        )
    if limit < 1:
        raise InvalidInputError("--limit must be at least 1")

    http = _client(client)
    # Over-fetch, because several activities collapse into one molecule and
    # records without a structure or a pChEMBL value are dropped.
    target_records = min(max(limit * 8, 200), 5000)

    path = (
        f"activity.json?target_chembl_id={target}"
        f"&pchembl_value__gte={pchembl_min}"
        f"&assay_type={assay}"
        f"&limit={min(target_records, _MAX_PAGE)}"
    )

    activities: list[dict[str, Any]] = []
    pages = 0
    malformed = 0
    while path and len(activities) < target_records and pages < 20:
        payload = http.fetch_json(path)
        if not isinstance(payload, dict):
            raise RemoteServiceError(
                f"ChEMBL returned {type(payload).__name__}, expected an object"
            )
        batch = payload.get("activities")
        if batch is None:
            raise RemoteServiceError(
                "ChEMBL response contained no 'activities' key.",
                hint=f"Keys present: {', '.join(sorted(payload)) or 'none'}",
            )
        activities.extend(batch)
        pages += 1
        meta = payload.get("page_meta") or {}
        next_path = meta.get("next")
        # 'next' arrives as an absolute path; strip the API prefix so it
        # resolves against the client's base URL.
        path = next_path.replace("/chembl/api/data/", "") if next_path else None

    # Aggregate per molecule: one row per compound, summarising its measurements.
    per_molecule: dict[str, dict[str, Any]] = {}
    for record in activities:
        smiles = record.get("canonical_smiles")
        pchembl = _as_float(record.get("pchembl_value"))
        molecule_id = record.get("molecule_chembl_id")
        if not smiles or pchembl is None or not molecule_id:
            malformed += 1
            continue
        entry = per_molecule.setdefault(molecule_id, {
            "molecule_chembl_id": molecule_id,
            "canonical_smiles": smiles,
            "pchembl_values": [],
            "standard_types": set(),
            "assay_ids": set(),
            "units": set(),
        })
        entry["pchembl_values"].append(pchembl)
        if record.get("standard_type"):
            entry["standard_types"].add(record["standard_type"])
        if record.get("assay_chembl_id"):
            entry["assay_ids"].add(record["assay_chembl_id"])
        if record.get("standard_units"):
            entry["units"].add(record["standard_units"])

    compounds = []
    for entry in per_molecule.values():
        values = entry["pchembl_values"]
        compounds.append({
            "molecule_chembl_id": entry["molecule_chembl_id"],
            "canonical_smiles": entry["canonical_smiles"],
            "pchembl_median": round(statistics.median(values), 2),
            "pchembl_max": round(max(values), 2),
            "pchembl_min": round(min(values), 2),
            "n_measurements": len(values),
            "standard_types": sorted(entry["standard_types"]),
            "n_assays": len(entry["assay_ids"]),
            "units": sorted(entry["units"]),
        })

    # Rank on the median, which is robust to one optimistic measurement.
    compounds.sort(key=lambda c: (-c["pchembl_median"], -c["n_measurements"]))
    selected = compounds[:limit]

    return {
        "source": "ChEMBL",
        "query": {
            "target_chembl_id": target,
            "pchembl_min": pchembl_min,
            "assay_type": f"{assay} ({ASSAY_TYPES[assay]})",
            "limit": limit,
        },
        "activities_examined": len(activities),
        "pages_fetched": pages,
        "records_skipped": malformed,
        "distinct_molecules_found": len(compounds),
        "returned": len(selected),
        "compounds": selected,
        "interpretation": (
            "pChEMBL is -log10 of the molar activity: 6 is 1 uM, 7 is 100 nM, "
            "9 is 1 nM. Higher is more potent. pchembl_median aggregates "
            "repeated measurements, which is more reliable than a single "
            "best value."
        ),
        "citation": "Zdrazil B et al. Nucleic Acids Res 2024;52:D1180-D1192",
    }
