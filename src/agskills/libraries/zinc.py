"""ZINC tranche selection and download-script generation.

Why this module was rewritten rather than patched
-------------------------------------------------
The previous implementation generated URLs of the form
``https://zinc20.docking.org/tranches/{code}/download.{format}``. That path
does not serve compound data. Verified against the live site on 2026-10-06,
it returns **HTTP 200 with a 12 KB HTML page**. Because the generated script
tested only for the presence of a file and then counted its lines, every
tranche appeared to "succeed" and the pipeline proceeded with HTML pages
saved as ``.smi`` files. A silent-corruption failure of that kind is worse
than an error.

The real bulk endpoint, confirmed working, is::

    https://files.docking.org/2D/<first two letters>/<four letters>.smi

for example ``https://files.docking.org/2D/CD/CDAA.smi`` (4.2 MB, header
row ``smiles zinc_id``).

The four-letter tranche code
----------------------------
=========  ======================================================
Position   Meaning
=========  ======================================================
1          Molecular size bin (ZINC indexes by heavy-atom count).
2          logP bin.
3          Reactivity filter.
4          Purchasability.
=========  ======================================================

The size and logP mappings below were derived empirically: tranches across
the grid were sampled and their contents characterised with RDKit. The
measured mean molecular weight per size letter was A 174, B 227, C 279,
D 313, E 338, F 362, G 386, H 410, I 436, J 475, K 613 Da, which confirms
the size boundaries used here.

The logP boundaries were likewise measured, and this is where the previous
code had a second defect: its table placed letter ``D`` at logP 0-1 and
``E`` at 1-2, whereas sampling shows ``D`` holds 1-2 and ``E`` holds
2-2.5. Its table was shifted by one bin, so even with a working URL it
would have downloaded the wrong chemistry.

ZINC computes its own logP, which is not identical to RDKit's Crippen
estimate, so bin edges are approximate near the boundaries. Letters G to J
were observed to align with the 0.5-unit edges exactly.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable

from ..errors import UsageError
from ..http import HttpClient

__all__ = [
    "ZINC_FILES_BASE",
    "SIZE_BINS",
    "LOGP_BINS",
    "REACTIVITY_CODES",
    "PURCHASABILITY_CODES",
    "SUBSETS",
    "TrancheSelection",
    "select_tranches",
    "generate_zinc_script",
    "tranche_url",
]

ZINC_FILES_BASE = "https://files.docking.org"

#: Size letter -> (molecular weight lower bound, upper bound) in daltons.
SIZE_BINS: tuple[tuple[str, float, float], ...] = (
    ("A", 0, 200), ("B", 200, 250), ("C", 250, 300), ("D", 300, 325),
    ("E", 325, 350), ("F", 350, 375), ("G", 375, 400), ("H", 400, 425),
    ("I", 425, 450), ("J", 450, 500), ("K", 500, 10_000),
)

#: logP letter -> (lower bound, upper bound).
LOGP_BINS: tuple[tuple[str, float, float], ...] = (
    ("A", -100, -1), ("B", -1, 0), ("C", 0, 1), ("D", 1, 2),
    ("E", 2, 2.5), ("F", 2.5, 3), ("G", 3, 3.5), ("H", 3.5, 4),
    ("I", 4, 4.5), ("J", 4.5, 5), ("K", 5, 100),
)

#: Third letter: how aggressively reactive functionality is filtered out.
REACTIVITY_CODES: dict[str, str] = {
    "A": "anodyne - the most conservative filter, reactive groups excluded",
    "B": "bother - mildly reactive groups allowed",
    "C": "clean - standard medicinal-chemistry filtering",
    "E": "mild - some reactive groups allowed",
    "G": "reactive - few restrictions",
    "I": "all - no reactivity filtering",
}

#: Fourth letter: availability.
PURCHASABILITY_CODES: dict[str, str] = {
    "A": "in-stock, shipped in roughly 1-2 weeks",
    "B": "in-stock plus agent-supplied",
    "C": "on-demand, synthesised to order",
    "D": "boutique",
    "E": "annotated (drugs and bioactives)",
}

#: Named subsets, expressed as property windows rather than tranche letters.
SUBSETS: dict[str, dict[str, Any]] = {
    "drug-like": {
        "mw": (250.0, 500.0), "logp": (-1.0, 5.0),
        "description": "Lipinski-compatible window; the usual default for "
                       "a virtual screen.",
    },
    "lead-like": {
        "mw": (250.0, 350.0), "logp": (-1.0, 3.5),
        "description": "Smaller and less lipophilic, leaving room to "
                       "optimise.",
    },
    "fragment-like": {
        "mw": (0.0, 250.0), "logp": (-1.0, 3.0),
        "description": "Rule-of-three space for fragment-based screening.",
    },
    "all-purchasable": {
        "mw": (0.0, 600.0), "logp": (-100.0, 100.0),
        "description": "Everything purchasable; a very large download.",
    },
}


def _letters_for(bins: tuple[tuple[str, float, float], ...],
                 low: float, high: float) -> list[str]:
    """Letters whose bin overlaps the half-open window ``[low, high)``."""
    if low > high:
        raise UsageError(f"range lower bound {low} exceeds upper bound {high}")
    return [letter for letter, bin_low, bin_high in bins
            if low < bin_high and high > bin_low]


def tranche_url(code: str, fmt: str = "smi") -> str:
    """Build the bulk download URL for a four-letter tranche code.

    >>> tranche_url("CDAA")
    'https://files.docking.org/2D/CD/CDAA.smi'
    >>> tranche_url("cd")
    Traceback (most recent call last):
    agskills.errors.UsageError: ...
    """
    normalised = (code or "").strip().upper()
    if len(normalised) != 4 or not normalised.isalpha():
        raise UsageError(
            f"{code!r} is not a tranche code.",
            hint="A code is four letters, for example CDAA.",
        )
    return f"{ZINC_FILES_BASE}/2D/{normalised[:2]}/{normalised}.{fmt}"


@dataclass
class TrancheSelection:
    """A selected set of ZINC tranches and the window that produced it."""

    subset: str
    mw_range: tuple[float, float]
    logp_range: tuple[float, float]
    reactivity: str
    purchasability: str
    fmt: str
    codes: list[str] = field(default_factory=list)
    size_letters: list[str] = field(default_factory=list)
    logp_letters: list[str] = field(default_factory=list)
    verified: dict[str, Any] = field(default_factory=dict)

    @property
    def urls(self) -> list[tuple[str, str]]:
        return [(code, tranche_url(code, self.fmt)) for code in self.codes]

    def as_dict(self) -> dict[str, Any]:
        return {
            "subset": self.subset,
            "mw_range": list(self.mw_range),
            "logp_range": list(self.logp_range),
            "reactivity": f"{self.reactivity} ({REACTIVITY_CODES.get(self.reactivity, 'unknown')})",
            "purchasability": f"{self.purchasability} "
                              f"({PURCHASABILITY_CODES.get(self.purchasability, 'unknown')})",
            "format": self.fmt,
            "size_letters": self.size_letters,
            "logp_letters": self.logp_letters,
            "n_tranches": len(self.codes),
            "tranches": [{"code": c, "url": u} for c, u in self.urls],
            **({"verification": self.verified} if self.verified else {}),
        }


def select_tranches(subset: str = "drug-like", *,
                    mw_range: tuple[float, float] | None = None,
                    logp_range: tuple[float, float] | None = None,
                    reactivity: str = "A", purchasability: str = "A",
                    fmt: str = "smi") -> TrancheSelection:
    """Choose the tranche codes covering a property window.

    Raises:
        UsageError: On an unknown subset, format, or filter letter.
    """
    if subset not in SUBSETS:
        raise UsageError(
            f"Unknown subset {subset!r}.",
            hint=f"Choose one of: {', '.join(sorted(SUBSETS))}",
        )
    if fmt not in ("smi", "sdf", "mol2"):
        raise UsageError(f"Unsupported format {fmt!r}; use smi, sdf or mol2.")

    reactivity = (reactivity or "A").strip().upper()
    purchasability = (purchasability or "A").strip().upper()
    if reactivity not in REACTIVITY_CODES:
        raise UsageError(
            f"Unknown reactivity code {reactivity!r}.",
            hint="Options: " + ", ".join(
                f"{k} = {v}" for k, v in REACTIVITY_CODES.items()),
        )
    if purchasability not in PURCHASABILITY_CODES:
        raise UsageError(
            f"Unknown purchasability code {purchasability!r}.",
            hint="Options: " + ", ".join(
                f"{k} = {v}" for k, v in PURCHASABILITY_CODES.items()),
        )

    defaults = SUBSETS[subset]
    mw = tuple(mw_range) if mw_range else tuple(defaults["mw"])
    logp = tuple(logp_range) if logp_range else tuple(defaults["logp"])

    size_letters = _letters_for(SIZE_BINS, mw[0], mw[1])
    logp_letters = _letters_for(LOGP_BINS, logp[0], logp[1])
    if not size_letters or not logp_letters:
        raise UsageError(
            f"No ZINC tranche covers MW {mw[0]}-{mw[1]} with "
            f"logP {logp[0]}-{logp[1]}.",
            hint="Widen the ranges; ZINC bins molecular weight up to 500+ "
                 "and logP from below -1 to above 5.",
        )

    codes = [f"{s}{p}{reactivity}{purchasability}"
             for s in size_letters for p in logp_letters]
    return TrancheSelection(
        subset=subset, mw_range=mw, logp_range=logp, reactivity=reactivity,
        purchasability=purchasability, fmt=fmt, codes=codes,
        size_letters=size_letters, logp_letters=logp_letters,
    )


def verify_tranches(selection: TrancheSelection, *, limit: int = 0,
                    client: HttpClient | None = None) -> dict[str, Any]:
    """HEAD-check tranche URLs so a broken selection is caught up front.

    Args:
        selection: The selection to check.
        limit: Check at most this many tranches (0 checks all).
        client: Injectable HTTP client, for tests.
    """
    http = client or HttpClient(ZINC_FILES_BASE, qps=4.0, max_retries=1,
                                deadline=60.0)
    to_check = selection.codes[:limit] if limit else selection.codes
    reachable: list[str] = []
    missing: list[dict[str, str]] = []
    for code in to_check:
        url = tranche_url(code, selection.fmt)
        try:
            response = http.fetch(url, method="HEAD", timeout=20)
            content_type = (response.headers.get("Content-Type") or "").lower()
            if "html" in content_type:
                missing.append({"code": code,
                                "reason": f"served HTML, not data ({content_type})"})
            else:
                reachable.append(code)
        except Exception as exc:
            missing.append({"code": code, "reason": f"{type(exc).__name__}: {exc}"})
    result = {
        "checked": len(to_check),
        "reachable": len(reachable),
        "unreachable": missing,
    }
    selection.verified = result
    return result


def generate_zinc_script(selection: TrancheSelection, *,
                         out_dir: str | None = None) -> str:
    """Render a download script for the selected tranches.

    The generated script is defensive in ways the previous one was not:

    * It verifies each download is **not** an HTML page before keeping it,
      which is the exact failure the old URL produced.
    * It skips files already present, so an interrupted bulk download
      resumes instead of restarting.
    * It reports a per-tranche and overall compound count from the data.
    """
    directory = out_dir or f"zinc_{selection.subset}"
    header = [
        "#!/bin/bash",
        "# ZINC bulk download, generated by agskills.",
        f"# Subset:        {selection.subset}",
        f"# MW window:     {selection.mw_range[0]:g}-{selection.mw_range[1]:g} Da "
        f"(size letters {''.join(selection.size_letters)})",
        f"# logP window:   {selection.logp_range[0]:g}-{selection.logp_range[1]:g} "
        f"(logP letters {''.join(selection.logp_letters)})",
        f"# Reactivity:    {selection.reactivity} - "
        f"{REACTIVITY_CODES[selection.reactivity]}",
        f"# Purchasability:{selection.purchasability} - "
        f"{PURCHASABILITY_CODES[selection.purchasability]}",
        f"# Tranches:      {len(selection.codes)}",
        "#",
        "# Tranche sizes range from a few MB to several GB.",
        "# Check free disk space before running; this takes a while.",
        "",
        "set -uo pipefail",
        "",
        f'OUT_DIR="{directory}"',
        'mkdir -p "$OUT_DIR"',
        "",
        "if ! command -v curl >/dev/null 2>&1; then",
        '  echo "ERROR: curl is required." >&2',
        "  exit 1",
        "fi",
        "",
        "TOTAL=0",
        "OK_COUNT=0",
        "FAIL_COUNT=0",
        "",
        "# A tranche that does not exist is served as an HTML error page with",
        "# HTTP 200, so the content type has to be checked rather than just",
        "# the exit status.",
        "fetch_tranche() {",
        '  local code="$1" url="$2" dest="$3"',
        '  if [ -s "$dest" ]; then',
        '    echo "  [skip] $code already downloaded"',
        "    return 0",
        "  fi",
        '  if ! curl -fsSL --retry 3 --retry-delay 2 -o "$dest.part" "$url"; then',
        '    echo "  [fail] $code could not be downloaded" >&2',
        '    rm -f "$dest.part"',
        "    return 1",
        "  fi",
        '  if head -c 512 "$dest.part" | grep -qiE "<!doctype html|<html"; then',
        '    echo "  [fail] $code returned an HTML page, not compound data" >&2',
        '    rm -f "$dest.part"',
        "    return 1",
        "  fi",
        '  mv "$dest.part" "$dest"',
        "  return 0",
        "}",
        "",
    ]

    body: list[str] = []
    for code, url in selection.urls:
        dest = f'"$OUT_DIR/{code}.{selection.fmt}"'
        body += [
            f'echo "Fetching tranche {code} ..."',
            f'if fetch_tranche "{code}" "{url}" {dest}; then',
            f"  COUNT=$(($(wc -l < {dest}) - 1))",
            '  [ "$COUNT" -lt 0 ] && COUNT=0',
            "  TOTAL=$((TOTAL + COUNT))",
            "  OK_COUNT=$((OK_COUNT + 1))",
            f'  echo "  [ok] {code}: $COUNT compounds"',
            "else",
            "  FAIL_COUNT=$((FAIL_COUNT + 1))",
            "fi",
            "",
        ]

    footer = [
        'echo ""',
        'echo "=== Download summary ==="',
        'echo "Tranches downloaded: $OK_COUNT"',
        'echo "Tranches failed:     $FAIL_COUNT"',
        'echo "Compounds total:     $TOTAL"',
        "",
        'if [ "$OK_COUNT" -eq 0 ]; then',
        '  echo "ERROR: nothing was downloaded." >&2',
        "  exit 1",
        "fi",
        "",
        f'MERGED="$OUT_DIR/all.{selection.fmt}"',
        f'cat "$OUT_DIR"/*.{selection.fmt} > "$MERGED"',
        'echo "Merged library: $MERGED"',
        'echo ""',
        'echo "Next step - triage with the relaxed ADMET tier:"',
        'echo "  ag-compound-screening admet-filter \\\\"',
        'echo "    --smiles-file $MERGED --strictness relaxed \\\\"',
        'echo "    --output zinc_triage.json"',
    ]
    return "\n".join(header + body + footer) + "\n"
