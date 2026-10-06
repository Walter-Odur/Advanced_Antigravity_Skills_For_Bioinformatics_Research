"""Offline ADMET screening with an explicit, two-tier check registry.

Design
------
A **check** is one named, independently-reported test. A **tier** is a named
selection of checks together with the policy limits that tier applies. Both
are data, so the number of checks a tier performs is *derived* from the
registry and reported in the output - it can never drift away from the
documentation the way the previous ``14 checks``/``34 checks`` claims did
(the code behind those claims ran three checks).

Two tiers ship, matching the pipeline's two decision points:

``relaxed`` (tier 1)
    Applied after docking, to choose seeds for generative chemistry. Wide
    windows, tolerant of alerts: the goal is a diverse seed set, and the
    generator can optimise away a liability later.

``strict`` (tier 2)
    Applied to final candidates before molecular dynamics and reporting.
    Literature thresholds, no structural alerts.

Accounting
----------
Every compound ends in exactly one of ``passed``, ``failed`` or ``errored``,
and the three always sum to the number of inputs. The code this replaces
counted passes with ``.get("overall_pass", False)`` and failures with
``not .get("overall_pass", True)``, so a compound with an unparseable SMILES
was counted in neither: three inputs were reported as
``total=3, passed=2, failed=0``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Literal

from .alerts import AlertReport, screen_alerts
from .descriptors import MolecularProperties, compute_properties
from .rules import evaluate_rule
from .smiles import CompoundRecord, parse_smiles

__all__ = [
    "CheckOutcome",
    "Check",
    "Tier",
    "TIERS",
    "CHECKS",
    "screen_compound",
    "screen_compounds",
    "tier_summary",
]

Status = Literal["pass", "fail", "skip"]


@dataclass
class CheckOutcome:
    """Result of one check on one compound.

    ``gating`` separates a check that decides the verdict from one that
    merely raises a flag. The distinction is scientific, not cosmetic. The
    Pfizer 3/75 criterion, for instance, identifies compounds at elevated
    risk of in-vivo toxicity (logP above 3 together with TPSA below 75); it
    was published as a signal to pay closer attention, not as a filter.
    Gating on it rejects most lead-like chemistry - ibuprofen and the whole
    4-anilinoquinazoline kinase-inhibitor series fail it - which is exactly
    the chemistry a permissive triage pass exists to retain.
    """

    check_id: str
    label: str
    category: str
    status: Status
    value: Any = None
    detail: str = ""
    limit: str = ""
    citation: str = ""
    gating: bool = True

    @property
    def passed(self) -> bool:
        return self.status == "pass"

    @property
    def blocks(self) -> bool:
        """Whether this outcome prevents the compound from passing."""
        return self.gating and self.status == "fail"

    @property
    def advises(self) -> bool:
        """Whether this outcome is a flag raised without blocking."""
        return not self.gating and self.status == "fail"

    def as_dict(self) -> dict[str, Any]:
        out = {
            "check": self.check_id,
            "label": self.label,
            "category": self.category,
            "status": self.status,
            "role": "gating" if self.gating else "advisory",
            "value": self.value,
        }
        if self.limit:
            out["limit"] = self.limit
        if self.detail:
            out["detail"] = self.detail
        if self.citation:
            out["citation"] = self.citation
        return out


@dataclass(frozen=True)
class Check:
    """A named check, independent of any tier's policy."""

    check_id: str
    label: str
    category: str
    evaluate: Callable[[MolecularProperties, AlertReport, dict], CheckOutcome]
    citation: str = ""


# ---------------------------------------------------------------------------
# Check implementations
# ---------------------------------------------------------------------------


def _numeric_check(check_id: str, label: str, category: str, attr: str,
                   *, unit: str = "", citation: str = "") -> Check:
    """Build a check comparing one descriptor against tier-supplied bounds.

    The tier policy supplies ``min`` and/or ``max`` under the check id. A
    descriptor whose value cannot be computed yields ``skip``, never a
    silent pass.
    """

    def _evaluate(p: MolecularProperties, _alerts: AlertReport,
                  policy: dict) -> CheckOutcome:
        value = getattr(p, attr, None)
        bounds = policy.get(check_id) or {}
        low, high = bounds.get("min"), bounds.get("max")
        limit = _format_bounds(low, high, unit)
        if value is None:
            return CheckOutcome(check_id, label, category, "skip",
                                detail=f"{label} could not be computed",
                                limit=limit, citation=citation)
        problems = []
        if low is not None and value < low:
            problems.append(f"{value:g}{unit} < {low:g}{unit}")
        if high is not None and value > high:
            problems.append(f"{value:g}{unit} > {high:g}{unit}")
        rounded = round(value, 3) if isinstance(value, float) else value
        return CheckOutcome(
            check_id, label, category,
            "fail" if problems else "pass",
            value=rounded,
            detail="; ".join(problems),
            limit=limit,
            citation=citation,
        )

    return Check(check_id, label, category, _evaluate, citation)


def _format_bounds(low: Any, high: Any, unit: str) -> str:
    if low is not None and high is not None:
        return f"{low:g}{unit} to {high:g}{unit}"
    if high is not None:
        return f"<= {high:g}{unit}"
    if low is not None:
        return f">= {low:g}{unit}"
    return "no limit"


def _rule_check(check_id: str, rule_name: str, label: str) -> Check:
    """Build a check that applies a published rule from :mod:`.rules`.

    The rule reports its own, fixed violation count. The tier may permit a
    number of violations via ``policy[check_id]["max_violations"]``, which
    changes the verdict but never the reported violations.
    """

    def _evaluate(p: MolecularProperties, _alerts: AlertReport,
                  policy: dict) -> CheckOutcome:
        result = evaluate_rule(rule_name, p)
        allowed = (policy.get(check_id) or {}).get(
            "max_violations", result.max_allowed_violations
        )
        ok = result.n_violations <= allowed
        return CheckOutcome(
            check_id, label, "drug-likeness rule",
            "pass" if ok else "fail",
            value=result.n_violations,
            detail="; ".join(result.violations),
            limit=f"<= {allowed} violation(s)",
            citation=result.citation,
        )

    return Check(check_id, label, "drug-likeness rule", _evaluate)


def _alert_check(check_id: str, catalog: str, label: str, citation: str) -> Check:
    """Build a check counting structural alerts from one catalogue."""

    def _evaluate(_p: MolecularProperties, alerts: AlertReport,
                  policy: dict) -> CheckOutcome:
        if catalog not in alerts.catalogs_screened:
            return CheckOutcome(check_id, label, "structural alert", "skip",
                                detail=f"{catalog} catalogue was not screened",
                                citation=citation)
        count = alerts.count(catalog)
        allowed = (policy.get(check_id) or {}).get("max", 0)
        hits = alerts.descriptions(catalog)
        return CheckOutcome(
            check_id, label, "structural alert",
            "pass" if count <= allowed else "fail",
            value=count,
            detail="; ".join(hits[:5]) + (f" (+{len(hits) - 5} more)"
                                          if len(hits) > 5 else ""),
            limit=f"<= {allowed} alert(s)",
            citation=citation,
        )

    return Check(check_id, label, "structural alert", _evaluate, citation)


#: Every available check, keyed by id.
CHECKS: dict[str, Check] = {
    c.check_id: c
    for c in (
        # --- physicochemical windows -----------------------------------
        _numeric_check("mw", "Molecular weight", "physicochemical", "mw",
                       unit=" Da"),
        _numeric_check("logp", "Lipophilicity (clogP)", "physicochemical",
                       "logp",
                       citation="Wildman & Crippen, J Chem Inf Comput Sci 1999"),
        _numeric_check("tpsa", "Topological polar surface area",
                       "physicochemical", "tpsa", unit=" A^2",
                       citation="Ertl et al. J Med Chem 2000;43:3714-3717"),
        _numeric_check("hbd", "H-bond donors", "physicochemical", "hbd"),
        _numeric_check("hba", "H-bond acceptors", "physicochemical", "hba"),
        _numeric_check("rotatable_bonds", "Rotatable bonds",
                       "physicochemical", "rotatable_bonds"),
        _numeric_check("heavy_atoms", "Heavy atoms", "physicochemical",
                       "heavy_atoms"),
        _numeric_check("molar_refractivity", "Molar refractivity",
                       "physicochemical", "molar_refractivity"),
        _numeric_check("fraction_csp3", "Fraction of sp3 carbons",
                       "physicochemical", "fraction_csp3"),
        # --- published rules -------------------------------------------
        _rule_check("lipinski", "lipinski", "Lipinski rule of five"),
        _rule_check("veber", "veber", "Veber oral bioavailability"),
        _rule_check("egan", "egan", "Egan absorption"),
        _rule_check("ghose", "ghose", "Ghose drug-like range"),
        _rule_check("muegge", "muegge", "Muegge drug-likeness"),
        _rule_check("pfizer_3_75", "pfizer_3_75", "Pfizer 3/75 toxicity risk"),
        _rule_check("gsk_4_400", "gsk_4_400", "GSK 4/400 developability"),
        _rule_check("golden_triangle", "golden_triangle", "Golden triangle"),
        # --- structural alerts -----------------------------------------
        _alert_check("pains", "pains", "PAINS alerts",
                     "Baell & Holloway, J Med Chem 2010;53:2719-2740"),
        _alert_check("brenk", "brenk", "Brenk alerts",
                     "Brenk et al. ChemMedChem 2008;3:435-444"),
        _alert_check("nih", "nih", "NIH alerts",
                     "Doveston et al. / NIH filter set"),
        # --- synthesis and composite -----------------------------------
        _numeric_check("sa_score", "Synthetic accessibility", "synthesis",
                       "sa_score",
                       citation="Ertl & Schuffenhauer, J Cheminform 2009;1:8"),
        _numeric_check("qed", "Quantitative estimate of drug-likeness",
                       "composite", "qed",
                       citation="Bickerton et al. Nat Chem 2012;4:90-98"),
        _numeric_check("formal_charge", "Formal charge", "physicochemical",
                       "formal_charge"),
        _numeric_check("stereocenters", "Stereocentres", "synthesis",
                       "stereocenters"),
    )
}


@dataclass(frozen=True)
class Tier:
    """A named screening policy: which checks run, and with what limits."""

    name: str
    label: str
    purpose: str
    checks: tuple[str, ...]
    policy: dict[str, dict[str, Any]]
    alert_catalogs: tuple[str, ...]
    #: Checks that are computed and reported but do not decide the verdict.
    #: Used for published *risk flags*, which were never intended as
    #: exclusion filters.
    advisory: tuple[str, ...] = ()

    @property
    def n_checks(self) -> int:
        return len(self.checks)

    @property
    def n_gating(self) -> int:
        return len([c for c in self.checks if c not in self.advisory])

    def is_gating(self, check_id: str) -> bool:
        return check_id not in self.advisory

    def describe(self) -> dict[str, Any]:
        """A JSON-serialisable description, used as output provenance."""
        return {
            "name": self.name,
            "label": self.label,
            "purpose": self.purpose,
            "n_checks": self.n_checks,
            "n_gating_checks": self.n_gating,
            "n_advisory_checks": len(self.advisory),
            "checks": [
                {
                    "check": cid,
                    "label": CHECKS[cid].label,
                    "category": CHECKS[cid].category,
                    "role": "gating" if self.is_gating(cid) else "advisory",
                    "policy": self.policy.get(cid, {}),
                }
                for cid in self.checks
            ],
        }


_RELAXED = Tier(
    name="relaxed",
    label="Tier 1 (relaxed)",
    purpose=(
        "Post-docking triage. Selects a diverse seed set for generative "
        "chemistry; liabilities that a generator can optimise away are "
        "tolerated here."
    ),
    checks=(
        "mw", "logp", "tpsa", "hbd", "hba", "rotatable_bonds",
        "heavy_atoms", "formal_charge",
        "lipinski",
        "pains", "brenk",
        "sa_score", "qed",
        "pfizer_3_75",
    ),
    policy={
        "mw": {"min": 100, "max": 700},
        "logp": {"min": -3.0, "max": 7.0},
        "tpsa": {"max": 200},
        "hbd": {"max": 7},
        "hba": {"max": 15},
        "rotatable_bonds": {"max": 15},
        "heavy_atoms": {"min": 8, "max": 90},
        "formal_charge": {"min": -2, "max": 2},
        # Two rule-of-five violations are tolerated at this tier. The
        # violation count itself is still the true, literature count.
        "lipinski": {"max_violations": 2},
        "pains": {"max": 1},
        "brenk": {"max": 3},
        "sa_score": {"max": 8.0},
        "qed": {"min": 0.10},
    },
    alert_catalogs=("pains", "brenk"),
    # Reported so the liability is visible, but not gating: this tier's
    # whole purpose is to keep a diverse seed set, and gating on a
    # toxicity-risk flag would discard most lead-like chemistry.
    advisory=("pfizer_3_75",),
)

_STRICT = Tier(
    name="strict",
    label="Tier 2 (strict)",
    purpose=(
        "Final candidate selection before molecular dynamics and reporting. "
        "Published thresholds applied at face value, no structural alerts."
    ),
    checks=(
        "mw", "logp", "tpsa", "hbd", "hba", "rotatable_bonds",
        "heavy_atoms", "molar_refractivity", "fraction_csp3",
        "formal_charge", "stereocenters",
        "lipinski", "veber", "egan", "ghose", "muegge",
        "pfizer_3_75", "gsk_4_400", "golden_triangle",
        "pains", "brenk", "nih",
        "sa_score", "qed",
    ),
    policy={
        "mw": {"min": 150, "max": 500},
        "logp": {"min": -0.4, "max": 5.0},
        "tpsa": {"min": 20, "max": 140},
        "hbd": {"max": 5},
        "hba": {"max": 10},
        "rotatable_bonds": {"max": 10},
        "heavy_atoms": {"min": 20, "max": 70},
        "molar_refractivity": {"min": 40, "max": 130},
        "fraction_csp3": {"min": 0.10},
        "formal_charge": {"min": -1, "max": 1},
        "stereocenters": {"max": 4},
        "lipinski": {"max_violations": 1},
        "pains": {"max": 0},
        "brenk": {"max": 0},
        "nih": {"max": 0},
        "sa_score": {"max": 6.0},
        "qed": {"min": 0.30},
    },
    alert_catalogs=("pains", "brenk", "nih"),
    # Even for final candidates these two are advisory. Pfizer 3/75 is a
    # risk signal rather than a disqualification, and the golden triangle
    # is evaluated here with Crippen logP standing in for logD at pH 7.4,
    # so it is too approximate to reject a compound on.
    advisory=("pfizer_3_75", "golden_triangle"),
)

#: The available tiers, keyed by the value of ``--strictness``.
TIERS: dict[str, Tier] = {"relaxed": _RELAXED, "strict": _STRICT}


def tier_summary() -> dict[str, int]:
    """Number of checks each tier performs. Used by docs and tests."""
    return {name: tier.n_checks for name, tier in TIERS.items()}


# ---------------------------------------------------------------------------
# Screening
# ---------------------------------------------------------------------------


@dataclass
class CompoundScreen:
    """Screening result for one compound."""

    name: str
    smiles: str
    status: Literal["pass", "fail", "error"]
    properties: MolecularProperties | None = None
    outcomes: list[CheckOutcome] = field(default_factory=list)
    error: str = ""
    source: str = ""

    @property
    def n_passed(self) -> int:
        return sum(1 for o in self.outcomes if o.status == "pass")

    @property
    def n_failed(self) -> int:
        return sum(1 for o in self.outcomes if o.status == "fail")

    @property
    def n_skipped(self) -> int:
        return sum(1 for o in self.outcomes if o.status == "skip")

    @property
    def failed_checks(self) -> list[str]:
        """Gating failures - the reasons the compound did not pass."""
        return [f"{o.label}: {o.detail or 'failed'}"
                for o in self.outcomes if o.blocks]

    @property
    def advisories(self) -> list[str]:
        """Flags raised without blocking. Worth reporting, not disqualifying."""
        return [f"{o.label}: {o.detail or 'flagged'}"
                for o in self.outcomes if o.advises]

    def as_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "name": self.name,
            "smiles": self.smiles,
            "status": self.status,
            "overall_pass": self.status == "pass",
        }
        if self.source:
            out["source"] = self.source
        if self.error:
            out["error"] = self.error
        if self.properties is not None:
            out["properties"] = self.properties.as_dict()
        if self.outcomes:
            out["checks_passed"] = self.n_passed
            out["checks_failed"] = self.n_failed
            out["checks_skipped"] = self.n_skipped
            out["failed_checks"] = self.failed_checks
            out["advisories"] = self.advisories
            out["checks"] = [o.as_dict() for o in self.outcomes]
        return out


def screen_compound(record: CompoundRecord, tier: Tier) -> CompoundScreen:
    """Screen one compound against one tier.

    An unparseable structure returns ``status="error"`` rather than raising,
    so a batch run reports every input.
    """
    mol = parse_smiles(record.smiles)
    if mol is None:
        return CompoundScreen(
            name=record.name, smiles=record.smiles, status="error",
            error="invalid SMILES: RDKit could not parse the structure",
            source=record.source,
        )

    from rdkit import Chem

    canonical = Chem.MolToSmiles(mol)
    try:
        properties = compute_properties(mol, smiles=canonical)
        alerts = screen_alerts(mol, tier.alert_catalogs)
    except Exception as exc:  # pragma: no cover - defensive
        return CompoundScreen(
            name=record.name, smiles=canonical, status="error",
            error=f"{type(exc).__name__}: {exc}", source=record.source,
        )

    outcomes = []
    for cid in tier.checks:
        outcome = CHECKS[cid].evaluate(properties, alerts, tier.policy)
        outcome.gating = tier.is_gating(cid)
        outcomes.append(outcome)

    # Only a gating failure blocks. A skipped check is not a failure, but it
    # is not evidence of a pass either: it is reported and left out of the
    # verdict. An advisory failure is reported as a flag.
    verdict = "fail" if any(o.blocks for o in outcomes) else "pass"
    return CompoundScreen(
        name=record.name, smiles=canonical, status=verdict,
        properties=properties, outcomes=outcomes, source=record.source,
    )


def screen_compounds(records: list[CompoundRecord],
                     strictness: str = "strict") -> dict[str, Any]:
    """Screen a batch of compounds and return a complete, balanced report.

    Raises:
        KeyError: If *strictness* is not a known tier.
    """
    try:
        tier = TIERS[strictness]
    except KeyError:
        raise KeyError(
            f"Unknown strictness {strictness!r}. "
            f"Available: {', '.join(sorted(TIERS))}"
        ) from None

    screens = [screen_compound(r, tier) for r in records]
    passed = [s for s in screens if s.status == "pass"]
    failed = [s for s in screens if s.status == "fail"]
    errored = [s for s in screens if s.status == "error"]

    # Rank survivors by how comfortably they passed, then by QED.
    def _rank_key(s: CompoundScreen):
        qed = (s.properties.qed if s.properties and s.properties.qed is not None
               else 0.0)
        return (-s.n_passed, -qed)

    ordered = sorted(passed, key=_rank_key) + sorted(failed, key=_rank_key) + errored

    flagged = [s for s in screens if s.advisories]
    return {
        "engine": "RDKit (offline, deterministic)",
        "strictness": tier.name,
        "tier": tier.label,
        "tier_purpose": tier.purpose,
        "checks_per_compound": tier.n_checks,
        "gating_checks_per_compound": tier.n_gating,
        "advisory_checks_per_compound": len(tier.advisory),
        "tier_definition": tier.describe(),
        "total_compounds": len(screens),
        "passed": len(passed),
        "failed": len(failed),
        "errored": len(errored),
        "passed_with_advisories": len(
            [s for s in passed if s.advisories]),
        "compounds_with_advisories": len(flagged),
        "advisory_note": (
            "An advisory check is computed and reported but does not decide "
            "the verdict. These are published risk *flags* rather than "
            "exclusion criteria: gating on them would reject most lead-like "
            "chemistry. Read the per-compound 'advisories' field and judge "
            "whether the flag matters for your series."
            if tier.advisory else ""
        ),
        "compounds": [s.as_dict() for s in ordered],
    }
