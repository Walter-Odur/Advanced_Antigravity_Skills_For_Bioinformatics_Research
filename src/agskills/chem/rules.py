"""Published drug-likeness and lead-likeness rule sets.

Every rule is evaluated against **its own authors' thresholds**, which is
the central correctness fix in this rewrite. The code this replaces folded
the screening tier into the rule itself: in "relaxed" mode the molecular
weight limit became 700 Da, so a 650 Da compound was reported as having
``lipinski_violations: 0``. That is simply false - 650 > 500 is a rule-of-five
violation whatever the screening policy is.

The separation is now explicit:

* A **rule** is a fixed, cited statement about a molecule. Lipinski's limit
  is 500 Da permanently, and :func:`lipinski` always reports the true
  violation count.
* A **tier** is a project policy about which rules must pass. Tiers live in
  :mod:`agskills.chem.admet` and may tolerate violations, but they never
  redefine what a violation is.

References
----------
Lipinski CA, Lombardo F, Dominy BW, Feeney PJ. *Adv Drug Deliv Rev*
    1997;23:3-25.
Ghose AK, Viswanadhan VN, Wendoloski JJ. *J Comb Chem* 1999;1:55-68.
Egan WJ, Merz KM, Baldwin JJ. *J Med Chem* 2000;43:3867-3877.
Veber DF, Johnson SR, Cheng HY, Smith BR, Ward KW, Kopple KD.
    *J Med Chem* 2002;45:2615-2623.
Muegge I, Heald SL, Brittelli D. *J Med Chem* 2001;44:1841-1846.
Congreve M, Carr R, Murray C, Jhoti H. *Drug Discov Today* 2003;8:876-877
    (rule of three).
Teague SJ, Davis AM, Leeson PD, Oprea T. *Angew Chem Int Ed* 1999;38:3743-3748
    (lead-likeness).
Hughes JD et al. *Bioorg Med Chem Lett* 2008;18:4872-4875 (Pfizer 3/75).
Gleeson MP. *J Med Chem* 2008;51:817-834 (GSK 4/400).
Johnson TW, Dress KR, Edwards M. *Bioorg Med Chem Lett* 2009;19:5560-5564
    (golden triangle).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

from .descriptors import MolecularProperties

__all__ = [
    "RuleResult",
    "RULES",
    "evaluate_rule",
    "evaluate_all",
    "lipinski",
    "veber",
    "ghose",
    "egan",
    "muegge",
    "rule_of_three",
    "lead_like",
    "pfizer_3_75",
    "gsk_4_400",
    "golden_triangle",
]


@dataclass
class RuleResult:
    """The outcome of applying one published rule to one molecule."""

    name: str
    passed: bool
    violations: list[str] = field(default_factory=list)
    max_allowed_violations: int = 0
    citation: str = ""
    description: str = ""

    @property
    def n_violations(self) -> int:
        return len(self.violations)

    def as_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "pass": self.passed,
            "violations": self.violations,
            "n_violations": self.n_violations,
            "max_allowed_violations": self.max_allowed_violations,
            "citation": self.citation,
            "description": self.description,
        }


def _cmp_max(label: str, value: float | int | None, limit: float | int,
             unit: str = "") -> str | None:
    """Return a violation string if *value* exceeds *limit*."""
    if value is None:
        return None
    if value > limit:
        return f"{label} {value:g}{unit} > {limit:g}{unit}"
    return None


def _cmp_min(label: str, value: float | int | None, limit: float | int,
             unit: str = "") -> str | None:
    if value is None:
        return None
    if value < limit:
        return f"{label} {value:g}{unit} < {limit:g}{unit}"
    return None


def _cmp_range(label: str, value: float | int | None, low: float | int,
               high: float | int, unit: str = "") -> str | None:
    if value is None:
        return None
    if not (low <= value <= high):
        return f"{label} {value:g}{unit} outside {low:g}-{high:g}{unit}"
    return None


def _collect(*checks: str | None) -> list[str]:
    return [c for c in checks if c]


# ---------------------------------------------------------------------------
# The rules
# ---------------------------------------------------------------------------


def lipinski(p: MolecularProperties) -> RuleResult:
    """Lipinski's rule of five.

    Thresholds are fixed: MW <= 500 Da, clogP <= 5, HBD <= 5, HBA <= 10.
    Donor and acceptor counts use Lipinski's own definitions (OH + NH
    groups; N + O atoms). Poor absorption is predicted when **two or more**
    criteria are violated, so up to one violation passes.
    """
    violations = _collect(
        _cmp_max("MW", p.mw, 500, " Da"),
        _cmp_max("clogP", p.logp, 5),
        _cmp_max("HBD", p.hbd_lipinski, 5),
        _cmp_max("HBA", p.hba_lipinski, 10),
    )
    return RuleResult(
        name="Lipinski Ro5",
        passed=len(violations) <= 1,
        violations=violations,
        max_allowed_violations=1,
        citation="Lipinski et al. Adv Drug Deliv Rev 1997;23:3-25",
        description="Oral bioavailability; <=1 violation permitted.",
    )


def veber(p: MolecularProperties) -> RuleResult:
    """Veber's oral-bioavailability criteria: <=10 rotatable bonds, TPSA <=140."""
    violations = _collect(
        _cmp_max("rotatable bonds", p.rotatable_bonds, 10),
        _cmp_max("TPSA", p.tpsa, 140, " A^2"),
    )
    return RuleResult(
        name="Veber",
        passed=not violations,
        violations=violations,
        citation="Veber et al. J Med Chem 2002;45:2615-2623",
        description="Oral bioavailability from flexibility and polar surface area.",
    )


def ghose(p: MolecularProperties) -> RuleResult:
    """Ghose drug-like filter on MW, logP, molar refractivity and atom count."""
    violations = _collect(
        _cmp_range("MW", p.mw, 160, 480, " Da"),
        _cmp_range("logP", p.logp, -0.4, 5.6),
        _cmp_range("molar refractivity", p.molar_refractivity, 40, 130),
        _cmp_range("heavy atoms", p.heavy_atoms, 20, 70),
    )
    return RuleResult(
        name="Ghose",
        passed=not violations,
        violations=violations,
        citation="Ghose et al. J Comb Chem 1999;1:55-68",
        description="Qualifying range for known drugs.",
    )


def egan(p: MolecularProperties) -> RuleResult:
    """Egan "egg" absorption model: TPSA <= 131.6 and logP <= 5.88."""
    violations = _collect(
        _cmp_max("TPSA", p.tpsa, 131.6, " A^2"),
        _cmp_max("logP", p.logp, 5.88),
    )
    return RuleResult(
        name="Egan",
        passed=not violations,
        violations=violations,
        citation="Egan et al. J Med Chem 2000;43:3867-3877",
        description="Passive intestinal absorption boundary.",
    )


def muegge(p: MolecularProperties) -> RuleResult:
    """Muegge pharmacophore-point drug-likeness filter."""
    violations = _collect(
        _cmp_range("MW", p.mw, 200, 600, " Da"),
        _cmp_range("logP", p.logp, -2, 5),
        _cmp_max("TPSA", p.tpsa, 150, " A^2"),
        _cmp_max("rings", p.rings, 7),
        _cmp_min("carbon atoms", p.carbons, 5),
        _cmp_min("heteroatoms", p.heteroatoms, 2),
        _cmp_max("rotatable bonds", p.rotatable_bonds, 15),
        _cmp_max("HBA", p.hba, 10),
        _cmp_max("HBD", p.hbd, 5),
    )
    return RuleResult(
        name="Muegge",
        passed=not violations,
        violations=violations,
        citation="Muegge et al. J Med Chem 2001;44:1841-1846",
        description="Pharmacophore-point drug-likeness filter.",
    )


def rule_of_three(p: MolecularProperties) -> RuleResult:
    """Astex rule of three for fragment screening libraries."""
    violations = _collect(
        _cmp_max("MW", p.mw, 300, " Da"),
        _cmp_max("logP", p.logp, 3),
        _cmp_max("HBD", p.hbd, 3),
        _cmp_max("HBA", p.hba, 3),
        _cmp_max("rotatable bonds", p.rotatable_bonds, 3),
    )
    return RuleResult(
        name="Rule of three",
        passed=not violations,
        violations=violations,
        citation="Congreve et al. Drug Discov Today 2003;8:876-877",
        description="Fragment-likeness for fragment-based screening.",
    )


def lead_like(p: MolecularProperties) -> RuleResult:
    """Lead-likeness: a smaller, less lipophilic starting point than a drug."""
    violations = _collect(
        _cmp_range("MW", p.mw, 200, 350, " Da"),
        _cmp_range("logP", p.logp, -1, 3),
        _cmp_max("rotatable bonds", p.rotatable_bonds, 7),
    )
    return RuleResult(
        name="Lead-like",
        passed=not violations,
        violations=violations,
        citation="Teague et al. Angew Chem Int Ed 1999;38:3743-3748",
        description="Headroom for optimisation during lead development.",
    )


def pfizer_3_75(p: MolecularProperties) -> RuleResult:
    """Pfizer 3/75 in-vivo toxicity risk flag.

    A compound is flagged when logP > 3 **and** TPSA < 75; such compounds
    were roughly 2.5-fold more likely to show in-vivo toxicity. Passing
    means "not flagged".
    """
    risky = p.logp > 3 and p.tpsa < 75
    violations = (
        [f"logP {p.logp:.2f} > 3 and TPSA {p.tpsa:.1f} A^2 < 75 "
         "(elevated in-vivo toxicity risk)"]
        if risky else []
    )
    return RuleResult(
        name="Pfizer 3/75",
        passed=not risky,
        violations=violations,
        citation="Hughes et al. Bioorg Med Chem Lett 2008;18:4872-4875",
        description="Toxicity risk flag for lipophilic, low-polarity compounds.",
    )


def gsk_4_400(p: MolecularProperties) -> RuleResult:
    """GSK 4/400: MW <= 400 Da and logP <= 4 for a favourable ADMET profile."""
    violations = _collect(
        _cmp_max("MW", p.mw, 400, " Da"),
        _cmp_max("logP", p.logp, 4),
    )
    return RuleResult(
        name="GSK 4/400",
        passed=not violations,
        violations=violations,
        citation="Gleeson MP. J Med Chem 2008;51:817-834",
        description="Developability guideline on size and lipophilicity.",
    )


def golden_triangle(p: MolecularProperties) -> RuleResult:
    """Golden-triangle region of good permeability and clearance.

    The original criterion is stated in terms of logD at pH 7.4. Crippen
    logP is used as a surrogate here, which is a deliberate approximation
    and is reported as such.
    """
    violations = _collect(
        _cmp_range("MW", p.mw, 200, 500, " Da"),
        _cmp_range("logP (surrogate for logD7.4)", p.logp, -2, 5),
    )
    return RuleResult(
        name="Golden triangle",
        passed=not violations,
        violations=violations,
        citation="Johnson et al. Bioorg Med Chem Lett 2009;19:5560-5564",
        description=("Permeability/clearance sweet spot. Uses Crippen logP in "
                     "place of logD7.4."),
    )


#: Every rule, keyed by the identifier used in CLI flags and JSON output.
RULES: dict[str, Callable[[MolecularProperties], RuleResult]] = {
    "lipinski": lipinski,
    "veber": veber,
    "ghose": ghose,
    "egan": egan,
    "muegge": muegge,
    "rule_of_three": rule_of_three,
    "lead_like": lead_like,
    "pfizer_3_75": pfizer_3_75,
    "gsk_4_400": gsk_4_400,
    "golden_triangle": golden_triangle,
}


def evaluate_rule(name: str, p: MolecularProperties) -> RuleResult:
    """Evaluate a single rule by name.

    Raises:
        KeyError: If *name* is not a known rule.
    """
    try:
        return RULES[name](p)
    except KeyError:
        raise KeyError(
            f"Unknown rule {name!r}. Available: {', '.join(sorted(RULES))}"
        ) from None


def evaluate_all(p: MolecularProperties,
                 names: list[str] | None = None) -> dict[str, RuleResult]:
    """Evaluate every rule (or the named subset) against *p*."""
    selected = names if names is not None else list(RULES)
    return {name: evaluate_rule(name, p) for name in selected}
