#!/usr/bin/env python3
"""Generate the ADMET threshold reference from the live check registry.

Hand-written threshold tables drift. The documentation this replaces
claimed "14 checks" and "34 checks" while the implementation behind those
claims ran three. Generating the table from
:mod:`agskills.chem.admet` means the document cannot disagree with the code.

Run from the workshop root:

    python tools/generate_admet_reference.py

A test (``test_the_skill_documentation_states_the_real_check_counts``)
fails if a SKILL.md states a check count the registry does not perform, so
forgetting to regenerate is caught rather than shipped.
"""

from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
if (ROOT / "src" / "agskills" / "__init__.py").is_file():
    sys.path.insert(0, str(ROOT / "src"))

from agskills.chem.admet import CHECKS, TIERS  # noqa: E402

#: Skills whose reference directory receives a copy.
TARGET_SKILLS = ("compound-screening", "drug-discovery-wizard")


def _limit_for(tier, check_id: str) -> str:
    policy = tier.policy.get(check_id, {})
    if "max_violations" in policy:
        return f"<= {policy['max_violations']} violation(s)"
    low, high = policy.get("min"), policy.get("max")
    if low is not None and high is not None:
        return f"{low:g} to {high:g}"
    if high is not None:
        return f"<= {high:g}"
    if low is not None:
        return f">= {low:g}"
    return "rule default"


def render() -> str:
    lines = [
        "# ADMET Screening Reference",
        "",
        "**This file is generated from the check registry in**",
        "`src/agskills/chem/admet.py`. Do not edit it by hand: regenerate it",
        "with `python tools/generate_admet_reference.py`. It exists so the",
        "documented thresholds cannot drift away from the ones the code",
        "applies, which is what happened to the previous '14 checks' and",
        "'34 checks' claims - the implementation behind them ran three.",
        "",
        "## Tiers",
        "",
        "| Tier | Checks | Gating | Advisory | Purpose |",
        "|---|---:|---:|---:|---|",
    ]
    for name, tier in TIERS.items():
        lines.append(
            f"| `{name}` ({tier.label}) | {tier.n_checks} | {tier.n_gating} "
            f"| {len(tier.advisory)} | {tier.purpose} |"
        )
    lines += [
        "",
        "A **gating** check decides the verdict. An **advisory** check is",
        "computed and reported but does not. Advisory checks are published",
        "risk *flags* rather than exclusion criteria: the Pfizer 3/75",
        "criterion flags logP above 3 together with TPSA below 75, which",
        "describes much of lead-like space, so gating on it would reject",
        "most real medicinal chemistry.",
        "",
    ]

    for name, tier in TIERS.items():
        lines += [
            f"## Tier `{name}`: {tier.label}",
            "",
            tier.purpose,
            "",
            f"{tier.n_checks} checks per compound ({tier.n_gating} gating, "
            f"{len(tier.advisory)} advisory).",
            "",
            "| Check | Role | Category | Limit at this tier | Source |",
            "|---|---|---|---|---|",
        ]
        for check_id in tier.checks:
            check = CHECKS[check_id]
            role = "gating" if tier.is_gating(check_id) else "advisory"
            lines.append(
                f"| `{check_id}` | {role} | {check.category} | "
                f"{_limit_for(tier, check_id)} | "
                f"{check.citation or 'project policy'} |"
            )
        lines.append("")

    lines += [
        "## Conventions",
        "",
        "**Molecular weight** is the average mass (RDKit `MolWt`), which is",
        "what Lipinski's rule means and what PubChem reports. The",
        "monoisotopic mass is reported separately as `exact_mw` and is",
        "never substituted: for aspirin those are 180.16 and 180.04 Da.",
        "",
        "**Donors and acceptors** are reported under both conventions,",
        "because they disagree and different rules were authored against",
        "different ones. `hbd`/`hba` use RDKit's refined definitions, which",
        "reproduce PubChem. `hbd_lipinski`/`hba_lipinski` use Lipinski's own",
        "(OH + NH groups; all N and O atoms). The rule-of-five",
        "implementation uses the Lipinski counts, as the 1997 paper",
        "specifies. For metformin these differ substantially: 3/1 refined",
        "against 5/5 Lipinski.",
        "",
        "**logP** is the Wildman-Crippen estimate. Lipinski's paper states",
        "the calculated-logP cutoff as 5; that is the form applied here.",
        "",
        "**A rule always reports its authors' thresholds.** A tier changes",
        "which rules must pass and how many violations are tolerated; it",
        "never redefines what a violation is. A 650 Da compound reports a",
        "molecular-weight violation at every tier.",
        "",
        "## Accounting",
        "",
        "Every compound ends in exactly one of `passed`, `failed` or",
        "`errored`, and the three always sum to the number of inputs. A",
        "check that could not be computed is reported as `skip` and",
        "excluded from the verdict - it is never treated as a pass.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    content = render()
    for skill in TARGET_SKILLS:
        directory = ROOT / ".agents" / "skills" / skill / "reference"
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / "admet_thresholds.md"
        path.write_text(content, encoding="utf-8", newline="\n")
        print(f"wrote {path.relative_to(ROOT)}")
    counts = ", ".join(f"{n}: {t.n_checks}" for n, t in TIERS.items())
    print(f"check counts: {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
