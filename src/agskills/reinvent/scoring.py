"""Scoring profiles for REINVENT 4 multi-parameter optimisation.

Component names, parameter names and transform syntax here were taken from
the REINVENT 4 checkout itself (``configs/SCORING.md`` and the worked
example in ``configs/scoring_components_example.toml``) rather than from
memory, so a generated config matches what the installed version parses.

How transforms work, and why they matter
----------------------------------------
REINVENT needs every component to return a score in [0, 1] so that the
aggregate is meaningful. A raw property such as molecular weight has to be
mapped onto that range by a *transform*:

``double_sigmoid``
    A window with soft edges: good inside ``low``-``high``, falling off
    outside. The right choice for a property with a preferred range, such
    as molecular weight or TPSA.
``reverse_sigmoid``
    Lower is better. Used for synthetic accessibility and for logP when the
    aim is to reduce lipophilicity.
``sigmoid``
    Higher is better.
``step`` / ``left_step`` / ``right_step``
    Hard thresholds. Useful for a count that must not be exceeded, but they
    give the optimiser no gradient, so prefer a sigmoid where possible.

A profile that simply demanded "MW under 500" with a hard step would tell
the generator nothing about whether 505 or 900 is worse. Every numeric
component in these profiles therefore uses a smooth transform.
"""

from __future__ import annotations

from typing import Any

from ..errors import UsageError

__all__ = [
    "SCORING_PROFILES",
    "AGGREGATION_TYPES",
    "build_scoring_section",
    "parse_component_spec",
    "UNWANTED_SMARTS",
]

#: Aggregation functions REINVENT accepts.
AGGREGATION_TYPES: dict[str, str] = {
    "geometric_mean": "Weighted geometric mean. One near-zero component "
                      "drags the total down, so every objective must be "
                      "satisfied. The usual choice.",
    "arithmetic_mean": "Weighted arithmetic mean. Components can compensate "
                       "for one another.",
}

#: The structural alerts shipped in REINVENT's own example configuration:
#: strained rings, peroxides, carbocations, disulfides, and heteroatom-
#: heteroatom motifs that are typically unstable or reactive.
UNWANTED_SMARTS: list[str] = [
    "[*;r{8-17}]",
    "[#8][#8]",
    "[#6;+]",
    "[#16][#16]",
    "[#7;!n][S;!$(S(=O)=O)]",
    "[#7;!n][#7;!n]",
    "C#C",
    "C(=[O,S])[O,S]",
    "[#7;!n][C;!$(C(=[O,N])[N,O])][#16;!s]",
    "[#7;!n][C;!$(C(=[O,N])[N,O])][#7;!n]",
    "[#7;!n][C;!$(C(=[O,N])[N,O])][#8;!o]",
    "[#8;!o][C;!$(C(=[O,N])[N,O])][#16;!s]",
    "[#8;!o][C;!$(C(=[O,N])[N,O])][#8;!o]",
    "[#16;!s][C;!$(C(=[O,N])[N,O])][#16;!s]",
]


def _window(component: str, name: str, low: float, high: float, *,
            weight: float = 1.0, steepness: float = 20.0) -> dict[str, Any]:
    """A soft preferred-range transform for a numeric component."""
    return {
        "component": component,
        "endpoint": {
            "name": name,
            "weight": weight,
            "transform": {
                "type": "double_sigmoid",
                "high": high,
                "low": low,
                "coef_div": high if high else 1.0,
                "coef_si": steepness,
                "coef_se": steepness,
            },
        },
    }


def _lower_is_better(component: str, name: str, low: float, high: float, *,
                     weight: float = 1.0, k: float = 0.5) -> dict[str, Any]:
    return {
        "component": component,
        "endpoint": {
            "name": name,
            "weight": weight,
            "transform": {"type": "reverse_sigmoid", "high": high,
                          "low": low, "k": k},
        },
    }


def _higher_is_better(component: str, name: str, low: float, high: float, *,
                      weight: float = 1.0, k: float = 0.5) -> dict[str, Any]:
    return {
        "component": component,
        "endpoint": {
            "name": name,
            "weight": weight,
            "transform": {"type": "sigmoid", "high": high, "low": low, "k": k},
        },
    }


def _plain(component: str, name: str, weight: float = 1.0) -> dict[str, Any]:
    """A component that already returns a value in [0, 1], such as QED."""
    return {"component": component,
            "endpoint": {"name": name, "weight": weight}}


def _alerts(weight: float = 1.0) -> dict[str, Any]:
    return {
        "component": "custom_alerts",
        "endpoint": {"name": "Unwanted substructures", "weight": weight,
                     "params": {"smarts": list(UNWANTED_SMARTS)}},
    }


#: Named multi-objective profiles. Each is a list of component definitions.
SCORING_PROFILES: dict[str, dict[str, Any]] = {
    "drug-like": {
        "description": "Oral small-molecule drug space: Lipinski-compatible "
                       "size and lipophilicity, good QED, synthesisable.",
        "aggregation": "geometric_mean",
        "components": [
            _plain("QED", "QED drug-likeness", weight=2.0),
            _window("MolecularWeight", "Molecular weight", 200.0, 500.0),
            _window("SlogP", "Lipophilicity (SlogP)", 1.0, 4.0),
            _window("TPSA", "Polar surface area", 20.0, 140.0),
            _lower_is_better("NumRotBond", "Rotatable bonds", 2.0, 10.0),
            _lower_is_better("HBondDonors", "H-bond donors", 0.0, 5.0),
            _lower_is_better("HBondAcceptors", "H-bond acceptors", 0.0, 10.0),
            _lower_is_better("SAScore", "Synthetic accessibility", 1.0, 5.0,
                             weight=1.5),
            _alerts(),
        ],
    },
    "lead-like": {
        "description": "Smaller and less lipophilic than a finished drug, "
                       "leaving room for optimisation.",
        "aggregation": "geometric_mean",
        "components": [
            _plain("QED", "QED drug-likeness", weight=1.5),
            _window("MolecularWeight", "Molecular weight", 200.0, 350.0),
            _window("SlogP", "Lipophilicity (SlogP)", 0.0, 3.0),
            _lower_is_better("NumRotBond", "Rotatable bonds", 1.0, 7.0),
            _lower_is_better("HBondDonors", "H-bond donors", 0.0, 3.0),
            _lower_is_better("HBondAcceptors", "H-bond acceptors", 0.0, 6.0),
            _lower_is_better("SAScore", "Synthetic accessibility", 1.0, 4.0,
                             weight=1.5),
            _alerts(),
        ],
    },
    "fragment-like": {
        "description": "Rule-of-three fragment space for fragment-based "
                       "screening.",
        "aggregation": "geometric_mean",
        "components": [
            _window("MolecularWeight", "Molecular weight", 100.0, 250.0),
            _window("SlogP", "Lipophilicity (SlogP)", -1.0, 3.0),
            _lower_is_better("HBondDonors", "H-bond donors", 0.0, 3.0),
            _lower_is_better("HBondAcceptors", "H-bond acceptors", 0.0, 3.0),
            _lower_is_better("NumRotBond", "Rotatable bonds", 0.0, 3.0),
            _lower_is_better("SAScore", "Synthetic accessibility", 1.0, 3.5,
                             weight=1.5),
        ],
    },
    "kinase-inhibitor": {
        "description": "Typical ATP-competitive kinase inhibitor property "
                       "space: larger and flatter, with aromatic ring "
                       "systems for the hinge.",
        "aggregation": "geometric_mean",
        "components": [
            _plain("QED", "QED drug-likeness"),
            _window("MolecularWeight", "Molecular weight", 300.0, 550.0),
            _window("SlogP", "Lipophilicity (SlogP)", 1.0, 5.0),
            _window("TPSA", "Polar surface area", 50.0, 120.0),
            _window("NumAromaticRings", "Aromatic rings", 2.0, 4.0),
            _lower_is_better("HBondDonors", "H-bond donors", 1.0, 3.0),
            _lower_is_better("NumRotBond", "Rotatable bonds", 2.0, 8.0),
            _lower_is_better("SAScore", "Synthetic accessibility", 1.0, 5.0),
            _alerts(),
        ],
    },
    "anti-tb": {
        "description": (
            "Property space for anti-tubercular agents. The mycobacterial "
            "cell envelope is unusually impermeable and lipid-rich, so "
            "known actives sit at higher lipophilicity and molecular "
            "weight than the oral-drug average - bedaquiline is 555 Da "
            "with clogP above 7. This profile therefore permits a higher "
            "logP than 'drug-like' while still demanding synthesisability "
            "and no reactive groups."
        ),
        "aggregation": "geometric_mean",
        "components": [
            _window("MolecularWeight", "Molecular weight", 250.0, 600.0),
            _window("SlogP", "Lipophilicity (SlogP)", 2.0, 6.0),
            _window("TPSA", "Polar surface area", 20.0, 120.0),
            _plain("QED", "QED drug-likeness", weight=0.5),
            _lower_is_better("NumRotBond", "Rotatable bonds", 2.0, 10.0),
            _lower_is_better("SAScore", "Synthetic accessibility", 1.0, 4.5,
                             weight=1.5),
            _alerts(),
        ],
        "caveat": (
            "A property profile is not an activity model. Whole-cell "
            "M. tuberculosis potency depends on envelope permeability and "
            "efflux, which no physicochemical descriptor captures. Add a "
            "target-based component (docking through DockStream, or a QSAR "
            "model through ChemProp or Qptuna) before treating output as "
            "anti-tubercular leads."
        ),
    },
    "cns": {
        "description": "CNS-penetrant space: small, low polar surface area, "
                       "few donors, moderate lipophilicity.",
        "aggregation": "geometric_mean",
        "components": [
            _plain("QED", "QED drug-likeness"),
            _window("MolecularWeight", "Molecular weight", 200.0, 400.0),
            _window("SlogP", "Lipophilicity (SlogP)", 1.0, 4.0),
            _window("TPSA", "Polar surface area", 20.0, 76.0),
            _lower_is_better("HBondDonors", "H-bond donors", 0.0, 2.0),
            _lower_is_better("NumRotBond", "Rotatable bonds", 0.0, 7.0),
            _lower_is_better("SAScore", "Synthetic accessibility", 1.0, 4.5),
            _alerts(),
        ],
        "caveat": "Based on the CNS MPO physicochemical consensus; it does "
                  "not model active efflux at the blood-brain barrier.",
    },
}


def parse_component_spec(spec: str) -> dict[str, Any]:
    """Parse a ``--component`` command-line specification.

    Syntax: ``Component[:name=VALUE][:weight=W][:transform=T][:low=L][:high=H][:k=K]``

    >>> parsed = parse_component_spec("SlogP:low=1:high=4:weight=2")
    >>> parsed["component"], parsed["endpoint"]["weight"]
    ('SlogP', 2.0)
    >>> parse_component_spec("QED")["endpoint"]["name"]
    'QED'

    Raises:
        UsageError: If the specification cannot be parsed.
    """
    parts = [p for p in (spec or "").split(":") if p != ""]
    if not parts:
        raise UsageError("An empty --component specification was given.")
    component = parts[0].strip()
    if not component:
        raise UsageError(f"--component {spec!r} has no component name.")

    options: dict[str, str] = {}
    for chunk in parts[1:]:
        if "=" not in chunk:
            raise UsageError(
                f"--component {spec!r}: expected key=value, got {chunk!r}.",
                hint="For example: SlogP:low=1:high=4:weight=2",
            )
        key, _, value = chunk.partition("=")
        options[key.strip().lower()] = value.strip()

    def number(key: str) -> float | None:
        if key not in options:
            return None
        try:
            return float(options[key])
        except ValueError:
            raise UsageError(
                f"--component {spec!r}: {key} must be a number, got "
                f"{options[key]!r}"
            ) from None

    endpoint: dict[str, Any] = {
        "name": options.get("name", component),
        "weight": number("weight") if "weight" in options else 1.0,
    }
    low, high, k = number("low"), number("high"), number("k")
    transform_type = options.get("transform")
    if transform_type is None and (low is not None or high is not None):
        # A bounded range implies a window; a single bound implies a
        # direction.
        transform_type = "double_sigmoid" if (low is not None and
                                              high is not None) else "sigmoid"
    if transform_type:
        transform: dict[str, Any] = {"type": transform_type}
        if low is not None:
            transform["low"] = low
        if high is not None:
            transform["high"] = high
        if transform_type == "double_sigmoid":
            transform.setdefault("coef_div", high if high else 1.0)
            transform.setdefault("coef_si", 20.0)
            transform.setdefault("coef_se", 20.0)
        elif k is not None:
            transform["k"] = k
        endpoint["transform"] = transform
    return {"component": component, "endpoint": endpoint}


def build_scoring_section(profile: str | None = None,
                          components: list[str] | None = None,
                          aggregation: str | None = None) -> dict[str, Any]:
    """Build the scoring definition from a profile and/or custom components.

    Args:
        profile: A key from :data:`SCORING_PROFILES`.
        components: Extra ``--component`` specifications, appended to the
            profile's own components.
        aggregation: Override the aggregation function.

    Raises:
        UsageError: If neither a profile nor components are given, or on an
            unknown profile or aggregation function.
    """
    if not profile and not components:
        raise UsageError(
            "Scoring needs at least one objective.",
            hint="Pass --scoring-profile (one of: "
                 f"{', '.join(sorted(SCORING_PROFILES))}) or one or more "
                 "--component specifications.",
        )

    chosen: list[dict[str, Any]] = []
    metadata: dict[str, Any] = {}
    default_aggregation = "geometric_mean"

    if profile:
        key = profile.strip().lower()
        if key not in SCORING_PROFILES:
            raise UsageError(
                f"Unknown scoring profile {profile!r}.",
                hint="Available: " + ", ".join(
                    f"{n} - {s['description'][:60]}"
                    for n, s in SCORING_PROFILES.items()
                ),
            )
        spec = SCORING_PROFILES[key]
        chosen.extend(spec["components"])
        default_aggregation = spec["aggregation"]
        metadata = {
            "profile": key,
            "description": spec["description"],
            **({"caveat": spec["caveat"]} if spec.get("caveat") else {}),
        }

    for raw in components or []:
        chosen.append(parse_component_spec(raw))

    aggregation_type = (aggregation or default_aggregation).strip()
    if aggregation_type not in AGGREGATION_TYPES:
        raise UsageError(
            f"Unknown aggregation {aggregation_type!r}.",
            hint="Options: " + ", ".join(
                f"{k} - {v}" for k, v in AGGREGATION_TYPES.items()),
        )

    return {
        "type": aggregation_type,
        "components": chosen,
        "metadata": metadata,
        "n_components": len(chosen),
    }
