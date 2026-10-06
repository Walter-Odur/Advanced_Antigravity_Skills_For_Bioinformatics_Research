"""REINVENT 4 TOML configuration generation.

TOML is emitted by a small writer rather than a library, because the shape
REINVENT expects uses arrays of tables in a specific nesting order
(``[[stage]]`` then ``[stage.scoring]`` then
``[[stage.scoring.component]]`` then
``[[stage.scoring.component.<Name>.endpoint]]``) and the generated file is
meant to be read and edited by the user. Keeping the writer here means the
output is ordered and commented.

Every generated config is validated before it is written: a missing seed
file for a generator that requires one, a transfer-learning run with no
input SMILES, or a sampling run with no model file are all caught at
generation time rather than minutes into a cluster job.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..errors import UsageError
from .priors import GENERATORS, find_priors_dir, find_reinvent_dir, resolve_prior
from .scoring import build_scoring_section

__all__ = [
    "RUN_MODES",
    "build_config",
    "render_toml",
    "toml_value",
]

#: REINVENT run modes and what each needs.
RUN_MODES: dict[str, dict[str, Any]] = {
    "staged_learning": {
        "description": "Reinforcement or curriculum learning. The main "
                       "optimisation mode.",
        "needs_scoring": True,
    },
    "transfer_learning": {
        "description": "Focus a prior on a set of known actives.",
        "needs_scoring": False,
    },
    "sampling": {
        "description": "Draw molecules from a model without optimising.",
        "needs_scoring": False,
    },
    "scoring": {
        "description": "Score an existing set of molecules.",
        "needs_scoring": True,
    },
}


def toml_value(value: Any) -> str:
    """Render a Python value as TOML.

    >>> toml_value(True)
    'true'
    >>> toml_value("priors/reinvent.prior")
    '"priors/reinvent.prior"'
    >>> toml_value([1, 2])
    '[1, 2]'
    >>> toml_value(0.0001)
    '0.0001'
    """
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int,)):
        return str(value)
    if isinstance(value, float):
        # Avoid scientific notation, which REINVENT's TOML parser accepts
        # but which is harder for a user to read and edit.
        text = f"{value:.10f}".rstrip("0").rstrip(".")
        return text or "0"
    if isinstance(value, Path):
        return toml_value(str(value))
    if isinstance(value, str):
        escaped = value.replace("\\", "\\\\").replace('"', '\\"')
        return f'"{escaped}"'
    if isinstance(value, (list, tuple)):
        if not value:
            return "[]"
        if all(isinstance(v, str) for v in value) and len(value) > 3:
            inner = ",\n    ".join(toml_value(v) for v in value)
            return f"[\n    {inner},\n]"
        return "[" + ", ".join(toml_value(v) for v in value) + "]"
    raise UsageError(f"Cannot render {type(value).__name__} as TOML")


@dataclass
class ConfigResult:
    """A generated configuration plus everything needed to act on it."""

    toml: str
    run_type: str
    generator: str
    settings: dict[str, Any] = field(default_factory=dict)
    expected_outputs: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    next_steps: list[str] = field(default_factory=list)
    scoring_metadata: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return {
            "run_type": self.run_type,
            "generator": self.generator,
            "settings": self.settings,
            "expected_outputs": self.expected_outputs,
            "scoring": self.scoring_metadata,
            "warnings": self.warnings,
            "next_steps": self.next_steps,
        }


def _emit_scoring(lines: list[str], scoring: dict[str, Any],
                  prefix: str) -> None:
    """Append a scoring section under *prefix* (``scoring`` or ``stage.scoring``)."""
    lines += [
        "",
        f"[{prefix}]",
        f"type = {toml_value(scoring['type'])}",
        "",
    ]
    for entry in scoring["components"]:
        component = entry["component"]
        endpoint = entry["endpoint"]
        lines += [
            f"[[{prefix}.component]]",
            f"[{prefix}.component.{component}]",
            f"[[{prefix}.component.{component}.endpoint]]",
            f"name = {toml_value(endpoint['name'])}",
            f"weight = {toml_value(endpoint.get('weight', 1.0))}",
        ]
        transform = endpoint.get("transform")
        if transform:
            for key in ("type", "high", "low", "k", "coef_div", "coef_si",
                        "coef_se"):
                if key in transform:
                    lines.append(
                        f"transform.{key} = {toml_value(transform[key])}"
                    )
        for key, value in (endpoint.get("params") or {}).items():
            lines.append(f"params.{key} = {toml_value(value)}")
        lines.append("")


def build_config(
    *,
    mode: str,
    generator: str = "reinvent",
    prior: str | None = None,
    device: str | None = None,
    scoring_profile: str | None = None,
    components: list[str] | None = None,
    aggregation: str | None = None,
    smiles_file: str | None = None,
    validation_smiles_file: str | None = None,
    model_file: str | None = None,
    num_steps: int | None = None,
    max_steps: int | None = None,
    min_steps: int | None = None,
    max_score: float | None = None,
    num_epochs: int | None = None,
    num_smiles: int | None = None,
    batch_size: int | None = None,
    sigma: int | None = None,
    learning_rate: float | None = None,
    csv_prefix: str | None = None,
    output_csv: str | None = None,
    output_model: str | None = None,
    checkpoint_file: str | None = None,
    diversity_filter: bool = True,
    inception_smiles_file: str | None = None,
    reinvent_dir: str | None = None,
    sample_strategy: str | None = None,
) -> ConfigResult:
    """Build a REINVENT 4 TOML configuration.

    Raises:
        UsageError: On an unknown mode or generator, a generator that needs
            seed structures without a ``smiles_file``, or a mode whose
            required parameters are missing.
    """
    run_type = (mode or "").strip().lower()
    if run_type not in RUN_MODES:
        raise UsageError(
            f"Unknown run mode {mode!r}.",
            hint="Options: " + ", ".join(
                f"{k} - {v['description']}" for k, v in RUN_MODES.items()),
        )
    generator_key = (generator or "reinvent").strip().lower()
    if generator_key not in GENERATORS:
        raise UsageError(
            f"Unknown generator {generator!r}.",
            hint="Options: " + ", ".join(GENERATORS),
        )

    spec = GENERATORS[generator_key]
    warnings: list[str] = []
    next_steps: list[str] = []

    directory = find_reinvent_dir(reinvent_dir)
    priors_dir = find_priors_dir(directory)
    prior_name, prior_path = resolve_prior(generator_key, prior,
                                           priors_dir=priors_dir)
    # Reference the prior by a path REINVENT can resolve. When the file was
    # found, use its real location; otherwise emit the conventional
    # relative path and warn.
    prior_reference = str(prior_path) if prior_path else f"priors/{Path(prior_name).name}"
    if prior_path is None:
        warnings.append(
            f"The prior model file {prior_name!r} was not found locally, so "
            f"the config references {prior_reference!r}. Download the priors "
            "or run 'check-setup' to see what is missing."
        )

    # A generator that consumes seed structures cannot run without them.
    if spec["needs_seeds"] and run_type in ("staged_learning", "sampling") \
            and not smiles_file:
        raise UsageError(
            f"The {generator_key} generator needs seed structures.",
            hint=f"Pass --smiles-file with {spec.get('seed_format', 'seed structures')}. "
                 "Use 'prepare-seeds' to validate and de-duplicate them first.",
        )
    if run_type == "transfer_learning" and not smiles_file:
        raise UsageError(
            "Transfer learning needs a set of molecules to learn from.",
            hint="Pass --smiles-file with the known actives.",
        )
    if run_type == "scoring" and not smiles_file:
        raise UsageError(
            "Scoring needs molecules to score.",
            hint="Pass --smiles-file.",
        )

    scoring: dict[str, Any] = {}
    if RUN_MODES[run_type]["needs_scoring"]:
        scoring = build_scoring_section(scoring_profile, components, aggregation)
    elif scoring_profile or components:
        warnings.append(
            f"A scoring profile was supplied but run mode {run_type!r} does "
            "not use one; it has been ignored."
        )

    resolved_device = device or "cuda:0"
    if resolved_device.startswith("cuda"):
        try:
            import torch
            if not torch.cuda.is_available():
                warnings.append(
                    f"device is {resolved_device!r} but no CUDA device is "
                    "visible here. If you will submit this to a GPU node that "
                    "is fine; to run locally, pass --device cpu."
                )
        except ImportError:
            warnings.append(
                f"device is {resolved_device!r} but PyTorch is not installed "
                "in this interpreter, so the setting could not be checked."
            )

    settings: dict[str, Any] = {
        "run_type": run_type,
        "generator": generator_key,
        "generator_architecture": spec["architecture"],
        "prior_file": prior_reference,
        "device": resolved_device,
    }
    expected_outputs: list[str] = []

    header = [
        f"# REINVENT 4 configuration - {RUN_MODES[run_type]['description']}",
        f"# Generator: {generator_key} ({spec['architecture']}) - {spec['use_for']}",
        "# Generated by agskills. Review before submitting a long run.",
        "",
        f"run_type = {toml_value(run_type)}",
        f"device = {toml_value(resolved_device)}",
        f'tb_logdir = "tb_logs"',
        f'json_out_config = "_{run_type}.json"',
        "",
    ]
    lines = list(header)

    if run_type == "staged_learning":
        steps = num_steps if num_steps is not None else 300
        resolved_batch = batch_size if batch_size is not None else 64
        prefix = csv_prefix or "staged_learning"
        checkpoint = checkpoint_file or f"{prefix}.chkpt"
        lines += [
            "[parameters]",
            "",
            f"summary_csv_prefix = {toml_value(prefix)}",
            "use_checkpoint = false",
            "purge_memories = false",
            "",
            f"prior_file = {toml_value(prior_reference)}",
            "# agent_file starts equal to the prior; replace it with a",
            "# checkpoint to resume an interrupted run.",
            f"agent_file = {toml_value(prior_reference)}",
        ]
        if smiles_file:
            lines.append(f"smiles_file = {toml_value(smiles_file)}")
        if sample_strategy and generator_key in ("mol2mol", "pepinvent"):
            lines.append(f"sample_strategy = {toml_value(sample_strategy)}")
            lines.append("distance_threshold = 100")
        lines += [
            f"batch_size = {toml_value(resolved_batch)}",
            "unique_sequences = true",
            "randomize_smiles = true",
            "",
            "[learning_strategy]",
            'type = "dap"',
            f"sigma = {toml_value(sigma if sigma is not None else 128)}",
            f"rate = {toml_value(learning_rate if learning_rate is not None else 0.0001)}",
            "",
        ]
        if diversity_filter:
            lines += [
                "# The diversity filter stops the agent collapsing onto one",
                "# scaffold, which is the most common failure of an",
                "# unconstrained RL run.",
                "[diversity_filter]",
                'type = "IdenticalMurckoScaffold"',
                "bucket_size = 25",
                "minscore = 0.4",
                "",
            ]
        if inception_smiles_file and generator_key == "reinvent":
            lines += [
                "[inception]",
                f"smiles_file = {toml_value(inception_smiles_file)}",
                "memory_size = 100",
                "sample_size = 10",
                "",
            ]
        lines += [
            "### Stage 1. Stages are a list, hence the double brackets.",
            "[[stage]]",
            "",
            f"chkpt_file = {toml_value(checkpoint)}",
            'termination = "simple"',
            f"max_score = {toml_value(max_score if max_score is not None else 0.7)}",
            f"min_steps = {toml_value(min_steps if min_steps is not None else max(10, steps // 10))}",
            f"max_steps = {toml_value(max_steps if max_steps is not None else steps)}",
        ]
        _emit_scoring(lines, scoring, "stage.scoring")
        settings.update({
            "max_steps": max_steps if max_steps is not None else steps,
            "batch_size": resolved_batch,
            "sigma": sigma if sigma is not None else 128,
            "summary_csv_prefix": prefix,
            "diversity_filter": diversity_filter,
        })
        expected_outputs = [f"{prefix}_1.csv", checkpoint, "tb_logs/"]
        molecules = (max_steps or steps) * resolved_batch
        next_steps = [
            f"This run will evaluate roughly {molecules:,} molecules "
            f"({max_steps or steps} steps x batch {resolved_batch}).",
            "Run it: ag-compound-synthesis run --config <this file> "
            "--output run_report.json",
            f"Then analyse: ag-compound-synthesis analyze-results --csv "
            f"{prefix}_1.csv --top-n 25 --output top_compounds.json",
        ]
        if (max_steps or steps) > 100 and resolved_device.startswith("cpu"):
            warnings.append(
                f"{max_steps or steps} RL steps on CPU will take many hours. "
                "Generate an HPC script with 'generate-hpc-script' and submit "
                "it to a GPU node."
            )

    elif run_type == "transfer_learning":
        epochs = num_epochs if num_epochs is not None else 10
        resolved_batch = batch_size if batch_size is not None else 50
        out_model = output_model or f"TL_{generator_key}.model"
        lines += [
            "[parameters]",
            "",
            f"num_epochs = {toml_value(epochs)}",
            f"save_every_n_epochs = {toml_value(max(1, epochs // 5))}",
            f"batch_size = {toml_value(resolved_batch)}",
            "# num_refs computes similarity against reference molecules.",
            "# Set it to 0 for datasets above ~200 molecules, or the run",
            "# becomes very slow.",
            "num_refs = 0",
            "sample_batch_size = 100",
            "",
            f"input_model_file = {toml_value(prior_reference)}",
            f"smiles_file = {toml_value(smiles_file)}",
            f"output_model_file = {toml_value(out_model)}",
        ]
        if validation_smiles_file:
            lines.append(
                f"validation_smiles_file = {toml_value(validation_smiles_file)}"
            )
        else:
            warnings.append(
                "No validation set was given, so there is no way to detect "
                "over-fitting to the training molecules. Pass "
                "--validation-smiles-file with a held-out subset."
            )
        if generator_key == "mol2mol":
            lines += [
                "",
                "# Mol2Mol learns from pairs of molecules; these thresholds",
                "# decide which pairs are similar enough to train on.",
                'pairs.type = "tanimoto"',
                "pairs.upper_threshold = 1.0",
                "pairs.lower_threshold = 0.7",
            ]
        lines.append("")
        settings.update({"num_epochs": epochs, "batch_size": resolved_batch,
                         "output_model_file": out_model})
        expected_outputs = [out_model, "tb_logs/"]
        next_steps = [
            "Run it: ag-compound-synthesis run --config <this file> "
            "--output tl_report.json",
            f"Then use {out_model} as the --prior of a staged_learning run "
            "to optimise within the focused chemical space.",
        ]

    elif run_type == "sampling":
        count = num_smiles if num_smiles is not None else 1000
        out_csv = output_csv or "sampled.csv"
        model = model_file or prior_reference
        lines += [
            "[parameters]",
            "",
            f"model_file = {toml_value(model)}",
            f"output_file = {toml_value(out_csv)}",
            f"num_smiles = {toml_value(count)}",
            "unique_molecules = true",
            "randomize_smiles = true",
        ]
        if smiles_file:
            lines.append(f"smiles_file = {toml_value(smiles_file)}")
        if sample_strategy and generator_key in ("mol2mol", "pepinvent"):
            lines.append(f"sample_strategy = {toml_value(sample_strategy)}")
        lines.append("")
        settings.update({"model_file": model, "num_smiles": count,
                         "output_file": out_csv})
        expected_outputs = [out_csv]
        next_steps = [
            f"Run it, then screen {out_csv} with: ag-compound-screening "
            f"admet-filter --csv {out_csv} --smiles-col SMILES "
            "--strictness relaxed --output triage.json",
        ]
        if spec["needs_seeds"] and not smiles_file:
            warnings.append(
                f"{generator_key} normally samples around input structures; "
                "without --smiles-file it will sample unconditionally."
            )

    else:  # scoring
        out_csv = output_csv or "scoring.csv"
        lines += [
            "[parameters]",
            "",
            f"smiles_file = {toml_value(smiles_file)}",
            f"output_csv = {toml_value(out_csv)}",
        ]
        _emit_scoring(lines, scoring, "scoring")
        settings.update({"smiles_file": smiles_file, "output_csv": out_csv})
        expected_outputs = [out_csv]
        next_steps = [
            f"Run it, then rank the results: ag-compound-synthesis "
            f"analyze-results --csv {out_csv} --output ranked.json",
        ]

    if scoring:
        settings["scoring_components"] = scoring["n_components"]
        settings["aggregation"] = scoring["type"]
        if scoring["metadata"].get("caveat"):
            warnings.append(scoring["metadata"]["caveat"])

    return ConfigResult(
        toml="\n".join(lines).rstrip() + "\n",
        run_type=run_type,
        generator=generator_key,
        settings=settings,
        expected_outputs=expected_outputs,
        warnings=warnings,
        next_steps=next_steps,
        scoring_metadata=scoring.get("metadata", {}),
    )


def render_toml(**kwargs) -> str:
    """Convenience wrapper returning only the TOML text."""
    return build_config(**kwargs).toml
