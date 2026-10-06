"""SLURM submission-script generation, shared by every skill.

Four skills generated SLURM scripts independently, with four slightly
different sets of ``#SBATCH`` directives and four separate bugs. They all
come through here now.

Generated scripts are deliberately defensive, because a silent failure at
hour three of a 24-hour allocation is expensive:

* ``set -euo pipefail`` so a failed step aborts rather than continuing.
* Inputs are checked for existence *before* the long-running step starts.
* The environment activation is checked, rather than assumed.
* ``--time``, ``--mem`` and GPU directives are validated here, so a typo is
  caught at generation time instead of at submission.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

from .errors import UsageError

__all__ = [
    "SlurmResources",
    "JOB_TYPES",
    "generate_slurm_script",
    "validate_walltime",
    "validate_memory",
]

_WALLTIME_RE = re.compile(
    r"^(?:(?P<days>\d+)-)?(?P<h>\d{1,3}):(?P<m>\d{2})(?::(?P<s>\d{2}))?$"
)
_MEMORY_RE = re.compile(r"^(?P<amount>\d+)(?P<unit>[KMGT])B?$", re.IGNORECASE)

#: Job types this module can generate, with a one-line description.
JOB_TYPES: dict[str, str] = {
    "gromacs-md": "GROMACS production molecular dynamics (GPU).",
    "reinvent": "REINVENT 4 generative chemistry run (GPU).",
    "vina-screen": "AutoDock Vina batch docking as a SLURM job array (CPU).",
    "colabfold": "ColabFold structure prediction (GPU).",
    "alphafold2": "AlphaFold2 structure prediction (GPU).",
}


def validate_walltime(value: str) -> str:
    """Validate a SLURM walltime.

    Accepts ``HH:MM``, ``HH:MM:SS`` and ``D-HH:MM:SS``.

    >>> validate_walltime("24:00:00")
    '24:00:00'
    >>> validate_walltime("2-12:00:00")
    '2-12:00:00'
    >>> validate_walltime("24h")
    Traceback (most recent call last):
    agskills.errors.UsageError: ...
    """
    match = _WALLTIME_RE.match((value or "").strip())
    if not match:
        raise UsageError(
            f"{value!r} is not a SLURM walltime.",
            hint="Use HH:MM:SS (for example 24:00:00) or D-HH:MM:SS.",
        )
    if int(match.group("m")) >= 60 or (match.group("s") and int(match.group("s")) >= 60):
        raise UsageError(f"{value!r} has a minutes or seconds field above 59.")
    return value.strip()


def validate_memory(value: str) -> str:
    """Validate a SLURM memory specification such as ``64G``.

    >>> validate_memory("128G")
    '128G'
    >>> validate_memory("lots")
    Traceback (most recent call last):
    agskills.errors.UsageError: ...
    """
    match = _MEMORY_RE.match((value or "").strip())
    if not match:
        raise UsageError(
            f"{value!r} is not a memory size.",
            hint="Use a number and a unit, for example 64G or 128000M.",
        )
    return value.strip().upper().removesuffix("B")


@dataclass
class SlurmResources:
    """Resource request for one SLURM job."""

    job_name: str = "agskills"
    partition: str | None = None
    account: str | None = None
    gpu: str | None = "a100"
    ngpu: int = 1
    ncpus: int = 8
    mem: str = "64G"
    time: str = "24:00:00"
    array: str | None = None
    conda_env: str | None = None
    module_loads: list[str] = field(default_factory=list)
    email: str | None = None
    extra_sbatch: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.time = validate_walltime(self.time)
        self.mem = validate_memory(self.mem)
        if self.ngpu < 0:
            raise UsageError("--ngpu cannot be negative")
        if self.ncpus < 1:
            raise UsageError("--ncpus must be at least 1")
        if self.array is not None and not re.match(
            r"^\d+(-\d+)?(:\d+)?(%\d+)?(,\d+(-\d+)?)*$", self.array.strip()
        ):
            raise UsageError(
                f"{self.array!r} is not a SLURM array specification.",
                hint="Use forms like 1-100, 1-100%10 or 1,3,5.",
            )

    def directives(self) -> list[str]:
        """Render the ``#SBATCH`` block."""
        lines = [
            f"#SBATCH --job-name={self.job_name}",
            f"#SBATCH --output={self.job_name}-%j.out",
            f"#SBATCH --error={self.job_name}-%j.err",
            f"#SBATCH --time={self.time}",
            f"#SBATCH --mem={self.mem}",
            f"#SBATCH --cpus-per-task={self.ncpus}",
        ]
        if self.array:
            # An array job writes one file per task, so the names must
            # include the task id or every task overwrites the same file.
            lines[1] = f"#SBATCH --output={self.job_name}-%A_%a.out"
            lines[2] = f"#SBATCH --error={self.job_name}-%A_%a.err"
            lines.append(f"#SBATCH --array={self.array}")
        if self.ngpu > 0 and self.gpu:
            lines.append(f"#SBATCH --gres=gpu:{self.gpu}:{self.ngpu}")
        elif self.ngpu > 0:
            lines.append(f"#SBATCH --gres=gpu:{self.ngpu}")
        if self.partition:
            lines.append(f"#SBATCH --partition={self.partition}")
        if self.account:
            lines.append(f"#SBATCH --account={self.account}")
        if self.email:
            lines.append(f"#SBATCH --mail-user={self.email}")
            lines.append("#SBATCH --mail-type=END,FAIL")
        lines.extend(self.extra_sbatch)
        return lines

    def preamble(self) -> list[str]:
        """Render module loads and environment activation, with checks."""
        lines: list[str] = []
        for module in self.module_loads:
            lines.append(f"module load {module}")
        if self.conda_env:
            lines += [
                "",
                '# Activate the environment. "conda activate" needs the shell',
                "# hook, which is not loaded in a non-interactive SLURM shell.",
                'if command -v conda >/dev/null 2>&1; then',
                '  eval "$(conda shell.bash hook)"',
                f'  conda activate {self.conda_env}',
                "else",
                f'  echo "ERROR: conda not found; cannot activate '
                f'{self.conda_env}" >&2',
                "  exit 1",
                "fi",
            ]
        return lines


def _require_file_block(label: str, variable: str) -> list[str]:
    return [
        f'if [ ! -f "{variable}" ]; then',
        f'  echo "ERROR: {label} not found: {variable}" >&2',
        "  exit 1",
        "fi",
    ]


def _header(job_type: str, resources: SlurmResources,
            settings: dict[str, Any]) -> list[str | None]:
    described = ", ".join(f"{k}={v}" for k, v in sorted(settings.items()) if v
                          not in (None, ""))
    return [
        "#!/bin/bash",
        *resources.directives(),
        "",
        f"# {JOB_TYPES.get(job_type, job_type)}",
        "# Generated by agskills. Review the directives above against your",
        "# cluster's partitions and GPU names before submitting.",
        f"# Settings: {described}" if described else None,
        "",
        "set -euo pipefail",
        "",
        'echo "Job $SLURM_JOB_ID starting on $(hostname) at $(date -Is)"',
        *resources.preamble(),
        "",
    ]


def _gromacs_md(resources: SlurmResources, settings: dict[str, Any]) -> list[str]:
    tpr = settings.get("tpr", "md.tpr")
    deffnm = settings.get("deffnm", "md")
    ns = settings.get("production_ns", 100)
    return [
        f'TPR="{tpr}"',
        f'DEFFNM="{deffnm}"',
        "",
        *_require_file_block("GROMACS run input", "$TPR"),
        "",
        'if ! command -v gmx >/dev/null 2>&1 && ! command -v gmx_mpi >/dev/null 2>&1; then',
        '  echo "ERROR: neither gmx nor gmx_mpi is on PATH" >&2',
        "  exit 1",
        "fi",
        'GMX=$(command -v gmx_mpi || command -v gmx)',
        "",
        f'echo "Running {ns} ns of production MD"',
        "",
        "# -maxh matches the walltime so GROMACS writes a final checkpoint",
        "# instead of being killed mid-write. -cpi resumes automatically if a",
        "# checkpoint from a previous attempt is present.",
        f'MAXH={settings.get("maxh", "23.5")}',
        'RESUME=""',
        'if [ -f "${DEFFNM}.cpt" ]; then',
        '  echo "Found ${DEFFNM}.cpt - resuming"',
        '  RESUME="-cpi ${DEFFNM}.cpt"',
        "fi",
        "",
        '"$GMX" mdrun \\',
        '  -s "$TPR" \\',
        '  -deffnm "$DEFFNM" \\',
        ('  -nb gpu -pme gpu -bonded gpu \\' if resources.ngpu > 0 else None),
        f'  -ntomp {resources.ncpus} \\',
        '  -maxh "$MAXH" \\',
        '  $RESUME',
        "",
        'echo "Production MD finished. Outputs: ${DEFFNM}.xtc ${DEFFNM}.edr ${DEFFNM}.gro"',
    ]


def _reinvent(resources: SlurmResources, settings: dict[str, Any]) -> list[str]:
    config = settings.get("config", "reinvent.toml")
    log = settings.get("log", "reinvent.log")
    return [
        f'CONFIG="{config}"',
        f'LOG="{log}"',
        "",
        *_require_file_block("REINVENT config", "$CONFIG"),
        "",
        'if ! command -v reinvent >/dev/null 2>&1; then',
        '  echo "ERROR: the reinvent executable is not on PATH." >&2',
        '  echo "Install it with: pip install -e /path/to/REINVENT4" >&2',
        "  exit 1",
        "fi",
        "",
        'nvidia-smi || echo "WARNING: nvidia-smi unavailable; is a GPU allocated?"',
        "",
        'reinvent --log-filename "$LOG" "$CONFIG"',
        "",
        'echo "REINVENT finished. Review $LOG and the CSV outputs it names."',
    ]


def _vina_screen(resources: SlurmResources, settings: dict[str, Any]) -> list[str]:
    receptor = settings.get("receptor", "receptor.pdbqt")
    ligand_dir = settings.get("ligand_dir", "ligands")
    site = settings.get("site_config", "binding_site.json")
    exhaustiveness = settings.get("exhaustiveness", 32)
    return [
        f'RECEPTOR="{receptor}"',
        f'LIGAND_DIR="{ligand_dir}"',
        f'SITE="{site}"',
        'OUT_DIR="docked"',
        "",
        *_require_file_block("receptor PDBQT", "$RECEPTOR"),
        *_require_file_block("binding site JSON", "$SITE"),
        'if [ ! -d "$LIGAND_DIR" ]; then',
        '  echo "ERROR: ligand directory not found: $LIGAND_DIR" >&2',
        "  exit 1",
        "fi",
        'mkdir -p "$OUT_DIR"',
        "",
        "# Read the box from the JSON written by 'target-preparation",
        "# define-site'. python3 is used rather than jq, which many clusters",
        "# do not install.",
        'read -r CX CY CZ SX SY SZ < <(python3 -c "',
        "import json,sys",
        "s=json.load(open(sys.argv[1]))",
        "print(s['center_x'],s['center_y'],s['center_z'],s['size_x'],s['size_y'],s['size_z'])",
        '" "$SITE")',
        'echo "Box centre ($CX $CY $CZ) size ($SX $SY $SZ)"',
        "",
        "# Each array task takes every Nth ligand, so the set is partitioned",
        "# with no overlap and no missed files.",
        'mapfile -t LIGANDS < <(find "$LIGAND_DIR" -name "*.pdbqt" | sort)',
        'TOTAL=${#LIGANDS[@]}',
        'if [ "$TOTAL" -eq 0 ]; then',
        '  echo "ERROR: no .pdbqt ligands in $LIGAND_DIR" >&2',
        "  exit 1",
        "fi",
        'TASK_ID=${SLURM_ARRAY_TASK_ID:-1}',
        'TASK_COUNT=${SLURM_ARRAY_TASK_COUNT:-1}',
        'echo "Task $TASK_ID of $TASK_COUNT over $TOTAL ligands"',
        "",
        'for (( i=TASK_ID-1; i<TOTAL; i+=TASK_COUNT )); do',
        '  LIG="${LIGANDS[$i]}"',
        '  BASE=$(basename "$LIG" .pdbqt)',
        '  if [ -s "$OUT_DIR/${BASE}_out.pdbqt" ]; then',
        '    echo "  skip $BASE (already docked)"',
        "    continue",
        "  fi",
        '  vina --receptor "$RECEPTOR" --ligand "$LIG" \\',
        '    --center_x "$CX" --center_y "$CY" --center_z "$CZ" \\',
        '    --size_x "$SX" --size_y "$SY" --size_z "$SZ" \\',
        f'    --exhaustiveness {exhaustiveness} --cpu {resources.ncpus} \\',
        '    --out "$OUT_DIR/${BASE}_out.pdbqt" \\',
        '    > "$OUT_DIR/${BASE}.log" 2>&1 || echo "  FAILED: $BASE" >&2',
        "done",
        "",
        'echo "Array task $TASK_ID complete."',
        'echo "Collect scores with: grep -H \'^   1\' docked/*.log"',
    ]


def _folding(job_type: str, resources: SlurmResources,
             settings: dict[str, Any]) -> list[str]:
    fasta = settings.get("fasta", "target.fasta")
    db = settings.get("db_path", "/data/colabfold_dbs")
    out = settings.get("out_dir", "predictions")
    if job_type == "colabfold":
        body = [
            'if ! command -v colabfold_batch >/dev/null 2>&1; then',
            '  echo "ERROR: colabfold_batch is not on PATH" >&2',
            "  exit 1",
            "fi",
            "",
            'colabfold_batch \\',
            '  --num-recycle 3 \\',
            '  --model-type auto \\',
            '  --amber --use-gpu-relax \\',
            '  "$FASTA" "$OUT_DIR"',
        ]
    else:
        body = [
            f'DB="{db}"',
            'if [ ! -d "$DB" ]; then',
            '  echo "ERROR: AlphaFold database directory not found: $DB" >&2',
            "  exit 1",
            "fi",
            "",
            "python3 /opt/alphafold/run_alphafold.py \\",
            '  --fasta_paths="$FASTA" \\',
            '  --output_dir="$OUT_DIR" \\',
            '  --data_dir="$DB" \\',
            "  --model_preset=monomer \\",
            "  --db_preset=full_dbs \\",
            "  --max_template_date=2024-01-01 \\",
            "  --use_gpu_relax=true",
        ]
    return [
        f'FASTA="{fasta}"',
        f'OUT_DIR="{out}"',
        "",
        *_require_file_block("input FASTA", "$FASTA"),
        'mkdir -p "$OUT_DIR"',
        "",
        *body,
        "",
        'echo "Prediction written to $OUT_DIR"',
        "echo 'Assess it with:'",
        "echo '  ag-target-preparation assess-structure --pdb "
        "<model>.pdb --output quality.json'",
    ]


_BUILDERS = {
    "gromacs-md": _gromacs_md,
    "reinvent": _reinvent,
    "vina-screen": _vina_screen,
}


def generate_slurm_script(job_type: str, resources: SlurmResources,
                          settings: dict[str, Any] | None = None) -> str:
    """Render a SLURM submission script.

    Args:
        job_type: One of :data:`JOB_TYPES`.
        resources: The resource request.
        settings: Job-specific paths and parameters.

    Raises:
        UsageError: If *job_type* is unknown, or a job type's required
            settings are missing.
    """
    if job_type not in JOB_TYPES:
        raise UsageError(
            f"Unknown job type {job_type!r}.",
            hint=f"Choose one of: {', '.join(sorted(JOB_TYPES))}",
        )
    options = dict(settings or {})

    if job_type == "vina-screen" and not resources.array:
        raise UsageError(
            "A vina-screen job needs --array to partition the ligand set.",
            hint="For example --array 1-100 to spread the ligands over 100 "
                 "tasks.",
        )
    if job_type == "reinvent" and not options.get("config"):
        raise UsageError("A reinvent job needs --config pointing at its TOML file.")

    if job_type in ("colabfold", "alphafold2"):
        body = _folding(job_type, resources, options)
    else:
        body = _BUILDERS[job_type](resources, options)

    lines: list[str | None] = [
        *_header(job_type, resources, options),
        *body,
        "",
        'echo "Job $SLURM_JOB_ID finished at $(date -Is)"',
    ]
    # None marks a line the job type chose to omit; "" is a real blank line.
    return "\n".join(line for line in lines if line is not None) + "\n"
