"""Submission instructions generated alongside each SLURM script.

A generated script is only half of what someone needs. The other half is
knowing which tools the cluster must provide, how to get them there, which
input files have to be in place, and what to type to submit and watch the
job. That is what this module writes: one Markdown file per script, named
after it.

Everything here is job-type specific and deliberately concrete. "Install the
dependencies" is not an instruction; ``conda install -c conda-forge gromacs``
is. Where a tool has no single install route - AlphaFold2 needs databases
measured in terabytes - the document says so rather than pretending.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .hpc import CPU_ONLY_JOB_TYPES, JOB_TYPES, SlurmResources, maxh_for

__all__ = ["ToolRequirement", "JobGuide", "JOB_GUIDES",
           "submission_instructions", "instructions_path"]


@dataclass(frozen=True)
class ToolRequirement:
    """One tool the cluster must provide, and how to get it there."""

    name: str
    #: The command that proves it is available.
    probe: str
    #: Most clusters provide scientific software as modules. Checked first
    #: because it needs no installation and is tuned for the hardware.
    module_hint: str | None = None
    #: A user-space install, for when no module exists.
    install: tuple[str, ...] = ()
    #: Anything that is not a package: databases, model weights, licences.
    extra: tuple[str, ...] = ()
    note: str | None = None


@dataclass(frozen=True)
class JobGuide:
    """Everything needed to document one job type."""

    summary: str
    tools: tuple[ToolRequirement, ...]
    #: Files that must sit beside the script, as (path-or-pattern, why).
    inputs: tuple[tuple[str, str], ...] = ()
    outputs: tuple[str, ...] = ()
    #: Job-type specific warnings worth reading before submitting.
    cautions: tuple[str, ...] = ()
    resumable: bool = False
    produced_by: str | None = None


GROMACS = ToolRequirement(
    name="GROMACS",
    probe="gmx --version   # or gmx_mpi --version",
    module_hint="module avail gromacs   # then e.g. module load gromacs/2024.2",
    install=(
        "# User-space install, if no module exists:",
        "conda install -c conda-forge gromacs",
    ),
    note=("Prefer the cluster's module. A site build is compiled against the "
          "local MPI, CUDA and hardware, and is usually several times faster "
          "than a generic conda build. Use the GPU-enabled build if you "
          "requested a GPU."),
)

VINA = ToolRequirement(
    name="AutoDock Vina",
    probe="vina --version",
    module_hint="module avail vina autodock",
    install=(
        "conda install -c conda-forge vina",
        "# or download a static binary:",
        "#   https://github.com/ccsb-scripps/AutoDock-Vina/releases",
    ),
    note=("Vina is CPU-only, so this job requests no GPU. The static binary "
          "release needs no compiler and is the simplest route on a cluster "
          "without conda."),
)

OPENBABEL = ToolRequirement(
    name="Open Babel or Meeko",
    probe="obabel -V   # or: python -c 'import meeko'",
    install=("conda install -c conda-forge openbabel",
             "pip install meeko"),
    note=("Only needed if the ligand files are not already PDBQT. Prepare "
          "them locally with the compound-screening skill and copy them "
          "over, and neither is required here."),
)

REINVENT = ToolRequirement(
    name="REINVENT 4",
    probe="reinvent --help",
    module_hint="module avail reinvent   # rarely provided; expect to install",
    install=(
        "git clone https://github.com/MolecularAI/REINVENT4",
        "pip install -e REINVENT4",
        "export REINVENT_DIR=$PWD/REINVENT4",
    ),
    extra=(
        "Prior model files, from https://doi.org/10.5281/zenodo.15641296",
        "Unpack them into $REINVENT_DIR/priors, or point "
        "REINVENT_PRIOR_BASE at wherever they live.",
    ),
    note="Do this inside a conda environment, not the cluster's system Python.",
)

PYTORCH = ToolRequirement(
    name="PyTorch",
    probe="python -c 'import torch; print(torch.__version__, "
          "torch.cuda.is_available())'",
    module_hint="module avail cuda   # load CUDA before installing PyTorch",
    install=(
        "# Pick the build matching the cluster's CUDA version:",
        "#   https://pytorch.org/get-started/locally/",
        "# e.g. for CUDA 12.1:",
        "pip install torch --index-url https://download.pytorch.org/whl/cu121",
    ),
    note=("The probe must print True for CUDA. A CPU-only PyTorch will run "
          "but take orders of magnitude longer, and the job requested a GPU."),
)

COLABFOLD = ToolRequirement(
    name="ColabFold",
    probe="colabfold_batch --help",
    module_hint="module avail colabfold alphafold",
    install=(
        "# localcolabfold installs colabfold_batch and its weights:",
        "#   https://github.com/YoshitakaMo/localcolabfold",
        "git clone https://github.com/YoshitakaMo/localcolabfold",
        "bash localcolabfold/install_colabbatch_linux.sh",
        'export PATH="$PWD/localcolabfold/colabfold-conda/bin:$PATH"',
    ),
    extra=(
        "Model weights download on first run, so the first job needs "
        "outbound network access from the compute node. If the nodes are "
        "offline, run one job on a login node first to populate the cache.",
    ),
    note=("ColabFold uses MMseqs2 against a remote server by default. On a "
          "cluster without outbound access you need a local database and "
          "--mmseqs flags, which is a much larger undertaking."),
)

ALPHAFOLD = ToolRequirement(
    name="AlphaFold2",
    probe="ls /opt/alphafold/run_alphafold.py",
    module_hint="module avail alphafold   # most sites provide this centrally",
    install=(
        "# From source, in a container:",
        "#   https://github.com/google-deepmind/alphafold",
        "# Most clusters provide AlphaFold as a module or Singularity image.",
        "# Ask your HPC support before building it yourself.",
    ),
    extra=(
        "The genetic databases are several terabytes. Do not download them "
        "into your home directory.",
        "Check for a shared copy first: ask support, or look under "
        "/scratch, /data or /shared for an alphafold directory.",
    ),
    note=("Of the five job types this is the one most worth asking your HPC "
          "support about rather than installing. If you only need a "
          "structure and not full AlphaFold2 accuracy, ColabFold is far "
          "easier to set up."),
)

CONDA = ToolRequirement(
    name="conda or mamba",
    probe="conda --version",
    module_hint="module avail anaconda miniconda mamba",
    install=(
        "# If the cluster has none, install Miniforge into your own space:",
        "curl -L -O https://github.com/conda-forge/miniforge/releases/"
        "latest/download/Miniforge3-Linux-x86_64.sh",
        "bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3",
        '$HOME/miniforge3/bin/conda init bash',
    ),
    note="Only needed if you are installing tools yourself rather than "
         "using modules.",
)


JOB_GUIDES: dict[str, JobGuide] = {
    "gromacs-md": JobGuide(
        summary="Runs production molecular dynamics from a prepared .tpr.",
        tools=(GROMACS,),
        inputs=(
            ("md.tpr", "the run input, produced by `gmx grompp`"),
        ),
        outputs=("md.xtc  - the trajectory",
                 "md.edr  - energies",
                 "md.gro  - the final frame",
                 "md.cpt  - the checkpoint, used to resume"),
        cautions=(
            "Run the preparation steps *before* submitting this. The "
            "md-simulation skill writes a `setup_md.sh` that does topology, "
            "box, solvation, ions, energy minimisation and equilibration, "
            "and ends by producing the .tpr this job needs.",
            "A production run rarely fits one allocation. This script stops "
            "cleanly before the walltime and resumes from the checkpoint, so "
            "resubmitting it continues rather than restarting.",
        ),
        resumable=True,
        produced_by="ag-md-simulation gmxapi-setup, then bash setup_md.sh",
    ),
    "reinvent": JobGuide(
        summary="Runs a REINVENT 4 generative chemistry job from a TOML "
                "config.",
        tools=(CONDA, PYTORCH, REINVENT),
        inputs=(
            ("reinvent.toml", "the run configuration"),
            ("the prior model named in the config", "e.g. priors/reinvent.prior"),
            ("any seed .smi file the config references", "for libinvent, "
             "linkinvent, mol2mol or transfer learning"),
        ),
        outputs=("<prefix>_1.csv  - sampled molecules with scores",
                 "<prefix>.chkpt  - the agent checkpoint",
                 "tb_logs/        - TensorBoard logs"),
        cautions=(
            "Run `ag-compound-synthesis preflight --config reinvent.toml` "
            "before submitting. It lists every file the config references "
            "that does not exist, which otherwise fails after the queue "
            "wait rather than before it.",
            "Paths inside the config are relative to the config file. Copy "
            "the whole directory to the cluster, not the config alone.",
        ),
        produced_by="ag-compound-synthesis generate-config",
    ),
    "vina-screen": JobGuide(
        summary="Docks a ligand library with AutoDock Vina as a SLURM job "
                "array, each task taking an even share of the ligands.",
        tools=(VINA, OPENBABEL),
        inputs=(
            ("receptor.pdbqt", "the prepared receptor"),
            ("binding_site.json", "the docking box; read at job time"),
            ("ligands/*.pdbqt", "one file per ligand"),
        ),
        outputs=("one log and one pose file per ligand",
                 "collect the scores afterwards with "
                 "`ag-compound-screening dock --parse-logs`"),
        cautions=(
            "This is a **CPU** job array and requests no GPU. Vina cannot "
            "use one.",
            "Size the array to the library: `--array 1-100` spreads the "
            "ligands over 100 tasks. Check your cluster's limit on "
            "concurrent array tasks (`scontrol show config | grep MaxArray`) "
            "and consider `--array 1-100%20` to cap how many run at once.",
            "Confirm the box is on the pocket you intend, not on a cofactor, "
            "before committing a large screen. The target-preparation skill "
            "flags cofactors for this reason.",
        ),
        produced_by="ag-target-preparation prepare-receptor (receptor and "
                    "site), then ag-compound-screening for the ligands",
    ),
    "colabfold": JobGuide(
        summary="Predicts structures for every sequence in a FASTA file "
                "using ColabFold.",
        tools=(COLABFOLD,),
        inputs=(("target.fasta", "one or more sequences"),),
        outputs=("<name>_relaxed_rank_001_*.pdb  - the best model",
                 "<name>_scores_rank_*.json      - pLDDT and PAE",
                 "<name>_coverage.png            - MSA depth"),
        cautions=(
            "Check pLDDT **in the binding site**, not just the global mean. "
            "`ag-target-preparation assess-structure` reports the "
            "low-confidence regions.",
            "The B-factor column holds pLDDT, where higher is better. Do not "
            "read it as a crystallographic B-factor.",
        ),
        produced_by="ag-target-preparation generate-af2-script",
    ),
    "alphafold2": JobGuide(
        summary="Predicts structures with the full AlphaFold2 pipeline, "
                "including its genetic databases.",
        tools=(ALPHAFOLD,),
        inputs=(("target.fasta", "one or more sequences"),
                ("the database directory", "set by --db in the script")),
        outputs=("ranked_0.pdb       - the best model",
                 "ranking_debug.json - pLDDT ranking",
                 "features.pkl       - MSA and template features"),
        cautions=(
            "The databases are several terabytes. Find the shared copy on "
            "your cluster rather than downloading your own.",
            "An MSA-based run takes hours per sequence, most of it on CPU "
            "building alignments rather than on the GPU.",
        ),
        produced_by="ag-target-preparation generate-af2-script "
                    "--method alphafold2",
    ),
}


def instructions_path(script_path: str | Path) -> Path:
    """Where the instructions for *script_path* belong.

    ``submit.sh`` gets ``submit_SUBMIT.md`` beside it, so the pair stays
    together when the directory is copied to a cluster.
    """
    path = Path(script_path)
    return path.with_name(f"{path.stem}_SUBMIT.md")


def _render_tool(tool: ToolRequirement, index: int) -> list[str]:
    lines = [f"### {index}. {tool.name}", ""]
    lines += ["Check whether it is already there:", "",
              "```bash", tool.probe, "```", ""]
    if tool.module_hint:
        lines += ["If that fails, look for a module - always prefer one over "
                  "installing it yourself:", "",
                  "```bash", tool.module_hint, "```", ""]
    if tool.install:
        lines += ["Otherwise install it:", "", "```bash"]
        lines += list(tool.install)
        lines += ["```", ""]
    if tool.extra:
        lines += ["You also need:", ""]
        lines += [f"-   {item}" for item in tool.extra]
        lines += [""]
    if tool.note:
        lines += [f"> {tool.note}", ""]
    return lines


def submission_instructions(
    job_type: str,
    resources: SlurmResources,
    script_path: str | Path,
    settings: dict[str, Any] | None = None,
) -> str:
    """Render the Markdown instructions for one generated script."""
    guide = JOB_GUIDES.get(job_type)
    options = dict(settings or {})
    script = Path(script_path).name
    cpu_only = job_type in CPU_ONLY_JOB_TYPES

    out: list[str] = [
        f"# How to run `{script}`",
        "",
        guide.summary if guide else JOB_TYPES.get(job_type, job_type),
        "",
        "Generated by agskills alongside the script. Work through it in "
        "order: install, stage the inputs, adjust the directives, submit.",
        "",
        "--------------------------------------------------------------",
        "",
        "## 1. What this job asks SLURM for",
        "",
        "| Resource | Value |",
        "|---|---|",
        f"| Walltime | `{resources.time}` |",
        f"| Memory | `{resources.mem}` |",
        f"| CPUs per task | `{resources.ncpus}` |",
    ]
    if cpu_only or resources.ngpu <= 0:
        out.append("| GPUs | none - this job is CPU-only |")
    else:
        out.append(f"| GPUs | `{resources.ngpu} x {resources.gpu or 'gpu'}` |")
    if resources.array:
        out.append(f"| Array | `{resources.array}` |")
    if resources.partition:
        out.append(f"| Partition | `{resources.partition}` |")
    if resources.account:
        out.append(f"| Account | `{resources.account}` |")
    out += [
        "",
        "**Check these against your cluster before submitting.** Partition "
        "names, GPU names and walltime limits are site-specific, and a job "
        "asking for something the cluster does not have sits in the queue "
        "indefinitely rather than failing.",
        "",
        "```bash",
        "sinfo -s                      # partitions and their limits",
        "sinfo -o '%P %G'              # which partitions have which GPUs",
        "scontrol show config | grep -i maxarray   # array size limit",
        "```",
        "",
    ]
    if not resources.partition:
        out += [
            "No `--partition` was set, so the job goes to your cluster's "
            "default. If that default has no GPU, or a short walltime limit, "
            "add the directive:",
            "",
            "```bash",
            "#SBATCH --partition=<name>",
            "```",
            "",
        ]

    out += ["--------------------------------------------------------------",
            "", "## 2. Tools the cluster must have", ""]
    if guide and guide.tools:
        out += ["Work through these on a **login node**, before submitting.",
                ""]
        for index, tool in enumerate(guide.tools, start=1):
            out += _render_tool(tool, index)
        if resources.module_loads:
            out += [
                "The script already loads these modules, because they were "
                "passed when it was generated:",
                "",
                "```bash",
                *(f"module load {m}" for m in resources.module_loads),
                "```",
                "",
                "If those names do not exist on your cluster, edit the "
                "script or regenerate it with the right `--module` values.",
                "",
            ]
        else:
            out += [
                "The script does **not** load any module, so whatever it "
                "needs must already be on `PATH` when the job starts. Two "
                "ways to arrange that:",
                "",
                "1.  Regenerate the script with the module names, which "
                "writes the `module load` lines into it:",
                "",
                "    ```bash",
                f"    ... generate-hpc-script --job-type {job_type} \\",
                "        --module gromacs/2024.2 --module cuda/12.1 ...",
                "    ```",
                "",
                "2.  Or add them to the script by hand, just below the "
                "`#SBATCH` block.",
                "",
            ]
        if resources.conda_env:
            out += [
                f"The script activates the conda environment "
                f"`{resources.conda_env}`, so create it first:",
                "",
                "```bash",
                f"conda create -n {resources.conda_env} python=3.11 -y",
                f"conda activate {resources.conda_env}",
                "```",
                "",
            ]
    else:
        out += ["No specific tool list is recorded for this job type.", ""]

    out += ["--------------------------------------------------------------",
            "", "## 3. Files that must be in place", ""]
    if guide and guide.inputs:
        out += ["Put these in the directory you submit from:", "",
                "| File | What it is |", "|---|---|"]
        for name, why in guide.inputs:
            resolved = options.get("tpr") if name == "md.tpr" else None
            label = resolved or name
            out.append(f"| `{label}` | {why} |")
        out += [""]
        if guide.produced_by:
            out += [f"Produced by: `{guide.produced_by}`.", ""]
        out += [
            "The script checks for them and exits before doing any work if "
            "one is missing, so a staging mistake costs seconds rather than "
            "a queue wait.",
            "",
        ]
    else:
        out += ["Nothing beyond the script itself.", ""]

    out += ["--------------------------------------------------------------",
            "", "## 4. Submit it", "",
            "```bash",
            f"sbatch {script}",
            "```",
            "",
            "That prints a job id. Use it to follow the job:",
            "",
            "```bash",
            "squeue --me                      # queued and running jobs",
            "squeue -j <jobid> -o '%T %R'     # state, and why it is waiting",
            f"tail -f {resources.job_name}-<jobid>.out    # live output",
            "sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS,ExitCode",
            "scancel <jobid>                  # give up",
            "```",
            ""]
    if resources.array:
        out += [
            "This is an array job, so each task has its own id "
            "(`<jobid>_<task>`) and its own output file. `squeue --me` shows "
            "them collapsed; `squeue -j <jobid> -r` expands them.",
            "",
        ]
    out += [
        "The job writes `"
        + (f"{resources.job_name}-%A_%a.out" if resources.array
           else f"{resources.job_name}-%j.out")
        + "` and the matching `.err` in the submission directory.",
        "",
    ]

    if guide and guide.outputs:
        out += ["## 5. What you get back", "", "```"]
        out += list(guide.outputs)
        out += ["```", ""]

    if guide and guide.resumable:
        out += [
            "## 6. If it hits the walltime",
            "",
            f"The script passes `-maxh {maxh_for(resources.time)}` - just "
            f"under the {resources.time} walltime - so the program stops "
            "cleanly and writes a checkpoint instead of being killed "
            "mid-write. To continue, submit the same script again in the "
            "same directory:",
            "",
            "```bash",
            f"sbatch {script}",
            "```",
            "",
            "It finds the checkpoint and resumes. Nothing needs editing.",
            "",
        ]

    if guide and guide.cautions:
        out += ["## Before you commit a large job", ""]
        out += [f"-   {c}" for c in guide.cautions]
        out += [""]

    out += [
        "## If it fails",
        "",
        "| Message | Means |",
        "|---|---|",
        "| `command not found` | Section 2: the tool is not on `PATH` "
        "inside the job. Loading a module on the login node does not carry "
        "into the job unless the script loads it too. |",
        "| `ERROR: ... not found` | Section 3: an input file is missing. |",
        "| `Invalid partition name` | Section 1: fix `--partition`. |",
        "| `Requested node configuration is not available` | Section 1: the "
        "GPU name, CPU count or memory exceeds anything the partition has. |",
        "| `bash^M: bad interpreter` | The file picked up CRLF endings in "
        "transfer. `dos2unix " + script + "`, or copy it again in binary "
        "mode. The generated file has none. |",
        "| `DUE TO TIME LIMIT` in `sacct` | The walltime was too short. "
        + ("Resubmit to resume." if guide and guide.resumable
           else "Raise `--time` and run again.") + " |",
        "| Job sits in `PD` forever | `squeue -j <jobid> -o '%R'` gives the "
        "reason. `(Resources)` means wait; `(PartitionConfig)` or "
        "`(QOSMaxGRESPerUser)` mean the request is impossible as written. |",
        "",
        "Nothing in this document is cluster-specific guesswork: everything "
        "above is either read from the script's own directives or a command "
        "you run to find out. When in doubt, your site's HPC support knows "
        "the partition names and the module names, which are the two things "
        "that most often need changing.",
    ]
    return "\n".join(out) + "\n"
