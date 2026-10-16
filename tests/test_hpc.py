"""SLURM script generation.

Generated scripts are checked three ways: that bash itself accepts them,
that they contain the directives and guards they claim to, and that invalid
resource requests are refused at generation time rather than at submission.

A generated script is also checked for LF line endings at the byte level. A
shell script written with CRLF fails on a Linux cluster with
``bash\\r: bad interpreter``, and that is precisely the sort of defect that
only appears after the file has been copied to the cluster.
"""

from __future__ import annotations

import subprocess

import pytest

from agskills.errors import UsageError
from agskills.hpc import (
    JOB_TYPES,
    SlurmResources,
    generate_slurm_script,
    validate_memory,
    validate_walltime,
)
from agskills.io_utils import write_text


def _resources(**kwargs) -> SlurmResources:
    defaults = dict(job_name="test", mem="64G", time="24:00:00")
    defaults.update(kwargs)
    return SlurmResources(**defaults)


JOB_SETTINGS = {
    "gromacs-md": {"production_ns": 100, "tpr": "md.tpr"},
    "reinvent": {"config": "run.toml", "log": "reinvent.log"},
    "vina-screen": {"receptor": "r.pdbqt", "ligand_dir": "ligands",
                    "site_config": "site.json"},
    "colabfold": {"fasta": "target.fasta"},
    "alphafold2": {"fasta": "target.fasta", "db_path": "/data/af"},
}


def _script(job_type: str, *, instructions: str | None = None,
            **resource_kwargs) -> str:
    """Render one script. ``instructions`` is a setting, not a resource."""
    if job_type == "vina-screen":
        resource_kwargs.setdefault("array", "1-100")
        resource_kwargs.setdefault("ngpu", 0)
        resource_kwargs.setdefault("gpu", None)
    settings = dict(JOB_SETTINGS[job_type])
    if instructions:
        settings["instructions"] = instructions
    return generate_slurm_script(job_type, _resources(**resource_kwargs),
                                 settings)


# ---------------------------------------------------------------------------
# Syntax
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("job_type", sorted(JOB_TYPES))
def test_generated_script_is_valid_bash(job_type, tmp_path, bash):
    """``bash -n`` must accept every generated script.

    This caught a real defect during development: an ``echo`` line ended
    with ``\\"``, which escaped its own closing quote and left the string
    unterminated.
    """
    path = write_text(_script(job_type), tmp_path / f"{job_type}.sh")
    from conftest import bash_check
    result = bash_check(bash, path)
    assert result.returncode == 0, (
        f"{job_type}.sh is not valid bash:\n{result.stderr}"
    )


@pytest.mark.parametrize("job_type", sorted(JOB_TYPES))
def test_generated_script_has_unix_line_endings(job_type, tmp_path):
    """A CRLF script fails on a Linux cluster before it runs a single step."""
    path = write_text(_script(job_type), tmp_path / f"{job_type}.sh")
    raw = path.read_bytes()
    assert b"\r" not in raw, f"{job_type}.sh contains carriage returns"
    assert raw.startswith(b"#!/bin/bash\n")


@pytest.mark.parametrize("job_type", sorted(JOB_TYPES))
def test_generated_script_aborts_on_error(job_type):
    """Without ``set -e`` a failed step lets the rest run against nothing."""
    assert "set -euo pipefail" in _script(job_type)


# ---------------------------------------------------------------------------
# Directives
# ---------------------------------------------------------------------------


def test_resource_request_appears_in_the_directives():
    script = _script("gromacs-md", ngpu=4, gpu="a100", ncpus=16,
                     mem="128G", time="48:00:00")
    assert "#SBATCH --gres=gpu:a100:4" in script
    assert "#SBATCH --mem=128G" in script
    assert "#SBATCH --cpus-per-task=16" in script
    assert "#SBATCH --time=48:00:00" in script


def test_cpu_only_job_requests_no_gpu():
    script = _script("vina-screen", ngpu=0, gpu=None, ncpus=16)
    assert "--gres=gpu" not in script


def test_array_job_names_its_output_per_task():
    """Without the task id every array task overwrites the same log.

    ``%A_%a`` expands to the array job id and the task index; ``%j`` alone
    is the same value for every task in the array.
    """
    script = _script("vina-screen", array="1-100%10")
    assert "#SBATCH --array=1-100%10" in script
    assert "%A_%a" in script
    assert "--output=test-%j.out" not in script


def test_non_array_job_uses_the_job_id():
    script = _script("reinvent")
    assert "--output=test-%j.out" in script
    assert "--array" not in script


def test_optional_directives_appear_only_when_asked():
    plain = _script("reinvent")
    assert "--partition" not in plain
    assert "--mail-user" not in plain

    detailed = _script("reinvent", partition="gpu", account="proj1",
                       email="me@example.org")
    assert "#SBATCH --partition=gpu" in detailed
    assert "#SBATCH --account=proj1" in detailed
    assert "#SBATCH --mail-user=me@example.org" in detailed
    assert "#SBATCH --mail-type=END,FAIL" in detailed


def test_conda_activation_uses_the_shell_hook():
    """``conda activate`` fails in a non-interactive shell without the hook.

    A SLURM job script is non-interactive, so a bare ``conda activate``
    aborts with "your shell has not been properly configured".
    """
    script = _script("reinvent", conda_env="reinvent4")
    assert "conda shell.bash hook" in script
    assert "conda activate reinvent4" in script
    # And it must fail loudly if conda is absent, not continue regardless.
    assert "ERROR: conda not found" in script


def test_no_environment_activation_when_none_is_requested():
    assert "conda activate" not in _script("reinvent", conda_env=None)


# ---------------------------------------------------------------------------
# Preconditions inside the script
# ---------------------------------------------------------------------------


def test_every_job_checks_its_inputs_before_the_long_step():
    """A missing input must fail in seconds, not after hours of allocation."""
    checks = {
        "gromacs-md": "$TPR",
        "reinvent": "$CONFIG",
        "vina-screen": "$RECEPTOR",
        "colabfold": "$FASTA",
        "alphafold2": "$FASTA",
    }
    for job_type, variable in checks.items():
        script = _script(job_type)
        assert f'if [ ! -f "{variable}" ]; then' in script, job_type
        assert "exit 1" in script


def test_every_job_checks_its_executable_is_present():
    expected = {
        "gromacs-md": "gmx",
        "reinvent": "reinvent",
        "colabfold": "colabfold_batch",
    }
    for job_type, executable in expected.items():
        script = _script(job_type)
        assert f"command -v {executable}" in script, job_type


def test_gromacs_job_resumes_from_a_checkpoint():
    """A 48-hour MD run that hits the walltime must be restartable.

    ``-maxh`` lets GROMACS write a final checkpoint rather than being
    killed mid-write, and ``-cpi`` picks it up on the next attempt.
    """
    script = _script("gromacs-md")
    assert "-maxh" in script
    assert "-cpi" in script
    assert ".cpt" in script


def test_gromacs_job_offloads_to_gpu_only_when_one_is_requested():
    with_gpu = _script("gromacs-md", ngpu=1)
    assert "-nb gpu" in with_gpu
    cpu_only = _script("gromacs-md", ngpu=0, gpu=None)
    assert "-nb gpu" not in cpu_only


def test_vina_array_partitions_the_ligands_without_overlap():
    """Each task must take a disjoint stride through the ligand list.

    A naive split by index range leaves tasks idle when the count does not
    divide evenly; a stride covers every ligand exactly once.
    """
    script = _script("vina-screen", array="1-50")
    assert "SLURM_ARRAY_TASK_ID" in script
    assert "SLURM_ARRAY_TASK_COUNT" in script
    assert "i+=TASK_COUNT" in script
    # Already-docked ligands are skipped, so a requeued task resumes.
    assert "already docked" in script


def test_vina_job_reads_the_box_from_the_site_file():
    """The box must come from the site JSON, not be hardcoded.

    ``jq`` is deliberately avoided: many clusters do not install it, while
    python3 is universally present.
    """
    script = _script("vina-screen")
    assert "center_x" in script
    assert "python3" in script
    # jq may be named in a comment explaining why it is avoided, but it
    # must never be invoked.
    code = [line for line in script.splitlines()
            if not line.lstrip().startswith("#")]
    assert not any("jq " in line or line.strip().endswith("jq")
                   for line in code)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("value", ["24:00:00", "2-12:00:00", "00:30:00",
                                   "168:00:00", "1-00:00:00", "12:30"])
def test_valid_walltimes_are_accepted(value):
    assert validate_walltime(value) == value


@pytest.mark.parametrize("value", ["24h", "1 day", "", "24:00:00:00",
                                   "24:70:00", "abc", "24:00:99"])
def test_invalid_walltimes_are_refused(value):
    with pytest.raises(UsageError):
        validate_walltime(value)


@pytest.mark.parametrize("value,expected", [
    ("64G", "64G"), ("128GB", "128G"), ("4000M", "4000M"), ("1T", "1T"),
    ("64g", "64G"),
])
def test_memory_specifications_are_normalised(value, expected):
    assert validate_memory(value) == expected


@pytest.mark.parametrize("value", ["lots", "64", "", "-4G", "G64"])
def test_invalid_memory_specifications_are_refused(value):
    with pytest.raises(UsageError):
        validate_memory(value)


def test_bad_resources_are_refused_at_construction():
    with pytest.raises(UsageError):
        _resources(ngpu=-1)
    with pytest.raises(UsageError):
        _resources(ncpus=0)
    with pytest.raises(UsageError):
        _resources(time="forever")


@pytest.mark.parametrize("array", ["1-100", "1-100%10", "1,3,5", "1-10:2",
                                   "5"])
def test_valid_array_specifications_are_accepted(array):
    assert _resources(array=array).array == array


@pytest.mark.parametrize("array", ["one to ten", "1..100", "-5", "1-"])
def test_invalid_array_specifications_are_refused(array):
    with pytest.raises(UsageError):
        _resources(array=array)


def test_vina_screen_without_an_array_is_refused():
    """A batch screen needs partitioning; one task would serialise it."""
    with pytest.raises(UsageError) as excinfo:
        generate_slurm_script("vina-screen", _resources(ngpu=0, gpu=None),
                              JOB_SETTINGS["vina-screen"])
    assert "--array" in str(excinfo.value)


def test_reinvent_without_a_config_is_refused():
    with pytest.raises(UsageError) as excinfo:
        generate_slurm_script("reinvent", _resources(), {})
    assert "--config" in str(excinfo.value)


def test_unknown_job_type_lists_the_options():
    with pytest.raises(UsageError) as excinfo:
        generate_slurm_script("mystery", _resources(), {})
    message = str(excinfo.value)
    for job_type in JOB_TYPES:
        assert job_type in message


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("job_type", sorted(JOB_TYPES))
def test_script_documents_what_it_is_and_what_to_review(job_type):
    script = _script(job_type)
    assert JOB_TYPES[job_type].rstrip(".") in script
    assert "agskills" in script
    # Partition and GPU names are site-specific, so the script has to say
    # the directives need checking rather than implying they are ready.
    assert "the directives above" in script
    assert "site-specific" in script


@pytest.mark.parametrize("job_type", sorted(JOB_TYPES))
def test_script_names_its_instructions_when_given_one(job_type):
    """A job log that says 'command not found' should say where to look."""
    script = _script(job_type, instructions="submit_SUBMIT.md")
    assert "submit_SUBMIT.md" in script
    # And it must not be listed as if it were a job parameter.
    settings_line = [line for line in script.splitlines()
                     if line.startswith("# Settings:")]
    assert not any("instructions=" in line for line in settings_line)


@pytest.mark.parametrize("job_type", sorted(JOB_TYPES))
def test_script_says_when_no_module_is_loaded(job_type):
    """Silence about modules is the commonest reason a job cannot find its
    tools: loading one on the login node does not carry into the job."""
    bare = _script(job_type)
    assert "No module is loaded" in bare

    loaded = _script(job_type, module_loads=["gromacs/2024.2", "cuda/12.1"])
    assert "module load gromacs/2024.2" in loaded
    assert "module load cuda/12.1" in loaded
    assert "No module is loaded" not in loaded


def test_script_records_the_settings_it_was_built_from():
    script = _script("gromacs-md", ngpu=2)
    assert "# Settings:" in script
    assert "production_ns=100" in script


def test_script_logs_start_and_finish():
    script = _script("reinvent")
    assert "starting on" in script
    assert "finished at" in script


def test_generation_is_deterministic():
    assert _script("gromacs-md") == _script("gromacs-md")


# ---------------------------------------------------------------------------
# Submission instructions
#
# A script on its own leaves four questions unanswered: which tools the
# cluster needs, how to install them, what to stage, and what to type.
# Every generated script gets a companion document that answers all four.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("job_type", sorted(JOB_TYPES))
def test_every_job_type_has_a_guide(job_type):
    from agskills.hpc_instructions import JOB_GUIDES

    assert job_type in JOB_GUIDES, f"no guide for {job_type}"
    guide = JOB_GUIDES[job_type]
    assert guide.summary
    assert guide.tools, f"{job_type} names no required tool"
    for tool in guide.tools:
        assert tool.probe, f"{job_type}/{tool.name} has no probe command"
        assert tool.module_hint or tool.install, (
            f"{job_type}/{tool.name} says neither how to find nor how to "
            "install it"
        )


@pytest.mark.parametrize("job_type", sorted(JOB_TYPES))
def test_instructions_answer_the_four_questions(job_type):
    from agskills.hpc_instructions import submission_instructions

    resources = _resources(
        **({"array": "1-100", "ngpu": 0, "gpu": None}
           if job_type == "vina-screen" else {}))
    text = submission_instructions(job_type, resources, "submit.sh",
                                   JOB_SETTINGS[job_type])

    # 1. what it asks for, 2. tools, 3. inputs, 4. how to submit
    for heading in ("## 1. What this job asks SLURM for",
                    "## 2. Tools the cluster must have",
                    "## 3. Files that must be in place",
                    "## 4. Submit it"):
        assert heading in text, f"{job_type} instructions lack {heading!r}"

    assert "sbatch submit.sh" in text
    assert "squeue" in text and "sacct" in text and "scancel" in text
    assert "## If it fails" in text
    # The resources quoted must be the ones actually requested.
    assert f"`{resources.time}`" in text
    assert f"`{resources.mem}`" in text


def test_instructions_do_not_promise_a_gpu_to_a_cpu_job():
    from agskills.hpc_instructions import submission_instructions

    resources = _resources(array="1-100", ngpu=0, gpu=None)
    text = submission_instructions("vina-screen", resources, "submit.sh",
                                   JOB_SETTINGS["vina-screen"])
    assert "CPU-only" in text
    assert "a100" not in text
    # An array job needs its own monitoring advice.
    assert "array job" in text
    assert "1-100" in text


def test_instructions_name_the_conda_environment_to_create():
    from agskills.hpc_instructions import submission_instructions

    resources = _resources(conda_env="reinvent4")
    text = submission_instructions("reinvent", resources, "submit.sh",
                                   JOB_SETTINGS["reinvent"])
    assert "conda create -n reinvent4" in text


def test_instructions_tell_a_moduleless_script_how_to_get_modules():
    from agskills.hpc_instructions import submission_instructions

    bare = submission_instructions("gromacs-md", _resources(), "submit.sh",
                                   JOB_SETTINGS["gromacs-md"])
    assert "does **not** load any module" in bare
    assert "--module" in bare

    loaded = submission_instructions(
        "gromacs-md", _resources(module_loads=["gromacs/2024.2"]),
        "submit.sh", JOB_SETTINGS["gromacs-md"])
    assert "module load gromacs/2024.2" in loaded


def test_instructions_contain_no_leaked_format_escapes():
    """A doubled %% belongs in argparse help, not in Markdown."""
    from agskills.hpc_instructions import submission_instructions

    for job_type in sorted(JOB_TYPES):
        resources = _resources(
            **({"array": "1-100", "ngpu": 0, "gpu": None}
               if job_type == "vina-screen" else {}))
        text = submission_instructions(job_type, resources, "submit.sh",
                                       JOB_SETTINGS[job_type])
        assert "%%" not in text, f"{job_type} instructions contain '%%'"


def test_instructions_path_sits_beside_the_script():
    from agskills.hpc_instructions import instructions_path

    assert instructions_path("submit.sh").name == "submit_SUBMIT.md"
    assert instructions_path("a/b/run_md.sh").as_posix() == \
        "a/b/run_md_SUBMIT.md"
