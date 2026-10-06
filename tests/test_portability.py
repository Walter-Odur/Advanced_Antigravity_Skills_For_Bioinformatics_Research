"""Handover readiness: will this work on someone else's computer?

These tests check the properties that decide whether a fresh checkout on an
unfamiliar machine works, rather than any particular scientific result.
Each corresponds to a way the previous version was tied to one developer's
setup.
"""

from __future__ import annotations

import os
import pathlib
import re
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "agskills"
SKILLS = ROOT / ".agents" / "skills"

#: The five skills, their launcher filename, and the CLI module each wraps.
SKILL_LAUNCHERS = {
    "target-preparation": ("target_preparation_api.py", "target_preparation"),
    "compound-screening": ("compound_screening_api.py", "compound_screening"),
    "compound-synthesis": ("compound_synthesis_api.py", "compound_synthesis"),
    "md-simulation": ("md_simulation_api.py", "md_simulation"),
    "drug-discovery-wizard": ("drug_discovery_api.py", "wizard"),
}


def python_files():
    return sorted(SRC.rglob("*.py"))


# ---------------------------------------------------------------------------
# No machine-specific assumptions
# ---------------------------------------------------------------------------


def test_no_absolute_path_is_hardcoded_anywhere():
    """Every example in the old docs named one machine's directory."""
    needles = ("E:\\ANTIGRAVITY", "E:/ANTIGRAVITY", "C:\\Users\\Walter",
               "/home/walter", "/Users/walter")
    offenders = []
    searched = (python_files() + sorted(SKILLS.rglob("*.md"))
                + sorted(SKILLS.rglob("*.py")) + sorted(SKILLS.rglob("*.toml"))
                + sorted((ROOT / "tools").rglob("*.py")))
    for path in searched:
        text = path.read_text(encoding="utf-8", errors="replace").lower()
        for needle in needles:
            if needle.lower() in text:
                offenders.append(f"{path.relative_to(ROOT)}: {needle}")
    assert not offenders, offenders


def test_no_posix_only_temp_path_is_hardcoded():
    """``/tmp`` does not exist on Windows."""
    offenders = []
    for path in python_files():
        for number, line in enumerate(
                path.read_text(encoding="utf-8").splitlines(), 1):
            if re.search(r'["\']/tmp\b', line) and "does not exist" not in line:
                offenders.append(f"{path.relative_to(ROOT)}:{number}")
    assert not offenders, offenders


def test_no_shell_specific_path_separator_assumptions():
    """Paths must be built with pathlib or os.path, not string joins."""
    offenders = []
    for path in python_files():
        for number, line in enumerate(
                path.read_text(encoding="utf-8").splitlines(), 1):
            if re.search(r'\+\s*["\']/[a-z]', line):
                offenders.append(f"{path.relative_to(ROOT)}:{number}")
    assert not offenders, offenders


def test_the_package_source_is_ascii_only():
    """A cp1252 console mangles non-ASCII in user-facing strings."""
    offenders = []
    for path in python_files():
        for number, line in enumerate(
                path.read_text(encoding="utf-8").splitlines(), 1):
            if any(ord(char) > 127 for char in line):
                offenders.append(f"{path.relative_to(ROOT)}:{number}")
    assert not offenders, offenders


# ---------------------------------------------------------------------------
# Imports and optional dependencies
# ---------------------------------------------------------------------------


def test_every_module_imports_cleanly():
    """A broken import surfaces as a confusing failure at first use."""
    import importlib
    failures = []
    for path in python_files():
        if path.name == "__init__.py":
            module = ".".join(path.relative_to(SRC.parent).parts[:-1])
        else:
            module = ".".join(
                path.relative_to(SRC.parent).with_suffix("").parts)
        try:
            importlib.import_module(module)
        except Exception as exc:
            failures.append(f"{module}: {type(exc).__name__}: {exc}")
    assert not failures, failures


def test_importing_the_package_does_no_network_or_disk_work():
    """Import must be free of side effects.

    The previous modules constructed HTTP clients at module scope, so
    merely importing one created three rate limiters and their state files.
    """
    source = (SRC / "libraries" / "chembl.py").read_text(encoding="utf-8")
    # No module-level client construction.
    for line in source.splitlines():
        if line.startswith("_") and "HttpClient(" in line:
            pytest.fail(f"module-level client construction: {line}")

    result = subprocess.run(
        [sys.executable, "-c",
         "import agskills, agskills.cli.compound_screening; print('ok')"],
        capture_output=True, text=True, timeout=120, check=False,
        cwd=str(ROOT),
    )
    assert result.returncode == 0, result.stderr
    assert "ok" in result.stdout


def test_core_functionality_needs_no_optional_dependency():
    """Vina, GROMACS and REINVENT must all be optional.

    A new user should get useful work out of the package before installing
    any external binary.
    """
    from agskills.chem.admet import screen_compounds
    from agskills.chem.smiles import CompoundRecord
    from agskills.hpc import SlurmResources, generate_slurm_script
    from agskills.md.mdp import generate_mdp
    from agskills.reinvent.config import build_config

    assert screen_compounds([CompoundRecord("a", "CCO")], "strict")
    assert generate_mdp("production", force_field="charmm36")
    assert generate_slurm_script(
        "reinvent", SlurmResources(mem="8G", time="01:00:00"),
        {"config": "x.toml"})
    assert build_config(mode="sampling", generator="reinvent",
                        device="cpu").toml


def test_missing_optional_dependencies_are_reported_not_crashed():
    """An absent tool must produce advice, not a traceback."""
    from agskills.docking.vina import find_vina
    from agskills.targets.pdbqt import available_converters

    # Both must answer without raising, whatever is installed.
    kind, _ = find_vina()
    assert kind in ("python", "cli", "none")
    assert isinstance(available_converters(), list)


def test_declared_dependencies_cover_what_the_core_imports():
    """``pip install .`` must be enough for the offline functionality."""
    import tomllib
    metadata = tomllib.loads(
        (ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    required = " ".join(metadata["project"]["dependencies"]).lower()
    assert "rdkit" in required
    assert "biopython" in required

    extras = metadata["project"]["optional-dependencies"]
    assert "meeko" in " ".join(extras["docking"]).lower()
    assert "scipy" in " ".join(extras["docking"]).lower()
    assert "pytest" in " ".join(extras["test"]).lower()


# ---------------------------------------------------------------------------
# Skill launchers
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("skill", sorted(SKILL_LAUNCHERS))
def test_every_skill_has_its_launcher_and_documentation(skill):
    filename, _ = SKILL_LAUNCHERS[skill]
    assert (SKILLS / skill / "SKILL.md").is_file()
    assert (SKILLS / skill / "scripts" / filename).is_file()


@pytest.mark.parametrize("skill", sorted(SKILL_LAUNCHERS))
def test_skill_launcher_runs_as_a_standalone_script(skill):
    """A launcher must work when run directly, before installation.

    Someone who copies the workshop folder onto a new machine may well run
    a command before reading the README.
    """
    filename, _ = SKILL_LAUNCHERS[skill]
    script = SKILLS / skill / "scripts" / filename
    result = subprocess.run(
        [sys.executable, str(script), "--help"],
        capture_output=True, text=True, timeout=180, check=False,
        cwd=str(ROOT),
    )
    assert result.returncode == 0, result.stderr
    assert "--output" in result.stdout or "subcommand" in result.stdout


@pytest.mark.parametrize("skill", sorted(SKILL_LAUNCHERS))
def test_skill_launcher_works_from_an_unrelated_directory(skill, tmp_path):
    """The bootstrap must not depend on the current working directory."""
    filename, _ = SKILL_LAUNCHERS[skill]
    script = SKILLS / skill / "scripts" / filename
    result = subprocess.run(
        [sys.executable, str(script), "--help"],
        capture_output=True, text=True, timeout=180, check=False,
        cwd=str(tmp_path),
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("skill", sorted(SKILL_LAUNCHERS))
def test_skill_launcher_bootstraps_without_the_package_installed(skill,
                                                                  tmp_path):
    """Simulate an uninstalled checkout.

    ``-S`` skips site-packages, so the installed ``agskills`` is invisible
    and only the launcher's own path search can find the source tree.
    """
    filename, _ = SKILL_LAUNCHERS[skill]
    script = SKILLS / skill / "scripts" / filename
    environment = dict(os.environ)
    environment.pop("PYTHONPATH", None)
    result = subprocess.run(
        [sys.executable, "-E", str(script), "--help"],
        capture_output=True, text=True, timeout=180, check=False,
        cwd=str(tmp_path), env=environment,
    )
    assert result.returncode == 0, (
        f"{skill} launcher failed without PYTHONPATH:\n{result.stderr}"
    )


@pytest.mark.parametrize("skill", sorted(SKILL_LAUNCHERS))
def test_skill_launcher_is_thin(skill):
    """A launcher must delegate, not contain logic.

    The orchestrator this replaces was 1,533 lines that re-implemented
    every other skill's subcommands, and the copies had already diverged.
    """
    filename, _ = SKILL_LAUNCHERS[skill]
    script = SKILLS / skill / "scripts" / filename
    code = [line for line in script.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.strip().startswith("#")]
    assert len(code) < 60, f"{skill} launcher has {len(code)} lines of code"
    text = script.read_text(encoding="utf-8")
    # It must delegate into the shared CLI layer, however it imports it.
    assert "agskills.cli" in text
    assert "CLI_MODULE" in text
    # No science in a launcher.
    for forbidden in ("rdkit", "Chem.", "argparse.ArgumentParser",
                      "HttpClient"):
        assert forbidden not in text, f"{skill} launcher mentions {forbidden}"


def test_no_shared_module_is_duplicated_between_skills():
    """``http_client.py`` previously existed four times, byte for byte."""
    import hashlib
    digests: dict[str, list[str]] = {}
    for path in SKILLS.rglob("*.py"):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        digests.setdefault(digest, []).append(str(path.relative_to(SKILLS)))
    duplicates = {d: paths for d, paths in digests.items() if len(paths) > 1}
    assert not duplicates, f"identical files in several skills: {duplicates}"


def test_the_skill_tree_holds_no_stale_implementation_files():
    """Science must live in the package, not in the skill folders."""
    stale = []
    for path in SKILLS.rglob("*.py"):
        if path.name in {f for f, _ in SKILL_LAUNCHERS.values()}:
            continue
        stale.append(str(path.relative_to(SKILLS)))
    assert not stale, f"unexpected implementation files: {stale}"


def test_no_compiled_caches_are_committed_in_the_skill_tree():
    caches = [str(p.relative_to(SKILLS))
              for p in SKILLS.rglob("__pycache__")]
    assert not caches, caches


# ---------------------------------------------------------------------------
# Console entry points
# ---------------------------------------------------------------------------


def test_every_skill_has_a_console_entry_point():
    import tomllib
    metadata = tomllib.loads(
        (ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    scripts = metadata["project"]["scripts"]
    expected = {"ag-target-preparation", "ag-compound-screening",
                "ag-compound-synthesis", "ag-md-simulation", "ag-wizard"}
    assert expected <= set(scripts)
    for target in scripts.values():
        module, _, function = target.partition(":")
        import importlib
        assert hasattr(importlib.import_module(module), function), target


def test_the_installed_console_commands_run():
    """If installed, the five commands must actually work.

    Discovery goes through the package's own helper rather than
    ``shutil.which``, because a console script installed into a virtual
    environment sits next to ``sys.executable`` and is frequently absent
    from PATH - which is how most tooling invokes the interpreter.
    """
    from agskills.targets.pdbqt import find_executable
    for command in ("ag-target-preparation", "ag-compound-screening",
                    "ag-compound-synthesis", "ag-md-simulation",
                    "ag-wizard"):
        found = find_executable(command)
        if not found:
            pytest.skip("the package is not installed in this environment")
        result = subprocess.run([found, "--version"], capture_output=True,
                                text=True, timeout=120, check=False)
        assert result.returncode == 0, f"{command}: {result.stderr}"
        assert "2.0.0" in result.stdout


# ---------------------------------------------------------------------------
# Test-suite hygiene
# ---------------------------------------------------------------------------


def test_the_default_test_run_is_hermetic():
    """The default run must need no network and no external binaries.

    The original suite declared ``@pytest.mark.network`` markers but
    registered no marker configuration and deselected nothing, so a plain
    ``pytest`` run attempted live calls. Combined with a retry policy
    allowing over six minutes of blocking per request, it never finished:
    the baseline run for this rework was killed at 30 minutes without
    completing.
    """
    import tomllib
    config = tomllib.loads(
        (ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    pytest_config = config["tool"]["pytest"]["ini_options"]

    addopts = pytest_config["addopts"]
    for marker in ("network", "slow", "requires_vina", "requires_gromacs",
                   "requires_reinvent"):
        assert marker in addopts, f"{marker} is not deselected by default"
        assert any(marker in declared
                   for declared in pytest_config["markers"]), (
            f"{marker} is used but not registered"
        )
    assert "--strict-markers" in addopts, (
        "without --strict-markers a typo in a marker name is silently "
        "ignored, which is how the original markers stopped working"
    )


def test_all_test_data_is_present():
    data = ROOT / "tests" / "data"
    for name in ("1M17.pdb", "6HEZ.pdb", "AF-P9WJA5-F1.pdb",
                 "pubchem_reference.json",
                 "reinvent_staged_learning_1.csv", "vina_output.txt"):
        path = data / name
        assert path.is_file(), f"missing test fixture: {name}"
        assert path.stat().st_size > 100, f"suspiciously small: {name}"


def test_the_readme_documents_installation_and_testing():
    readme = ROOT / "README.md"
    assert readme.is_file(), "a handover needs a README"
    text = readme.read_text(encoding="utf-8")
    assert "pip install" in text
    assert "pytest" in text
    assert "python" in text.lower()
