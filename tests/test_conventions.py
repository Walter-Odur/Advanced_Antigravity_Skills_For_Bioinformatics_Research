"""Packaging conventions, checked against the reference bundle's own rules.

The conventions come from the Google DeepMind science-skills bundle vendored
at ``science-skills/`` -- specifically its ``workflow_skill_creator`` skill,
which is its authoritative spec, plus what its 40 skills do in practice.
These tests pin the ones a future edit could silently break.

They also enforce the cost budget. A ``SKILL.md`` is loaded into context
every time its skill triggers, and the JSON a subcommand writes is read back
by the caller, so both are paid for in tokens on every use. Sizes are
asserted rather than hoped for.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SKILLS = ROOT / ".agents" / "skills"
REFERENCE = ROOT / "science-skills" / "skills"

#: Skill folder -> the launcher inside it.
SKILL_SCRIPTS = {
    "target_preparation": "target_preparation_api.py",
    "compound_screening": "compound_screening_api.py",
    "compound_synthesis": "compound_synthesis_api.py",
    "md_simulation": "md_simulation_api.py",
    "drug_discovery_wizard": "drug_discovery_api.py",
}

#: A rough but stable token estimate. Good enough for a budget assertion.
CHARS_PER_TOKEN = 4

#: Ceiling per SKILL.md, in tokens. The reference bundle's median is ~2,139
#: and its largest is ~4,351, so this keeps us at or below their typical
#: cost rather than merely inside their range.
SKILL_MD_TOKEN_BUDGET = 2_700


def skill_dirs() -> list[Path]:
    return sorted(d for d in SKILLS.iterdir() if d.is_dir())


def frontmatter(path: Path) -> dict[str, str]:
    """Parse the YAML frontmatter without requiring a YAML library."""
    text = path.read_text(encoding="utf-8")
    assert text.startswith("---\n"), f"{path} has no frontmatter"
    body = text.split("---\n", 2)[1]
    fields: dict[str, str] = {}
    key = None
    for line in body.splitlines():
        match = re.match(r"^(\w[\w-]*):\s*(.*)$", line)
        if match:
            key = match.group(1)
            fields[key] = match.group(2).strip()
        elif key and line.strip():
            fields[key] = (fields[key] + " " + line.strip()).strip()
    return fields


# ---------------------------------------------------------------------------
# Naming
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("skill", sorted(SKILL_SCRIPTS))
def test_folder_name_is_snake_case(skill):
    """Folders are snake_case, as all 40 reference skills are."""
    directory = SKILLS / skill
    assert directory.is_dir(), f"missing skill folder: {skill}"
    assert re.fullmatch(r"[a-z][a-z0-9_]*", skill), skill
    assert "-" not in skill


@pytest.mark.parametrize("skill", sorted(SKILL_SCRIPTS))
def test_yaml_name_is_kebab_case_and_matches_the_folder(skill):
    """The YAML ``name:`` is the folder name with hyphens.

    That is the exact relationship every reference skill holds, e.g. folder
    ``chembl_database`` with ``name: chembl-database``.
    """
    name = frontmatter(SKILLS / skill / "SKILL.md")["name"]
    assert re.fullmatch(r"[a-z][a-z0-9-]*", name), name
    assert name == skill.replace("_", "-")


@pytest.mark.parametrize("skill", sorted(SKILL_SCRIPTS))
def test_frontmatter_respects_the_documented_limits(skill):
    """name <= 64 chars, description <= 1024, per the skill-creator spec."""
    fields = frontmatter(SKILLS / skill / "SKILL.md")
    assert len(fields["name"]) <= 64
    description = fields.get("description", "")
    assert description, f"{skill} has no description"
    assert len(description) <= 1024, len(description)
    # The description decides whether the skill triggers, so it must say
    # both what it does and when not to use it.
    assert len(description) >= 120, f"{skill}: description too thin to route on"


def test_the_reference_bundle_really_uses_this_convention():
    """Guard the premise: if the reference changed, these tests should know.

    Skipped rather than failed when the bundle is absent, so the suite still
    runs on a checkout without it.
    """
    if not REFERENCE.is_dir():
        pytest.skip("the science-skills reference bundle is not present")
    checked = 0
    for directory in sorted(REFERENCE.iterdir()):
        skill_md = directory / "SKILL.md"
        if not skill_md.is_file():
            continue
        name = frontmatter(skill_md).get("name", "")
        assert name == directory.name.replace("_", "-"), (
            f"{directory.name}: reference convention differs from assumed"
        )
        checked += 1
    assert checked >= 30, f"only inspected {checked} reference skills"


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("skill", sorted(SKILL_SCRIPTS))
def test_skill_has_the_expected_layout(skill):
    """SKILL.md, scripts/ and references/ -- the latter plural."""
    directory = SKILLS / skill
    assert (directory / "SKILL.md").is_file()
    assert (directory / "scripts" / SKILL_SCRIPTS[skill]).is_file()
    assert (directory / "references").is_dir()
    # 'reference' singular was the old spelling and must not come back.
    assert not (directory / "reference").exists()


@pytest.mark.parametrize("skill", sorted(SKILL_SCRIPTS))
def test_skill_ships_citations(skill):
    """Every reference skill that touches a data source ships citation.bib."""
    bib = SKILLS / skill / "references" / "citation.bib"
    assert bib.is_file(), f"{skill} has no references/citation.bib"
    text = bib.read_text(encoding="utf-8")
    entries = re.findall(r"@(\w+)\{([^,]+),", text)
    assert len(entries) >= 5, f"{skill}: only {len(entries)} citations"
    keys = [k for _, k in entries]
    assert len(keys) == len(set(keys)), f"{skill}: duplicate citation keys"
    assert text.count("{") == text.count("}"), f"{skill}: unbalanced braces"
    # A citation without a resolvable identifier is not much use.
    dois = len(re.findall(r"^\s*DOI\s*=", text, re.M))
    urls = len(re.findall(r"^\s*URL\s*=", text, re.M))
    assert dois + urls >= len(entries), f"{skill}: citations lack DOI/URL"


def test_bundle_has_its_manifest_and_licence_table():
    manifest = json.loads((SKILLS / "plugin.json").read_text(encoding="utf-8"))
    for field in ("name", "version", "description"):
        assert manifest.get(field), f"plugin.json lacks {field}"
    assert re.fullmatch(r"\d+\.\d+\.\d+", manifest["version"])

    licences = (SKILLS / "SKILL_LICENSES.md").read_text(encoding="utf-8")
    for skill in SKILL_SCRIPTS:
        assert skill.replace("_", "-") in licences, (
            f"{skill} is missing from SKILL_LICENSES.md"
        )


def test_no_tests_are_shipped_inside_the_skills():
    """Tests confirm the skills work; they are not part of them.

    Anything under ``.agents/skills`` is shipped, so a stray test file there
    is dead weight a user downloads and an agent may read.
    """
    strays = [
        str(p.relative_to(SKILLS))
        for p in SKILLS.rglob("*")
        if p.is_file() and (
            p.name.startswith("test_")
            or p.name in {"conftest.py", "pytest.ini", "panel.py"}
            or "__pycache__" in p.parts
        )
    ]
    assert not strays, f"test or cache artefacts shipped in the skills: {strays}"


def test_only_launchers_live_in_the_skill_scripts_directories():
    """Implementation lives in the package, not duplicated per skill."""
    expected = set(SKILL_SCRIPTS.values())
    found = {p.name for p in SKILLS.rglob("scripts/*.py")}
    assert found == expected, f"unexpected scripts: {found - expected}"


# ---------------------------------------------------------------------------
# PEP 723 and uv
# ---------------------------------------------------------------------------


def pep723(path: Path) -> dict:
    """Extract and parse the PEP 723 inline metadata block."""
    import tomllib

    text = path.read_text(encoding="utf-8")
    match = re.search(r"^# /// script\s*$(.+?)^# ///\s*$", text,
                      re.M | re.S)
    assert match, f"{path} has no PEP 723 block"
    stripped = "".join(
        line[2:] if line.startswith("# ") else line[1:] + "\n"
        for line in match.group(1).splitlines(keepends=True)
    )
    return tomllib.loads(stripped)


@pytest.mark.parametrize("skill", sorted(SKILL_SCRIPTS))
def test_launcher_declares_valid_pep723_metadata(skill):
    """63 of the reference bundle's 66 scripts use this; it is how uv resolves.

    The block must be valid TOML, because ``uv run`` parses it and a
    malformed header fails the skill before it starts.
    """
    meta = pep723(SKILLS / skill / "scripts" / SKILL_SCRIPTS[skill])
    assert meta.get("requires-python"), f"{skill} declares no Python floor"
    for dep in meta.get("dependencies", []):
        assert re.fullmatch(r"[A-Za-z0-9._-]+([<>=!~].+)?", dep), dep


def test_declared_dependencies_are_installable_and_bounded():
    """Nothing that cannot be installed may be declared.

    ``vina`` builds from source against Boost and fails on Windows;
    ``reinvent`` is a git checkout and pulls PyTorch. Declaring either would
    break every invocation of the skill rather than just the one subcommand
    that needs it, so both stay optional prerequisites.
    """
    forbidden = {"vina", "reinvent", "torch", "gromacs", "openbabel"}
    for skill, script in SKILL_SCRIPTS.items():
        deps = pep723(SKILLS / skill / "scripts" / script).get(
            "dependencies", [])
        roots = {re.split(r"[<>=!~\[]", d)[0].lower() for d in deps}
        assert not (roots & forbidden), (
            f"{skill} declares an uninstallable dependency: "
            f"{roots & forbidden}"
        )


def test_md_simulation_needs_no_third_party_dependency():
    """It runs on the standard library, which is why it works anywhere.

    The reference spec says to prefer stdlib and avoid libraries that need
    extra installation where possible; this skill is the case where that is
    fully achievable.
    """
    meta = pep723(SKILLS / "md_simulation" / "scripts" /
                  SKILL_SCRIPTS["md_simulation"])
    assert not meta.get("dependencies"), meta.get("dependencies")


@pytest.mark.parametrize("skill", sorted(SKILL_SCRIPTS))
@pytest.mark.slow
def test_launcher_runs_under_uv(skill):
    """``uv run`` must resolve the header and start the skill.

    This is the invocation every SKILL.md documents, so it is the one that
    has to work. Marked slow because the first run populates uv's cache.
    """
    if not shutil.which("uv"):
        pytest.skip("uv is not installed")
    script = SKILLS / skill / "scripts" / SKILL_SCRIPTS[skill]
    result = subprocess.run(
        ["uv", "run", "--no-project", str(script), "--help"],
        capture_output=True, text=True, timeout=900, check=False,
        cwd=str(ROOT),
    )
    assert result.returncode == 0, (
        f"uv run failed for {skill}:\n{result.stderr[-2000:]}"
    )
    assert "--output" in result.stdout or "subcommand" in result.stdout


# ---------------------------------------------------------------------------
# SKILL.md structure
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("skill", sorted(SKILL_SCRIPTS))
def test_skill_md_has_the_canonical_sections(skill):
    """The section set the reference spec's Rule 6 prescribes."""
    text = (SKILLS / skill / "SKILL.md").read_text(encoding="utf-8")
    for heading in ("## Prerequisites", "## Core Rules", "## Utility Script",
                    "## Common Options", "## Common Mistakes",
                    "## Reference", "## Workflow"):
        assert heading in text, f"{skill} lacks '{heading}'"


@pytest.mark.parametrize("skill", sorted(SKILL_SCRIPTS))
def test_skill_md_names_the_uv_prerequisite(skill):
    """Every reference skill that runs code points at the ``uv`` skill."""
    text = (SKILLS / skill / "SKILL.md").read_text(encoding="utf-8")
    assert "`uv` skill" in text
    assert "uv run scripts/" in text


@pytest.mark.parametrize("skill", sorted(SKILL_SCRIPTS))
def test_skill_md_follows_the_licence_notification_pattern(skill):
    """``.licenses/<skill>_LICENSE.txt``, as chembl_database does."""
    text = (SKILLS / skill / "SKILL.md").read_text(encoding="utf-8")
    assert f".licenses/{skill}_LICENSE.txt" in text
    assert "timestamp" in text


@pytest.mark.parametrize("skill", sorted(SKILL_SCRIPTS))
def test_skill_md_requires_output_and_says_so(skill):
    text = (SKILLS / skill / "SKILL.md").read_text(encoding="utf-8")
    assert "**Required.**" in text or "is required for every" in text
    assert "--output" in text


@pytest.mark.parametrize("skill", sorted(SKILL_SCRIPTS))
def test_skill_md_asks_to_be_mentioned_in_the_output(skill):
    """The reference skills' own notification rule."""
    text = (SKILLS / skill / "SKILL.md").read_text(encoding="utf-8")
    assert "mention" in text.lower() and "output" in text.lower()


@pytest.mark.parametrize("skill", sorted(SKILL_SCRIPTS))
def test_skill_md_links_only_to_files_that_exist(skill):
    """A broken reference link wastes a read and misleads the agent."""
    directory = SKILLS / skill
    text = (directory / "SKILL.md").read_text(encoding="utf-8")
    missing = []
    for target in re.findall(r"\]\((references/[^)#]+|\.\./[^)#]+)\)", text):
        if not (directory / target).resolve().is_file():
            missing.append(target)
    assert not missing, f"{skill} links to missing files: {missing}"


# ---------------------------------------------------------------------------
# Token economy
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("skill", sorted(SKILL_SCRIPTS))
def test_skill_md_stays_within_the_token_budget(skill):
    """A SKILL.md is loaded on every trigger, so its size is a per-use cost.

    Depth belongs in ``references/``, which is read only when needed.
    """
    size = (SKILLS / skill / "SKILL.md").stat().st_size
    tokens = size // CHARS_PER_TOKEN
    assert tokens <= SKILL_MD_TOKEN_BUDGET, (
        f"{skill}: SKILL.md is ~{tokens} tokens, over the "
        f"{SKILL_MD_TOKEN_BUDGET} budget. Move detail into references/."
    )


def test_compact_output_is_an_order_of_magnitude_cheaper_than_full():
    """The default detail level must actually save context.

    Measured on a 50-compound strict screen: the per-check array alone is
    ~84% of the full output, which comes to roughly 1,400 tokens per
    compound. ``compact`` is the default because the caller reads this JSON
    back.
    """
    from agskills.chem.admet import screen_compounds
    from agskills.chem.smiles import CompoundRecord

    from panel import EGFR_ACTIVES

    records = [CompoundRecord(f"c{i}", EGFR_ACTIVES[i % len(EGFR_ACTIVES)])
               for i in range(50)]
    compact = len(json.dumps(screen_compounds(records, "strict", "compact")))
    full = len(json.dumps(screen_compounds(records, "strict", "full")))

    assert compact * 8 < full, (
        f"compact is only {full / compact:.1f}x cheaper than full; "
        "the saving should be around an order of magnitude"
    )
    per_compound = compact / CHARS_PER_TOKEN / len(records)
    assert per_compound < 200, (
        f"compact costs ~{per_compound:.0f} tokens per compound"
    )


def test_every_detail_level_reports_the_same_verdicts():
    """Cheaper output must not mean different answers.

    Only the per-compound payload changes between levels; counts, verdicts
    and ordering are identical.
    """
    from agskills.chem.admet import screen_compounds
    from agskills.chem.smiles import CompoundRecord

    from panel import EGFR_ACTIVES

    records = [CompoundRecord(f"c{i}", s) for i, s in enumerate(EGFR_ACTIVES)]
    records.append(CompoundRecord("bad", "NOT_A_MOLECULE"))

    reports = {level: screen_compounds(records, "strict", level)
               for level in ("compact", "standard", "full")}
    reference = reports["compact"]
    for level, report in reports.items():
        for key in ("total_compounds", "passed", "failed", "errored",
                    "checks_per_compound"):
            assert report[key] == reference[key], (level, key)
        assert [c["name"] for c in report["compounds"]] == \
               [c["name"] for c in reference["compounds"]], level
        assert [c["overall_pass"] for c in report["compounds"]] == \
               [c["overall_pass"] for c in reference["compounds"]], level


def test_the_fixed_tier_definition_is_not_embedded_by_default():
    """~830 tokens of identical overhead in every report is worth removing.

    It is reachable on demand through ``describe-tiers``, which the compact
    output names.
    """
    from agskills.chem.admet import screen_compounds
    from agskills.chem.smiles import CompoundRecord

    compact = screen_compounds([CompoundRecord("x", "CCO")], "strict",
                               "compact")
    assert "tier_definition" not in compact
    assert "describe-tiers" in compact["tier_definition_command"]

    full = screen_compounds([CompoundRecord("x", "CCO")], "strict", "full")
    assert "tier_definition" in full


@pytest.mark.parametrize("skill", sorted(SKILL_SCRIPTS))
def test_skill_md_warns_against_printing_whole_files(skill):
    """The reference spec's Rule 4: stdout stays short, results go to a file."""
    text = (SKILLS / skill / "SKILL.md").read_text(encoding="utf-8")
    assert "--output" in text
    # Each skill must document the detail control that bounds its output.
    assert "--detail" in text, f"{skill} does not document --detail"


def test_stdout_summary_is_bounded(tmp_path, capsys):
    """A status summary must stay short however many compounds are screened.

    stdout goes straight into the caller's context, so a 1,000-compound run
    must not print 1,000 lines.
    """
    from agskills.cli import compound_screening

    from panel import EGFR_ACTIVES

    smiles = [EGFR_ACTIVES[i % len(EGFR_ACTIVES)] for i in range(300)]
    code = compound_screening.main([
        "admet-filter", "--smiles", *smiles,
        "--strictness", "relaxed",
        "--output", str(tmp_path / "out.json"),
    ])
    assert code == 0
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) < 60, f"stdout printed {len(lines)} lines for 300 compounds"
