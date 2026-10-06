#!/usr/bin/env python3
"""Launcher for the compound-synthesis skill.

This file is a thin entry point. All behaviour lives in the ``agskills``
package, which every skill shares, so a fix reaches all of them at once
rather than being applied to five diverged copies.

The skills this replaces each carried their own copy of the shared
modules - ``http_client.py`` existed four times byte for byte, and the
orchestrator re-implemented every other skill's subcommands. The copies
had already drifted: a docking argument bug was present in one and absent
in another.

Run it either way:

    # Installed (recommended):
    pip install -e <workshop root>
    ag-compound-synthesis --help

    # Or directly from a checkout, with no installation:
    python .agents/skills/compound-synthesis/scripts/compound_synthesis_api.py --help
"""

from __future__ import annotations

import sys
from pathlib import Path

#: The CLI module in ``agskills.cli`` that this skill exposes.
CLI_MODULE = "compound_synthesis"


def _bootstrap() -> None:
    """Make ``agskills`` importable whether or not it is installed.

    Walks up from this file looking for a ``src/agskills`` package and puts
    it on ``sys.path``. This keeps the skill runnable straight from a
    checkout, which matters when someone copies the workshop folder onto a
    new machine and runs a command before installing anything.
    """
    try:
        import agskills  # noqa: F401
        return
    except ImportError:
        pass

    for parent in Path(__file__).resolve().parents:
        candidate = parent / "src"
        if (candidate / "agskills" / "__init__.py").is_file():
            sys.path.insert(0, str(candidate))
            return

    sys.exit(
        "ERROR: the 'agskills' package could not be found.\n"
        "Install it from the workshop root:\n"
        "    pip install -e .\n"
        "or set PYTHONPATH to the directory that contains 'agskills'."
    )


def main() -> int:
    _bootstrap()
    from importlib import import_module
    return import_module(f"agskills.cli.{CLI_MODULE}").main()


if __name__ == "__main__":
    raise SystemExit(main())
