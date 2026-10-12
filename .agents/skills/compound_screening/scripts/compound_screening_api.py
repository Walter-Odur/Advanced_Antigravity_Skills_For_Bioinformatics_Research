#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "polite-http",
#   "rdkit>=2023.9",
#   "biopython>=1.81",
#   "meeko>=0.5",
#   "scipy>=1.10",
#   "gemmi>=0.6",
# ]
# ///
"""Entry point for the compound-screening skill.

Run it with ``uv run``, which reads the PEP 723 header above and resolves
the dependencies into a cached environment. Nothing is installed into the
user's interpreter, nothing needs elevated privilege, and no manual install
step is required::

    uv run scripts/compound_screening_api.py --help

The behaviour lives in the ``agskills`` package, which all five skills
share. The skills this replaced each carried their own copy of the common
modules -- ``http_client.py`` existed four times byte for byte -- and the
copies had already drifted apart, so a bug fixed in one persisted in the
others.

``agskills`` is located by walking up from this file to the repository's
``src`` directory, so the skill runs straight from a checkout. If the
package is installed instead, that installation is used.

Docking additionally needs an AutoDock Vina binary, which is
deliberately not declared above: the ``vina`` wheel builds from
source against Boost and fails on Windows, so declaring it would
break every invocation of the skill rather than just the docking
subcommand. Install it separately with
``conda install -c conda-forge vina`` when you need to dock.
"""

from __future__ import annotations

import sys
from pathlib import Path

#: The module in ``agskills.cli`` that implements this skill.
CLI_MODULE = "compound_screening"


def _bootstrap() -> None:
    """Make ``agskills`` importable, installed or not."""
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
        "Run this script with 'uv run', or install the package from the "
        "workspace root with 'pip install -e .'."
    )


def main() -> int:
    _bootstrap()
    from importlib import import_module
    return import_module(f"agskills.cli.{CLI_MODULE}").main()


if __name__ == "__main__":
    raise SystemExit(main())
