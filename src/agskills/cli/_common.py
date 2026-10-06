"""Shared CLI plumbing: error handling, exit codes and output contracts.

Every subcommand of every skill honours the same contract:

* ``--output`` is required and receives the machine-readable result.
* stdout carries a short human summary; the file carries the data.
* A deliberate failure prints a one-line message plus a hint and exits with
  a documented status code (see :mod:`agskills.errors`). It does not print
  a traceback, because the exception types are part of the interface.
* ``--quiet`` suppresses the summary; ``--verbose`` enables library logging.

The code this replaces called ``sys.exit(1)`` from inside the science
functions after printing to stderr, which made every failure mode
indistinguishable to a caller and made the functions untestable without
capturing output and catching ``SystemExit``.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Any, Callable, Sequence

from ..errors import AgSkillsError, UsageError
from ..io_utils import write_json, write_text

__all__ = [
    "build_root_parser",
    "add_output",
    "add_compound_inputs",
    "add_common_flags",
    "compounds_from_args",
    "emit_json",
    "emit_text",
    "run_cli",
    "summarise",
]


class _HelpfulParser(argparse.ArgumentParser):
    """An ArgumentParser whose errors point at the right command."""

    def error(self, message: str):  # noqa: D102
        self.print_usage(sys.stderr)
        sys.stderr.write(f"\nerror: {message}\n")
        sys.stderr.write(
            f"\nRun '{self.prog} --help' to see the available options.\n"
        )
        raise SystemExit(UsageError.exit_code)


def build_root_parser(prog: str, description: str, epilog: str = ""
                      ) -> tuple[argparse.ArgumentParser, Any]:
    """Create the root parser and its subcommand registry."""
    parser = _HelpfulParser(
        prog=prog, description=description, epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    from .. import __version__
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s (agskills {__version__})")
    sub = parser.add_subparsers(dest="command", required=True,
                                metavar="<subcommand>")
    return parser, sub


def add_output(parser: argparse.ArgumentParser, *, help_text: str,
               required: bool = True) -> None:
    """Add the mandatory ``--output`` argument."""
    parser.add_argument("--output", required=required, metavar="FILE",
                        help=help_text)


def add_common_flags(parser: argparse.ArgumentParser) -> None:
    """Add flags every subcommand shares."""
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-q", "--quiet", action="store_true",
                       help="Suppress the human-readable summary on stdout.")
    group.add_argument("-v", "--verbose", action="store_true",
                       help="Log library activity, including HTTP retries.")


def add_compound_inputs(parser: argparse.ArgumentParser, *,
                        required: bool = True) -> None:
    """Add the standard trio of compound input options."""
    group = parser.add_mutually_exclusive_group(required=required)
    group.add_argument("--smiles", nargs="+", metavar="SMILES",
                       help="One or more SMILES strings.")
    group.add_argument("--csv", metavar="FILE",
                       help="CSV file with a SMILES column.")
    group.add_argument("--smiles-file", metavar="FILE",
                       help="SMILES file: one structure per line, with an "
                            "optional name in the second column.")
    parser.add_argument("--names", nargs="+", metavar="NAME",
                        help="Compound names for --smiles. Must match the "
                             "number of SMILES exactly.")
    parser.add_argument("--smiles-col", default="SMILES", metavar="COL",
                        help="SMILES column name in --csv (default: SMILES).")
    parser.add_argument("--name-col", default="Name", metavar="COL",
                        help="Name column in --csv (default: Name).")


def compounds_from_args(args: argparse.Namespace):
    """Build the compound list from parsed arguments."""
    from ..chem.smiles import load_compounds
    return load_compounds(
        smiles=getattr(args, "smiles", None),
        names=getattr(args, "names", None),
        csv_path=getattr(args, "csv", None),
        smiles_file=getattr(args, "smiles_file", None),
        smiles_col=getattr(args, "smiles_col", "SMILES"),
        name_col=getattr(args, "name_col", "Name"),
    )


def emit_json(data: Any, args: argparse.Namespace,
              summary: Sequence[str] = ()) -> Path:
    """Write the JSON result and print the summary."""
    path = write_json(data, args.output)
    _print_summary(path, summary, args)
    return path


def emit_text(text: str, args: argparse.Namespace,
              summary: Sequence[str] = (), *, executable: bool = False) -> Path:
    """Write a text artefact (script, TOML, SMILES) and print the summary."""
    path = write_text(text, args.output, executable=executable)
    _print_summary(path, summary, args)
    return path


def _print_summary(path: Path, summary: Sequence[str],
                   args: argparse.Namespace) -> None:
    if getattr(args, "quiet", False):
        return
    for line in summary:
        print(line)
    print(f"Written: {path}")


def summarise(result: dict[str, Any], keys: Sequence[tuple[str, str]]
              ) -> list[str]:
    """Build summary lines from ``(label, key)`` pairs present in *result*."""
    lines = []
    for label, key in keys:
        if key in result and result[key] is not None:
            lines.append(f"  {label}: {result[key]}")
    return lines


def run_cli(parser: argparse.ArgumentParser,
            argv: Sequence[str] | None = None) -> int:
    """Parse arguments, dispatch, and translate exceptions into exit codes.

    Returns the process exit status rather than calling ``sys.exit``, so
    tests can assert on it directly.
    """
    args = parser.parse_args(argv)

    if getattr(args, "verbose", False):
        logging.basicConfig(
            level=logging.INFO, stream=sys.stderr,
            format="%(levelname)s %(name)s: %(message)s",
        )

    handler: Callable[[argparse.Namespace], Any] | None = getattr(
        args, "handler", None)
    if handler is None:  # pragma: no cover - argparse enforces a subcommand
        parser.error("no subcommand selected")
        return UsageError.exit_code

    try:
        handler(args)
    except AgSkillsError as exc:
        print(f"error: {exc.message}", file=sys.stderr)
        if exc.hint:
            print(f"hint: {exc.hint}", file=sys.stderr)
        return exc.exit_code
    except (KeyError, ValueError) as exc:
        # The science layer raises these for bad enum values and bad data.
        print(f"error: {exc}", file=sys.stderr)
        return UsageError.exit_code
    except BrokenPipeError:  # pragma: no cover
        return 0
    except KeyboardInterrupt:
        print("\ninterrupted", file=sys.stderr)
        return 130
    return 0
