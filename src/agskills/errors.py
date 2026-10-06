"""Typed exceptions shared by every skill.

Each exception carries an ``exit_code`` so the CLI layer can translate a
failure into a stable, scriptable process status instead of a bare
traceback. Codes are deliberately distinct and documented:

====  ===========================================================
Code  Meaning
====  ===========================================================
0     Success.
1     Unexpected internal error (a traceback the user should report).
2     Bad usage: argparse-level error, conflicting or missing flags.
3     Invalid input data (unparseable SMILES, malformed PDB, bad FASTA).
4     A required optional dependency or external binary is missing.
5     A remote service failed or was unreachable.
6     A required local file or directory does not exist.
====  ===========================================================
"""

from __future__ import annotations

__all__ = [
    "AgSkillsError",
    "UsageError",
    "InvalidInputError",
    "MissingDependencyError",
    "RemoteServiceError",
    "ResourceNotFoundError",
]


class AgSkillsError(Exception):
    """Base class for every error this package raises deliberately."""

    exit_code: int = 1

    def __init__(self, message: str, *, hint: str | None = None):
        super().__init__(message)
        self.message = message
        self.hint = hint

    def __str__(self) -> str:
        if self.hint:
            return f"{self.message}\nHint: {self.hint}"
        return self.message


class UsageError(AgSkillsError):
    """The caller combined or omitted arguments in an unsupported way."""

    exit_code = 2


class InvalidInputError(AgSkillsError):
    """Input data could not be parsed or failed validation."""

    exit_code = 3


class MissingDependencyError(AgSkillsError):
    """An optional Python package or external executable is unavailable."""

    exit_code = 4

    def __init__(self, what: str, *, install: str | None = None):
        hint = f"Install it with: {install}" if install else None
        super().__init__(f"Required dependency not available: {what}", hint=hint)
        self.what = what
        self.install = install


class RemoteServiceError(AgSkillsError):
    """A remote API returned an error or could not be reached."""

    exit_code = 5


class ResourceNotFoundError(AgSkillsError):
    """A local file or directory that the command needs does not exist."""

    exit_code = 6
