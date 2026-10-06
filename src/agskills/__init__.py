"""agskills - structure-based drug discovery skills.

Five skills share this one implementation:

``target-preparation``
    Fetch, clean, assess and prepare protein targets; define docking sites.
``compound-screening``
    Search compound libraries, dock, and screen with a two-tier ADMET
    battery.
``compound-synthesis``
    Generate and optimise molecules with REINVENT 4.
``md-simulation``
    Generate validated GROMACS workflows and cluster submission scripts.
``drug-discovery-wizard``
    Orchestrate the four above end to end.

Everything is importable as a library as well as runnable from the command
line, and no function writes to ``sys.exit`` or prints in place of
returning a value.
"""

from __future__ import annotations

__version__ = "2.0.0"

from .errors import (
    AgSkillsError,
    InvalidInputError,
    MissingDependencyError,
    RemoteServiceError,
    ResourceNotFoundError,
    UsageError,
)

__all__ = [
    "__version__",
    "AgSkillsError",
    "UsageError",
    "InvalidInputError",
    "MissingDependencyError",
    "RemoteServiceError",
    "ResourceNotFoundError",
]
