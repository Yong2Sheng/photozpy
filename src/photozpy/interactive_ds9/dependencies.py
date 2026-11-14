"""Dependency discovery helpers for the DS9/XPA integration.

This module keeps external-program checks separate from the DS9 adapter
itself.  The interactive reviewer currently talks to DS9 through the command
line XPA tools instead of importing a Python DS9 binding, so these helpers
verify that the required executables are available on ``PATH``.
"""

from dataclasses import dataclass
from shutil import which


@dataclass(frozen=True)
class XPACommands:
    """Absolute command paths for the XPA executables used by Photozpy."""

    xpaset: str
    xpaget: str


class XPANotAvailableError(RuntimeError):
    """Raised when one or more required XPA commands cannot be found."""


def find_xpa_commands() -> XPACommands:
    """Find the XPA command line tools required to communicate with DS9.

    Returns
    -------
    XPACommands
        The resolved executable paths for ``xpaset`` and ``xpaget``.

    Raises
    ------
    XPANotAvailableError
        If either command is missing from the current ``PATH``.
    """

    commands = {
        "xpaset": which("xpaset"),
        "xpaget": which("xpaget"),
    }
    missing = [name for name, path in commands.items() if path is None]

    if missing:
        raise XPANotAvailableError(
            "XPA commands not found on PATH: " + ", ".join(missing)
        )

    return XPACommands(
        xpaset=commands["xpaset"],
        xpaget=commands["xpaget"],
    )
