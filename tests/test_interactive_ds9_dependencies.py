"""Unit tests for XPA executable discovery."""

from collections.abc import Mapping

import pytest

from photozpy.interactive_ds9 import dependencies
from photozpy.interactive_ds9.dependencies import (
    XPACommands,
    XPANotAvailableError,
    find_xpa_commands,
)


def test_find_xpa_commands_returns_resolved_executables(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Both discovered executable paths are returned in a typed container."""

    command_paths = {
        "xpaset": "/opt/xpa/bin/xpaset",
        "xpaget": "/opt/xpa/bin/xpaget",
    }
    monkeypatch.setattr(dependencies, "which", command_paths.get)

    result = find_xpa_commands()

    assert result == XPACommands(
        xpaset="/opt/xpa/bin/xpaset",
        xpaget="/opt/xpa/bin/xpaget",
    )


@pytest.mark.parametrize(
    ("available", "expected_missing"),
    [
        ({"xpaget": "/opt/xpa/bin/xpaget"}, "xpaset"),
        ({"xpaset": "/opt/xpa/bin/xpaset"}, "xpaget"),
        ({}, "xpaset, xpaget"),
    ],
)
def test_find_xpa_commands_reports_every_missing_executable(
    monkeypatch: pytest.MonkeyPatch,
    available: Mapping[str, str],
    expected_missing: str,
) -> None:
    """The error identifies exactly which required commands are absent."""

    monkeypatch.setattr(dependencies, "which", available.get)

    with pytest.raises(
        XPANotAvailableError,
        match=rf"XPA commands not found on PATH: {expected_missing}",
    ):
        find_xpa_commands()
