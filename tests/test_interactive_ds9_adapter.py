"""Unit tests for the subprocess-based DS9/XPA adapter."""

from pathlib import Path
import subprocess
from unittest.mock import Mock, call

import pytest

from photozpy.interactive_ds9 import adapter as adapter_module
from photozpy.interactive_ds9.adapter import (
    DS9Adapter,
    DS9CommunicationError,
)
from photozpy.interactive_ds9.dependencies import XPACommands


@pytest.fixture
def xpa_commands() -> XPACommands:
    """Return deterministic executable paths for adapter tests."""

    return XPACommands(
        xpaset="/opt/xpa/bin/xpaset",
        xpaget="/opt/xpa/bin/xpaget",
    )


@pytest.fixture
def ds9_adapter(
    monkeypatch: pytest.MonkeyPatch,
    xpa_commands: XPACommands,
) -> DS9Adapter:
    """Build an adapter without depending on the host XPA installation."""

    monkeypatch.setattr(
        adapter_module,
        "find_xpa_commands",
        lambda: xpa_commands,
    )
    return DS9Adapter()


def test_adapter_discovers_commands_during_initialization(
    monkeypatch: pytest.MonkeyPatch,
    xpa_commands: XPACommands,
) -> None:
    """Initialization records the access-point name and discovered commands."""

    find_commands = Mock(return_value=xpa_commands)
    monkeypatch.setattr(adapter_module, "find_xpa_commands", find_commands)

    adapter = DS9Adapter(
        software_name="photozpy-ds9",
        reset_scale_on_load=False,
    )

    assert adapter.software_name == "photozpy-ds9"
    assert adapter.commands == xpa_commands
    assert adapter.reset_scale_on_load is False
    find_commands.assert_called_once_with()


def test_ensure_running_reuses_existing_target(ds9_adapter, monkeypatch):
    monkeypatch.setattr(ds9_adapter, 'get', Mock(return_value='8.6'))
    start = Mock()
    monkeypatch.setattr(adapter_module.subprocess, 'Popen', start)
    ds9_adapter.ensure_running()
    start.assert_not_called()


def test_ensure_running_starts_named_target(ds9_adapter, monkeypatch):
    ds9_adapter.software_name = 'photozpy-cr-preview'
    monkeypatch.setattr(ds9_adapter, 'get', Mock(side_effect=[DS9CommunicationError('no target'), '8.6']))
    monkeypatch.setattr(adapter_module, 'which', lambda name: '/opt/ds9/bin/ds9')
    start = Mock()
    monkeypatch.setattr(adapter_module.subprocess, 'Popen', start)
    ds9_adapter.ensure_running()
    assert start.call_args.args[0] == ['/opt/ds9/bin/ds9', '-title', 'photozpy-cr-preview']


def test_ensure_running_reports_launch_failure(ds9_adapter, monkeypatch):
    monkeypatch.setattr(ds9_adapter, 'get', Mock(side_effect=DS9CommunicationError('no target')))
    monkeypatch.setattr(adapter_module, 'which', lambda name: '/opt/ds9/bin/ds9')
    process = Mock()
    process.poll.return_value = 1
    monkeypatch.setattr(adapter_module.subprocess, 'Popen', Mock(return_value=process))
    with pytest.raises(DS9CommunicationError, match='did not become available'):
        ds9_adapter.ensure_running()


@pytest.mark.parametrize('response', ['0', 'no', '', 'false'])
def test_empty_frame_is_not_reported_as_success(ds9_adapter, monkeypatch, response):
    query = Mock(return_value=response)
    monkeypatch.setattr(ds9_adapter, 'get', query)
    with pytest.raises(DS9CommunicationError, match='no FITS image'):
        ds9_adapter.verify_image_loaded()
    query.assert_called_once_with('frame', 'has', 'fits')


@pytest.mark.parametrize('response', ['1', 'yes', 'true'])
def test_image_presence_check_accepts_loaded_frame(ds9_adapter, monkeypatch, response):
    monkeypatch.setattr(ds9_adapter, 'get', Mock(return_value=response))
    ds9_adapter.verify_image_loaded()


def test_set_runs_xpaset_with_separate_command_tokens(
    monkeypatch: pytest.MonkeyPatch,
    ds9_adapter: DS9Adapter,
) -> None:
    """Adapter arguments are passed directly to subprocess without a shell."""

    completed = subprocess.CompletedProcess(
        args=[],
        returncode=0,
        stdout="",
        stderr="",
    )
    run = Mock(return_value=completed)
    monkeypatch.setattr(adapter_module.subprocess, "run", run)

    ds9_adapter.set("scale", "mode", "zscale")

    run.assert_called_once_with(
        [
            "/opt/xpa/bin/xpaset",
            "-p",
            "ds9",
            "scale",
            "mode",
            "zscale",
        ],
        check=False,
        text=True,
        capture_output=True,
    )


def test_set_raises_communication_error_with_stderr(
    monkeypatch: pytest.MonkeyPatch,
    ds9_adapter: DS9Adapter,
) -> None:
    """A failed xpaset process is translated into the adapter exception."""

    completed = subprocess.CompletedProcess(
        args=[],
        returncode=1,
        stdout="",
        stderr="  no xpaset access points match template: ds9  \n",
    )
    monkeypatch.setattr(
        adapter_module.subprocess,
        "run",
        Mock(return_value=completed),
    )

    with pytest.raises(
        DS9CommunicationError,
        match="no xpaset access points match template: ds9",
    ):
        ds9_adapter.set("fits", "image.fits")


def test_get_returns_stripped_xpaget_output(
    monkeypatch: pytest.MonkeyPatch,
    ds9_adapter: DS9Adapter,
) -> None:
    """Successful xpaget output is stripped before it reaches callers."""

    completed = subprocess.CompletedProcess(
        args=[],
        returncode=0,
        stdout="  ds9 8.6  \n",
        stderr="",
    )
    run = Mock(return_value=completed)
    monkeypatch.setattr(adapter_module.subprocess, "run", run)

    result = ds9_adapter.get("version")

    assert result == "ds9 8.6"
    run.assert_called_once_with(
        ["/opt/xpa/bin/xpaget", "ds9", "version"],
        check=False,
        text=True,
        capture_output=True,
    )


def test_get_raises_communication_error_with_stderr(
    monkeypatch: pytest.MonkeyPatch,
    ds9_adapter: DS9Adapter,
) -> None:
    """A failed xpaget process is translated into the adapter exception."""

    completed = subprocess.CompletedProcess(
        args=[],
        returncode=1,
        stdout="",
        stderr="  invalid host name specified  \n",
    )
    monkeypatch.setattr(
        adapter_module.subprocess,
        "run",
        Mock(return_value=completed),
    )

    with pytest.raises(
        DS9CommunicationError,
        match="invalid host name specified",
    ):
        ds9_adapter.get("version")


def test_load_fits_quotes_path_and_applies_display_defaults(
    tmp_path: Path,
    ds9_adapter: DS9Adapter,
) -> None:
    """Loading uses a DS9-safe path before scale and zoom commands."""

    image_path = tmp_path / 'night with spaces and "quotes"' / "image.fits"
    set_command = Mock()
    ds9_adapter.set = set_command

    ds9_adapter.load_fits(
        image_path,
        scale_mode="minmax",
        scale_function="sqrt",
    )

    escaped_path = str(image_path.resolve()).replace("\\", "\\\\").replace(
        '"', '\\"'
    )
    assert set_command.call_args_list == [
        call("fits", f'"{escaped_path}"'),
        call("scale", "mode", "minmax"),
        call("scale", "sqrt"),
        call("zoom", "to", "fit"),
    ]


def test_load_fits_keeps_manual_scale_when_advancing_or_going_back(
    tmp_path: Path,
    ds9_adapter: DS9Adapter,
) -> None:
    """Only the first display initializes scale; later loads retain UI choices."""
    set_command = Mock()
    ds9_adapter.set = set_command
    first, second = tmp_path / 'first.fits', tmp_path / 'second.fits'

    ds9_adapter.load_fits(first)
    assert call('scale', 'mode', 'zscale') in set_command.call_args_list
    assert call('scale', 'asinh') in set_command.call_args_list

    # Represent the settings selected manually in DS9 before pressing Good.
    ds9_adapter.set_scale_function('log')
    ds9_adapter.set_scale_mode('minmax')
    set_command.reset_mock()

    ds9_adapter.load_fits(second)
    ds9_adapter.load_fits(first)

    assert set_command.call_args_list == [
        call('fits', f'"{second.resolve()}"'),
        call('zoom', 'to', 'fit'),
        call('fits', f'"{first.resolve()}"'),
        call('zoom', 'to', 'fit'),
    ]


def test_comparison_load_only_sends_fits_command_and_selects_mask(tmp_path, ds9_adapter):
    ds9_adapter.set = Mock()
    image = tmp_path / 'night with spaces' / 'image.fits'
    ds9_adapter.load_fits(image, configure_display=False, hdu='MASK')
    ds9_adapter.set.assert_called_once_with('fits', f'"{image.resolve()}[MASK]"')
    assert not ds9_adapter._has_loaded_image


def test_load_only_does_not_accept_conflicting_scale_arguments(tmp_path, ds9_adapter):
    ds9_adapter.set = Mock()
    with pytest.raises(ValueError, match='configure_display'):
        ds9_adapter.load_fits(tmp_path / 'image.fits', scale_function='log', configure_display=False)
    ds9_adapter.set.assert_not_called()


@pytest.mark.parametrize('hdu', [-1, True, 'MASK]; frame delete all', 'MASK\n'])
def test_invalid_hdu_selector_is_rejected_before_communication(tmp_path, ds9_adapter, hdu):
    ds9_adapter.set = Mock()
    with pytest.raises(ValueError, match='hdu'):
        ds9_adapter.load_fits(tmp_path / 'image.fits', hdu=hdu)
    ds9_adapter.set.assert_not_called()


def test_reset_scale_option_reapplies_defaults_on_each_load(
    tmp_path: Path,
    ds9_adapter: DS9Adapter,
) -> None:
    """Callers can explicitly request the former reset-on-every-image policy."""
    ds9_adapter.reset_scale_on_load = True
    set_command = Mock()
    ds9_adapter.set = set_command

    ds9_adapter.load_fits(tmp_path / 'first.fits')
    ds9_adapter.load_fits(tmp_path / 'second.fits')

    assert set_command.call_args_list.count(call('scale', 'asinh')) == 2
    assert set_command.call_args_list.count(call('scale', 'mode', 'zscale')) == 2


def test_explicit_stretch_after_initialization_does_not_reset_scale_mode(
    tmp_path: Path,
    ds9_adapter: DS9Adapter,
) -> None:
    """An explicit display argument still applies after the first image."""
    ds9_adapter.set = Mock()
    ds9_adapter.load_fits(tmp_path / 'first.fits')
    ds9_adapter.set.reset_mock()

    ds9_adapter.load_fits(tmp_path / 'second.fits', scale_function='log')

    assert ds9_adapter.set.call_args_list == [
        call('fits', f'"{(tmp_path / "second.fits").resolve()}"'),
        call('scale', 'log'),
        call('zoom', 'to', 'fit'),
    ]


def test_failed_first_display_retries_initialization(
    tmp_path: Path,
    ds9_adapter: DS9Adapter,
) -> None:
    """An unsuccessful display must not suppress initialization on Retry."""
    ds9_adapter.set = Mock(side_effect=DS9CommunicationError('DS9 unavailable'))
    with pytest.raises(DS9CommunicationError):
        ds9_adapter.load_fits(tmp_path / 'first.fits')

    ds9_adapter.set = Mock()
    ds9_adapter.load_fits(tmp_path / 'first.fits')

    assert call('scale', 'asinh') in ds9_adapter.set.call_args_list
    assert call('scale', 'mode', 'zscale') in ds9_adapter.set.call_args_list


@pytest.mark.parametrize(
    ("method_name", "arguments", "expected_call"),
    [
        ("set_zscale", (), call("scale", "mode", "zscale")),
        ("set_cmap", ("heat",), call("cmap", "heat")),
        ("clear_regions", (), call("region", "delete")),
        ("set_scale_mode", ("minmax",), call("scale", "mode", "minmax")),
        ("set_scale_function", ("log",), call("scale", "log")),
        (
            "set_scale_limits",
            (10.5, 42.0),
            call("scale", "limits", "10.5", "42.0"),
        ),
        ("zoom_to_fit", (), call("zoom", "to", "fit")),
    ],
)
def test_command_helpers_delegate_to_set(
    ds9_adapter: DS9Adapter,
    method_name: str,
    arguments: tuple[object, ...],
    expected_call: object,
) -> None:
    """Small convenience methods preserve their intended XPA token order."""

    set_command = Mock()
    ds9_adapter.set = set_command

    getattr(ds9_adapter, method_name)(*arguments)

    assert set_command.call_args == expected_call


@pytest.mark.parametrize("method_name", ["load_region", "save_region"])
def test_region_paths_are_resolved_and_quoted(
    tmp_path: Path,
    ds9_adapter: DS9Adapter,
    method_name: str,
) -> None:
    """Region helpers protect paths containing spaces from DS9 parsing."""

    region_path = tmp_path / "region files" / "source.reg"
    set_command = Mock()
    ds9_adapter.set = set_command

    getattr(ds9_adapter, method_name)(region_path)

    action = "load" if method_name == "load_region" else "save"
    set_command.assert_called_once_with(
        "region",
        action,
        f'"{region_path.resolve()}"',
    )


def test_set_cmap_parameters_delegates_valid_values(
    ds9_adapter: DS9Adapter,
) -> None:
    """Valid contrast and bias values are converted to XPA strings."""

    set_command = Mock()
    ds9_adapter.set = set_command

    ds9_adapter.set_cmap_parameters(contrast=2.5, bias=0.25)

    set_command.assert_called_once_with("cmap", "2.5", "0.25")


@pytest.mark.parametrize(
    ("contrast", "bias", "message"),
    [
        (-0.1, 0.5, "contrast must be between 0 and 10"),
        (10.1, 0.5, "contrast must be between 0 and 10"),
        (1.0, -0.1, "bias must be between 0 and 1"),
        (1.0, 1.1, "bias must be between 0 and 1"),
    ],
)
def test_set_cmap_parameters_rejects_out_of_range_values(
    ds9_adapter: DS9Adapter,
    contrast: float,
    bias: float,
    message: str,
) -> None:
    """Invalid display-control values fail before sending an XPA command."""

    set_command = Mock()
    ds9_adapter.set = set_command

    with pytest.raises(ValueError, match=message):
        ds9_adapter.set_cmap_parameters(contrast=contrast, bias=bias)

    set_command.assert_not_called()
