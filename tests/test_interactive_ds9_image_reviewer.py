"""Unit tests for the Tkinter-controlled DS9 image reviewer."""

from collections.abc import Callable
from pathlib import Path
from typing import cast
from unittest.mock import Mock

import pytest

from photozpy.interactive_ds9 import image_reviewer as reviewer_module
from photozpy.interactive_ds9.adapter import DS9Adapter
from photozpy.interactive_ds9.image_reviewer import ImageReviewer


class FakeVariable:
    """Small replacement for ``tk.StringVar`` used in headless tests."""

    def __init__(self) -> None:
        self.value = ""

    def set(self, value: str) -> None:
        """Record the latest value assigned by the reviewer."""

        self.value = value

    def get(self) -> str:
        """Return the current value through the normal Tk variable API."""

        return self.value


class FakeWidget:
    """Store Tk widget options without requiring a display server."""

    def __init__(self, *args: object, **kwargs: object) -> None:
        self.options = dict(kwargs)
        self.pack_options: dict[str, object] = {}

    def pack(self, **kwargs: object) -> None:
        """Record layout options supplied by the implementation."""

        self.pack_options = dict(kwargs)

    def config(self, **kwargs: object) -> None:
        """Apply state changes in the same style as a Tk widget."""

        self.options.update(kwargs)

    def cget(self, key: str) -> object:
        """Return one configured option through the normal Tk widget API."""

        return self.options[key]


class FakeRoot(FakeWidget):
    """Controllable replacement for the Tk root and its update loop."""

    def __init__(self) -> None:
        super().__init__()
        self.window_title = ""
        self.protocol_handlers: dict[str, Callable[[], None]] = {}
        self.on_update: Callable[[], None] | None = None
        self.quit_called = False
        self.destroy_called = False
        self.update_idletasks_calls = 0
        self.update_calls = 0

    def title(self, value: str) -> None:
        """Record the configured window title."""

        self.window_title = value

    def resizable(self, width: bool, height: bool) -> None:
        """Record the configured resize policy."""

        self.options["resizable"] = (width, height)

    def protocol(self, name: str, callback: Callable[[], None]) -> None:
        """Record a window-manager callback."""

        self.protocol_handlers[name] = callback

    def quit(self) -> None:
        """Record that Tk's event loop was asked to stop."""

        self.quit_called = True

    def destroy(self) -> None:
        """Record that the Tk root was destroyed."""

        self.destroy_called = True

    def update_idletasks(self) -> None:
        """Count idle-task updates made by ``ImageReviewer.run``."""

        self.update_idletasks_calls += 1

    def update(self) -> None:
        """Run an optional callback to control the test event loop."""

        self.update_calls += 1
        if self.on_update is not None:
            callback = self.on_update
            self.on_update = None
            callback()


@pytest.fixture
def fake_tkinter(
    monkeypatch: pytest.MonkeyPatch,
) -> list[tuple[str, str]]:
    """Replace Tk widgets and collect message-box errors for each test."""

    errors: list[tuple[str, str]] = []
    monkeypatch.setattr(reviewer_module.tk, "Tk", FakeRoot)
    monkeypatch.setattr(reviewer_module.tk, "StringVar", FakeVariable)
    monkeypatch.setattr(reviewer_module.tk, "Label", FakeWidget)
    monkeypatch.setattr(reviewer_module.tk, "Frame", FakeWidget)
    monkeypatch.setattr(reviewer_module.tk, "Button", FakeWidget)
    monkeypatch.setattr(
        reviewer_module.messagebox,
        "showerror",
        lambda title, message: errors.append((title, message)),
    )
    # Most tests exercise state transitions, not diagnostic file output.
    monkeypatch.setattr(ImageReviewer, "_debug", lambda self, message: None)
    return errors


def _create_reviewer(
    images: list[Path],
    rejected_dir: Path,
    adapter: DS9Adapter | None = None,
) -> ImageReviewer:
    """Construct a reviewer after the headless Tk fixture is active."""

    if adapter is None:
        adapter = Mock(spec=DS9Adapter)
    return ImageReviewer(images, rejected_dir, adapter=adapter)


def _write_placeholder_fits(path: Path) -> Path:
    """Create a disposable file; unit tests do not parse FITS contents."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"temporary test image")
    return path


def test_initialization_normalizes_paths_without_creating_rejected_dir(
    tmp_path: Path,
    fake_tkinter: list[tuple[str, str]],
) -> None:
    """Opening the reviewer alone has no rejected-directory side effect."""

    image = tmp_path / "night one" / "image.fits"
    rejected_dir = tmp_path / "rejected"

    reviewer = _create_reviewer([image], rejected_dir)
    root = cast(FakeRoot, reviewer.root)

    assert reviewer.images == [image.resolve()]
    assert reviewer.rejected_dir == rejected_dir.resolve()
    assert not rejected_dir.exists()
    assert reviewer.index == 0
    assert reviewer.history == []
    assert reviewer.rejected == []
    assert root.window_title == "Photozpy DS9 Image Review"
    assert "WM_DELETE_WINDOW" in root.protocol_handlers


def test_show_current_loads_image_and_sets_ready_button_states(
    tmp_path: Path,
    fake_tkinter: list[tuple[str, str]],
) -> None:
    """A successful DS9 load enables decisions and disables Retry."""

    image = tmp_path / "image.fits"
    adapter = Mock(spec=DS9Adapter)
    reviewer = _create_reviewer([image], tmp_path / "rejected", adapter)

    reviewer.show_current()

    adapter.load_fits.assert_called_once_with(image.resolve())
    assert reviewer.status.get() == "1/1  image.fits"
    assert reviewer.good_button.cget("state") == reviewer_module.tk.NORMAL
    assert reviewer.reject_button.cget("state") == reviewer_module.tk.NORMAL
    assert reviewer.back_button.cget("state") == reviewer_module.tk.DISABLED
    assert reviewer.retry_button.cget("state") == reviewer_module.tk.DISABLED
    assert fake_tkinter == []


def test_show_current_keeps_position_and_enables_retry_after_ds9_error(
    tmp_path: Path,
    fake_tkinter: list[tuple[str, str]],
) -> None:
    """A display failure cannot accidentally advance or accept an image."""

    image = tmp_path / "image.fits"
    adapter = Mock(spec=DS9Adapter)
    adapter.load_fits.side_effect = RuntimeError("DS9 is unavailable")
    reviewer = _create_reviewer([image], tmp_path / "rejected", adapter)

    reviewer.show_current()

    assert reviewer.index == 0
    assert reviewer.good_button.cget("state") == reviewer_module.tk.DISABLED
    assert reviewer.reject_button.cget("state") == reviewer_module.tk.DISABLED
    assert reviewer.retry_button.cget("state") == reviewer_module.tk.NORMAL
    assert fake_tkinter == [
        (
            "DS9 error",
            f"Could not display:\n{image.resolve()}\n\nDS9 is unavailable",
        )
    ]


def test_mark_good_records_decision_and_advances(
    tmp_path: Path,
    fake_tkinter: list[tuple[str, str]],
) -> None:
    """Good records only review state and leaves the image in place."""

    image = _write_placeholder_fits(tmp_path / "image.fits")
    reviewer = _create_reviewer([image], tmp_path / "rejected")
    reviewer.show_current = Mock()

    reviewer.mark_good()

    assert reviewer.index == 1
    assert reviewer.history == ["good"]
    assert image.exists()
    reviewer.show_current.assert_called_once_with()


def test_mark_reject_creates_directory_lazily_and_moves_image(
    tmp_path: Path,
    fake_tkinter: list[tuple[str, str]],
) -> None:
    """Reject creates its directory on demand and records the new path."""

    image = _write_placeholder_fits(tmp_path / "source" / "image.fits")
    rejected_dir = tmp_path / "rejected"
    reviewer = _create_reviewer([image], rejected_dir)
    reviewer.show_current = Mock()

    reviewer.mark_reject()

    destination = rejected_dir / image.name
    assert not image.exists()
    assert destination.exists()
    assert reviewer.rejected == [destination]
    assert reviewer.history == ["reject"]
    assert reviewer.index == 1
    reviewer.show_current.assert_called_once_with()
    assert fake_tkinter == []


def test_mark_reject_refuses_to_overwrite_existing_destination(
    tmp_path: Path,
    fake_tkinter: list[tuple[str, str]],
) -> None:
    """A filename collision preserves both files and the review position."""

    image = _write_placeholder_fits(tmp_path / "source" / "image.fits")
    destination = _write_placeholder_fits(
        tmp_path / "rejected" / "image.fits"
    )
    reviewer = _create_reviewer([image], destination.parent)
    reviewer.show_current = Mock()

    reviewer.mark_reject()

    assert image.exists()
    assert destination.exists()
    assert reviewer.index == 0
    assert reviewer.history == []
    assert reviewer.rejected == []
    reviewer.show_current.assert_not_called()
    assert fake_tkinter == [
        (
            "Reject error",
            f"Rejected file already exists:\n{destination}",
        )
    ]


def test_mark_reject_preserves_state_when_move_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fake_tkinter: list[tuple[str, str]],
) -> None:
    """A filesystem failure leaves the image available for another action."""

    image = _write_placeholder_fits(tmp_path / "source" / "image.fits")
    rejected_dir = tmp_path / "rejected"
    reviewer = _create_reviewer([image], rejected_dir)
    reviewer.show_current = Mock()

    def fail_move(source: str, destination: str) -> None:
        raise OSError("disk is read-only")

    monkeypatch.setattr(reviewer_module.shutil, "move", fail_move)

    reviewer.mark_reject()

    assert image.exists()
    assert reviewer.index == 0
    assert reviewer.history == []
    assert reviewer.rejected == []
    reviewer.show_current.assert_not_called()
    assert fake_tkinter == [
        (
            "Reject error",
            f"Could not reject:\n{image}\n\ndisk is read-only",
        )
    ]


def test_go_back_undoes_good_decision_without_moving_file(
    tmp_path: Path,
    fake_tkinter: list[tuple[str, str]],
) -> None:
    """Back returns to a previously accepted image and removes its history."""

    image = _write_placeholder_fits(tmp_path / "image.fits")
    reviewer = _create_reviewer([image], tmp_path / "rejected")
    reviewer.index = 1
    reviewer.history = ["good"]
    reviewer.show_current = Mock()

    reviewer.go_back()

    assert reviewer.index == 0
    assert reviewer.history == []
    assert image.exists()
    reviewer.show_current.assert_called_once_with()


def test_go_back_restores_previously_rejected_image(
    tmp_path: Path,
    fake_tkinter: list[tuple[str, str]],
) -> None:
    """Back reverses the file move before revisiting a rejected image."""

    image = _write_placeholder_fits(tmp_path / "source" / "image.fits")
    rejected_dir = tmp_path / "rejected"
    reviewer = _create_reviewer([image], rejected_dir)
    reviewer.show_current = Mock()
    reviewer.mark_reject()
    destination = rejected_dir / image.name
    reviewer.show_current.reset_mock()

    reviewer.go_back()

    assert image.exists()
    assert not destination.exists()
    assert reviewer.index == 0
    assert reviewer.history == []
    assert reviewer.rejected == []
    reviewer.show_current.assert_called_once_with()
    assert fake_tkinter == []


def test_go_back_preserves_state_when_rejected_file_is_missing(
    tmp_path: Path,
    fake_tkinter: list[tuple[str, str]],
) -> None:
    """A missing rejected file is reported without corrupting history."""

    image = tmp_path / "source" / "image.fits"
    rejected_path = tmp_path / "rejected" / image.name
    reviewer = _create_reviewer([image], rejected_path.parent)
    reviewer.index = 1
    reviewer.history = ["reject"]
    reviewer.rejected = [rejected_path]
    reviewer.show_current = Mock()

    reviewer.go_back()

    assert reviewer.index == 1
    assert reviewer.history == ["reject"]
    assert reviewer.rejected == [rejected_path]
    reviewer.show_current.assert_not_called()
    assert fake_tkinter == [
        (
            "Restore error",
            f"Rejected file does not exist:\n{rejected_path}",
        )
    ]


def test_finish_disables_decisions_but_allows_back_after_progress(
    tmp_path: Path,
    fake_tkinter: list[tuple[str, str]],
) -> None:
    """Completion is visible while Back remains available for correction."""

    image = tmp_path / "image.fits"
    reviewer = _create_reviewer([image], tmp_path / "rejected")
    reviewer.index = 1
    reviewer.rejected = [tmp_path / "rejected" / image.name]

    reviewer.finish()

    assert reviewer.status.get() == "Review complete. Rejected: 1"
    assert reviewer.good_button.cget("state") == reviewer_module.tk.DISABLED
    assert reviewer.reject_button.cget("state") == reviewer_module.tk.DISABLED
    assert reviewer.back_button.cget("state") == reviewer_module.tk.NORMAL
    assert reviewer.retry_button.cget("state") == reviewer_module.tk.DISABLED


def test_quit_review_stops_loop_and_destroys_window(
    tmp_path: Path,
    fake_tkinter: list[tuple[str, str]],
) -> None:
    """Quit performs both pieces required for notebook-friendly shutdown."""

    reviewer = _create_reviewer([], tmp_path / "rejected")
    root = cast(FakeRoot, reviewer.root)
    reviewer._is_running = True

    reviewer.quit_review()

    assert reviewer._is_running is False
    assert root.quit_called is True
    assert root.destroy_called is True


def test_run_returns_after_quit_callback_in_explicit_update_loop(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fake_tkinter: list[tuple[str, str]],
) -> None:
    """The custom event loop returns cleanly in the notebook use case."""

    image = tmp_path / "image.fits"
    adapter = Mock(spec=DS9Adapter)
    reviewer = _create_reviewer([image], tmp_path / "rejected", adapter)
    root = cast(FakeRoot, reviewer.root)
    root.on_update = reviewer.quit_review
    monkeypatch.setattr(reviewer_module.time, "sleep", lambda seconds: None)

    result = reviewer.run()

    assert result == []
    adapter.load_fits.assert_called_once_with(image.resolve())
    assert root.update_idletasks_calls == 1
    assert root.update_calls == 1
    assert root.quit_called is True
    assert root.destroy_called is True
    assert reviewer._is_running is False
