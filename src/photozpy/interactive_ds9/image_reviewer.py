"""Tkinter-based DS9 image quality reviewer.

The reviewer displays FITS files in DS9 one at a time and keeps the workflow
controls in a small Tkinter window.  It is intentionally narrow: DS9 handles
image display and manual visual inspection, while this module manages the file
queue and moves rejected files to a user-provided directory.
"""

from datetime import datetime
import json
from pathlib import Path
import shutil
import sys
import time
import tkinter as tk
from typing import Iterable

from .adapter import DS9Adapter


class ImageReviewer:
    """Review FITS images in DS9 and move rejected files to a directory.

    The GUI exposes simple ``Good``, ``Reject``, ``Back``, ``Retry``, and
    ``Quit`` actions.  Rejected images are moved, not deleted, so a mistaken
    rejection can be undone with ``Back`` while the session is still active.

    Moving files does not update an existing ImageFileCollection or other
    collection snapshot. Rebuild the collection after the review completes.
    The first display connects to DS9 or starts it automatically. Subsequent
    images reuse that instance; Retry rechecks readiness after a load failure.
    DS9 remains open when the controls are closed.

    Parameters
    ----------
    images
        Iterable of FITS paths to inspect in order.
    rejected_dir
        Directory where rejected images are moved.  The directory is created
        lazily on the first rejection.
    adapter
        Optional DS9 adapter.  Supplying one is mainly useful for tests or for
        targeting a non-default DS9/XPA access point.
    isolated
        Run the GUI in a separate Python process. Defaults to True on macOS
        for the standard DS9 adapter, keeping Tk out of the notebook kernel.
        Custom adapters default to in-process operation. GUI widget attributes
        are available only in-process. Hidden filenames are excluded.
    """

    def __init__(
        self,
        images: Iterable[str | Path],
        rejected_dir: str | Path,
        adapter: DS9Adapter | None = None,
        *,
        isolated: bool | None = None,
    ) -> None:
        """Initialize image paths, review state, and the Tk control window."""

        self.images: list[Path] = [
            Path(path).expanduser().resolve()
            for path in images
            if not Path(path).name.startswith(".")
        ]

        self.rejected_dir: Path = Path(rejected_dir).expanduser().resolve()

        # Subclasses may carry custom behavior that cannot be sent to a worker.
        adapter_type = type(adapter)
        standard_adapter = adapter is None or adapter_type is DS9Adapter
        self._isolated = (
            sys.platform == "darwin" and standard_adapter
            if isolated is None else isolated
        )
        if self._isolated and not standard_adapter:
            raise TypeError("Custom adapters require isolated=False.")
        self.adapter = adapter
        self._ds9_ready: bool = False
        self._display_verified: bool = False

        self.index: int = 0
        self.rejected: list[Path] = []
        self.history: list[str] = []
        self._is_running: bool = False
        self._closed = False
        self._state_path: Path | None = None
        self.error_message: str = ""
        self.debug_log_path: Path = Path("ds9_reviewer_debug.log").resolve()
        self._debug("--- new ImageReviewer session ---")
        self._debug("init: created ImageReviewer")

        if self._isolated:
            return

        self.adapter = adapter if adapter is not None else DS9Adapter()

        self.root: tk.Tk = tk.Tk()
        self.root.title("Photozpy DS9 Image Review")
        self.root.resizable(False, False)
        self.root.protocol("WM_DELETE_WINDOW", self.quit_review)
        self._debug("init: created Tk root")

        self.status: tk.StringVar = tk.StringVar(master=self.root)
        tk.Label(
            self.root,
            textvariable=self.status,
            width=60,
            anchor="w",
        ).pack(padx=12, pady=(12, 8))

        self.error_label = tk.Label(
            self.root, text="", wraplength=600, justify=tk.LEFT,
            anchor="w", fg="firebrick",
        )
        self.error_label.pack(padx=12, pady=(0, 8))

        button_frame: tk.Frame = tk.Frame(self.root)
        button_frame.pack(padx=12, pady=(0, 12))

        self.good_button: tk.Button = tk.Button(
            button_frame,
            text="Good",
            width=12,
            command=self.mark_good,
        )
        self.good_button.pack(side=tk.LEFT, padx=4)

        self.reject_button: tk.Button = tk.Button(
            button_frame,
            text="Reject",
            width=12,
            command=self.mark_reject,
        )
        self.reject_button.pack(side=tk.LEFT, padx=4)

        self.back_button: tk.Button = tk.Button(
            button_frame,
            text="Back",
            width=12,
            command=self.go_back,
            state=tk.DISABLED,
        )
        self.back_button.pack(side=tk.LEFT, padx=4)

        self.retry_button: tk.Button = tk.Button(
            button_frame,
            text="Retry",
            width=12,
            command=self.show_current,
            state=tk.DISABLED,
        )
        self.retry_button.pack(side=tk.LEFT, padx=4)

        self.quit_button: tk.Button = tk.Button(
            button_frame,
            text="Quit",
            width=12,
            command=self.quit_review,
        )
        self.quit_button.pack(side=tk.LEFT, padx=4)

    def _debug(self, message: str) -> None:
        """Append a debug message to the reviewer log file.

        Debug logging must never interrupt the visual review workflow, so file
        system errors are intentionally ignored.
        """

        try:
            timestamp = datetime.now().isoformat(timespec="milliseconds")
            with self.debug_log_path.open("a", encoding="utf-8") as file:
                file.write(f"{timestamp} {message}\n")
                file.flush()
        except OSError:
            pass

    def _show_error(self, title: str, message: str) -> None:
        """Report errors in the control window, without a native modal dialog."""
        self.error_message = f"{title}: {message}"
        self.error_label.config(text=self.error_message)
        self._debug(self.error_message)

    def _require_local_gui(self) -> None:
        if self._isolated:
            raise RuntimeError(
                "Use the buttons in the isolated review window; "
                "notebook-side GUI actions are unavailable."
            )

    def _save_state(self) -> None:
        """Checkpoint decisions before the next external display operation."""
        if self._state_path is None:
            return
        state = {
            "index": self.index,
            "history": self.history,
            "rejected": [str(path) for path in self.rejected],
        }
        temporary = self._state_path.with_suffix(".tmp")
        temporary.write_text(json.dumps(state), encoding="utf-8")
        temporary.replace(self._state_path)

    def _restore_state(self, state: dict) -> None:
        self.index = state["index"]
        self.history = list(state["history"])
        self.rejected = [Path(path) for path in state["rejected"]]

    @property
    def current_image(self) -> Path:
        """Return the image at the current review position.

        The caller is responsible for ensuring ``self.index`` is within the
        image list bounds before accessing this property.
        """

        return self.images[self.index]

    def show_current(self) -> None:
        """Display the current image in DS9 and update button states.

        If DS9 cannot display the current image, decision buttons are disabled
        and ``Retry`` remains available so transient XPA failures can be tried
        again without advancing the queue.
        """

        self._require_local_gui()
        self._debug(f"show_current: entered index={self.index}")
        self._save_state()
        self.error_message = ""
        self.error_label.config(text="")

        if self.index >= len(self.images):
            self._debug("show_current: index past end, finishing")
            self.finish()
            return

        path = self.current_image
        self._debug(f"show_current: loading {path}")
        self.status.set(
            f"{self.index + 1}/{len(self.images)}  {path.name}"
        )

        try:
            if not self._ds9_ready:
                self.adapter.ensure_running()
                self._ds9_ready = True
                self._display_verified = False
            self.adapter.load_fits(path)
            if not self._display_verified:
                self.adapter.verify_image_loaded()
                self._display_verified = True
        except Exception as error:
            self._ds9_ready = False
            self._display_verified = False
            self._debug(f"show_current: load failed: {error}")
            self.good_button.config(state=tk.DISABLED)
            self.reject_button.config(state=tk.DISABLED)
            self.back_button.config(
                state=tk.NORMAL if self.index > 0 else tk.DISABLED
            )
            self.retry_button.config(state=tk.NORMAL)

            self._show_error(
                "DS9 error",
                f"Could not display:\n{path}\n\n{error}",
            )
            return

        self.good_button.config(state=tk.NORMAL)
        self.reject_button.config(state=tk.NORMAL)
        self.back_button.config(
            state=tk.NORMAL if self.index > 0 else tk.DISABLED
        )
        self.retry_button.config(state=tk.DISABLED)
        self._debug("show_current: load succeeded")

    def mark_good(self) -> None:
        """Accept the current image and advance to the next image."""

        self._require_local_gui()
        self._debug(f"mark_good: index={self.index}")
        self.history.append("good")
        self.index += 1
        self.show_current()

    def mark_reject(self) -> None:
        """Move the current image to the rejected directory and advance.

        The rejected directory is created only when it is first needed.  This
        avoids leaving an empty directory behind when a user opens the reviewer
        but does not reject any images.
        """

        self._require_local_gui()
        source: Path = self.current_image
        self._debug(f"mark_reject: source={source}")

        try:
            self.rejected_dir.mkdir(parents=True, exist_ok=True)
        except OSError as error:
            self._debug(f"mark_reject: mkdir failed: {error}")
            self._show_error(
                "Reject error",
                f"Could not create rejected directory:\n"
                f"{self.rejected_dir}\n\n{error}",
            )
            return

        destination: Path = self.rejected_dir / source.name

        if destination.exists():
            self._debug(f"mark_reject: destination exists: {destination}")
            self._show_error(
                "Reject error",
                f"Rejected file already exists:\n{destination}",
            )
            return

        try:
            shutil.move(str(source), str(destination))
        except OSError as error:
            self._debug(f"mark_reject: move failed: {error}")
            self._show_error(
                "Reject error",
                f"Could not reject:\n{source}\n\n{error}",
            )
            return

        self._debug(f"mark_reject: moved to {destination}")
        self.rejected.append(destination)
        self.history.append("reject")
        self.index += 1
        self.show_current()

    def go_back(self) -> None:
        """Undo the previous decision and display the previous image.

        This method keeps the history deliberately simple.  A ``"good"`` entry
        only moves the index backward.  A ``"reject"`` entry also moves the
        rejected file back to its original path before changing the index.
        """

        self._require_local_gui()
        self._debug(f"go_back: entered index={self.index}")

        if self.index == 0 or not self.history:
            self._debug("go_back: nothing to undo")
            return

        previous_index: int = self.index - 1
        previous_decision: str = self.history[-1]
        self._debug(
            f"go_back: previous_index={previous_index}, "
            f"previous_decision={previous_decision}"
        )

        if previous_decision == "reject":
            source: Path = self.images[previous_index]
            rejected_path: Path = self.rejected_dir / source.name

            if source.exists():
                self._debug(f"go_back: source already exists: {source}")
                self._show_error(
                    "Restore error",
                    f"Original path already exists:\n{source}",
                )
                return

            if not rejected_path.exists():
                self._debug(
                    f"go_back: rejected file missing: {rejected_path}"
                )
                self._show_error(
                    "Restore error",
                    f"Rejected file does not exist:\n{rejected_path}",
                )
                return

            # Undoing a rejection restores the file before changing review state.
            try:
                shutil.move(str(rejected_path), str(source))
            except OSError as error:
                self._debug(f"go_back: restore failed: {error}")
                self._show_error(
                    "Restore error",
                    f"Could not restore:\n{rejected_path}\n\n{error}",
                )
                return

            self._debug(f"go_back: restored {rejected_path} to {source}")
            self.rejected.remove(rejected_path)

        self.history.pop()
        self.index = previous_index
        self.show_current()

    def quit_review(self) -> None:
        """Close the review window and stop the reviewer event loop."""

        self._require_local_gui()
        self._debug("quit_review: entered")
        self._is_running = False
        if self._closed:
            return
        self._closed = True

        try:
            self._debug("quit_review: before quit")
            self.root.quit()
            self._debug("quit_review: after quit")
        except tk.TclError as error:
            self._debug(f"quit_review: quit TclError: {error}")

        try:
            self._debug("quit_review: before destroy")
            self.root.destroy()
            self._debug("quit_review: after destroy")
        except tk.TclError as error:
            self._debug(f"quit_review: TclError: {error}")
            pass

    def finish(self) -> None:
        """Mark the review as complete and disable decision buttons."""

        self._require_local_gui()
        self._debug("finish: entered")
        self.status.set(
            f"Review complete. Rejected: {len(self.rejected)}"
        )
        self.good_button.config(state=tk.DISABLED)
        self.reject_button.config(state=tk.DISABLED)
        self.back_button.config(
            state=tk.NORMAL if self.index > 0 else tk.DISABLED
        )
        self.retry_button.config(state=tk.DISABLED)

    def run(self) -> list[Path]:
        """Run the GUI event loop and return successfully rejected paths.

        macOS uses a disposable GUI process by default. Other platforms and
        custom in-process adapters use the explicit Tk update loop. Interrupts
        close the GUI; completed file decisions are not rolled back.
        """

        if self._closed:
            raise RuntimeError("This review session is closed; create a new reviewer.")
        if self._isolated:
            from ._review_process import run_isolated

            config = None if self.adapter is None else {
                "software_name": self.adapter.software_name,
                "reset_scale_on_load": self.adapter.reset_scale_on_load,
            }
            self._is_running = True
            try:
                run_isolated(
                    self.images, self.rejected_dir, config, self._restore_state,
                )
            finally:
                self._is_running = False
                self._closed = True
            return self.rejected

        self._debug("run: entered")
        self._is_running = True

        try:
            # Map the controls before attempting external DS9 communication.
            self.root.update_idletasks()
            self.show_current()
            self._debug("run: before update loop")
            while self._is_running:
                try:
                    self.root.update_idletasks()
                    self.root.update()
                except tk.TclError as error:
                    self._debug(f"run: update loop TclError: {error}")
                    break
                if self._is_running:
                    time.sleep(0.01)
        finally:
            self.quit_review()
        self._debug("run: after update loop")

        self._debug("run: returning")
        return self.rejected
