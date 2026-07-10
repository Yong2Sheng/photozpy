"""Tkinter-based DS9 image quality reviewer.

The reviewer displays FITS files in DS9 one at a time and keeps the workflow
controls in a small Tkinter window.  It is intentionally narrow: DS9 handles
image display and manual visual inspection, while this module manages the file
queue and moves rejected files to a user-provided directory.
"""

from datetime import datetime
from pathlib import Path
import shutil
import time
import tkinter as tk
from tkinter import messagebox
from typing import Iterable

from .adapter import DS9Adapter


class ImageReviewer:
    """Review FITS images in DS9 and move rejected files to a directory.

    The GUI exposes simple ``Good``, ``Reject``, ``Back``, ``Retry``, and
    ``Quit`` actions.  Rejected images are moved, not deleted, so a mistaken
    rejection can be undone with ``Back`` while the session is still active.

    Moving files does not update an existing ImageFileCollection or other
    collection snapshot. Rebuild the collection after the review completes.

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
    """

    def __init__(
        self,
        images: Iterable[str | Path],
        rejected_dir: str | Path,
        adapter: DS9Adapter | None = None,
    ) -> None:
        """Initialize image paths, review state, and the Tk control window."""

        self.images: list[Path] = [
            Path(path).expanduser().resolve()
            for path in images
        ]

        self.rejected_dir: Path = Path(rejected_dir).expanduser().resolve()

        self.adapter: DS9Adapter = adapter or DS9Adapter()

        self.index: int = 0
        self.rejected: list[Path] = []
        self.history: list[str] = []
        self._is_running: bool = False
        self.debug_log_path: Path = Path("ds9_reviewer_debug.log").resolve()
        self._debug("--- new ImageReviewer session ---")
        self._debug("init: created ImageReviewer")

        self.root: tk.Tk = tk.Tk()
        self.root.title("Photozpy DS9 Image Review")
        self.root.resizable(False, False)
        self.root.protocol("WM_DELETE_WINDOW", self.quit_review)
        self._debug("init: created Tk root")

        self.status: tk.StringVar = tk.StringVar()
        tk.Label(
            self.root,
            textvariable=self.status,
            width=60,
            anchor="w",
        ).pack(padx=12, pady=(12, 8))

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

        self._debug(f"show_current: entered index={self.index}")

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
            self.adapter.load_fits(path)
        except Exception as error:
            self._debug(f"show_current: load failed: {error}")
            self.good_button.config(state=tk.DISABLED)
            self.reject_button.config(state=tk.DISABLED)
            self.back_button.config(
                state=tk.NORMAL if self.index > 0 else tk.DISABLED
            )
            self.retry_button.config(state=tk.NORMAL)

            messagebox.showerror(
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

        source: Path = self.current_image
        self._debug(f"mark_reject: source={source}")

        try:
            self.rejected_dir.mkdir(parents=True, exist_ok=True)
        except OSError as error:
            self._debug(f"mark_reject: mkdir failed: {error}")
            messagebox.showerror(
                "Reject error",
                f"Could not create rejected directory:\n"
                f"{self.rejected_dir}\n\n{error}",
            )
            return

        destination: Path = self.rejected_dir / source.name

        if destination.exists():
            self._debug(f"mark_reject: destination exists: {destination}")
            messagebox.showerror(
                "Reject error",
                f"Rejected file already exists:\n{destination}",
            )
            return

        try:
            shutil.move(str(source), str(destination))
        except OSError as error:
            self._debug(f"mark_reject: move failed: {error}")
            messagebox.showerror(
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
                messagebox.showerror(
                    "Restore error",
                    f"Original path already exists:\n{source}",
                )
                return

            if not rejected_path.exists():
                self._debug(
                    f"go_back: rejected file missing: {rejected_path}"
                )
                messagebox.showerror(
                    "Restore error",
                    f"Rejected file does not exist:\n{rejected_path}",
                )
                return

            # Undoing a rejection restores the file before changing review state.
            try:
                shutil.move(str(rejected_path), str(source))
            except OSError as error:
                self._debug(f"go_back: restore failed: {error}")
                messagebox.showerror(
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

        self._debug("quit_review: entered")
        self._is_running = False

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

        The loop uses explicit ``update()`` calls instead of Tk's blocking
        ``mainloop()``.  This is more stable when the reviewer is launched from
        a Jupyter notebook, where Tk's native event loop can leave the cell in a
        stuck state even after the window is destroyed.
        """

        self._debug("run: entered")
        self._is_running = True

        if not self.images:
            self._debug("run: no images")
            self.finish()
        else:
            self.show_current()

        self._debug("run: before update loop")
        while self._is_running:
            try:
                # Keep control in Python instead of handing the whole thread to
                # Tk.  This makes notebook shutdown behavior much easier to
                # control and diagnose.
                self.root.update_idletasks()
                self.root.update()
            except tk.TclError as error:
                self._debug(f"run: update loop TclError: {error}")
                self._is_running = False
                break
            time.sleep(0.01)
        self._debug("run: after update loop")

        self._debug("run: returning")
        return self.rejected
