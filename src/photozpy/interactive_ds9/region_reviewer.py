"""Edit paired photometry regions using the existing DS9 preview navigation."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
import tkinter as tk

from .adapter import DS9Adapter
from .image_previewer import ImagePreviewer
from ..photometry.regions import PhotometryRegions


class RegionReviewer(ImagePreviewer):
    """Review one or more images with per-image aperture/background files.

    Save & Next validates and persists the current pair before advancing. The
    last image is saved by the same button. Previous reloads saved files and
    discards unconfirmed current edits; Quit never saves unconfirmed edits.
    Initial catalog positions remain fixed blue crosses. macOS uses the
    existing disposable GUI worker, keeping Tk out of the notebook kernel.
    No science image is moved or written, and no source magnitude is changed.
    Same-source/group filter transitions preserve zoom and image-coordinate
    pan, including Previous. First display and target/field changes zoom to fit.
    """

    target_name = "photozpy-region-review"

    def __init__(
        self,
        sessions: Iterable[PhotometryRegions],
        adapter: DS9Adapter | None = None,
        *,
        isolated: bool | None = None,
    ) -> None:
        self.sessions = tuple(sessions)
        self._display_group: tuple | None = None
        if not self.sessions:
            raise ValueError("At least one region session is required.")
        super().__init__(
            [session.image_path for session in self.sessions],
            adapter=adapter,
            isolated=isolated,
        )
        if not self._isolated:
            self.root.title("Photozpy Aperture / Annulus Review")
            self.good_button.config(text="Save & Next")

    def _load_image(self, path: Path) -> None:
        """Preserve the user's view only within the same field/source group."""
        records = self.sessions[self.index].records
        group = (
            records[0]["field_name"],
            tuple(sorted((row["source_role"], row["source_id"]) for row in records)),
        )
        try:
            self.adapter.load_fits(path, preserve_view=group == self._display_group)
        except Exception:
            self._display_group = None
            raise
        self._display_group = group

    def show_current(self) -> None:
        """Display the image, both editable files, and immutable initial crosses."""
        super().show_current()
        if self.index >= len(self.sessions) or self.error_message:
            return
        try:
            session = self.sessions[self.index]
            self.adapter.clear_regions()
            self.adapter.load_region(session.aperture_path)
            self.adapter.add_regions(session.background_display_regions())
            self.adapter.add_regions(session.initial_regions())
        except Exception as error:
            self.good_button.config(state=tk.DISABLED)
            self.retry_button.config(state=tk.NORMAL)
            self._show_error("Region load error", str(error))

    def mark_good(self) -> None:
        """Read DS9 image-coordinate edits, validate/save, and only then advance."""
        if self.index >= len(self.sessions) or self.error_message:
            return
        try:
            text = self.adapter.get("regions", "-system", "image", "-format", "ds9")
            image_name = self.current_image.name
            review_warnings = self.sessions[self.index].save_review(text)
        except Exception as error:
            self._show_error("Region save error", str(error))
            # Keep editing and allow another save attempt without reloading.
            self.error_message = ""
            return
        self.error_label.config(text="")
        super().mark_good()
        # show_current/finish may clear the error label; show saved warnings last.
        # Do not use _show_error: a warning must not block subsequent navigation.
        if review_warnings:
            self.error_label.config(
                text=(
                    f"Warning — saved {image_name}:\n"
                    + "\n".join(review_warnings)
                    + (f"\n{self.error_message}" if self.error_message else "")
                ),
                fg="firebrick",
            )

    def run(self) -> list:
        """Wait for review; saved manifests remain usable after interruption."""
        if self._closed:
            raise RuntimeError("This region review is closed; create another session.")
        if not self._isolated:
            return super().run()
        from ._review_process import run_isolated

        config = {
            "software_name": (
                self.target_name if self.adapter is None else self.adapter.software_name
            ),
            "reset_scale_on_load": False,
        }
        try:
            run_isolated(
                self.images,
                self.rejected_dir,
                config,
                self._restore_state,
                region_sessions=[
                    {
                        "root": str(session.root),
                        "manifest": session.manifest_path.as_posix(),
                    }
                    for session in self.sessions
                ],
            )
        finally:
            self._closed = True
            for session in self.sessions:
                session._state = PhotometryRegions.load(
                    session.root, session.manifest_path
                )._state
        return []
