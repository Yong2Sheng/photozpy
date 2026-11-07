"""Stage exposure-quality decisions while comparing Before/After/Mask in DS9."""

from __future__ import annotations

import tkinter as tk
import time
from collections.abc import Iterable
from pathlib import Path

from .adapter import DS9Adapter
from .comparison_reviewer import ImageComparisonReviewer
from .image_reviewer import ImageReviewer


class CosmicRayReviewer(ImageComparisonReviewer):
    """Stage Reject decisions without moving files in the GUI process.

    Previous undoes the previous Keep/Reject decision. Normal Quit returns
    the staged cleaned-image paths to Artifacts, which applies the decisions.
    Closing early leaves unvisited images unchanged. An interrupted or failed
    GUI does not apply staged decisions. No scientific calculation occurs.
    """

    target_name = 'photozpy-cr-review'

    def __init__(
        self, pairs: Iterable[tuple[str | Path, str | Path]],
        adapter: DS9Adapter | None = None, *, isolated: bool | None = None,
    ) -> None:
        super().__init__(pairs, adapter=adapter, isolated=isolated)
        if self._isolated:
            return
        self.root.title('Photozpy Cosmic-ray Quality Review')
        self.good_button.config(text='Keep')
        self.reject_button.config(text='Reject', command=self.mark_reject)
        self.sync_button = tk.Button(self.root, text='Sync views', command=self.sync_views)
        self.sync_button.pack(pady=(0, 8))
        tk.Label(
            self.root, wraplength=600,
            text='Reject stages removal of the original exposure. On Quit, the raw '
            'file moves to rejected/ and derived files are deleted. Previous undoes '
            'pending decisions; unvisited images are kept.',
        ).pack(padx=12, pady=(0, 8))

    def mark_good(self) -> None:
        """Keep the current exposure only when its comparison is ready."""
        self._require_local_gui()
        if self.index < len(self.images) and self.good_button.cget('state') == tk.NORMAL:
            super().mark_good()

    def show_current(self) -> None:
        """Allow quality decisions only after a successful comparison load."""
        super().show_current()
        self.reject_button.config(state=self.good_button.cget('state'))

    def mark_reject(self) -> None:
        """Stage rejection of the current original exposure and advance."""
        self._require_local_gui()
        if self.index >= len(self.images) or self.good_button.cget('state') != tk.NORMAL:
            return
        self.rejected.append(self.current_image)
        self.history.append('reject')
        self.index += 1
        self.show_current()

    def go_back(self) -> None:
        """Undo a staged decision; files have not been moved or deleted."""
        self._require_local_gui()
        if not self.history:
            return
        self.index -= 1
        if self.history.pop() == 'reject':
            self.rejected.remove(self.current_image)
        self.show_current()

    def finish(self) -> None:
        """Keep staged decisions reversible until normal Quit applies them."""
        ImageReviewer.finish(self)
        self.status.set(
            f'Review complete. Pending rejects: {len(self.rejected)}. '
            'Previous undoes; Quit applies decisions.'
        )
        self._save_state()

    def run(self) -> list[Path]:
        """Return staged decisions after normal completion, without file writes."""
        if not self._isolated:
            if self._closed:
                raise RuntimeError('This review is closed; create another session.')
            self._is_running = True
            try:
                self.root.update_idletasks()
                self.show_current()
                # Propagate event-loop failures instead of treating them as Quit.
                while self._is_running:
                    self.root.update_idletasks()
                    self.root.update()
                    if self._is_running:
                        time.sleep(0.01)
            finally:
                self.quit_review()
            return self.rejected.copy()
        if self._closed:
            raise RuntimeError('This review is closed; create another session.')
        from ._review_process import run_isolated

        config = {
            'software_name': self.target_name if self.adapter is None else self.adapter.software_name,
            'reset_scale_on_load': False,
        }
        try:
            run_isolated(
                self.images, self.rejected_dir, config, self._restore_state,
                comparison_pairs=self.pairs, cosmic_ray_review=True,
            )
        finally:
            self._closed = True
        return self.rejected.copy()
