"""Browse completed FITS products without file-moving quality decisions."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
import tkinter as tk

from .adapter import DS9Adapter
from .image_reviewer import ImageReviewer


class ImagePreviewer(ImageReviewer):
    """Reuse the reviewer display/navigation with no Reject operation.

    FITS image values are displayed directly, without modifying files or
    applying their MASK/CMASK extensions to the science display. DS9 starts
    automatically in a dedicated window and remains open after Quit.
    macOS uses the existing disposable GUI process to protect notebook Tk state.
    """

    target_name = 'photozpy-master-preview'

    def __init__(
        self, images: Iterable[str | Path], adapter: DS9Adapter | None = None,
        *, isolated: bool | None = None,
    ) -> None:
        selected: list[Path] = [
            Path(path).expanduser().resolve() for path in images
            if not Path(path).name.startswith('.')
        ]
        if not selected:
            raise ValueError('At least one FITS product is required for preview.')
        for path in selected:
            if not path.is_file():
                raise FileNotFoundError(path)

        # The inherited queue needs a directory, but preview never uses it to
        # create folders or move files. mark_reject is blocked independently.
        super().__init__(selected, selected[0].parent, adapter=adapter, isolated=isolated)
        if self._isolated:
            return
        if adapter is None:
            self.adapter = DS9Adapter(software_name=self.target_name)
        self.root.title('Photozpy Master Image Preview')
        self.good_button.config(text='Next')
        self.back_button.config(text='Previous')
        self.reject_button.pack_forget()

    def mark_reject(self) -> None:
        """Disallow file moves, including direct calls outside the GUI."""
        raise RuntimeError('Product preview cannot reject or move files.')

    def show_current(self) -> None:
        """Use the original display/error handling, keeping Reject disabled."""
        super().show_current()
        self.reject_button.config(state=tk.DISABLED)

    def finish(self) -> None:
        """Keep Previous available at the end until the user chooses Quit."""
        super().finish()
        self.status.set('Preview complete. Use Previous or Quit.')

    def run(self) -> list[Path]:
        """Wait for browsing to finish; no FITS files are written or moved."""
        if self._closed:
            raise RuntimeError('This preview is closed; create another session.')
        if self._isolated:
            from ._review_process import run_isolated

            config = {
                'software_name': self.target_name if self.adapter is None else self.adapter.software_name,
                'reset_scale_on_load': False,
            }
            try:
                run_isolated(
                    self.images, self.rejected_dir, config, self._restore_state,
                    preview_only=True,
                )
            finally:
                self._closed = True
            return []
        return super().run()
