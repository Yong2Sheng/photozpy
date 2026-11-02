"""User-paced Before/After/Mask comparisons using the existing review controls."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
import tkinter as tk

from .adapter import DS9Adapter
from .image_reviewer import ImageReviewer


class ImageComparisonReviewer(ImageReviewer):
    """Inspect paired FITS without Good/Reject decisions or file moves.

    Reuse ImageReviewer's queue, Back/Next navigation, error display, and macOS
    process isolation. A dedicated DS9 target keeps quality-review frames intact.
    Frames 1, 2, and 3 display Before, After, and the After CMASK respectively.
    A manual three-column, one-row grid keeps them side by side. The initial
    zoom-to-fit is computed after all three frames have their final tile sizes.
    DS9 reads the CMASK extension directly, without making a temporary image.
    Backups are never deleted. DS9 stays open. No scientific calculation occurs.

    Before and After start with identical stretch/limits and view geometry.
    Adjust any frame in DS9, then click Sync views to copy its zoom/position
    to all three. A science-frame source also synchronizes science display
    settings; a Mask source leaves both science scales unchanged. The mask's
    black/white 0..1 scale is never copied to science images or overwritten.
    Next/Previous explicitly carry the selected frame's current zoom/position
    to the new pair, without re-fitting or re-synchronizing intensity settings.
    Do not rely on DS9's per-file geometry history. Use Sync views to match scales.
    Do not use
    DS9's global scale lock: it would also apply science settings to the mask.
    """

    target_name = 'photozpy-cr-preview'

    def __init__(
        self,
        pairs: Iterable[tuple[str | Path, str | Path]],
        adapter: DS9Adapter | None = None,
        *,
        isolated: bool | None = None,
    ) -> None:
        self.pairs: list[tuple[Path, Path]] = [
            (Path(before).expanduser().resolve(), Path(after).expanduser().resolve())
            for before, after in pairs
        ]
        if not self.pairs:
            raise ValueError('At least one Before/After pair is required.')
        for pair in self.pairs:
            for path in pair:
                if not path.is_file():
                    raise FileNotFoundError(path)
            if pair[0] == pair[1]:
                raise ValueError('Before and After must be separate retained files.')

        self._display_ready = False
        self._loaded = False
        super().__init__(
            [after for _, after in self.pairs], self.pairs[0][0].parent,
            adapter=adapter, isolated=isolated,
        )
        if self._isolated:
            return
        if adapter is None:
            self.adapter = DS9Adapter(software_name=self.target_name)

        self.root.title('Photozpy Cosmic-ray Comparison')
        self.good_button.config(text='Next')
        self.back_button.config(text='Previous')
        # Reuse this button's position, replacing the file-moving action.
        self.reject_button.config(text='Sync views', command=self.sync_views)

    def mark_reject(self) -> None:
        """Comparison sessions never move input files, even via direct calls."""
        raise RuntimeError('Comparison preview cannot reject or move files.')

    def _capture_view(self, *, geometry_only: bool = False) -> dict[str, tuple[str, ...]]:
        """Read the selected frame; Mask supplies geometry, never science scale."""
        source_frame: str = self.adapter.get('frame')
        if source_frame not in ('1', '2', '3'):
            raise ValueError('Select Before, After, or Mask in DS9 before synchronizing.')

        queries = {
            'zoom': ('zoom',), 'pan': ('pan', 'image'),
        }
        if source_frame != '3' and not geometry_only:
            queries.update({
                'stretch': ('scale',), 'limits': ('scale', 'limits'),
                'cmap': ('cmap',), 'contrast': ('cmap', 'value'),
                'invert': ('cmap', 'invert'),
                'log_exponent': ('scale', 'log', 'exp'),
            })

        view = {name: tuple(self.adapter.get(*query).split()) for name, query in queries.items()}
        view['source_frame'] = (source_frame,)
        return view

    def _apply_view(self, view: dict[str, tuple[str, ...]]) -> None:
        """Match geometry from any frame; science scales only from science."""
        source_frame = view['source_frame'][0]
        for frame in ('1', '2', '3'):
            self.adapter.set('frame', frame)
            self.adapter.set('zoom', 'to', *view['zoom'])
            self.adapter.set('pan', 'to', *view['pan'], 'image')
            if frame != '3' and source_frame != '3':
                self.adapter.set('scale', 'log', 'exp', *view['log_exponent'])
                self.adapter.set('scale', *view['stretch'])
                self.adapter.set('scale', 'limits', *view['limits'])
                self.adapter.set('cmap', *view['cmap'])
                self.adapter.set('cmap', *view['contrast'])
                self.adapter.set('cmap', 'invert', *view['invert'])
        self.adapter.set('frame', source_frame)

    def sync_views(self) -> None:
        """Use the selected Before/After/Mask frame as the view reference.

        Mask synchronizes zoom and image-coordinate position only. Its binary
        scale must not replace the intensity scaling of the science frames.
        Keep the source frame selected so users can continue locating features.
        """
        if not self._loaded:
            return
        try:
            self._apply_view(self._capture_view())
            self.error_label.config(text='')
        except Exception as error:
            self._show_error('Could not synchronize views', str(error))

    def show_current(self) -> None:
        """Load only the pair selected by the user, with a separate binary mask."""
        self._require_local_gui()
        if self.index >= len(self.pairs):
            self.finish()
            return

        before, after = self.pairs[self.index]
        self.status.set(
            f'{self.index + 1}/{len(self.pairs)}: {after.name}\n'
            'Frame 1: Before | Frame 2: After | Frame 3: Mask'
        )
        self.good_button.config(state=tk.DISABLED)
        try:
            if not self._display_ready:
                self.adapter.ensure_running()
                self.adapter.set('tile', 'yes')
                self.adapter.set('tile', 'mode', 'grid')
                self.adapter.set('tile', 'grid', 'mode', 'manual')
                self.adapter.set('tile', 'grid', 'direction', 'x')
                self.adapter.set('tile', 'grid', 'layout', '3', '1')
                self.adapter.set('scale', 'lock', 'no')
                self.adapter.set('scale', 'lock', 'limits', 'no')
                self.adapter.set('lock', 'colorbar', 'no')
                self.adapter.set('preserve', 'pan', 'yes')
                self._display_ready = True

            # Read the live outgoing view before any frame changes or loads.
            # DS9 can restore older per-file geometry despite preserve pan.
            geometry = self._capture_view(geometry_only=True) if self._loaded else None

            # Restore geometry immediately after each load, rather than load
            # all three first and then visibly resize the complete layout.
            for frame, path, hdu in (('1', before, None), ('2', after, None), ('3', after, 'CMASK')):
                self.adapter.set('frame', frame)
                self.adapter.load_fits(path, configure_display=False, hdu=hdu)
                if not self._loaded:
                    self.adapter.verify_image_loaded()
                if geometry is not None:
                    self.adapter.set('zoom', 'to', *geometry['zoom'])
                    self.adapter.set('pan', 'to', *geometry['pan'], 'image')
                if frame == '3' and not self._loaded:
                    self.adapter.set_scale_function('linear')
                    self.adapter.set_scale_limits(0., 1.)
                    self.adapter.set_cmap('gray')
                    self.adapter.set_cmap_parameters(1., 0.5)
                    self.adapter.set('cmap', 'invert', 'no')
            if geometry is not None:
                self.adapter.set('frame', geometry['source_frame'][0])

            if not self._loaded:
                # Adding frames changes tile dimensions. Fit only once the
                # complete layout exists, then copy this view to all frames.
                self.adapter.set('frame', '1')
                self.adapter.set_scale_mode('zscale')
                self.adapter.set_scale_function('asinh')
                self.adapter.set('update', 'now')
                self.adapter.zoom_to_fit()
                view = self._capture_view()
                self._apply_view(view)

            self._loaded = True
            self.error_label.config(text='')
            self.good_button.config(state=tk.NORMAL)
            self.retry_button.config(state=tk.DISABLED)
        except Exception as error:
            self._show_error('Could not load comparison', str(error))
            self.retry_button.config(state=tk.NORMAL)
            self._display_ready = False

        self.back_button.config(state=tk.NORMAL if self.index else tk.DISABLED)
        self.reject_button.config(state=tk.NORMAL if self._loaded else tk.DISABLED)
        self._save_state()

    def finish(self) -> None:
        """Keep Previous available after reaching the last pair."""
        super().finish()
        self.status.set('Comparison complete. Backups retained; use Previous or Quit.')
        self.reject_button.config(state=tk.NORMAL if self._loaded else tk.DISABLED)
        self._save_state()

    def run(self) -> list[Path]:
        """Wait for user-paced comparison; no decision moves or deletes files."""
        if self._closed:
            raise RuntimeError('This comparison is closed; create another session.')
        if self._isolated:
            from ._review_process import run_isolated

            config = {
                'software_name': self.target_name if self.adapter is None else self.adapter.software_name,
                'reset_scale_on_load': False,
            }
            try:
                run_isolated(
                    self.images, self.rejected_dir, config, self._restore_state,
                    comparison_pairs=self.pairs,
                )
            finally:
                self._closed = True
            return []

        return super().run()
