"""Headless comparison checks; never launch DS9 or touch observing data."""

from pathlib import Path

from astropy.io import fits
from astropy.nddata import CCDData
import astropy.units as u
import numpy as np
import pytest

from photozpy.interactive_ds9.comparison_reviewer import ImageComparisonReviewer
from test_interactive_ds9_image_reviewer import fake_tkinter as fake_tkinter


class FakeDS9:
    """Record per-frame display state independently of image loading."""

    def __init__(self):
        self.frame = '1'
        self.frames = {}
        self.loaded = []
        self.starts = 0
        self.commands = []
        self.fit_after_loads = []
        self.queries = []
        self.restore_geometry_on_load = False

    def ensure_running(self):
        self.starts += 1

    def verify_image_loaded(self):
        if not any(frame == self.frame for frame, _, _ in self.loaded):
            raise RuntimeError('DS9 has no FITS image in the active frame; loading failed.')

    def state(self):
        return self.frames.setdefault(self.frame, {
            ('scale',): 'asinh', ('scale', 'limits'): '0 1000',
            ('zoom',): '1 1', ('pan', 'image'): '8 8',
            ('cmap',): 'gray', ('cmap', 'value'): '1 0.5', ('cmap', 'invert'): 'no',
            ('scale', 'log', 'exp'): '1000',
        })

    def get(self, *args):
        self.queries.append(args)
        return self.frame if args == ('frame',) else self.state()[args]

    def set(self, *args):
        self.commands.append(args)
        if args[0] == 'frame':
            self.frame = args[1]
        elif args[:3] == ('scale', 'log', 'exp'):
            self.state()[('scale', 'log', 'exp')] = args[3]
        elif args[0] == 'scale' and args[1] not in ('mode', 'lock'):
            key = ('scale', 'limits') if args[1] == 'limits' else ('scale',)
            self.state()[key] = ' '.join(args[2:] if key[-1] == 'limits' else args[1:])
        elif args[:2] == ('zoom', 'to'):
            self.state()[('zoom',)] = ' '.join(args[2:])
        elif args[:2] == ('pan', 'to'):
            self.state()[('pan', 'image')] = ' '.join(args[2:4])
        elif args[0] == 'cmap':
            if args[1] == 'invert':
                self.state()[('cmap', 'invert')] = args[2]
            elif len(args) == 3:
                self.state()[('cmap', 'value')] = ' '.join(args[1:])
            else:
                self.state()[('cmap',)] = args[1]

    def load_fits(self, path, *, configure_display=True, hdu=None):
        assert configure_display is False
        data = fits.getdata(path) if hdu is None else fits.getdata(path, hdu)
        if self.restore_geometry_on_load:
            # Simulate the native viewer recalling an older per-file view.
            self.state()[('zoom',)] = '1 1'
            self.state()[('pan', 'image')] = '9 9'
        self.loaded.append((self.frame, Path(path), data.copy()))

    def set_scale_mode(self, mode):
        self.set('scale', 'mode', mode)

    def set_scale_function(self, function):
        self.set('scale', function)

    def set_scale_limits(self, low, high):
        self.set('scale', 'limits', str(low), str(high))

    def set_cmap(self, name):
        self.set('cmap', name)

    def set_cmap_parameters(self, contrast, bias):
        self.set('cmap', str(contrast), str(bias))

    def zoom_to_fit(self):
        self.fit_after_loads.append(len(self.loaded))
        self.set('zoom', 'to', '0.5', '0.5')


def pair(directory, name):
    backup = directory / 'before_cosmic_ray'
    backup.mkdir(exist_ok=True)
    before, after = backup / name, directory / name
    mask = np.zeros((16, 16), dtype=bool)
    mask[4, 5] = True
    CCDData(np.full((16, 16), 100.), unit=u.adu).write(before)
    CCDData(np.full((16, 16), 100.), unit=u.adu).write(after)
    with fits.open(after, mode='update') as hdus:
        hdus.append(fits.ImageHDU(mask.astype(np.uint8), name='CMASK'))
    return before, after


@pytest.mark.usefixtures('fake_tkinter')
def test_next_previous_sync_and_fixed_mask_without_file_moves(tmp_path):
    pairs = [pair(tmp_path, 'a.fits'), pair(tmp_path, 'b.fits')]
    originals = {path: path.read_bytes() for item in pairs for path in item}
    adapter = FakeDS9()
    reviewer = ImageComparisonReviewer(pairs, adapter=adapter, isolated=False)
    reviewer.show_current()
    assert [frame for frame, _, _ in adapter.loaded] == ['1', '2', '3']
    assert set(np.unique(adapter.loaded[2][2])) == {0, 1}
    assert adapter.starts == 1
    assert ('tile', 'grid', 'layout', '3', '1') in adapter.commands
    assert ('tile', 'grid', 'mode', 'manual') in adapter.commands
    assert ('tile', 'mode', 'row') not in adapter.commands
    assert adapter.fit_after_loads == [3]
    assert all(adapter.frames[frame][('zoom',)] == '0.5 0.5' for frame in ('1', '2', '3'))

    adapter.set('frame', '2')
    adapter.set('scale', 'log')
    adapter.set('scale', 'limits', '20', '200')
    adapter.set('zoom', 'to', '3', '3')
    adapter.set('pan', 'to', '6', '7', 'image')
    reviewer.sync_views()
    assert adapter.frames['1'] == adapter.frames['2']
    assert adapter.frames['3'][('scale',)] == 'linear'
    assert adapter.frames['3'][('scale', 'limits')] == '0.0 1.0'

    commands_before = len(adapter.commands)
    queries_before = len(adapter.queries)
    reviewer.mark_good()
    assert reviewer.index == 1
    assert adapter.frames['2'][('scale',)] == 'log'
    assert adapter.frames['2'][('scale', 'limits')] == '20 200'
    assert adapter.fit_after_loads == [3]
    assert adapter.frames['2'][('zoom',)] == '3 3'
    assert adapter.queries[queries_before:] == [('frame',), ('zoom',), ('pan', 'image')]
    assert all(command[0] in ('frame', 'zoom', 'pan') for command in adapter.commands[commands_before:])
    reviewer.go_back()
    assert reviewer.index == 0
    with pytest.raises(RuntimeError, match='cannot reject'):
        reviewer.mark_reject()
    assert {path: path.read_bytes() for path in originals} == originals


@pytest.mark.usefixtures('fake_tkinter')
def test_run_reads_mask_extension_without_temporary_image_or_file_moves(tmp_path):
    pairs = [pair(tmp_path, 'a.fits')]
    adapter = FakeDS9()
    reviewer = ImageComparisonReviewer(pairs, adapter=adapter, isolated=False)
    reviewer.root.on_update = reviewer.quit_review
    assert reviewer.run() == []
    assert adapter.loaded[2][1] == pairs[0][1]
    assert not list(tmp_path.glob('mask.fits'))
    assert all(path.exists() for path in pairs[0])
    assert reviewer.root.destroy_called


@pytest.mark.usefixtures('fake_tkinter')
@pytest.mark.parametrize('source_frame', ['1', '2', '3'])
def test_navigation_restores_live_view_instead_of_old_file_geometry(tmp_path, source_frame):
    adapter = FakeDS9()
    pairs = [pair(tmp_path, 'a.fits'), pair(tmp_path, 'b.fits')]
    reviewer = ImageComparisonReviewer(pairs, adapter=adapter, isolated=False)
    reviewer.show_current()
    adapter.restore_geometry_on_load = True

    # No Sync click: navigation must read the user's most recent manual view.
    adapter.set('frame', source_frame)
    adapter.set('zoom', 'to', '0.4', '0.4')
    adapter.set('pan', 'to', '6', '7', 'image')
    reviewer.mark_good()
    for frame in ('1', '2', '3'):
        assert adapter.frames[frame][('zoom',)] == '0.4 0.4'
        assert adapter.frames[frame][('pan', 'image')] == '6 7'
    assert adapter.frame == source_frame

    # Previous must carry the new view, not a cached Sync view or old file view.
    adapter.set('zoom', 'to', '0.7', '0.7')
    adapter.set('pan', 'to', '3', '4', 'image')
    reviewer.go_back()
    for frame in ('1', '2', '3'):
        assert adapter.frames[frame][('zoom',)] == '0.7 0.7'
        assert adapter.frames[frame][('pan', 'image')] == '3 4'
    assert adapter.frame == source_frame
    assert adapter.fit_after_loads == [3]
    assert adapter.frames['3'][('scale',)] == 'linear'
    assert adapter.frames['3'][('scale', 'limits')] == '0.0 1.0'


@pytest.mark.usefixtures('fake_tkinter')
@pytest.mark.parametrize('source_frame', ['1', '2', '3'])
def test_sync_uses_selected_frame_and_separates_mask_scale(tmp_path, source_frame):
    adapter = FakeDS9()
    reviewer = ImageComparisonReviewer([pair(tmp_path, 'a.fits')], adapter=adapter, isolated=False)
    reviewer.show_current()

    # Make science displays distinct so Mask-driven synchronization cannot
    # silently replace either with the other science frame's display settings.
    adapter.frames['1'][('scale',)] = 'log'
    adapter.frames['1'][('scale', 'limits')] = '10 100'
    adapter.frames['2'][('scale',)] = 'sqrt'
    adapter.frames['2'][('scale', 'limits')] = '20 200'
    adapter.set('frame', source_frame)
    adapter.set('zoom', 'to', '5', '5')
    adapter.set('pan', 'to', '4', '5', 'image')
    displays = {
        frame: {key: value for key, value in state.items() if key not in [('zoom',), ('pan', 'image')]}
        for frame, state in adapter.frames.items()
    }

    reviewer.sync_views()

    assert adapter.frame == source_frame
    for frame in ('1', '2', '3'):
        assert adapter.frames[frame][('zoom',)] == '5 5'
        assert adapter.frames[frame][('pan', 'image')] == '4 5'
        actual_display = {
            key: value for key, value in adapter.frames[frame].items()
            if key not in [('zoom',), ('pan', 'image')]
        }
        expected = displays[frame] if source_frame == '3' or frame == '3' else displays[source_frame]
        assert actual_display == expected


def test_missing_before_fails_without_creating_a_window(tmp_path):
    with pytest.raises(FileNotFoundError):
        ImageComparisonReviewer([(tmp_path / 'missing.fits', tmp_path / 'also_missing.fits')])


def test_missing_mask_disables_next_and_allows_retry(tmp_path, request):
    errors = request.getfixturevalue('fake_tkinter')
    before, after = pair(tmp_path, 'a.fits')
    with fits.open(after, mode='update') as hdus:
        del hdus['CMASK']
    reviewer = ImageComparisonReviewer([(before, after)], adapter=FakeDS9(), isolated=False)
    reviewer.show_current()
    assert reviewer.good_button.options['state'] == 'disabled'
    assert reviewer.retry_button.options['state'] == 'normal'
    assert 'CMASK' in errors[0][1]


def test_silently_ignored_load_keeps_comparison_in_error_state(tmp_path, request):
    errors = request.getfixturevalue('fake_tkinter')
    adapter = FakeDS9()
    adapter.load_fits = lambda *args, **kwargs: None
    reviewer = ImageComparisonReviewer([pair(tmp_path, 'a.fits')], adapter=adapter, isolated=False)
    reviewer.show_current()
    assert not reviewer._loaded
    assert reviewer.good_button.options['state'] == 'disabled'
    assert reviewer.retry_button.options['state'] == 'normal'
    assert 'no FITS image' in errors[0][1]
