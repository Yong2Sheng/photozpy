"""Post-cleaning review safety checks on temporary synthetic exposures."""

import json
from pathlib import Path

from astropy.io import fits
from ccdproc import ImageFileCollection
import pytest

from photozpy import Artifacts, ImageCollectionGroup
from photozpy.calibration import _artifact_review as file_review
from photozpy.interactive_ds9 import cosmic_ray_reviewer as gui
from photozpy.interactive_ds9.cosmic_ray_reviewer import CosmicRayReviewer
from test_artifacts import image, fake_lacosmic as fake_lacosmic
from test_interactive_ds9_comparison import FakeDS9
from test_interactive_ds9_image_reviewer import fake_tkinter as fake_tkinter


pytestmark = pytest.mark.usefixtures('fake_lacosmic')


def batch(root, *, names=('a.fits', 'b.fits')):
    """Create raw originals and matching uncleaned working copies."""
    analysis = root / 'analysis'
    for name in names:
        raw = image(root, name, EXPTIME=100., **{'DATE-OBS': '2024-02-16T01:40:05'})
        with fits.open(raw, mode='update') as hdul:
            for key in ('BIASCORR', 'DARKCORR', 'FLATCORR'):
                del hdul[0].header[key]
        image(analysis, name, EXPTIME=100., **{'DATE-OBS': '2024-02-16T01:40:05'})
    return analysis


def decisions(monkeypatch, selected=(), error=None):
    """Supply GUI decisions without opening native windows or DS9."""
    class Reviewer:
        def __init__(self, pairs):
            self.pairs = pairs

        def run(self):
            if error is not None:
                raise error
            return [after for _, after in self.pairs if after.name in selected]

    monkeypatch.setattr(gui, 'CosmicRayReviewer', Reviewer)


def test_reject_moves_original_cleans_products_and_updates_selection(tmp_path, monkeypatch):
    analysis = batch(tmp_path)
    bias = image(analysis, 'Master_Bias.fits', IMTYPE='Master Bias')
    raw = tmp_path / 'a.fits'
    original = raw.read_bytes()
    cleaner = Artifacts(ImageCollectionGroup(analysis), overwrite=True)
    cleaner.remove_cosmic_ray(plot_mask=True)
    unrelated = analysis / 'a.reg'
    unrelated.write_text('manual regions')
    decisions(monkeypatch, ['a.fits'])

    result = cleaner.review_in_ds9()

    assert result is cleaner.image_collection
    assert isinstance(result, ImageCollectionGroup)
    assert result[0].files == ['Master_Bias.fits', 'b.fits']
    assert (tmp_path / 'rejected' / raw.name).read_bytes() == original
    assert not raw.exists()
    assert not (analysis / raw.name).exists()
    assert not (analysis / 'before_cosmic_ray' / raw.name).exists()
    assert not (analysis / 'a_cr_mask.png').exists()
    assert unrelated.read_text() == 'manual regions'
    assert bias.exists()
    assert list(cleaner.cosmic_ray_counts) == [analysis / 'b.fits']
    assert list(cleaner.backup_paths) == [analysis / 'b.fits']
    record = json.loads((tmp_path / 'rejected' / 'a.fits.rejection.json').read_text())
    assert record['stage'] == 'post_cosmic_ray_review'
    assert record['original'] == str(raw)
    # Root discovery must not recursively reintroduce rejected exposures.
    rerun = ImageFileCollection(tmp_path, keywords='*')
    assert rerun.files == ['b.fits']


@pytest.mark.parametrize('error', [None, RuntimeError('GUI failed'), KeyboardInterrupt()])
def test_no_rejection_or_gui_failure_preserves_all_files(tmp_path, monkeypatch, error):
    analysis = batch(tmp_path)
    cleaner = Artifacts(ImageCollectionGroup(analysis), overwrite=True)
    cleaner.remove_cosmic_ray()
    originals = {path: path.read_bytes() for path in tmp_path.rglob('*.fits')}
    decisions(monkeypatch, ['a.fits'] if error else [], error)
    if error is None:
        assert cleaner.review_in_ds9() is cleaner.image_collection
    else:
        with pytest.raises(type(error)):
            cleaner.review_in_ds9()
    assert all(path.read_bytes() == content for path, content in originals.items())
    assert not (tmp_path / 'rejected').exists()


@pytest.mark.parametrize('problem', ['missing', 'date', 'exposure', 'archive', 'record', 'registered', 'raw_alias'])
def test_preflight_blocks_unsafe_rejection(tmp_path, monkeypatch, problem):
    analysis = batch(tmp_path)
    cleaner = Artifacts(ImageCollectionGroup(analysis), overwrite=True)
    cleaner.remove_cosmic_ray()
    raw = tmp_path / 'a.fits'
    if problem == 'missing':
        raw.unlink()
    elif problem in ('date', 'exposure'):
        with fits.open(raw, mode='update') as hdul:
            hdul[0].header['DATE-OBS' if problem == 'date' else 'EXPTIME'] = 'wrong'
    elif problem in ('archive', 'record'):
        rejected = tmp_path / 'rejected'
        rejected.mkdir()
        (rejected / ('a.fits' if problem == 'archive' else 'a.fits.rejection.json')).write_bytes(b'keep')
    elif problem == 'registered':
        with fits.open(analysis / 'a.fits', mode='update') as hdul:
            hdul[0].header['ALIGN'] = 'Registered'
    else:
        raw.unlink()
        raw.symlink_to(analysis / 'a.fits')
    originals = {path: path.read_bytes() for path in tmp_path.rglob('*') if path.is_file()}

    def unexpected_gui(*args, **kwargs):
        pytest.fail('Unsafe paths reached the GUI.')

    monkeypatch.setattr(gui, 'CosmicRayReviewer', unexpected_gui)
    with pytest.raises((ValueError, FileExistsError)):
        cleaner.review_in_ds9()
    assert all(path.read_bytes() == content for path, content in originals.items())


def test_reject_all_is_empty_and_stays_empty_after_refresh(tmp_path, monkeypatch):
    analysis = batch(tmp_path, names=('a.fits',))
    cleaner = Artifacts(ImageCollectionGroup(analysis), overwrite=True)
    cleaner.remove_cosmic_ray()
    decisions(monkeypatch, ['a.fits'])
    result = cleaner.review_in_ds9()
    assert result[0].files == []
    image(analysis, 'unselected.fits')
    result[0].refresh()
    assert result[0].files == []
    from photozpy.calibration._calibration_io import select
    assert select(result[0], IMTYPE='Light') == []
    assert cleaner.backup_paths == cleaner.cosmic_ray_counts == {}


def test_multi_root_inference_preserves_order(tmp_path, monkeypatch):
    directories = [batch(tmp_path / name) for name in ('night_b', 'night_a')]
    cleaner = Artifacts(ImageCollectionGroup(directories), overwrite=True)
    cleaner.remove_cosmic_ray()
    decisions(monkeypatch, ['a.fits'])
    result = cleaner.review_in_ds9()
    assert [child.files for child in result] == [['b.fits'], ['b.fits']]
    assert [Path(child.location) for child in result] == directories
    assert all((directory.parent / 'rejected' / 'a.fits').exists() for directory in directories)


def test_failed_archive_publication_does_not_delete_products(tmp_path, monkeypatch):
    analysis = batch(tmp_path)
    cleaner = Artifacts(ImageCollectionGroup(analysis), overwrite=True)
    cleaner.remove_cosmic_ray()
    decisions(monkeypatch, ['a.fits'])
    original_link = file_review.os.link

    def fail_record(source, target):
        if str(target).endswith('.json'):
            raise OSError('record disk full')
        return original_link(source, target)

    monkeypatch.setattr(file_review.os, 'link', fail_record)
    with pytest.raises(OSError, match='disk full'):
        cleaner.review_in_ds9()
    assert (tmp_path / 'a.fits').exists()
    assert (analysis / 'a.fits').exists()
    assert (analysis / 'before_cosmic_ray' / 'a.fits').exists()
    assert not (tmp_path / 'rejected' / 'a.fits').exists()
    assert len(cleaner.backup_paths) == 2


def test_cleanup_failure_keeps_raw_and_excludes_rejected_result(tmp_path, monkeypatch):
    analysis = batch(tmp_path)
    cleaner = Artifacts(ImageCollectionGroup(analysis), overwrite=True)
    cleaner.remove_cosmic_ray()
    decisions(monkeypatch, ['a.fits'])
    original_unlink = Path.unlink

    def fail_cleaned(path, *args, **kwargs):
        if path == analysis / 'a.fits':
            raise PermissionError('cannot remove cleaned file')
        return original_unlink(path, *args, **kwargs)

    monkeypatch.setattr(Path, 'unlink', fail_cleaned)
    with pytest.raises(PermissionError):
        cleaner.review_in_ds9()
    assert (tmp_path / 'rejected' / 'a.fits').exists()
    assert cleaner.image_collection[0].files == ['b.fits']
    assert list(cleaner.backup_paths) == [analysis / 'b.fits']


@pytest.mark.usefixtures('fake_tkinter')
def test_gui_stages_reject_and_previous_undoes_without_file_moves(tmp_path):
    from test_interactive_ds9_comparison import pair
    pairs = [pair(tmp_path, name) for name in ('a.fits', 'b.fits')]
    originals = {path: path.read_bytes() for item in pairs for path in item}
    reviewer = CosmicRayReviewer(pairs, adapter=FakeDS9(), isolated=False)
    reviewer.show_current()
    reviewer.mark_reject()
    assert reviewer.rejected == [pairs[0][1]]
    reviewer.go_back()
    assert reviewer.rejected == []
    reviewer.mark_good()
    reviewer.mark_reject()
    assert reviewer.rejected == [pairs[1][1]]
    assert reviewer.good_button.cget('state') == 'disabled'
    assert reviewer.reject_button.cget('state') == 'disabled'
    reviewer.go_back()
    assert reviewer.index == 1 and not reviewer.rejected
    assert all(path.read_bytes() == content for path, content in originals.items())


def test_isolated_gui_passes_review_mode_and_returns_decisions(tmp_path, monkeypatch):
    from test_interactive_ds9_comparison import pair
    from photozpy.interactive_ds9 import _review_process
    from photozpy.interactive_ds9.image_reviewer import ImageReviewer
    monkeypatch.setattr(ImageReviewer, '_debug', lambda *args: None)
    pairs = [pair(tmp_path, 'a.fits')]
    reviewer = CosmicRayReviewer(pairs, isolated=True)

    def run(images, rejected_dir, config, restore, **kwargs):
        assert kwargs['cosmic_ray_review'] is True
        assert kwargs['comparison_pairs'] == pairs
        assert config['software_name'] == 'photozpy-cr-review'
        restore({'index': 1, 'history': ['reject'], 'rejected': [str(pairs[0][1])]})

    monkeypatch.setattr(_review_process, 'run_isolated', run)
    assert reviewer.run() == [pairs[0][1]]
    assert all(path.exists() for path in pairs[0])


@pytest.mark.parametrize('after_unlink', [False, True])
def test_interrupt_at_raw_unlink_never_loses_original(tmp_path, monkeypatch, after_unlink):
    analysis = batch(tmp_path)
    cleaner = Artifacts(ImageCollectionGroup(analysis), overwrite=True)
    cleaner.remove_cosmic_ray()
    raw = tmp_path / 'a.fits'
    original = raw.read_bytes()
    decisions(monkeypatch, ['a.fits'])
    unlink = Path.unlink

    def interrupt(path, *args, **kwargs):
        if path == raw:
            if after_unlink:
                unlink(path, *args, **kwargs)
            raise KeyboardInterrupt()
        return unlink(path, *args, **kwargs)

    monkeypatch.setattr(Path, 'unlink', interrupt)
    with pytest.raises(KeyboardInterrupt):
        cleaner.review_in_ds9()
    retained = tmp_path / 'rejected' / 'a.fits'
    assert (retained if after_unlink else raw).read_bytes() == original
    assert (analysis / 'a.fits').exists()
    assert (analysis / 'a.fits' in cleaner.backup_paths) is not after_unlink
    assert ('a.fits' in cleaner.image_collection[0].files) is not after_unlink


@pytest.mark.parametrize('problem', ['layout', 'missing_date', 'raw_is_calibrated'])
def test_layout_and_raw_identity_are_required(tmp_path, monkeypatch, problem):
    analysis = batch(tmp_path)
    cleaner = Artifacts(ImageCollectionGroup(analysis), overwrite=True)
    cleaner.remove_cosmic_ray()
    decisions(monkeypatch, ['a.fits'])
    if problem == 'layout':
        output = tmp_path / 'other_output'
        analysis.rename(output)
        cleaner._backup_paths = {output / 'a.fits': output / 'before_cosmic_ray' / 'a.fits'}
    else:
        with fits.open(tmp_path / 'a.fits', mode='update') as hdul:
            if problem == 'missing_date':
                del hdul[0].header['DATE-OBS']
            else:
                hdul[0].header['BIASCORR'] = 'Yes'
    with pytest.raises(ValueError):
        cleaner.review_in_ds9()
    assert (tmp_path / 'a.fits').exists()
    assert not (tmp_path / 'rejected').exists()


@pytest.mark.usefixtures('fake_tkinter')
@pytest.mark.parametrize('failure', [RuntimeError('window failed'), KeyboardInterrupt()])
def test_gui_run_failure_propagates_without_file_changes(tmp_path, failure):
    from test_interactive_ds9_comparison import pair
    pairs = [pair(tmp_path, 'a.fits')]
    reviewer = CosmicRayReviewer(pairs, adapter=FakeDS9(), isolated=False)

    def fail():
        reviewer.mark_reject()
        raise failure

    reviewer.root.on_update = fail
    with pytest.raises(type(failure)):
        reviewer.run()
    assert all(path.exists() for path in pairs[0])
    assert reviewer.root.destroy_called


@pytest.mark.usefixtures('fake_tkinter')
def test_review_child_dispatch_uses_quality_reviewer(tmp_path, monkeypatch):
    from photozpy.interactive_ds9 import _review_process
    from test_interactive_ds9_comparison import pair
    import sys
    pairs = [pair(tmp_path, 'a.fits')]
    request = tmp_path / 'request.json'
    request.write_text(json.dumps({
        'adapter': {}, 'comparison_pairs': [[str(path) for path in pairs[0]]],
        'cosmic_ray_review': True,
    }))
    monkeypatch.setattr(sys, 'argv', ['review', str(request), str(tmp_path / 'state.json')])
    from photozpy.interactive_ds9 import adapter as adapter_module
    monkeypatch.setattr(adapter_module, 'DS9Adapter', lambda **kwargs: FakeDS9())
    calls = []

    def run(reviewer):
        calls.append(type(reviewer))
        assert reviewer._isolated is False
        assert reviewer._state_path == tmp_path / 'state.json'

    monkeypatch.setattr(CosmicRayReviewer, 'run', run)
    _review_process.main()
    assert calls == [CosmicRayReviewer]


@pytest.mark.usefixtures('fake_tkinter')
def test_failed_new_pair_load_disables_reject_even_after_success(tmp_path):
    from test_interactive_ds9_comparison import pair
    pairs = [pair(tmp_path, name) for name in ('a.fits', 'b.fits')]
    reviewer = CosmicRayReviewer(pairs, adapter=FakeDS9(), isolated=False)
    reviewer.show_current()
    with fits.open(pairs[1][1], mode='update') as hdul:
        del hdul['CMASK']
    reviewer.mark_good()
    assert reviewer.reject_button.cget('state') == 'disabled'
    reviewer.mark_reject()
    assert reviewer.rejected == []
    assert reviewer.index == 1


def test_conflicts_created_during_review_are_rechecked(tmp_path, monkeypatch):
    analysis = batch(tmp_path)
    cleaner = Artifacts(ImageCollectionGroup(analysis), overwrite=True)
    cleaner.remove_cosmic_ray()

    class Reviewer:
        def __init__(self, pairs):
            self.pairs = pairs

        def run(self):
            destination = tmp_path / 'rejected' / 'a.fits'
            destination.parent.mkdir()
            destination.write_bytes(b'previous rejection')
            return [analysis / 'a.fits']

    monkeypatch.setattr(gui, 'CosmicRayReviewer', Reviewer)
    with pytest.raises(FileExistsError):
        cleaner.review_in_ds9()
    assert (tmp_path / 'a.fits').exists()
    assert (analysis / 'a.fits').exists()
    assert (tmp_path / 'rejected' / 'a.fits').read_bytes() == b'previous rejection'


@pytest.mark.usefixtures('fake_tkinter')
def test_normal_early_quit_returns_pending_rejects_without_moving_files(tmp_path):
    from test_interactive_ds9_comparison import pair
    pairs = [pair(tmp_path, name) for name in ('a.fits', 'b.fits')]
    reviewer = CosmicRayReviewer(pairs, adapter=FakeDS9(), isolated=False)

    def choose():
        reviewer.mark_reject()
        reviewer.quit_review()

    reviewer.root.on_update = choose
    assert reviewer.run() == [pairs[0][1]]
    assert reviewer.index == 1
    assert all(path.exists() for item in pairs for path in item)
