"""Preview navigation and isolation without DS9, Tk, or observing data."""

from unittest.mock import Mock

from astropy.io import fits
from ccdproc import ImageFileCollection
import numpy as np
import pytest

from photozpy import Combine
from photozpy.interactive_ds9.adapter import DS9Adapter
from photozpy.interactive_ds9.image_previewer import ImagePreviewer
from test_interactive_ds9_image_reviewer import fake_tkinter as fake_tkinter


def products(tmp_path):
    paths = [tmp_path / 'Master_a.fits', tmp_path / 'Master_b.fits']
    for path in paths:
        fits.PrimaryHDU(np.ones((8, 8))).writeto(path)
    return paths


@pytest.mark.usefixtures('fake_tkinter')
def test_navigation_never_moves_products_and_preserves_previous_at_end(tmp_path):
    paths = products(tmp_path)
    originals = {path: path.read_bytes() for path in paths}
    adapter = Mock(spec=DS9Adapter)
    preview = ImagePreviewer(paths, adapter=adapter, isolated=False)
    preview.show_current()
    assert preview.reject_button.pack_options['hidden']
    assert preview.reject_button.options['state'] == 'disabled'
    preview.mark_good()
    assert preview.current_image == paths[1]
    preview.mark_good()
    assert 'Preview complete' in preview.status.get()
    preview.go_back()
    assert preview.current_image == paths[1]
    preview.go_back()
    assert preview.current_image == paths[0]
    with pytest.raises(RuntimeError, match='cannot reject'):
        preview.mark_reject()
    preview.quit_review()
    adapter.ensure_running.assert_called_once()
    assert {path: path.read_bytes() for path in originals} == originals
    assert set(tmp_path.iterdir()) == set(paths)


@pytest.mark.usefixtures('fake_tkinter')
def test_isolated_preview_passes_read_only_mode_to_worker(tmp_path, monkeypatch):
    from photozpy.interactive_ds9 import _review_process

    paths = products(tmp_path)
    worker = Mock()
    monkeypatch.setattr(_review_process, 'run_isolated', worker)
    preview = ImagePreviewer(paths, isolated=True)
    assert preview.run() == []
    assert worker.call_args.kwargs == {'preview_only': True}
    assert worker.call_args.args[2]['software_name'] == 'photozpy-master-preview'
    with pytest.raises(RuntimeError, match='closed'):
        preview.run()


def test_preview_worker_dispatches_to_nonmoving_previewer(tmp_path, monkeypatch):
    import json
    from photozpy.interactive_ds9 import _review_process, adapter, image_previewer

    request, state = tmp_path / 'request.json', tmp_path / 'state.json'
    request.write_text(json.dumps({
        'images': ['Master_a.fits'], 'adapter': None, 'preview_only': True,
    }))
    constructor = Mock()
    monkeypatch.setattr(adapter, 'DS9Adapter', Mock())
    monkeypatch.setattr(image_previewer, 'ImagePreviewer', constructor)
    monkeypatch.setattr(_review_process.sys, 'argv', ['worker', str(request), str(state)])
    _review_process.main()
    assert constructor.call_args.args == (['Master_a.fits'],)
    assert constructor.call_args.kwargs['isolated'] is False
    constructor.return_value.run.assert_called_once()
    assert constructor.return_value._state_path == state


def test_combine_preview_selects_only_last_masters_without_recombining(tmp_path, monkeypatch):
    from photozpy.interactive_ds9 import image_previewer

    for index in range(2):
        fits.PrimaryHDU(
            np.full((8, 8), 10. + index),
            fits.Header({'IMTYPE': 'Bias', 'BUNIT': 'ADU'}),
        ).writeto(tmp_path / f'bias_{index}.fits')
    combiner = Combine(ImageFileCollection(tmp_path, keywords='*'))
    with pytest.raises(ValueError, match='No completed combine'):
        combiner.preview_in_ds9()
    result = combiner.combine_bias_or_dark('Bias', quite=True)
    originals = {path: path.read_bytes() for path in tmp_path.glob('*.fits')}
    constructor = Mock()
    monkeypatch.setattr(image_previewer, 'ImagePreviewer', constructor)
    combiner.preview_in_ds9()
    constructor.assert_called_once_with((tmp_path / 'Master_Bias.fits',))
    constructor.return_value.run.assert_called_once()
    assert combiner.image_collection is result
    assert {path: path.read_bytes() for path in originals} == originals


def test_missing_products_fail_before_creating_gui(tmp_path):
    with pytest.raises(FileNotFoundError):
        ImagePreviewer([tmp_path / 'missing.fits'])
    with pytest.raises(ValueError, match='At least one'):
        ImagePreviewer([])
