"""Cosmic-ray cleaning on temporary synthetic CCDs, never observing data."""

import astropy.units as u
from astropy.io import fits
from astropy.nddata import CCDData
from ccdproc import ImageFileCollection
import numpy as np
import pytest

from photozpy import Artifacts, ImageCollectionGroup, Registration
from photozpy.calibration import artifacts as artifacts_module
from photozpy.calibration import registration as registration_module


def image(directory, name='light.fits', *, mask=None, data=None, **cards):
    directory.mkdir(parents=True, exist_ok=True)
    header = fits.Header({
        'IMTYPE': 'Light', 'FIELD': 'field_a', 'FILTER': "SDSS_g'",
        'BUNIT': 'ADU', 'GAIN': 2., 'RDNOISE': 7.1,
        'BIASCORR': 'Yes', 'DARKCORR': 'Yes', 'FLATCORR': 'Yes',
    })
    header.update(cards)
    data = np.full((16, 16), 100.) if data is None else data
    hdus = [fits.PrimaryHDU(data, header)]
    if mask is not None:
        hdus.append(fits.ImageHDU(np.uint8(mask), name='MASK'))
    path = directory / name
    fits.HDUList(hdus).writeto(path)
    return path


def collection(directory):
    return ImageFileCollection(location=directory, keywords='*', glob_exclude='.*')


@pytest.fixture
def fake_lacosmic(monkeypatch):
    calls = []

    def clean(ccd, *, gain_apply, gain, readnoise, **kwargs):
        calls.append((gain_apply, gain, readnoise, kwargs))
        assert ccd.data.dtype.isnative
        result = ccd.copy()
        result.mask = np.zeros(ccd.shape, dtype=bool)
        result.mask[5, 6] = True
        result.data[5, 6] = 101.
        # Reproduce CCDProc 2.5.1's unit label with gain_apply=False.
        result.unit = u.electron
        return result

    monkeypatch.setattr(artifacts_module, 'lacosmic', clean)
    return calls


def test_header_parameters_mask_union_and_retained_calibration(tmp_path, fake_lacosmic):
    mask = np.zeros((16, 16), dtype=bool)
    mask[2, 3] = True
    light = image(tmp_path, mask=mask)
    bias = image(tmp_path, 'Master_Bias.fits', IMTYPE='Master Bias')
    bias_bytes = bias.read_bytes()
    cleaner = Artifacts(ImageCollectionGroup(tmp_path), overwrite=True)

    result = cleaner.remove_cosmic_ray(sigclip=6., objlim=8.)

    assert isinstance(result, ImageCollectionGroup)
    assert result is cleaner.image_collection
    assert fake_lacosmic == [(False, 2., 7.1, {'sigclip': 6., 'objlim': 8.})]
    assert bias.read_bytes() == bias_bytes
    cleaned = CCDData.read(light)
    assert cleaned.unit == u.adu
    assert cleaned.mask[2, 3] and not cleaned.mask[5, 6]
    assert cleaned.mask.sum() == 1
    cmask = fits.getdata(light, 'CMASK').astype(bool)
    assert cmask[5, 6] and not cmask[2, 3]
    assert cmask.sum() == 1
    assert cleaned.meta['CRCORR'] == 'Yes'
    assert cleaned.meta['CRNPIX'] == 1
    assert cleaned.meta['CRMETHOD'] == 'LACosmic'
    assert cleaned.data[0, 0] == 100.
    assert cleaner.cosmic_ray_counts == {light: 1}


def test_separate_output_preserves_source_and_collection_members(tmp_path, fake_lacosmic):
    source, output = tmp_path / 'analysis', tmp_path / 'output'
    light = image(source)
    flat = image(source, 'Master_Flat.fits', IMTYPE='Master Flat')
    originals = {path: path.read_bytes() for path in (light, flat)}
    result = Artifacts(collection(source)).remove_cosmic_ray(save_location=output)

    assert {path: path.read_bytes() for path in originals} == originals
    assert (output / light.name).exists()
    assert not (output / flat.name).exists()
    assert len(result.files) == 2
    assert fits.getheader(output / light.name)['CRCORR'] == 'Yes'


def test_in_place_backups_are_exact_and_excluded_from_working_collection(tmp_path, fake_lacosmic):
    light = image(tmp_path)
    original = light.read_bytes()
    cleaner = Artifacts(collection(tmp_path), overwrite=True)
    result = cleaner.remove_cosmic_ray()

    backup = tmp_path / 'before_cosmic_ray' / light.name
    assert backup.read_bytes() == original
    assert cleaner.backup_paths == {light: backup}
    assert result.files == [light.name]
    assert fits.getheader(light)['CRCORR'] == 'Yes'


def test_existing_backup_is_never_overwritten(tmp_path, fake_lacosmic):
    light = image(tmp_path)
    before = light.read_bytes()
    backup = image(tmp_path / 'before_cosmic_ray')
    backup_before = backup.read_bytes()

    with pytest.raises(FileExistsError):
        Artifacts(collection(tmp_path), overwrite=True).remove_cosmic_ray()

    assert not fake_lacosmic
    assert light.read_bytes() == before
    assert backup.read_bytes() == backup_before


def test_backup_failure_precedes_any_cleaning(tmp_path, fake_lacosmic, monkeypatch):
    light = image(tmp_path)
    before = light.read_bytes()

    def fail(*args):
        raise OSError('No space left')

    monkeypatch.setattr(Artifacts, '_copy_backup', staticmethod(fail))
    with pytest.raises(OSError, match='No space left'):
        Artifacts(collection(tmp_path), overwrite=True).remove_cosmic_ray()
    assert not fake_lacosmic
    assert light.read_bytes() == before


def test_explicit_backup_collision_does_not_touch_inputs(tmp_path, fake_lacosmic):
    first = image(tmp_path / 'a')
    second = image(tmp_path / 'b')
    before = {path: path.read_bytes() for path in (first, second)}
    with pytest.raises(ValueError, match='Multiple images'):
        Artifacts(ImageCollectionGroup([first.parent, second.parent]), overwrite=True).remove_cosmic_ray(
            backup_dir=tmp_path / 'backups',
        )
    assert not fake_lacosmic
    assert {path: path.read_bytes() for path in before} == before


def test_preview_opens_only_after_batch_and_can_be_reopened(tmp_path, fake_lacosmic, monkeypatch):
    light = image(tmp_path)
    cleaner = Artifacts(collection(tmp_path), overwrite=True)
    reviews = []

    class Reviewer:
        def __init__(self, pairs):
            reviews.append(pairs)
            assert fits.getheader(light)['CRCORR'] == 'Yes'
            assert cleaner.image_collection.files == [light.name]

        def run(self):
            assert len(fake_lacosmic) == 1

    from photozpy.interactive_ds9 import comparison_reviewer
    monkeypatch.setattr(comparison_reviewer, 'ImageComparisonReviewer', Reviewer)

    cleaner.remove_cosmic_ray(preview=True)
    cleaner.preview_in_ds9()
    backup = tmp_path / 'before_cosmic_ray' / light.name
    assert reviews == [[(backup, light)], [(backup, light)]]


def test_preview_failure_keeps_cleaned_collection_and_backup(tmp_path, fake_lacosmic, monkeypatch):
    light = image(tmp_path)
    original = light.read_bytes()
    cleaner = Artifacts(collection(tmp_path), overwrite=True)

    def fail():
        raise RuntimeError('DS9 unavailable')

    monkeypatch.setattr(cleaner, 'preview_in_ds9', fail)
    with pytest.raises(RuntimeError, match='DS9 unavailable'):
        cleaner.remove_cosmic_ray(preview=True)
    assert cleaner.backup_paths[light].read_bytes() == original
    assert fits.getheader(light)['CRCORR'] == 'Yes'
    assert cleaner.image_collection.files == [light.name]


def test_saturation_adu_uses_each_images_gain(tmp_path, fake_lacosmic):
    first = image(tmp_path, 'a.fits', GAIN=1.)
    second = image(tmp_path, 'b.fits', GAIN=2.)

    Artifacts(collection(tmp_path), overwrite=True).remove_cosmic_ray(
        saturation_adu=65000.,
    )

    assert {call[1]: call[3]['satlevel'] for call in fake_lacosmic} == {
        1.: 65000., 2.: 130000.,
    }
    for path in (first, second):
        assert CCDData.read(path).unit == u.adu
        assert fits.getdata(path)[0, 0] == 100.


@pytest.mark.parametrize('threshold', [0., -1., np.nan, np.inf, True, '65000'])
def test_invalid_saturation_adu_fails_before_writes(tmp_path, fake_lacosmic, threshold):
    light = image(tmp_path)
    before = light.read_bytes()

    with pytest.raises(ValueError, match='saturation_adu'):
        Artifacts(collection(tmp_path), overwrite=True).remove_cosmic_ray(
            saturation_adu=threshold,
        )

    assert not fake_lacosmic
    assert light.read_bytes() == before


def test_saturation_unit_conflict_and_legacy_satlevel(tmp_path, fake_lacosmic):
    light = image(tmp_path)
    before = light.read_bytes()
    cleaner = Artifacts(collection(tmp_path), overwrite=True)

    with pytest.raises(ValueError, match='not both'):
        cleaner.remove_cosmic_ray(saturation_adu=65000., satlevel=100000.)
    assert not fake_lacosmic
    assert light.read_bytes() == before

    cleaner.remove_cosmic_ray(satlevel=100000.)
    assert fake_lacosmic[0][3]['satlevel'] == 100000.


def test_single_file_saturation_conversion(tmp_path, fake_lacosmic):
    light = image(tmp_path, GAIN=2.)
    Artifacts._remove_cosmic_ray(light, overwrite=True, saturation_adu=65000.)
    assert fake_lacosmic[0][3]['satlevel'] == 130000.


def test_saturation_conversion_overflow_is_preflighted(tmp_path, fake_lacosmic):
    first = image(tmp_path, 'a.fits', GAIN=1.)
    second = image(tmp_path, 'b.fits', GAIN=1.e308)
    before = {path: path.read_bytes() for path in (first, second)}

    with pytest.raises(ValueError, match='multiplied by GAIN'):
        Artifacts(collection(tmp_path), overwrite=True).remove_cosmic_ray(
            saturation_adu=65000.,
        )

    assert not fake_lacosmic
    assert {path: path.read_bytes() for path in before} == before


@pytest.mark.parametrize('cards', [
    {'BIASCORR': 'No'}, {'DARKCORR': 'No'}, {'FLATCORR': 'No'},
    {'ALIGN': 'Registered'}, {'CRCORR': 'Yes'}, {'BUNIT': 'electron'},
    {'GAIN': 0.}, {'RDNOISE': -1.},
])
def test_batch_preflight_prevents_partial_cleaning(tmp_path, fake_lacosmic, cards):
    first = image(tmp_path, 'a.fits')
    invalid = image(tmp_path, 'z.fits', **cards)
    before = {path: path.read_bytes() for path in (first, invalid)}

    with pytest.raises(ValueError):
        Artifacts(collection(tmp_path), overwrite=True).remove_cosmic_ray()

    assert not fake_lacosmic
    assert {path: path.read_bytes() for path in before} == before


def test_in_place_requires_overwrite_and_repeat_is_rejected(tmp_path, fake_lacosmic):
    light = image(tmp_path)
    cleaner = Artifacts(collection(tmp_path))
    with pytest.raises(FileExistsError):
        cleaner.remove_cosmic_ray()
    assert not fake_lacosmic

    cleaner.remove_cosmic_ray(overwrite=True)
    before = light.read_bytes()
    with pytest.raises(ValueError, match='already Yes'):
        cleaner.remove_cosmic_ray(overwrite=True)
    assert light.read_bytes() == before
    assert len(fake_lacosmic) == 1


def test_filter_selection_and_calibration_selection_guard(tmp_path, fake_lacosmic):
    light = image(tmp_path)
    other = image(tmp_path, 'other.fits', FIELD='field_b')
    before = other.read_bytes()
    Artifacts(collection(tmp_path), overwrite=True).remove_cosmic_ray(
        headers_values={'IMTYPE': 'Light', 'FIELD': 'field_a'},
    )
    assert fits.getheader(light)['CRCORR'] == 'Yes'
    assert other.read_bytes() == before

    with pytest.raises(ValueError, match='selection must'):
        Artifacts(collection(tmp_path)).remove_cosmic_ray(headers_values={'IMTYPE': 'Flat'})


def test_unsupported_options_and_extensions_fail_before_writes(tmp_path):
    light = image(tmp_path)
    cleaner = Artifacts(collection(tmp_path), overwrite=True)
    before = light.read_bytes()
    with pytest.raises(TypeError):
        cleaner.remove_cosmic_ray(nonexistent_parameter=True)
    with pytest.raises(ValueError, match='controlled by headers'):
        cleaner.remove_cosmic_ray(gain_apply=True)
    assert light.read_bytes() == before

    with fits.open(light, mode='update') as hdul:
        hdul.append(fits.ImageHDU(np.ones((16, 16)), name='UNCERT'))
    before = light.read_bytes()
    with pytest.raises(ValueError, match='unsupported additional'):
        cleaner.remove_cosmic_ray()
    assert light.read_bytes() == before


def test_plot_mask_creates_only_binary_mask_png(tmp_path, fake_lacosmic, monkeypatch):
    light = image(tmp_path)
    plots = []
    original_plot = artifacts_module.plot_image

    def capture_plot(**kwargs):
        plots.append(kwargs)
        return original_plot(**kwargs)

    monkeypatch.setattr(artifacts_module, 'plot_image', capture_plot)
    Artifacts(collection(tmp_path), overwrite=True).remove_cosmic_ray(
        plot_mask=True,
    )
    assert list(tmp_path.glob('*.png')) == [tmp_path / f'{light.stem}_cr_mask.png']
    assert plots[0]['norm'].vmin == 0.
    assert plots[0]['norm'].vmax == 1.
    assert plots[0]['cmap'] == 'gray'
    assert set(np.unique(plots[0]['ccd_data'].data)) == {0., 1.}


@pytest.mark.parametrize('option', ['plot_before_after', 'norm_percent'])
def test_removed_plot_options_report_actionable_error(tmp_path, fake_lacosmic, option):
    light = image(tmp_path)
    before = light.read_bytes()
    with pytest.raises(TypeError, match='use DS9'):
        Artifacts(collection(tmp_path), overwrite=True).remove_cosmic_ray(**{option: False})
    assert not fake_lacosmic
    assert light.read_bytes() == before


def test_cmask_guides_detection_but_cleaned_values_reach_combine(tmp_path, fake_lacosmic, monkeypatch):
    from photozpy import Combine
    from skimage.transform import SimilarityTransform

    reference = image(tmp_path, 'a.fits')
    target = image(tmp_path, 'b.fits')
    cleaned = Artifacts(collection(tmp_path), overwrite=True).remove_cosmic_ray()
    calls = []

    def find(source, reference, **kwargs):
        assert source.mask[5, 6] and reference.mask[5, 6]
        calls.append(True)
        return SimilarityTransform(), (np.zeros((3, 2)), np.zeros((3, 2)))

    # Different repaired values prove the marked inputs contribute to the mean.
    with fits.open(target, mode='update') as hdul:
        hdul[0].data[5, 6] = 103.
    monkeypatch.setattr(registration_module, 'find_transform', find)
    registered = Registration(cleaned, overwrite=True).register_images("SDSS_g'")
    assert len(calls) == 1
    assert CCDData.read(reference).mask is None
    assert not CCDData.read(target).mask[5, 6]
    assert fits.getdata(reference, 'CMASK')[5, 6]
    assert fits.getdata(target, 'CMASK')[5, 6]
    Combine(registered).combine_light_or_flat(image_type='Light')
    master = CCDData.read(tmp_path / "Master_field_a_SDSS_g'.fits")
    assert master.data[5, 6] == pytest.approx(102.)
    assert not master.mask[5, 6]


def test_real_lacosmic_removes_impulse_without_gain_scaling_star(tmp_path):
    rng = np.random.default_rng(20)
    y, x = np.mgrid[:96, :96]
    data = 1000. + rng.normal(0., 5., (96, 96))
    data += 1500. * np.exp(-((x - 30.) ** 2 + (y - 32.) ** 2) / (2 * 2. ** 2))
    aperture = (x - 30.) ** 2 + (y - 32.) ** 2 < 7. ** 2
    original_flux = np.sum(data[aperture] - 1000.)
    data[65, 67] += 20000.
    light = image(tmp_path, data=data)

    Artifacts(collection(tmp_path), overwrite=True).remove_cosmic_ray()

    cleaned = CCDData.read(light)
    assert cleaned.unit == u.adu
    cmask = fits.getdata(light, 'CMASK').astype(bool)
    assert cleaned.mask is None
    assert cmask[65, 67]
    assert abs(cleaned.data[65, 67] - 1000.) < 50.
    assert not cmask[32, 30]
    assert np.median(cleaned.data) == pytest.approx(1000., abs=2.)
    assert np.sum(cleaned.data[aperture] - 1000.) / original_flux == pytest.approx(1., abs=0.01)
