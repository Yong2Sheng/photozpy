"""Registration checks on synthetic CCDs, including a real astroalign fit."""

from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.wcs import WCS
from ccdproc import ImageFileCollection
from scipy.ndimage import shift

from photozpy import Combine, ImageCollectionGroup, Registration
from photozpy.calibration import registration as registration_module
from photozpy.fits_utilities import has_wcs
from photozpy.telescope import SARA_CT


def image(directory, name, *, field='field_a', band="SDSS_g'", value=100., mask=None, **cards):
    directory.mkdir(parents=True, exist_ok=True)
    data = np.full((8, 8), value) if np.ndim(value) == 0 else np.asarray(value)
    header = fits.Header({
        'IMTYPE': 'Light', 'FIELD': field, 'FILTER': band,
        'BUNIT': 'ADU', 'OBJECT': 'unchanged_object',
    })
    header.update(cards)
    hdus = [fits.PrimaryHDU(data, header)]
    if mask is not None:
        hdus.append(fits.ImageHDU(np.asarray(mask, dtype=np.uint8), name='MASK'))
    path = directory / name
    fits.HDUList(hdus).writeto(path)
    return path


def collection(directory):
    return ImageFileCollection(location=directory, keywords='*', glob_exclude='.*')


@pytest.fixture
def fake_register(monkeypatch):
    calls = []

    def align(source, target, *, min_area, detection_sigma, propagate_mask):
        calls.append((source, target, min_area, detection_sigma, propagate_mask))
        assert source.data.dtype.isnative and target.data.dtype.isnative
        footprint = np.zeros(target.shape, dtype=bool)
        footprint[:, 0] = True
        if source.mask is not None:
            footprint |= source.mask
        return np.array(source.data, copy=True), footprint

    monkeypatch.setattr(registration_module, 'register', align)
    return calls


def test_field_grouping_cross_filter_alignment_and_retained_calibrations(tmp_path, fake_register):
    reference = image(tmp_path, 'a_g_001.fits', value=100.)
    target = image(tmp_path, 'a_r_001.fits', band="SDSS_r'", value=200.)
    reference_b = image(tmp_path, 'b_g_001.fits', field='field_b', value=300.)
    target_b = image(tmp_path, 'b_i_001.fits', field='field_b', band="SDSS_i'", value=400.)
    calibration = image(tmp_path, 'Master_Bias.fits', IMTYPE='Master Bias')
    before_calibration = calibration.read_bytes()
    reference_data = fits.getdata(reference).copy()

    processor = Registration(ImageCollectionGroup(tmp_path), SARA_CT, overwrite=True)
    result = processor.Register_images(filter='SDSS_g', min_area=7, detection_sigma=8)

    assert result is processor.image_collection
    assert isinstance(result, ImageCollectionGroup)
    assert set(processor.reference_paths) == {reference, reference_b}
    assert len(fake_register) == 2
    assert all(call[2:] == (7, 8, True) for call in fake_register)
    assert {float(call[1].data[0, 0]) for call in fake_register} == {100., 300.}
    assert calibration.read_bytes() == before_calibration
    np.testing.assert_array_equal(fits.getdata(reference), reference_data)

    for path, ref in ((target, reference), (target_b, reference_b)):
        header = fits.getheader(path)
        assert header['ALIGN'] == 'Registered'
        assert header['REGREF'] == ref.name
        assert len(header['REFHASH']) == 64
        assert header['OBJECT'] == 'unchanged_object'
        output = CCDData.read(path)
        assert output.mask[:, 0].all()
        assert not output.mask[:, 1:].any()
    assert fits.getheader(reference)['ALIGN'] == 'Reference'


def test_stable_reference_choice_and_explicit_reference(tmp_path, fake_register):
    a = image(tmp_path, 'g_001.fits')
    b = image(tmp_path, 'g_002.fits')
    files = ImageFileCollection(location=tmp_path, filenames=[b.name, a.name], keywords='*')
    processor = Registration(files, overwrite=True)
    plans = processor._plan_groups("SDSS_g'", None, None, True)
    assert plans[0].reference_image_path == a

    processor.register_images(reference_images={'field_a': b})
    assert processor.reference_paths == (b,)
    assert len(fake_register) == 1


def test_manual_references_per_field_and_progress_direction(tmp_path, fake_register, monkeypatch):
    target_a = image(tmp_path, 'a_g_001.fits', value=10.)
    reference_a = image(tmp_path, 'a_i_003.fits', band="SDSS_i'", value=20.)
    target_b = image(tmp_path, 'b_g_001.fits', field='field_b', value=30.)
    reference_b = image(tmp_path, 'b_z_009.fits', field='field_b', band="SDSS_z'", value=40.)
    reference_pixels = {path: fits.getdata(path).copy() for path in (reference_a, reference_b)}
    messages = []

    class Progress:
        def __init__(self, paths, **kwargs):
            self.paths = paths

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def __iter__(self):
            return iter(self.paths)

        def set_postfix_str(self, message):
            messages.append(message)

    monkeypatch.setattr(registration_module, 'tqdm', Progress)
    processor = Registration(ImageCollectionGroup(tmp_path), SARA_CT, overwrite=True)
    processor.register_images(
        reference_images={'field_a': reference_a, 'field_b': reference_b},
        detection_sigma=2, min_area=8,
    )

    assert set(processor.reference_paths) == {reference_a, reference_b}
    assert {float(call[1].data[0, 0]) for call in fake_register} == {20., 40.}
    assert all(call[2:] == (8, 2, True) for call in fake_register)
    assert messages == [
        f'Aligning {target_a.name} to {reference_a.name}',
        f'Aligning {target_b.name} to {reference_b.name}',
    ]
    for reference, pixels in reference_pixels.items():
        np.testing.assert_array_equal(fits.getdata(reference), pixels)

    fake_register.clear()
    messages.clear()
    processor.register_images(reference_images={'field_a': reference_a, 'field_b': reference_b})
    assert not fake_register
    assert all(message.startswith('Reusing ') for message in messages)


@pytest.mark.parametrize('problem', [
    'missing', 'wrong_field', 'backup', 'incomplete', 'unknown_field', 'bad_key',
])
def test_invalid_manual_selection_fails_before_any_writes(tmp_path, fake_register, problem):
    reference_a = image(tmp_path, 'a_g_001.fits')
    reference_b = image(tmp_path, 'b_g_001.fits', field='field_b')
    references = {'field_a': reference_a, 'field_b': reference_b}
    if problem == 'missing':
        references['field_b'] = tmp_path / 'missing.fits'
    elif problem == 'wrong_field':
        references['field_a'] = reference_b
    elif problem == 'backup':
        references['field_a'] = image(tmp_path / 'before_cosmic_ray', 'a_g_001.fits')
    elif problem == 'incomplete':
        del references['field_b']
    elif problem == 'unknown_field':
        references['typo_field'] = reference_a
    elif problem == 'bad_key':
        references[''] = reference_a

    before = {path: path.read_bytes() for path in (reference_a, reference_b)}
    with pytest.raises((ValueError, FileNotFoundError)):
        Registration(collection(tmp_path), overwrite=True).register_images(reference_images=references)
    assert not fake_register
    assert {path: path.read_bytes() for path in before} == before


def test_manual_selection_can_override_filter_shortcut(tmp_path, fake_register):
    first = image(tmp_path, 'a_g_001.fits')
    chosen = image(tmp_path, 'a_i_003.fits', band="SDSS_i'")
    processor = Registration(collection(tmp_path), SARA_CT, overwrite=True)
    processor.register_images(reference_filter='SDSS_g', reference_images={'field_a': chosen})
    assert processor.reference_paths == (chosen,)
    assert fits.getheader(first)['REGREF'] == chosen.name


def test_matching_registration_is_reused_without_resampling(tmp_path, fake_register):
    reference = image(tmp_path, 'g_001.fits')
    target = image(tmp_path, 'g_002.fits')
    processor = Registration(collection(tmp_path), overwrite=True)
    processor.register_images("SDSS_g'")
    before = {path: path.read_bytes() for path in (reference, target)}
    fake_register.clear()

    processor.register_images("SDSS_g'", overwrite=False)

    assert not fake_register
    assert {path: path.read_bytes() for path in before} == before


def test_changed_reference_pixels_are_rejected_before_writes(tmp_path, fake_register):
    reference = image(tmp_path, 'g_001.fits')
    target = image(tmp_path, 'g_002.fits')
    processor = Registration(collection(tmp_path), overwrite=True)
    processor.register_images("SDSS_g'")
    with fits.open(reference, mode='update') as hdul:
        hdul[0].data += 1
    before = target.read_bytes()
    with pytest.raises(ValueError, match='reference is missing or changed'):
        processor.register_images("SDSS_g'")
    assert target.read_bytes() == before


def test_changing_to_a_resampled_reference_is_rejected(tmp_path, fake_register):
    image(tmp_path, 'g_001.fits')
    target = image(tmp_path, 'g_002.fits')
    processor = Registration(collection(tmp_path), overwrite=True)
    processor.register_images("SDSS_g'")
    with pytest.raises(ValueError, match='new reference'):
        processor.register_images(reference_images={'field_a': target})


def test_separate_outputs_preserve_original_headers_and_pixels(tmp_path, fake_register):
    raw, output = tmp_path / 'analysis', tmp_path / 'separate'
    reference = image(raw, 'g_001.fits')
    target = image(raw, 'r_001.fits', band="SDSS_r'")
    originals = {path: path.read_bytes() for path in (reference, target)}
    processor = Registration(collection(raw))
    result = processor.register_images("SDSS_g'", save_location=output)
    assert {path: path.read_bytes() for path in originals} == originals
    assert Path(result.location) == output
    assert fits.getheader(output / reference.name)['ALIGN'] == 'Reference'
    assert fits.getheader(output / target.name)['ALIGN'] == 'Registered'


def test_original_files_need_explicit_overwrite(tmp_path, fake_register):
    reference = image(tmp_path, 'g_001.fits')
    target = image(tmp_path, 'g_002.fits')
    before = {path: path.read_bytes() for path in (reference, target)}
    with pytest.raises(FileExistsError):
        Registration(collection(tmp_path)).register_images("SDSS_g'")
    assert not fake_register
    assert {path: path.read_bytes() for path in before} == before


@pytest.mark.parametrize('problem', ['field', 'filter', 'legacy_marker', 'uncertainty', 'shape', 'wcs'])
def test_invalid_input_is_rejected_before_reference_is_marked(tmp_path, fake_register, problem):
    reference = image(tmp_path, 'g_001.fits')
    value = np.ones((9, 9)) if problem == 'shape' else 100.
    target = image(tmp_path, 'r_001.fits', band="SDSS_r'", value=value)
    with fits.open(target, mode='update') as hdul:
        if problem == 'field':
            del hdul[0].header['FIELD']
        elif problem == 'filter':
            del hdul[0].header['FILTER']
        elif problem == 'legacy_marker':
            hdul[0].header['ALIGN'] = 'Registered'
        elif problem == 'uncertainty':
            hdul.append(fits.ImageHDU(np.ones((8, 8)), name='UNCERT'))
        elif problem == 'wcs':
            hdul[0].header.update(celestial_header())
    before = {path: path.read_bytes() for path in (reference, target)}
    with pytest.raises(ValueError):
        Registration(collection(tmp_path), overwrite=True).register_images("SDSS_g'")
    assert not fake_register
    assert {path: path.read_bytes() for path in before} == before


def celestial_header():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    wcs.wcs.crval = [150, -30]
    wcs.wcs.crpix = [4, 4]
    wcs.wcs.cdelt = [-0.0001, 0.0001]
    return wcs.to_header()


def test_unresampled_reference_can_keep_its_wcs(tmp_path, fake_register):
    reference = image(tmp_path, 'g_001.fits')
    target = image(tmp_path, 'r_001.fits', band="SDSS_r'")
    with fits.open(reference, mode='update') as hdul:
        hdul[0].header.update(celestial_header())
    Registration(collection(tmp_path), overwrite=True).register_images("SDSS_g'")
    assert has_wcs(reference)
    assert not has_wcs(target)


def test_missing_reference_band_fails_before_any_write(tmp_path, fake_register):
    image(tmp_path, 'r_001.fits', band="SDSS_r'")
    with pytest.raises(ValueError, match='no reference image'):
        Registration(collection(tmp_path), overwrite=True).register_images("SDSS_g'")
    assert not fake_register


def test_groups_with_repeated_field_names_use_local_references(tmp_path, fake_register):
    one, two = tmp_path / 'night_one', tmp_path / 'night_two'
    for root, value in ((one, 100.), (two, 200.)):
        image(root, 'g_001.fits', value=value)
        image(root, 'r_001.fits', band="SDSS_r'", value=value)
    processor = Registration(ImageCollectionGroup((one, two)), overwrite=True)
    processor.register_images("SDSS_g'")
    assert {path.parent for path in processor.reference_paths} == {one, two}
    assert {float(call[1].data[0, 0]) for call in fake_register} == {100., 200.}


def test_shared_output_collisions_are_rejected(tmp_path, fake_register):
    one, two = tmp_path / 'night_one', tmp_path / 'night_two'
    image(one, 'g_001.fits')
    image(two, 'g_001.fits')
    out = tmp_path / 'out'
    with pytest.raises(ValueError, match='same output'):
        Registration(ImageCollectionGroup((one, two))).register_images(
            "SDSS_g'", save_location=out,
        )
    assert not out.exists()


def test_failed_transform_does_not_mark_or_overwrite_target(tmp_path, monkeypatch):
    image(tmp_path, 'g_001.fits')
    target = image(tmp_path, 'r_001.fits', band="SDSS_r'")
    before = target.read_bytes()

    def fail(*args, **kwargs):
        raise ValueError('no matching stars')

    monkeypatch.setattr(registration_module, 'register', fail)
    with pytest.raises(RuntimeError, match='registration to'):
        Registration(collection(tmp_path), overwrite=True).register_images("SDSS_g'")
    assert target.read_bytes() == before
    assert 'ALIGN' not in fits.getheader(target)


def test_registration_mask_survives_light_combination(tmp_path, fake_register):
    mask = np.zeros((8, 8), dtype=bool)
    mask[:, 0] = True
    image(tmp_path, 'g_001.fits', value=10., mask=mask)
    image(tmp_path, 'g_002.fits', value=30., mask=mask)
    working = Registration(collection(tmp_path), overwrite=True).register_images("SDSS_g'")
    combiner = Combine(working)
    combiner.combine_light_or_flat('Light', quite=True)
    master = CCDData.read(combiner.master_paths[0])
    assert master.mask[:, 0].all()
    assert not master.mask[:, 1:].any()
    np.testing.assert_allclose(master.data[:, 1:], 20.)


def test_light_combine_rejects_unregistered_or_different_grids(tmp_path):
    a = image(tmp_path, 'g_001.fits')
    b = image(tmp_path, 'g_002.fits')
    with pytest.raises(ValueError, match='register Light'):
        Combine(collection(tmp_path)).combine_light_or_flat('Light', quite=True)
    for path, digest in ((a, 'one'), (b, 'two')):
        with fits.open(path, mode='update') as hdul:
            hdul[0].header['ALIGN'] = 'Registered'
            hdul[0].header['REGREF'] = 'reference.fits'
            hdul[0].header['REFHASH'] = digest
    with pytest.raises(ValueError, match='different registration'):
        Combine(collection(tmp_path)).combine_light_or_flat('Light', quite=True)


def synthetic_star_field():
    rng = np.random.default_rng(741)
    yy, xx = np.mgrid[:128, :128]
    data = rng.normal(100., 0.4, size=(128, 128))
    stars = [
        (19., 21., 1100.), (45., 18., 1700.), (80., 20., 900.),
        (107., 17., 1300.), (25., 48., 800.), (63., 41., 1600.),
        (99., 46., 1900.), (16., 80., 1400.), (48., 75., 2000.),
        (78., 79., 1000.), (111., 87., 1500.), (31., 106., 2100.),
        (62., 111., 1200.), (94., 111., 1800.),
    ]
    for x, y, peak in stars:
        data += peak * np.exp(-((xx - x)**2 + (yy - y)**2) / (2 * 1.7**2))
    return data, stars


@pytest.mark.parametrize('with_cmask', [False, True])
def test_real_astroalign_recovers_subpixel_shift_and_aperture_flux(tmp_path, with_cmask):
    """Exercise the installed engine, not a mock or transformation oracle."""
    data, stars = synthetic_star_field()
    reference = image(tmp_path, 'g_001.fits', value=data)
    translated = shift(data, shift=(2.4, -3.1), order=3, mode='constant', cval=100.)
    target = image(tmp_path, 'r_001.fits', band="SDSS_r'", value=translated)
    if with_cmask:
        cmask = np.zeros(data.shape, dtype=np.uint8)
        cmask[60, 60] = 1
        for path in (reference, target):
            with fits.open(path, mode='update') as hdul:
                hdul.append(fits.ImageHDU(cmask, name='CMASK'))
    original_reference = fits.getdata(reference).copy()

    Registration(collection(tmp_path), overwrite=True).register_images("SDSS_g'")

    aligned = CCDData.read(target)
    if with_cmask:
        transformed_cmask = fits.getdata(target, 'CMASK').astype(bool)
        assert transformed_cmask[58, 63]
        assert not transformed_cmask[60, 60]
        assert not aligned.mask[transformed_cmask].any()
        assert fits.getdata(reference, 'CMASK')[60, 60]
    np.testing.assert_array_equal(fits.getdata(reference), original_reference)
    assert aligned.mask.any()
    yy, xx = np.mgrid[:128, :128]
    for x, y, _ in stars:
        aperture = (xx - x)**2 + (yy - y)**2 <= 6**2
        assert not aligned.mask[aperture].any()
        signal = np.clip(aligned.data - 100., 0, None)
        flux = np.sum(signal[aperture])
        reference_flux = np.sum((data - 100.)[aperture])
        assert flux / reference_flux == pytest.approx(1.0, abs=0.03)
        center_x = np.sum(xx[aperture] * signal[aperture]) / flux
        center_y = np.sum(yy[aperture] * signal[aperture]) / flux
        assert abs(center_x - x) < 0.15
        assert abs(center_y - y) < 0.15


def test_legacy_cosmic_mask_rejected_before_registration_or_combination(tmp_path, fake_register):
    """Do not silently reuse old outputs whose two mask meanings are mixed."""
    mask = np.zeros((8, 8), dtype=bool)
    mask[3, 3] = True
    reference = image(tmp_path, 'a.fits', mask=mask, CRCORR='Yes')
    target = image(tmp_path, 'b.fits', mask=mask, CRCORR='Yes')
    originals = {path: path.read_bytes() for path in (reference, target)}
    with pytest.raises(ValueError, match='legacy cosmic-ray MASK'):
        Registration(collection(tmp_path), overwrite=True).register_images("SDSS_g'")
    with pytest.raises(ValueError, match='legacy cosmic-ray MASK'):
        Combine(collection(tmp_path)).combine_light_or_flat('Light')
    assert {path: path.read_bytes() for path in originals} == originals
    assert not fake_register


def test_reusing_registration_does_not_reject_a_later_valid_wcs(tmp_path, fake_register):
    image(tmp_path, 'g_001.fits')
    target = image(tmp_path, 'g_002.fits')
    processor = Registration(collection(tmp_path), overwrite=True)
    processor.register_images("SDSS_g'")
    with fits.open(target, mode='update') as hdul:
        hdul[0].header.update(celestial_header())
    before = target.read_bytes()
    fake_register.clear()
    processor.register_images("SDSS_g'", overwrite=False)
    assert not fake_register
    assert target.read_bytes() == before
