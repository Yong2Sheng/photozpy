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


@pytest.mark.parametrize('with_cmask', [False, True])
@pytest.mark.parametrize('entry_point', ['register_images', 'Register_images'])
def test_sigma_candidates_retry_in_order_per_target_and_reuse(
    tmp_path, monkeypatch, with_cmask, entry_point, caplog,
):
    image(tmp_path, 'a.fits')
    first = image(tmp_path, 'b.fits', value=200.)
    second = image(tmp_path, 'c.fits', value=300.)
    if with_cmask:
        for path in tmp_path.glob('*.fits'):
            with fits.open(path, mode='update') as hdul:
                hdul.append(fits.ImageHDU(np.zeros((8, 8), dtype=np.uint8), name='CMASK'))
    calls = []

    def align(source, reference, *masks, min_area, detection_sigma, **kwargs):
        value = int(source.data[0, 0])
        calls.append((value, detection_sigma))
        assert min_area == 8
        if value == 200 and detection_sigma != 7:
            raise registration_module.MaxIterError('no matching triangles')
        result = (source.data.copy(), np.zeros(source.shape, dtype=bool))
        return (*result, np.zeros(source.shape, dtype=bool)) if with_cmask else result

    monkeypatch.setattr(
        registration_module.Registration if with_cmask else registration_module,
        '_align_with_cmask' if with_cmask else 'register', staticmethod(align) if with_cmask else align,
    )
    processor = Registration(collection(tmp_path), overwrite=True)
    with caplog.at_level('INFO', logger=registration_module.__name__):
        getattr(processor, entry_point)(
            "SDSS_g'", min_area=8, detection_sigmas=[5, 6, 7, 8],
        )
    assert calls == [(200, 5), (200, 6), (200, 7), (300, 5)]
    assert fits.getheader(first)['REGSIG'] == 7
    assert fits.getheader(second)['REGSIG'] == 5
    assert 'matched at detection_sigma=7' in caplog.text
    before = {path: path.read_bytes() for path in tmp_path.glob('*.fits')}
    calls.clear()
    processor.register_images("SDSS_g'", detection_sigmas=[10, 12], overwrite=False)
    assert not calls
    assert all(path.read_bytes() == data for path, data in before.items())


@pytest.mark.parametrize('failure', [
    registration_module.MaxIterError('no matching triangles'),
    ValueError('Reference stars in source image are less than the minimum value (3).'),
])
def test_sigma_candidates_exhaustion_preserves_target(tmp_path, monkeypatch, failure):
    image(tmp_path, 'a.fits')
    target = image(tmp_path, 'b.fits')
    before = target.read_bytes()
    calls = []

    def fail(*args, detection_sigma, **kwargs):
        calls.append(detection_sigma)
        raise failure

    monkeypatch.setattr(registration_module, 'register', fail)
    with pytest.raises(RuntimeError, match=r'attempted detection_sigmas=\[5, 6, 7\]') as error:
        Registration(collection(tmp_path), overwrite=True).register_images(
            "SDSS_g'", detection_sigmas=[5, 6, 7],
        )
    assert error.value.__cause__ is failure
    assert calls == [5, 6, 7]
    assert target.read_bytes() == before


@pytest.mark.parametrize('failure', [OSError('read error'), ValueError('invalid output shape')])
def test_sigma_candidates_do_not_retry_unrelated_errors(tmp_path, monkeypatch, failure):
    image(tmp_path, 'a.fits')
    target = image(tmp_path, 'b.fits')
    before = target.read_bytes()
    calls = []

    def fail(*args, detection_sigma, **kwargs):
        calls.append(detection_sigma)
        raise failure

    monkeypatch.setattr(registration_module, 'register', fail)
    with pytest.raises(RuntimeError) as error:
        Registration(collection(tmp_path), overwrite=True).register_images(
            "SDSS_g'", detection_sigmas=[5, 6, 7],
        )
    assert error.value.__cause__ is failure
    assert calls == [5]
    assert target.read_bytes() == before


@pytest.mark.parametrize('sigmas', [[], '5,12', 5, [0, 5], [-1], [float('nan')],
                                    [float('inf')], [True], ['5'], [None]])
def test_invalid_sigma_candidates_fail_before_writes(tmp_path, fake_register, sigmas):
    reference = image(tmp_path, 'a.fits')
    target = image(tmp_path, 'b.fits')
    before = {p: p.read_bytes() for p in (reference, target)}
    with pytest.raises(ValueError, match='detection_sigmas'):
        Registration(collection(tmp_path), overwrite=True).register_images(
            "SDSS_g'", detection_sigmas=sigmas,
        )
    assert not fake_register
    assert all(p.read_bytes() == data for p, data in before.items())


def test_sigma_candidates_preserve_explicit_order_and_override_single_sigma(tmp_path, monkeypatch):
    image(tmp_path, 'a.fits')
    target = image(tmp_path, 'b.fits')
    calls = []

    def align(source, reference, *, detection_sigma, **kwargs):
        calls.append(detection_sigma)
        if detection_sigma == 8.5:
            raise registration_module.MaxIterError('no matching triangles')
        return source.data.copy(), np.zeros(source.shape, dtype=bool)

    monkeypatch.setattr(registration_module, 'register', align)
    Registration(collection(tmp_path), overwrite=True).register_images(
        "SDSS_g'", detection_sigma=99, detection_sigmas=[8.5, 5.5, 7.5],
    )
    assert calls == [8.5, 5.5]
    assert fits.getheader(target)['REGSIG'] == 5.5


@pytest.mark.parametrize('entrypoint', ['register_images', 'Register_images'])
def test_auto_references_follow_each_fields_collection_order(tmp_path, fake_register, entrypoint):
    # Deliberately put a calibration first and choose a non-alphabetical order.
    bias = image(tmp_path, 'bias.fits', IMTYPE='Master Bias')
    a_g = image(tmp_path, 'a_g.fits')
    a_i = image(tmp_path, 'a_i.fits', band="SDSS_i'", value=200.)
    b_g = image(tmp_path, 'b_g.fits', field='field_b')
    b_z = image(tmp_path, 'b_z.fits', field='field_b', band="SDSS_z'", value=300.)
    files = ImageFileCollection(
        location=tmp_path,
        filenames=[p.name for p in (bias, b_z, a_i, b_g, a_g)], keywords='*',
    )
    processor = Registration(files, SARA_CT, overwrite=True)
    getattr(processor, entrypoint)(reference_images='auto')
    assert processor.reference_paths == (b_z, a_i)
    assert len(fake_register) == 2
    assert fits.getheader(a_g)['REGREF'] == a_i.name
    assert fits.getheader(b_g)['REGREF'] == b_z.name
    assert 'ALIGN' not in fits.getheader(bias)

    # Repeating auto with the same order reuses targets without resampling.
    fake_register.clear()
    processor.register_images(reference_images='auto')
    assert not fake_register


def test_auto_references_respect_child_boundaries(tmp_path, fake_register):
    first = image(tmp_path / 'one', 'a.fits')
    image(tmp_path / 'one', 'b.fits', band="SDSS_r'")
    second = image(tmp_path / 'two', 'a.fits', band="SDSS_i'")
    image(tmp_path / 'two', 'b.fits', band="SDSS_z'")
    processor = Registration(ImageCollectionGroup([first.parent, second.parent]), overwrite=True)
    processor.register_images(reference_images='auto')
    assert set(processor.reference_paths) == {first, second}
    assert len(fake_register) == 2


def test_auto_reordered_registered_images_fail_before_writes(tmp_path, fake_register):
    first = image(tmp_path, 'a.fits')
    second = image(tmp_path, 'b.fits')
    Registration(collection(tmp_path), overwrite=True).register_images(reference_images='auto')
    before = {p: p.read_bytes() for p in (first, second)}
    reordered = ImageFileCollection(
        location=tmp_path, filenames=[second.name, first.name], keywords='*',
    )
    fake_register.clear()
    with pytest.raises(ValueError, match='resampled target'):
        Registration(reordered, overwrite=True).register_images(reference_images='auto')
    assert not fake_register
    assert {p: p.read_bytes() for p in before} == before


@pytest.mark.parametrize('setting', ['automatic', '', True])
def test_invalid_reference_mode_fails_before_writes(tmp_path, fake_register, setting):
    first = image(tmp_path, 'a.fits')
    before = first.read_bytes()
    with pytest.raises(TypeError, match='reference_images'):
        Registration(collection(tmp_path), overwrite=True).register_images(reference_images=setting)
    assert first.read_bytes() == before
    assert not fake_register


@pytest.mark.parametrize('entrypoint', ['register_images', 'Register_images'])
def test_discard_input_wcs_only_on_newly_resampled_outputs(tmp_path, fake_register, entrypoint):
    reference = image(tmp_path / 'input', 'a.fits')
    target = image(tmp_path / 'input', 'b.fits', band="SDSS_r'")
    for path in (reference, target):
        with fits.open(path, mode='update') as hdul:
            hdul[0].header.update(celestial_header())
            hdul[0].header['CTYPE1A'] = 'RA---TAN'
            hdul[0].header['CRPIX1A'] = 4.
            hdul[0].header['A_ORDER'] = 2
            hdul[0].header['B_ORDER'] = 2
            hdul[0].header['B_0_2'] = 1e-5
            hdul[0].header['A_2_0'] = 1e-5
    before = {p: p.read_bytes() for p in (reference, target)}
    output = tmp_path / 'output'
    processor = Registration(collection(reference.parent), overwrite=True)
    getattr(processor, entrypoint)(
        reference_images='auto', discard_input_wcs=True, save_location=output,
    )
    assert {p: p.read_bytes() for p in before} == before
    assert has_wcs(output / reference.name)
    h = fits.getheader(output / target.name)
    assert not has_wcs(output / target.name)
    assert not any(k.startswith(('CTYPE', 'CRPIX', 'CRVAL', 'CDELT', 'PC', 'A_')) for k in h)
    assert h['FILTER'] == "SDSS_r'"
    assert h['ALIGN'] == 'Registered'
    assert h['REGREF'] == reference.name
    assert 'removed input spatial WCS' in str(h['HISTORY'])
    # Newly solved WCS survives a repeated run which reuses the existing grid.
    with fits.open(output / target.name, mode='update') as hdul:
        hdul[0].header.update(celestial_header())
    fake_register.clear()
    Registration(collection(output), overwrite=True).register_images(
        reference_images='auto', discard_input_wcs=True,
    )
    assert not fake_register
    assert has_wcs(output / target.name)


@pytest.mark.parametrize("entry_point", ["register_images", "Register_images"])
def test_field_smoothing_fallback_restarts_earlier_success_from_native_pixels(tmp_path, monkeypatch, entry_point):
    from skimage.transform import SimilarityTransform

    data = np.zeros((8, 8))
    data[3, 3] = 100.
    reference = image(tmp_path, 'a.fits', value=data)
    first = image(tmp_path, 'b.fits', value=data + 20)
    last = image(tmp_path, 'c.fits', value=data + 40)
    originals = {p: p.read_bytes() for p in (reference, first, last)}
    calls = []

    def ordinary(source, target, **kwargs):
        calls.append(('ordinary', source.header['FILENAME'] if 'FILENAME' in source.header else
                      float(source.data[0, 0]), kwargs['detection_sigma']))
        # No reference/target has been published, even after b succeeded.
        assert all(p.read_bytes() == before for p, before in originals.items())
        if source.data[0, 0] == 40:
            raise registration_module.MaxIterError('no triangles')
        return source.data + 1000, np.zeros(source.shape, bool)

    def smoothed(source, target, **kwargs):
        calls.append(('smoothed', kwargs['detection_sigma']))
        assert source.data[3, 3] < 100
        assert all(p.read_bytes() == before for p, before in originals.items())
        return SimilarityTransform(), (np.zeros((5, 2)), np.zeros((5, 2)))

    def apply(transform, source, target, **kwargs):
        # The actual resampling sees the original sharp source, not ordinary
        # pass output (+1000) or the smoothed detection image.
        assert source.data[0, 0] in (20., 40.)
        assert source.data[3, 3] == source.data[0, 0] + 100
        return source.data.copy(), np.zeros(source.shape, bool)

    monkeypatch.setattr(registration_module, 'register', ordinary)
    monkeypatch.setattr(registration_module, 'find_transform', smoothed)
    monkeypatch.setattr(registration_module, 'apply_transform', apply)
    processor = Registration(collection(tmp_path), overwrite=True)
    getattr(processor, entry_point)("SDSS_g'", detection_sigmas=[3, 4], fallback_smoothing_sigma=2.)
    assert calls == [('ordinary', 20., 3), ('ordinary', 40., 3), ('ordinary', 40., 4),
                     ('smoothed', 3), ('smoothed', 3)]
    for path in (reference, first, last):
        assert fits.getheader(path)['REGSMTH'] == 2.
    np.testing.assert_array_equal(fits.getdata(first), data + 20)
    np.testing.assert_array_equal(fits.getdata(last), data + 40)
    assert fits.getheader(first)['REGSIG'] == 3
    # A repeat call uses the complete, consistent field without a new fit.
    processor.register_images("SDSS_g'", fallback_smoothing_sigma=2., overwrite=False)
    assert len(calls) == 5


def test_failed_smoothing_leaves_entire_field_unchanged(tmp_path, monkeypatch):
    a = image(tmp_path, 'a.fits')
    b = image(tmp_path, 'b.fits', value=200.)
    before = {p: p.read_bytes() for p in (a, b)}

    def fail(*args, **kwargs):
        raise registration_module.MaxIterError('no triangles')

    monkeypatch.setattr(registration_module, 'register', fail)
    monkeypatch.setattr(registration_module, 'find_transform', fail)
    with pytest.raises(RuntimeError):
        Registration(collection(tmp_path), overwrite=True).register_images(
            "SDSS_g'", detection_sigmas=[3, 4], fallback_smoothing_sigma=2.)
    assert all(p.read_bytes() == value for p, value in before.items())


@pytest.mark.parametrize('setting', [0, -1, float('nan'), float('inf'), True, '2'])
def test_invalid_smoothing_setting_rejected_without_writes(tmp_path, setting):
    a = image(tmp_path, 'a.fits')
    before = a.read_bytes()
    with pytest.raises(ValueError, match='fallback_smoothing_sigma'):
        Registration(collection(tmp_path), overwrite=True).register_images(
            "SDSS_g'", fallback_smoothing_sigma=setting)
    assert a.read_bytes() == before


def test_fallback_ordinary_success_and_mixed_resume_guard(tmp_path, fake_register):
    a = image(tmp_path, 'a.fits')
    b = image(tmp_path, 'b.fits', value=200.)
    processor = Registration(collection(tmp_path), overwrite=True)
    processor.register_images("SDSS_g'", fallback_smoothing_sigma=2.)
    assert all(fits.getheader(p)['REGSMTH'] == 0 for p in (a, b))
    c = image(tmp_path, 'c.fits', value=300.)
    before = {p: p.read_bytes() for p in (a, b, c)}
    with pytest.raises(ValueError, match='unresampled inputs'):
        Registration(collection(tmp_path), overwrite=True).register_images(
            "SDSS_g'", fallback_smoothing_sigma=2.)
    assert all(p.read_bytes() == value for p, value in before.items())


def test_smoothing_ignores_masked_spike_and_keeps_science_and_cmask_separate(monkeypatch):
    from skimage.transform import SimilarityTransform
    data, _ = synthetic_star_field()
    mask = np.zeros(data.shape, bool)
    mask[60, 60] = True
    data[60, 60] = 1e10
    source = CCDData(data.copy(), unit='adu')
    reference = CCDData(data.copy(), unit='adu')

    def find(s, r, **kwargs):
        assert s.mask[60, 60]
        assert s.data[60, 60] < 1e5
        assert s.data[60, 61] < 1e5
        return SimilarityTransform(), (None, None)

    def apply(t, s, r, **kwargs):
        assert s.data[60, 60] == 1e10
        assert s.mask is None
        return s.data.copy(), np.zeros(s.shape, bool)

    monkeypatch.setattr(registration_module, 'find_transform', find)
    monkeypatch.setattr(registration_module, 'apply_transform', apply)
    aligned, footprint, cmask = Registration._align_with_cmask(
        source, reference, mask, mask, min_area=8, detection_sigma=3, smoothing_sigma=2.)
    np.testing.assert_array_equal(aligned, data)
    assert not footprint.any()
    np.testing.assert_array_equal(cmask, mask)


def test_real_smoothed_detection_resamples_original_flux_once(tmp_path, monkeypatch):
    data, stars = synthetic_star_field()
    reference = image(tmp_path, 'g_001.fits', value=data)
    translated = shift(data, shift=(2.4, -3.1), order=3, mode='constant', cval=100.)
    target = image(tmp_path, 'r_001.fits', band="SDSS_r'", value=translated)

    def fail_raw(*args, **kwargs):
        raise registration_module.MaxIterError('exercise smoothing fallback')

    monkeypatch.setattr(registration_module, 'register', fail_raw)
    Registration(collection(tmp_path), overwrite=True).register_images(
        "SDSS_g'", fallback_smoothing_sigma=2.)
    aligned = CCDData.read(target)
    yy, xx = np.mgrid[:128, :128]
    for x, y, _ in stars:
        aperture = (xx - x)**2 + (yy - y)**2 <= 6**2
        flux = np.sum((aligned.data - 100.)[aperture])
        assert flux / np.sum((data - 100.)[aperture]) == pytest.approx(1., abs=.03)
        cx = np.sum(xx[aperture] * (aligned.data - 100.)[aperture]) / flux
        cy = np.sum(yy[aperture] * (aligned.data - 100.)[aperture]) / flux
        assert abs(cx - x) < .15 and abs(cy - y) < .15
    np.testing.assert_array_equal(fits.getdata(reference), data)
    assert fits.getheader(target)['REGSMTH'] == 2.


def test_fallback_does_not_retry_processing_errors(tmp_path, monkeypatch):
    a = image(tmp_path, 'a.fits')
    b = image(tmp_path, 'b.fits')
    before = {p: p.read_bytes() for p in (a, b)}

    def fail(*args, **kwargs):
        raise OSError('disk read failed')

    def forbidden(*args, **kwargs):
        pytest.fail('processing errors must not trigger smooth matching')
    monkeypatch.setattr(registration_module, 'register', fail)
    monkeypatch.setattr(registration_module, 'find_transform', forbidden)
    with pytest.raises(RuntimeError) as error:
        Registration(collection(tmp_path), overwrite=True).register_images(
            "SDSS_g'", fallback_smoothing_sigma=2.)
    assert isinstance(error.value.__cause__, OSError)
    assert all(p.read_bytes() == data for p, data in before.items())


# Review regressions: recorded field mode must remain authoritative even when
# the caller omits the fallback option on a later run.
def _complete_smoothed_field(tmp_path, monkeypatch):
    reference = image(tmp_path, 'a.fits')
    target = image(tmp_path, 'b.fits', value=200.)

    def ordinary_failure(*args, **kwargs):
        raise registration_module.MaxIterError('force field fallback')

    def smoothed(source, reference, *masks, **kwargs):
        return source.data.copy(), np.zeros(source.shape, bool), None

    monkeypatch.setattr(registration_module, 'register', ordinary_failure)
    monkeypatch.setattr(Registration, '_align_with_cmask', staticmethod(smoothed))
    Registration(collection(tmp_path), overwrite=True).register_images(
        "SDSS_g'", fallback_smoothing_sigma=2.)
    return reference, target


def test_default_rerun_rejects_new_native_target_in_recorded_smoothed_field(tmp_path, monkeypatch):
    reference, target = _complete_smoothed_field(tmp_path, monkeypatch)
    added = image(tmp_path, 'c.fits', value=300.)
    before = {p: p.read_bytes() for p in (reference, target, added)}

    def ordinary(source, reference, **kwargs):
        return source.data.copy(), np.zeros(source.shape, bool)

    monkeypatch.setattr(registration_module, 'register', ordinary)
    with pytest.raises(ValueError, match='REGSMTH'):
        Registration(collection(tmp_path), overwrite=True).register_images("SDSS_g'")
    assert all(p.read_bytes() == value for p, value in before.items())


def test_single_light_fallback_repeat_is_read_only_with_overwrite_false(tmp_path, fake_register):
    reference = image(tmp_path, 'a.fits')
    processor = Registration(collection(tmp_path), overwrite=True)
    processor.register_images("SDSS_g'", fallback_smoothing_sigma=2.)
    before = reference.read_bytes()
    processor.register_images("SDSS_g'", fallback_smoothing_sigma=2., overwrite=False)
    assert reference.read_bytes() == before
    assert not fake_register


@pytest.mark.parametrize('fallback', [None, 2.])
def test_cleared_align_with_registration_provenance_cannot_resample_again(tmp_path, monkeypatch, fallback):
    reference, target = _complete_smoothed_field(tmp_path, monkeypatch)
    with fits.open(target, mode='update') as hdul:
        del hdul[0].header['ALIGN']
    before = {p: p.read_bytes() for p in (reference, target)}
    with pytest.raises(ValueError, match='registration markers'):
        Registration(collection(tmp_path), overwrite=True).register_images(
            "SDSS_g'", fallback_smoothing_sigma=fallback)
    assert all(p.read_bytes() == value for p, value in before.items())


@pytest.mark.parametrize('fallback', [None, 2.])
@pytest.mark.parametrize('overwrite', [False, True])
def test_reference_only_selection_preserves_recorded_smoothed_mode(tmp_path, monkeypatch, fallback, overwrite):
    reference, _ = _complete_smoothed_field(tmp_path, monkeypatch)
    before = reference.read_bytes()
    selected = ImageFileCollection(tmp_path, filenames=[reference.name], keywords='*')
    Registration(selected, overwrite=overwrite).register_images(
        "SDSS_g'", fallback_smoothing_sigma=fallback)
    assert reference.read_bytes() == before
    assert fits.getheader(reference)['REGSMTH'] == 2.


def test_default_complete_recorded_field_is_reused_without_fitting(tmp_path, monkeypatch):
    reference, target = _complete_smoothed_field(tmp_path, monkeypatch)
    before = {p: p.read_bytes() for p in (reference, target)}

    def forbidden(*args, **kwargs):
        pytest.fail('a complete field must not be fitted again')
    monkeypatch.setattr(registration_module, 'register', forbidden)
    monkeypatch.setattr(Registration, '_align_with_cmask', staticmethod(forbidden))
    Registration(collection(tmp_path)).register_images("SDSS_g'")
    assert all(p.read_bytes() == value for p, value in before.items())


@pytest.mark.parametrize('mode', [None, 0., -1., 'unknown', True])
@pytest.mark.parametrize('fallback', [None, 2.])
def test_mixed_or_invalid_recorded_modes_fail_on_every_entry_path(tmp_path, monkeypatch, mode, fallback):
    reference, target = _complete_smoothed_field(tmp_path, monkeypatch)
    with fits.open(target, mode='update') as hdul:
        if mode is None:
            del hdul[0].header['REGSMTH']
        else:
            hdul[0].header['REGSMTH'] = mode
    before = {p: p.read_bytes() for p in (reference, target)}
    with pytest.raises(ValueError, match='REGSMTH'):
        Registration(collection(tmp_path), overwrite=True).register_images(
            "SDSS_g'", fallback_smoothing_sigma=fallback)
    assert all(p.read_bytes() == value for p, value in before.items())


def test_changing_recorded_smoothing_sigma_requires_native_copies(tmp_path, monkeypatch):
    reference, target = _complete_smoothed_field(tmp_path, monkeypatch)
    before = {p: p.read_bytes() for p in (reference, target)}
    with pytest.raises(ValueError, match='REGSMTH'):
        Registration(collection(tmp_path), overwrite=True).register_images(
            "SDSS_g'", fallback_smoothing_sigma=3.)
    assert all(p.read_bytes() == value for p, value in before.items())


def test_legacy_complete_field_can_default_reuse_but_not_opt_into_fallback(tmp_path, fake_register):
    reference = image(tmp_path, 'a.fits')
    target = image(tmp_path, 'b.fits')
    Registration(collection(tmp_path), overwrite=True).register_images("SDSS_g'")
    before = {p: p.read_bytes() for p in (reference, target)}
    fake_register.clear()
    Registration(collection(tmp_path)).register_images("SDSS_g'")
    assert not fake_register
    with pytest.raises(ValueError, match='REGSMTH'):
        Registration(collection(tmp_path), overwrite=True).register_images(
            "SDSS_g'", fallback_smoothing_sigma=2.)
    assert all(p.read_bytes() == value for p, value in before.items())


@pytest.mark.parametrize('marker,value', [('REGSMTH', 2.), ('REGSIG', 3),
                                          ('REGREF', 'a.fits'), ('REFHASH', 'some-grid')])
@pytest.mark.parametrize('align', [None, ''])
def test_any_orphan_registration_marker_is_rejected(tmp_path, fake_register, marker, value, align):
    reference = image(tmp_path, 'a.fits')
    target = image(tmp_path, 'b.fits', **{marker: value})
    if align is not None:
        with fits.open(target, mode='update') as hdul:
            hdul[0].header['ALIGN'] = align
    before = {p: p.read_bytes() for p in (reference, target)}
    with pytest.raises(ValueError, match='registration markers'):
        Registration(collection(tmp_path), overwrite=True).register_images("SDSS_g'")
    assert not fake_register
    assert all(p.read_bytes() == data for p, data in before.items())


def test_partial_publication_has_context_and_default_retry_is_blocked(tmp_path, monkeypatch):
    reference = image(tmp_path, 'a.fits')
    first = image(tmp_path, 'b.fits', value=200.)
    last = image(tmp_path, 'c.fits', value=300.)

    def fail(*args, **kwargs):
        raise registration_module.MaxIterError('force smooth mode')

    def smooth(source, reference, *masks, **kwargs):
        return source.data.copy(), np.zeros(source.shape, bool), None
    monkeypatch.setattr(registration_module, 'register', fail)
    monkeypatch.setattr(Registration, '_align_with_cmask', staticmethod(smooth))
    original_update = registration_module.update_fits_header

    def publish(path, cards, **kwargs):
        if not cards and kwargs.get('output_path') == last:
            raise OSError('simulated publication failure')
        return original_update(path, cards, **kwargs)
    monkeypatch.setattr(registration_module, 'update_fits_header', publish)
    last_before = last.read_bytes()
    with pytest.raises(OSError, match='FIELD=field_a: publication failed.*after 2/3'):
        Registration(collection(tmp_path), overwrite=True).register_images(
            "SDSS_g'", fallback_smoothing_sigma=2.)
    assert fits.getheader(reference)['REGSMTH'] == 2.
    assert fits.getheader(first)['REGSMTH'] == 2.
    assert last.read_bytes() == last_before
    before_retry = {p: p.read_bytes() for p in (reference, first, last)}
    with pytest.raises(ValueError, match='REGSMTH'):
        Registration(collection(tmp_path), overwrite=True).register_images("SDSS_g'")
    assert all(p.read_bytes() == value for p, value in before_retry.items())


def test_field_mode_selection_is_independent_between_fields(tmp_path, monkeypatch):
    a = image(tmp_path, 'a.fits')
    b = image(tmp_path, 'b.fits', value=200.)
    c = image(tmp_path, 'c.fits', field='field_b', value=300.)
    d = image(tmp_path, 'd.fits', field='field_b', value=400.)
    smooth_inputs = []

    def ordinary(source, reference, **kwargs):
        if source.header['FIELD'] == 'field_b':
            raise registration_module.MaxIterError('field b needs smoothing')
        return source.data.copy(), np.zeros(source.shape, bool)

    def smoothed(source, reference, *masks, **kwargs):
        smooth_inputs.append(source.header['FIELD'])
        return source.data.copy(), np.zeros(source.shape, bool), None
    monkeypatch.setattr(registration_module, 'register', ordinary)
    monkeypatch.setattr(Registration, '_align_with_cmask', staticmethod(smoothed))
    Registration(collection(tmp_path), overwrite=True).register_images(
        "SDSS_g'", fallback_smoothing_sigma=2.)
    assert [fits.getheader(p)['REGSMTH'] for p in (a, b, c, d)] == [0., 0., 2., 2.]
    assert smooth_inputs == ['field_b']


def test_exhausted_smoothing_error_reports_phase_width_and_thresholds(tmp_path, monkeypatch):
    reference = image(tmp_path, 'a.fits')
    target = image(tmp_path, 'b.fits')
    before = {p: p.read_bytes() for p in (reference, target)}

    def fail(*args, **kwargs):
        raise registration_module.MaxIterError('no triangles')
    monkeypatch.setattr(registration_module, 'register', fail)
    monkeypatch.setattr(registration_module, 'find_transform', fail)
    with pytest.raises(RuntimeError, match=r'FIELD=field_a.*phases_attempted=ordinary \+ smoothed.*smoothing_sigma=2.0') as error:
        Registration(collection(tmp_path), overwrite=True).register_images(
            "SDSS_g'", detection_sigmas=[3, 4], fallback_smoothing_sigma=2.)
    assert 'detection_sigmas=[3, 4]' in str(error.value)
    assert 'No results published' in str(error.value)
    assert isinstance(error.value.__cause__, registration_module.MaxIterError)
    assert all(p.read_bytes() == value for p, value in before.items())


@pytest.mark.parametrize('overwrite', [False, True])
def test_real_fallback_separate_output_matches_single_apply_with_wcs_and_masks(tmp_path, monkeypatch, overwrite):
    """Independent pixel-level oracle through staging AND final FITS publication."""
    from astroalign import find_transform as actual_find, apply_transform as actual_apply
    from scipy.ndimage import gaussian_filter
    from skimage.transform import warp

    native = tmp_path / 'native'
    output = tmp_path / 'registered'
    data, _ = synthetic_star_field()
    validity = np.zeros(data.shape, bool)
    validity[120, 120] = True
    reference = image(native, 'a.fits', value=data, mask=validity)
    moved = shift(data, (2.4, -3.1), order=3, mode='constant', cval=100.)
    target = image(native, 'b.fits', value=moved, band="SDSS_r'", mask=validity)
    cmask = np.zeros(data.shape, np.uint8)
    cmask[60, 60] = 1
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [64., 64.]
    wcs.wcs.crval = [163., -5.]
    wcs.wcs.cdelt = [-.0001, .0001]
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    for path in (reference, target):
        with fits.open(path, mode='update') as hdul:
            hdul[0].header.update(wcs.to_header())
            hdul.append(fits.ImageHDU(cmask, name='CMASK'))
    originals = {p: p.read_bytes() for p in (reference, target)}
    src = Registration._read_ccd(target)
    ref = Registration._read_ccd(reference)

    def detection(ccd):
        result = ccd.copy()
        invalid = ccd.mask | cmask.astype(bool)
        weights = gaussian_filter((~invalid).astype(float), 2.)
        result.data = gaussian_filter(np.where(invalid, 0., ccd.data), 2.) / weights
        result.mask = invalid
        return result

    transform, _ = actual_find(detection(src), detection(ref), min_area=5, detection_sigma=3)
    expected, footprint = actual_apply(transform, src, ref, propagate_mask=True)
    expected_cmask = warp(cmask, inverse_map=transform.inverse, output_shape=ref.shape,
                          order=0, mode='constant', cval=0, preserve_range=True).astype(bool)
    calls = []

    def matching(source, destination, **kwargs):
        calls.append(kwargs['detection_sigma'])
        if np.array_equal(source.data, src.data):
            raise registration_module.MaxIterError('force ordinary-mode exhaustion')
        return actual_find(source, destination, **kwargs)
    monkeypatch.setattr(registration_module, 'find_transform', matching)
    registration = Registration(collection(native), overwrite=overwrite)
    registration.register_images(
        "SDSS_g'", detection_sigmas=[3, 4], fallback_smoothing_sigma=2.,
        save_location=output, discard_input_wcs=True)
    assert calls == [3, 4, 3]
    assert all(p.read_bytes() == value for p, value in originals.items())
    np.testing.assert_allclose(fits.getdata(output / target.name), expected, rtol=1e-12, atol=1e-10)
    np.testing.assert_array_equal(fits.getdata(output / target.name, 'MASK').astype(bool),
                                  footprint | ~np.isfinite(expected))
    np.testing.assert_array_equal(fits.getdata(output / target.name, 'CMASK').astype(bool), expected_cmask)
    np.testing.assert_array_equal(fits.getdata(output / reference.name), data)
    np.testing.assert_array_equal(fits.getdata(output / reference.name, 'MASK').astype(bool), validity)
    np.testing.assert_array_equal(fits.getdata(output / reference.name, 'CMASK'), cmask)
    assert has_wcs(output / reference.name) and not has_wcs(output / target.name)
    assert all(fits.getheader(output / p.name)['REGSMTH'] == 2. for p in (reference, target))
    count = len(calls)
    # Reusing native inputs + the same output folder is an overwrite decision,
    # not implicit reuse of destination pixels.
    with pytest.raises(FileExistsError):
        Registration(collection(native)).register_images(
            "SDSS_g'", save_location=output, fallback_smoothing_sigma=2., discard_input_wcs=True)
    assert len(calls) == count
    # A successful processor now selects its returned output collection. To
    # intentionally recompute native pixels, construct from the native inputs.
    Registration(collection(native)).register_images(
        "SDSS_g'", detection_sigmas=[3, 4], fallback_smoothing_sigma=2.,
        save_location=output, overwrite=True, discard_input_wcs=True)
    assert calls == [3, 4, 3, 3, 4, 3]
    assert all(p.read_bytes() == value for p, value in originals.items())
    np.testing.assert_allclose(fits.getdata(output / target.name), expected, rtol=1e-12, atol=1e-10)


def test_failed_separate_publication_can_recover_from_unchanged_native_inputs(tmp_path, monkeypatch):
    native = tmp_path / 'native'
    output = tmp_path / 'registered'
    reference = image(native, 'a.fits')
    target = image(native, 'b.fits', value=200.)
    original = {p: p.read_bytes() for p in (reference, target)}

    def fail(*args, **kwargs):
        raise registration_module.MaxIterError('force smooth mode')

    def smooth(source, reference, *masks, **kwargs):
        return source.data.copy(), np.zeros(source.shape, bool), None
    monkeypatch.setattr(registration_module, 'register', fail)
    monkeypatch.setattr(Registration, '_align_with_cmask', staticmethod(smooth))
    update = registration_module.update_fits_header

    def publish(path, cards, **kwargs):
        if not cards and kwargs.get('output_path') == output / target.name:
            raise OSError('simulated output disk error')
        return update(path, cards, **kwargs)
    monkeypatch.setattr(registration_module, 'update_fits_header', publish)
    with pytest.raises(OSError, match='publication failed'):
        Registration(collection(native)).register_images(
            "SDSS_g'", fallback_smoothing_sigma=2., save_location=output)
    assert (output / reference.name).exists() and not (output / target.name).exists()
    assert all(p.read_bytes() == value for p, value in original.items())
    monkeypatch.setattr(registration_module, 'update_fits_header', update)
    Registration(collection(native), overwrite=True).register_images(
        "SDSS_g'", fallback_smoothing_sigma=2., save_location=output)
    assert all(p.read_bytes() == value for p, value in original.items())
    assert fits.getheader(output / target.name)['REGSMTH'] == 2.
