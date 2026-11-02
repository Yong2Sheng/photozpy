"""Scientific and workflow regression checks using disposable synthetic CCDs."""

from pathlib import Path

import astropy.units as u
from astropy.io import fits
from astropy.nddata import CCDData, StdDevUncertainty
from ccdproc import ImageFileCollection
import numpy as np
import pytest

from photozpy import Combine, ImageCollectionGroup, Reduction
from photozpy.telescope import SARA_CT


def image(directory, name, kind, value, **cards):
    directory.mkdir(parents=True, exist_ok=True)
    data = np.asarray(value) if np.ndim(value) else np.full((2, 2), value, dtype=float)
    header = fits.Header({'IMTYPE': kind, 'BUNIT': 'ADU', 'OBJECT': 'original'})
    header.update(cards)
    path = directory / name
    fits.PrimaryHDU(data=data, header=header).writeto(path)
    return path


def collection(directory):
    return ImageFileCollection(location=directory, keywords='*', glob_exclude='.*')


def contents(path: Path) -> np.ndarray:
    return fits.getdata(path)


def test_complete_bias_dark_flat_chain_has_analytic_result(tmp_path):
    """Recover source * mean(response), including distinct flat exposures."""
    response = np.array([[1.0, 1.0], [2.0, 4.0]])
    image(tmp_path, 'bias_1.fits', 'Bias', 9)
    image(tmp_path, 'bias_2.fits', 'Bias', 11)
    image(tmp_path, 'dark_1.fits', 'Dark', 30, EXPTIME=10)
    image(tmp_path, 'dark_2.fits', 'Dark', 30, EXPTIME=10)
    image(tmp_path, 'flat_1.fits', 'Flat', 10 + 2 + 100 * response, EXPTIME=1, FILTER="SDSS_g'")
    image(tmp_path, 'flat_2.fits', 'Flat', 10 + 4 + 200 * response, EXPTIME=2, FILTER="SDSS_g'")
    light = image(tmp_path, 'target.fits', 'Light', 10 + 10 + 20 * response, EXPTIME=5, FILTER="SDSS_g'", FIELD='field_a')

    working = Combine(ImageCollectionGroup(tmp_path)).combine_bias_or_dark('Bias', quite=True)
    assert any('Master_Bias.fits' in child.files for child in working)
    working = Reduction(working, overwrite=True).apply_bias_correction()
    np.testing.assert_allclose(contents(tmp_path / 'dark_1.fits'), 20)

    working = Combine(working).combine_bias_or_dark('Dark', quite=True)
    working = Reduction(working, overwrite=True).apply_dark_correction()
    np.testing.assert_allclose(contents(light), 20 * response)

    working = Combine(working).combine_light_or_flat('Flat', quite=True)
    master = tmp_path / "Master_Flat_SDSS_g'.fits"
    np.testing.assert_allclose(contents(master), response / np.median(response), rtol=1e-6)
    working = Reduction(working, overwrite=True).apply_flat_correction()

    np.testing.assert_allclose(contents(light), 20 * np.mean(response), rtol=1e-6)
    header = fits.getheader(light)
    assert header['BIASCORR'] == header['DARKCORR'] == header['FLATCORR'] == 'Yes'
    assert header['FIELD'] == 'field_a'
    assert header['OBJECT'] == 'original'
    assert contents(light).dtype.kind == 'f' and contents(light).dtype.itemsize == 4
    assert len(working[0].files) == 10  # Seven inputs and three masters.


def test_bias_sigma_clipping_retains_original_median_mad_recipe(tmp_path):
    for index, value in enumerate((10, 10, 10, 10, 1000)):
        image(tmp_path, f'bias_{index}.fits', 'Bias', value)
    combiner = Combine(collection(tmp_path))
    output = combiner.combine_bias_or_dark('Bias', quite=True)
    master = combiner.master_paths[0]
    np.testing.assert_array_equal(contents(master), np.full((2, 2), 10))
    assert fits.getheader(master)['NCOMBINE'] == 5
    assert fits.getheader(master)['IMTYPE'] == 'Master Bias'
    assert master.name in output.files
    assert contents(master).dtype.itemsize == 4


def test_single_flat_copy_and_multiple_filters_do_not_depend_on_object(tmp_path):
    data = np.array([[10., 20.], [30., 40.]])
    image(tmp_path, 'flat_g.fits', 'Flat', data, FILTER='SDSS_g', OBJECT='FLAT FRAME')
    image(tmp_path, 'flat_r.fits', 'Flat', 3 * data, FILTER='SDSS_r', OBJECT='ANOTHER LABEL')
    combiner = Combine(collection(tmp_path), SARA_CT)
    combiner.combine_light_or_flat('Flat', quite=True)
    assert len(combiner.master_paths) == 2
    np.testing.assert_array_equal(contents(tmp_path / "Master_Flat_SDSS_g'.fits"), data)
    assert fits.getheader(tmp_path / "Master_Flat_SDSS_r'.fits")['NCOMBINE'] == 1


def test_custom_flat_scale_is_used_instead_of_unbound_auto_variable(tmp_path):
    image(tmp_path, 'flat_1.fits', 'Flat', 2, FILTER='g')
    image(tmp_path, 'flat_2.fits', 'Flat', 4, FILTER='g')
    combiner = Combine(collection(tmp_path))
    combiner.combine_light_or_flat('Flat', scale_function=lambda array: 2.0, quite=True)
    np.testing.assert_allclose(contents(combiner.master_paths[0]), 6.0)


@pytest.mark.parametrize('problem', ['exposure', 'state', 'shape', 'unit'])
def test_invalid_combine_inputs_fail_before_master_creation(tmp_path, problem):
    image(tmp_path, 'dark_1.fits', 'Dark', 30, EXPTIME=10)
    cards = {'EXPTIME': 20 if problem == 'exposure' else 10}
    if problem == 'state':
        cards['BIASCORR'] = 'Yes'
    if problem == 'unit':
        cards['BUNIT'] = 'electron'
    data = np.ones((3, 3)) if problem == 'shape' else 30
    image(tmp_path, 'dark_2.fits', 'Dark', data, **cards)
    with pytest.raises(ValueError):
        Combine(collection(tmp_path)).combine_bias_or_dark('Dark', quite=True)
    assert not (tmp_path / 'Master_Dark.fits').exists()


def test_master_overwrite_is_explicit(tmp_path):
    image(tmp_path, 'bias.fits', 'Bias', 10)
    master = image(tmp_path, 'Master_Bias.fits', 'Master Bias', 20)
    before = master.read_bytes()
    with pytest.raises(FileExistsError):
        Combine(collection(tmp_path)).combine_bias_or_dark('Bias', quite=True)
    assert master.read_bytes() == before


def test_registered_light_grouping_uses_field_not_object(tmp_path):
    image(tmp_path, 'a.fits', 'Light', 5, FILTER='g', FIELD='field_a')
    image(tmp_path, 'b.fits', 'Light', 9, FILTER='g', FIELD='field_b')
    combiner = Combine(collection(tmp_path))
    combiner.combine_light_or_flat('Light', quite=True)
    assert len(combiner.master_paths) == 2
    np.testing.assert_allclose(contents(tmp_path / 'Master_field_a_g.fits'), 5)
    np.testing.assert_allclose(contents(tmp_path / 'Master_field_b_g.fits'), 9)


def test_bias_flag_is_written_only_to_corrected_copy(tmp_path):
    raw, out = tmp_path / 'analysis', tmp_path / 'alternate_output'
    master = image(raw, 'Master_Bias.fits', 'Master Bias', 10)
    target = image(raw, 'light.fits', 'Light', 25)
    before = target.read_bytes()
    reducer = Reduction(collection(raw))
    working = reducer.apply_bias_correction(save_location=out)
    assert target.read_bytes() == before
    assert 'BIASCORR' not in fits.getheader(target)
    np.testing.assert_allclose(contents(out / target.name), 15)
    assert fits.getheader(out / target.name)['BIASCORR'] == 'Yes'
    assert master.is_file()
    assert any(str(out / target.name) in str(name) for name in working.files)


def test_scaled_dark_flag_and_value_only_change_the_output(tmp_path):
    raw, out = tmp_path / 'analysis', tmp_path / 'alternate_output'
    image(raw, 'Master_Dark.fits', 'Master Dark', 20, EXPTIME=10, BIASCORR='Yes')
    target = image(raw, 'flat.fits', 'Flat', 100, EXPTIME=2, BIASCORR='Yes')
    before = target.read_bytes()
    Reduction(collection(raw)).apply_dark_correction(save_location=out)
    assert target.read_bytes() == before
    np.testing.assert_allclose(contents(out / target.name), 96)
    assert fits.getheader(out / target.name)['DARKCORR'] == 'Yes'


def test_flat_uses_mean_normalization_not_median(tmp_path):
    image(tmp_path, 'Master_Flat.fits', 'Master Flat', [[1., 1.], [2., 8.]], FILTER='g')
    target = image(tmp_path, 'light.fits', 'Light', [[10., 10.], [20., 80.]], FILTER='g')
    Reduction(collection(tmp_path), overwrite=True).apply_flat_correction()
    np.testing.assert_allclose(contents(target), 30)
    assert fits.getheader(target)['FLATCORR'] == 'Yes'


@pytest.mark.parametrize('kind,master_type,target_type,flag', [
    ('bias', 'Master Bias', 'Dark', 'BIASCORR'),
    ('dark', 'Master Dark', 'Flat', 'DARKCORR'),
    ('flat', 'Master Flat', 'Light', 'FLATCORR'),
])
def test_repeated_correction_is_rejected_even_with_overwrite(tmp_path, kind, master_type, target_type, flag):
    image(tmp_path, 'master.fits', master_type, 1, EXPTIME=10, FILTER='g', BIASCORR='Yes')
    target = image(tmp_path, 'target.fits', target_type, 10, EXPTIME=10, FILTER='g', **{flag: 'Yes'})
    before = target.read_bytes()
    reducer = Reduction(collection(tmp_path), overwrite=True)
    with pytest.raises(ValueError, match='already Yes'):
        getattr(reducer, f'apply_{kind}_correction')()
    assert target.read_bytes() == before


@pytest.mark.parametrize('problem', ['missing', 'duplicate', 'shape', 'time', 'state'])
def test_bad_master_or_exposure_does_not_partially_correct_batch(tmp_path, problem):
    extra = {'EXPTIME': 10, 'BIASCORR': 'Yes'}
    if problem != 'missing':
        image(tmp_path, 'master.fits', 'Master Dark', 20, **extra)
    if problem == 'duplicate':
        image(tmp_path, 'master_2.fits', 'Master Dark', 20, **extra)
    first = image(tmp_path, 'a.fits', 'Flat', 100, **extra)
    if problem == 'time':
        extra['EXPTIME'] = 0
    if problem == 'state':
        extra['BIASCORR'] = 'No'
    data = np.ones((3, 3)) if problem == 'shape' else 100
    second = image(tmp_path, 'b.fits', 'Light', data, **extra)
    before = {path: path.read_bytes() for path in (first, second)}
    with pytest.raises(ValueError):
        Reduction(collection(tmp_path), overwrite=True).apply_dark_correction()
    assert {path: path.read_bytes() for path in before} == before


def test_missing_flat_in_later_filter_fails_before_any_write(tmp_path):
    image(tmp_path, 'Master_Flat_g.fits', 'Master Flat', 2, FILTER='g')
    a = image(tmp_path, 'a.fits', 'Light', 20, FILTER='g')
    b = image(tmp_path, 'b.fits', 'Light', 20, FILTER='r')
    before = {path: path.read_bytes() for path in (a, b)}
    with pytest.raises(ValueError):
        Reduction(collection(tmp_path), overwrite=True).apply_flat_correction()
    assert {path: path.read_bytes() for path in before} == before


def test_groups_use_their_own_masters(tmp_path):
    one, two = tmp_path / 'night_1', tmp_path / 'night_2'
    image(one, 'Master_Bias.fits', 'Master Bias', 10)
    image(two, 'Master_Bias.fits', 'Master Bias', 30)
    a = image(one, 'light.fits', 'Light', 100)
    b = image(two, 'light.fits', 'Light', 100)
    result = Reduction(ImageCollectionGroup((one, two)), overwrite=True).apply_bias_correction()
    assert isinstance(result, ImageCollectionGroup)
    np.testing.assert_allclose(contents(a), 90)
    np.testing.assert_allclose(contents(b), 70)


def test_reduction_output_collisions_are_detected_before_writing(tmp_path):
    one, two = tmp_path / 'night_1', tmp_path / 'night_2'
    for directory in (one, two):
        image(directory, 'Master_Bias.fits', 'Master Bias', 10)
        image(directory, 'light.fits', 'Light', 100)
    output = tmp_path / 'output'
    with pytest.raises(ValueError, match='same output'):
        Reduction(ImageCollectionGroup((one, two))).apply_bias_correction(save_location=output)
    assert not output.exists()


def test_masks_uncertainties_removed_as_in_original_reduction(tmp_path):
    image(tmp_path, 'Master_Bias.fits', 'Master Bias', 10)
    target = tmp_path / 'light.fits'
    ccd = CCDData(np.full((2, 2), 100.), unit=u.adu, meta={'IMTYPE': 'Light'})
    ccd.mask = np.zeros((2, 2), dtype=bool)
    ccd.uncertainty = StdDevUncertainty(np.ones((2, 2)))
    ccd.write(target)
    Reduction(collection(tmp_path), overwrite=True).apply_bias_correction()
    result = CCDData.read(target)
    assert result.mask is None
    assert result.uncertainty is None
    assert len(fits.getheader(target)) > 0


def test_unsigned_input_does_not_retain_stale_bzero_after_subtraction(tmp_path):
    image(tmp_path, 'Master_Bias.fits', 'Master Bias', 10)
    target = image(tmp_path, 'light.fits', 'Light', np.full((2, 2), 100, dtype=np.uint16))
    Reduction(collection(tmp_path), overwrite=True).apply_bias_correction()
    np.testing.assert_allclose(contents(target), 90)
    assert 'BZERO' not in fits.getheader(target)


def test_no_target_and_misleading_type_checks_are_explicit(tmp_path):
    image(tmp_path, 'bias.fits', 'Bias', 10)
    assert Combine.check_imtype(collection(tmp_path), 'Bias')
    assert not Combine.check_imtype(collection(tmp_path), 'Dark')
    with pytest.raises(ValueError, match='No images'):
        Reduction(collection(tmp_path)).apply_bias_correction()


def test_combine_keeps_night_boundaries_and_preflights_shared_output(tmp_path):
    one, two = tmp_path / 'night_1', tmp_path / 'night_2'
    a = image(one, 'bias.fits', 'Bias', 10)
    b = image(two, 'bias.fits', 'Bias', 30)
    before = {path: path.read_bytes() for path in (a, b)}
    group = ImageCollectionGroup((one, two))
    output = tmp_path / 'masters'
    with pytest.raises(ValueError, match='same output'):
        Combine(group).combine_bias_or_dark('Bias', save_location=output, quite=True)
    assert not output.exists()

    result = Combine(group).combine_bias_or_dark('Bias', quite=True)
    assert len(result) == 2
    np.testing.assert_allclose(contents(one / 'Master_Bias.fits'), 10)
    np.testing.assert_allclose(contents(two / 'Master_Bias.fits'), 30)
    assert {path: path.read_bytes() for path in before} == before


def test_reduced_output_has_valid_checksums(tmp_path):
    image(tmp_path, 'Master_Bias.fits', 'Master Bias', 10)
    target = tmp_path / 'light.fits'
    header = fits.Header({'IMTYPE': 'Light', 'BUNIT': 'ADU'})
    fits.PrimaryHDU(np.full((2, 2), 100.), header).writeto(target, checksum=True)
    Reduction(collection(tmp_path), overwrite=True).apply_bias_correction()
    with fits.open(target, checksum=True) as hdul:
        assert hdul[0].verify_checksum() == 1
        assert hdul[0].verify_datasum() == 1
        np.testing.assert_allclose(hdul[0].data, 90)


def test_twilight_flat_clips_transient_star_after_sky_normalization(tmp_path):
    """A star hidden by changing raw sky levels is rejected after scaling."""
    originals = {}
    for index, sky in enumerate((1000., 2000., 3000., 4000., 5000.)):
        data = np.full((5, 5), sky)
        if index == 2:
            data[2, 2] += 1000.
        path = image(
            tmp_path, f'flat_{index}.fits', 'Flat', data,
            FILTER='g', BIASCORR='Yes', DARKCORR='Yes',
        )
        originals[path] = path.read_bytes()

    combiner = Combine(collection(tmp_path))
    combiner.combine_light_or_flat('Flat', quite=True)
    master = combiner.master_paths[0]

    np.testing.assert_allclose(contents(master), 1.0, atol=1e-6)
    assert {path: path.read_bytes() for path in originals} == originals
    header = fits.getheader(master)
    assert header['BIASCORR'] == header['DARKCORR'] == 'Yes'
    assert header['FILTER'] == 'g'
    assert header['NCOMBINE'] == 5


def test_flat_rejection_does_not_erase_persistent_spatial_structure(tmp_path):
    """A fixed feature is not a temporal outlier; normalization must retain it."""
    for index, sky in enumerate((1000., 2000., 3000., 4000., 5000.)):
        data = np.full((5, 5), sky)
        data[2, 2] *= 1.2
        image(tmp_path, f'flat_{index}.fits', 'Flat', data, FILTER='g')

    combiner = Combine(collection(tmp_path))
    combiner.combine_light_or_flat('Flat', quite=True)
    expected = np.ones((5, 5))
    expected[2, 2] = 1.2
    np.testing.assert_allclose(contents(combiner.master_paths[0]), expected, atol=1e-6)


def test_explicitly_disabled_flat_scaling_stays_unscaled(tmp_path):
    image(tmp_path, 'flat_1.fits', 'Flat', 2, FILTER='g')
    image(tmp_path, 'flat_2.fits', 'Flat', 4, FILTER='g')
    combiner = Combine(collection(tmp_path))
    combiner.combine_light_or_flat('Flat', scale_function=None, quite=True)
    np.testing.assert_allclose(contents(combiner.master_paths[0]), 3)


@pytest.mark.parametrize('factor', [0, -1, np.nan, np.inf])
def test_invalid_flat_normalization_cannot_write_a_master(tmp_path, factor):
    image(tmp_path, 'flat_1.fits', 'Flat', 2, FILTER='g')
    image(tmp_path, 'flat_2.fits', 'Flat', 4, FILTER='g')
    with pytest.raises(ValueError, match='positive finite scalar'):
        Combine(collection(tmp_path)).combine_light_or_flat(
            'Flat', scale_function=lambda array: factor, quite=True,
        )
    assert not (tmp_path / 'Master_Flat_g.fits').exists()
