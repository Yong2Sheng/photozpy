"""Synthetic FITS checks; no observatory data are modified by these tests."""

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from ccdproc import ImageFileCollection

from photozpy import HeaderCorrection
from photozpy.collections import ImageCollectionGroup
from photozpy.sources import ObservationField
from photozpy.telescope import SARA_CT


@pytest.fixture
def field():
    return ObservationField(
        field_name='science_field', telescope=SARA_CT,
        nominal_pointing=SkyCoord(150, -30, unit='deg'), file_pattern='target_',
    )


def write_image(directory, name, **cards):
    directory.mkdir(parents=True, exist_ok=True)
    header = fits.Header({
        'DATE-OBS': '2024-02-16T00:00:00', 'TIMESYS': 'UTC',
        'OBJECT': 'original_object',
    })
    header.update(cards)
    path = directory / name
    fits.HDUList([
        fits.PrimaryHDU(np.arange(6, dtype=np.int16).reshape(2, 3), header),
        fits.ImageHDU(np.full((2, 3), 7.0), name='EXTRA'),
    ]).writeto(path)
    return path


def collection(directory):
    return ImageFileCollection(location=directory, glob_exclude='._*', keywords='*')


def test_field_ccd_time_and_hdus_with_separate_output(tmp_path, field):
    raw = tmp_path / 'raw'
    path = write_image(raw, 'target_SDSS_g_001.fits', RA='10:00:01', DEC='-30:00:01')
    before = path.read_bytes()
    output = tmp_path / 'processed'
    processor = HeaderCorrection(
        collection(raw), observation_fields=field, save_location=output, telescope_name='SARA-CT',
    )
    result = processor.correct_headers_by_filename()
    assert result is processor.image_collection
    assert Path(result.location) == output
    assert path.read_bytes() == before
    with fits.open(output / path.name) as hdul:
        header = hdul[0].header
        assert header['FIELD'] == 'science_field'
        assert header['OBJECT'] == 'original_object'
        assert header['RA'] == '10:00:01'
        assert header['DEC'] == '-30:00:01'
        assert header['GAIN'] == 1.0
        assert header['RDNOISE'] == 7.1
        assert header['TELESCOP'] == 'SARA-CT'
        assert header['TIMESYS'] == 'UTC'
        assert header['MJD-OBS'] == pytest.approx(60356.0, abs=1e-9)
        assert len(hdul) == 2
        np.testing.assert_array_equal(hdul[0].data, np.arange(6).reshape(2, 3))
        np.testing.assert_array_equal(hdul[1].data, np.full((2, 3), 7.0))


def test_standard_field_needs_no_placeholder_source(tmp_path, field):
    standard = replace(field, field_name='standard_field', file_pattern='DLS_')
    path = write_image(tmp_path / 'raw', 'DLS_SDSS_g_001.fits')
    processor = HeaderCorrection(
        collection(path.parent), observation_fields=(field, standard), overwrite=True,
    )
    processor.correct_headers_by_filename()
    header = fits.getheader(path)
    assert header['FIELD'] == 'standard_field'
    assert header['RA'] == '10:00:00.000'
    assert header['DEC'] == '-30:00:00.000'


def test_each_file_selects_ccd_by_its_own_observation_date(tmp_path):
    directory = tmp_path / 'raw'
    old = write_image(directory, 'Bias_old.fits', **{'DATE-OBS': '2018-01-02T00:00:00'})
    new = write_image(directory, 'Bias_new.fits')
    HeaderCorrection(collection(directory), telescope=SARA_CT, overwrite=True).correct_headers_by_filename()
    assert fits.getheader(old)['GAIN'] == 2.0
    assert fits.getheader(old)['RDNOISE'] == pytest.approx(19.4)
    assert fits.getheader(new)['GAIN'] == 1.0
    assert fits.getheader(new)['RDNOISE'] == 7.1


@pytest.mark.parametrize('problem', ['unmatched', 'missing_time', 'invalid_date', 'outside_ccd', 'incomplete_pointing', 'wrong_unit'])
def test_bad_input_fails_before_any_batch_write(tmp_path, field, problem):
    directory = tmp_path / 'raw'
    first = write_image(directory, 'target_001.fits')
    cards = {}
    second_name = 'target_002.fits'
    if problem == 'unmatched':
        second_name = 'unknown_002.fits'
    elif problem == 'invalid_date':
        cards['DATE-OBS'] = 'not-a-date'
    elif problem == 'outside_ccd':
        cards['DATE-OBS'] = '2010-01-01T00:00:00'
    elif problem == 'incomplete_pointing':
        cards['RA'] = '10:00:00'
    elif problem == 'wrong_unit':
        cards['BUNIT'] = 'electron'
    second = write_image(directory, second_name, **cards)
    if problem == 'missing_time':
        with fits.open(second, mode='update') as hdul:
            del hdul[0].header['TIMESYS']
    before = {p: p.read_bytes() for p in (first, second)}
    with pytest.raises(ValueError):
        HeaderCorrection(collection(directory), observation_fields=field, overwrite=True).correct_headers_by_filename()
    assert {p: p.read_bytes() for p in (first, second)} == before


def test_missing_timesys_requires_explicit_fallback(tmp_path, field):
    path = write_image(tmp_path / 'raw', 'target_001.fits')
    with fits.open(path, mode='update') as hdul:
        del hdul[0].header['TIMESYS']
    HeaderCorrection(
        collection(path.parent), observation_fields=field, time_scale='utc', overwrite=True,
    ).correct_headers_by_filename()
    assert fits.getheader(path)['TIMESYS'] == 'UTC'


def test_explicit_time_scale_cannot_silently_override_timesys(tmp_path, field):
    path = write_image(tmp_path / 'raw', 'target_001.fits')
    before = path.read_bytes()
    with pytest.raises(ValueError, match='conflicts'):
        HeaderCorrection(
            collection(path.parent), observation_fields=field, time_scale='ut1', overwrite=True,
        ).correct_headers_by_filename()
    assert path.read_bytes() == before


def test_ambiguous_field_rules_fail_without_writing(tmp_path, field):
    path = write_image(tmp_path / 'raw', 'target_001.fits')
    before = path.read_bytes()
    other = replace(field, field_name='other_field')
    with pytest.raises(ValueError, match='ambiguous'):
        HeaderCorrection(
            collection(path.parent), observation_fields=(field, other), overwrite=True,
        ).correct_headers_by_filename()
    assert path.read_bytes() == before


def test_skipped_files_are_not_claimed_as_written_outputs(tmp_path, field):
    directory = tmp_path / 'raw'
    good = write_image(directory, 'target_001.fits')
    unknown = write_image(directory, 'unknown_001.fits')
    output = tmp_path / 'processed'
    processor = HeaderCorrection(collection(directory), observation_fields=field, save_location=output)
    selected = processor.correct_headers_by_filename(unmatched='skip')
    assert selected.files == [good.name]
    assert processor.skipped_files == (unknown,)
    assert not (output / unknown.name).exists()


def test_existing_outputs_require_explicit_overwrite(tmp_path, field):
    path = write_image(tmp_path / 'raw', 'target_001.fits')
    before = path.read_bytes()
    with pytest.raises(FileExistsError):
        HeaderCorrection(collection(path.parent), observation_fields=field).correct_headers_by_filename()
    assert path.read_bytes() == before


@pytest.mark.parametrize('image_type', ['Flat', 'Light'])
def test_filters_need_not_have_both_flat_and_light(tmp_path, image_type):
    path = write_image(
        tmp_path / 'raw', f'{image_type}_SDSS_g_001.fits',
        IMTYPE=image_type, FILTER='SDSS_g',
    )
    HeaderCorrection(
        collection(path.parent), telescope=SARA_CT, overwrite=True,
    ).correct_filter_headers()
    assert fits.getheader(path)['FILTER'] == "SDSS_g'"


def test_filter_output_contains_calibration_files_too(tmp_path):
    directory = tmp_path / 'raw'
    bias = write_image(directory, 'Bias_001.fits', IMTYPE='Bias')
    light = write_image(directory, 'target_001.fits', IMTYPE='Light', FILTER='SDSS r')
    output = tmp_path / 'processed'
    result = HeaderCorrection(
        collection(directory), telescope=SARA_CT, save_location=output,
    ).correct_filter_headers()
    assert set(result.files) == {bias.name, light.name}
    assert (output / bias.name).is_file()
    assert fits.getheader(output / light.name)['FILTER'] == "SDSS_r'"
    assert fits.getheader(light)['FILTER'] == 'SDSS r'


def test_unknown_filter_aborts_without_partial_changes(tmp_path):
    directory = tmp_path / 'raw'
    first = write_image(directory, 'a.fits', IMTYPE='Light', FILTER='SDSS_g')
    second = write_image(directory, 'b.fits', IMTYPE='Light', FILTER='unknown')
    before = {p: p.read_bytes() for p in (first, second)}
    with pytest.raises(ValueError):
        HeaderCorrection(collection(directory), telescope=SARA_CT, overwrite=True).correct_filter_headers()
    assert {p: p.read_bytes() for p in (first, second)} == before


def test_filename_filter_lookup_and_custom_alias_mapping(tmp_path):
    directory = tmp_path / 'raw'
    path = write_image(directory, 'target_SDSS_i_001.fits', IMTYPE='Light')
    processor = HeaderCorrection(collection(directory), telescope=SARA_CT, overwrite=True)
    processor.correct_filter_headers(from_file_name=True)
    assert fits.getheader(path)['FILTER'] == "SDSS_i'"
    processor.correct_filter_headers(filter_mapping={"SDSS_i'": 'SDSS z'})
    assert fits.getheader(path)['FILTER'] == "SDSS_z'"


def test_group_output_and_filename_collision_preflight(tmp_path):
    one, two = tmp_path / 'one', tmp_path / 'two'
    a = write_image(one, 'Bias_001.fits')
    b = write_image(two, 'Bias_001.fits')
    group = ImageCollectionGroup((one, two))
    output = tmp_path / 'processed'
    with pytest.raises(ValueError, match='Duplicate output'):
        HeaderCorrection(group, telescope=SARA_CT, save_location=output).correct_headers_by_filename()
    assert not output.exists()
    before = (a.read_bytes(), b.read_bytes())
    result = HeaderCorrection(group, telescope=SARA_CT, overwrite=True).correct_headers_by_filename()
    assert isinstance(result, ImageCollectionGroup)
    assert len(result) == 2
    assert (a.read_bytes(), b.read_bytes()) != before
    assert all(child.summary['imtype'][0] == 'Bias' for child in result)


def test_refactored_aliases_preserve_legacy_sdss_mapping():
    for band in 'ugriz':
        assert SARA_CT.resolve_filter_alias(f'SDSS_{band}') == f"SDSS_{band}'"
