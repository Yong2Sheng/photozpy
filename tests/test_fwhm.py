"""Scientific FWHM regression checks on synthetic stars, never observing data."""

from astropy.io import fits
from astropy.modeling.models import Moffat2D
from astropy.wcs import WCS
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
import pytest

from photozpy import FWHM, ImageCollectionGroup


def star_field(width=5., beta=3., *, sky=500., gradient=False, outlier=False, count=9):
    rng = np.random.default_rng(1041)
    yy, xx = np.mgrid[:384, :384]
    data = sky + rng.normal(0., 2., xx.shape)
    if gradient:
        data += 0.08 * xx + 0.04 * yy
    positions = [(x, y) for y in (65.3, 180.1, 307.2) for x in (66.2, 180.4, 306.3)][:count]
    for index, (x, y) in enumerate(positions):
        fwhm = 10. if outlier and index == 0 else width
        gamma = fwhm / (2 * np.sqrt(2**(1 / beta) - 1))
        data += Moffat2D(1000. + 100 * index, x, y, gamma, beta)(xx, yy)
    return data, positions


def image(directory, name='master.fits', *, width=5., beta=3., data=None, **cards):
    directory.mkdir(parents=True, exist_ok=True)
    header = fits.Header({
        'IMTYPE': 'Master Light', 'FIELD': 'field_a', 'FILTER': "SDSS_g'",
        'BUNIT': 'ADU', 'NCOMBINE': 10,
    })
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [192., 192.]
    wcs.wcs.crval = [100., -30.]
    wcs.wcs.cdelt = [-0.0001, 0.0001]
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    header.update(wcs.to_header())
    header.update(cards)
    if data is None:
        data, _ = star_field(width, beta)
    validity = np.zeros(data.shape, dtype=np.uint8)
    validity[:, :4] = 1
    # Deliberately mark every pixel as a past cosmic detection. CMASK must not
    # silently become an exclusion array during this measurement.
    cmask = np.ones(data.shape, dtype=np.uint8)
    path = directory / name
    fits.HDUList([
        fits.PrimaryHDU(data, header),
        fits.ImageHDU(validity, name='MASK'),
        fits.ImageHDU(cmask, name='CMASK'),
    ]).writeto(path, checksum=True)
    return path


@pytest.mark.parametrize('width,beta', [(3., 2.), (5., 3.), (8., 4.5)])
def test_recovers_analytic_moffat_width_and_preserves_all_science_arrays(tmp_path, width, beta):
    path = image(tmp_path, width=width, beta=beta)
    with fits.open(path) as hdul:
        arrays = [hdu.data.copy() for hdu in hdul]
        old_wcs = WCS(hdul[0].header).to_header()
    estimator = FWHM(ImageCollectionGroup(tmp_path), overwrite=True)
    working = estimator.find_fwhm()
    assert working is estimator.image_collection
    assert not estimator.failed_paths
    assert len(estimator.results) == 1
    result = estimator.results[0]
    assert result['fwhm'] == pytest.approx(width, rel=0.03)
    assert result['nstars'] >= 5
    with fits.open(path, checksum=True) as hdul:
        assert hdul[0].header['FWHM'] == pytest.approx(width, rel=0.03)
        assert hdul[0].header['FWHMMET'] == 'Moffat'
        assert hdul[0].header['FWHMN'] == result['nstars']
        assert WCS(hdul[0].header).to_header() == old_wcs
        for hdu, original in zip(hdul, arrays):
            np.testing.assert_array_equal(hdu.data, original)
            assert hdu.verify_checksum() == 1
            assert hdu.verify_datasum() == 1


def test_background_offset_and_gradient_do_not_change_width(tmp_path):
    first, _ = star_field(sky=50.)
    second, _ = star_field(sky=5000., gradient=True)
    image(tmp_path, 'a.fits', data=first)
    image(tmp_path, 'b.fits', data=second)
    estimator = FWHM(ImageCollectionGroup(tmp_path))
    originals = {path: path.read_bytes() for path in tmp_path.glob('*.fits')}
    estimator.find_fwhm(write_fwhm=False)
    assert not estimator.failed_paths
    widths = np.asarray(estimator.results['fwhm'])
    assert len(widths) == 2
    assert np.ptp(widths) < 0.1
    assert {path: path.read_bytes() for path in originals} == originals


def test_per_image_widths_are_not_pooled_and_plot_shows_counts(tmp_path):
    a = image(tmp_path, 'a.fits', width=4., FIELD='field_a')
    b = image(tmp_path, 'b.fits', width=8., FIELD='field_b', FILTER="SDSS_i'")
    estimator = FWHM(ImageCollectionGroup(tmp_path), overwrite=True)
    estimator.find_fwhm()
    assert not estimator.failed_paths
    assert fits.getheader(a)['FWHM'] == pytest.approx(4., rel=0.03)
    assert fits.getheader(b)['FWHM'] == pytest.approx(8., rel=0.03)
    figure = estimator.plot_fwhm()
    assert isinstance(figure, Figure)
    axes = figure.axes[0]
    assert [tick.get_text() for tick in axes.get_xticklabels()] == ['field_a', 'field_b']
    assert all(text.get_text().startswith('n=') for text in axes.texts)
    assert 'scatter' in axes.get_title()
    plt.close(figure)


def test_width_outlier_is_rejected_within_image(tmp_path):
    data, positions = star_field(outlier=True)
    image(tmp_path, data=data)
    estimator = FWHM(ImageCollectionGroup(tmp_path))
    estimator.find_fwhm(write_fwhm=False)
    assert not estimator.failed_paths
    assert estimator.results['fwhm'][0] == pytest.approx(5., rel=0.03)
    rejected = estimator.star_measurements[estimator.star_measurements['reason'] == 'fwhm_outlier']
    assert len(rejected) >= 1
    assert np.min(np.hypot(rejected['x'] - positions[0][0], rejected['y'] - positions[0][1])) < 1.


def test_neighbor_saturated_and_elongated_candidates_do_not_define_width(tmp_path):
    data, _ = star_field()
    yy, xx = np.indices(data.shape)
    data += Moffat2D(2000., 120., 120., 5., 3.)(xx, yy)
    data += Moffat2D(1600., 133., 120., 5., 3.)(xx, yy)
    data += Moffat2D(20000., 250., 120., 5., 3.)(xx, yy)
    # A strongly elongated object with a distinctly non-circular profile.
    data += 1800. * (1 + ((xx - 120.) / 12.)**2 + ((yy - 250.) / 3.)**2)**(-3.)
    image(tmp_path, data=data)
    estimator = FWHM(ImageCollectionGroup(tmp_path))
    estimator.find_fwhm(write_fwhm=False, saturation_adu=10000.)
    assert not estimator.failed_paths
    assert estimator.results['fwhm'][0] == pytest.approx(5., rel=0.03)
    stars = estimator.star_measurements
    assert 'neighbor' in stars['reason']
    assert 'saturated' in stars['reason']
    elongated = stars[np.hypot(stars['x'] - 120., stars['y'] - 250.) < 5.]
    assert len(elongated) and not np.any(elongated['accepted'])


def test_insufficient_stars_preserve_existing_header_and_report_failure(tmp_path):
    data, _ = star_field(count=2)
    path = image(tmp_path, data=data, FWHM=99.)
    original = path.read_bytes()
    estimator = FWHM(ImageCollectionGroup(tmp_path), overwrite=True)
    estimator.find_fwhm()
    assert path in estimator.failed_paths
    assert not len(estimator.results)
    assert path.read_bytes() == original
    with pytest.raises(ValueError, match='No successful'):
        estimator.plot_fwhm()


def test_empty_validity_area_is_not_used_as_background_or_sources(tmp_path):
    path = image(tmp_path)
    with fits.open(path, mode='update') as hdul:
        hdul['MASK'].data[:] = 1
    estimator = FWHM(ImageCollectionGroup(tmp_path))
    estimator.find_fwhm(write_fwhm=False)
    assert 'valid image area' in estimator.failed_paths[path]


@pytest.mark.parametrize('options', [
    {'fit_box_size': 40}, {'fit_box_size': 9}, {'min_sources': 1},
    {'background_inner_radius': 10}, {'max_axis_ratio': 1},
    {'detect_threshold': 0}, {'max_sources': 3, 'min_sources': 5},
])
def test_invalid_controls_fail_without_writes(tmp_path, options):
    path = image(tmp_path)
    original = path.read_bytes()
    with pytest.raises(ValueError):
        FWHM(ImageCollectionGroup(tmp_path), overwrite=True).find_fwhm(**options)
    assert path.read_bytes() == original


def test_header_updates_require_permission(tmp_path):
    image(tmp_path)
    with pytest.raises(ValueError, match='overwrite=True'):
        FWHM(ImageCollectionGroup(tmp_path)).find_fwhm()
