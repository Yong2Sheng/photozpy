"""Numerical regressions for the original count/error equations and new wiring."""

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from ccdproc import ImageFileCollection
import numpy as np
import pytest

from photozpy.photometry.photometry import Photometry
from photozpy.photometry.regions import PhotometryRegions, PhotometryRegionCollection
from photozpy.sources import ObservationField, SourceInfo, Sources, FilterMagnitude
from photozpy.telescope import SARA_CT


def setup_image(tmp_path, flux=1000.0, *, mask_source=False, cmask=False):
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [49, 49]
    wcs.wcs.crval = [163, -5]
    wcs.wcs.cdelt = [-0.0001, 0.0001]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    data = np.full((96, 96), 60.0)
    data[30, 30] += flux
    header = wcs.to_header()
    header.update(
        IMTYPE="Master Light",
        FIELD="DLS",
        FILTER="SDSS_g'",
        FWHM=3.0,
        BUNIT="adu",
        GAIN=2.0,
        RDNOISE=5.0,
        NCOMBINE=4,
    )
    invalid = np.zeros(data.shape, dtype=np.uint8)
    if mask_source:
        invalid[30, 30] = 1
    path = tmp_path / "Master_DLS_g.fits"
    hdus = [fits.PrimaryHDU(data, header), fits.ImageHDU(invalid, name="MASK")]
    if cmask:
        hdus.append(fits.ImageHDU(np.ones(data.shape, dtype=np.uint8), name="CMASK"))
    fits.HDUList(hdus).writeto(path)
    field = ObservationField(
        field_name="DLS",
        telescope=SARA_CT,
        nominal_pointing=SkyCoord(163, -5, unit="deg"),
    )
    sources = Sources(
        SourceInfo(
            source_name="standard",
            source_role="standard",
            source_coordinate=wcs.pixel_to_world(30, 30),
            observation_field=field,
        )
    )
    reference = FilterMagnitude(
        telescope=SARA_CT,
        filter_names=["SDSS g"],
        magnitudes=[15],
        errors=[0.02],
        significance=None,
        magnitude_system="AB",
    )
    sources.set_magnitudes(kind="reference", magnitudes=reference)
    regions = PhotometryRegions.generate(
        path,
        sources,
        root=tmp_path,
        aperture_scaling=1,
        inner_annulus_scaling=2,
        outer_annulus_scaling=3,
    )
    return path, sources, reference, regions


def review(regions):
    # Background deliberately offset from the source in clean constant sky.
    regions.save_review(
        "image\ncircle(31,31,3) # tag={id:standard} tag={aperture}\n"
        "annulus(66,66,6,9) # tag={id:standard} tag={background}\n"
    )


def test_original_error_formula_numeric_values_and_average_stack_scaling():
    arguments = dict(
        gain=2 * u.electron / u.adu,
        readout_noise=5 * u.electron,
        n_ap=10 * u.pix,
        n_ann=100 * u.pix,
        mean_bkg=5 * u.adu / u.pix,
        aperture_counts=250 * u.adu,
        n_com=4,
    )
    error, background_error = Photometry.estimate_source_region_error(**arguments)
    assert error.value[0] ** 2 == pytest.approx(49.0625)
    assert background_error.value[0] ** 2 == pytest.approx(24.0625)
    arguments["n_com"] = 1
    single, background_single = Photometry.estimate_source_region_error(**arguments)
    np.testing.assert_allclose(single.value, error.value * 2)
    np.testing.assert_allclose(background_single.value, background_error.value * 2)


def test_reviewed_pixel_regions_recover_known_net_counts_and_preserve_other_slots(
    tmp_path,
):
    path, sources, reference, regions = setup_image(tmp_path, cmask=True)
    review(regions)
    originals = [
        path.read_bytes(),
        regions.aperture_path.read_bytes(),
        regions.background_path.read_bytes(),
    ]
    photometry = Photometry(ImageFileCollection(tmp_path, filenames=[path.name]))
    table = photometry.run_photometry(sources, regions=regions, verbose=False)
    assert table["src"][0].to_value(u.adu) == pytest.approx(1000)
    assert table["bkg"][0].to_value(u.adu) == pytest.approx(60 * np.pi * 3**2)
    assert table["mean_bkg"][0].value == 60
    assert table["background_x"][0] == 65
    assert table["aperture_radius"][0] == 3
    assert table["mag_inst"][0].value == pytest.approx(-7.5)
    assert table["status"][0] == "measured"
    assert not table["is_upper_limit"][0]
    instrumental = sources.get_magnitudes(kind="instrumental")
    assert instrumental.magnitudes.value[1] == pytest.approx(-7.5)
    assert sources.get_magnitudes(kind="reference") is reference
    assert sources.get_magnitudes(kind="calibrated") is None
    # Compare pipeline errors to direct calls of the user's unchanged functions.
    expected, expected_ul = Photometry.estimate_source_region_error(
        gain=2 * u.electron / u.adu,
        readout_noise=5 * u.electron,
        n_ap=table["n_ap"][0],
        n_ann=table["n_ann"][0],
        mean_bkg=table["mean_bkg"][0],
        aperture_counts=table["src+bkg"][0],
        n_com=4,
    )
    assert table["error_counts"][0].value == expected.value[0]
    assert table["error_counts_ul"][0].value == expected_ul.value[0]
    assert [
        path.read_bytes(),
        regions.aperture_path.read_bytes(),
        regions.background_path.read_bytes(),
    ] == originals
    table["src"][0] = 0 * u.adu
    assert photometry.results["src"][0].value == pytest.approx(1000)


@pytest.mark.parametrize("flux", [1.0, -2.0])
def test_nondetection_retains_original_limit_and_storage_flag(tmp_path, flux):
    path, sources, _, regions = setup_image(tmp_path, flux=flux)
    review(regions)
    table = Photometry(ImageFileCollection(tmp_path)).run_photometry(
        sources, regions=regions, verbose=False
    )
    assert table["status"][0] == "upper_limit"
    assert np.isnan(table["mag_inst_error"][0].value)
    assert table["mag_inst"][0].value == pytest.approx(
        -2.5 * np.log10(3 * table["error_counts_ul"][0].value)
    )
    assert sources.get_magnitudes(kind="instrumental").errors.value[1] == -99


def test_source_mask_is_not_silently_area_corrected(tmp_path):
    path, sources, reference, regions = setup_image(tmp_path, mask_source=True)
    review(regions)
    table = Photometry(ImageFileCollection(tmp_path)).run_photometry(
        sources, regions=regions, verbose=False
    )
    assert table["status"][0] == "invalid_aperture"
    assert np.isnan(table["src"][0].value)
    assert sources.get_magnitudes(kind="instrumental") is None
    assert sources.get_magnitudes(kind="reference") is reference


def test_review_required_rejection_and_duplicate_measurements(tmp_path):
    path, sources, _, regions = setup_image(tmp_path)
    photometry = Photometry(ImageFileCollection(tmp_path))
    with pytest.raises(ValueError, match="Complete region review"):
        photometry.run_photometry(sources, regions=regions, verbose=False)
    review(regions)
    with pytest.raises(ValueError, match="Duplicate source/filter"):
        photometry.run_photometry(
            sources,
            regions=PhotometryRegionCollection([regions, regions]),
            verbose=False,
        )
    assert sources.get_magnitudes(kind="instrumental") is None
    regions.save_review(regions.initial_regions())
    table = photometry.run_photometry(sources, regions=regions, verbose=False)
    assert table["status"][0] == "rejected"
    assert sources.get_magnitudes(kind="instrumental") is None


def test_invalid_detector_metadata_is_rejected_before_attaching_results(tmp_path):
    path, sources, reference, regions = setup_image(tmp_path)
    review(regions)
    with fits.open(path, mode="update") as hdul:
        hdul[0].header["GAIN"] = 0
    with pytest.raises(ValueError, match="GAIN"):
        Photometry(ImageFileCollection(tmp_path)).run_photometry(
            sources, regions=regions, verbose=False
        )
    assert sources.get_magnitudes(kind="reference") is reference
    assert sources.get_magnitudes(kind="instrumental") is None


def test_masked_background_pixel_changes_effective_area_not_mean(tmp_path):
    path, sources, _, regions = setup_image(tmp_path)
    review(regions)
    measurer = Photometry(ImageFileCollection(tmp_path))
    before = measurer.run_photometry(sources, regions=regions, verbose=False)
    with fits.open(path, mode="update") as hdul:
        hdul[0].data[65, 72] = 1e8
        hdul["MASK"].data[65, 72] = 1
    after = measurer.run_photometry(sources, regions=regions, verbose=False)
    assert after["mean_bkg"][0].value == 60
    assert after["src"][0].value == pytest.approx(1000)
    assert after["n_ann"][0].value == pytest.approx(before["n_ann"][0].value - 1)
    assert after["error_counts"][0].value > before["error_counts"][0].value


def test_multiple_stars_use_each_ids_actual_radius_not_file_order(tmp_path):
    path, original, _, _ = setup_image(tmp_path)
    with fits.open(path, mode="update") as hdul:
        hdul[0].data[30, 65] += 800
        wcs = WCS(hdul[0].header)
    first = original.source_info
    second = SourceInfo(
        source_name="second",
        source_role="standard",
        observation_field=first.observation_field,
        source_coordinate=wcs.pixel_to_world(65, 30),
    )
    sources = Sources(first, second)
    regions = PhotometryRegions.generate(path, sources, root=tmp_path, overwrite=True)
    regions.save_review(
        "image\ncircle(66,31,5) # tag={id:second} tag={aperture}\n"
        "circle(31,31,3) # tag={id:standard} tag={aperture}\n"
        "annulus(66,66,6,9) # tag={id:standard} tag={background}\n"
        "annulus(31,66,8,12) # tag={id:second} tag={background}\n"
    )
    table = Photometry(ImageFileCollection(tmp_path)).run_photometry(
        sources, regions=regions, verbose=False
    )
    assert list(table["source_id"]) == ["standard", "second"]
    np.testing.assert_allclose(table["src"].value, [1000, 800], atol=1e-8)
    np.testing.assert_allclose(table["aperture_radius"], [3, 5])
    assert table["n_ap"][1].value == pytest.approx(25 * np.pi)


def test_unified_photometry_handles_both_roles_in_one_frame_and_updates_originals(
    tmp_path, monkeypatch
):
    path, standards, reference, standard_regions = setup_image(tmp_path)
    with fits.open(path, mode="update") as hdul:
        hdul[0].data[30, 65] += 800
        wcs = WCS(hdul[0].header)
    science = Sources(
        SourceInfo(
            source_name="science-target",
            source_role="science",
            observation_field=standards.source_info.observation_field,
            source_coordinate=wcs.pixel_to_world(65, 30),
        )
    )
    science_regions = PhotometryRegions.generate(path, science, root=tmp_path)
    review(standard_regions)
    science_regions.save_review(
        "image\ncircle(66,31,3) # tag={id:science-target} tag={aperture}\n"
        "annulus(31,66,6,9) # tag={id:science-target} tag={background}\n"
    )
    measurer = Photometry(ImageFileCollection(tmp_path))
    original_open = fits.open
    reads = []

    def counting_open(filename, *args, **kwargs):
        if str(filename) == str(path):
            reads.append(filename)
        return original_open(filename, *args, **kwargs)

    monkeypatch.setattr(fits, "open", counting_open)
    table = measurer.run_photometry(
        sources=(standards, science),
        regions=(standard_regions, science_regions),
        verbose=False,
    )
    assert len(reads) == 1
    assert list(table["source_role"]) == ["standard", "science"]
    assert len(set(table["image_path"])) == 1
    np.testing.assert_allclose(table["src"].value, [1000, 800], atol=1e-8)
    assert standards.get_magnitudes(kind="reference") is reference
    assert standards.get_magnitudes(kind="instrumental") is not None
    assert science.get_magnitudes(kind="instrumental") is not None
    assert all(
        group.get_magnitudes(kind="calibrated") is None
        for group in (standards, science)
    )
    # A single combined Sources object is equally valid, preserving shared slots.
    combined = Sources.combine(standards, science)
    again = measurer.run_photometry(
        combined,
        regions=PhotometryRegionCollection([standard_regions, science_regions]),
        verbose=False,
    )
    np.testing.assert_allclose(again["mag_inst"].value, table["mag_inst"].value)

    # Carry real measured ADU counts through the new zero-point/calibration entry points.
    zero_points = standards.calculate_zero_points()
    assert zero_points.set_index("filter").loc["SDSS_g'", "zero_point"] == pytest.approx(22.5)
    calibrated_table = science.calibrate_targets(standard_sources=standards)
    expected_magnitude = 15 - 2.5 * np.log10(800 / 1000)
    calibrated = science.get_magnitudes(kind="calibrated")
    assert calibrated.magnitudes[1].value == pytest.approx(expected_magnitude)
    np.testing.assert_array_equal(
        calibrated.errors.value, science.get_magnitudes(kind="instrumental").errors.value
    )
    assert calibrated_table.set_index("filter").loc[
        "SDSS_g'", "calibrated_magnitude"
    ] == pytest.approx(expected_magnitude)
    assert standards.get_magnitudes(kind="reference") is reference
    assert standards.get_magnitudes(kind="calibrated") is None
    assert len(reads) == 2  # Calibration reads stored measurements, not FITS images.


def test_conflicting_source_owners_fail_without_replacing_measurements(tmp_path):
    path, sources, reference, regions = setup_image(tmp_path)
    review(regions)
    with pytest.raises(ValueError, match="Conflicting independent"):
        Photometry(ImageFileCollection(tmp_path)).run_photometry(
            (sources, Sources(sources.source_info)),
            regions=regions,
            verbose=False,
        )
    assert sources.get_magnitudes(kind="reference") is reference
    assert sources.get_magnitudes(kind="instrumental") is None


@pytest.mark.parametrize('kwargs', [
    {'mask_neighbors': 'yes'}, {'detection_sigma': 0}, {'min_area': True},
    {'mask_growth_fwhm': float('nan')}, {'min_pixels': 0},
    {'min_area_fraction': 0}, {'min_area_fraction': 1.1},
    {'stability_growth_fwhms': [1, -2]}, {'stability_growth_fwhms': '1,2'},
    {'stability_tolerance': float('inf')},
])
def test_local_background_options_validate(kwargs):
    from photozpy.photometry import LocalBackgroundOptions
    with pytest.raises(ValueError):
        LocalBackgroundOptions(**kwargs)


def test_masked_local_sky_preserves_aperture_sum_and_improves_injected_sky(tmp_path):
    from photozpy.photometry import LocalBackgroundOptions
    path, sources, _, regions = setup_image(tmp_path)
    rng = np.random.default_rng(834)
    with fits.open(path, mode='update') as hdul:
        yy, xx = np.indices(hdul[0].data.shape)
        hdul[0].data += rng.normal(0, 1.5, xx.shape)
        hdul[0].data += 80 * np.exp(-((xx - 82)**2 + (yy - 65)**2) / (2 * 2.6**2))
    regions.save_review(
        'image\ncircle(31,31,3) # tag={id:standard} tag={aperture}\n'
        'annulus(66,66,12,20) # tag={id:standard} tag={background}\n'
    )
    original = path.read_bytes()
    saved = regions.aperture_path.read_bytes(), regions.background_path.read_bytes()
    phot = Photometry(ImageFileCollection(tmp_path))
    plain = phot.run_photometry(sources, regions=regions, verbose=False)
    options = LocalBackgroundOptions()
    masked = phot.run_photometry(sources, regions=regions, background_options=options, verbose=False)
    assert masked['src+bkg'][0] == plain['src+bkg'][0]
    assert abs(masked['mean_bkg'][0].value - 60) < abs(plain['mean_bkg'][0].value - 60)
    assert .5 < masked['background_area_fraction'][0] < 1
    assert masked['background_stability'][0] == 'stable'
    assert masked.meta['local_background_options']['mask_growth_fwhm'] == 1.5
    assert len(__import__('json').loads(masked['background_trials'][0])) == 3
    assert path.read_bytes() == original
    assert saved == (regions.aperture_path.read_bytes(), regions.background_path.read_bytes())


def test_local_sky_area_rejection_and_area_only_mode(tmp_path):
    from photozpy.photometry import LocalBackgroundOptions
    path, sources, _, regions = setup_image(tmp_path)
    review(regions)
    phot = Photometry(ImageFileCollection(tmp_path))
    rejected = phot.run_photometry(
        sources, regions=regions, verbose=False,
        background_options=LocalBackgroundOptions(mask_neighbors=False, min_pixels=10000),
    )
    assert rejected['status'][0] == 'invalid_background'
    assert np.isnan(rejected['mag_inst'][0])
    assert 'effective background area' in rejected['reason'][0]
    with fits.open(path, mode='update') as hdul:
        del hdul[0].header['FWHM']
    accepted = phot.run_photometry(
        sources, regions=regions, verbose=False,
        background_options=LocalBackgroundOptions(mask_neighbors=False),
    )
    assert accepted['status'][0] == 'measured'
    assert accepted['background_stability'][0] == 'not_checked'
    with pytest.raises(ValueError, match='FWHM'):
        phot.run_photometry(sources, regions=regions, background_options=LocalBackgroundOptions())


def test_no_sources_does_not_create_edge_mask_or_modify_float32_data():
    from photozpy.photometry.local_background import LocalBackgroundOptions, _LocalBackground
    from photutils.aperture import CircularAnnulus
    data = np.full((50, 50), 60, dtype=np.float32)
    invalid = np.zeros(data.shape, bool)
    invalid[0, 0] = True
    data[0, 0] = np.nan
    original = data.copy()
    estimator = _LocalBackground(data, invalid, 3., LocalBackgroundOptions(min_pixels=5))
    stats, result = estimator.statistics(CircularAnnulus((20, 20), 5, 10), 3, 1.5)
    assert result['valid']
    assert result['area_fraction'] == pytest.approx(1)
    assert result['mean'] == 60
    np.testing.assert_array_equal(data, original)


@pytest.mark.parametrize('trial_valid,expected', [(True, 'unstable'), (False, 'insufficient_area')])
def test_stability_flags_do_not_select_best_snr_or_change_primary(tmp_path, monkeypatch, trial_valid, expected):
    from photozpy.photometry import LocalBackgroundOptions
    from photozpy.photometry import local_background as sky
    path, sources, _, regions = setup_image(tmp_path)
    review(regions)
    original = sky._LocalBackground.statistics

    def statistics(self, annulus, clip_sigma, growth):
        stats, result = original(self, annulus, clip_sigma, growth)
        if growth == 2.:
            result['mean'] += 100
            result['valid'] = trial_valid
        return stats, result

    monkeypatch.setattr(sky._LocalBackground, 'statistics', statistics)
    phot = Photometry(ImageFileCollection(tmp_path))
    result = phot.run_photometry(
        sources, regions=regions, verbose=False, background_options=LocalBackgroundOptions(),
    )
    assert result['background_stability'][0] == expected
    assert result['background_warning'][0]
    assert result['mean_bkg'][0].value == 60
    assert result['src'][0].value == pytest.approx(1000)
    assert result['status'][0] == 'measured'


@pytest.mark.parametrize('center', [(20., 20.), (2., 2.)])
def test_diagnostic_masks_match_actual_mean_and_exact_area(tmp_path, center):
    from photozpy.photometry.local_background import LocalBackgroundOptions, _LocalBackground
    from photozpy.photometry.background_diagnostics import BackgroundDiagnostic
    from photutils.aperture import CircularAperture, CircularAnnulus
    rng = np.random.default_rng(671)
    data = rng.normal(60, 1, (40, 40))
    invalid = np.zeros(data.shape, bool)
    invalid[20, 16] = True
    data[20, 26] = 500
    options = LocalBackgroundOptions(mask_neighbors=False, min_pixels=1)
    estimator = _LocalBackground(data, invalid, 3., options)
    estimator.distance = np.full(data.shape, np.inf)
    estimator.distance[25, 20] = 0
    ann = CircularAnnulus(center, 3, 7)
    stats, result = estimator.statistics(ann, 3, 1.5)
    diagnostic = BackgroundDiagnostic.capture(
        estimator, CircularAperture(center, 2), ann, stats, 1.5,
        {'source_id': 'test', 'field': 'DLS', 'filter': "SDSS_g'", 'image_path': 'test.fits'},
    )
    assert diagnostic.image[diagnostic.mean_labels == 1].mean() == pytest.approx(float(stats.mean))
    assert diagnostic.area_weights[diagnostic.area_labels == 1].sum() == pytest.approx(result['n_ann'])
    assert diagnostic.areas['after_clipping'] == pytest.approx(result['n_ann'])
    areas = list(diagnostic.areas.values())
    assert all(a >= b - 1e-8 for a, b in zip(areas, areas[1:]))
    if center == (20., 20.):
        assert set(np.unique(diagnostic.mean_labels)) == {0, 1, 2, 3, 4}
    else:
        assert diagnostic.areas['on_image'] < diagnostic.areas['geometric']
    export = tmp_path / 'diagnostic.fits'
    diagnostic.save(export)
    restored = BackgroundDiagnostic.load(export)
    np.testing.assert_array_equal(restored.mean_labels, diagnostic.mean_labels)
    np.testing.assert_array_equal(restored.area_labels, diagnostic.area_labels)
    assert restored.areas == diagnostic.areas
    assert restored.metadata == diagnostic.metadata
    with pytest.raises(OSError):
        restored.save(export)
    import matplotlib.pyplot as plt
    fig = restored.plot()
    assert len(fig.axes) == 3
    fig.savefig(tmp_path / 'diagnostic.png')
    plt.close(fig)


def test_public_diagnostic_is_a_snapshot_and_resets_after_plain_run(tmp_path):
    from photozpy.photometry import LocalBackgroundOptions
    path, sources, _, regions = setup_image(tmp_path)
    review(regions)
    phot = Photometry(ImageFileCollection(tmp_path))
    table = phot.run_photometry(sources, regions=regions, verbose=False, background_options=LocalBackgroundOptions())
    diagnostic = phot.get_background_diagnostic('standard', "SDSS_g'")
    assert diagnostic.areas['after_clipping'] == pytest.approx(table['n_ann'][0].value)
    assert diagnostic.areas['after_source_mask'] == pytest.approx(table['background_source_masked_area'][0].value)
    original = diagnostic.image.copy()
    # Plotting never reruns measurements or reads subsequently changed FITS pixels.
    with fits.open(path, mode='update') as hdul:
        hdul[0].data[:] = 500
    diagnostic.image[:] = 0
    snapshot = phot.get_background_diagnostic('standard', "SDSS_g'")
    np.testing.assert_array_equal(snapshot.image, original)
    import matplotlib.pyplot as plt
    fig = phot.plot_background_diagnostic('standard', "SDSS_g'")
    plt.close(fig)
    with pytest.raises(ValueError, match='No background diagnostic'):
        phot.get_background_diagnostic('absent', "SDSS_g'")
    phot.run_photometry(sources, regions=regions, verbose=False)
    with pytest.raises(ValueError, match='No background diagnostic'):
        phot.get_background_diagnostic('standard', "SDSS_g'")


def test_insufficient_background_still_has_inspectable_diagnostic(tmp_path):
    from photozpy.photometry import LocalBackgroundOptions
    path, sources, _, regions = setup_image(tmp_path)
    review(regions)
    phot = Photometry(ImageFileCollection(tmp_path))
    table = phot.run_photometry(
        sources, regions=regions, verbose=False,
        background_options=LocalBackgroundOptions(min_pixels=10000),
    )
    assert table['status'][0] == 'invalid_background'
    diagnostic = phot.get_background_diagnostic('standard', "SDSS_g'")
    assert diagnostic.areas['after_clipping'] < 10000


def test_empty_background_mask_can_be_exported_and_reloaded(tmp_path):
    from photozpy.photometry import LocalBackgroundOptions, BackgroundDiagnostic
    path, sources, _, regions = setup_image(tmp_path)
    review(regions)
    with fits.open(path, mode='update') as hdul:
        hdul['MASK'].data[50:80, 50:80] = 1
    phot = Photometry(ImageFileCollection(tmp_path))
    table = phot.run_photometry(sources, regions=regions, verbose=False,
                                background_options=LocalBackgroundOptions(mask_neighbors=False))
    assert table['status'][0] == 'invalid_background'
    diagnostic = phot.get_background_diagnostic('standard', "SDSS_g'")
    assert diagnostic.areas['after_clipping'] == 0
    diagnostic.save(tmp_path / 'empty.fits')
    restored = BackgroundDiagnostic.load(tmp_path / 'empty.fits')
    assert np.isnan(restored.metadata['mean'])
    assert restored.metadata['mean_pixel_count'] == 0


def test_source_object_plots_all_available_filters_in_telescope_order(tmp_path):
    from photozpy.photometry import LocalBackgroundOptions
    from dataclasses import replace
    path, sources, _, green = setup_image(tmp_path)
    review(green)
    red_path = tmp_path / 'Master_DLS_r.fits'
    with fits.open(path) as hdul:
        hdul[0].header['FILTER'] = "SDSS_r'"
        hdul.writeto(red_path)
    red = PhotometryRegions.generate(red_path, sources, root=tmp_path)
    review(red)
    # Reverse the measurement order to test canonical ordering in the figure.
    phot = Photometry(ImageFileCollection(tmp_path))
    phot.run_photometry(sources, regions=(red, green), verbose=False,
                        background_options=LocalBackgroundOptions())
    diagnostics = phot.get_background_diagnostics(sources)
    assert list(diagnostics) == ["SDSS_g'", "SDSS_r'"]
    assert list(phot.get_background_diagnostics(sources.source_info)) == list(diagnostics)
    assert "SDSS_u'" not in diagnostics
    import matplotlib.pyplot as plt
    fig = phot.plot_background_diagnostics(sources)
    assert len(fig.axes) == 8  # Three panels + one annotation axis per filter.
    assert "SDSS_g'" in fig.axes[0].get_title()
    assert "SDSS_r'" in fig.axes[4].get_title()
    fig.savefig(tmp_path / 'all_filters.png')
    plt.close(fig)
    diagnostics["SDSS_g'"].image[:] = -999
    assert not np.all(phot.get_background_diagnostics(sources)["SDSS_g'"].image == -999)
    # A same-name entry from another field must never leak into the selection.
    other = phot.get_background_diagnostics(sources)["SDSS_g'"]
    other.metadata['field'] = 'other'
    phot._background_diagnostics[('other', 'standard', "SDSS_g'")] = other
    assert all(d.metadata['field'] == 'DLS' for d in phot.get_background_diagnostics(sources).values())
    with pytest.raises(TypeError, match='source must'):
        phot.plot_background_diagnostics('standard')
    multiple = Sources(sources.source_info, replace(sources.source_info, source_name='another'))
    with pytest.raises(ValueError, match='exactly one'):
        phot.get_background_diagnostics(multiple)
    with pytest.raises(ValueError, match='No background diagnostics'):
        phot.get_background_diagnostics(Sources(replace(sources.source_info, source_name='absent')))


def test_batch_mask_plots_keep_same_name_fields_separate_and_report_missing(tmp_path, monkeypatch):
    from photozpy.photometry import LocalBackgroundOptions
    from photozpy.photometry import background_diagnostics as diagnostic_module
    from dataclasses import replace
    path, sources, _, regions = setup_image(tmp_path)
    review(regions)
    phot = Photometry(ImageFileCollection(tmp_path))
    phot.run_photometry(sources, regions=regions, verbose=False,
                        background_options=LocalBackgroundOptions())
    first = sources.source_info
    second = replace(first, observation_field=replace(first.observation_field, field_name='other'))
    missing = replace(first, source_name='rejected source')
    diagnostic = phot.get_background_diagnostic('standard', "SDSS_g'")
    diagnostic.metadata['field'] = 'other'
    phot._background_diagnostics[('other', 'standard', "SDSS_g'")] = diagnostic
    calls = []

    def plot(items):
        calls.append([(d.metadata['field'], d.metadata['source_id']) for d in items])
        return object()

    monkeypatch.setattr(diagnostic_module, 'plot_background_collection', plot)
    with pytest.warns(UserWarning, match='DLS/rejected source'):
        figures = phot.plot_background_diagnostics(Sources(first, missing, second))
    assert list(figures) == [('DLS', 'standard'), ('other', 'standard')]
    assert calls == [[('DLS', 'standard')], [('other', 'standard')]]
    assert figures[('DLS', 'standard')] is not figures[('other', 'standard')]
    calls.clear()
    with pytest.raises(ValueError, match='No background diagnostics'):
        phot.plot_background_diagnostics(Sources(missing, replace(missing, source_name='absent')))
    with pytest.raises(ValueError, match='No background diagnostics'):
        phot.plot_background_diagnostics(Sources())
    assert calls == []
