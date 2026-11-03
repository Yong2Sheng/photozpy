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
