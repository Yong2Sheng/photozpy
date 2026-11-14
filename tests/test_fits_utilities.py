from photozpy.fits_utilities import (
    MissingWCSError,
    create_image_norm,
    get_fov_center_skycoord,
    get_wcs,
    has_wcs,
    plot_image,
    update_fits_header,
)
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from astropy.wcs import WCS
from astropy.visualization import AsinhStretch, ImageNormalize, SqrtStretch
from astropy.nddata import CCDData
from astropy.io import fits
from astropy.coordinates import SkyCoord
import pytest
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import matplotlib

matplotlib.use("Agg")


def test_update_header_preserves_all_hdus_and_original_file(tmp_path):
    source = tmp_path / "raw.fits"
    output = tmp_path / "processed" / "raw.fits"
    primary = np.arange(12, dtype=np.int16).reshape(3, 4)
    extension = np.arange(6, dtype=np.float32).reshape(2, 3)
    fits.HDUList([
        fits.PrimaryHDU(primary),
        fits.ImageHDU(extension, name="SCI"),
    ]).writeto(source)
    before = source.read_bytes()

    written = update_fits_header(
        source, {"OBJECT": ("field_a", "Observation field")},
        output_path=output,
    )

    assert written == output
    assert source.read_bytes() == before
    with fits.open(output) as hdul:
        assert len(hdul) == 2
        assert hdul[0].header["OBJECT"] == "field_a"
        assert hdul[0].header.comments["OBJECT"] == "Observation field"
        np.testing.assert_array_equal(hdul[0].data, primary)
        np.testing.assert_array_equal(hdul[1].data, extension)
        assert hdul[1].name == "SCI"


def test_update_header_requires_explicit_overwrite(tmp_path):
    source = _write_fits_image(tmp_path / "raw.fits")
    before = source.read_bytes()
    with pytest.raises(FileExistsError):
        update_fits_header(source, {"OBJECT": "field_a"})
    assert source.read_bytes() == before
    update_fits_header(source, {"OBJECT": "field_a"}, overwrite=True)
    assert fits.getheader(source)["OBJECT"] == "field_a"


def test_failed_header_update_does_not_damage_input(tmp_path):
    source = _write_fits_image(tmp_path / "raw.fits")
    before = source.read_bytes()
    with pytest.raises(ValueError):
        update_fits_header(source, {"OBJECT": object()}, overwrite=True)
    assert source.read_bytes() == before
    assert sorted(p.name for p in tmp_path.iterdir()) == ["raw.fits"]


def test_update_header_edits_only_requested_extension(tmp_path):
    source = tmp_path / "raw.fits"
    fits.HDUList([
        fits.PrimaryHDU(), fits.ImageHDU(np.ones((2, 2)), name="SCI")
    ]).writeto(source)
    update_fits_header(
        source, {"TESTKEY": 42}, hdu_index=1, overwrite=True
    )
    with fits.open(source) as hdul:
        assert "TESTKEY" not in hdul[0].header
        assert hdul[1].header["TESTKEY"] == 42


def test_update_header_preserves_unsigned_scaling_and_checksums(tmp_path):
    source = tmp_path / "unsigned.fits"
    data = np.array([[0, 32768, 65535]], dtype=np.uint16)
    fits.PrimaryHDU(data).writeto(source, checksum=True)
    update_fits_header(source, {"FIELD": "science_field"}, overwrite=True)
    with fits.open(source, checksum=True) as hdul:
        np.testing.assert_array_equal(hdul[0].data, data)
        assert hdul[0].header['BZERO'] == 32768
        assert hdul[0].verify_checksum() == 1
        assert hdul[0].verify_datasum() == 1


def _write_fits_image(
    path: Path,
    data: np.ndarray | None = None,
    *,
    header: fits.Header | None = None,
) -> Path:
    if data is None:
        data = np.arange(100, dtype=float).reshape(10, 10)

    fits.PrimaryHDU(data=data, header=header).writeto(path)

    return path


def _celestial_wcs_header() -> fits.Header:
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [10.0, 20.0]
    wcs.wcs.crpix = [5.0, 5.0]
    wcs.wcs.cdelt = [-0.0002777778, 0.0002777778]

    return wcs.to_header()


def test_create_image_norm_returns_default_asinh_normalization() -> None:
    data = np.arange(100, dtype=float).reshape(10, 10)

    norm = create_image_norm(data)

    assert isinstance(norm, ImageNormalize)
    assert isinstance(norm.stretch, AsinhStretch)
    assert norm.vmin < norm.vmax


def test_create_image_norm_accepts_supported_stretch_name() -> None:
    data = np.arange(100, dtype=float).reshape(10, 10)

    norm = create_image_norm(data, stretch="sqrt")

    assert isinstance(norm.stretch, SqrtStretch)


@pytest.mark.parametrize(
    ("data", "match"),
    [
        (np.arange(5, dtype=float), "two-dimensional"),
        (np.full((3, 3), np.nan), "no finite pixels"),
        (np.ones((3, 3)), "more than one finite display value"),
    ],
)
def test_create_image_norm_rejects_invalid_data(
    data: np.ndarray,
    match: str,
) -> None:
    with pytest.raises(ValueError, match=match):
        create_image_norm(data)


def test_create_image_norm_rejects_invalid_percentile_range() -> None:
    data = np.arange(100, dtype=float).reshape(10, 10)

    with pytest.raises(ValueError, match="min_percent.*max_percent"):
        create_image_norm(data, min_percent=99.0, max_percent=1.0)


def test_create_image_norm_rejects_unknown_stretch() -> None:
    data = np.arange(100, dtype=float).reshape(10, 10)

    with pytest.raises(ValueError, match="stretch"):
        create_image_norm(data, stretch="invalid")  # type: ignore[arg-type]


def test_plot_image_from_fits_returns_figure_without_saving_by_default(
    tmp_path: Path,
) -> None:
    fits_path = _write_fits_image(tmp_path / "preview.fits")

    fig, ax = plot_image(fits_path=fits_path, show_colorbar=False)

    assert isinstance(fig, Figure)
    assert isinstance(ax, Axes)
    assert ax.get_title() == "preview"
    assert not (tmp_path / "preview.png").exists()

    plt.close(fig)


def test_plot_image_from_ccddata_uses_metadata_for_title() -> None:
    ccd_data = CCDData(
        np.arange(100, dtype=float).reshape(10, 10),
        unit="adu",
        meta={"OBJECT": "target", "FILTER": "g"},
    )

    fig, ax = plot_image(ccd_data=ccd_data, show_colorbar=False)

    assert ax.get_title() == "target_g"

    plt.close(fig)


def test_plot_image_saves_when_save_dir_is_provided_and_closes(
    tmp_path: Path,
) -> None:
    fits_path = _write_fits_image(tmp_path / "science.fits")
    save_dir = tmp_path / "plots"

    result = plot_image(
        fits_path=fits_path,
        save_dir=save_dir,
        filename_suffix="_quicklook",
        show_colorbar=False,
        close=True,
    )

    assert result is None
    assert (save_dir / "science_quicklook.png").exists()


def test_plot_image_saves_to_explicit_save_path(tmp_path: Path) -> None:
    fits_path = _write_fits_image(tmp_path / "science.fits")
    save_path = tmp_path / "custom" / "image.png"

    result = plot_image(
        fits_path=fits_path,
        save_path=save_path,
        show_colorbar=False,
        close=True,
    )

    assert result is None
    assert save_path.exists()


def test_plot_image_forwards_norm_kwargs(tmp_path: Path) -> None:
    fits_path = _write_fits_image(tmp_path / "science.fits")

    fig, ax = plot_image(
        fits_path=fits_path,
        show_colorbar=False,
        norm_kwargs={"stretch": "sqrt", "max_percent": 99.0},
    )

    assert isinstance(ax.images[0].norm.stretch, SqrtStretch)

    plt.close(fig)


def test_plot_image_requires_exactly_one_image_input(tmp_path: Path) -> None:
    fits_path = _write_fits_image(tmp_path / "science.fits")
    ccd_data = CCDData(
        np.arange(100, dtype=float).reshape(10, 10),
        unit="adu",
    )

    with pytest.raises(ValueError, match="exactly one"):
        plot_image()

    with pytest.raises(ValueError, match="exactly one"):
        plot_image(fits_path=fits_path, ccd_data=ccd_data)


def test_plot_image_rejects_non_2d_image(tmp_path: Path) -> None:
    fits_path = _write_fits_image(
        tmp_path / "cube.fits",
        data=np.ones((2, 3, 4), dtype=float),
    )

    with pytest.raises(ValueError, match="2D image"):
        plot_image(fits_path=fits_path)


def test_has_wcs_and_get_wcs_for_celestial_header(tmp_path: Path) -> None:
    fits_path = _write_fits_image(
        tmp_path / "wcs.fits",
        header=_celestial_wcs_header(),
    )

    assert has_wcs(fits_path)
    assert isinstance(get_wcs(fits_path), WCS)


def test_get_wcs_raises_for_missing_celestial_wcs(tmp_path: Path) -> None:
    fits_path = _write_fits_image(tmp_path / "no_wcs.fits")

    assert not has_wcs(fits_path)

    with pytest.raises(MissingWCSError, match="No celestial WCS"):
        get_wcs(fits_path)


def test_get_fov_center_skycoord_returns_center_coordinate(
    tmp_path: Path,
) -> None:
    fits_path = _write_fits_image(
        tmp_path / "wcs.fits",
        header=_celestial_wcs_header(),
    )

    center = get_fov_center_skycoord(fits_path)

    assert isinstance(center, SkyCoord)


def test_get_fov_center_skycoord_rejects_non_2d_hdu(
    tmp_path: Path,
) -> None:
    fits_path = _write_fits_image(
        tmp_path / "cube_wcs.fits",
        data=np.ones((2, 3, 4), dtype=float),
        header=_celestial_wcs_header(),
    )

    with pytest.raises(ValueError, match="does not contain a 2D image"):
        get_fov_center_skycoord(fits_path)
