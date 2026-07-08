from photozpy.fits_utilities import (
    MissingWCSError,
    create_image_norm,
    get_fov_center_skycoord,
    get_wcs,
    has_wcs,
    plot_image,
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
