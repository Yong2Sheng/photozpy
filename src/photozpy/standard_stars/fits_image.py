# Standard library
from pathlib import Path

# Third-party: Astropy coordinates / FITS / WCS
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord


def has_wcs(
    file: str | Path,
    hdu_index: int = 0,
    require_celestial: bool = True,
) -> bool:
    """
    Check whether a FITS HDU contains usable WCS information.

    This method attempts to construct an ``astropy.wcs.WCS`` object from the
    specified FITS header. By default, it requires that the WCS include a
    celestial component, such as RA/Dec axes.

    Parameters
    ----------
    file : str | Path
        Path to the FITS file.
    hdu_index : int, optional
        HDU index to inspect. Default is 0.
    require_celestial : bool, optional
        If True, require a celestial WCS component. If False, accept any parsed
        WCS with at least one axis. Default is True.

    Returns
    -------
    bool
        True if usable WCS information is present, otherwise False.

    Notes
    -----
    Any exception encountered while opening the FITS file or parsing the WCS is
    caught and treated as a False result.
    """

    try:
        with fits.open(file, memmap=False) as hdul:
            header = hdul[hdu_index].header
            w = WCS(header)

            if require_celestial:
                return bool(w.has_celestial)

            return w.naxis > 0

    except Exception:
        return False


def get_wcs(
    fits_path: Path | str,
    hdu_index: int = 0,
) -> WCS:

    # check is WCS is already solved
    if not has_wcs(fits_path, hdu_index=hdu_index):
        raise MissingWCSError(f"No celestial WCS found in {fits_path}")

    # access the data and wcs header
    with fits.open(fits_path, memmap=False) as hdul:
        header = hdul[hdu_index].header

    return WCS(header)


def get_fov_center_skycoord(
    fits_path: Path | str,
    hdu_index: int = 0,
) -> SkyCoord:
    """
    Compute the sky coordinate of the image center pixel using WCS.

    The geometric center of the image array is computed in pixel coordinates and
    then transformed into a celestial coordinate using the FITS WCS solution.

    Parameters
    ----------
    fits_path : Path | str
        FITS image file.
    hdu_index : int, optional
        HDU index containing the image data and WCS header. Default is 0.

    Returns
    -------
    SkyCoord
        Sky coordinate corresponding to the center pixel of the image.

    Raises
    ------
    MissingWCSError
        If the FITS file does not contain usable celestial WCS information.
    """

    # check is WCS is already solved
    if not has_wcs(fits_path, hdu_index=hdu_index):
        raise MissingWCSError(f"No celestial WCS found in {fits_path}")

    # access the data
    with fits.open(fits_path, memmap=False) as hdul:
        data = hdul[hdu_index].data

    # get wcs object
    wcs = get_wcs(
        fits_path=fits_path,
        hdu_index=hdu_index,
    )

    # get center pixel coordinate
    ny, nx = data.shape
    x_center = (nx - 1) / 2
    y_center = (ny - 1) / 2

    # transform the center coordinate to the SkyCoord
    center_coord = pixel_to_skycoord(x_center, y_center, wcs, origin=0)

    return center_coord


class MissingWCSError(ValueError):
    """
    Raised when a FITS file does not contain usable celestial WCS information.
    """
