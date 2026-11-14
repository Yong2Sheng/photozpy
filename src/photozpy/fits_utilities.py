# Standard library
from collections.abc import Iterable, Mapping
import os
from pathlib import Path
import tempfile
from typing import Any, Literal

# Third-party
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import sigma_clip
from astropy.time import Time
from ccdproc import ImageFileCollection
from astropy.visualization import (
    AsinhStretch,
    ImageNormalize,
    LinearStretch,
    LogStretch,
    SqrtStretch,
)
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import ArrayLike


ImageStretchName = Literal["asinh", "sqrt", "log", "linear"]


def update_fits_header(
    fits_path: str | Path,
    header_values: Mapping[str, Any],
    *,
    output_path: str | Path | None = None,
    hdu_index: int | str | tuple[str, int] = 0,
    overwrite: bool = False,
    remove_keys: Iterable[str] = (),
) -> Path:
    """Write header cards in one HDU while preserving the complete FITS file.

    Values may be ordinary FITS values or ``(value, comment)`` tuples.
    remove_keys are deleted from the selected header before applying new values;
    missing keys are ignored. Use this for replacing an obsolete WCS completely.
    An empty mapping copies the complete file without changing header cards.
    No time conversion, filter mapping, or image-data calculation is performed.
    Only the selected HDU's header is edited; all other HDUs are retained.
    If the input contains checksum cards, checksums are regenerated for every
    output HDU so header edits cannot leave stale integrity metadata.

    If output_path is None, replace the input file (requires overwrite=True).
    Otherwise write a separate output file. The complete output is first
    written to a temporary file in the destination directory, so a failed FITS
    serialization does not leave a partially overwritten input. Existing output
    files are replaced only when overwrite=True. Return the written path.
    """
    source = Path(fits_path).expanduser().resolve()
    destination = (
        source if output_path is None
        else Path(output_path).expanduser().resolve()
    )
    if not source.is_file():
        raise FileNotFoundError(source)
    if destination.exists() and not overwrite:
        raise FileExistsError(destination)

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = None
    try:
        with fits.open(
            source, mode="readonly", memmap=False, do_not_scale_image_data=True
        ) as hdul:
            for key in remove_keys:
                hdul[hdu_index].header.remove(key, ignore_missing=True, remove_all=True)
            for key, value in header_values.items():
                hdul[hdu_index].header[key] = value
            # Preserve the output suffix so Astropy retains compressed formats.
            descriptor, temporary_name = tempfile.mkstemp(
                prefix=f".{destination.name}.",
                suffix="".join(destination.suffixes),
                dir=destination.parent,
            )
            os.close(descriptor)
            temporary_path = Path(temporary_name)
            checksums_present = any(
                'CHECKSUM' in hdu.header or 'DATASUM' in hdu.header
                for hdu in hdul
            )
            hdul.writeto(
                temporary_path, overwrite=True, checksum=checksums_present
            )

        if overwrite:
            os.replace(temporary_path, destination)
        else:
            # A hard link publishes the complete file without overwriting an
            # output that another process may have created since our check.
            os.link(temporary_path, destination)
        return destination
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


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


class MissingWCSError(ValueError):
    """
    Raised when a FITS file does not contain usable celestial WCS information.
    """


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
    if data is None or data.ndim != 2:
        raise ValueError(
            f"HDU {hdu_index} in {fits_path} does not contain a 2D image."
        )
    ny, nx = data.shape
    x_center = (nx - 1) / 2
    y_center = (ny - 1) / 2

    # transform the center coordinate to the SkyCoord
    center_coord = pixel_to_skycoord(x_center, y_center, wcs, origin=0)

    return center_coord


def create_image_norm(
    data: ArrayLike,
    *,
    stretch: ImageStretchName = "asinh",
    min_percent: float = 0.5,
    max_percent: float = 99.5,
    clip_sigma: float = 3.0,
    clip_maxiters: int | None = 5,
    clip: bool = True,
) -> ImageNormalize:
    """
    Create a display normalization for a 2D FITS image.

    The default settings are tuned for quick-look plots that should make faint
    sources easier to see. The cut levels are estimated from finite,
    sigma-clipped pixels so that saturated stars, cosmic rays, and other bright
    outliers do not dominate the contrast.

    Parameters
    ----------
    data : array-like
        Two-dimensional image data.
    stretch : {"asinh", "sqrt", "log", "linear"}, optional
        Stretch function used by `astropy.visualization.ImageNormalize`.
        The default is ``"asinh"``, which is usually a good compromise for
        faint structure and bright compact sources.
    min_percent : float, optional
        Lower percentile used to define ``vmin`` after sigma clipping.
        Default is 0.5.
    max_percent : float, optional
        Upper percentile used to define ``vmax`` after sigma clipping.
        Default is 99.5.
    clip_sigma : float, optional
        Sigma threshold passed to `astropy.stats.sigma_clip`. Default is 3.0.
    clip_maxiters : int or None, optional
        Maximum number of sigma-clipping iterations. Default is 5.
    clip : bool, optional
        Whether values outside ``vmin`` and ``vmax`` are clipped by the
        returned normalization. Default is True.

    Returns
    -------
    astropy.visualization.ImageNormalize
        Normalization object that can be passed to ``imshow(..., norm=norm)``.

    Raises
    ------
    ValueError
        If ``data`` is not two-dimensional, contains no finite pixels, has an
        invalid percentile range, or uses an unsupported stretch name.

    Notes
    -----
    This function only controls display scaling. It does not modify the input
    image data and should not be used for photometric measurements.
    """

    image_data = np.asarray(data)

    if image_data.ndim != 2:
        raise ValueError(
            "'data' must be a two-dimensional image array."
        )

    if not 0 <= min_percent < max_percent <= 100:
        raise ValueError(
            "'min_percent' and 'max_percent' must satisfy "
            "0 <= min_percent < max_percent <= 100."
        )

    finite_pixels = image_data[np.isfinite(image_data)]
    if finite_pixels.size == 0:
        raise ValueError("Image data contains no finite pixels.")

    clipped_pixels = sigma_clip(
        finite_pixels,
        sigma=clip_sigma,
        maxiters=clip_maxiters,
        masked=True,
    ).compressed()

    # Extremely small or pathological images can be fully clipped. Fall back to
    # all finite pixels so the function remains usable for quick-look plotting.
    if clipped_pixels.size == 0:
        clipped_pixels = finite_pixels

    vmin = float(np.percentile(clipped_pixels, min_percent))
    vmax = float(np.percentile(clipped_pixels, max_percent))

    # A constant image has no useful contrast range. Use the full finite range
    # as a final fallback before raising a clear error.
    if vmax <= vmin:
        vmin = float(np.min(finite_pixels))
        vmax = float(np.max(finite_pixels))

    if vmax <= vmin:
        raise ValueError(
            "Image data must contain more than one finite display value."
        )

    stretches = {
        "asinh": AsinhStretch(),
        "sqrt": SqrtStretch(),
        "log": LogStretch(),
        "linear": LinearStretch(),
    }

    try:
        stretch_object = stretches[stretch]
    except KeyError as exc:
        raise ValueError(
            "'stretch' must be one of: 'asinh', 'sqrt', 'log', 'linear'."
        ) from exc

    return ImageNormalize(
        vmin=vmin,
        vmax=vmax,
        stretch=stretch_object,
        clip=clip,
    )


def plot_image(
    *,
    fits_path: str | Path | None = None,
    ccd_data: CCDData | None = None,
    hdu_index: int = 0,
    save_path: str | Path | None = None,
    save_dir: str | Path | None = None,
    filename_suffix: str = "",
    cmap: str = "gray",
    origin: str = "lower",
    interpolation: str = "nearest",
    figsize: tuple[float, float] = (8.0, 8.0),
    dpi: int = 300,
    show_colorbar: bool = True,
    close: bool = False,
    norm: ImageNormalize | None = None,
    norm_kwargs: dict[str, Any] | None = None,
) -> tuple[Figure, Axes] | None:
    """
    Plot a simple quick-look image from a FITS file or CCDData object.

    This function intentionally does not draw coordinates, markers, apertures,
    regions, labels, or zoomed cutouts. Domain-specific overlays should be
    handled by the caller or by module-specific plotting functions.

    Parameters
    ----------
    fits_path : str, pathlib.Path, or None, optional
        Path to the FITS image. Provide either ``fits_path`` or ``ccd_data``,
        but not both.
    ccd_data : astropy.nddata.CCDData or None, optional
        CCDData object containing a 2D image. Provide either ``fits_path`` or
        ``ccd_data``, but not both.
    hdu_index : int, optional
        HDU index used when reading ``fits_path``. Default is 0.
    save_path : str, pathlib.Path, or None, optional
        Full output path for the saved figure. If provided, this takes priority
        over ``save_dir`` and ``filename_suffix``.
    save_dir : str, pathlib.Path, or None, optional
        Directory where the figure is saved. If provided and ``save_path`` is
        None, the output filename is generated from the input image stem and
        ``filename_suffix``. If both ``save_path`` and ``save_dir`` are None,
        the figure is not saved.
    filename_suffix : str, optional
        Suffix appended to the default output filename before ``.png`` when
        ``save_dir`` is provided and ``save_path`` is None.
    cmap : str, optional
        Matplotlib colormap. Default is ``"gray"``.
    origin : str, optional
        Image origin passed to ``imshow``. Default is ``"lower"``.
    interpolation : str, optional
        Interpolation mode passed to ``imshow``. Default is ``"nearest"``.
    figsize : tuple of float, optional
        Figure size in inches. Default is ``(8.0, 8.0)``.
    dpi : int, optional
        Resolution used when saving the figure. Default is 300.
    show_colorbar : bool, optional
        Whether to draw a colorbar. Default is True.
    close : bool, optional
        If True, close the figure before returning. This is useful for batch
        plotting after saving figures. Default is False.
    norm : astropy.visualization.ImageNormalize or None, optional
        Precomputed normalization. If None, ``create_image_norm`` is called.
    norm_kwargs : dict or None, optional
        Keyword arguments forwarded to ``create_image_norm`` when ``norm`` is
        None. For example: ``{"stretch": "asinh", "max_percent": 99.5}``.

    Returns
    -------
    tuple[matplotlib.figure.Figure, matplotlib.axes.Axes] or None
        Returns ``(fig, ax)`` when ``close`` is False. Returns None when
        ``close`` is True.

    Raises
    ------
    ValueError
        If neither or both of ``fits_path`` and ``ccd_data`` are provided, or
        if the image data is not two-dimensional.
    """

    if (fits_path is None) == (ccd_data is None):
        raise ValueError(
            "Provide exactly one of 'fits_path' or 'ccd_data'."
        )

    if fits_path is not None:
        fits_path = Path(fits_path)
        image_data = np.asarray(fits.getdata(fits_path, ext=hdu_index))
        default_stem = fits_path.stem
    else:
        image_data = np.asarray(ccd_data.data)
        header = getattr(ccd_data, "header", None)
        if header is None:
            header = ccd_data.meta

        object_name = header.get("OBJECT", "image")
        filter_name = header.get("FILTER", "unknown_filter")
        default_stem = f"{object_name}_{filter_name}"

    if image_data.ndim != 2:
        raise ValueError("plot_image() expects a 2D image array.")

    if norm is None:
        if norm_kwargs is None:
            norm_kwargs = {}

        norm = create_image_norm(
            image_data,
            **norm_kwargs,
        )

    fig, ax = plt.subplots(figsize=figsize)
    image_artist = ax.imshow(
        image_data,
        cmap=cmap,
        origin=origin,
        interpolation=interpolation,
        norm=norm,
    )

    if show_colorbar:
        fig.colorbar(image_artist, ax=ax)

    ax.set_title(default_stem)

    should_save = save_path is not None or save_dir is not None

    if should_save:
        if save_path is None:
            if save_dir is None:
                raise ValueError(
                    "'save_dir' must be provided when 'save_path' is None."
                )
            save_dir = Path(save_dir)
            save_dir.mkdir(parents=True, exist_ok=True)
            save_path = save_dir / f"{default_stem}{filename_suffix}.png"
        else:
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)

        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")

    if close:
        plt.close(fig)
        return None

    return fig, ax


# Header operations migrated from HeaderManipulation. Collection-based
# signatures and existing scientific behavior are retained in this migration.
def check_headers(image_collection, **headers_values):
    """
    Check if the headers and their values are consistent in the image collection

    Parameters
    ----------
    image_collection: ccdproc.ImageFileCollection; the collection of images
    headers_values: dict; the header names and corresponding values. The values should be in tuple.
                          If there is only one value for the header, use (value,)

    Returns
    -------
    final : boolean; True or False
    """

    results = []
    for header, value in headers_values.items():
        # get the header values from the image collection
        values_in_collection = image_collection.values(header)
        # remove the repeated values to speed up
        values_in_collection = [*set(values_in_collection)]

        result = []
        for val in value:
            if val in values_in_collection:
                result += [True]
            else:
                result += [False]

        results += [all(result)]

    final = all(results)

    return final


def correct_headers(input_collection, save_location="",
                    overwrite=True, time_transformation=False, **headers_values):
    """
    Change or add headers and values.

    Parameters
    ----------
    input_collection: ccdproc.ImageFileCollection; the collection of images.
    save_location: str; the location to save the new files.
    overwrite: boolean; overwrite the files or not if the files already exist.
    time_transformation: boolean; transformation the time format from UT1 to MJD and add a MJD-OBS header to meet the fits standard
    headers_values: dict; the dictionary of the headers and their values.

    Returns
    ------
    new_collection: ccdproc.ImageFileCollection; The new image collection.
    """

    if not headers_values:  # python evaluate empty object to be False
        raise ValueError(
            "Please provide at least one header and one header values!")

    # create the save_location if not exist
    if save_location != "":
        Path(save_location).mkdir(parents=False, exist_ok=True)

    for hdu in input_collection.hdus(
            overwrite=overwrite, save_location=save_location):
        for key, value in headers_values.items():
            hdu.header[key] = value
            if time_transformation:
                date_time = Time(
                    hdu.header["DATE-OBS"], format="isot", scale="ut1")
                mjd_time = date_time.mjd
                hdu.header["MJD-OBS"] = mjd_time

    # Here I want to return a new collection of the modified files because
    # the summary of the input image collection won't be updated by itself
    if save_location == "":
        new_collection = ImageFileCollection(
            location=input_collection.location,
            filenames=input_collection.files)
    else:
        new_collection = ImageFileCollection(
            location=save_location, filenames=input_collection.files)

    return new_collection


def map_filters(name, **filter_dict):
    """
    Map filter names to standard filter names.

    Parameters
    ----------
    name: the input filter name to be matched to a standard one.
    filter_dict: customized filter mapping dictionary.

    Return
    ------
    mapped_filter: str; the name of the mapped standard filter.
    """

    if filter_dict == {}:
        # use the default mappers
        mappers = {"Sloan g'2": "SDSS_g'",
                   "Sloan r'2": "SDSS_r'",
                   "Sloan i'2": "SDSS_i'",
                   "Sloan z'2": "SDSS_z'",
                   "Sloan_g'2": "SDSS_g'",
                   "Sloan_r'2": "SDSS_r'",
                   "Sloan_i'2": "SDSS_i'",
                   "Sloan_z'2": "SDSS_z'",
                   "SDSS g": "SDSS_g'",
                   "SDSS r": "SDSS_r'",
                   "SDSS i": "SDSS_i'",
                   "SDSS z": "SDSS_z'",
                   "SDSS_g": "SDSS_g'",
                   "SDSS_r": "SDSS_r'",
                   "SDSS_i": "SDSS_i'",
                   "SDSS_z": "SDSS_z'",
                   "BESSEL B": "Bessel_B",
                   "BESSEL V": "Bessel_V",
                   "BESSEL R": "Bessel_R",
                   "BESSEL I": "Bessel_I",
                   "B": "ubb",
                   "UVM2": "um2",
                   "U": "uuu",
                   "V": "uvv",
                   "UVW1": "uw1",
                   "UVW2": "uw2"}
    else:
        mappers = filter_dict

    mapped_filters = [
        value for key,
        value in mappers.items() if key == name]
    # remove the duplicated filter names
    mapped_filters = [*set(mapped_filters)]

    if len(mapped_filters) == 0:
        raise ValueError(
            f"Mapping the input filter failed! No mapped filter {name} was found!")
    elif len(mapped_filters) > 1:
        raise ValueError(
            f"Mapping the input filter failed! More than two mapped filters are found: {mapped_filters}")
    else:
        return mapped_filters[0]


def get_header_values(image_collection, header):
    """
    Get the header values from an image collection. The repeated values will be removed.

    Parameters
    ----------
    image_collection
    headers

    Returns
    -------
    header_values: list
    """

    return image_collection.values(header, unique=True)


def edit_headers_with_fitsio(image_collection, save_location="",
                             additional_name=None, overwrite=True, **headers_values, ):
    """
    This edits the header of multi-extension fits files since the iteration of ImageCollection doesn't work on
    multi-extension fits files.

    Parameters
    ----------
    image_collection
    save_location
    headers_values: dict; the dictionary of the headers you want to add/edit
    additional_name: str; the additional name appending to the original file name

    Returns
    -------
    None
    """

    file_paths = image_collection.files_filtered(include_path=True)

    for file_path in file_paths:

        file_path = Path(file_path)
        if save_location != "":
            save_location = Path(save_location)
            if additional_name is None:
                save_location = save_location / file_path.name
            else:
                file_name = file_path.stem + "_master" + file_path.suffix
                save_location = save_location / file_name
        else:
            if additional_name is None:
                save_location = file_path
            else:
                file_name = file_path.stem + "_master" + file_path.suffix
                save_location = file_path.parent / file_name

        hdul = fits.open(file_path)
        for key, value in headers_values.items():
            hdul[0].header[key] = value
        hdul.writeto(save_location, overwrite=overwrite)
        hdul.close()
