# Standard library
from collections.abc import Mapping
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal

# Third-party: NumPy / pandas
import numpy as np
import numpy.typing as npt
import pandas as pd

# Third-party: table / text formatting
import tables as tb
from tabulate import tabulate

# Third-party: plotting
import matplotlib.pyplot as plt

# Third-party: Astropy core
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.units import Quantity
from astropy.visualization import simple_norm
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel

# Third-party: HEALPix
from astropy_healpix import HEALPix

# Local imports
from ..fits_utilities import get_fov_center_skycoord, get_wcs
from ..sources import FilterMagnitude, ObservationField, SourceInfo, Sources
from ..telescope import Telescope
from ..telescope.filterset import UnknownFilterAliasError
from ..interactive_ds9.adapter import DS9Adapter


@dataclass(slots=True, frozen=True, kw_only=True)
class CatalogMetaData:
    """
    Store structural metadata for a compiled standard-star catalog table.

    This data object records the catalog file location, the resolved PyTables
    table node information, and the schema-level metadata needed to interpret
    the catalog content. It is intended to describe the catalog table itself,
    not to hold the table rows in memory.

    Attributes
    ----------
    catalog_path : Path
        Path to the compiled HDF5 catalog file.
    table_node_path : str
        Full PyTables node path of the catalog table inside the HDF5 file.
    table_node_name : str
        Name of the table node used to locate the catalog table.
    table_title : str
        Human-readable table title stored in the PyTables table metadata.
    table_columns : tuple[str, ...]
        Column names defined in the catalog table.
    bucket_size : int
        Number of sources stored in each spatial bucket of the compiled catalog.
    ra_unit : str
        Unit string used for right ascension values in the catalog.
    dec_unit : str
        Unit string used for declination values in the catalog.
    mag_system : str
        Photometric magnitude system used by the catalog.
    nside : int
        HEALPix NSIDE value used when building the spatial index.
    order : int
        Ordering or hierarchy parameter stored in the compiled catalog metadata.
    source : str
        Source survey or catalog name from which the standard stars were compiled.
    version : str
        Version string of the compiled catalog schema or dataset.
    """

    catalog_path: Path
    table_node_path: str
    table_node_name: str
    table_title: str
    table_columns: tuple[str, ...]
    bucket_size: int
    ra_unit: str
    dec_unit: str
    mag_system: str
    nside: int
    order: int
    source: str
    version: str


class StandardStarCatalog:
    """
    Represent a compiled standard-star catalog and its table-level metadata.

    This class resolves the target table node inside a compiled HDF5 catalog,
    reads the table metadata, and stores the result as a ``CatalogMetaData``
    instance for later use by higher-level matching and calibration routines.

    Parameters
    ----------
    catalog_path : str | Path
        Path to the compiled HDF5 standard-star catalog.
    table_node_name : str, optional
        Name of the PyTables table node to resolve. Default is ``"std"``.

    Attributes
    ----------
    catalog_path : Path
        Path to the compiled HDF5 catalog file.
    table_node_name : str
        Name of the target PyTables table node.
    metadata : CatalogMetaData
        Parsed metadata describing the resolved catalog table.
    table_node_path : str
        Full PyTables node path of the resolved catalog table.
    table_title : str
        Human-readable title of the resolved catalog table.
    """

    def __init__(
        self,
        catalog_path: str | Path,
        table_node_name: str = "std",
    ) -> None:
        """
        Initialize the standard-star catalog interface from a compiled HDF5 file.

        Parameters
        ----------
        catalog_path : str | Path
            Path to the compiled HDF5 standard-star catalog.
        table_node_name : str, optional
            Name of the PyTables table node to resolve. Default is ``"std"``.
        """

        self.catalog_path = Path(catalog_path)
        self.table_node_name: str = table_node_name

        # read the metadata
        self.metadata = self._read_catalog_metadata()

        # store some useful strings
        self.table_node_path: str = self.metadata.table_node_path
        self.table_title: str = self.metadata.table_title
        self._matched_stars: pd.DataFrame | None = None

    @property
    def matched_stars(self) -> pd.DataFrame:
        """Return the last completed query's fully filtered candidate table.

        Includes all original catalog and optional pixel-position columns, in
        query row order with the original index. Source IDs are decoded from
        UTF-8 bytes for notebook display. Magnitudes and errors are unchanged.
        Both find_standard_stars and find_standard_sources populate this table;
        to_sources alone does not perform a query or replace it.

        A defensive copy prevents notebook edits from changing the stored query
        result. An empty successful query returns an empty DataFrame with its
        columns intact. Before the first successful query, raises RuntimeError;
        a failed query leaves the previous snapshot available.
        """
        if self._matched_stars is None:
            raise RuntimeError("No matched-star table is available; run a catalog query first.")
        stars = self._matched_stars.copy(deep=True)
        if "id_name" in stars:
            stars["id_name"] = stars["id_name"].map(
                lambda value: value.decode("utf-8") if isinstance(value, bytes) else value
            )
        return stars

    @staticmethod
    def resolve_magnitude_columns(
        columns: tuple[str, ...],
        telescope: Telescope,
    ) -> tuple[dict[str, str], dict[str, str]]:
        """Resolve catalog columns through Telescope's existing filter aliases.

        Return magnitude and uncertainty mappings in telescope filter order.
        Uncertainty columns use the compiled catalogs' ``<magnitude>_err``
        convention. Unrecognized columns are ignored; ambiguous aliases and
        catalogs with no matching filters raise ValueError. No approximate
        string matching or photometric-system transformation is performed.
        """
        matched: dict[str, str] = {}
        for column in columns:
            try:
                canonical = telescope.resolve_filter_alias(column)
            except UnknownFilterAliasError:
                continue
            if canonical in matched:
                raise ValueError(
                    f"Multiple catalog columns match {canonical!r}: "
                    f"{matched[canonical]!r}, {column!r}. "
                    "Provide an explicit magnitude_columns override."
                )
            matched[canonical] = column
        if not matched:
            raise ValueError("No catalog magnitude columns match the telescope filters.")
        magnitudes = {
            name: matched[name] for name in telescope.filters() if name in matched
        }
        errors = {
            name: column + "_err"
            for name, column in magnitudes.items() if column + "_err" in columns
        }
        return magnitudes, errors

    def to_sources(
        self,
        stars: pd.DataFrame,
        *,
        observation_field: ObservationField,
        magnitude_columns: Mapping[str, str] | None = None,
        magnitude_system: Literal["AB", "Vega"],
        error_columns: Mapping[str, str] | None = None,
        coordinate_frame: str = "icrs",
    ) -> Sources:
        """Attach catalog reference magnitudes to individual standard sources.

        Parameters
        ----------
        stars : pandas.DataFrame
            Selected catalog rows, normally returned by ``find_standard_stars``.
            Must contain ``id_name``, ``ra``, and ``dec``. Byte-string source
            identifiers are decoded as UTF-8; row order is preserved.
        observation_field : ObservationField
            Shared field definition. Its Telescope instance is reused by every
            reference FilterMagnitude.
        magnitude_columns : Mapping[str, str] or None, optional
            Telescope filter names or aliases mapped to catalog magnitude columns.
            By default, resolve catalog columns through Telescope aliases and
            discover corresponding _err columns. An explicit mapping overrides
            discovery. At least one filter is required; aliases must be unique.
        magnitude_system : {"AB", "Vega"}
            Explicit system of the selected columns. Catalogs may contain mixed
            systems, so this is not inferred from table-level metadata.
        error_columns : Mapping[str, str] or None, optional
            Filter names or aliases mapped to catalog uncertainty columns. Only
            filters in magnitude_columns may be supplied. Missing errors remain
            NaN; no uncertainties are estimated here.
        coordinate_frame : str, optional
            Frame of the catalog RA/Dec, default ICRS as in the existing catalog
            footprint query. Angular units come from catalog metadata. No frame
            is inferred from the observed image's WCS.

        Returns
        -------
        Sources
            One SourceInfo per row, with reference magnitudes attached. Both
            instrumental and calibrated slots remain None. Empty input returns
            an empty collection after column and configuration validation.

        Notes
        -----
        This adapts the original notebook's coordinate and magnitude extraction;
        it does not centroid, calibrate, convert photometric systems, or write
        files. Catalog values (including NaN and any numeric sentinels) are
        preserved. Catalog-specific missing-value rules require a separate,
        explicit preprocessing step. Duplicate field/source identities fail
        rather than silently replacing a star.
        """
        if not isinstance(stars, pd.DataFrame):
            raise TypeError("stars must be a pandas DataFrame.")
        if not isinstance(observation_field, ObservationField):
            raise TypeError("observation_field must be an ObservationField.")
        if magnitude_system not in ("AB", "Vega"):
            raise ValueError("magnitude_system must be 'AB' or 'Vega'.")
        if magnitude_columns is None:
            magnitude_columns, discovered_errors = self.resolve_magnitude_columns(
                tuple(stars.columns), observation_field.telescope,
            )
            if error_columns is None:
                error_columns = discovered_errors
        if not magnitude_columns:
            raise ValueError("magnitude_columns must contain at least one filter.")

        telescope = observation_field.telescope
        mag_columns: dict[str, str] = {}
        for filter_name, column in magnitude_columns.items():
            canonical = telescope.resolve_filter_alias(filter_name)
            if canonical in mag_columns:
                raise ValueError(f"Duplicate canonical filter: {canonical!r}.")
            mag_columns[canonical] = column

        err_columns: dict[str, str] = {}
        for filter_name, column in (error_columns or {}).items():
            canonical = telescope.resolve_filter_alias(filter_name)
            if canonical not in mag_columns:
                raise ValueError(f"Error filter {canonical!r} has no magnitude column.")
            if canonical in err_columns:
                raise ValueError(f"Duplicate error filter: {canonical!r}.")
            err_columns[canonical] = column

        required = {"id_name", "ra", "dec", *mag_columns.values(), *err_columns.values()}
        missing = required.difference(stars.columns)
        if missing:
            raise KeyError(f"Missing catalog columns: {sorted(missing)}")

        infos: list[SourceInfo] = []
        references: list[FilterMagnitude] = []
        for _, row in stars.iterrows():
            identifier = row["id_name"]
            if isinstance(identifier, bytes):
                name = identifier.decode("utf-8")
            elif isinstance(identifier, str):
                name = identifier
            else:
                raise TypeError("Catalog id_name values must be strings or bytes.")

            ra, dec = float(row["ra"]), float(row["dec"])
            if not np.all(np.isfinite([ra, dec])):
                raise ValueError(f"Non-finite catalog coordinate for {name!r}.")
            infos.append(SourceInfo(
                source_name=name,
                source_role="standard",
                source_coordinate=SkyCoord(
                    ra=ra, dec=dec,
                    unit=(u.Unit(self.metadata.ra_unit), u.Unit(self.metadata.dec_unit)),
                    frame=coordinate_frame,
                ),
                observation_field=observation_field,
            ))
            references.append(FilterMagnitude(
                telescope=telescope,
                filter_names=tuple(mag_columns),
                magnitudes=[float(row[column]) for column in mag_columns.values()],
                errors=[
                    float(row[err_columns[name]]) if name in err_columns else np.nan
                    for name in mag_columns
                ],
                significance=None,
                magnitude_system=magnitude_system,
            ))

        sources = Sources.from_list(infos)
        for single, reference in zip(sources, references):
            single.set_magnitudes(kind="reference", magnitudes=reference)
        return sources

    def find_standard_sources(
        self,
        fits_path: str | Path,
        *,
        observation_field: ObservationField,
        magnitude_columns: Mapping[str, str] | None = None,
        magnitude_system: Literal["AB", "Vega"],
        error_columns: Mapping[str, str] | None = None,
        hdu_index: int = 0,
        edge_margin_pix: float = 20.0,
        mag_limit: float | None = None,
    ) -> Sources:
        """Query one field image and return standards with reference results.

        Retains ``find_standard_stars`` and its original footprint, edge, and
        optional all-selected-filters magnitude cuts. The selected FITS HDU must
        have FIELD matching observation_field.field_name. Magnitude/error column
        mappings and the explicit AB/Vega system follow ``to_sources``.
        By default, columns and errors are discovered using Telescope aliases.
        If mag_limit is supplied, it applies to every discovered magnitude
        column, including telescope filters not present in the query image.

        Pixel positions in ``self.final_stars`` are initial WCS projections,
        not centroids. Returned sources use catalog ICRS coordinates, matching
        the existing query. This method does not read FWHM, create regions,
        modify an existing Sources collection, or write files.
        """
        if not isinstance(observation_field, ObservationField):
            raise TypeError("observation_field must be an ObservationField.")
        if magnitude_columns is None:
            magnitude_columns, discovered_errors = self.resolve_magnitude_columns(
                self.metadata.table_columns, observation_field.telescope,
            )
            if error_columns is None:
                error_columns = discovered_errors
        # Validate mappings on an empty table before querying the catalog.
        columns = {"id_name", "ra", "dec", *magnitude_columns.values()}
        columns.update((error_columns or {}).values())
        self.to_sources(
            pd.DataFrame(columns=sorted(columns)),
            observation_field=observation_field,
            magnitude_columns=magnitude_columns,
            error_columns=error_columns,
            magnitude_system=magnitude_system,
        )
        header = fits.getheader(fits_path, ext=hdu_index)
        if header.get("FIELD") != observation_field.field_name:
            raise ValueError("Image FIELD does not match observation_field.field_name.")

        stars = self.find_standard_stars(
            fits_path=fits_path,
            hdu_index=hdu_index,
            edge_margin_pix=edge_margin_pix,
            add_pixel_columns=True,
            mag_limit=mag_limit,
            mag_columns=list(magnitude_columns.values()),
        )
        sources = self.to_sources(
            stars,
            observation_field=observation_field,
            magnitude_columns=magnitude_columns,
            error_columns=error_columns,
            magnitude_system=magnitude_system,
        )
        self._matched_stars = stars.copy(deep=True)
        return sources

    def preview_standard_stars(
        self,
        fits_path: str | Path,
        sources: Sources,
        *,
        adapter: DS9Adapter | None = None,
    ) -> DS9Adapter:
        """Inspect catalog candidates interactively in a dedicated DS9 window.

        Display the primary image and label each standard's initial catalog
        position with a blue cross. All sources must belong to the image's
        FIELD. No centroiding, acceptance decision, or region-file write takes
        place. DS9 remains open and this method returns immediately after load.

        A supplied adapter selects an existing XPA target; otherwise use
        photozpy-standard-preview. This target's current regions are replaced.
        Later aperture review will record edits separately; moving these
        inspection crosses does not update SourceInfo or reference results.
        The legacy matplotlib plot method remains available for compatibility.
        """
        header = fits.getheader(fits_path, ext=0)
        lines = ["# Region file format: DS9 version 4.1", "icrs"]
        for info in sources.source_infos:
            if info.source_role != "standard":
                raise ValueError("Standard preview requires standard sources.")
            if info.observation_field.field_name != header.get("FIELD"):
                raise ValueError("Source field does not match image FIELD.")
            coordinate = info.source_coordinate.icrs
            # Labels are display-only; keep region syntax out of catalog names.
            label = info.source_name.replace("{", "(").replace("}", ")")
            label = label.replace("\\", " ").replace("\n", " ").replace("\r", " ")
            lines.append(
                f"point({coordinate.ra.deg:.10f},{coordinate.dec.deg:.10f}) "
                f"# point=cross 15 color=blue width=2 text={{{label}}} "
                "select=0 edit=0 move=0 delete=0"
            )

        viewer = adapter if adapter is not None else DS9Adapter(
            software_name="photozpy-standard-preview",
        )
        viewer.ensure_running()
        viewer.load_fits(fits_path, hdu=0)
        viewer.verify_image_loaded()
        viewer.clear_regions()
        viewer.add_regions("\n".join(lines) + "\n")
        return viewer

    def print_metadata(
        self,
    ) -> None:

        metadata_dict = asdict(self.metadata)

        # 对超长字段做特殊处理
        metadata_dict["table_columns"] = format_columns_multiline(
            metadata_dict["table_columns"],
            n_per_line=4,
        )

        rows = list(metadata_dict.items())

        print(
            tabulate(
                rows,
                headers=["Field", "Value"],
                tablefmt="fancy_grid",
            )
        )

        return

    def _read_catalog_metadata(
        self,
    ) -> CatalogMetaData:
        """
        Read table-level metadata from the compiled standard-star catalog.

        This internal helper resolves the target table node, opens the HDF5
        catalog, reads the table title, column names, and required table
        attributes, and packages them into a ``CatalogMetaData`` instance.

        Returns
        -------
        CatalogMetaData
            Metadata describing the resolved catalog table.

        Raises
        ------
        ValueError
            Raised if the target table node cannot be uniquely resolved.
        KeyError
            Raised if a required table attribute is missing.
        tables.NoSuchNodeError
            Raised if the resolved node path cannot be opened from the file.
        """

        # get standard star table node path
        node_path: str = resolve_table_path_by_name(
            h5_path=self.catalog_path,
            table_node_name=self.table_node_name,
        )

        # read the basic metadata
        with tb.open_file(self.catalog_path, mode="r") as h5:
            table = h5.get_node(node_path)

            # construct metedata class
            return CatalogMetaData(
                catalog_path=self.catalog_path,
                table_node_path=node_path,
                table_node_name=self.table_node_name,
                table_title=table.title,
                table_columns=tuple(table.colnames),
                bucket_size=table.attrs["bucket_size"],
                ra_unit=table.attrs["ra_unit"],
                dec_unit=table.attrs["dec_unit"],
                mag_system=table.attrs["mag_system"],
                nside=table.attrs["nside"],
                order=table.attrs["order"],
                source=table.attrs["source"],
                version=table.attrs["version"],
            )

    def _get_searching_radius(
        self,
        fits_path: str | Path,
        hdu_index: int = 0,
    ) -> Quantity:
        """
        Compute a coarse search radius that fully encloses the image FOV.

        This method estimates the radius of the circumcircle of the image footprint
        on the sky. It first computes the sky coordinate of the image center, then
        obtains the sky coordinates of the image footprint corners from the WCS, and
        finally returns the maximum angular separation between the center and the
        corners.

        The returned radius is intended for coarse candidate querying only. It may
        include sources that lie outside the true image footprint but inside the
        circumscribed search circle. Such sources should be removed in later geometric
        filtering steps.

        Parameters
        ----------
        fits_path : str | Path
            Path to the FITS image file.
        hdu_index : int, optional
            HDU index containing the image data and WCS header. Default is 0.

        Returns
        -------
        Quantity
            Angular search radius that encloses the full image footprint.

        Raises
        ------
        MissingWCSError
            If the FITS file does not contain usable celestial WCS information.
        KeyError
            If required WCS header keywords are missing.
        ValueError
            If the WCS footprint cannot be computed from the FITS header.
        """

        center_coord: SkyCoord = get_fov_center_skycoord(
            fits_path,
            hdu_index=hdu_index,
        )

        wcs = get_wcs(
            fits_path=fits_path,
            hdu_index=hdu_index,
        )

        # for a square FOV, the footprint is the ra and dec
        # of the image corners
        footprint: npt.NDArray[np.float64] = wcs.calc_footprint()
        corner_coords = SkyCoord(
            ra=footprint[:, 0],
            dec=footprint[:, 1],
            unit=(u.deg, u.deg),
            frame="icrs",
        )

        return max(center_coord.separation(corner_coords))

    def _query_standard_from_catalog(
        self,
        fits_path: str | Path,
        hdu_index: int = 0,
    ) -> dict[str, Any]:
        """
        Query candidate standard stars from the compiled catalog using a coarse
        HEALPix-based spatial search.

        This method performs the catalog-side preselection step for standard-star
        matching. It first computes the image center sky coordinate and a coarse search
        radius that encloses the full image footprint. It then uses a HEALPix cone
        search to identify candidate pixels, converts these pixels into bucket IDs for
        fast PyTables lookup, and finally refines the result in memory by keeping only
        rows whose ``ipix`` values are contained in the queried HEALPix pixel set.

        This method does not yet perform exact image-footprint filtering, edge-margin
        filtering, or filter-availability filtering. Those steps should be handled by
        later geometric and photometric selection functions.

        Parameters
        ----------
        fits_path : str | Path
            Path to the FITS image file whose field of view will be used for the query.
        hdu_index : int, optional
            HDU index containing the image data and WCS header. Default is 0.

        Returns
        -------
        dict[str, Any]
            Dictionary containing query metadata, timing information, and the filtered
            candidate rows from the catalog. The returned dictionary currently includes:

            - ``nside`` : HEALPix NSIDE used for the catalog.
            - ``order`` : HEALPix ordering scheme used for the catalog.
            - ``bucket_size`` : Number of HEALPix pixels grouped into one catalog bucket.
            - ``n_ipix`` : Number of HEALPix pixels returned by the cone search.
            - ``n_bucket`` : Number of unique catalog buckets queried.
            - ``n_candidates_bucket`` : Number of rows returned by the bucket query.
            - ``n_final`` : Number of rows remaining after exact ``ipix`` filtering.
            - ``dt_ipix_sec`` : Time spent in HEALPix cone search.
            - ``dt_bucket_query_sec`` : Time spent reading candidate rows from the table.
            - ``dt_total_sec`` : Total runtime of this method.
            - ``rows`` : Structured NumPy array of candidate catalog rows.
            - ``cond`` : PyTables condition string used for the bucket query.

        Raises
        ------
        MissingWCSError
            If the FITS file does not contain usable celestial WCS information.
        ValueError
            If the catalog table path cannot be resolved or if required metadata are
            invalid.
        tables.NoSuchNodeError
            If the catalog table node cannot be opened from the HDF5 file.

        Notes
        -----
        This method is intentionally a coarse candidate query stage. It may return
        sources that lie outside the true image footprint but inside the enclosing
        search region. Such sources should be removed by later filtering steps such as
        ``filter_stars_inside_fov(...)`` and ``filter_stars_near_edge(...)``.
        """

        # get fov center sky coordinate
        center_coord = get_fov_center_skycoord(
            fits_path=fits_path,
            hdu_index=hdu_index,
        )

        # get searching radius
        radius = self._get_searching_radius(
            fits_path=fits_path,
            hdu_index=hdu_index,
        )

        t0 = time.perf_counter()
        # open catalog file and start searching
        with tb.open_file(self.catalog_path, mode="r") as h5:
            table = h5.get_node(self.table_node_path)

            # compute healpix map
            hp = HEALPix(
                nside=self.metadata.nside,
                order=self.metadata.order,
                frame="icrs",
            )

            # 1) HEALPix coarse selection: ipix list within radius
            t1 = time.perf_counter()
            ipix = hp.cone_search_lonlat(
                center_coord.ra,
                center_coord.dec,
                radius)
            ipix_array = np.asarray(ipix, dtype=np.int64)
            buckets = np.unique(ipix_array // self.metadata.bucket_size).astype(np.int64)
            t2 = time.perf_counter()

            # 2) query by bucket using PyTables in-kernel condition (Numexpr)
            # Build condition: (bucket==b0) | (bucket==b1) | ...
            # NOTE: this is fine because radius is small, so bucket count is small.
            cond = "(" + ") | (".join([f"(bucket == {b})" for b in buckets]) + ")"
            t3 = time.perf_counter()
            rows_bucket = table.read_where(cond)
            t4 = time.perf_counter()

            # 3) get the final filtered rows
            mask_ipix = np.isin(rows_bucket["ipix"], ipix_array)
            rows = rows_bucket[mask_ipix]
            t5 = time.perf_counter()

            out = {
                "nside": self.metadata.nside,
                "order": self.metadata.order,
                "bucket_size": self.metadata.bucket_size,
                "n_ipix": int(len(ipix_array)),
                "n_bucket": int(len(buckets)),
                "n_candidates_bucket": int(len(rows_bucket)),
                "n_final": int(len(rows)),
                "dt_ipix_sec": t2 - t1,
                "dt_bucket_query_sec": t4 - t3,
                "dt_total_sec": (t5 - t0),
                "rows": rows,  # numpy structured array
                "cond": cond,
            }

            return out

    def _filter_star_in_fov(
        self,
        fits_path: str | Path,
        hdu_index: int = 0
    ) -> pd.DataFrame:
        """
        Filter queried candidate standard stars to those contained within the image FOV.

        This method takes the coarse candidate stars returned by
        ``_query_standard_from_catalog(...)``, converts the queried catalog rows into a
        pandas ``DataFrame``, constructs ``SkyCoord`` objects from the candidate
        RA/Dec columns, and keeps only those stars whose sky positions are contained
        within the FITS image footprint defined by the image WCS.

        This step performs image-footprint filtering only. It does not yet exclude
        stars near the image edge or stars missing magnitudes in requested filters.

        Parameters
        ----------
        fits_path : str | Path
            Path to the FITS image file whose field of view defines the valid sky
            region for standard-star selection.
        hdu_index : int, optional
            HDU index containing the image data and WCS header. Default is 0.

        Returns
        -------
        pd.DataFrame
            DataFrame containing only the queried candidate stars that lie inside the
            image field of view.

        Raises
        ------
        MissingWCSError
            If the FITS file does not contain usable celestial WCS information.
        KeyError
            If required catalog columns such as ``ra`` or ``dec`` are missing from the
            queried candidate table.
        ValueError
            If WCS-based footprint filtering cannot be completed.

        Notes
        -----
        This method relies on ``SkyCoord.contained_by(...)`` to test whether candidate
        sky coordinates fall inside the image footprint. Depending on the available WCS
        shape information, it may later be necessary to pass the image explicitly to the
        containment check.
        """

        # query stars from the catalog
        queried_stars = self._query_standard_from_catalog(
            fits_path=fits_path,
            hdu_index=hdu_index,
        )

        # convert the queried stars to a DataFrame for further operations
        df_star = pd.DataFrame.from_records(queried_stars["rows"])

        star_coords = SkyCoord(
            ra=df_star.loc[:, "ra"].to_numpy(),
            dec=df_star.loc[:, "dec"].to_numpy(),
            unit=(u.deg, u.deg),
            frame="icrs",
        )

        wcs = get_wcs(fits_path=fits_path, hdu_index=hdu_index)

        # mask out stars not in the fov
        mask: npt.NDArray[np.bool_] = star_coords.contained_by(wcs)
        df_star_in_fov = df_star.loc[mask].copy()

        return df_star_in_fov

    def _filter_stars_near_edge(
        self,
        df_star_in_fov: pd.DataFrame,
        fits_path: str | Path,
        hdu_index: int = 0,
        edge_margin_pix: float = 20.0,
        add_pixel_columns: bool = True,
    ) -> pd.DataFrame:
        """
        Remove stars that are too close to the image edge in pixel space.

        This method assumes the input stars have already been filtered to lie inside
        the image FOV. It converts the stars' sky coordinates to pixel coordinates
        using the FITS WCS, computes each star's minimum distance to the four image
        edges, and keeps only stars whose minimum edge distance is larger than the
        requested pixel margin.

        Parameters
        ----------
        df_star_in_fov : pd.DataFrame
            DataFrame of candidate stars already confirmed to be inside the image
            footprint. It must contain ``ra`` and ``dec`` columns in degrees.
        fits_path : str | Path
            Path to the FITS image file.
        hdu_index : int, optional
            HDU index containing the image data and WCS header. Default is 0.
        edge_margin_pix : float, optional
            Minimum allowed distance to the nearest image edge, in pixels.
            Default is 20.0.
        add_pixel_columns : bool, optional
            If True, add ``x_pix``, ``y_pix``, and ``edge_dist_pix`` columns to the
            returned DataFrame. Default is True.

        Returns
        -------
        pd.DataFrame
            DataFrame containing only stars farther than ``edge_margin_pix`` from
            the nearest image edge.

        Raises
        ------
        KeyError
            If ``ra`` or ``dec`` columns are missing.
        MissingWCSError
            If the FITS file does not contain usable celestial WCS information.
        ValueError
            If the FITS image is not two-dimensional.
        """
        # The edge cut is defined in pixel space, so the input table must at least
        # contain sky coordinates that can be projected onto the detector.
        if "ra" not in df_star_in_fov.columns or "dec" not in df_star_in_fov.columns:
            raise KeyError("Input DataFrame must contain 'ra' and 'dec' columns.")

        # A negative edge margin is physically meaningless.
        if edge_margin_pix < 0:
            raise ValueError("edge_margin_pix must be non-negative.")

        # Load the WCS used to convert sky coordinates to detector pixel coordinates.
        wcs: WCS = get_wcs(fits_path=fits_path, hdu_index=hdu_index)

        # Read the image so we can determine its 2D shape.
        with fits.open(fits_path, memmap=False) as hdul:
            data = hdul[hdu_index].data

        # This method assumes a normal 2D imaging HDU.
        if data is None or data.ndim != 2:
            raise ValueError("Expected a 2D FITS image in the requested HDU.")

        # NumPy image shape is (ny, nx) = (rows, columns).
        ny: int
        nx: int
        ny, nx = data.shape

        # Build sky coordinates for the stars that already passed the FOV filter.
        star_coords = SkyCoord(
            ra=df_star_in_fov.loc[:, "ra"].to_numpy(),
            dec=df_star_in_fov.loc[:, "dec"].to_numpy(),
            unit=(u.deg, u.deg),
            frame="icrs",
        )

        # Convert world coordinates to 0-based pixel coordinates.
        # x increases left -> right, y increases bottom -> top in the WCS transform.
        x_pix, y_pix = skycoord_to_pixel(star_coords, wcs, origin=0, mode="all")

        # Compute each star's distance to the four detector edges in pixels.
        # For a valid in-FOV star:
        #   left   edge distance = x
        #   right  edge distance = (nx - 1) - x
        #   bottom edge distance = y
        #   top    edge distance = (ny - 1) - y
        d_left = x_pix
        d_right = (nx - 1) - x_pix
        d_bottom = y_pix
        d_top = (ny - 1) - y_pix

        # The nearest edge controls whether the star is safe for photometry.
        edge_dist_pix = np.minimum.reduce([d_left, d_right, d_bottom, d_top])

        # Keep only stars whose nearest-edge distance is larger than the requested margin.
        mask = edge_dist_pix > edge_margin_pix

        # Copy the filtered subset so later column edits are explicit and safe.
        df_safe = df_star_in_fov.loc[mask].copy()

        # Optionally keep the derived pixel-space quantities for debugging,
        # plotting, and downstream visualization.
        if add_pixel_columns:
            df_safe.loc[:, "x_pix"] = x_pix[mask]
            df_safe.loc[:, "y_pix"] = y_pix[mask]
            df_safe.loc[:, "edge_dist_pix"] = edge_dist_pix[mask]

        return df_safe

    def _filter_stars_by_mag_limit(
        self,
        df_star: pd.DataFrame,
        mag_limit: float,
        mag_columns: list[str],
    ) -> pd.DataFrame:
        """
        Filter stars by an upper magnitude limit across a set of catalog columns.

        This internal helper keeps only stars whose magnitudes are not fainter than the
        requested limit in all specified magnitude columns. In the astronomical
        magnitude system, larger magnitude values correspond to fainter sources, so
        this filter applies an upper limit using ``<= mag_limit``.

        Parameters
        ----------
        df_star : pd.DataFrame
            Input star table.
        mag_limit : float
            Maximum allowed magnitude value. Stars fainter than this limit in any of
            the selected columns are removed.
        mag_columns : list[str]
            List of catalog magnitude column names to check.

        Returns
        -------
        pd.DataFrame
            Filtered DataFrame containing only stars that satisfy the magnitude limit
            in all requested magnitude columns.

        Raises
        ------
        KeyError
            If any requested magnitude column is missing from the input DataFrame.

        Notes
        -----
        This method currently requires explicit catalog magnitude column names. In a
        future update, the magnitude-column selection should be derived from a catalog-
        level mapping between standard filter names and catalog magnitude columns. That
        mapping is intended to live in the compiled catalog metadata rather than in a
        telescope-specific filter class.
        """

        df_filtered = df_star.loc[(df_star[mag_columns] <= mag_limit).all(axis=1)].copy()

        return df_filtered

    def find_standard_stars(
        self,
        fits_path: str | Path,
        hdu_index: int = 0,
        edge_margin_pix: float = 20.0,
        add_pixel_columns: bool = True,
        mag_limit: float | None = None,
        mag_columns: list[str] | None = None,
    ) -> pd.DataFrame:
        """
        Find usable standard stars for a given FITS image.

        This public method runs the standard-star selection pipeline for a single FITS
        image. It first queries candidate stars from the compiled catalog and filters
        them to those contained within the image field of view. It then removes stars
        that are too close to the image edge in pixel space. An optional magnitude
        limit can also be applied across a specified set of catalog magnitude columns.

        The returned table is intended to contain only stars that are geometrically
        usable for later photometric calibration. Optional pixel-space columns can be
        added for debugging, plotting, and downstream visualization.

        Parameters
        ----------
        fits_path : str | Path
            Path to the FITS image file.
        hdu_index : int, optional
            HDU index containing the image data and WCS header. Default is 0.
        edge_margin_pix : float, optional
            Minimum allowed distance from the nearest image edge, in pixels.
            Default is 20.0.
        add_pixel_columns : bool, optional
            If True, add pixel-coordinate and edge-distance columns to the returned
            DataFrame. Default is True.
        mag_limit : float | None, optional
            Optional upper magnitude limit. If provided, stars are kept only if all
            selected magnitude columns satisfy ``magnitude <= mag_limit``. If None, no
            magnitude filtering is applied.
        mag_columns : list[str] | None, optional
            Catalog magnitude column names to use for the magnitude-limit filter.
            This must be provided when ``mag_limit`` is not None. If ``mag_limit`` is
            None, this parameter is ignored.

        Returns
        -------
        pd.DataFrame
            DataFrame containing the final set of standard stars that lie inside the
            image FOV, pass the edge-distance cut, and optionally satisfy the magnitude
            limit.

        Raises
        ------
        MissingWCSError
            If the FITS file does not contain usable celestial WCS information.
        KeyError
            If required catalog columns such as ``ra`` or ``dec`` are missing, or if a
            requested magnitude column is not present.
        ValueError
            If geometric filtering cannot be completed because of invalid image or WCS
            information, or if ``mag_limit`` is provided without ``mag_columns``.

        Notes
        -----
        This method is the public API for standard-star selection. The intermediate
        results are stored on the instance as ``self.stars_in_fov`` and
        ``self.final_stars`` for convenience in debugging and comparison.
        ``matched_stars`` exposes a snapshot after all cuts, including the
        optional magnitude limit; final_stars retains its legacy pre-limit
        meaning. The returned DataFrame and its byte-string IDs are unchanged.

        The current magnitude-filter interface requires explicit catalog magnitude
        column names. In a future update, this should be replaced by a catalog-driven
        mapping from standard filter names to catalog magnitude columns. That mapping
        should be stored in the compiled catalog metadata and loaded into
        ``self.metadata`` when the catalog is read, rather than being delegated to a
        telescope-specific filter class.
        """

        # query the catalog and find standard stars
        # inside the image FOV
        self.stars_in_fov = self._filter_star_in_fov(
            fits_path=fits_path,
            hdu_index=hdu_index,
        )

        self.final_stars = self._filter_stars_near_edge(
            df_star_in_fov=self.stars_in_fov,
            fits_path=fits_path,
            hdu_index=hdu_index,
            edge_margin_pix=edge_margin_pix,
            add_pixel_columns=add_pixel_columns,
        )

        if mag_limit is None:
            selected_stars = self.final_stars
        else:
            if mag_columns is None:
                raise ValueError("You must provide the magnitude names.")
            else:
                selected_stars = self._filter_stars_by_mag_limit(
                    df_star=self.final_stars,
                    mag_limit=mag_limit,
                    mag_columns=mag_columns,
                )

        self._matched_stars = selected_stars.copy(deep=True)
        return selected_stars

    def plot_standard_stars_on_fits(
        self,
        fits_path: str | Path,
        df_stars: pd.DataFrame | None = None,
        hdu_index: int = 0,
        marker_size: float = 80.0,
        marker_edgecolor: str = "cyan",
        marker_facecolor: str = "none",
        marker_linewidth: float = 1.5,
        cmap: str = "gray",
        stretch: str = "sqrt",
        percent: float = 99.0,
        figsize: tuple[float, float] = (8, 8),
        title: str | None = None,
    ):
        """
        Plot selected standard stars on a FITS image in world coordinates.

        The FITS image is displayed with WCSAxes so that the x and y axes are shown
        in celestial coordinates (typically RA/Dec). The selected stars are marked
        as small circles, and RA/Dec grid lines are overlaid.

        Parameters
        ----------
        fits_path : str | Path
            Path to the FITS image file.
        df_stars : pd.DataFrame | None, optional
            DataFrame containing standard stars to plot. If None, this method will
            try to use ``self.final_stars``.
            The table must contain ``ra`` and ``dec`` columns in degrees.
        hdu_index : int, optional
            HDU index containing the image data and WCS header. Default is 0.
        marker_size : float, optional
            Marker size passed to ``scatter_coord``. Default is 80.0.
        marker_edgecolor : str, optional
            Marker edge color. Default is ``"cyan"``.
        marker_facecolor : str, optional
            Marker face color. Default is ``"none"``.
        marker_linewidth : float, optional
            Marker line width. Default is 1.5.
        cmap : str, optional
            Colormap for the FITS image. Default is ``"gray"``.
        stretch : str, optional
            Stretch passed to ``simple_norm``. Default is ``"sqrt"``.
        percent : float, optional
            Percentile passed to ``simple_norm``. Default is 99.0.
        figsize : tuple[float, float], optional
            Figure size in inches. Default is (8, 8).
        title : str | None, optional
            Plot title. If None, a default title is used.

        Returns
        -------
        tuple[matplotlib.figure.Figure, astropy.visualization.wcsaxes.WCSAxes]
            The created figure and axes.
        """
        if df_stars is None:
            if not hasattr(self, "final_stars"):
                raise ValueError(
                    "No star table was provided and self.final_stars does not exist."
                )
            df_stars = self.final_stars

        if "ra" not in df_stars.columns or "dec" not in df_stars.columns:
            raise KeyError("df_stars must contain 'ra' and 'dec' columns.")

        # Read image data and WCS.
        with fits.open(fits_path, memmap=False) as hdul:
            data = hdul[hdu_index].data
            header = hdul[hdu_index].header

        if data is None or data.ndim != 2:
            raise ValueError("Expected a 2D FITS image in the requested HDU.")

        wcs = WCS(header).celestial

        # Build star coordinates for overplotting.
        star_coords = SkyCoord(
            ra=df_stars["ra"].to_numpy(),
            dec=df_stars["dec"].to_numpy(),
            unit=(u.deg, u.deg),
            frame="icrs",
        )

        # Create WCSAxes so the axes are in RA/Dec.
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection=wcs)

        # Display image.
        norm = simple_norm(data, stretch=stretch, percent=percent)
        ax.imshow(data, origin="lower", cmap=cmap, norm=norm)

        # Plot RA/Dec grid.
        ax.grid(color="white", alpha=0.5, linestyle=":")

        # Label axes in world coordinates.
        ax.set_xlabel("RA")
        ax.set_ylabel("Dec")

        # Prevent autoscale from changing when overlays are added.
        ax.set_autoscale_on(False)

        # Plot stars as small circles in sky coordinates.
        ax.scatter_coord(
            star_coords,
            s=marker_size,
            edgecolor=marker_edgecolor,
            facecolor=marker_facecolor,
            linewidth=marker_linewidth,
        )

        if title is None:
            title = "Standard stars on FITS image"
        ax.set_title(title)

        return fig, ax


# module level functions
def resolve_table_path_by_name(
    h5_path: str | Path,
    table_node_name: str
) -> str:
    """
    Resolve the unique PyTables table path for a given table node name.

    The function walks the HDF5 node tree, searches all ``Table`` nodes, and
    returns the full node path of the unique table whose node name matches the
    requested value.

    Parameters
    ----------
    h5_path : str | Path
        Path to the HDF5 catalog file.
    table_node_name : str
        Name of the target PyTables table node.

    Returns
    -------
    str
        Full PyTables node path of the matched table.

    Raises
    ------
    ValueError
        Raised if no table with the requested node name exists, or if more than
        one matching table is found.
    """

    with tb.open_file(h5_path, mode="r") as h5:
        node_paths = [
            node._v_pathname
            for node in h5.walk_nodes("/", classname="Table")
            if node._v_name == table_node_name
        ]

        if len(node_paths) == 0:
            raise ValueError(f"Catalog has no table named '{table_node_name}'.")
        elif len(node_paths) > 1:
            raise ValueError(f"Catalog has more than one table named '{table_node_name}'.")
        else:
            return node_paths[0]


def format_columns_multiline(columns: tuple[str, ...], n_per_line: int = 4) -> str:
    lines = []
    for i in range(0, len(columns), n_per_line):
        chunk = columns[i:i + n_per_line]
        lines.append(", ".join(chunk))
    return "\n".join(lines)
