"""Connect reviewed region files to the original optical photometry equations.

Numerical error and magnitude calculations remain in Photometry. This module
only handles source/image association, masks, per-source files and result slots.
"""

from __future__ import annotations

from typing import TYPE_CHECKING
from collections.abc import Sequence
from pathlib import Path
import warnings

import astropy.units as u
from astropy.io import fits
from astropy.table import QTable
import numpy as np
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
from regions import CirclePixelRegion, CircleAnnulusPixelRegion, Regions

from ..calibration._calibration_io import children, paths
from ..sources import FilterMagnitude, Sources
from ..sources.filter_magnitude import UPPER_LIMIT_ERROR
from .regions import PhotometryRegions, PhotometryRegionCollection

if TYPE_CHECKING:
    from .photometry import Photometry


def _pairs(session: PhotometryRegions) -> dict[str, dict]:
    """Read actual final file geometry, preserving IDs and independent radii."""
    pairs: dict[str, dict] = {}
    known = {row["source_id"] for row in session.records if row["accepted"]}
    for path, kind, expected in (
        (session.aperture_path, "aperture", CirclePixelRegion),
        (session.background_path, "background", CircleAnnulusPixelRegion),
    ):
        with warnings.catch_warnings(record=True) as messages:
            warnings.simplefilter("always")
            parsed = Regions.read(path, format="ds9")
        if messages:
            raise ValueError(f"Region file could not be parsed completely: {path}.")
        for region in parsed:
            tags = region.meta.get("tag", [])
            ids = [tag[3:] for tag in tags if tag.startswith("id:")]
            if len(ids) != 1 or ids[0] not in known or kind not in tags:
                raise ValueError(f"Unexpected region ID/role in {path}.")
            if not isinstance(region, expected) or not region.meta.get("include", True):
                raise ValueError(f"Expected included circular pixel regions in {path}.")
            pair = pairs.setdefault(ids[0], {})
            if kind in pair:
                raise ValueError(f"Duplicate {kind} region for {ids[0]!r}.")
            pair[kind] = region
    if set(pairs) != known or any(
        set(pair) != {"aperture", "background"} for pair in pairs.values()
    ):
        raise ValueError(
            "Saved files do not match the accepted source IDs and paired region records."
        )
    return pairs


def run_aperture_photometry(
    measurer: Photometry,
    sources: Sources | Sequence[Sources],
    regions: (
        PhotometryRegions
        | PhotometryRegionCollection
        | Sequence[PhotometryRegions | PhotometryRegionCollection]
    ),
    *,
    bkg_clip_sigma: float,
    src_detection_sigma: float,
    hdu: int,
    require_review: bool,
) -> QTable:
    """Run the old measurement sequence with the new reviewed-file associations.

    Source-slot updates are deferred until every input has been processed.
    Rejected/unusable filters stay NaN in an otherwise measured source; a source
    with no usable measurements has an undefined instrumental slot. Duplicate
    source/filter measurements are rejected rather than silently overwritten.
    """
    if not isinstance(sources, Sources):
        if not isinstance(sources, Sequence) or isinstance(sources, (str, bytes)):
            raise TypeError("sources must be Sources or a sequence of Sources.")
        sources = Sources.combine(*sources)
    if not len(sources):
        raise ValueError("A non-empty Sources collection is required.")
    if hdu != 0:
        raise ValueError("Reviewed photometry requires primary image hdu=0.")
    for name, value in (
        ("bkg_clip_sigma", bkg_clip_sigma),
        ("src_detection_sigma", src_detection_sigma),
    ):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be positive and finite.")
    groups = (
        (regions,)
        if isinstance(regions, (PhotometryRegions, PhotometryRegionCollection))
        else regions
    )
    if not isinstance(groups, Sequence) or isinstance(groups, (str, bytes)):
        raise TypeError(
            "regions must be region sessions, collections, or a sequence of them."
        )
    sessions: list[PhotometryRegions] = []
    for group in groups:
        if isinstance(group, PhotometryRegions):
            sessions.append(group)
        elif isinstance(group, PhotometryRegionCollection):
            sessions.extend(group.sessions)
        else:
            raise TypeError(
                "Each regions entry must be a region session or collection."
            )
    if not sessions:
        raise ValueError("No region sessions were supplied.")
    by_image: dict[Path, list[PhotometryRegions]] = {}
    for session in sessions:
        by_image.setdefault(session.image_path, []).append(session)
    sessions = [
        session for image_sessions in by_image.values() for session in image_sessions
    ]
    selected_images = {
        path for child in children(measurer._image_collection) for path in paths(child)
    }
    measurements: list[dict] = []
    seen: set[tuple[str, str, str]] = set()
    attempted: dict[tuple[str, str], list[dict]] = {}
    # Reuse data/header/mask for the same frame even when role-specific region
    # files are separate. Release it when changing images to bound memory use.
    last_image: Path | None = None

    for supplied in sessions:
        # Refresh only this known session's records, not the image directory.
        session = PhotometryRegions.load(supplied.root, supplied.manifest_path)
        if session.image_path not in selected_images:
            raise ValueError(
                f"Region image is not in the supplied image collection: {session.image_path}"
            )
        records = session.records
        if require_review and any(not row["reviewed"] for row in records):
            raise ValueError(
                f"Complete region review before photometry: {session.image_path.name}"
            )
        pairs = _pairs(session)
        if session.image_path != last_image:
            with fits.open(session.image_path, mode="readonly", memmap=False) as hdul:
                data = np.asarray(hdul[0].data, dtype=float)
                header = hdul[0].header.copy()
                if data.ndim != 2:
                    raise ValueError("Expected a primary 2D science image.")
                invalid = ~np.isfinite(data)
                if "MASK" in hdul:
                    flags = hdul["MASK"].data
                    if flags is None or flags.shape != data.shape:
                        raise ValueError("MASK shape must match the primary image.")
                    invalid |= np.asarray(flags, dtype=bool)
            last_image = session.image_path
        if u.Unit(header["BUNIT"]) != u.adu:
            raise u.UnitsError(
                "The original photometry error model requires image counts in ADU."
            )
        gain_value = float(header["GAIN"])
        read_noise = float(header["RDNOISE"])
        n_com = float(header["NCOMBINE"])
        if (
            not np.all(np.isfinite([gain_value, read_noise, n_com]))
            or gain_value <= 0
            or read_noise < 0
            or n_com < 1
            or not n_com.is_integer()
        ):
            raise ValueError(
                "Require finite GAIN>0, RDNOISE>=0 and positive integer NCOMBINE."
            )
        gain = gain_value * u.electron / u.adu
        readout_noise = read_noise * u.electron

        for record in records:
            info = sources.get_info(
                record["source_id"], field_name=record["field_name"]
            )
            if (
                header.get("FIELD") != info.observation_field.field_name
                or header.get("IMTYPE") != "Master Light"
            ):
                raise ValueError("Region source must match a Master Light FIELD.")
            filter_name = info.observation_field.telescope.resolve_filter_alias(
                str(header["FILTER"])
            )
            if session._state.get("filter_name", filter_name) != filter_name:
                raise ValueError(
                    "Region record FILTER does not match the image header."
                )
            identity = (record["field_name"], record["source_id"])
            key = (*identity, filter_name)
            if key in seen:
                raise ValueError(f"Duplicate source/filter measurement: {key!r}.")
            seen.add(key)
            row = {
                "id": len(measurements) + 1,
                "source_id": info.source_name,
                "source_role": info.source_role,
                "field": identity[0],
                "filter": filter_name,
                "image_path": session.image_path.relative_to(session.root).as_posix(),
                "aperture_path": session.aperture_path.relative_to(
                    session.root
                ).as_posix(),
                "background_path": session.background_path.relative_to(
                    session.root
                ).as_posix(),
                **{
                    name: np.nan
                    for name in (
                        "xcenter",
                        "ycenter",
                        "background_x",
                        "background_y",
                        "aperture_radius",
                        "inner_radius",
                        "outer_radius",
                        "n_ap",
                        "n_ann",
                        "mean_bkg",
                        "src+bkg",
                        "bkg",
                        "src",
                        "error_counts",
                        "error_counts_ul",
                        "mag_inst",
                        "mag_inst_error",
                        "src_significance",
                        "aperture_valid_fraction",
                        "annulus_valid_fraction",
                    )
                },
                "is_upper_limit": False,
                "status": "rejected",
                "reason": "Rejected during region review.",
            }
            measurements.append(row)
            attempted.setdefault(identity, []).append(row)
            if not record["accepted"]:
                continue
            pair = pairs[record["source_id"]]
            src_region, bkg_region = pair["aperture"], pair["background"]
            aperture = CircularAperture(
                (src_region.center.x, src_region.center.y), src_region.radius
            )
            annulus = CircularAnnulus(
                (bkg_region.center.x, bkg_region.center.y),
                bkg_region.inner_radius,
                bkg_region.outer_radius,
            )
            row.update(
                xcenter=float(src_region.center.x),
                ycenter=float(src_region.center.y),
                background_x=float(bkg_region.center.x),
                background_y=float(bkg_region.center.y),
                aperture_radius=float(src_region.radius),
                inner_radius=float(bkg_region.inner_radius),
                outer_radius=float(bkg_region.outer_radius),
                aperture_valid_fraction=PhotometryRegions._valid_fraction(
                    aperture, invalid
                ),
                annulus_valid_fraction=PhotometryRegions._valid_fraction(
                    annulus, invalid
                ),
            )
            if row["aperture_valid_fraction"] < 1 - 1e-8:
                row.update(
                    status="invalid_aperture",
                    reason="Source aperture contains invalid pixels or crosses the image edge; no missing-flux correction applied.",
                )
                continue

            # Retain the original sigma-clipped mean, geometric source area,
            # clipped annulus sum_aper_area, exact sum and background subtraction.
            mean_bkg, bkg_stats = measurer.estimate_mean_background(
                data,
                annulus,
                clip_sigma=bkg_clip_sigma,
                mask=invalid,
            )
            n_ap = aperture.area * u.pix
            n_ann_eff = bkg_stats.sum_aper_area / u.pix
            if (
                not np.all(np.isfinite(mean_bkg.value))
                or not np.all(np.isfinite(n_ann_eff.value))
                or np.any(n_ann_eff.value <= 0)
            ):
                row.update(
                    status="invalid_background",
                    reason="No finite background mean or positive effective annulus area.",
                )
                continue
            total_bkg = mean_bkg * aperture.area * u.pix
            phot_table = aperture_photometry(
                data, aperture, method="exact", mask=invalid
            )
            phot_table["aperture_sum"].name = "src+bkg"
            phot_table["src+bkg"].unit = u.adu
            phot_src = phot_table["src+bkg"] - total_bkg

            # Call the user's original error expressions without reimplementing
            # them, and retain both measured and background-only count errors.
            error_adu, error_ul_adu = measurer.estimate_source_region_error(
                gain=gain,
                readout_noise=readout_noise,
                n_ap=n_ap,
                n_ann=n_ann_eff,
                aperture_counts=phot_table["src+bkg"],
                mean_bkg=mean_bkg,
                n_com=n_com,
            )
            row.update(
                **{
                    "n_ap": float(n_ap.value),
                    "n_ann": float(np.atleast_1d(n_ann_eff.value)[0]),
                    "mean_bkg": float(mean_bkg.value[0]),
                    "src+bkg": float(phot_table["src+bkg"].value[0]),
                    "bkg": float(total_bkg.value[0]),
                    "src": float(phot_src.value[0]),
                    "error_counts": float(error_adu.value[0]),
                    "error_counts_ul": float(error_ul_adu.value[0]),
                }
            )
            if (
                not np.all(np.isfinite([row["error_counts"], row["error_counts_ul"]]))
                or min(row["error_counts"], row["error_counts_ul"]) <= 0
            ):
                row.update(
                    status="invalid_noise",
                    reason="Original noise model returned non-finite or non-positive count uncertainty.",
                )
                continue
            m_inst, m_inst_error, significance = measurer.counts2mag(
                src_counts=phot_src,
                error_counts=error_adu,
                error_counts_ul=error_ul_adu,
                detection_significance=src_detection_sigma,
            )
            upper_limit = bool(
                significance[0] < src_detection_sigma or phot_src.value[0] <= 0
            )
            row.update(
                mag_inst=float(m_inst.value[0]),
                mag_inst_error=float(m_inst_error.value[0]),
                src_significance=float(significance[0]),
                is_upper_limit=upper_limit,
                status="upper_limit" if upper_limit else "measured",
                reason="",
            )

    # Prepare all results before mutating any source's shared measurement slots.
    pending: list[tuple[tuple[str, str], FilterMagnitude | None]] = []
    for identity, rows in attempted.items():
        usable = [row for row in rows if row["status"] in ("measured", "upper_limit")]
        info = sources.get_info(identity[1], field_name=identity[0])
        result = (
            None
            if not usable
            else FilterMagnitude(
                telescope=info.observation_field.telescope,
                filter_names=[row["filter"] for row in usable],
                magnitudes=[row["mag_inst"] for row in usable],
                errors=[
                    (
                        UPPER_LIMIT_ERROR.to_value(u.mag)
                        if row["is_upper_limit"]
                        else row["mag_inst_error"]
                    )
                    for row in usable
                ],
                significance=[row["src_significance"] for row in usable],
                magnitude_system="instrumental",
            )
        )
        pending.append((identity, result))
    table = QTable(rows=measurements)
    for name in ("src+bkg", "bkg", "src", "error_counts", "error_counts_ul"):
        table[name].unit = u.adu
    for name in ("mag_inst", "mag_inst_error"):
        table[name].unit = u.mag
    for name in ("n_ap", "n_ann"):
        table[name].unit = u.pix
    table["mean_bkg"].unit = u.adu / u.pix
    table.meta.update(
        magnitude_system="instrumental",
        count_convention="average-stack ADU counts; no exposure-time normalization",
        detection_significance=src_detection_sigma,
    )
    for identity, result in pending:
        sources.set_magnitudes(
            kind="instrumental",
            magnitudes=result,
            source_name=identity[1],
            field_name=identity[0],
        )
    return table
