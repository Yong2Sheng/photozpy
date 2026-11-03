"""Generate and persist per-image apertures using the legacy centroid workflow.

Retains CCD_Regions' background subtraction, repeated center-of-mass estimates,
and final quadratic centroid. Legacy Swift helpers remain unchanged.
"""

from __future__ import annotations

from copy import deepcopy
import json
import os
from pathlib import Path
import re
import tempfile
import warnings

from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
import numpy as np
from photutils.aperture import CircularAperture, CircularAnnulus
from photutils.background import Background2D, MedianBackground
from photutils.centroids import centroid_com, centroid_quadratic, centroid_sources
from regions import (
    CirclePixelRegion,
    CircleAnnulusPixelRegion,
    PointPixelRegion,
    Regions,
)

from ..sources import Sources
from ..calibration._calibration_io import ImageInput, children, select


class PhotometryRegions:
    """One image's shared aperture/background files and source review records.

    Paths in the manifest are relative to ``root`` (normally analysis_dir).
    Source names, including standard catalog id_name, identify paired regions.
    Coordinates in records are zero-based pixels; DS9 image coordinates are
    converted explicitly at the serialization boundary. SourceInfo catalog
    coordinates and all magnitude slots are left untouched.

    Use generate for one image/group, or generate_regions to process all
    supplied sources and filters automatically. Regeneration requires
    overwrite=True and discards saved manual edits; calling review again does
    not regenerate files. Internal manifests support the isolated GUI worker.
    """

    def __init__(self, root: Path, manifest: Path, state: dict) -> None:
        self.root = root.expanduser().resolve()
        self.manifest_path = manifest
        self._state = state

    def _resolve(self, relative: str | Path) -> Path:
        """Resolve a relative session path without escaping the analysis root."""
        path = Path(relative)
        resolved = (self.root / path).resolve()
        if path.is_absolute() or not resolved.is_relative_to(self.root):
            raise ValueError(
                f"Session path must be relative to the analysis root: {path}"
            )
        return resolved

    @property
    def image_path(self) -> Path:
        """Return the current absolute image path, derived from the root."""
        return self._resolve(self._state["image_path"])

    @property
    def records(self) -> list[dict]:
        """Return a copy of source IDs, relative paths, positions and decisions."""
        return deepcopy(self._state["sources"])

    @property
    def aperture_path(self) -> Path:
        """Return the current absolute aperture-file path."""
        return self._resolve(self._state["aperture_path"])

    @property
    def background_path(self) -> Path:
        """Return the current absolute background-file path."""
        return self._resolve(self._state["background_path"])

    @classmethod
    def load(cls, root: str | Path, manifest_path: str | Path) -> PhotometryRegions:
        """Restore a session from a root-relative JSON manifest, without writing."""
        session = cls(Path(root), Path(manifest_path), {})
        state = json.loads(session._resolve(manifest_path).read_text(encoding="utf-8"))
        if state.get("version") != 1:
            raise ValueError("Unsupported photometry-region manifest version.")
        session._state = state
        for key in ("image_path", "aperture_path", "background_path"):
            path = session._resolve(state[key])
            if not path.is_file():
                raise FileNotFoundError(path)
        return session

    @staticmethod
    def _valid_fraction(aperture, invalid: np.ndarray) -> float:
        """Measure valid geometric coverage, without inferring stack depth."""
        mask = aperture.to_mask(method="exact")
        large, small = mask.get_overlap_slices(invalid.shape)
        if large is None:
            return 0.0
        return float(np.sum(mask.data[small] * ~invalid[large]) / aperture.area)

    @staticmethod
    def _filename_part(value: str) -> str:
        """Keep readable names while excluding path separators and traversal."""
        part = re.sub(r"[^A-Za-z0-9_.+'~-]+", "_", value).strip(".")
        if not part:
            raise ValueError(f"Cannot construct a region filename from {value!r}.")
        return part

    @classmethod
    def _output_layout(
        cls, header: fits.Header, sources: Sources
    ) -> tuple[Path, str, str]:
        """Group standards by field, science sources by name; files mark filters."""
        infos = sources.source_infos
        if not infos:
            raise ValueError("At least one source is required.")
        field = infos[0].observation_field
        if any(
            info.observation_field.field_name != header.get("FIELD") for info in infos
        ):
            raise ValueError("Source field does not match image FIELD.")
        filter_name = field.telescope.resolve_filter_alias(str(header["FILTER"]))
        field_part = cls._filename_part(field.field_name)
        filter_part = cls._filename_part(filter_name)
        if all(info.source_role == "standard" for info in infos):
            group, stem = f"{field_part}_standards", filter_part
        elif len(infos) == 1 and infos[0].source_role == "science":
            group = cls._filename_part(infos[0].source_name)
            stem = f"{field_part}_{filter_part}"
        else:
            raise ValueError(
                "Use generate_regions to split science targets into individual files."
            )
        return Path("photometry/regions") / group, stem, filter_name

    @classmethod
    def _plan_regions(
        cls,
        image_collection: ImageInput,
        sources: Sources,
        root: str | Path,
        *,
        overwrite: bool = False,
        load_existing: bool = False,
    ) -> list[tuple[Path, Sources]]:
        """Resolve image/source groups and validate all expected paths before I/O."""
        if not isinstance(sources, Sources) or not len(sources):
            raise ValueError("A non-empty Sources collection is required.")
        root_path = Path(root).expanduser().resolve()
        collections = children(image_collection)
        if any(child.ext != 0 for child in collections):
            raise ValueError("Region generation requires collection ext=0.")
        fields = dict.fromkeys(
            info.observation_field.field_name for info in sources.source_infos
        )
        plans: list[tuple[Path, Sources]] = []
        destinations: set[Path] = set()
        for field_name in fields:
            in_field = sources.for_field(field_name)
            telescope = in_field.source_infos[0].observation_field.telescope
            paths = {
                path
                for child in collections
                for path in select(child, IMTYPE="Master Light", FIELD=field_name)
            }
            if not paths:
                raise ValueError(f"No Master Light images for FIELD={field_name!r}.")
            headers = {path: fits.getheader(path) for path in paths}
            order = {name: index for index, name in enumerate(telescope.filters())}
            paths = sorted(
                paths,
                key=lambda path: (
                    order[telescope.resolve_filter_alias(str(headers[path]["FILTER"]))],
                    str(path),
                ),
            )
            standards = in_field.by_role("standard")
            groups = ([standards] if len(standards) else []) + list(
                in_field.by_role("science")
            )
            for group in groups:
                for path in paths:
                    path.relative_to(root_path)
                    directory, stem, _ = cls._output_layout(headers[path], group)
                    for suffix in ("aperture.reg", "background.reg", "regions.json"):
                        destination = (
                            root_path / directory / f"{stem}_{suffix}"
                        ).resolve()
                        if not destination.is_relative_to(root_path):
                            raise ValueError("Output path escapes analysis root.")
                        if destination in destinations:
                            raise ValueError(
                                f"Multiple images or source groups map to {destination}."
                            )
                        if load_existing and not destination.is_file():
                            raise FileNotFoundError(
                                f"Missing saved region file: {destination}; nothing was regenerated."
                            )
                        if not load_existing and destination.exists() and not overwrite:
                            raise FileExistsError(destination)
                        destinations.add(destination)
                    plans.append((path, group))
        return plans

    @classmethod
    def load_regions(
        cls,
        image_collection: ImageInput,
        sources: Sources,
        *,
        root: str | Path,
    ) -> PhotometryRegionCollection:
        """Read existing per-source/filter regions without generating anything.

        Uses the same FIELD, source-role grouping, filter ordering and naming
        rules as generate_regions. The caller supplies the existing image
        collection, current Sources and analysis root; no manual manifest paths
        are needed. Retains saved centers, radii, rejection and review records.

        Missing files, image associations or source identities raise an explicit
        error. No centroiding, scaling, GUI, file writes or automatic replacement
        occurs. Sources/reference magnitudes are not restored or modified: these
        are supplied by the caller. Returned records/review() and photometry
        inputs work exactly as with newly generated regions.
        """
        root_path = Path(root).expanduser().resolve()
        plans = cls._plan_regions(
            image_collection, sources, root_path, load_existing=True
        )
        sessions: list[PhotometryRegions] = []
        for image, group in plans:
            directory, stem, _ = cls._output_layout(fits.getheader(image), group)
            session = cls.load(root_path, directory / f"{stem}_regions.json")
            if session.image_path != image:
                raise ValueError(
                    f"Saved regions refer to a different image: {session.image_path}."
                )
            expected = {
                (info.observation_field.field_name, info.source_name, info.source_role)
                for info in group.source_infos
            }
            recorded = {
                (row["field_name"], row["source_id"], row["source_role"])
                for row in session.records
            }
            if recorded != expected:
                raise ValueError(
                    f"Saved region source identities do not match Sources for {image.name}."
                )
            sessions.append(session)
        return PhotometryRegionCollection(sessions)

    @classmethod
    def generate_regions(
        cls,
        image_collection: ImageInput,
        sources: Sources,
        *,
        root: str | Path,
        aperture_scaling: float = 3.0,
        inner_annulus_scaling: float = 5.0,
        outer_annulus_scaling: float = 8.0,
        box_size: int | None = None,
        max_centroid_shift: float | None = None,
        overwrite: bool = False,
    ) -> PhotometryRegionCollection:
        """Generate all supplied sources' Master Light regions in one operation.

        Select existing collection members by IMTYPE='Master Light' and each
        source's FIELD, without rescanning directories or requiring FITS paths
        in ObservationField. Both ImageFileCollection and ImageCollectionGroup
        are supported. Standards share one pair per field/filter; each science
        source gets its own pair per field/filter, even in a shared field.

        Order follows source field order, source group, then Telescope filter
        order and image path, keeping each target's filters together. Every
        image uses its own FWHM and centroiding.
        box_size=None retains the old notebook's 19-pixel standard and 25-pixel
        science windows. Explicit box_size applies to all groups. Scaling and
        maximum-displacement controls follow generate.

        Reject missing field images, filename collisions and existing outputs
        before writing any group. A later calculation/I/O failure can leave
        earlier groups written; this operation is not a batch transaction.
        No inputs or source results are changed. Return a collection with one
        review() call and a combined records list for notebook table display.
        """
        root_path = Path(root).expanduser().resolve()
        plans = cls._plan_regions(
            image_collection, sources, root_path, overwrite=overwrite
        )
        sessions = [
            cls.generate(
                path,
                group,
                root=root_path,
                aperture_scaling=aperture_scaling,
                inner_annulus_scaling=inner_annulus_scaling,
                outer_annulus_scaling=outer_annulus_scaling,
                box_size=(
                    (19 if group.source_infos[0].source_role == "standard" else 25)
                    if box_size is None
                    else box_size
                ),
                max_centroid_shift=max_centroid_shift,
                overwrite=overwrite,
            )
            for path, group in plans
        ]
        return PhotometryRegionCollection(sessions)

    @classmethod
    def generate(
        cls,
        image_path: str | Path,
        sources: Sources,
        *,
        root: str | Path,
        aperture_scaling: float = 3.0,
        inner_annulus_scaling: float = 5.0,
        outer_annulus_scaling: float = 8.0,
        box_size: int = 19,
        max_centroid_shift: float | None = None,
        overwrite: bool = False,
    ) -> PhotometryRegions:
        """Centroid sources and write one circle/annulus file pair for an image.

        All three scaling factors multiply this image's positive header FWHM
        in pixels. They must satisfy 0 < aperture < inner annulus < outer
        annulus. box_size is an odd local centroid window, default 19 pixels as
        in the standard-field notebook. max_centroid_shift is measured from the
        initial position in pixels; None permits at most box_size // 2 pixels.

        Retains the old Background2D recipe (20-pixel boxes, 5-box filtering,
        MedianBackground, 3-sigma clipping), two centroid_com iterations, and a
        final centroid_quadratic. A failed quadratic tries centroid_com; failed,
        non-finite, masked or excessively shifted results revert to the initial
        WCS projection. Background-estimation failure also falls back and is
        recorded. MASK and non-finite values are excluded; CMASK is ignored.

        Both standard and science Sources are supported, within one image FIELD.
        All standards for the field may share a pair; a single science source
        gets its own pair. Use generate_regions to split multiple science
        sources automatically. Standard files share a field_standards directory;
        science files share a source-name directory and include FIELD in their
        filenames. All filenames identify FILTER. Comments record field, source
        IDs and the root-relative image path. No file path is stored in SourceInfo.
        No flux, uncertainty, or detection significance is measured. Valid-area
        fractions are diagnostics only: partial regions remain for manual review.
        They do not indicate whether all stack exposures contributed to a pixel.
        Image data, SourceInfo, and magnitude results are never modified.
        """
        scales = np.asarray(
            [aperture_scaling, inner_annulus_scaling, outer_annulus_scaling], float
        )
        if not np.all(np.isfinite(scales)) or not (
            0 < scales[0] < scales[1] < scales[2]
        ):
            raise ValueError(
                "Require 0 < aperture_scaling < inner_annulus_scaling < outer_annulus_scaling."
            )
        if (
            isinstance(box_size, bool)
            or not isinstance(box_size, int)
            or box_size < 3
            or box_size % 2 != 1
        ):
            raise ValueError("box_size must be an odd integer >= 3.")
        shift_limit = (
            box_size // 2 if max_centroid_shift is None else float(max_centroid_shift)
        )
        if not np.isfinite(shift_limit) or shift_limit <= 0:
            raise ValueError("max_centroid_shift must be positive and finite.")
        if not len(sources):
            raise ValueError("At least one source is required.")

        root_path = Path(root).expanduser().resolve()
        image = Path(image_path).expanduser().resolve()
        image_relative = image.relative_to(root_path)
        header = fits.getheader(image)
        directory, stem, filter_name = cls._output_layout(header, sources)
        manifest = directory / f"{stem}_regions.json"
        state = {
            "version": 1,
            "image_path": image_relative.as_posix(),
            "field_name": sources.source_infos[0].observation_field.field_name,
            "filter_name": filter_name,
            "aperture_path": (directory / f"{stem}_aperture.reg").as_posix(),
            "background_path": (directory / f"{stem}_background.reg").as_posix(),
            "scaling_factors": {
                "aperture": float(scales[0]),
                "annulus_inner": float(scales[1]),
                "annulus_outer": float(scales[2]),
            },
            "sources": [],
            "reviewed": False,
        }
        session = cls(root_path, manifest, state)
        for path in (
            session.aperture_path,
            session.background_path,
            session._resolve(manifest),
        ):
            if path.exists() and not overwrite:
                raise FileExistsError(
                    f"{path}; use the existing session.review() to revisit edits, "
                    "or overwrite=True to regenerate explicitly."
                )

        with fits.open(image, mode="readonly", memmap=False) as hdul:
            data = np.asarray(hdul[0].data, dtype=float)
            header = hdul[0].header
            if data.ndim != 2:
                raise ValueError("Expected a primary 2D image.")
            invalid = ~np.isfinite(data)
            if "MASK" in hdul:
                flags = hdul["MASK"].data
                if flags is None or flags.shape != data.shape:
                    raise ValueError("MASK shape must match the primary image.")
                invalid |= np.asarray(flags, dtype=bool)
        fwhm = float(header["FWHM"])
        if not np.isfinite(fwhm) or fwhm <= 0:
            raise ValueError("Image FWHM must be positive and finite in pixels.")
        wcs = WCS(header)
        if not wcs.has_celestial:
            raise ValueError("Image requires celestial WCS.")
        state["fwhm"] = fwhm
        state["centroid_controls"] = {"box_size": box_size, "max_shift": shift_limit}
        radii = scales * fwhm
        background_error = ""
        try:
            background = Background2D(
                data,
                box_size=(20, 20),
                filter_size=(5, 5),
                mask=invalid,
                sigma_clip=SigmaClip(sigma=3),
                bkg_estimator=MedianBackground(),
            )
            signal = data - background.background
        except (ValueError, TypeError) as error:
            signal = None
            background_error = f"Background estimation failed: {error}"

        for info in sources.source_infos:
            if info.observation_field.field_name != header.get("FIELD"):
                raise ValueError("Source field does not match image FIELD.")
            name = info.source_name
            if any(character in name for character in "{}\\\n\r"):
                raise ValueError("Source ID contains unsupported DS9 tag characters.")
            x, y = skycoord_to_pixel(info.source_coordinate, wcs, origin=0, mode="all")
            initial = np.array([float(x), float(y)])
            if not np.all(np.isfinite(initial)) or not (
                0 <= x < data.shape[1] and 0 <= y < data.shape[0]
            ):
                raise ValueError(
                    f"Initial source coordinate is outside the image: {name!r}."
                )
            center = initial.copy()
            method, reason = "initial_fallback", background_error
            if signal is not None:
                try:

                    def locate(position, function):
                        with warnings.catch_warnings(record=True) as messages:
                            warnings.simplefilter("always")
                            xx, yy = centroid_sources(
                                signal,
                                [position[0]],
                                [position[1]],
                                box_size=box_size,
                                centroid_func=function,
                                mask=invalid,
                            )
                        result = np.array([xx[0], yy[0]])
                        if messages or not np.all(np.isfinite(result)):
                            raise ValueError(
                                "Centroid fit warned or returned a non-finite position."
                            )
                        if np.linalg.norm(result - initial) > shift_limit:
                            raise ValueError(
                                "Centroid exceeds maximum displacement from initial position."
                            )
                        ix, iy = np.rint(result).astype(int)
                        if (
                            not (0 <= ix < data.shape[1] and 0 <= iy < data.shape[0])
                            or invalid[iy, ix]
                        ):
                            raise ValueError(
                                "Centroid is outside the valid image area."
                            )
                        return result

                    for _ in range(2):
                        center = locate(center, centroid_com)
                    try:
                        center = locate(center, centroid_quadratic)
                        method = "quadratic"
                    except (ValueError, TypeError) as error:
                        center = locate(center, centroid_com)
                        method = "com_fallback"
                        reason = f"Quadratic fit failed; used COM: {error}"
                except (ValueError, TypeError) as error:
                    center = initial.copy()
                    reason = str(error)
            aperture = CircularAperture(center, r=radii[0])
            annulus = CircularAnnulus(center, r_in=radii[1], r_out=radii[2])
            state["sources"].append(
                {
                    "source_id": name,
                    "field_name": info.observation_field.field_name,
                    "source_role": info.source_role,
                    "aperture_path": state["aperture_path"],
                    "background_path": state["background_path"],
                    "initial_x": initial[0],
                    "initial_y": initial[1],
                    "auto_x": center[0],
                    "auto_y": center[1],
                    "x": center[0],
                    "y": center[1],
                    "background_x": center[0],
                    "background_y": center[1],
                    "center_offset": 0.0,
                    "centroid_method": method,
                    "centroid_reason": reason,
                    "aperture_radius": float(radii[0]),
                    "inner_radius": float(radii[1]),
                    "outer_radius": float(radii[2]),
                    "generated_radii": radii.tolist(),
                    "accepted": True,
                    "reviewed": False,
                    "aperture_valid_fraction": cls._valid_fraction(aperture, invalid),
                    "annulus_valid_fraction": cls._valid_fraction(annulus, invalid),
                }
            )
        session._publish(state)
        return session

    @staticmethod
    def _region_text(
        records: list[dict], kind: str, metadata: dict | None = None
    ) -> str:
        """Serialize circles or annuli with stable tags and DS9 1-based pixels."""
        lines = ["# Region file format: DS9 version 4.1", "image"]
        if metadata is not None:
            for key in ("field_name", "filter_name", "image_path"):
                if key in metadata:
                    lines.append(
                        f"# {key.upper()}={json.dumps(metadata[key], ensure_ascii=False)}"
                    )
            lines.append(
                "# SOURCE_NAMES="
                + json.dumps(
                    [row["source_id"] for row in records if row["accepted"]],
                    ensure_ascii=False,
                )
            )
        for row in records:
            if not row["accepted"]:
                continue
            x, y = row["x"] + 1, row["y"] + 1
            if kind == "aperture":
                shape = f"circle({x:.10f},{y:.10f},{row['aperture_radius']:.10f})"
                color = "green"
            else:
                x = row.get("background_x", row["x"]) + 1
                y = row.get("background_y", row["y"]) + 1
                shape = f"annulus({x:.10f},{y:.10f},{row['inner_radius']:.10f},{row['outer_radius']:.10f})"
                color = "red"
            label = f" text={{{row['source_id']}}}" if kind == "aperture" else ""
            lines.append(
                f"{shape} # color={color} tag={{id:{row['source_id']}}} tag={{{kind}}}{label}"
            )
        return "\n".join(lines) + "\n"

    def _publish(self, state: dict) -> None:
        """Stage all files before replacement; publish the manifest last.

        Individual replacements are atomic. The three-file operation is not a
        filesystem transaction; an I/O failure is raised and navigation stops.
        """
        payloads = {
            self.aperture_path: self._region_text(state["sources"], "aperture", state),
            self.background_path: self._region_text(
                state["sources"], "background", state
            ),
            self._resolve(self.manifest_path): json.dumps(
                state, indent=2, allow_nan=False
            ),
        }
        staged: list[tuple[Path, Path]] = []
        try:
            for destination, text in payloads.items():
                destination.parent.mkdir(parents=True, exist_ok=True)
                descriptor, temporary = tempfile.mkstemp(
                    dir=destination.parent, prefix=".regions-"
                )
                path = Path(temporary)
                staged.append((path, destination))
                with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
                    stream.write(text)
            for temporary, destination in staged:
                os.replace(temporary, destination)
        finally:
            for temporary, _ in staged:
                temporary.unlink(missing_ok=True)
        self._state = state

    def initial_regions(self) -> str:
        """Return fixed initial-position crosses, excluded from measured files."""
        lines = ["image"]
        for row in self._state["sources"]:
            lines.append(
                f"point({row['initial_x'] + 1:.10f},{row['initial_y'] + 1:.10f}) "
                "# point=cross 15 color=blue width=2 tag={initial} select=0 move=0 edit=0 delete=0"
            )
        return "\n".join(lines) + "\n"

    def background_display_regions(self) -> str:
        """Load annuli without text labels, including previously generated files.

        Strip display text in memory only; source ID tags and geometry remain
        unchanged. Only the aperture displays the star name in the DS9 frame.
        """
        regions = Regions.read(self.background_path, format="ds9")
        for region in regions:
            region.meta.pop("text", None)
        return regions.serialize(format="ds9") if regions else "image\n"

    def save_review(self, region_text: str) -> list[str]:
        """Validate DS9 edits and save paired regions plus per-source decisions.

        DS9 exports must use image coordinates and preserve id/type tags. Both
        regions absent means rejection for this image. A missing counterpart,
        duplicate/unknown ID, excluded shape, non-circular geometry or invalid
        radius blocks saving. Independent finite centers and manual radius
        changes are retained. Return non-blocking warnings for offsets above
        0.01 pixel and aperture/sky overlap. Initial crosses are never written
        to files. Saved records contain both centers and the warning messages.
        """
        pairs: dict[str, dict] = {}
        known = {row["source_id"] for row in self._state["sources"]}
        with warnings.catch_warnings(record=True) as messages:
            warnings.simplefilter("always")
            parsed = Regions.parse(region_text, format="ds9")
        if messages:
            raise ValueError(
                "DS9 regions could not be parsed completely; export image coordinates."
            )
        for region in parsed:
            tags = region.meta.get("tag", [])
            if isinstance(region, PointPixelRegion) and "initial" in tags:
                continue
            ids = [tag[3:] for tag in tags if tag.startswith("id:")]
            kinds = [tag for tag in tags if tag in ("aperture", "background")]
            if len(ids) != 1 or ids[0] not in known or len(kinds) != 1:
                raise ValueError(
                    "Every measured region requires a known source ID and one role tag."
                )
            kind = kinds[0]
            expected = (
                CirclePixelRegion if kind == "aperture" else CircleAnnulusPixelRegion
            )
            if not isinstance(region, expected) or not region.meta.get("include", True):
                raise ValueError(
                    "Expected included circular aperture/annulus regions in image coordinates."
                )
            pair = pairs.setdefault(ids[0], {})
            if kind in pair:
                raise ValueError(f"Duplicate {kind} for {ids[0]!r}.")
            pair[kind] = region

        state = deepcopy(self._state)
        review_warnings: list[str] = []
        for row in state["sources"]:
            pair = pairs.get(row["source_id"], {})
            if pair and set(pair) != {"aperture", "background"}:
                raise ValueError(
                    f"Incomplete aperture/annulus pair for {row['source_id']!r}."
                )
            if pair:
                aperture, background = pair["aperture"], pair["background"]
                xy = np.array([aperture.center.x, aperture.center.y], float)
                bxy = np.array([background.center.x, background.center.y], float)
                radii = np.array(
                    [aperture.radius, background.inner_radius, background.outer_radius],
                    float,
                )
                if not np.all(np.isfinite(np.r_[xy, bxy, radii])):
                    raise ValueError(
                        "Aperture and annulus centers/radii must be finite."
                    )
                if not (0 < radii[0] < radii[1] < radii[2]):
                    raise ValueError(
                        "Require 0 < aperture radius < annulus inner radius < outer radius."
                    )
                row.update(
                    x=float(xy[0]),
                    y=float(xy[1]),
                    background_x=float(bxy[0]),
                    background_y=float(bxy[1]),
                    center_offset=float(np.linalg.norm(xy - bxy)),
                    aperture_radius=float(radii[0]),
                    inner_radius=float(radii[1]),
                    outer_radius=float(radii[2]),
                )
                if row["center_offset"] > 0.01:
                    review_warnings.append(
                        f"{row['source_id']}: background center offset "
                        f"{row['center_offset']:.2f} pixels; saved both positions."
                    )
                distance = row["center_offset"]
                if distance + radii[0] > radii[1] and distance - radii[0] < radii[2]:
                    review_warnings.append(
                        f"{row['source_id']}: background annulus overlaps the source aperture; "
                        "check sky contamination before photometry."
                    )
            row.update(accepted=bool(pair), reviewed=True)
        # Recompute coverage diagnostics at manually edited positions/radii.
        with fits.open(self.image_path, mode="readonly", memmap=False) as hdul:
            invalid = ~np.isfinite(hdul[0].data)
            if "MASK" in hdul:
                invalid |= np.asarray(hdul["MASK"].data, dtype=bool)
        for row in state["sources"]:
            center = (row["x"], row["y"])
            row["aperture_valid_fraction"] = self._valid_fraction(
                CircularAperture(center, row["aperture_radius"]), invalid
            )
            row["annulus_valid_fraction"] = self._valid_fraction(
                CircularAnnulus(
                    (
                        row.get("background_x", row["x"]),
                        row.get("background_y", row["y"]),
                    ),
                    row["inner_radius"],
                    row["outer_radius"],
                ),
                invalid,
            )
        state["reviewed"] = True
        state["review_warnings"] = review_warnings
        self._publish(state)
        return review_warnings

    def review(self, *, isolated: bool | None = None) -> None:
        """Open the editable DS9 reviewer; Next saves, Quit discards unsaved edits."""
        from ..interactive_ds9.region_reviewer import RegionReviewer

        try:
            RegionReviewer([self], isolated=isolated).run()
        finally:
            self._state = self.load(self.root, self.manifest_path)._state


class PhotometryRegionCollection:
    """Group automatically generated per-image/source sessions for one review.

    Standards remain grouped by field; science targets have individual files.
    sessions preserves generation order, and records provides one combined list
    with image/filter context suitable for pandas.DataFrame. No magnitude
    measurements are combined or overwritten by this container.
    """

    def __init__(self, sessions: list[PhotometryRegions]) -> None:
        self.sessions = tuple(sessions)

    def __len__(self) -> int:
        """Return the number of image/source-group review entries."""
        return len(self.sessions)

    @property
    def records(self) -> list[dict]:
        """Return copied source records with relative image and filter names."""
        return [
            {
                **row,
                "image_path": session._state["image_path"],
                "filter_name": session._state["filter_name"],
            }
            for session in self.sessions
            for row in session.records
        ]

    def review(self, *, isolated: bool | None = None) -> None:
        """Review the entire queue; Save & Next saves the current file pair."""
        from ..interactive_ds9.region_reviewer import RegionReviewer

        try:
            RegionReviewer(self.sessions, isolated=isolated).run()
        finally:
            for session in self.sessions:
                session._state = PhotometryRegions.load(
                    session.root,
                    session.manifest_path,
                )._state
