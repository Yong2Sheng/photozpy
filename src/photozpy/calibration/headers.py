"""Standardize FITS headers from observation-field and telescope metadata.

This module coordinates filename classification, field matching, and header
updates. It uses collection helpers for refresh and ``update_fits_header`` for
file writing; it does not require identified sources or magnitude results.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from datetime import date
from fnmatch import fnmatchcase
from pathlib import Path
import re
from typing import Any, Literal, TypeAlias

import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from ccdproc import ImageFileCollection
from tqdm.auto import tqdm

from ..collections import ImageCollectionGroup
from ..collections.image_collection_helpers import (
    refresh_image_collection,
    refresh_image_collection_group,
)
from ..fits_utilities import update_fits_header
from ..sources.observation_field import ObservationField
from ..telescope import Telescope


ImageInput: TypeAlias = ImageFileCollection | ImageCollectionGroup
ImageType: TypeAlias = Literal["Bias", "Dark", "Flat", "Light"]

# Astropy accepts both plain card values and (value, comment) tuples. Keep
# their value type flexible and let fits.Header validate each planned card.
HeaderCards: TypeAlias = dict[str, Any]
HeaderPlan: TypeAlias = tuple[ImageFileCollection, Path, HeaderCards]


class HeaderCorrection:
    """Coordinate FITS-header standardization for one telescope session.

    Parameters
    ----------
    image_collection : ImageFileCollection or ImageCollectionGroup
        Input images. Group children are processed in their original order.
    observation_fields : ObservationField or iterable of ObservationField, optional
        Science and standard observation fields used to match light filenames.
        A standard field does not need identified stars. Repeated definitions
        of the same field name must reuse the same ObservationField instance.
        Default is an empty tuple, suitable for calibration-only processing.
    telescope : Telescope or None, optional
        Telescope configuration for the session. If omitted, use the first
        observation field's telescope. All fields must reference this same
        instance. Required when no observation fields are supplied.
    telescope_name : str or None, optional
        Explicit value to write to TELESCOP during filename correction.
        If None, leave the existing card unchanged, including its absence.
    time_scale : str or None, optional
        Astropy time scale used only when a file has no TIMESYS declaration.
        Must agree with TIMESYS when both are supplied. If None, each file
        must declare a supported TIMESYS. No scale conversion is performed.
    save_location : str, Path, or None, optional
        Output directory. None or an empty string means in-place output.
        Processing methods may override this setting for an individual call.
    overwrite : bool, optional
        Allow replacement of existing destination files. Default is False;
        in-place header changes therefore require explicit overwrite=True.

    Raises
    ------
    TypeError
        An input collection, field, or telescope has an unsupported type, or
        neither telescope nor observation fields provide a Telescope instance.
    ValueError
        Field definitions conflict, fields use different telescope instances,
        telescope_name is empty, or time_scale is unsupported.

    Notes
    -----
    Construction does not write files. Processing methods build and validate
    header updates before writing, then return the updated collection and store
    it in image_collection. Output collections contain fixed file selections.

    FIELD identifies the observation field; OBJECT is preserved. Missing
    pointing may be supplied from nominal_pointing in ICRS, which is distinct
    from an image center measured from WCS. CCD parameters are selected using
    each image's observation date. MJD-OBS records the exposure start without
    a midpoint adjustment or time-scale conversion.

    Image arrays and all HDUs are retained. Individual output files are written
    atomically, but a batch is not a transaction if a later disk operation fails.
    Shared output directories require unique filenames across group children.
    """

    def __init__(
        self,
        image_collection: ImageInput,
        *,
        observation_fields: ObservationField | Iterable[ObservationField] = (),
        telescope: Telescope | None = None,
        telescope_name: str | None = None,
        time_scale: str | None = None,
        save_location: str | Path | None = None,
        overwrite: bool = False,
    ) -> None:
        """Validate session metadata and retain the input collection."""
        if not isinstance(
            image_collection, (ImageFileCollection, ImageCollectionGroup)
        ):
            raise TypeError(
                "image_collection must be an ImageFileCollection "
                "or ImageCollectionGroup."
            )

        # Normalize fields and reject inconsistent definitions of one name.
        if isinstance(observation_fields, ObservationField):
            observation_fields = (observation_fields,)

        fields: list[ObservationField] = list(observation_fields)
        unique_fields: dict[str, ObservationField] = {}

        for field in fields:
            if not isinstance(field, ObservationField):
                raise TypeError(
                    "Every observation field must be an ObservationField."
                )

            existing: ObservationField | None = unique_fields.get(field.field_name)
            if existing is not None and existing is not field:
                raise ValueError(
                    f"Conflicting field definitions for {field.field_name!r}; "
                    "reuse the same ObservationField instance."
                )

            unique_fields[field.field_name] = field

        self._fields: tuple[ObservationField, ...] = tuple(unique_fields.values())

        # One session uses one telescope configuration, including calibration
        # frames that cannot be associated with a particular light field.
        if telescope is None and self._fields:
            telescope = self._fields[0].telescope

        if not isinstance(telescope, Telescope):
            raise TypeError("Supply a Telescope or at least one observation field.")

        if any(field.telescope is not telescope for field in self._fields):
            raise ValueError("All fields must use the session's Telescope instance.")

        if telescope_name is not None and (
            not isinstance(telescope_name, str) or not telescope_name.strip()
        ):
            raise ValueError("telescope_name must be a non-empty string or None.")

        if time_scale is not None:
            time_scale = time_scale.lower()
            if time_scale not in Time.SCALES:
                raise ValueError(f"Unsupported time scale: {time_scale!r}.")

        self._image_collection: ImageInput = image_collection
        self._telescope: Telescope = telescope
        self._telescope_name: str | None = telescope_name
        self._time_scale: str | None = time_scale
        self._save_location: str | Path | None = save_location
        self._overwrite: bool = overwrite
        self._skipped_files: tuple[Path, ...] = ()

    @property
    def image_collection(self) -> ImageInput:
        """Return the input collection or the latest successful output.

        A successful processing call replaces this reference with its returned
        collection. Merely accessing the property does not refresh or write it.
        """
        return self._image_collection

    @property
    def telescope(self) -> Telescope:
        """Return the telescope configuration shared by all session fields."""
        return self._telescope

    @property
    def skipped_files(self) -> tuple[Path, ...]:
        """Return absolute paths skipped by the last successful filename call.

        This tuple is updated only after correct_headers_by_filename succeeds.
        Filter correction does not replace the previous skip record.
        """
        return self._skipped_files

    @staticmethod
    def _children(collection: ImageInput) -> tuple[ImageFileCollection, ...]:
        """Expose either input form as an ordered tuple of child collections.

        Existing collections are returned by reference, without reconstruction.
        Raise TypeError for unsupported inputs.
        """
        if isinstance(collection, ImageCollectionGroup):
            return collection.collections

        if isinstance(collection, ImageFileCollection):
            return (collection,)

        raise TypeError("Expected an ImageFileCollection or ImageCollectionGroup.")

    @staticmethod
    def _refresh(collection: ImageInput) -> None:
        """Refresh input collections in place using the collection helpers.

        Directory discovery remains governed by each child's configuration;
        refreshing an explicit file selection does not broaden that selection.
        """
        if isinstance(collection, ImageCollectionGroup):
            refresh_image_collection_group(collection)

        elif isinstance(collection, ImageFileCollection):
            refresh_image_collection(collection)

        else:
            raise TypeError("Expected an ImageFileCollection or ImageCollectionGroup.")

    @classmethod
    def _files(
        cls,
        collection: ImageInput,
    ) -> Iterator[tuple[ImageFileCollection, Path]]:
        """Yield each owning collection and absolute file path in input order.

        Resolve relative filenames against their collection location and omit
        macOS AppleDouble sidecars beginning with ``._``. No files are written.
        """
        for child in cls._children(collection):
            for filename in child.files:
                path: Path = Path(filename)

                if not path.is_absolute():
                    path = Path(child.location or '.') / path

                if not path.name.startswith('._'):
                    yield child, path.resolve()

    def _observation_time(self, header: fits.Header, path: Path) -> Time:
        """Interpret DATE-OBS in its declared or explicitly supplied time scale.

        path is used in diagnostic messages. Raise ValueError for missing or
        invalid DATE-OBS, unsupported or absent TIMESYS without a fallback, or
        a conflict between TIMESYS and the session's time_scale.
        """
        declared: str = str(header.get("TIMESYS", "")).strip().lower()
        if declared and self._time_scale and declared != self._time_scale:
            raise ValueError(f"{path.name}: TIMESYS conflicts with time_scale.")

        scale: str | None = declared or self._time_scale
        if scale not in Time.SCALES:
            raise ValueError(
                f"{path.name}: missing or unsupported TIMESYS; "
                "supply time_scale explicitly when TIMESYS is absent."
            )

        if "DATE-OBS" not in header:
            raise ValueError(f"{path.name}: missing DATE-OBS.")

        try:
            return Time(str(header['DATE-OBS']).strip(), format="isot", scale=scale)

        except ValueError as exc:
            raise ValueError(f"{path.name}: invalid DATE-OBS.") from exc

    def _match_field(self, filename: str) -> ObservationField | None:
        """Find one field by its case-sensitive filename substring or glob.

        Fields with file_pattern=None are ignored. Return None for no match;
        raise ValueError when more than one field matches the same basename.
        """
        matches: list[ObservationField] = []

        for field in self._fields:
            pattern: str | None = field.file_pattern
            if pattern is None:
                continue

            # A glob applies to the complete basename; a plain pattern keeps
            # the legacy substring meaning used by the existing field setup.
            matched: bool = (
                fnmatchcase(filename, pattern) if any(c in pattern for c in '*?[')
                else pattern in filename
            )
            if matched:
                matches.append(field)

        if len(matches) > 1:
            raise ValueError(f"{filename}: ambiguous field patterns.")

        return matches[0] if matches else None

    @staticmethod
    def _image_type(filename: str, forced: ImageType | None) -> ImageType:
        """Use an explicit type or infer a calibration token from the basename.

        Bias, Dark, and Flat tokens are matched case-insensitively with filename
        boundaries. Other names are classified as Light and later require field
        matching. The public caller validates any explicitly forced type.
        """
        if forced is not None:
            return forced

        match: re.Match[str] | None = re.search(
            r"(?:^|[_\-\s])(bias|dark|flat)(?=[_\-\s\d.]|$)",
            filename, re.IGNORECASE,
        )

        return match.group(1).capitalize() if match else "Light"

    def correct_headers_by_filename(
        self,
        save_location: str | Path | None = None,
        overwrite: bool | None = None,
        *,
        image_type: ImageType | None = None,
        time_transformation: bool = True,
        overwrite_pointing: bool = False,
        unmatched: Literal["raise", "skip"] = "raise",
    ) -> ImageInput:
        """Classify images and write field, CCD, and optional time metadata.

        Parameters
        ----------
        save_location : str, Path, or None, optional
            Override the session output directory. None uses the constructor
            setting; an empty string explicitly selects in-place output.
        overwrite : bool or None, optional
            Override the session overwrite setting. None uses that setting.
        image_type : {"Bias", "Dark", "Flat", "Light"} or None, optional
            Force one type for every image. Otherwise infer calibration types
            from filenames and treat remaining files as Light. Light images
            must match exactly one observation field in either case.
        time_transformation : bool, optional
            Write exposure-start MJD-OBS and its TIMESYS. Default is True.
            DATE-OBS and its time scale are still required when False because
            the observation date selects the valid CCD configuration.
        overwrite_pointing : bool, optional
            Replace RA and DEC with nominal field pointing in ICRS. Default is
            False: retain both existing cards, or fill both when absent. Having
            only one card is an error unless replacement is explicitly enabled.
        unmatched : {"raise", "skip"}, optional
            How to handle light filenames that match no field. Default is
            "raise". Skipping omits those files from the output and records
            their paths in skipped_files. Ambiguous matches always raise.

        Returns
        -------
        ImageFileCollection or ImageCollectionGroup
            Written files as explicit selections, with the same single/group
            form as the session input. Also stored in image_collection.

        Raises
        ------
        ValueError
            Classification or field matching fails, dates or time scales are
            invalid, no valid CCD exists, BUNIT is present but not ADU, pointing
            is incomplete, no files remain, or output paths collide.
        FileExistsError
            A destination exists and overwriting was not enabled.
        OSError
            A file cannot be read or written.

        Notes
        -----
        OBJECT is never replaced. FIELD is written only for matched light
        fields. GAIN and RDNOISE are selected per observation date; no image
        array is converted to a different unit. TELESCOP is changed only when
        telescope_name was supplied. Filter normalization is a separate call.

        File matches, card values, and destination paths are checked before
        writing. A later I/O failure may still leave earlier outputs written.
        """
        if image_type not in (None, "Bias", "Dark", "Flat", "Light"):
            raise ValueError("image_type must be Bias, Dark, Flat, Light, or None.")

        if unmatched not in ("raise", "skip"):
            raise ValueError("unmatched must be 'raise' or 'skip'.")

        self._refresh(self._image_collection)
        plans: list[HeaderPlan] = []
        skipped: list[Path] = []

        # Build every update before _apply is allowed to write any file.
        for child, path in self._files(self._image_collection):
            kind: ImageType = self._image_type(path.name, image_type)
            field: ObservationField | None = (
                self._match_field(path.name) if kind == "Light" else None
            )

            if kind == "Light" and field is None:
                if unmatched == "skip":
                    skipped.append(path)
                    continue

                raise ValueError(f"{path.name}: no observation field matches.")

            header: fits.Header = fits.getheader(path, ext=child.ext)
            if 'BUNIT' in header and str(header['BUNIT']).strip().lower() != 'adu':
                raise ValueError(
                    f"{path.name}: BUNIT is not ADU; header correction does not "
                    "convert image data to ADU."
                )

            # Resolve the CCD on this image's date, not the directory's date.
            observed: Time = self._observation_time(header, path)
            obs_date: date = observed.to_value('datetime').date()
            gain = self._telescope.gain(obs_date).to_value(u.electron / u.adu)
            rdnoise = self._telescope.rdnoise(obs_date).to_value(u.electron / u.pix)

            cards: HeaderCards = {
                "IMTYPE": kind,
                "GAIN": (float(gain), "CCD gain [electron/ADU]"),
                "RDNOISE": (float(rdnoise), "Per-pixel RMS read noise [electron/pixel]"),
                "BUNIT": "ADU",
            }

            if self._telescope_name is not None:
                cards['TELESCOP'] = self._telescope_name

            if time_transformation:
                cards['MJD-OBS'] = (float(observed.mjd), "Start of observation in TIMESYS")
                cards['TIMESYS'] = observed.scale.upper()

            if field is not None:
                cards['FIELD'] = (field.field_name, "ObservationField identifier")
                has_ra, has_dec = 'RA' in header, 'DEC' in header

                if not overwrite_pointing and has_ra != has_dec:
                    raise ValueError(
                        f"{path.name}: incomplete RA/DEC; set overwrite_pointing=True "
                        "to replace both with nominal field pointing."
                    )

                if overwrite_pointing or not has_ra:
                    ra, dec = field.nominal_pointing.icrs.to_string(
                        'hmsdms', sep=':', precision=3, pad=True
                    ).split()
                    cards['RA'] = (ra, "Nominal field pointing RA (ICRS)")
                    cards['DEC'] = (dec, "Nominal field pointing Dec (ICRS)")

            plans.append((child, path, cards))

        result: ImageInput = self._apply(
            self._image_collection, plans, save_location, overwrite
        )
        self._skipped_files = tuple(skipped)

        return result

    def correct_filter_headers(
        self,
        input_collection: ImageInput | None = None,
        overwrite: bool | None = None,
        save_location: str | Path | None = None,
        *,
        from_file_name: bool = False,
        filter_mapping: Mapping[str, str] | None = None,
    ) -> ImageInput:
        """Normalize FILTER on Flat and Light images using telescope aliases.

        Parameters
        ----------
        input_collection : ImageFileCollection, ImageCollectionGroup, or None, optional
            Images to process. None uses the current session collection.
            IMTYPE must already be populated, normally by filename correction.
        overwrite : bool or None, optional
            Override the session overwrite setting. None uses that setting.
        save_location : str, Path, or None, optional
            Override the session output directory. None uses the constructor
            setting; an empty string explicitly selects in-place output.
        from_file_name : bool, optional
            If False, resolve the existing FILTER value. If True, identify one
            filter in each Flat/Light basename. Default is False.
        filter_mapping : mapping of str to str or None, optional
            Header mode: remap an exact FILTER value before alias resolution;
            values absent from the mapping are resolved unchanged.
            Filename mode: match mapping keys as substrings of the basename
            and resolve their target values. Without a mapping, filename mode
            searches canonical names, allowing omission of a final apostrophe.
            Every target must resolve through the session Telescope.

        Returns
        -------
        ImageFileCollection or ImageCollectionGroup
            All selected files, retaining non-Flat/Light images unchanged.
            The result follows input_collection's single/group form and also
            becomes the session's image_collection.

        Raises
        ------
        ValueError
            IMTYPE or a required FILTER is absent, a filter cannot be resolved,
            filename matching is ambiguous, no files remain, or output paths
            collide. A batch containing only Flat or only Light is valid.
        FileExistsError
            A destination exists and overwriting was not enabled.
        OSError
            A file cannot be read or written.

        Notes
        -----
        Filter labels are standardized without converting magnitudes, fluxes,
        or image arrays. Other image types are copied unchanged when writing
        to a different directory. In-place processing leaves them untouched.
        """
        collection: ImageInput = (
            self._image_collection if input_collection is None else input_collection
        )
        self._refresh(collection)
        plans: list[HeaderPlan] = []

        for child, path in self._files(collection):
            header: fits.Header = fits.getheader(path, ext=child.ext)
            kind: str = str(header.get('IMTYPE', '')).strip().lower()

            if not kind:
                raise ValueError(
                    f"{path.name}: missing IMTYPE; correct image types first."
                )

            cards: HeaderCards = {}
            if kind in ('flat', 'light'):
                if from_file_name:
                    standard = self._filter_from_filename(path.name, filter_mapping)

                else:
                    if 'FILTER' not in header:
                        raise ValueError(f"{path.name}: missing FILTER.")

                    alias: str = str(header['FILTER']).strip()
                    if filter_mapping is not None:
                        alias = filter_mapping.get(alias, alias)

                    standard = self._telescope.resolve_filter_alias(alias)

                cards['FILTER'] = standard

            plans.append((child, path, cards))

        return self._apply(collection, plans, save_location, overwrite)

    def _filter_from_filename(
        self,
        filename: str,
        mapping: Mapping[str, str] | None,
    ) -> str:
        """Resolve exactly one canonical filter from a filename.

        Mapping keys use substring matching. Without a mapping, canonical
        names require non-alphanumeric boundaries and allow an omitted final
        apostrophe. For example, SDSS_g matches SDSS_g_001, not SDSS_green_001.
        Raise ValueError for no matches, conflicting matches, or unknown aliases.
        """
        if mapping is not None:
            candidates: set[str] = {
                self._telescope.resolve_filter_alias(target)
                for token, target in mapping.items() if token in filename
            }

        else:
            # Do not infer physical passbands from a single letter appearing
            # anywhere in a filename.
            candidates = set()
            for standard in self._telescope.filters():
                token: str = standard.rstrip("'")
                pattern: str = (
                    r"(?<![A-Za-z0-9])" + re.escape(token) + r"'?(?![A-Za-z0-9])"
                )
                if re.search(pattern, filename):
                    candidates.add(standard)

        if len(candidates) != 1:
            raise ValueError(
                f"{filename}: expected exactly one filter match, got {candidates!r}."
            )

        return candidates.pop()

    def _apply(
        self,
        source: ImageInput,
        plans: list[HeaderPlan],
        save_location: str | Path | None,
        overwrite: bool | None,
    ) -> ImageInput:
        """Validate destinations, write planned cards, and rebuild collections.

        Parameters
        ----------
        source : ImageFileCollection or ImageCollectionGroup
            Collection that owns the planned files and determines result form.
        plans : list of HeaderPlan
            Ordered (owning collection, absolute input path, header updates)
            tuples. Empty updates copy a file for separate output or leave an
            in-place file untouched.
        save_location : str, Path, or None
            Output override. None inherits the constructor setting; an empty
            string selects in-place output.
        overwrite : bool or None
            Overwrite override. None inherits the constructor setting.

        Returns
        -------
        ImageFileCollection or ImageCollectionGroup
            Written file selections, also assigned to image_collection.
            Original summary keywords and newly written card names are included.

        Notes
        -----
        Validate all destination paths, header values, and per-child directory
        boundaries before writing. Each file is published atomically, but a
        later I/O failure does not roll back earlier files in the same batch.
        Collection summaries are rebuilt only after writing completes.
        """
        if not plans:
            raise ValueError("No files selected for header correction.")

        location: str | Path | None = (
            self._save_location if save_location is None else save_location
        )
        output_dir: Path | None = (
            None if location in (None, '') else Path(location).expanduser().resolve()
        )
        replace: bool = self._overwrite if overwrite is None else overwrite
        destinations: list[Path] = []
        seen: set[Path] = set()

        # Preflight every destination and card before the first write.
        for child, path, cards in plans:
            destination: Path = path if output_dir is None else output_dir / path.name
            if destination in seen:
                raise ValueError(f"Duplicate output path: {destination}.")

            seen.add(destination)
            changes_file: bool = bool(cards) or destination != path
            if changes_file and destination.exists() and not replace:
                raise FileExistsError(destination)

            header: fits.Header = fits.getheader(path, ext=child.ext).copy()
            for key, value in cards.items():
                header[key] = value

            destinations.append(destination)

        # Each reconstructed child must still represent one output directory.
        child_paths: list[tuple[ImageFileCollection, Path, list[Path]]] = []
        for child in self._children(source):
            paths: list[Path] = [
                destination for (owner, _, _), destination in zip(plans, destinations)
                if owner is child
            ]
            if paths:
                parent: Path = paths[0].parent
                if any(path.parent != parent for path in paths):
                    raise ValueError("Each child must have one output directory.")

                child_paths.append((child, parent, paths))

        # Keep writing separate from planning and validation.
        for (child, path, cards), destination in tqdm(
            zip(plans, destinations),
            total=len(plans),
            desc="Updating FITS headers",
            unit="file",
        ):
            if cards or destination != path:
                update_fits_header(
                    path, cards, output_path=destination,
                    hdu_index=child.ext, overwrite=replace,
                )

        # Include new card names so FIELD/IMTYPE/FILTER appear in the summary.
        children: list[ImageFileCollection] = [
            ImageFileCollection(
                location=parent,
                filenames=[path.name for path in paths],
                keywords=list(dict.fromkeys([
                    *(str(key).lower() for key in child.keywords),
                    *(key.lower() for owner, _, cards in plans
                      if owner is child for key in cards),
                ])),
                ext=child.ext,
            )
            for child, parent, paths in child_paths
        ]
        result: ImageInput = (
            ImageCollectionGroup.from_collections(children)
            if isinstance(source, ImageCollectionGroup) else children[0]
        )
        self._image_collection = result

        return result
