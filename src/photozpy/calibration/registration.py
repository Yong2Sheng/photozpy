"""Register Light images to a common pixel grid within each observation field.

Originally written by Yong Sheng at Clemson University, 2023, for Photozpy.
Advisor: Dr. Marco Ajello.

The original workflow is retained: select Light images, choose a reference in
one filter, call astroalign.register for the remaining images, and mark ALIGN.
The refactoring updates collection/field APIs, reference tracking, and masks.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
import hashlib
import logging
import re
from pathlib import Path
from typing import Literal
from tempfile import TemporaryDirectory

from astroalign import MaxIterError, apply_transform, find_transform, register
from astropy.io import fits
from astropy.nddata import CCDData
import numpy as np
from scipy.ndimage import gaussian_filter
from skimage.transform import warp
from tqdm.auto import tqdm

from ..fits_utilities import has_wcs, update_fits_header
from ..telescope import Telescope
from ._calibration_io import (
    ImageInput, check_destinations, children, rebuild, refresh, select,
    validate_compatible, write_ccd,
)


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class _RegistrationPlan:
    """One field's validated inputs and output paths; not an image-data model."""

    field_name: str
    image_paths: tuple[Path, ...]
    reference_image_path: Path
    reference_hash: str
    destinations: dict[Path, Path]
    discard_input_wcs: bool = False


class Registration:
    """Align reduced Light images, retaining the original astroalign algorithm.

    Parameters
    ----------
    image_collection : ImageFileCollection or ImageCollectionGroup
        Working images. Only IMTYPE='Light' members are transformed. Each child
        is processed separately, with Light frames grouped by FIELD, not OBJECT.
    telescope : Telescope or None, optional
        Refactored telescope for reference-filter alias resolution. Without a
        telescope, FILTER labels must match the requested reference filter.
    overwrite : bool, optional
        Allow replacement of existing destination files. Default is False;
        explicitly enable it when updating analysis working copies in place.

    Notes
    -----
    All filters of a FIELD share one reference grid, as in the original code.
    Prefer reference_images to choose a particular FITS for every FIELD. Its
    keys are FIELD header values, not SourceInfo names. Other images, including
    other filters, are aligned to that field's selected image. There is no
    automatic quality ranking or sensitivity-based filter selection.
    Default source-detection settings remain min_area=5 and detection_sigma=5.
    Astroalign's interpolation, transformation, and flux behavior are unchanged;
    this class does not add a photometric rescaling or claim flux conservation.

    CMASK is used for star detection only and resampled separately for inspection.
    It never becomes the output validity MASK. Footprints and existing validity
    masks are saved as MASK. Unrelated
    FITS extensions, uncertainty arrays, and (by default) pre-existing spatial
    WCS on targets are rejected before writing. discard_input_wcs=True removes
    old WCS on newly resampled outputs for workflows solving Master WCS later.
    An unresampled reference may retain its own valid WCS.

    ALIGN retains the original Reference/Registered values. REGREF and REFHASH
    identify the reference filename and a digest of its pixel grid and mask.
    Already registered targets are reused only when these identifiers agree;
    changing references requires pre-registration working copies. Any REGSMTH
    records additionally require a complete, consistent selected field on every
    call, including the default path. Orphan markers without ALIGN are rejected.
    Headers do not authenticate target arrays if pixels/markers are manually
    replaced or all provenance is stripped.

    Methods return updated collections and retain calibration/master members.
    Each output file is written atomically. By default a later transform failure
    may leave earlier targets completed. The optional smoothing fallback stages
    a whole field before publication and retries from native inputs. Publication
    I/O errors are still not a multi-file transaction. No writes or registrations
    occur in the constructor.
    """

    def __init__(
        self,
        image_collection: ImageInput,
        telescope: Telescope | None = None,
        *,
        overwrite: bool = False,
    ) -> None:
        children(image_collection)
        if telescope is not None and not isinstance(telescope, Telescope):
            raise TypeError('telescope must be a refactored Telescope or None.')

        self._image_collection: ImageInput = image_collection
        self._telescope: Telescope | None = telescope
        self._overwrite: bool = overwrite
        self._reference_paths: tuple[Path, ...] = ()

    @property
    def image_collection(self) -> ImageInput:
        """Return the current selection after the latest successful registration."""
        return self._image_collection

    @property
    def reference_paths(self) -> tuple[Path, ...]:
        """Return one output reference path per field from the latest success."""
        return self._reference_paths

    def register_images(
        self,
        reference_filter: str | None = None,
        min_area: int = 5,
        detection_sigma: float = 5,
        *,
        detection_sigmas: Sequence[float] | None = None,
        fallback_smoothing_sigma: float | None = None,
        reference_images: Mapping[str, str | Path] | Literal['auto'] | None = None,
        save_location: str | Path | None = None,
        overwrite: bool | None = None,
        discard_input_wcs: bool = False,
    ) -> ImageInput:
        """Register each FIELD, returning working images with masks and markers.

        Parameters
        ----------
        reference_filter : str or None, optional
            Optional fallback for fields omitted from reference_images. Resolve
            aliases via telescope when supplied; choose the first candidate by
            filename, preferring a valid existing Reference on repeat calls.
            This is a deterministic shortcut, not a quality recommendation.
            Leave None when explicitly choosing a file for every field.
        min_area : int, optional
            Minimum connected detection-region area, in pixels, for accepting
            a candidate star in each image. Default 5 retains the original API.
            This is neither the number of stars nor a radius/FWHM. Increasing
            it rejects small structures but can lose compact or faint stars;
            decreasing it admits smaller sources and potentially artifacts.
            SEP performs detection filtering, so this is not a count of raw
            pixels individually exceeding detection_sigma times the noise.
        detection_sigma : float, optional
            Positive detection threshold multiplier; default 5, unchanged.
            Astroalign estimates each image's background separately and passes
            detection_sigma * background.globalrms to SEP for detection on the
            background-subtracted image. Lower values admit fainter candidates
            and more noise; higher values may leave too few common stars.
            This is not an integrated stellar S/N or a matching tolerance.
        detection_sigmas : sequence of float or None, optional
            Explicit positive, finite detection thresholds to try in order for
            each target, e.g. [5, 6, 7, 8, 9, 10, 11, 12]. Overrides
            detection_sigma when supplied. Stop at the first successful fit;
            retry only failed star matching, including too few detected stars.
            Exhausting the list raises RuntimeError without writing the target.
            None retains the single detection_sigma behavior. Successful targets
            record the actual threshold in the FITS header REGSIG.
        fallback_smoothing_sigma : float or None, optional
            Opt in with e.g. 2.0 (Gaussian sigma in pixels). First exhaust the
            detection_sigmas on original detection images. If any target still
            fails matching, discard staged results and restart ALL targets in
            that FIELD using smoothed detection copies and the same thresholds.
            Science pixels are never smoothed; the fitted transform resamples
            the original science image once. MASK/CMASK remain separate.
            Results are published only after every target in the field succeeds.
            REGSMTH records 0 for ordinary detection or the smoothing sigma on
            every output, including the unchanged reference. Repeated complete
            fields (including reference-only selections) with compatible
            REGSMTH are reused. Recorded mode consistency is checked even when
            this option is omitted. Old/mixed registrations require restoring
            pre-registration copies. Orphan registration markers without ALIGN
            are rejected on all paths. Default None keeps ordinary matching for
            new inputs and reference-based reuse for pure legacy selections.
        reference_images : mapping of str to path, "auto", or None, optional
            Use 'auto' to select the first selected Light image of each FIELD
            in collection order, separately within each child directory. All
            filters in that field align to this reference; no quality ranking
            or reference_filter restriction is applied in auto mode. The chosen
            filename is logged and available in reference_paths. Existing
            registration must still match it; auto never changes an old grid
            silently if collection order or selection has changed.
            Alternatively, supply a FIELD-to-specific-FITS mapping, for example
            {'science_field': analysis_dir / 'science_i_003.fits'}.
            Pick one selected Light image per field; all other exposures and
            filters in that field align to its grid. The reference's pixels are
            unchanged. A different field may use a different reference filter.
            Paths must exist and belong to the field's current selection, not
            its raw or before_cosmic_ray backups. Relative paths are resolved
            against the process working directory; explicit full paths or Paths
            built from analysis_dir avoid ambiguity. Explicit choices take
            precedence over reference_filter. Unspecified fields require that
            fallback; they are not silently assigned an unrelated reference.
            Repeated field names across directories are resolved separately;
            an explicit reference must belong to the corresponding child.
        save_location : str, Path, or None, optional
            Output directory. None or '' updates working files in place.
            A shared output directory must not cause filename collisions.
        overwrite : bool or None, optional
            Override the constructor's overwrite policy. This does not permit
            silently resampling an image already tied to another reference.

        discard_input_wcs : bool, optional
            Default False rejects WCS-bearing targets before writing. Set True
            to remove old spatial WCS (including alternate solutions and SIP)
            from newly resampled output headers. Original input files are not
            separately edited; in-place registration replaces working copies.
            The unresampled reference and reused targets retain their WCS.
            Solve WCS on the registered/stacked images before sky-coordinate use.

        Returns
        -------
        ImageFileCollection or ImageCollectionGroup
            The same input form, with registered output paths and retained
            non-Light members. Also assigned to image_collection.

        Raises
        ------
        ValueError
            Field/filter/reference selection fails, shapes/units disagree,
            prior registration metadata is missing or incompatible, or an
            unsupported WCS, extension, uncertainty, or invalid setting is found.
        FileExistsError
            An output exists and overwriting has not been enabled.
        FileNotFoundError
            A user-selected reference path is not an existing file.
        RuntimeError
            Astroalign fails or produces invalid output. A failed target is
            never marked Registered or written as a successful result.

        Notes
        -----
        Progress titles stay fixed per FIELD; the current file and reference
        appear at the end as 'Aligning image.fits to reference.fits'. Counts
        exclude the retained reference. Existing matching registration is
        reported as Reusing and is not resampled a second time.
        Star detection runs independently on each image/reference pair; the
        algorithm looks for common patterns, not every star in the reference.
        Detection settings are forwarded unchanged to astroalign. Optional
        detection_sigmas enables ordered per-target retries using only the
        supplied values. Completed matching registrations are still reused.
        The legacy notebook's 2/8 combination is not a universal optimum.
        Registration is a resampling operation, not plate solving.
        """
        return self.Register_images(
            filter=reference_filter,
            min_area=min_area,
            detection_sigma=detection_sigma,
            detection_sigmas=detection_sigmas,
            fallback_smoothing_sigma=fallback_smoothing_sigma,
            reference_images=reference_images,
            save_location=save_location,
            overwrite=overwrite,
            discard_input_wcs=discard_input_wcs,
        )

    def Register_images(
        self,
        filter: str | None = None,
        min_area: int = 5,
        detection_sigma: float = 5,
        *,
        detection_sigmas: Sequence[float] | None = None,
        fallback_smoothing_sigma: float | None = None,
        reference_images: Mapping[str, str | Path] | Literal['auto'] | None = None,
        save_location: str | Path | None = None,
        overwrite: bool | None = None,
        discard_input_wcs: bool = False,
    ) -> ImageInput:
        """Retain the original entry point; see register_images for its contract.

        The old filter argument selects the reference band. The lowercase
        register_images entry point uses the clearer name reference_filter.
        """
        if not isinstance(min_area, int) or isinstance(min_area, bool) or min_area <= 0:
            raise ValueError('min_area must be a positive integer.')
        if detection_sigmas is not None and (
            not isinstance(detection_sigmas, Sequence)
            or isinstance(detection_sigmas, (str, bytes))
            or len(detection_sigmas) == 0
        ):
            raise ValueError('detection_sigmas must be a nonempty sequence of positive finite numbers.')
        candidates = (detection_sigma,) if detection_sigmas is None else tuple(detection_sigmas)
        for sigma in candidates:
            if (
                isinstance(sigma, (bool, np.bool_))
                or not isinstance(sigma, (int, float, np.integer, np.floating))
                or not np.isfinite(sigma) or sigma <= 0
            ):
                name = 'detection_sigma' if detection_sigmas is None else 'detection_sigmas'
                raise ValueError(f'{name} must contain only positive finite numbers.')

        if fallback_smoothing_sigma is not None and (
            isinstance(fallback_smoothing_sigma, (bool, np.bool_))
            or not isinstance(fallback_smoothing_sigma, (int, float, np.integer, np.floating))
            or not np.isfinite(fallback_smoothing_sigma) or fallback_smoothing_sigma <= 0
        ):
            raise ValueError('fallback_smoothing_sigma must be a positive finite number or None.')

        if not isinstance(discard_input_wcs, bool):
            raise TypeError('discard_input_wcs must be a bool.')
        replace: bool = self._overwrite if overwrite is None else overwrite
        plans: list[_RegistrationPlan] = self._plan_groups(
            filter, reference_images, save_location, replace, discard_input_wcs,
        )
        # Recorded field recipes remain authoritative on every entry path,
        # including calls that omit the fallback option. Validate all fields
        # before any publication; a default call cannot bypass this guard.
        reuse_modes = [self._fallback_reuse_mode(plan, fallback_smoothing_sigma)
                       for plan in plans]
        replacements: dict[Path, Path] = {}

        for plan, reuse_mode in zip(plans, reuse_modes):
            logger.info(
                'Registration | FIELD=%s\n%s\nReference: %s (%d Light images)',
                plan.field_name, '-' * 64,
                plan.reference_image_path.name, len(plan.image_paths),
            )
            if fallback_smoothing_sigma is None:
                self._Register_images(
                    plan, min_area=min_area, detection_sigma=detection_sigma,
                    detection_sigmas=candidates, overwrite=replace,
                    **({'smoothing_sigma': reuse_mode} if reuse_mode is not None else {}),
                )
            else:
                self._register_with_fallback(
                    plan, min_area=min_area, candidates=candidates,
                    smoothing_sigma=float(fallback_smoothing_sigma), overwrite=replace,
                )
            replacements.update(plan.destinations)

        self._image_collection = rebuild(self._image_collection, replacements=replacements)
        self._reference_paths = tuple(
            plan.destinations[plan.reference_image_path] for plan in plans
        )
        return self._image_collection

    @staticmethod
    def _matching_failure(exc: Exception) -> bool:
        return isinstance(exc, MaxIterError) or (
            isinstance(exc, ValueError)
            and str(exc).startswith('Reference stars in ')
            and 'minimum value' in str(exc)
        )

    @staticmethod
    def _fallback_reuse_mode(
        plan: _RegistrationPlan, smoothing_sigma: float | None,
    ) -> float | None:
        """Recognize a complete recorded field, including a reference-only one.

        Pure legacy inputs retain default behavior. Once any member carries
        REGSMTH, every selected member must form a complete, consistent field,
        even if the current call does not request fallback. Native additions or
        partial publication require restoration instead of silently mixing modes.
        """
        headers = [fits.getheader(path) for path in plan.image_paths]
        targets = [h for path, h in zip(plan.image_paths, headers)
                   if path != plan.reference_image_path]
        has_recorded_mode = any('REGSMTH' in h for h in headers)
        if not has_recorded_mode:
            if smoothing_sigma is None or not any(h.get('ALIGN') == 'Registered' for h in targets):
                return None
        else:
            complete = all(
                h.get('ALIGN') == ('Reference' if path == plan.reference_image_path else 'Registered')
                for path, h in zip(plan.image_paths, headers)
            )
            values = [h.get('REGSMTH') for h in headers]
            valid_modes = all(
                not isinstance(mode, (bool, np.bool_))
                and isinstance(mode, (int, float, np.integer, np.floating))
                and np.isfinite(mode) and mode >= 0 for mode in values
            )
            if complete and valid_modes and len(set(values)) == 1:
                mode = float(values[0])
                if smoothing_sigma is None or mode in (0., smoothing_sigma):
                    return mode
        raise ValueError(
            f'{plan.field_name}: registration requires unresampled inputs '
            'or a complete field with one compatible REGSMTH mode; restore '
            'pre-registration working copies and rerun cosmic-ray cleaning if needed. '
            'Do not clear ALIGN markers on resampled files.'
        )

    @staticmethod
    def _register_with_fallback(
        plan: _RegistrationPlan, *, min_area: int, candidates: Sequence[float],
        smoothing_sigma: float, overwrite: bool,
    ) -> None:
        mode = Registration._fallback_reuse_mode(plan, smoothing_sigma)
        if mode is not None:
            Registration._Register_images(
                plan, min_area=min_area, detection_sigmas=candidates,
                overwrite=overwrite, smoothing_sigma=mode,
            )
            return

        # Stage a complete field on disk, keeping memory bounded and all input
        # pixels native until a single detection recipe succeeds for the field.
        # A failed match never publishes a partially registered field.
        with TemporaryDirectory(prefix='photozpy-registration-') as temporary:
            for mode in (0., smoothing_sigma):
                directory = Path(temporary) / str(mode)
                staged = replace(plan, destinations={
                    path: directory / path.name for path in plan.image_paths
                })
                try:
                    Registration._Register_images(
                        staged, min_area=min_area, detection_sigmas=candidates,
                        overwrite=False, smoothing_sigma=mode,
                    )
                except (RuntimeError, OSError) as exc:
                    if mode == 0. and Registration._matching_failure(exc.__cause__):
                        logger.info(
                            'FIELD=%s: raw detection thresholds exhausted; restarting ALL '
                            'targets from unresampled inputs with detection smoothing sigma=%s pix.',
                            plan.field_name, smoothing_sigma,
                        )
                        continue
                    phase = 'ordinary' if mode == 0. else 'smoothed'
                    phases = 'ordinary' if mode == 0. else 'ordinary + smoothed'
                    raise RuntimeError(
                        f'FIELD={plan.field_name}: {phase} registration failed; '
                        f'phases_attempted={phases}, detection_sigmas={list(candidates)}, '
                        f'smoothing_sigma={mode} pix. No results published for this field. {exc}'
                    ) from (exc.__cause__ if exc.__cause__ is not None else exc)
                # Complete recipe chosen. Publish each file atomically using
                # the original overwrite policy, including new output folders.
                published = 0
                for path in plan.image_paths:
                    try:
                        update_fits_header(
                            staged.destinations[path], {},
                            output_path=plan.destinations[path], overwrite=overwrite,
                        )
                    except OSError as exc:
                        raise type(exc)(
                            f'FIELD={plan.field_name}: publication failed at '
                            f'{plan.destinations[path]} after {published}/{len(plan.image_paths)} '
                            f'files; REGSMTH={mode}. Earlier outputs may already be replaced. '
                            'Restore pre-registration copies before in-place retry, or rerun '
                            'from unchanged native inputs into a separate output directory. '
                            f'{exc}'
                        ) from exc
                    published += 1
                return

    def _plan_groups(
        self,
        reference_filter: str | None,
        reference_images: Mapping[str, str | Path] | Literal['auto'] | None,
        save_location: str | Path | None,
        overwrite: bool,
        discard_input_wcs: bool = False,
    ) -> list[_RegistrationPlan]:
        """Validate all fields and destinations before the first output write."""
        refresh(self._image_collection)
        automatic = isinstance(reference_images, str) and reference_images == 'auto'
        if reference_filter is not None and not automatic:
            reference_filter = self._filter_name(reference_filter)

        if reference_images is not None and not automatic and not isinstance(reference_images, Mapping):
            raise TypeError("reference_images must be 'auto', a FIELD-to-FITS mapping, or None.")

        explicit: dict[str, Path] = {}
        references = reference_images if isinstance(reference_images, Mapping) else {}
        for field_name, supplied_path in references.items():
            if not isinstance(field_name, str) or not field_name.strip():
                raise ValueError('Reference keys must be nonempty FIELD names.')

            reference_path: Path = Path(supplied_path).expanduser().resolve()
            if not reference_path.is_file():
                raise FileNotFoundError(f'{field_name}: reference FITS does not exist: {reference_path}')

            explicit[field_name] = reference_path

        plans: list[_RegistrationPlan] = []
        all_destinations: list[Path] = []
        write_destinations: list[Path] = []
        seen_fields: set[str] = set()
        output_dir = (
            None if save_location in (None, '')
            else Path(save_location).expanduser().resolve()
        )

        for child in children(self._image_collection):
            light_paths: list[Path] = select(child, IMTYPE='Light')
            groups: dict[str, list[Path]] = {}

            for path in light_paths:
                header = fits.getheader(path)
                field_name = str(header.get('FIELD', '')).strip()
                if not field_name:
                    raise ValueError(f'{path.name}: missing FIELD; correct headers first.')
                groups.setdefault(field_name, []).append(path)

            for field_name, image_paths in groups.items():
                seen_fields.add(field_name)
                validate_compatible(image_paths)
                headers: dict[Path, fits.Header] = {
                    path: fits.getheader(path) for path in image_paths
                }

                for path, header in headers.items():
                    if not header.get('FILTER'):
                        raise ValueError(f'{path.name}: missing FILTER.')
                    self._filter_name(str(header['FILTER']))
                    self._validate_extensions(path)
                    align = header.get('ALIGN')
                    if align not in (None, '', 'Reference', 'Registered'):
                        raise ValueError(f'{path.name}: unsupported ALIGN state {align!r}.')
                    if align in (None, '') and any(
                        key in header for key in ('REGSMTH', 'REGSIG', 'REGREF', 'REFHASH')
                    ):
                        raise ValueError(
                            f'{path.name}: registration markers remain without ALIGN; '
                            'restore pre-registration working copies instead of resampling '
                            'an array with unknown registration history.'
                        )

                # Auto follows collection order; explicit/filter choices retain
                # their existing behavior and all reference-identity checks.
                if automatic:
                    reference_image_path = image_paths[0]
                elif field_name in explicit:
                    reference_image_path = explicit[field_name]
                    if reference_image_path not in image_paths:
                        raise ValueError(
                            f'{field_name}: reference {reference_image_path.name} '
                            'must be a selected Light image in this field/directory.'
                        )
                else:
                    if reference_filter is None:
                        raise ValueError(
                            f'{field_name}: choose a FITS in reference_images '
                            'or supply reference_filter as a fallback.'
                        )
                    candidates = [
                        path for path in image_paths
                        if self._filter_name(str(headers[path]['FILTER'])) == reference_filter
                    ]
                    existing = [path for path in candidates if headers[path].get('ALIGN') == 'Reference']
                    if len(existing) > 1:
                        raise ValueError(f'{field_name}: multiple existing reference images.')
                    if not candidates:
                        raise ValueError(f'{field_name}: no reference image in {reference_filter}.')
                    reference_image_path = existing[0] if existing else sorted(candidates)[0]

                if headers[reference_image_path].get('ALIGN') == 'Registered':
                    raise ValueError('A resampled target cannot silently become a new reference.')

                reference_ccd = self._read_ccd(reference_image_path)
                reference_hash = self._fingerprint(
                    reference_ccd, self._read_cmask(reference_image_path),
                )
                destinations: dict[Path, Path] = {}

                for path in image_paths:
                    header = headers[path]
                    if header.get('ALIGN') in ('Reference', 'Registered'):
                        if not self._matches_reference(header, reference_image_path, reference_hash):
                            raise ValueError(
                                f'{path.name}: registration reference is missing or changed; '
                                'restore pre-registration working copies.'
                            )
                        if header.get('ALIGN') == 'Reference' and path != reference_image_path:
                            raise ValueError('Changing the reference requires pre-registration copies.')

                    if path != reference_image_path and not self._matches_reference(
                        header, reference_image_path, reference_hash,
                    ):
                        # A transformed array cannot retain an old spatial WCS.
                        # This first pipeline stage expects solving after stacking.
                        if not discard_input_wcs and (has_wcs(path) or any(
                            key.startswith(('CTYPE', 'CRPIX', 'CRVAL', 'CDELT', 'CROTA'))
                            or key in ('CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
                                       'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2')
                            for key in header
                        )):
                            raise ValueError(
                                f'{path.name}: spatial WCS needs explicit transformation; '
                                'use pre-WCS images or set discard_input_wcs=True and solve '
                                'WCS again after registration/stacking.'
                            )

                    destination = path if output_dir is None else output_dir / path.name
                    destinations[path] = destination
                    all_destinations.append(destination)
                    if destination != path or not self._matches_reference(
                        header, reference_image_path, reference_hash,
                    ):
                        write_destinations.append(destination)

                plans.append(_RegistrationPlan(
                    field_name=field_name,
                    image_paths=tuple(image_paths),
                    reference_image_path=reference_image_path,
                    reference_hash=reference_hash,
                    destinations=destinations,
                    discard_input_wcs=discard_input_wcs,
                ))

        if not plans:
            raise ValueError('No Light images selected for registration.')
        if explicit.keys() - seen_fields:
            raise ValueError('Explicit references contain unselected FIELD names.')

        # Check all paths, including reused members, for cross-field collisions.
        if len(all_destinations) != len(set(all_destinations)):
            raise ValueError('Multiple registration inputs map to the same output path.')
        check_destinations(write_destinations, overwrite)
        return plans

    @staticmethod
    def _without_spatial_wcs(header: fits.Header) -> fits.Header:
        """Remove primary/alternate spatial WCS and SIP from a resampled output.

        Keep observation metadata (e.g. RA/DEC pointing), photometric headers,
        and the original input file. The output deliberately needs plate solving.
        """
        result = header.copy()
        pattern = re.compile(
            r'(?:(?:WCSAXES|WCSNAME|LONPOLE|LATPOLE|RADESYS|RADECSYS|EQUINOX)[A-Z]?'
            r'|(?:CTYPE|CUNIT|CRPIX|CRVAL|CDELT|CROTA|CNAME|CRDER|CSYER)\d+[A-Z]?'
            r'|(?:PC|CD|PV|PS)\d+_\d+[A-Z]?'
            r'|(?:A|B|AP|BP)_(?:ORDER|DMAX|\d+_\d+)'
            r'|(?:CPDIS|CQDIS|D2IMDIS)\d+|(?:DP|DQ|D2IM)\d+\..+)'
        )
        removed = False
        for key in list(result):
            if pattern.fullmatch(key):
                del result[key]
                removed = True
        if removed:
            result.add_history('Registration removed input spatial WCS after resampling; solve WCS again.')
        return result

    def _filter_name(self, name: str) -> str:
        """Resolve the reference band without depending on legacy Telescope."""
        name = name.strip()
        return self._telescope.resolve_filter_alias(name) if self._telescope else name

    @staticmethod
    def _validate_extensions(path: Path) -> None:
        """Accept the primary image, validity MASK, and detection-only CMASK.

        Uncertainty and unrelated arrays cannot be copied unchanged across a
        coordinate transformation. Reject them instead of silently discarding
        or claiming to propagate them correctly.
        """
        with fits.open(path, memmap=True) as hdul:
            unsupported = [hdu.name for hdu in hdul[1:] if hdu.name not in ('MASK', 'CMASK')]
            if hdul[0].header.get('CRCORR') == 'Yes' and 'MASK' in hdul and 'CMASK' not in hdul:
                raise ValueError(
                    f'{path.name}: legacy cosmic-ray MASK; restore pre-registration '
                    'copies and separate cosmic-ray detections into CMASK first.'
                )
            if 'CMASK' in hdul:
                cmask = hdul['CMASK'].data
                if cmask is None or cmask.shape != hdul[0].data.shape or not np.all(np.isin(cmask, (0, 1))):
                    raise ValueError(f'{path.name}: CMASK must be binary and match the image shape.')
        if unsupported:
            raise ValueError(f'{path.name}: unsupported registration extensions {unsupported}.')

    @staticmethod
    def _read_ccd(path: Path) -> CCDData:
        """Read native-endian image data without altering its numerical values."""
        ccd: CCDData = CCDData.read(path)
        if ccd.uncertainty is not None:
            raise ValueError(f'{path.name}: uncertainty resampling is not implemented.')
        ccd.data = np.asarray(ccd.data, dtype=ccd.data.dtype.newbyteorder('='))
        if not np.all(np.isfinite(ccd.data)):
            raise ValueError(f'{path.name}: registration requires finite image data.')
        return ccd

    @staticmethod
    def _read_cmask(path: Path) -> np.ndarray | None:
        """Read detection history explicitly; CCDData.mask remains validity-only."""
        with fits.open(path, memmap=True) as hdul:
            return np.asarray(hdul['CMASK'].data, dtype=bool).copy() if 'CMASK' in hdul else None

    @staticmethod
    def _fingerprint(ccd: CCDData, cmask: np.ndarray | None = None) -> str:
        """Identify reference pixels, validity, and optional detection history."""
        digest = hashlib.sha256()
        data = np.ascontiguousarray(ccd.data)
        digest.update(str((data.shape, data.dtype.str)).encode())
        digest.update(data.tobytes())
        mask = np.zeros(data.shape, dtype=bool) if ccd.mask is None else ccd.mask
        digest.update(np.ascontiguousarray(mask, dtype=bool).tobytes())
        if cmask is not None:
            digest.update(b'CMASK')
            digest.update(np.ascontiguousarray(cmask, dtype=bool).tobytes())
        return digest.hexdigest()

    @staticmethod
    def _align_with_cmask(
        source: CCDData, reference: CCDData,
        source_cmask: np.ndarray | None, reference_cmask: np.ndarray | None,
        *, min_area: int, detection_sigma: float, smoothing_sigma: float = 0.,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
        """Use cosmic detections for star finding, never for output validity.

        Split astroalign.register into its original find/apply calls to reuse
        the same transform for the independent CMASK. Science interpolation
        stays unchanged; nearest-neighbor CMASK resampling retains binary labels.
        """
        detection_images = []
        for image, cmask in ((source, source_cmask), (reference, reference_cmask)):
            detection = image.copy()
            if cmask is not None:
                detection.mask = cmask if image.mask is None else image.mask | cmask
            if smoothing_sigma > 0:
                # Normalized convolution prevents invalid pixels/cosmic detections
                # from bleeding into neighboring stellar detection pixels.
                valid = (np.ones(image.shape, dtype=bool) if detection.mask is None
                         else ~detection.mask)
                weights = gaussian_filter(valid.astype(float), smoothing_sigma)
                values = gaussian_filter(np.where(valid, image.data, 0.).astype(float), smoothing_sigma)
                detection.data = np.divide(
                    values, weights, out=np.full(image.shape, np.median(image.data)),
                    where=weights > 0,
                )
            detection_images.append(detection)

        transform, _ = find_transform(
            *detection_images, min_area=min_area, detection_sigma=detection_sigma,
        )
        aligned, footprint = apply_transform(
            transform, source, reference, propagate_mask=True,
        )
        aligned_cmask = None
        if source_cmask is not None:
            aligned_cmask = warp(
                source_cmask.astype(np.uint8), inverse_map=transform.inverse,
                output_shape=reference.shape, order=0, mode='constant',
                cval=0, preserve_range=True,
            ).astype(bool)
        return aligned, footprint, aligned_cmask

    @staticmethod
    def _matches_reference(header: fits.Header, reference: Path, digest: str) -> bool:
        """Recognize only registration markers tied to this reference version."""
        return (
            header.get('ALIGN') in ('Reference', 'Registered')
            and header.get('REGREF') == reference.name
            and header.get('REFHASH') == digest
        )

    @staticmethod
    def _Register_images(
        plan: _RegistrationPlan,
        min_area: int = 5,
        detection_sigma: float = 5,
        *,
        overwrite: bool,
        detection_sigmas: Sequence[float] | None = None,
        smoothing_sigma: float | None = None,
    ) -> None:
        """Register one FIELD, keeping the original per-target astroalign call."""
        reference_image_path: Path = plan.reference_image_path
        reference_data: CCDData = Registration._read_ccd(reference_image_path)
        reference_cmask = Registration._read_cmask(reference_image_path)
        if Registration._fingerprint(reference_data, reference_cmask) != plan.reference_hash:
            raise ValueError('Reference pixels changed after registration planning.')

        # Mark the unresampled reference without rewriting its image arrays.
        reference_header = fits.getheader(reference_image_path)
        reference_destination = plan.destinations[reference_image_path]
        markers = {'REGREF': reference_image_path.name, 'REFHASH': plan.reference_hash}
        if smoothing_sigma is not None:
            markers['REGSMTH'] = (smoothing_sigma, 'Detection-only Gaussian sigma (pix); 0=none')
        if reference_destination != reference_image_path or not Registration._matches_reference(
            reference_header, reference_image_path, plan.reference_hash,
        ):
            update_fits_header(
                reference_image_path, {'ALIGN': 'Reference', **markers},
                output_path=reference_destination, overwrite=overwrite,
            )

        register_image_paths: list[Path] = [
            path for path in plan.image_paths if path != reference_image_path
        ]
        registered_count: int = 0
        reused_count: int = 0

        with tqdm(
            register_image_paths, desc=f'Registration | {plan.field_name}',
            unit='image', disable=not register_image_paths,
        ) as progress:
            for path in progress:
                header_target: fits.Header = fits.getheader(path)
                destination: Path = plan.destinations[path]

                if Registration._matches_reference(
                    header_target, reference_image_path, plan.reference_hash,
                ):
                    progress.set_postfix_str(
                        f'Reusing {path.name} (reference: {reference_image_path.name})'
                    )
                    if destination != path:
                        update_fits_header(path, markers, output_path=destination, overwrite=overwrite)
                    reused_count += 1
                    continue

                progress.set_postfix_str(
                    f'Aligning {path.name} to {reference_image_path.name}'
                )
                target_data: CCDData = Registration._read_ccd(path)
                target_cmask = Registration._read_cmask(path)

                # Each attempt uses the original, unresampled input. Only
                # matching failures can benefit from another detection threshold.
                candidates = (detection_sigma,) if detection_sigmas is None else detection_sigmas
                attempted = []
                for sigma in candidates:
                    attempted.append(sigma)
                    try:
                        aligned_cmask = None
                        if target_cmask is not None or reference_cmask is not None or smoothing_sigma:
                            smoothing_options = ({'smoothing_sigma': smoothing_sigma}
                                                 if smoothing_sigma else {})
                            img_aligned, footprint, aligned_cmask = Registration._align_with_cmask(
                                target_data, reference_data, target_cmask, reference_cmask,
                                min_area=min_area, detection_sigma=sigma, **smoothing_options,
                            )
                        else:
                            img_aligned, footprint = register(
                                target_data, reference_data,
                                min_area=min_area, detection_sigma=sigma,
                                propagate_mask=True,
                            )
                    except Exception as exc:
                        matching_failure = Registration._matching_failure(exc)
                        if matching_failure and len(attempted) < len(candidates):
                            logger.info(
                                '%s: no match at detection_sigma=%s; trying detection_sigma=%s.',
                                path.name, sigma, candidates[len(attempted)],
                            )
                            continue
                        raise RuntimeError(
                            f'{path.name}: registration to {reference_image_path.name} failed '
                            f'(attempted detection_sigmas={attempted}, min_area={min_area}).'
                        ) from exc
                    break
                if len(attempted) > 1:
                    logger.info('%s: matched at detection_sigma=%s.', path.name, sigma)

                if img_aligned.shape != reference_data.shape or footprint.shape != reference_data.shape:
                    raise RuntimeError(f'{path.name}: registration returned an invalid output shape.')
                mask = np.asarray(footprint, dtype=bool) | ~np.isfinite(img_aligned)
                if np.all(mask):
                    raise RuntimeError(f'{path.name}: registration has no valid coverage.')

                if plan.discard_input_wcs:
                    header_target = Registration._without_spatial_wcs(header_target)

                if smoothing_sigma is not None:
                    header_target['REGSMTH'] = markers['REGSMTH']
                header_target['REGSIG'] = (sigma, 'Detection sigma used for registration')
                header_target['ALIGN'] = 'Registered'
                header_target['REGREF'] = reference_image_path.name
                header_target['REFHASH'] = plan.reference_hash
                aligned_ccd = CCDData(
                    img_aligned, unit=target_data.unit, meta=header_target,
                    mask=mask,
                )
                write_ccd(aligned_ccd, destination, overwrite=overwrite, cmask=aligned_cmask)
                registered_count += 1

        logger.info(
            'FIELD=%s: registered %d images; reused %d; reference retained.\n',
            plan.field_name, registered_count, reused_count,
        )
