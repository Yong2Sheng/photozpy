"""Register Light images to a common pixel grid within each observation field.

Originally written by Yong Sheng at Clemson University, 2023, for Photozpy.
Advisor: Dr. Marco Ajello.

The original workflow is retained: select Light images, choose a reference in
one filter, call astroalign.register for the remaining images, and mark ALIGN.
The refactoring updates collection/field APIs, reference tracking, and masks.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import hashlib
import logging
from pathlib import Path

from astroalign import apply_transform, find_transform, register
from astropy.io import fits
from astropy.nddata import CCDData
import numpy as np
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
    FITS extensions, uncertainty arrays, and pre-existing spatial WCS on images
    to be resampled are rejected before writing: they need their own explicit
    transformation model. For this pipeline, solve WCS on Master Light later.
    An unresampled reference may retain its own valid WCS.

    ALIGN retains the original Reference/Registered values. REGREF and REFHASH
    identify the reference filename and a digest of its pixel grid and mask.
    Already registered targets are reused only when these identifiers agree;
    changing references requires pre-registration working copies.

    Methods return updated collections and retain calibration/master members.
    Each output file is written atomically; the batch is not a transaction if a
    later transform fails. No writes or registrations occur in the constructor.
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
        reference_images: Mapping[str, str | Path] | None = None,
        save_location: str | Path | None = None,
        overwrite: bool | None = None,
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
        reference_images : mapping of str to path, or None, optional
            FIELD-to-specific-FITS mapping, for example
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
        Both detection settings are forwarded unchanged to astroalign; neither
        is adjusted automatically after failure. The legacy notebook's 2/8
        combination is an explicit starting choice, not a universal optimum.
        Registration is a resampling operation, not plate solving.
        """
        return self.Register_images(
            filter=reference_filter,
            min_area=min_area,
            detection_sigma=detection_sigma,
            reference_images=reference_images,
            save_location=save_location,
            overwrite=overwrite,
        )

    def Register_images(
        self,
        filter: str | None = None,
        min_area: int = 5,
        detection_sigma: float = 5,
        *,
        reference_images: Mapping[str, str | Path] | None = None,
        save_location: str | Path | None = None,
        overwrite: bool | None = None,
    ) -> ImageInput:
        """Retain the original entry point; see register_images for its contract.

        The old filter argument selects the reference band. The lowercase
        register_images entry point uses the clearer name reference_filter.
        """
        if not isinstance(min_area, int) or isinstance(min_area, bool) or min_area <= 0:
            raise ValueError('min_area must be a positive integer.')
        if not np.isfinite(detection_sigma) or detection_sigma <= 0:
            raise ValueError('detection_sigma must be positive and finite.')

        replace: bool = self._overwrite if overwrite is None else overwrite
        plans: list[_RegistrationPlan] = self._plan_groups(
            filter, reference_images, save_location, replace,
        )
        replacements: dict[Path, Path] = {}

        for plan in plans:
            logger.info(
                'Registration | FIELD=%s\n%s\nReference: %s (%d Light images)',
                plan.field_name, '-' * 64,
                plan.reference_image_path.name, len(plan.image_paths),
            )
            self._Register_images(
                plan, min_area=min_area, detection_sigma=detection_sigma,
                overwrite=replace,
            )
            replacements.update(plan.destinations)

        self._image_collection = rebuild(self._image_collection, replacements=replacements)
        self._reference_paths = tuple(
            plan.destinations[plan.reference_image_path] for plan in plans
        )
        return self._image_collection

    def _plan_groups(
        self,
        reference_filter: str | None,
        reference_images: Mapping[str, str | Path] | None,
        save_location: str | Path | None,
        overwrite: bool,
    ) -> list[_RegistrationPlan]:
        """Validate all fields and destinations before the first output write."""
        refresh(self._image_collection)
        if reference_filter is not None:
            reference_filter = self._filter_name(reference_filter)

        if reference_images is not None and not isinstance(reference_images, Mapping):
            raise TypeError('reference_images must map FIELD names to FITS paths.')

        explicit: dict[str, Path] = {}
        for field_name, supplied_path in (reference_images or {}).items():
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

                # Retain the original reference-filter strategy, but make the
                # default order deterministic and respect an existing reference.
                if field_name in explicit:
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
                        if has_wcs(path) or any(
                            key.startswith(('CTYPE', 'CRPIX', 'CRVAL', 'CDELT', 'CROTA'))
                            or key in ('CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
                                       'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2')
                            for key in header
                        ):
                            raise ValueError(
                                f'{path.name}: spatial WCS needs explicit transformation; '
                                'use pre-WCS images for this registration stage.'
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
        *, min_area: int, detection_sigma: float,
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

                # Original star-pattern registration, including its default
                # interpolation and source-detection parameters. CCDData allows
                # astroalign to propagate validity masks onto the footprint.
                # Cosmic detections affect star finding only, not pixel validity.
                # Each image supplies its own background/noise estimate. These
                # are detection controls, not a list of required reference stars
                # or a position-matching tolerance; do not tune them silently.
                try:
                    aligned_cmask = None
                    if target_cmask is not None or reference_cmask is not None:
                        img_aligned, footprint, aligned_cmask = Registration._align_with_cmask(
                            target_data, reference_data, target_cmask, reference_cmask,
                            min_area=min_area, detection_sigma=detection_sigma,
                        )
                    else:
                        img_aligned, footprint = register(
                            target_data, reference_data,
                            min_area=min_area, detection_sigma=detection_sigma,
                            propagate_mask=True,
                        )
                except Exception as exc:
                    raise RuntimeError(
                        f'{path.name}: registration to {reference_image_path.name} failed '
                        f'(detection_sigma={detection_sigma}, min_area={min_area}).'
                    ) from exc

                if img_aligned.shape != reference_data.shape or footprint.shape != reference_data.shape:
                    raise RuntimeError(f'{path.name}: registration returned an invalid output shape.')
                mask = np.asarray(footprint, dtype=bool) | ~np.isfinite(img_aligned)
                if np.all(mask):
                    raise RuntimeError(f'{path.name}: registration has no valid coverage.')

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
