"""Build master calibration frames using the existing ccdproc combine recipe."""

from __future__ import annotations

from collections.abc import Callable
import logging
from numbers import Real
from pathlib import Path
import re
from textwrap import fill
from typing import Literal, TypeAlias

from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import mad_std
from ccdproc import ImageFileCollection, combine as ccdpro_combine
import numpy as np
from numpy.typing import NDArray

from ..telescope import Telescope
from ._calibration_io import (
    ImageInput, check_destinations, children, exposure, is_corrected, paths,
    rebuild, refresh, select, validate_compatible, write_ccd,
)


logger = logging.getLogger(__name__)


ScaleFunction: TypeAlias = Callable[[NDArray], float | np.floating]
CombinePlan: TypeAlias = tuple[
    ImageFileCollection, list[Path], Path, str | None, str | None,
]


class Combine:
    """Combine calibration frames separately within each child collection.

    Parameters
    ----------
    image_collection : ImageFileCollection or ImageCollectionGroup
        Working images with standardized IMTYPE, FILTER, and BUNIT headers.
        Master files are selected by metadata, not by filename or OBJECT.
    telescope : Telescope or None, optional
        Refactored telescope used to resolve filter aliases. If omitted, FILTER
        values must already use exactly matching canonical names.
    overwrite : bool, optional
        Allow replacement of existing master files. Default is False.

    Notes
    -----
    The original combination recipe is retained: sigma-clipped average,
    np.ma.median center, mad_std deviation, default 5-sigma thresholds, and
    float32 calculation. Flat inputs are scaled by inverse median by default,
    before sigma clipping so differing sky levels do not hide outliers.
    A single input is copied without this multi-image scaling, as before.

    Bias/Dark grouping uses IMTYPE, Flat uses FILTER, and registered Light uses
    FIELD plus FILTER. Directory boundaries are never merged automatically.
    Mixed correction states and different dark exposure times are rejected.

    Methods return explicit selections including original images and new
    masters, and update image_collection. Passing that result to Reduction
    avoids relying on a snapshot refresh to discover newly created files.
    Output is written as CCDData, not as a copy of unrelated FITS extensions.
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
        self._master_paths: tuple[Path, ...] = ()

    @property
    def image_collection(self) -> ImageInput:
        """Return working images plus masters from the latest successful call."""
        return self._image_collection

    @property
    def master_paths(self) -> tuple[Path, ...]:
        """Return master files written by the latest successful combine call."""
        return self._master_paths

    def preview_in_ds9(self) -> None:
        """Browse masters from the latest successful call without moving files.

        Start or reuse a dedicated DS9 window. Next/Previous browse one master
        at a time; Quit closes the controls and leaves DS9 open. The notebook
        waits for Quit. On macOS the controls use an isolated GUI process.
        This works for Bias, Dark, Flat, and Light masters and does not repeat
        combination or refresh/change the working collection.
        """
        from ..interactive_ds9.image_previewer import ImagePreviewer

        if not self._master_paths:
            raise ValueError('No completed combine batch is available for preview.')
        ImagePreviewer(self._master_paths).run()

    @staticmethod
    def inv_median(a: NDArray) -> float | np.floating:
        """Return the original inverse-median flat scale, rejecting invalid data."""
        median = np.median(a)
        if not np.isfinite(median) or median <= 0:
            raise ValueError('Flat median must be positive and finite.')
        return 1 / median

    @staticmethod
    def check_imtype(image_collection: ImageInput, imtype: str) -> bool:
        """Check that every selected image has the requested IMTYPE.

        Return False for empty selections. Unlike the old character-by-character
        membership check, this checks the complete type string on every image.
        """
        image_paths = [path for child in children(image_collection) for path in paths(child)]
        return bool(image_paths) and all(
            fits.getheader(path).get('IMTYPE') == imtype for path in image_paths
        )

    def combine_bias_or_dark(
        self,
        image_type: Literal['Bias', 'Dark'],
        method: str = 'sigma clip',
        sigma_clip_low: float = 5,
        sigma_clip_high: float = 5,
        save_location: str | Path | None = None,
        quite: bool = False,
        *,
        overwrite: bool | None = None,
    ) -> ImageInput:
        """Create one Master Bias or Master Dark per input child directory.

        Parameters
        ----------
        image_type : {'Bias', 'Dark'}
            Raw image type to combine. Existing masters are not input frames.
        method : str, optional
            Only 'sigma clip' is supported, matching the implemented old recipe.
        sigma_clip_low, sigma_clip_high : float, optional
            Positive lower and upper clipping thresholds. Defaults are both 5.
        save_location : str, Path, or None, optional
            Destination directory; None or '' writes beside each child input.
        quite : bool, optional
            Suppress INFO output. The existing parameter name is retained.
        overwrite : bool or None, optional
            Override the constructor's overwrite policy for this call.

        Returns
        -------
        ImageFileCollection or ImageCollectionGroup
            Original working images and newly generated master files.

        Raises
        ------
        ValueError
            No inputs exist, shapes/units or correction states differ, dark
            exposure times differ, or destination paths collide.
        FileExistsError
            A master exists and overwriting has not been enabled.

        Notes
        -----
        Dark frames must have one common positive EXPTIME. This method does
        not scale individual dark inputs before combining or subtract bias.
        Apply bias correction before building a bias-subtracted Master Dark.
        """
        if image_type not in ('Bias', 'Dark'):
            raise ValueError('image_type must be Bias or Dark.')

        self._validate_method(method, sigma_clip_low, sigma_clip_high)
        refresh(self._image_collection)
        plans: list[CombinePlan] = []

        for child in children(self._image_collection):
            image_list_to_combine: list[Path] = select(child, IMTYPE=image_type)
            if not image_list_to_combine:
                continue

            self._validate_inputs(image_list_to_combine)
            if image_type == 'Dark':
                exposures = [exposure(fits.getheader(path), path) for path in image_list_to_combine]
                if len(set(exposures)) != 1:
                    raise ValueError('Combine dark frames with one common EXPTIME.')

            destination = self._output_dir(child, save_location) / f'Master_{image_type}.fits'
            plans.append((child, image_list_to_combine, destination, None, None))

        return self._execute(
            plans, image_type, sigma_clip_low, sigma_clip_high,
            scale_function=None, quite=quite, overwrite=overwrite,
        )

    def combine_light_or_flat(
        self,
        image_type: Literal['Flat', 'Light'],
        method: str = 'sigma clip',
        scale_function: Literal['auto'] | ScaleFunction | None = 'auto',
        sigma_clip_low: float = 5,
        sigma_clip_high: float = 5,
        save_location: str | Path | None = None,
        quite: bool = False,
        *,
        overwrite: bool | None = None,
    ) -> ImageInput:
        """Combine Flat by FILTER or already registered Light by FIELD/FILTER.

        Parameters
        ----------
        image_type : {'Flat', 'Light'}
            Type to combine. This method does not register light images.
        method : str, optional
            Only 'sigma clip' is supported.
        scale_function : 'auto', callable, or None, optional
            For Flat, 'auto' retains inverse-median scaling; a callable supplies
            an explicit positive finite scale, and None disables input scaling.
            Scaling is applied to each image before clipping. Single-image
            groups are copied without scaling, as in the old implementation.
            Light combinations are unscaled, preserving the original behavior.
        sigma_clip_low, sigma_clip_high : float, optional
            Positive clipping thresholds; defaults are both 5.
        save_location : str, Path, or None, optional
            Output directory; None or '' writes beside each child input.
        quite : bool, optional
            Suppress INFO output; existing parameter spelling retained.
        overwrite : bool or None, optional
            Override the constructor's overwrite policy.

        Returns
        -------
        ImageFileCollection or ImageCollectionGroup
            Working images with new masters included explicitly.

        Notes
        -----
        Multi-image Flat/Light products retain float32 conversion and removal
        of uncertainty. Flat masks are removed as before; Light masks are kept
        so registration footprints remain excluded and uncovered output pixels
        stay identifiable. No new uncertainty model is added.
        Flat correction later uses ccdproc's mean normalization independently
        of this per-input stacking scale. OBJECT is not used for grouping.

        Multi-image flats are normalized in memory one filter group at a time
        before being passed to ccdproc. Unlike passing scale= to ccdproc 2.5.1,
        this makes rejection operate on normalized pixels. No normalized
        intermediate FITS files are written or raw inputs modified.
        """
        if image_type not in ('Flat', 'Light'):
            raise ValueError('image_type must be Flat or Light.')
        if scale_function != 'auto' and scale_function is not None and not callable(scale_function):
            raise TypeError("scale_function must be 'auto', callable, or None.")

        self._validate_method(method, sigma_clip_low, sigma_clip_high)
        refresh(self._image_collection)
        plans: list[CombinePlan] = []

        for child in children(self._image_collection):
            selected = select(child, IMTYPE=image_type)
            groups: dict[tuple[str | None, str], list[Path]] = {}

            for path in selected:
                header = fits.getheader(path)
                if not header.get('FILTER'):
                    raise ValueError(f'{path.name}: missing FILTER.')
                filter_name = str(header['FILTER']).strip()
                if self._telescope is not None:
                    filter_name = self._telescope.resolve_filter_alias(filter_name)

                field_name = None
                if image_type == 'Light':
                    if not header.get('FIELD'):
                        raise ValueError(f'{path.name}: missing FIELD.')
                    field_name = str(header['FIELD'])

                groups.setdefault((field_name, filter_name), []).append(path)

            for (field_name, filter_name), image_list_to_combine in groups.items():
                self._validate_inputs(image_list_to_combine)
                if image_type == 'Light' and len(image_list_to_combine) > 1:
                    self._validate_registration(image_list_to_combine)
                label = 'Flat' if field_name is None else self._safe_name(field_name)
                destination = (
                    self._output_dir(child, save_location)
                    / f'Master_{label}_{self._safe_name(filter_name)}.fits'
                )
                plans.append((child, image_list_to_combine, destination, field_name, filter_name))

        scale_function_ = self.inv_median if scale_function == 'auto' else scale_function
        return self._execute(
            plans, image_type, sigma_clip_low, sigma_clip_high,
            scale_function=scale_function_, quite=quite, overwrite=overwrite,
        )

    @staticmethod
    def _validate_method(method: str, low: float, high: float) -> None:
        """Reject unimplemented methods and invalid thresholds before writing."""
        if method != 'sigma clip':
            raise ValueError("Only method='sigma clip' is implemented.")
        if not np.isfinite(low) or not np.isfinite(high) or low <= 0 or high <= 0:
            raise ValueError('Sigma-clipping thresholds must be positive and finite.')

    @staticmethod
    def _validate_inputs(image_paths: list[Path]) -> None:
        """Check geometry, units, and consistent prior calibration states."""
        validate_compatible(image_paths)
        for keyword in ('BIASCORR', 'DARKCORR', 'FLATCORR'):
            states = {is_corrected(fits.getheader(path), keyword) for path in image_paths}
            if len(states) > 1:
                raise ValueError(f'Cannot combine mixed {keyword} states.')

    @staticmethod
    def _validate_registration(image_paths: list[Path]) -> None:
        """Require a shared, explicitly tracked grid for multi-Light stacking.

        Registration outputs carry ALIGN, REGREF, and REFHASH. Reject a mix
        of unregistered images or different references rather than averaging
        incompatible pixel coordinates. Single-image copies need no alignment.
        """
        reference_ids: set[tuple[str, str]] = set()
        for path in image_paths:
            header = fits.getheader(path)
            with fits.open(path, memmap=True) as hdul:
                if header.get('CRCORR') == 'Yes' and 'MASK' in hdul and 'CMASK' not in hdul:
                    raise ValueError(
                        f'{path.name}: legacy cosmic-ray MASK; restore pre-registration '
                        'copies and use separate CMASK detections before combining.'
                    )
            if (
                header.get('ALIGN') not in ('Reference', 'Registered')
                or not header.get('REGREF')
                or not header.get('REFHASH')
            ):
                raise ValueError(f'{path.name}: register Light images before combining.')
            reference_ids.add((header['REGREF'], header['REFHASH']))

        if len(reference_ids) != 1:
            raise ValueError('Light images use different registration references.')

    @staticmethod
    def _safe_name(value: str) -> str:
        """Keep metadata labels from becoming filesystem paths."""
        return re.sub(r'[^A-Za-z0-9_.\-\x27]', '_', value)

    @staticmethod
    def _output_dir(
        child: ImageFileCollection, save_location: str | Path | None,
    ) -> Path:
        """Use a requested output directory or the child's existing location."""
        if save_location in (None, ''):
            return Path(child.location or '.').resolve()
        return Path(save_location).expanduser().resolve()

    @staticmethod
    def _normalize_flats(
        image_paths: list[Path],
        scale_function: ScaleFunction,
    ) -> list[CCDData]:
        """Scale one filter's calibrated flats before outlier rejection.

        Read each CCDData and evaluate the same caller-selected scaling
        function on its complete image array. Apply that dimensionless scalar
        through CCDData arithmetic, preserving metadata and propagating any
        existing mask/uncertainty into the temporary normalized inputs.

        Only in-memory copies are changed. The final multi-flat output still
        follows the existing float32/no-mask/no-uncertainty policy. Memory must
        accommodate the normalized inputs for one filter group.
        """
        normalized_flats: list[CCDData] = []

        for path in image_paths:
            ccd: CCDData = CCDData.read(path)
            scale_factor = scale_function(ccd.data)

            if (
                not isinstance(scale_factor, Real)
                or not np.isfinite(scale_factor)
                or scale_factor <= 0
            ):
                raise ValueError(
                    f'{path.name}: flat scale must be a positive finite scalar.'
                )

            normalized: CCDData = ccd.multiply(scale_factor)
            normalized.meta = ccd.meta.copy()
            normalized_flats.append(normalized)

        return normalized_flats

    def _execute(
        self,
        plans: list[CombinePlan],
        image_type: str,
        sigma_clip_low: float,
        sigma_clip_high: float,
        *,
        scale_function: ScaleFunction | None,
        quite: bool,
        overwrite: bool | None,
    ) -> ImageInput:
        """Report exact inputs, run the original recipe, and publish masters.

        Unless quite=True, log the count, wrapped input filenames, output
        filename, and completion at INFO level. Group titles identify FILTER
        and FIELD when present; separators and blank lines distinguish groups.
        Logging handlers and levels are configured by the caller.
        """
        if not plans:
            raise ValueError(f'No {image_type} frames selected.')

        replace = self._overwrite if overwrite is None else overwrite
        check_destinations([plan[2] for plan in plans], replace)
        additions: dict[int, list[Path]] = {}

        for child, image_list_to_combine, destination, field_name, filter_name in plans:
            image_number: int = len(image_list_to_combine)
            combined: CCDData

            if not quite:
                action = 'Copying' if image_number == 1 else 'Combining'
                frame_label = 'frame' if image_number == 1 else 'frames'
                group_label: str = image_type
                if filter_name is not None:
                    group_label += f' | FILTER={filter_name}'
                if field_name is not None:
                    group_label += f' | FIELD={field_name}'

                input_names: str = fill(
                    ', '.join(path.name for path in image_list_to_combine),
                    width=88,
                    initial_indent='  ',
                    subsequent_indent='  ',
                    break_long_words=False,
                    break_on_hyphens=False,
                )
                logger.info(
                    '%s\n%s\n%s %d %s %s → %s\nInput files:\n%s',
                    group_label, '-' * 64,
                    action, image_number, image_type, frame_label,
                    destination.name, input_names,
                )

            # Keep the single-image copy path. For multiple flats, apply the
            # chosen scale before clipping; other recipes retain their order.
            if image_number == 1:
                combined = CCDData.read(image_list_to_combine[0])

            elif image_type == 'Flat':
                flat_inputs: list[Path] | list[CCDData] = image_list_to_combine
                if scale_function is not None:
                    flat_inputs = self._normalize_flats(
                        image_list_to_combine, scale_function,
                    )

                combined = ccdpro_combine(
                    flat_inputs,
                    method='average',
                    sigma_clip=True,
                    sigma_clip_low_thresh=sigma_clip_low,
                    sigma_clip_high_thresh=sigma_clip_high,
                    sigma_clip_func=np.ma.median,
                    sigma_clip_dev_func=mad_std,
                    dtype=np.float32,
                )
                # Release this filter's in-memory inputs before the next group.
                del flat_inputs

            else:
                combined = ccdpro_combine(
                    image_list_to_combine,
                    method='average',
                    sigma_clip=True,
                    sigma_clip_low_thresh=sigma_clip_low,
                    sigma_clip_high_thresh=sigma_clip_high,
                    sigma_clip_func=np.ma.median,
                    sigma_clip_dev_func=mad_std,
                    dtype=np.float32,
                )

            combined.meta['IMTYPE'] = f'Master {image_type}'
            combined.meta['NCOMBINE'] = image_number
            if filter_name is not None:
                combined.meta['FILTER'] = filter_name
            if field_name is not None:
                combined.meta['FIELD'] = field_name

            if image_number > 1 and image_type in ('Flat', 'Light'):
                combined.data = combined.data.astype(np.float32)
                combined.uncertainty = None
                if image_type == 'Flat':
                    combined.mask = None

            write_ccd(combined, destination, overwrite=replace)
            additions.setdefault(id(child), []).append(destination)

            if not quite:
                logger.info('Master written successfully.\n')

        self._image_collection = rebuild(self._image_collection, additions=additions)
        self._master_paths = tuple(plan[2] for plan in plans)
        return self._image_collection
