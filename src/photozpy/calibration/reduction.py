"""Apply master bias, dark, and flat corrections using the existing recipes."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, TypeAlias

import astropy.units as u
from astropy.io import fits
from astropy.nddata import CCDData
from ccdproc import subtract_bias, subtract_dark, flat_correct
import numpy as np
from tqdm.auto import tqdm

from ..telescope import Telescope
from ._calibration_io import (
    ImageInput, check_destinations, children, exposure, is_corrected, rebuild,
    refresh, select, validate_compatible, write_ccd,
)


CorrectionKind: TypeAlias = Literal['bias', 'dark', 'flat']
ReductionPlan: TypeAlias = tuple[Path, Path, Path]


class Reduction:
    """Apply existing master frames to working images without legacy APIs.

    Parameters
    ----------
    image_collection : ImageFileCollection or ImageCollectionGroup
        Working images and masters. Each child's file selection must explicitly
        include its masters. Pass the collection returned by Combine, or use a
        directory-backed collection to discover new files during refresh.
    telescope : Telescope or None, optional
        Refactored telescope for FILTER alias resolution. When omitted, master
        and light FILTER values must already match exactly. Gain/read noise are
        not recalculated here; these corrections use image and master arrays.
    overwrite : bool, optional
        Allow replacement of existing output files. Default is False. Set True
        explicitly when updating the reviewed working copies in analysis/.

    Notes
    -----
    The original scientific recipes are retained: subtract_bias, subtract_dark
    with EXPTIME in seconds and scale=True, and flat_correct with its default
    mean normalization. Results are float32 with mask and uncertainty removed,
    as in the old Reduction. This does not introduce uncertainty propagation.

    Masters are resolved within each child, never across nights automatically.
    Scaled dark subtraction requires bias-corrected inputs and Master Dark.
    Flat and light must have matching prior bias/dark correction states.
    A correction already marked Yes cannot be applied again, even if overwrite
    is True: overwrite controls files, not the scientific processing state.

    Methods return the entire working selection, replacing corrected file paths
    when saving elsewhere and retaining masters/uncorrected members. CCDData
    controls the result HDUs; unrelated input extensions are not copied.

    Each correction displays one progress bar with a fixed stage description.
    Its postfix shows the current master and target filenames after the stats.
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

    @property
    def image_collection(self) -> ImageInput:
        """Return the latest successful working selection, including masters."""
        return self._image_collection

    def apply_bias_correction(
        self,
        save_location: str | Path | None = None,
        *,
        overwrite: bool | None = None,
    ) -> ImageInput:
        """Subtract one Master Bias from each Dark, Flat, and Light image.

        Parameters
        ----------
        save_location : str, Path, or None, optional
            Output directory. None or '' updates files beside their inputs.
        overwrite : bool or None, optional
            Override the constructor policy; in-place correction requires True.

        Returns
        -------
        ImageFileCollection or ImageCollectionGroup
            Updated working selection, also stored in image_collection.

        Raises
        ------
        ValueError
            No targets exist, a child lacks exactly one Master Bias, geometry
            or units differ, or a target was already calibrated out of order.
        FileExistsError
            An output exists and overwriting has not been enabled.

        Notes
        -----
        Retains subtract_bias and the original float32/no-mask/no-uncertainty
        output. BIASCORR='Yes' is written into the result before serialization,
        never into an untouched source when saving to a different directory.
        """
        plans = self._plan('bias', save_location, overwrite)
        masters: dict[Path, CCDData] = {}

        with tqdm(plans, desc='Bias correction', unit='image') as progress:
            for source, master_bias_path, destination in progress:
                progress.set_postfix_str(
                    f'Applying {master_bias_path.name} to {source.name}'
                )

                if master_bias_path not in masters:
                    masters[master_bias_path] = CCDData.read(master_bias_path)

                master_bias_ccd: CCDData = masters[master_bias_path]
                ccd: CCDData = CCDData.read(source)

                # Original bias correction and output policy.
                ccd = subtract_bias(ccd, master_bias_ccd)
                ccd.data = ccd.data.astype(np.float32)
                ccd.uncertainty = None
                ccd.mask = None

                self._write_result(ccd, destination, 'BIASCORR', overwrite)

        return self._finish(plans)

    def apply_dark_correction(
        self,
        save_location: str | Path | None = None,
        *,
        overwrite: bool | None = None,
    ) -> ImageInput:
        """Subtract exposure-scaled Master Dark from Flat and Light images.

        Parameters
        ----------
        save_location : str, Path, or None, optional
            Output directory. None or '' updates the working images in place.
        overwrite : bool or None, optional
            Override the constructor's file-overwrite policy.

        Returns
        -------
        ImageFileCollection or ImageCollectionGroup
            Working images with corrected paths and retained master members.

        Raises
        ------
        ValueError
            Master selection fails, EXPTIME is missing/nonpositive, geometry or
            units differ, bias correction is absent, or DARKCORR is already Yes.
        FileExistsError
            An output exists and overwriting is disabled.

        Notes
        -----
        Preserves subtract_dark(exposure_time='EXPTIME', exposure_unit=u.second,
        scale=True). Both Master Dark and targets must first be bias-corrected
        so scaling does not also scale a bias offset. No new exposure scaling
        or read-noise model is introduced. Results receive DARKCORR='Yes'.
        """
        plans = self._plan('dark', save_location, overwrite)
        masters: dict[Path, CCDData] = {}

        with tqdm(plans, desc='Dark correction', unit='image') as progress:
            for source, master_dark_path, destination in progress:
                progress.set_postfix_str(
                    f'Applying {master_dark_path.name} to {source.name}'
                )

                if master_dark_path not in masters:
                    masters[master_dark_path] = CCDData.read(master_dark_path)

                master_dark_ccd: CCDData = masters[master_dark_path]
                ccd: CCDData = CCDData.read(source)

                # Original dark scaling, including the exposure unit and flag.
                ccd = subtract_dark(
                    ccd, master_dark_ccd,
                    exposure_time='EXPTIME', exposure_unit=u.second,
                    scale=True,
                )
                ccd.data = ccd.data.astype(np.float32)
                ccd.uncertainty = None
                ccd.mask = None

                self._write_result(ccd, destination, 'DARKCORR', overwrite)

        return self._finish(plans)

    def apply_flat_correction(
        self,
        save_location: str | Path | None = None,
        *,
        overwrite: bool | None = None,
    ) -> ImageInput:
        """Correct each Light image using its filter's unique Master Flat.

        Parameters
        ----------
        save_location : str, Path, or None, optional
            Output directory. None or '' updates working images in place.
        overwrite : bool or None, optional
            Override the constructor's file-overwrite policy.

        Returns
        -------
        ImageFileCollection or ImageCollectionGroup
            Updated working images, including all original master members.

        Raises
        ------
        ValueError
            No Light frames exist, FILTER/master selection fails, image shapes
            or units differ, prior correction states disagree, a Light was
            already flat-corrected, or the numerical output is non-finite.
        FileExistsError
            A destination exists and overwriting is disabled.

        Notes
        -----
        Uses flat_correct(object_ccd, master_flat_ccd), retaining ccdproc's
        default mean normalization. No min_value floor, extra clipping, or
        alternate normalization is silently introduced. Output is float32,
        without mask/uncertainty, and carries FLATCORR='Yes'.
        """
        plans = self._plan('flat', save_location, overwrite)
        masters: dict[Path, CCDData] = {}

        with tqdm(plans, desc='Flat correction', unit='image') as progress:
            for source, master_flat_path, destination in progress:
                progress.set_postfix_str(
                    f'Applying {master_flat_path.name} to {source.name}'
                )

                if master_flat_path not in masters:
                    masters[master_flat_path] = CCDData.read(master_flat_path)

                master_flat_ccd: CCDData = masters[master_flat_path]
                object_ccd: CCDData = CCDData.read(source)

                # Original flat correction and output policy.
                object_ccd_corrected = flat_correct(object_ccd, master_flat_ccd)
                object_ccd_corrected.data = object_ccd_corrected.data.astype(np.float32)
                object_ccd_corrected.uncertainty = None
                object_ccd_corrected.mask = None

                self._write_result(
                    object_ccd_corrected, destination, 'FLATCORR', overwrite,
                )

        return self._finish(plans)

    def _filter(self, header: fits.Header, path: Path) -> str:
        """Resolve a required FILTER using the optional refactored telescope."""
        if not header.get('FILTER'):
            raise ValueError(f'{path.name}: missing FILTER.')
        name: str = str(header['FILTER']).strip()
        return self._telescope.resolve_filter_alias(name) if self._telescope else name

    def _plan(
        self,
        kind: CorrectionKind,
        save_location: str | Path | None,
        overwrite: bool | None,
    ) -> list[ReductionPlan]:
        """Validate all target/master pairs and output paths before writing."""
        refresh(self._image_collection)
        target_types: dict[CorrectionKind, list[str]] = {
            'bias': ['Dark', 'Flat', 'Light'],
            'dark': ['Flat', 'Light'],
            'flat': ['Light'],
        }
        master_type = {'bias': 'Master Bias', 'dark': 'Master Dark', 'flat': 'Master Flat'}[kind]
        flag = {'bias': 'BIASCORR', 'dark': 'DARKCORR', 'flat': 'FLATCORR'}[kind]
        plans: list[ReductionPlan] = []
        master_paths: set[Path] = set()

        for child in children(self._image_collection):
            correct_paths: list[Path] = select(child, IMTYPE=target_types[kind])
            if not correct_paths:
                continue

            candidates = select(child, IMTYPE=master_type)
            master_paths.update(candidates)

            for source in correct_paths:
                header: fits.Header = fits.getheader(source)
                if is_corrected(header, flag):
                    raise ValueError(f'{source.name}: {flag} is already Yes; refusing double correction.')

                matching = candidates
                if kind == 'flat':
                    filter_name = self._filter(header, source)
                    matching = [
                        path for path in candidates
                        if self._filter(fits.getheader(path), path) == filter_name
                    ]

                if len(matching) != 1:
                    raise ValueError(
                        f'{source.name}: expected exactly one applicable {master_type}, '
                        f'found {len(matching)}.'
                    )

                master_path: Path = matching[0]
                master_header = fits.getheader(master_path)
                validate_compatible([source, master_path])
                self._validate_state(kind, source, header, master_header)

                if kind == 'dark':
                    exposure(header, source)
                    exposure(master_header, master_path)

                destination = (
                    source if save_location in (None, '')
                    else Path(save_location).expanduser().resolve() / source.name
                )
                plans.append((source, master_path, destination))

        if not plans:
            raise ValueError(f'No images selected for {kind} correction.')
        destinations = [plan[2] for plan in plans]
        if master_paths.intersection(destinations):
            raise ValueError('A correction destination would overwrite a master frame.')

        replace = self._overwrite if overwrite is None else overwrite
        check_destinations(destinations, replace)
        return plans

    @staticmethod
    def _validate_state(
        kind: CorrectionKind,
        source: Path,
        header: fits.Header,
        master_header: fits.Header,
    ) -> None:
        """Reject double subtraction and mismatched calibration histories."""
        if kind == 'bias':
            if is_corrected(header, 'DARKCORR') or is_corrected(header, 'FLATCORR'):
                raise ValueError(f'{source.name}: bias correction must precede dark/flat correction.')

        elif kind == 'dark':
            if not is_corrected(header, 'BIASCORR') or not is_corrected(master_header, 'BIASCORR'):
                raise ValueError(
                    f'{source.name}: bias-correct the target and Master Dark '
                    'before scaling dark.'
                )
            if is_corrected(header, 'FLATCORR'):
                raise ValueError(f'{source.name}: dark correction must precede flat correction.')

        else:
            for keyword in ('BIASCORR', 'DARKCORR'):
                if is_corrected(header, keyword) != is_corrected(master_header, keyword):
                    raise ValueError(f'{source.name}: Master Flat has a different {keyword} state.')

    def _write_result(
        self, ccd: CCDData, destination: Path, flag: str, overwrite: bool | None,
    ) -> None:
        """Put the correction marker in the output, then write it atomically."""
        if not np.all(np.isfinite(ccd.data)):
            raise ValueError(f'{destination.name}: correction produced non-finite data.')

        ccd.meta[flag] = 'Yes'
        replace = self._overwrite if overwrite is None else overwrite
        write_ccd(ccd, destination, overwrite=replace)

    def _finish(self, plans: list[ReductionPlan]) -> ImageInput:
        """Keep masters/other images and follow corrected copies when relocated."""
        replacements = {source: destination for source, _, destination in plans}
        self._image_collection = rebuild(self._image_collection, replacements=replacements)
        return self._image_collection
