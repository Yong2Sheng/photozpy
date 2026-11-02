"""Remove cosmic rays from reduced Light frames before image registration.

Originally written by Yong Sheng at Clemson University, 2023, for Photozpy.
Advisor: Dr. Marco Ajello.

The original L.A.Cosmic call, header-derived gain/read noise, and
gain_apply=False are retained. Collection and output handling use the refactored workflow.
Only binary masks are plotted; image previews belong in DS9.
"""

from __future__ import annotations

from collections.abc import Mapping
import inspect
import logging
from numbers import Real
import os
from pathlib import Path
import shutil
import tempfile
from typing import Any

import astropy.units as u
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.visualization import ImageNormalize
from ccdproc import cosmicray_lacosmic as lacosmic
import numpy as np
from tqdm.auto import tqdm

from ..fits_utilities import plot_image
from ._calibration_io import (
    ImageInput, check_destinations, children, is_corrected, rebuild,
    refresh, select, write_ccd,
)


logger = logging.getLogger(__name__)


class Artifacts:
    """Clean cosmic rays using the existing CCDProc L.A.Cosmic implementation.

    Parameters
    ----------
    image_collection : ImageFileCollection or ImageCollectionGroup
        Reviewed working images after bias, dark, and flat correction. Only
        IMTYPE='Light' is selected by default, including standard fields.
    overwrite : bool, optional
        Permit replacement of existing destination files. Default is False;
        use True explicitly to update working images in analysis/.

    Notes
    -----
    Keep this as an explicit step before Registration, not an automatic hidden
    part of alignment. GAIN is interpreted as electron/ADU and RDNOISE as the
    numerical per-pixel RMS read noise in electrons. gain_apply=False retains
    ADU values and units; no extra gain multiplication is performed.

    Original L.A.Cosmic defaults are used unless keyword arguments override
    them (currently sigclip=4.5, sigfrac=0.3, objlim=5.0, niter=4 in CCDProc 2.5.1).
    Thresholds can mistake sharp stellar cores for cosmic rays: inspect masks
    on representative frames rather than treating detections as ground truth.

    Cleaned pixel values are saved with a separate CMASK detection extension.
    CMASK supports inspection and registration source detection, not exclusion
    from Light combination. Any pre-existing validity MASK remains separate.
    As in the original code, uncertainty is omitted from the output.
    Unrelated FITS extensions are not copied; inputs must have a primary image
    and optionally MASK. CRCORR prevents accidental repeat processing. ALIGN
    markers prevent running this stage after registration has already occurred.
    """

    def __init__(self, image_collection: ImageInput, *, overwrite: bool = False) -> None:
        children(image_collection)
        self._image_collection: ImageInput = image_collection
        self._overwrite: bool = overwrite
        self._cosmic_ray_counts: dict[Path, int] = {}
        self._backup_paths: dict[Path, Path] = {}

    @property
    def image_collection(self) -> ImageInput:
        """Return working images with cleaned output paths after success."""
        return self._image_collection

    @property
    def cosmic_ray_counts(self) -> dict[Path, int]:
        """Return new masked-pixel counts per output from the latest success."""
        return self._cosmic_ray_counts.copy()

    @property
    def backup_paths(self) -> dict[Path, Path]:
        """Map cleaned outputs to retained pre-cleaning FITS for comparison."""
        return self._backup_paths.copy()

    @staticmethod
    def _copy_backup(source: Path, destination: Path) -> None:
        """Publish an exact pre-cleaning copy without replacing older backups."""
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary: Path | None = None
        try:
            with tempfile.NamedTemporaryFile(
                prefix='.', suffix='.fits', dir=destination.parent, delete=False,
            ) as output:
                temporary = Path(output.name)
                with source.open('rb') as input_file:
                    shutil.copyfileobj(input_file, output)
            shutil.copystat(source, temporary)
            os.link(temporary, destination)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)

    def preview_in_ds9(self) -> None:
        """Browse the last completed batch without cleaning or moving files.

        Start or reuse a dedicated DS9 comparison window after processing.
        Previous/Next control navigation; Quit closes the controls without
        deleting backups. Before/After/Mask occupy frames 1/2/3. Requires the
        retained pre-cleaning paths; errors do not undo completed cleaning.
        """
        from ..interactive_ds9.comparison_reviewer import ImageComparisonReviewer

        if not self._backup_paths:
            raise ValueError('No completed cleaning batch is available for comparison.')

        pairs = [(before, after) for after, before in self._backup_paths.items()]
        ImageComparisonReviewer(pairs).run()

    @staticmethod
    def _validate_header(header: fits.Header, path: Path) -> tuple[float, float]:
        """Require the original processing stage and usable detector metadata."""
        if header.get('IMTYPE') != 'Light':
            raise ValueError(f'{path.name}: cosmic-ray removal requires Light frames.')

        for keyword in ('BIASCORR', 'DARKCORR', 'FLATCORR'):
            if not is_corrected(header, keyword):
                raise ValueError(f'{path.name}: complete {keyword} before cosmic-ray removal.')

        if header.get('ALIGN') in ('Reference', 'Registered'):
            raise ValueError(f'{path.name}: remove cosmic rays before registration.')
        if is_corrected(header, 'CRCORR'):
            raise ValueError(f'{path.name}: CRCORR is already Yes; refusing repeat cleaning.')
        if str(header.get('BUNIT', '')).strip().lower() != 'adu':
            raise ValueError(f'{path.name}: cosmic-ray removal expects ADU input.')
        if header.get('NAXIS') != 2:
            raise ValueError(f'{path.name}: expected a 2D primary image.')

        try:
            gain = float(header['GAIN'])
            rdnoise = float(header['RDNOISE'])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f'{path.name}: missing or invalid GAIN/RDNOISE.') from exc
        if not np.isfinite(gain) or gain <= 0:
            raise ValueError(f'{path.name}: GAIN must be positive and finite.')
        if not np.isfinite(rdnoise) or rdnoise < 0:
            raise ValueError(f'{path.name}: RDNOISE must be nonnegative and finite.')

        return gain, rdnoise

    @staticmethod
    def _saturation_electrons(saturation_adu: float | None, gain: float) -> float | None:
        """Convert a positive finite ADU threshold using this frame's gain."""
        if saturation_adu is None:
            return None

        if (
            isinstance(saturation_adu, bool)
            or not isinstance(saturation_adu, Real)
            or not np.isfinite(saturation_adu)
            or saturation_adu <= 0
        ):
            raise ValueError('saturation_adu must be a positive finite number or None.')

        satlevel: float = float(saturation_adu) * gain
        if not np.isfinite(satlevel) or satlevel <= 0:
            raise ValueError('saturation_adu multiplied by GAIN must be positive and finite.')

        return satlevel

    @staticmethod
    def _validate_kwargs(
        kwargs: Mapping[str, Any], saturation_adu: float | None = None,
    ) -> None:
        """Keep header-based noise parameters and ADU output explicit."""
        removed_options = {'plot_before_after', 'norm_percent'}.intersection(kwargs)
        if removed_options:
            raise TypeError(
                f'{sorted(removed_options)} were removed; use DS9 for image '
                'previews and plot_mask=True for a binary mask PNG.'
            )

        Artifacts._saturation_electrons(saturation_adu, gain=1.)
        if saturation_adu is not None and 'satlevel' in kwargs:
            raise ValueError('Use saturation_adu or legacy satlevel, not both.')

        reserved = {'gain', 'readnoise', 'gain_apply'}.intersection(kwargs)
        if reserved:
            raise ValueError(
                f'{sorted(reserved)} are controlled by headers and gain_apply=False.'
            )
        # Detect unsupported options before beginning a batch, without calling
        # the numerical algorithm or importing an alternative implementation.
        inspect.signature(lacosmic).bind(
            None, gain_apply=False, gain=1., readnoise=1., **kwargs,
        )

    @staticmethod
    def _remove_cosmic_ray(
        fits_path: str | Path,
        *,
        plot_mask: bool = False,
        save_location: str | Path | None = None,
        overwrite: bool = False,
        saturation_adu: float | None = None,
        **kwargs: Any,
    ) -> Path:
        """Clean one reduced Light frame, keeping the original L.A.Cosmic call.

        GAIN and RDNOISE come from the image header. Additional keyword options
        are forwarded unchanged, except gain/readnoise/gain_apply are reserved.
        saturation_adu is converted to electrons using this frame's GAIN.
        None retains the library default or a legacy satlevel keyword; supplying
        saturation_adu and satlevel together is an error.
        Write the cleaned ADU image, separate CMASK, CRCORR='Yes', and CRNPIX
        (newly masked candidate pixels). Return the written path.

        plot_mask saves *_cr_mask.png with a fixed black/white range of 0 to 1.
        Image previews and before/after comparisons are left to DS9. No plots
        are generated by default. Only the requested working file is changed.
        """
        fits_path = Path(fits_path).expanduser().resolve()
        destination = (
            fits_path if save_location in (None, '')
            else Path(save_location).expanduser().resolve() / fits_path.name
        )
        Artifacts._validate_kwargs(kwargs, saturation_adu)
        gain, rdnoise = Artifacts._validate_header(fits.getheader(fits_path), fits_path)
        satlevel: float | None = Artifacts._saturation_electrons(saturation_adu, gain)
        lacosmic_kwargs: dict[str, Any] = dict(kwargs)
        if satlevel is not None:
            lacosmic_kwargs['satlevel'] = satlevel

        check_destinations([destination], overwrite)

        with fits.open(fits_path, memmap=True) as hdul:
            if any(hdu.name != 'MASK' for hdu in hdul[1:]):
                raise ValueError(f'{fits_path.name}: unsupported additional FITS extensions.')

        # Original CCDData input, with native byte order for the C extension.
        ccd: CCDData = CCDData.read(fits_path)
        ccd.data = np.asarray(ccd.data, dtype=ccd.data.dtype.newbyteorder('='))
        if ccd.unit != u.adu or not np.all(np.isfinite(ccd.data)):
            raise ValueError(f'{fits_path.name}: expected finite ADU image data.')
        input_mask = (
            np.zeros(ccd.shape, dtype=bool) if ccd.mask is None
            else np.asarray(ccd.mask, dtype=bool).copy()
        )

        # Retain Yong's original algorithm, detector parameters, and ADU policy.
        newccd: CCDData = lacosmic(
            ccd,
            gain_apply=False,
            gain=gain,
            readnoise=rdnoise,
            **lacosmic_kwargs,
        )
        # CCDProc 2.5.1 labels CCDData as electrons even with gain_apply=False.
        # Retain the original explicit BUNIT='adu' policy: the returned numbers
        # have not been gain multiplied, so restore the input unit, not pixels.
        newccd.unit = ccd.unit

        if newccd.shape != ccd.shape:
            raise RuntimeError('L.A.Cosmic returned incompatible image geometry.')
        if not np.all(np.isfinite(newccd.data)):
            raise RuntimeError('L.A.Cosmic returned non-finite image values.')

        output_mask = (
            np.zeros(newccd.shape, dtype=bool) if newccd.mask is None
            else np.asarray(newccd.mask, dtype=bool)
        )
        cosmic_mask: np.ndarray = output_mask & ~input_mask
        newccd.mask = input_mask if np.any(input_mask) else None
        newccd.uncertainty = None
        newccd.meta['CRCORR'] = 'Yes'
        newccd.meta['CRNPIX'] = int(np.count_nonzero(cosmic_mask))
        newccd.meta['CRMETHOD'] = 'LACosmic'
        write_ccd(newccd, destination, overwrite=overwrite, cmask=cosmic_mask)

        if plot_mask:
            plot_image(
                ccd_data=CCDData(cosmic_mask.astype(float), unit=u.dimensionless_unscaled),
                norm=ImageNormalize(vmin=0., vmax=1.),
                cmap='gray',
                save_path=destination.parent / f'{fits_path.stem}_cr_mask.png',
                close=True,
            )

        return destination

    def remove_cosmic_ray(
        self,
        *,
        plot_mask: bool = False,
        save_location: str | Path | None = None,
        headers_values: Mapping[str, object] | None = None,
        overwrite: bool | None = None,
        saturation_adu: float | None = None,
        backup_dir: str | Path | None = None,
        preview: bool = False,
        **kwargs: Any,
    ) -> ImageInput:
        """Clean selected Light frames and return the updated working collection.

        Parameters
        ----------
        plot_mask : bool, optional
            Save *_cr_mask.png beside each output, using a fixed grayscale
            range: black=0 (unmasked), white=1 (masked). Default is False.
            This PNG setting is independent of the DS9 preview option.
        save_location : str, Path, or None, optional
            Output directory; None or '' updates working images in place.
        headers_values : mapping or None, optional
            Extra header criteria, such as FIELD or FILTER. Selection always
            remains restricted to Light; calibration frames cannot be included.
        overwrite : bool or None, optional
            Override the constructor file-overwrite policy. CRCORR still blocks
            repeat cleaning; overwriting is not permission to change science twice.
        saturation_adu : float or None, optional
            Positive finite saturation threshold in the input image's ADU units,
            for example 65000. Convert separately for each frame using its GAIN;
            users do not need to look up or multiply the gain themselves.
            None preserves the previous CCDProc default. This setting neither
            clips pixel values nor changes output units, which remain ADU.
        backup_dir : str, Path, or None, optional
            Retain exact pre-cleaning copies here. For in-place processing,
            None uses before_cosmic_ray/ under each input directory. Existing
            backup files are never overwritten. With separate outputs and no
            explicit backup_dir, the unchanged input files serve as Before.
            Backups are excluded from the returned collection and never deleted
            automatically. Delete the backup directory yourself after checking.
        preview : bool, optional
            Open a user-paced Before/After/Mask comparison in DS9 only after
            the whole batch succeeds. Default False. The notebook cell waits
            until Quit; on macOS the Tk controls run in a separate process.
        **kwargs
            Original CCDProc L.A.Cosmic options, e.g. sigclip, objlim, niter,
            retaining library defaults when omitted. Prefer saturation_adu for
            saturation. Legacy satlevel (electrons) is accepted only when
            saturation_adu is None; conflicting units are never resolved silently.

        Returns
        -------
        ImageFileCollection or ImageCollectionGroup
            Working selection with cleaned paths, retaining calibration members.
            Also stored in image_collection. cosmic_ray_counts records the new
            masked-pixel counts, not a guarantee that every detection was cosmic.

        Notes
        -----
        Header prerequisites and destination collisions are checked across the
        batch before writing. A later numerical or plotting failure can leave
        earlier outputs written. Run on a representative subset and inspect the
        mask before cleaning the full night, especially for undersampled stars.
        """
        Artifacts._validate_kwargs(kwargs, saturation_adu)
        refresh(self._image_collection)
        criteria: dict[str, object] = {
            key.upper(): value for key, value in (headers_values or {}).items()
        }
        if 'IMTYPE' in criteria and criteria['IMTYPE'] != 'Light':
            raise ValueError("Cosmic-ray selection must use IMTYPE='Light'.")
        criteria['IMTYPE'] = 'Light'

        selected: list[Path] = []
        for child in children(self._image_collection):
            selected.extend(select(child, **criteria))
        if not selected:
            raise ValueError('No Light images selected for cosmic-ray removal.')

        replace: bool = self._overwrite if overwrite is None else overwrite
        destinations = [
            path if save_location in (None, '')
            else Path(save_location).expanduser().resolve() / path.name
            for path in selected
        ]
        for path in selected:
            gain, _ = Artifacts._validate_header(fits.getheader(path), path)
            Artifacts._saturation_electrons(saturation_adu, gain)
            with fits.open(path, memmap=True) as hdul:
                if any(hdu.name != 'MASK' for hdu in hdul[1:]):
                    raise ValueError(f'{path.name}: unsupported additional FITS extensions.')
        check_destinations(destinations, replace)

        before_paths: dict[Path, Path] = {}
        copies: dict[Path, Path] = {}
        for source, destination in zip(selected, destinations):
            if backup_dir is not None:
                before = Path(backup_dir).expanduser().resolve() / source.name
            elif source == destination:
                before = source.parent / 'before_cosmic_ray' / source.name
            else:
                before = source

            before_paths[destination] = before
            if before != source:
                copies[source] = before
            elif source == destination:
                raise ValueError('backup_dir must not replace the input image.')

        if set(copies.values()).intersection(set(selected) | set(destinations)):
            raise ValueError('Backup paths must not overlap input or output images.')
        check_destinations(list(copies.values()), overwrite=False)

        # Complete every backup before modifying the first working image.
        for source, backup in copies.items():
            self._copy_backup(source, backup)

        replacements: dict[Path, Path] = {}
        counts: dict[Path, int] = {}
        with tqdm(selected, desc='Cosmic-ray removal', unit='image') as progress:
            for path in progress:
                progress.set_postfix_str(path.name)
                destination = self._remove_cosmic_ray(
                    path, plot_mask=plot_mask, save_location=save_location,
                    overwrite=replace, saturation_adu=saturation_adu, **kwargs,
                )
                replacements[path] = destination
                counts[destination] = fits.getheader(destination)['CRNPIX']

        self._image_collection = rebuild(self._image_collection, replacements=replacements)
        self._cosmic_ray_counts = counts
        self._backup_paths = before_paths
        logger.info(
            'Cosmic-ray removal: %d Light images; %d newly masked pixels.\n',
            len(counts), sum(counts.values()),
        )
        if copies:
            logger.info('Pre-cleaning backups retained in: %s', ', '.join(
                str(path) for path in sorted({path.parent for path in copies.values()})
            ))
        if preview:
            self.preview_in_ds9()

        return self._image_collection
