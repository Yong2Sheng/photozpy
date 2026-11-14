"""Add celestial WCS using the locally installed Astrometry.net solver.

Originally written by Yong Sheng at Clemson University, 2023, for Photozpy.
Advisor: Dr. Marco Ajello.

Retains the original per-file solve-field call and batch failure reporting.
Collection selection, output validation, and FITS publication use the refactored
workflow. Online solving and the legacy Sources dependency have been removed.
"""

from __future__ import annotations

from collections.abc import Sequence
import logging
from numbers import Real
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile

import astropy.units as u
from astropy.coordinates import Angle
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from tqdm.auto import tqdm

from ..fits_utilities import update_fits_header
from ._calibration_io import (
    ImageInput, check_destinations, children, rebuild, refresh, select,
)


logger = logging.getLogger(__name__)

# Replace old coordinate transforms, including alternate WCS and SIP, rather
# than mixing obsolete CD/PC or distortion terms with the new solution.
_WCS_KEY = re.compile(
    r'^(?:WCSAXES|WCSNAME|CTYPE\d+|CUNIT\d+|CRPIX\d+|CRVAL\d+|'
    r'CDELT\d+|CROTA\d+|(?:CD|PC)\d+_\d+|(?:PV|PS)\d+_\d+|'
    r'LONPOLE|LATPOLE|RADESYS|RADECSYS|EQUINOX)[A-Z]?$'
    r'|^(?:A|B|AP|BP)_(?:ORDER|DMAX|\d+_\d+)$'
    r'|^(?:CPDIS|CQDIS|D2IMDIS|DET2IM)\d+$|^(?:DP|DQ|D2IM)\d+(?:\..*)?$'
)


class PlateSolving:
    """Solve selected images locally, without Sources, FWHM, or network uploads.

    Parameters
    ----------
    image_collection : ImageFileCollection or ImageCollectionGroup
        Working collection containing Master Light products. Directory
        boundaries and explicit selections are retained.
    overwrite : bool, optional
        Permit replacing destination FITS files, including in-place WCS updates.
        Default False. No process is launched by the constructor.

    Notes
    -----
    Astrometry.net and suitable local index files must already be installed.
    The Python package remains importable when solve-field is unavailable.
    Images are not resampled; only primary-header WCS cards are replaced.
    All science data, masks, units, and unrelated HDUs remain intact. CMASK and
    MASK are not applied to source extraction by this image-based solver.
    """

    def __init__(self, image_collection: ImageInput, *, overwrite: bool = False) -> None:
        children(image_collection)
        self._image_collection: ImageInput = image_collection
        self._overwrite: bool = overwrite
        self._solved_paths: tuple[Path, ...] = ()
        self._failed_paths: dict[Path, str] = {}

    @property
    def image_collection(self) -> ImageInput:
        """Return the latest working selection, including successfully solved paths."""
        return self._image_collection

    @property
    def solved_paths(self) -> tuple[Path, ...]:
        """Return outputs successfully written during the latest batch."""
        return self._solved_paths

    @property
    def failed_paths(self) -> dict[Path, str]:
        """Return input paths and diagnostics for failures in the latest batch."""
        return self._failed_paths.copy()

    @staticmethod
    def _positive(value: float, name: str) -> float:
        """Reject invalid numerical controls before launching an external process."""
        if isinstance(value, bool) or not isinstance(value, Real) or not np.isfinite(value) or value <= 0:
            raise ValueError(f'{name} must be positive and finite.')
        return float(value)

    @staticmethod
    def _executable(executable: str | Path) -> str:
        """Search PATH, then common Homebrew locations for the default solver.

        An explicit executable remains authoritative: a miss is an error rather
        than silently switching to a different installation.
        """
        found = shutil.which(str(Path(executable).expanduser()))
        if found is None and str(executable) == 'solve-field':
            for candidate in ('/opt/homebrew/bin/solve-field', '/usr/local/bin/solve-field'):
                found = shutil.which(candidate)
                if found is not None:
                    break
        if found is None:
            raise FileNotFoundError(
                f'Local solver not found: {executable}. Searched PATH and, for '
                'the default solve-field, common Homebrew locations. Install '
                'Astrometry.net or specify executable for a custom installation.'
            )
        return found

    @staticmethod
    def _extra_options(extra_args: str | Sequence[str] | None) -> list[str]:
        """Allow solver tuning but keep input/output ownership inside this call."""
        args = shlex.split(extra_args) if isinstance(extra_args, str) else list(extra_args or ())
        # An explicit tuning list prevents solver output redirection or extra
        # input files, including abbreviated aliases accepted by getopt_long.
        valued = {
            '--config', '--index-dir', '--index-file', '--scale-units',
            '--scale-low', '--scale-high', '--downsample', '--cpulimit',
            '--depth', '--objs', '--parity', '--tweak-order', '--pixel-error',
            '--odds-to-solve', '--odds-to-tune-up', '--sigma', '--uniformize',
        }
        flags = {'--no-verify', '--no-tweak', '--guess-scale', '--resort',
                 '--no-remove-lines', '--no-background-subtraction'}
        index = 0
        while index < len(args):
            token = args[index]
            if not isinstance(token, str):
                raise TypeError('extra_args must contain strings.')
            option, separator, value = token.partition('=')
            if option in flags and not separator:
                index += 1
                continue
            if option not in valued:
                raise ValueError(f'Unsupported extra_args option: {token}. Use named parameters or supported long tuning options.')
            if not separator:
                index += 1
                if index >= len(args) or not isinstance(args[index], str) or args[index].startswith('--'):
                    raise ValueError(f'Missing value for extra_args option: {option}.')
            elif not value:
                raise ValueError(f'Missing value for extra_args option: {option}.')
            if option in ('--config', '--index-dir', '--index-file'):
                # The solver runs inside its private work directory; resolve
                # user paths while still in the notebook's original directory.
                path = str(Path(value if separator else args[index]).expanduser().resolve())
                args[index] = f'{option}={path}' if separator else path
            index += 1
        return args

    @staticmethod
    def _pointing(
        header: fits.Header, ra: str | float | None, dec: str | float | None,
        ra_unit: str,
    ) -> tuple[float, float] | None:
        """Convert pointing hints to degrees without silently guessing numeric RA.

        Sexagesimal RA (colon/space separated or with h units) is hours in
        auto mode; a plain numeric RA is degrees. ra_unit='hourangle' explicitly
        supports decimal hours. DEC is always interpreted in degrees.
        """
        if ra_unit not in ('auto', 'deg', 'hourangle'):
            raise ValueError("ra_unit must be 'auto', 'deg', or 'hourangle'.")
        if (ra is None) != (dec is None):
            raise ValueError('Supply both ra and dec, or neither to use FITS pointing.')
        if ra is None:
            ra, dec = header.get('RA'), header.get('DEC')
            if ra is None or dec is None:
                return None
        unit = ra_unit
        if unit == 'auto':
            text = str(ra).strip()
            unit = 'hourangle' if ':' in text or len(text.split()) > 1 or 'h' in text.lower() else 'deg'
        ra_deg = float(Angle(ra, unit=u.hourangle if unit == 'hourangle' else u.deg).degree)
        dec_deg = float(Angle(dec, unit=u.deg).degree)
        if not np.all(np.isfinite([ra_deg, dec_deg])) or not -90 <= dec_deg <= 90:
            raise ValueError('Pointing must be finite, with DEC between -90 and 90 degrees.')
        return ra_deg % 360, dec_deg

    @staticmethod
    def _validated_wcs(path: Path, shape: tuple[int, int]) -> fits.Header:
        """Require a finite, nonsingular 2D celestial solution, retaining SIP."""
        header = fits.getheader(path)
        # Solver .wcs files may be header-only (NAXIS=0). Supply the known
        # science geometry for validation without modifying that file.
        header['NAXIS'] = 2
        header['NAXIS2'], header['NAXIS1'] = shape
        try:
            wcs = WCS(header, relax=True, fix=False)
        except ValueError as exc:
            raise RuntimeError(f'Solver produced an invalid celestial WCS: {exc}') from exc
        if not wcs.has_celestial or wcs.pixel_n_dim != 2 or wcs.world_n_dim != 2:
            raise RuntimeError('Solver did not produce a two-dimensional celestial WCS.')
        matrix = wcs.pixel_scale_matrix
        if not np.all(np.isfinite(matrix)) or np.linalg.det(matrix) == 0:
            raise RuntimeError('Solver produced a non-finite or singular WCS scale.')
        ny, nx = shape
        positions = np.array([[0., 0.], [nx - 1., ny - 1.], [(nx - 1.) / 2, (ny - 1.) / 2]])
        if not np.all(np.isfinite(wcs.all_pix2world(positions, 0))):
            raise RuntimeError('Solver WCS gives non-finite image coordinates.')
        return wcs.to_header(relax=True)

    @staticmethod
    def add_wcs_locally_for_file(
        file_path: str | Path, ra: str | float | None = None,
        dec: str | float | None = None, radius: float = 10,
        detect_threshold: float = 8, extra_args: str | Sequence[str] | None = None,
        clean: bool = True, *, output_path: str | Path | None = None,
        overwrite: bool = False, executable: str | Path = 'solve-field',
        timeout: float | None = None, ra_unit: str = 'auto',
    ) -> subprocess.CompletedProcess[str]:
        """Solve one FITS locally and atomically publish only its new WCS.

        Parameters
        ----------
        file_path : str or Path
            FITS with a two-dimensional primary image.
        ra, dec : str, float, or None
            Optional search-center pair. If omitted, read RA/DEC headers; if
            those are incomplete, solve without a pointing hint. Numeric RA is
            degrees in auto mode, sexagesimal RA hours; DEC is degrees.
        radius : float, optional
            Search radius in degrees when a pointing hint is available. Default 10.
        detect_threshold : float, optional
            Source-detection significance passed to --nsigma. Default 8 retains
            the old local recipe. This is not a FWHM or coordinate tolerance.
        extra_args : str, sequence of str, or None
            Long solver tuning options, e.g. ['--downsample', '2'] or
            ['--scale-units', 'arcsecperpix', '--scale-low', '0.3',
            '--scale-high', '0.4']. Strings are parsed with shlex, never a shell.
            Input/output redirection and named-parameter duplicates are rejected.
        clean : bool, optional
            Retained legacy entry: True defaults to updating the input, False
            to <stem>solved.fits. Explicit output_path overrides that destination.
            Temporary solver files are always confined to a disposable directory.
        output_path : str, Path, or None
            Explicit output FITS path. All original HDUs and pixel values remain.
        overwrite : bool, optional
            Required to replace any existing output, including the input itself.
        executable : str or Path, optional
            Automatically search PATH then common Homebrew locations by default.
            An explicit name/path overrides discovery and must be executable.
        timeout : float or None, optional
            Wall-clock limit in seconds for the solver process. None adds no limit.
        ra_unit : {'auto', 'deg', 'hourangle'}
            RA hint interpretation; use hourangle explicitly for decimal hours.

        Returns
        -------
        subprocess.CompletedProcess[str]
            Captured solver command, status, stdout and stderr, as before.
            Success requires exit status zero, a successful .solved marker, and
            a usable .wcs file. FWHM is neither required nor calculated.

        Raises
        ------
        FileNotFoundError, ValueError, FileExistsError
            Missing input/executable, invalid controls, or forbidden replacement.
        RuntimeError
            Unsolved field, process timeout/failure, or invalid WCS. No FITS
            destination is published until a valid solution has been obtained.
        """
        source = Path(file_path).expanduser().resolve()
        if not source.is_file():
            raise FileNotFoundError(source)
        header = fits.getheader(source)
        if header.get('NAXIS') != 2:
            raise ValueError(f'{source.name}: expected a 2D primary image.')
        threshold = PlateSolving._positive(detect_threshold, 'detect_threshold')
        radius = PlateSolving._positive(radius, 'radius')
        if timeout is not None:
            timeout = PlateSolving._positive(timeout, 'timeout')
        options = PlateSolving._extra_options(extra_args)
        pointing = PlateSolving._pointing(header, ra, dec, ra_unit)
        solver = PlateSolving._executable(executable)
        destination = (
            Path(output_path).expanduser().resolve() if output_path is not None
            else source if clean else source.with_name(f'{source.stem}solved.fits')
        )
        check_destinations([destination], overwrite)

        with tempfile.TemporaryDirectory(prefix='photozpy-solve-') as directory:
            work = Path(directory)
            # A private byte-for-byte copy also protects the input from external
            # tools and keeps apostrophes/spaces in product filenames out of their
            # intermediate command construction.
            staged = work / 'input.fits'
            shutil.copyfile(source, staged)
            wcs_path, solved_path = work / 'solution.wcs', work / 'solution.solved'
            cmd: list[str] = [
                solver, str(staged), '--overwrite', '--no-plots',
                '--dir', str(work), '--out', 'solution', '--new-fits', 'none',
                '--wcs', str(wcs_path), '--solved', str(solved_path),
                '--nsigma', str(threshold),
            ]
            if pointing is not None:
                cmd += ['--ra', str(pointing[0]), '--dec', str(pointing[1]), '--radius', str(radius)]
            cmd.extend(options)
            try:
                result: subprocess.CompletedProcess[str] = subprocess.run(
                    cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    text=True, cwd=work, timeout=timeout,
                )
            except subprocess.TimeoutExpired as exc:
                raise RuntimeError(f'{source.name}: solve-field timed out after {timeout} seconds.') from exc
            if result.returncode != 0:
                raise RuntimeError(
                    f'{source.name}: solve-field exited with {result.returncode}. '
                    f'stdout={result.stdout[-2000:]!r}; stderr={result.stderr[-2000:]!r}'
                )
            if not solved_path.is_file() or solved_path.read_bytes()[:1] != b'\x01' or not wcs_path.is_file():
                raise RuntimeError(
                    f'{source.name}: solve-field did not solve the image. '
                    f'stdout={result.stdout[-2000:]!r}; stderr={result.stderr[-2000:]!r}'
                )
            solved_header = PlateSolving._validated_wcs(
                wcs_path, (header['NAXIS2'], header['NAXIS1']),
            )
            # WCS.to_header can also include observation-time metadata. Keep
            # those original scientific headers; replace coordinate cards only.
            values = {
                card.keyword: (card.value, card.comment)
                for card in solved_header.cards if _WCS_KEY.fullmatch(card.keyword)
            }
            values['HISTORY'] = 'Celestial WCS solved locally with Astrometry.net.'
            update_fits_header(
                source, values, output_path=destination, overwrite=overwrite,
                remove_keys=[key for key in header if _WCS_KEY.fullmatch(key)],
            )
        return result

    def add_wcs_locally(
        self, image_type: str = 'Master Light', detect_threshold: float = 8, *,
        radius: float = 10, extra_args: str | Sequence[str] | None = None,
        save_location: str | Path | None = None, overwrite: bool | None = None,
        executable: str | Path = 'solve-field', timeout: float | None = None,
        ra_unit: str = 'auto',
    ) -> ImageInput:
        """Solve a batch locally, retaining the original continue-on-failure loop.

        image_type selects IMTYPE (default Master Light). detect_threshold is
        forwarded unchanged to every per-file --nsigma call; radius is degrees.
        extra_args, executable, timeout, and ra_unit follow the per-file method.
        save_location=None updates selected working images in place; otherwise
        successful copies go to that directory. overwrite=None uses the
        constructor policy. Replacing files always requires explicit permission
        through overwrite=True, including when replacing an existing WCS.

        Return the updated collection, retaining non-selected images. Failed
        inputs remain in the selection unchanged; inspect failed_paths and logs
        before proceeding to photometry. solved_paths lists only successful
        outputs. The progress bar keeps its title fixed and current filename at
        the end. Results are committed one file at a time, not as a transaction
        for the whole batch. No online fallback or automatic threshold retries.
        """
        refresh(self._image_collection)
        selected = [
            path for child in children(self._image_collection)
            for path in select(child, IMTYPE=image_type)
        ]
        if not selected:
            raise ValueError(f'No {image_type} images selected for local plate solving.')
        replace = self._overwrite if overwrite is None else overwrite
        output = None if save_location in (None, '') else Path(save_location).expanduser().resolve()
        destinations = [path if output is None else output / path.name for path in selected]
        check_destinations(destinations, replace)
        solver = self._executable(executable)
        self._positive(detect_threshold, 'detect_threshold')
        self._positive(radius, 'radius')
        if timeout is not None:
            self._positive(timeout, 'timeout')
        options = self._extra_options(extra_args)
        # Catch invalid metadata/settings before any image is published.
        for path in selected:
            header = fits.getheader(path)
            if header.get('NAXIS') != 2:
                raise ValueError(f'{path.name}: expected a 2D primary image.')
            self._pointing(header, None, None, ra_unit)

        self._solved_paths, self._failed_paths = (), {}
        replacements: dict[Path, Path] = {}
        with tqdm(zip(selected, destinations), total=len(selected), desc='Plate solving', unit='image') as progress:
            for source, destination in progress:
                progress.set_postfix_str(f'Solving {source.name}')
                try:
                    self.add_wcs_locally_for_file(
                        source, radius=radius, detect_threshold=detect_threshold,
                        extra_args=options, output_path=destination, overwrite=replace,
                        executable=solver, timeout=timeout, ra_unit=ra_unit,
                    )
                except (OSError, ValueError, RuntimeError) as exc:
                    self._failed_paths[source] = str(exc)
                    logger.error('%s: local WCS solving failed: %s', source.name, exc)
                else:
                    replacements[source] = destination
                    self._solved_paths += (destination,)
                    logger.info('WCS solved: %s', source.name)

        self._image_collection = rebuild(self._image_collection, replacements=replacements)
        logger.info('Local plate solving: %d solved, %d failed.', len(self._solved_paths), len(self._failed_paths))
        return self._image_collection
