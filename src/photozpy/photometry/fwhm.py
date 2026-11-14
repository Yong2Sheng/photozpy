"""Measure a representative stellar FWHM independently for each Master Light.

Originally written by Yong Sheng at Clemson University, 2023, for Photozpy.
Advisor: Dr. Marco Ajello.

Retains the original detect-peaks, centroid, measure-stars, and sigma-clip
workflow. The old radial-profile half-height estimate seeds a bounded circular
Moffat fit; it is no longer the final width measurement. No image is resampled.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from pathlib import Path
import logging
from numbers import Integral, Real

from astropy.io import fits
from astropy.stats import mad_std, sigma_clip, sigma_clipped_stats
from astropy.table import Table
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
from photutils.centroids import centroid_com, centroid_sources
from photutils.detection import find_peaks
from photutils.profiles import RadialProfile
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
import sep
from tqdm.auto import tqdm

from ..calibration._calibration_io import ImageInput, children, rebuild, refresh, select
from ..fits_utilities import update_fits_header


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class _Star:
    """One candidate's fit diagnostics, including explicit rejection reasons."""

    image: str
    x: float
    y: float
    fwhm: float = np.nan
    gamma: float = np.nan
    beta: float = np.nan
    peak: float = np.nan
    background: float = np.nan
    residual_rms: float = np.nan
    axis_ratio: float = np.nan
    accepted: bool = False
    reason: str = ''


@dataclass(frozen=True)
class _Frame:
    """An image's clipped mean and between-star scatter, both in pixels."""

    image: str
    field: str
    filter: str
    fwhm: float
    scatter: float
    nstars: int
    ncandidates: int


class FWHM:
    """Measure bright, isolated star widths with circular 2D Moffat profiles.

    Parameters
    ----------
    image_collection : ImageFileCollection or ImageCollectionGroup
        Combined working images. Each selected Master Light is measured
        independently; different fields, filters, and directories are not pooled.
    overwrite : bool, optional
        Permit updating FITS headers. Default False; required when write_fwhm=True.

    Notes
    -----
    Moffat uses I(r) = background + amplitude * (1 + r**2/gamma**2)**(-beta).
    FWHM = 2 * gamma * sqrt(2**(1/beta) - 1) is a full width in pixels,
    not a radius or an angular seeing estimate. A circular model supplies one
    unambiguous width; markedly elongated stars are rejected, not circularized.
    No Gaussian fallback is used. Fits and background estimation may use the
    validity MASK; CMASK is not applied to the cleaned image. No pixel array,
    WCS, uncertainty, aperture, or calibration value is modified.
    """

    def __init__(self, image_collection: ImageInput, *, overwrite: bool = False) -> None:
        children(image_collection)
        self._image_collection: ImageInput = image_collection
        self._overwrite: bool = overwrite
        self._frames: list[_Frame] = []
        self._stars: list[_Star] = []
        self._failed_paths: dict[Path, str] = {}

    @property
    def image_collection(self) -> ImageInput:
        """Return the latest working selection with refreshed header summaries."""
        return self._image_collection

    @property
    def results(self) -> Table:
        """Return per-image clipped means, stellar scatter, and retained counts.

        scatter describes variation between accepted stars, not uncertainty of
        the mean. Failed frames are listed separately in failed_paths.
        """
        table = Table(rows=[asdict(frame) for frame in self._frames])
        if self._frames:
            table['fwhm'].unit = table['scatter'].unit = 'pix'
        return table

    @property
    def star_measurements(self) -> Table:
        """Return candidate positions, fit diagnostics, acceptance, and reasons."""
        table = Table(rows=[asdict(star) for star in self._stars])
        if self._stars:
            for name in ('x', 'y', 'fwhm', 'gamma'):
                table[name].unit = 'pix'
        return table

    @property
    def failed_paths(self) -> dict[Path, str]:
        """Return failures; their FITS headers have not been updated this run."""
        return self._failed_paths.copy()

    @staticmethod
    def _positive(value: float, name: str) -> float:
        if isinstance(value, bool) or not isinstance(value, Real) or not np.isfinite(value) or value <= 0:
            raise ValueError(f'{name} must be positive and finite.')
        return float(value)

    @staticmethod
    def _initial_width(data: np.ndarray, center: tuple[float, float], invalid: np.ndarray) -> float:
        """Reuse radial half-height estimation as a seed, interpolating the crossing."""
        edges = np.arange(0., min(data.shape) / 2, 0.5)
        profile = RadialProfile(data, center, edges, mask=invalid)
        radius, values = profile.radius, profile.profile
        finite = np.isfinite(values)
        radius, values = radius[finite], values[finite]
        if len(values) < 3 or values[0] <= 0:
            return 4.
        half = values[0] / 2
        crossings = np.flatnonzero(values[1:] <= half)
        if not len(crossings):
            return 4.
        index = int(crossings[0] + 1)
        r0, r1 = radius[index - 1:index + 1]
        v0, v1 = values[index - 1:index + 1]
        return float(2 * (r0 + (half - v0) * (r1 - r0) / (v1 - v0)))

    @staticmethod
    def _fit_star(
        data: np.ndarray, signal: np.ndarray, invalid: np.ndarray,
        x: float, y: float, *, image: Path, half_box: int,
        sky_inner: float, sky_outer: float, max_axis_ratio: float,
        saturation_adu: float | None,
    ) -> _Star:
        """Fit one small stamp after robust local-annulus background subtraction."""
        star = _Star(str(image), x, y)
        xi, yi = int(round(x)), int(round(y))
        radius = int(np.ceil(sky_outer))
        sky = signal[yi - radius:yi + radius + 1, xi - radius:xi + radius + 1]
        sky_invalid = invalid[yi - radius:yi + radius + 1, xi - radius:xi + radius + 1]
        sy, sx = np.indices(sky.shape, dtype=float)
        distance = np.hypot(sx + xi - radius - x, sy + yi - radius - y)
        annulus = (distance >= sky_inner) & (distance <= sky_outer) & ~sky_invalid
        if np.count_nonzero(annulus) < 30:
            return replace(star, reason='insufficient_background')
        _, background, noise = sigma_clipped_stats(sky[annulus], sigma=3, maxiters=5)
        if not np.all(np.isfinite([background, noise])) or noise <= 0:
            return replace(star, reason='invalid_background')

        ys = slice(yi - half_box, yi + half_box + 1)
        xs = slice(xi - half_box, xi + half_box + 1)
        stamp = signal[ys, xs] - background
        valid = ~invalid[ys, xs]
        if np.mean(valid) < 0.95:
            return replace(star, reason='invalid_stamp')
        if saturation_adu is not None and np.any(data[ys, xs][valid] >= saturation_adu):
            return replace(star, reason='saturated')
        yy, xx = np.indices(stamp.shape, dtype=float)
        xx += xi - half_box
        yy += yi - half_box
        peak = float(np.max(stamp[valid]))
        if peak <= 0:
            return replace(star, reason='nonpositive_peak')

        width = FWHM._initial_width(
            stamp, (x - xi + half_box, y - yi + half_box), ~valid,
        )
        width = float(np.clip(width, 2., half_box / 2))
        gamma = width / (2 * np.sqrt(2**(1 / 3.) - 1))
        initial = [peak, x, y, gamma, 3., 0.]
        lower = [0., x - 2.5, y - 2.5, 0.4, 1.1, -5 * noise]
        upper = [np.inf, x + 2.5, y + 2.5, float(half_box), 15., 5 * noise]

        def model(parameters: np.ndarray) -> np.ndarray:
            amplitude, cx, cy, scale, beta, offset = parameters
            rr = ((xx - cx)**2 + (yy - cy)**2) / scale**2
            return offset + amplitude * (1 + rr)**(-beta)

        def residual(parameters: np.ndarray) -> np.ndarray:
            return (model(parameters)[valid] - stamp[valid]) / noise

        fit = least_squares(
            residual, initial, bounds=(lower, upper), loss='soft_l1',
            x_scale='jac', max_nfev=500,
        )
        if not fit.success or not np.all(np.isfinite(fit.x)):
            return replace(star, reason='fit_failed')
        amplitude, cx, cy, gamma, beta, offset = map(float, fit.x)
        fwhm = 2 * gamma * np.sqrt(np.expm1(np.log(2) / beta))
        star = replace(
            star, x=cx, y=cy, fwhm=float(fwhm), gamma=gamma, beta=beta,
            peak=amplitude, background=float(background + offset),
        )
        # Reject bound-limited/truncated fits instead of accepting a numerical
        # answer just because the optimizer reported convergence.
        if np.any(fit.active_mask) or beta <= 1.101 or beta >= 14.99:
            return replace(star, reason='parameter_bound')
        if not 1. <= fwhm < half_box:
            return replace(star, reason='width_outside_stamp')
        if np.hypot(cx - x, cy - y) > 2.:
            return replace(star, reason='centroid_shift')

        rr = np.hypot(xx - cx, yy - cy)
        core = (rr <= min(1.5 * fwhm, half_box)) & valid
        rms = float(np.sqrt(np.mean((model(fit.x)[core] - stamp[core])**2)))
        star = replace(star, residual_rms=rms)
        if rms > max(3 * noise, 0.05 * amplitude):
            return replace(star, reason='poor_residual')

        # Positive high-significance core pixels give an independent elongation
        # check; these moments are diagnostics, not the reported Moffat width.
        weights = np.where((rr <= fwhm) & valid & (stamp > 3 * noise), stamp, 0.)
        total = float(np.sum(weights))
        if total <= 0 or np.count_nonzero(weights) < 6:
            return replace(star, reason='insufficient_core')
        mx, my = float(np.sum(weights * xx) / total), float(np.sum(weights * yy) / total)
        dx, dy = xx - mx, yy - my
        covariance = np.array([
            [np.sum(weights * dx * dx), np.sum(weights * dx * dy)],
            [np.sum(weights * dx * dy), np.sum(weights * dy * dy)],
        ]) / total
        eigenvalues = np.linalg.eigvalsh(covariance)
        if eigenvalues[0] <= 0:
            return replace(star, reason='invalid_shape')
        axis_ratio = float(np.sqrt(eigenvalues[1] / eigenvalues[0]))
        star = replace(star, axis_ratio=axis_ratio)
        if axis_ratio > max_axis_ratio:
            return replace(star, reason='elongated')
        return replace(star, accepted=True, reason='accepted')

    @staticmethod
    def _measure_image(
        path: Path, *, detect_threshold: float, fit_box_size: int,
        sky_inner: float, sky_outer: float, isolation_radius: float,
        max_sources: int, min_sources: int, clip_sigma: float,
        saturation_adu: float | None, max_axis_ratio: float,
    ) -> tuple[_Frame | None, list[_Star], str]:
        """Read once, select isolated peaks, fit stars, and clip within this image."""
        with fits.open(path, memmap=False) as hdul:
            data = np.ascontiguousarray(hdul[0].data, dtype=float)
            header = hdul[0].header
            invalid = ~np.isfinite(data)
            if 'MASK' in hdul:
                flags = hdul['MASK'].data
                if flags is None or flags.shape != data.shape:
                    raise ValueError('MASK shape does not match the primary image.')
                invalid |= np.asarray(flags, dtype=bool)
        if np.mean(~invalid) < 0.5:
            return None, [], 'Insufficient valid image area.'
        # SEP handles broad sky gradients; each candidate still gets an
        # independent local-annulus residual background estimate below.
        safe_data = data.copy()
        safe_data[invalid] = np.median(data[~invalid])
        background = sep.Background(safe_data, mask=invalid)
        signal = safe_data - background.back()
        noise = background.rms()
        if not np.all(np.isfinite(noise)) or np.max(noise) <= 0:
            return None, [], 'Cannot estimate positive background noise.'
        border = int(np.ceil(sky_outer)) + 3
        if min(data.shape) <= 2 * border:
            return None, [], 'Image too small for the background annulus.'
        peaks = find_peaks(
            signal, threshold=detect_threshold * noise, box_size=5,
            mask=invalid, border_width=border,
        )
        if peaks is None or not len(peaks):
            return None, [], 'No bright candidate sources detected.'
        coordinates = np.column_stack((peaks['x_peak'], peaks['y_peak'])).astype(float)
        distances = cKDTree(coordinates).query(coordinates, k=2)[0][:, 1]
        order = np.argsort(np.asarray(peaks['peak_value']))[::-1]
        centroids = centroid_sources(
            signal, coordinates[:, 0], coordinates[:, 1], box_size=7,
            centroid_func=centroid_com, mask=invalid,
        )
        stars: list[_Star] = []
        accepted = 0
        for index in order[:max_sources * 10]:
            x, y = float(centroids[0][index]), float(centroids[1][index])
            if distances[index] < isolation_radius:
                stars.append(_Star(str(path), *coordinates[index], reason='neighbor'))
                continue
            if not np.all(np.isfinite([x, y])) or np.hypot(x - coordinates[index, 0], y - coordinates[index, 1]) > 2.5:
                stars.append(_Star(str(path), *coordinates[index], reason='invalid_centroid'))
                continue
            try:
                star = FWHM._fit_star(
                    data, signal, invalid, x, y, image=path,
                    half_box=fit_box_size // 2, sky_inner=sky_inner,
                    sky_outer=sky_outer, max_axis_ratio=max_axis_ratio,
                    saturation_adu=saturation_adu,
                )
            except (ValueError, RuntimeError, FloatingPointError) as exc:
                star = _Star(str(path), x, y, reason=f'fit_error: {exc}')
            stars.append(star)
            accepted += int(star.accepted)
            if accepted >= max_sources:
                break

        good = [index for index, star in enumerate(stars) if star.accepted]
        if len(good) < min_sources:
            return None, stars, f'Only {len(good)} reliable stars; require {min_sources}.'
        widths = np.array([stars[index].fwhm for index in good])
        clipped = sigma_clip(widths, sigma=clip_sigma, maxiters=5, cenfunc='median', stdfunc=mad_std)
        for index, rejected in zip(good, np.ma.getmaskarray(clipped)):
            if rejected:
                stars[index] = replace(stars[index], accepted=False, reason='fwhm_outlier')
        retained = clipped.compressed()
        if len(retained) < min_sources:
            return None, stars, f'Only {len(retained)} stars remain after clipping; require {min_sources}.'
        frame = _Frame(
            str(path), str(header['FIELD']), str(header['FILTER']),
            float(np.mean(retained)), float(np.std(retained, ddof=1)),
            len(retained), len(stars),
        )
        return frame, stars, ''

    def find_fwhm(
        self, *, image_type: str = 'Master Light', detect_threshold: float = 10,
        fit_box_size: int = 41, background_inner_radius: float = 25,
        background_outer_radius: float = 35, isolation_radius: float | None = None,
        max_sources: int = 50, min_sources: int = 5, clip_sigma: float = 3,
        saturation_adu: float | None = None, max_axis_ratio: float = 1.5,
        write_fwhm: bool = True, overwrite: bool | None = None,
    ) -> ImageInput:
        """Measure each image independently and optionally update its FWHM header.

        Parameters
        ----------
        image_type : str, optional
            IMTYPE selection, default Master Light.
        detect_threshold : float, optional
            Local background-RMS multiple for unsmoothed peak detection. Default
            10 retains the old FWHM candidate threshold; it is not aperture S/N
            or Astrometry.net's smoothed detection significance.
        fit_box_size : int, optional
            Odd Moffat fitting-stamp width in pixels, default 41. Fitted full
            widths must be at least 1 pixel and smaller than half this stamp.
            Candidate centroiding uses a separate 7-pixel box, not this stamp.
        background_inner_radius, background_outer_radius : float, optional
            Local background-annulus radii in pixels, defaults 25 and 35.
            Inner radius must exceed the stamp half-width; outer must be larger.
        isolation_radius : float or None, optional
            Minimum separation from other detected peaks, in pixels. None uses
            fit_box_size. Fainter undetected neighbors can still contaminate a fit.
        max_sources, min_sources : int, optional
            Maximum successful fits to collect and minimum retained after
            clipping, defaults 50 and 5. At most 10*max_sources bright candidates
            are considered; actual rejection reasons are retained for inspection.
        clip_sigma : float, optional
            Within-image FWHM clipping threshold, default 3, median/MAD center
            and scale, at most five iterations. No between-image clipping.
        saturation_adu : float or None, optional
            Optional ceiling in the current image's ADU values, before local
            sky subtraction. None imposes no numerical saturation ceiling.
            Raw CCD saturation (e.g. 65535) cannot blindly be transferred to
            calibrated, averaged frames; saturation may already be diluted.
        max_axis_ratio : float, optional
            Core moment major/minor ratio limit, default 1.5; must exceed 1.
            This rejects elongation, not a guarantee that every object is stellar.
        write_fwhm : bool, optional
            Default True writes FWHM (full width in pixels), FWHMSTD (between-star
            sample scatter), FWHMN (retained stars), and FWHMMET='Moffat'.
            False measures without changing any FITS file.
        overwrite : bool or None, optional
            Header-update permission; None uses the constructor policy.

        Returns
        -------
        ImageFileCollection or ImageCollectionGroup
            Updated working selection. Results and star_measurements hold
            diagnostics; failed_paths records images with unchanged headers.

        Notes
        -----
        Numerical bounds beta=(1.1,15), gamma=(0.4,stamp half-width) prevent
        unconstrained wing/width solutions; bound-limited fits are rejected.
        Positive amplitude, bounded center and residual-background fitting use
        soft-L1 least squares scaled by the local sky RMS. A residual check
        rejects RMS above max(3*sky RMS, 5% of fitted peak) inside 1.5 FWHM.
        These are quality controls, not uncertainty estimates. Wider PSFs need
        larger stamps and annuli. Failed frames retain any pre-existing FWHM;
        inspect failed_paths before using old headers downstream. No Gaussian
        fallback or replacement by a nightly average is performed.
        """
        for name, value in (
            ('detect_threshold', detect_threshold), ('background_inner_radius', background_inner_radius),
            ('background_outer_radius', background_outer_radius), ('clip_sigma', clip_sigma),
            ('max_axis_ratio', max_axis_ratio),
        ):
            self._positive(value, name)
        for name, value in (('fit_box_size', fit_box_size), ('max_sources', max_sources), ('min_sources', min_sources)):
            if isinstance(value, bool) or not isinstance(value, Integral):
                raise ValueError(f'{name} must be an integer.')
        if fit_box_size < 11 or fit_box_size % 2 == 0:
            raise ValueError('fit_box_size must be odd and at least 11.')
        if not 3 <= min_sources <= max_sources:
            raise ValueError('Require 3 <= min_sources <= max_sources.')
        if not fit_box_size // 2 < background_inner_radius < background_outer_radius:
            raise ValueError('Background annulus must lie outside the fitting stamp.')
        if max_axis_ratio <= 1:
            raise ValueError('max_axis_ratio must exceed 1.')
        isolation_radius = float(fit_box_size) if isolation_radius is None else self._positive(isolation_radius, 'isolation_radius')
        if saturation_adu is not None:
            self._positive(saturation_adu, 'saturation_adu')
        replace_headers = self._overwrite if overwrite is None else overwrite
        if write_fwhm and not replace_headers:
            raise ValueError('Writing FWHM headers requires overwrite=True; use write_fwhm=False for measurement only.')
        refresh(self._image_collection)
        selected = [path for child in children(self._image_collection) for path in select(child, IMTYPE=image_type)]
        if not selected:
            raise ValueError(f'No {image_type} images selected for FWHM measurement.')
        for path in selected:
            header = fits.getheader(path)
            if header.get('NAXIS') != 2 or not header.get('FIELD') or not header.get('FILTER'):
                raise ValueError(f'{path.name}: require a 2D image and FIELD/FILTER headers.')
            if str(header.get('BUNIT', '')).strip().lower() != 'adu':
                raise ValueError(f'{path.name}: expected ADU image data.')
        self._frames, self._stars, self._failed_paths = [], [], {}
        with tqdm(selected, desc='Measuring FWHM', unit='image') as progress:
            for path in progress:
                progress.set_postfix_str(path.name)
                try:
                    frame, stars, error = self._measure_image(
                        path, detect_threshold=detect_threshold, fit_box_size=fit_box_size,
                        sky_inner=background_inner_radius, sky_outer=background_outer_radius,
                        isolation_radius=isolation_radius, max_sources=max_sources,
                        min_sources=min_sources, clip_sigma=clip_sigma,
                        saturation_adu=saturation_adu, max_axis_ratio=max_axis_ratio,
                    )
                    self._stars.extend(stars)
                    if frame is None:
                        raise ValueError(error)
                    if write_fwhm:
                        update_fits_header(path, {
                            'FWHM': (frame.fwhm, 'Clipped mean circular Moffat full width [pixel]'),
                            'FWHMSTD': (frame.scatter, 'Between-star FWHM sample scatter [pixel]'),
                            'FWHMN': (frame.nstars, 'Stars retained for FWHM measurement'),
                            'FWHMMET': ('Moffat', 'Circular 2D Moffat; per-image clipped mean'),
                        }, overwrite=True)
                except (OSError, ValueError, RuntimeError) as exc:
                    self._failed_paths[path] = str(exc)
                    logger.error('%s: FWHM measurement failed: %s', path.name, exc)
                else:
                    self._frames.append(frame)
                    logger.info('%s: FWHM=%.3f pix; scatter=%.3f pix; %d stars.', path.name, frame.fwhm, frame.scatter, frame.nstars)
        self._image_collection = rebuild(self._image_collection)
        return self._image_collection

    def plot_fwhm(self) -> Figure:
        """Compare successful images by field/filter; return a notebook-ready figure.

        Points show per-image clipped means. Error bars show between-star sample
        scatter, not uncertainty of the mean; labels give retained star counts.
        There is no time-trend line or nightly averaging. Failed images are not
        plotted; inspect failed_paths. No file is saved or figure shown here.
        """
        if not self._frames:
            raise ValueError('No successful FWHM measurements to plot.')
        fields = list(dict.fromkeys(frame.field for frame in self._frames))
        filters = list(dict.fromkeys(frame.filter for frame in self._frames))
        figure, axes = plt.subplots(figsize=(max(7., 2.2 * len(fields)), 4.8), constrained_layout=True)
        offsets = np.linspace(-0.25, 0.25, len(filters)) if len(filters) > 1 else [0.]
        for index, band in enumerate(filters):
            frames = [frame for frame in self._frames if frame.filter == band]
            x = [fields.index(frame.field) + offsets[index] for frame in frames]
            axes.errorbar(x, [frame.fwhm for frame in frames], yerr=[frame.scatter for frame in frames],
                          fmt='o', capsize=4, label=band, color=f'C{index % 10}')
            for position, frame in zip(x, frames):
                axes.annotate(f'n={frame.nstars}', (position, frame.fwhm), xytext=(0, 8), textcoords='offset points', ha='center', fontsize=8)
        axes.set_xticks(range(len(fields)), fields)
        axes.set_xlabel('Observation field')
        axes.set_ylabel('Moffat FWHM [pixel]')
        axes.set_title('Combined-image FWHM; bars show between-star scatter')
        axes.grid(axis='y', alpha=0.25)
        axes.legend(title='Filter')
        return figure
