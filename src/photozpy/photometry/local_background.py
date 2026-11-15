"""Optional source masking and diagnostics for local annulus photometry."""

from dataclasses import asdict, dataclass
from numbers import Real

import numpy as np
from astropy.stats import SigmaClip
from photutils.aperture import ApertureStats
from scipy.ndimage import distance_transform_edt


@dataclass(frozen=True)
class LocalBackgroundOptions:
    """Controls for local sky only; source-aperture pixels remain unchanged.

    Growth is measured beyond SEP's detected source footprints, in each image's
    own header FWHM. Effective area includes validity, neighbor masks and sigma
    clipping. Stability compares induced source-count shifts to the existing
    count-error model; it is a diagnostic, not an empirical noise calibration.
    """

    mask_neighbors: bool = True
    detection_sigma: float = 3.0
    min_area: int = 5
    mask_growth_fwhm: float = 1.5
    min_area_fraction: float = 0.5
    min_pixels: int = 50
    stability_growth_fwhms: tuple[float, ...] = (1.0, 1.5, 2.0)
    stability_tolerance: float = 1.0

    def __post_init__(self):
        if not isinstance(self.mask_neighbors, bool):
            raise ValueError('mask_neighbors must be boolean.')
        for name in ('detection_sigma', 'mask_growth_fwhm', 'stability_tolerance'):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Real) or not np.isfinite(value) or value <= 0:
                raise ValueError(f'{name} must be positive and finite.')
        for name in ('min_area', 'min_pixels'):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ValueError(f'{name} must be a positive integer.')
        value = self.min_area_fraction
        if isinstance(value, bool) or not isinstance(value, Real) or not np.isfinite(value) or not 0 < value <= 1:
            raise ValueError('min_area_fraction must be in (0, 1].')
        values = self.stability_growth_fwhms
        if not isinstance(values, (tuple, list)):
            raise ValueError('stability_growth_fwhms must be a sequence of positive finite numbers.')
        if any(isinstance(v, bool) or not isinstance(v, Real) or not np.isfinite(v) or v <= 0 for v in values):
            raise ValueError('stability_growth_fwhms must contain positive finite numbers.')
        object.__setattr__(self, 'stability_growth_fwhms', tuple(values))

    def as_dict(self):
        """Serializable provenance stored with the measurement table."""
        return asdict(self)


class _LocalBackground:
    """Detect once per image, then compute independent local annulus statistics."""

    def __init__(self, data, invalid, fwhm, options):
        self.data = data
        self.invalid = invalid
        self.options = options
        self.fwhm = fwhm
        self.distance = None
        if options.mask_neighbors:
            if not isinstance(fwhm, Real) or not np.isfinite(fwhm) or fwhm <= 0:
                raise ValueError('Neighbor masking requires this image\'s positive finite FWHM in pixels.')
            if not np.any(~invalid):
                self.distance = np.full(data.shape, np.inf)
                return
            import sep

            image = np.array(data, dtype=np.float32, order='C', copy=True)
            image[invalid] = np.median(image[~invalid]) if np.any(~invalid) else 0
            background = sep.Background(image, mask=invalid)
            rms = max(float(background.globalrms), float(np.finfo(np.float32).eps))
            _, segmentation = sep.extract(
                image - background.back(), options.detection_sigma * rms,
                minarea=options.min_area, mask=invalid, segmentation_map=True,
            )
            detected = segmentation > 0
            self.distance = (
                distance_transform_edt(~detected) if np.any(detected)
                else np.full(data.shape, np.inf)
            )

    def statistics(self, annulus, clip_sigma, growth):
        mask = self.invalid
        if self.distance is not None:
            mask = mask | (self.distance <= growth * self.fwhm)
        stats = ApertureStats(
            self.data, annulus, mask=mask,
            sigma_clip=SigmaClip(sigma=clip_sigma, maxiters=10), sum_method='exact',
        )
        n_ann = float(stats.sum_aper_area.value)
        mean = float(stats.mean)
        fraction = n_ann / annulus.area
        valid = (
            np.isfinite(mean) and np.isfinite(n_ann)
            and n_ann >= self.options.min_pixels
            and fraction >= self.options.min_area_fraction
        )
        return stats, {
            'growth_fwhm': growth, 'mean': mean, 'n_ann': n_ann,
            'area_fraction': fraction, 'valid': bool(valid),
        }
