"""Snapshots of the actual local-sky masks used for a measurement."""

from dataclasses import dataclass
from pathlib import Path
import json

import numpy as np
from astropy.io import fits


@dataclass
class BackgroundDiagnostic:
    """One primary annulus snapshot, with distinct mean and area selections.

    Labels: 0 outside annulus, 1 used, 2 invalid, 3 source mask, 4 sigma clipped.
    Mean statistics use pixel-center inclusion. Effective areas instead use
    exact fractional overlap and their own clipping mask, as in ApertureStats.
    These are diagnostic copies; editing them cannot change measured results.
    """

    metadata: dict
    image: np.ndarray
    mean_labels: np.ndarray
    area_labels: np.ndarray
    area_weights: np.ndarray
    detected: np.ndarray
    source_mask: np.ndarray
    areas: dict

    @classmethod
    def capture(cls, estimator, aperture, annulus, stats, growth, metadata):
        exact = annulus.to_mask(method='exact')
        large, small = exact.get_overlap_slices(estimator.data.shape)
        if large is None:
            return None
        weights = exact.data[small].copy()
        center = annulus.to_mask(method='center').data[small] > 0
        invalid = estimator.invalid[large]
        detected = np.zeros(weights.shape, bool)
        masked = np.zeros(weights.shape, bool)
        if estimator.distance is not None:
            detected = estimator.distance[large] == 0
            masked = estimator.distance[large] <= growth * estimator.fwhm

        def labels(inside, final_mask):
            result = np.zeros(weights.shape, np.uint8)
            result[inside] = 1
            result[inside & invalid] = 2
            result[inside & ~invalid & masked] = 3
            result[inside & ~invalid & ~masked & final_mask] = 4
            return result

        mean_labels = labels(center, np.ma.getmaskarray(stats.data_cutout))
        area_labels = labels(weights > 0, np.ma.getmaskarray(stats.data_sumcutout))
        after_valid = float(np.sum(weights[~invalid]))
        after_sources = float(np.sum(weights[~invalid & ~masked]))
        final_area = float(np.sum(weights[area_labels == 1]))
        areas = {
            'geometric': float(annulus.area),
            'on_image': float(weights.sum()),
            'after_validity': after_valid,
            'after_source_mask': after_sources,
            'after_clipping': final_area,
        }
        # Enforce the link between visualization and the actual measurement.
        measured_area = float(stats.sum_aper_area.value)
        if np.isfinite(measured_area) and not np.isclose(final_area, measured_area):
            raise RuntimeError('Background diagnostic area disagrees with photometry.')
        details = {
            **metadata,
            'x_origin': large[1].start, 'y_origin': large[0].start,
            'aperture_x': float(aperture.positions[0]),
            'aperture_y': float(aperture.positions[1]),
            'aperture_radius': float(aperture.r),
            'background_x': float(annulus.positions[0]),
            'background_y': float(annulus.positions[1]),
            'inner_radius': float(annulus.r_in), 'outer_radius': float(annulus.r_out),
            'fwhm': estimator.fwhm, 'growth_fwhm': growth,
            'growth_pixels': growth * estimator.fwhm if estimator.distance is not None else 0.,
            'mean': float(stats.mean),
            'mean_pixel_count': int(np.count_nonzero(mean_labels == 1)),
            'options': estimator.options.as_dict(),
        }
        return cls(details, estimator.data[large].copy(), mean_labels, area_labels,
                   weights, detected.copy(), masked.copy(), areas)

    def _draw_panels(self, axes):
        """Draw one snapshot on supplied axes, shared by single/batch views."""
        from matplotlib.colors import ListedColormap, BoundaryNorm
        from matplotlib.patches import Circle
        from astropy.visualization import simple_norm

        m = self.metadata
        extent = (m['x_origin'] - .5, m['x_origin'] + self.image.shape[1] - .5,
                  m['y_origin'] - .5, m['y_origin'] + self.image.shape[0] - .5)
        finite = self.image[np.isfinite(self.image)]
        norm = simple_norm(finite, 'asinh', percent=99) if finite.size else None
        colors = ['#f2f2f2', '#70bb88', '#747474', '#ef9b45', '#9d77b5']
        cmap = ListedColormap(colors)
        bounds = BoundaryNorm(np.arange(-.5, 5.5), cmap.N)
        axes[0].imshow(self.image, origin='lower', extent=extent, cmap='gray', norm=norm)
        axes[0].set_title('Image + reviewed regions')
        for ax, data, title in zip(axes[1:], [self.mean_labels, self.area_labels],
                                   ['Pixels used for sky mean', 'Pixels used for exact area']):
            ax.imshow(data, origin='lower', extent=extent, cmap=cmap, norm=bounds,
                      interpolation='nearest')
            ax.set_title(title)
        for ax in axes:
            ax.add_patch(Circle((m['aperture_x'], m['aperture_y']), m['aperture_radius'],
                                fill=False, edgecolor='#08c838', linewidth=1.3))
            for radius in [m['inner_radius'], m['outer_radius']]:
                ax.add_patch(Circle((m['background_x'], m['background_y']), radius,
                                    fill=False, edgecolor='#008bc4', linewidth=1.3))
            ax.set_xlim(extent[:2])
            ax.set_ylim(extent[2:])
            ax.set_xlabel('x [zero-based image pixel]')
        axes[0].set_ylabel('y [zero-based image pixel]')
        return colors

    def _caption(self):
        m = self.metadata
        a = self.areas
        text = (f"Exact area [pixel²]: geometric {a['geometric']:.1f} → on image {a['on_image']:.1f}"
                f" → valid {a['after_validity']:.1f} → source-masked {a['after_source_mask']:.1f}"
                f" → clipped {a['after_clipping']:.1f}\n"
                f"Final effective fraction: {a['after_clipping']/a['geometric']:.1%}; "
                f"sky mean uses {m['mean_pixel_count']} pixel centers.\n"
                f"SEP detection: {m['options']['detection_sigma']:g}σ, min area {m['options']['min_area']} pixels; "
                f"clipping: {m.get('clip_sigma', 3):g}σ.\n"
                "Green circle: source aperture; blue circles: sky annulus. "
                "Mean selection and exact-area selection can differ at boundaries and after clipping.")
        return text

    def plot(self):
        """Return one filter's image, masks and area stages."""
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch

        m = self.metadata
        fig, axes = plt.subplots(1, 3, figsize=(14, 6.4))
        fig.subplots_adjust(left=.055, right=.985, bottom=.30, top=.82, wspace=.20)
        colors = self._draw_panels(axes)
        fig.suptitle(
            f"{m['source_id']} | {m['field']} | {m['filter']}\n"
            f"Sky mean = {m['mean']:.5f} ADU/pixel; "
            f"mask growth = {m['growth_fwhm']:g} FWHM ({m['growth_pixels']:.2f} pixels)",
            fontsize=13,
        )
        fig.legend(handles=[Patch(facecolor=colors[i], label=label) for i, label in
                            [(1, 'Used'), (2, 'Invalid'), (3, 'Source mask'), (4, 'Sigma clipped')]],
                   loc='lower center', bbox_to_anchor=(.5, .20), ncol=4, frameon=False)
        fig.text(.5, .095, self._caption(), ha='center', va='center', fontsize=9)
        return fig

    def save(self, path, *, overwrite=False):
        """Export the snapshot to FITS for pixel-level inspection after restart."""
        m = self.metadata
        header = fits.Header()
        header['DIAGVER'] = 1
        for key, value in [('SOURCE', m['source_id']), ('FIELD', m['field']),
                           ('FILTER', m['filter']), ('XORIGIN', m['x_origin']),
                           ('YORIGIN', m['y_origin']), ('GROWPIX', m['growth_pixels'])]:
            header[key] = value
        for index, name in enumerate(['Outside annulus', 'Used', 'Invalid', 'Source mask', 'Sigma clipped']):
            header[f'LABEL{index}'] = name

        def serializable(value):
            if isinstance(value, dict):
                return {k: serializable(v) for k, v in value.items()}
            if isinstance(value, (tuple, list)):
                return [serializable(v) for v in value]
            if isinstance(value, float) and not np.isfinite(value):
                return None
            return value

        payload = json.dumps(serializable({'metadata': m, 'areas': self.areas}), allow_nan=False)
        hdus = fits.HDUList([
            fits.PrimaryHDU(self.image, header),
            fits.ImageHDU(self.mean_labels, name='MEANLABEL'),
            fits.ImageHDU(self.area_labels, name='AREALABEL'),
            fits.ImageHDU(self.area_weights, name='AREAWEIGHT'),
            fits.ImageHDU(self.detected.astype(np.uint8), name='DETECTED'),
            fits.ImageHDU(self.source_mask.astype(np.uint8), name='SOURCEMASK'),
            fits.BinTableHDU.from_columns([
                fits.Column(name='JSON', format=f'{len(payload)}A', array=[payload]),
            ], name='METADATA'),
        ])
        hdus.writeto(Path(path), overwrite=overwrite)

    @classmethod
    def load(cls, path):
        """Restore an exported snapshot without accessing source images or regions."""
        with fits.open(path, memmap=False) as hdul:
            if hdul[0].header.get('DIAGVER') != 1:
                raise ValueError('Unsupported background diagnostic version.')
            payload = json.loads(hdul['METADATA'].data['JSON'][0])
            arrays = [np.array(hdul[name].data, copy=True) for name in
                      ['PRIMARY', 'MEANLABEL', 'AREALABEL', 'AREAWEIGHT', 'DETECTED', 'SOURCEMASK']]
        if any(a.shape != arrays[0].shape for a in arrays) or arrays[0].ndim != 2:
            raise ValueError('Incompatible diagnostic array shapes.')
        for a in arrays[1:3]:
            if not np.all(np.isin(a, [0, 1, 2, 3, 4])):
                raise ValueError('Unknown diagnostic pixel label.')
        arrays[4] = arrays[4].astype(bool)
        arrays[5] = arrays[5].astype(bool)
        # JSON null represents an unavailable statistic in an invalid annulus.
        if payload['metadata']['mean'] is None:
            payload['metadata']['mean'] = np.nan
        from .local_background import LocalBackgroundOptions
        payload['metadata']['options'] = LocalBackgroundOptions(**payload['metadata']['options']).as_dict()
        return cls(payload['metadata'], *arrays, payload['areas'])


def plot_background_collection(diagnostics):
    """Draw every available filter of one source in one figure, using snapshots."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    count = len(diagnostics)
    fig = plt.figure(figsize=(14, 5.0 * count + .8))
    grid = fig.add_gridspec(
        2 * count, 3, height_ratios=[3.6, .9] * count,
        left=.055, right=.985, bottom=.055 / count, top=1 - 1.0 / (5 * count + .8),
        hspace=.42, wspace=.20,
    )
    for index, diagnostic in enumerate(diagnostics):
        axes = [fig.add_subplot(grid[2 * index, column]) for column in range(3)]
        colors = diagnostic._draw_panels(axes)
        m = diagnostic.metadata
        axes[0].set_title(f"{m['filter']} — Image + reviewed regions")
        caption = fig.add_subplot(grid[2 * index + 1, :])
        caption.axis('off')
        caption.text(.5, .5,
                     f"{m['filter']} | sky = {m['mean']:.5f} ADU/pixel | "
                     f"growth = {m['growth_fwhm']:g} FWHM ({m['growth_pixels']:.2f} pixels)\n"
                     + diagnostic._caption(), ha='center', va='center', fontsize=8.5)
    first = diagnostics[0].metadata
    fig.suptitle(f"{first['source_id']} | Field: {first['field']} | All available background masks",
                 fontsize=14, y=1 - .2 / (5 * count + .8))
    fig.legend(handles=[Patch(facecolor=colors[i], label=label) for i, label in
                        [(1, 'Used'), (2, 'Invalid'), (3, 'Source mask'), (4, 'Sigma clipped')]],
               loc='lower center', ncol=4, frameon=False)
    return fig
