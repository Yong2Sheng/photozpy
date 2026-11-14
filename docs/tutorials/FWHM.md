# Measure FWHM in combined images

Run this after Master Light combination and local plate solving. WCS is retained
but is not required by the pixel-based width measurement. Each field/filter
product gets its own FWHM; there is no nightly averaging.

```python
from photozpy import FWHM

# Initialize independent FWHM measurements for the combined images.
fwhm_measurement: FWHM = FWHM(
    image_collection=analysis_images,
    overwrite=True,
)

# Fit bright, isolated stars and save each Master Light's own FWHM.
analysis_images: ImageCollectionGroup = fwhm_measurement.find_fwhm()

# Inspect per-image widths, stellar scatter, and retained star counts.
fwhm_measurement.results

# Check failures before using FWHM to generate photometric apertures.
fwhm_measurement.failed_paths

# Compare fields and filters; bars show between-star scatter, not mean errors.
figure = fwhm_measurement.plot_fwhm()
figure
```

For measurement without any FITS changes, use `find_fwhm(write_fwhm=False)`.
`plot_fwhm()` returns a Matplotlib figure and does not save files automatically.
Use `figure.savefig(...)` if you want an exported scientific plot.

## What is retained from the original code?

The original workflow detects local peaks, refines centroids, measures individual
stars, and sigma-clips their widths. This implementation retains Photutils
`find_peaks`, `centroid_sources`/`centroid_com`, and the radial-profile half-height
idea as a fit initializer. The final measurement is now a bounded circular 2D
Moffat fit, rather than a discretely sampled radial-profile crossing.

The existing `fwhm.py` placeholder is now the measurement module. Superseded
`SourceDetection.find_fwhm`, `find_fwhm_from_image`, `get_radial_bkg`, and
`get_mad_clipped_average` were removed from `source_detection.py`, so there is
only one FWHM implementation. General legacy peak detection and centroid
functions remain for their later refactoring. The photometry package
exports the refactored FWHM interface without importing the legacy Photometry
and SourceDetection implementations. The remaining legacy functions stay
available for later refactoring; the aperture-generation and photometry
algorithms are not changed. The historical `docs/Tutorial.ipynb` still contains
an old `SourceDetection.find_fwhm` example; use the cell above instead.

The new measurement uses ImageFileCollection/ImageCollectionGroup helpers and
the atomic FITS-header utility, not legacy CollectionManager or Sources.
It works with the existing Photutils 1.13 environment; no upgrade is required.

## Selection and background

1. Read each Master Light once, excluding non-finite and validity-MASK pixels.
   CMASK is not used to exclude cleaned values. No array in the file is changed.
2. SEP estimates a broad background map and local RMS; Photutils finds peaks
   above `detect_threshold * RMS`, default 10 as in the old FWHM candidate code.
   This unsmoothed threshold is not aperture S/N or Astrometry.net's smoothed
   source-detection significance.
3. Search only far enough from the edge to fit a full background annulus.
   Rank candidates by background-subtracted peak brightness. Reject peaks with
   another detected peak inside `isolation_radius` (default 41 pixels).
4. Refine initial centers in a separate 7-pixel window. Fit in a 41-pixel stamp;
   a 301-pixel centroid window is not reused.
5. Measure a sigma-clipped residual sky level/RMS in an annulus of radii 25 and
   35 pixels after broad background subtraction. The Moffat fit also permits a
   small residual constant background offset.

These defaults are starting values, not a guarantee of isolation or optimality
for every image. Neighbors below the candidate detection threshold, round
galaxies, and saturation diluted by stacking can escape automatic selection.
Inspect the candidate table and actual images if results look unexpected.

`saturation_adu=None` does not impose an arbitrary ceiling. Set it only if you
know the saturation threshold in these processed image values. Raw 16-bit limits
cannot automatically be applied to bias/dark/flat-corrected, averaged pixels.
An explicit ceiling rejects any fitting stamp containing a value at or above it.

## Moffat definition and quality checks

The fitted circular profile is:

```text
I(r) = background + amplitude * (1 + r^2 / gamma^2)^(-beta)
FWHM = 2 * gamma * sqrt(2^(1/beta) - 1)
```

FWHM is a full width in pixels. It is not a radius, a major/minor-axis average,
or a measurement of atmospheric seeing alone. Combined images also reflect
focus, tracking, registration interpolation, and stacking.

The center, positive amplitude, gamma, beta, and residual background are fitted
with bounded soft-L1 least squares. Beta is free within (1.1, 15), initialized at
3; it is not fixed to a theoretical atmospheric value. Gamma ranges from 0.4
pixels to the stamp half-width. A fit that reaches a numerical bound is rejected.
There is no automatic Gaussian fallback.

Additional checks reject failed convergence, non-finite parameters, full widths
below 1 pixel or beyond the stamp half-width, center displacement above 2 pixels,
and large residuals. Inside 1.5 FWHM, residual RMS must not exceed the larger of
three sky RMS and five percent of the fitted amplitude. Positive core moments
provide an independent elongation diagnostic; major/minor ratio above 1.5 is
rejected. These thresholds are quality heuristics, not statistical confidence
limits or photometric uncertainties.

For a much broader PSF, increase `fit_box_size` and both annulus radii explicitly;
the annulus must remain outside the fitting stamp. `max_axis_ratio` and detection,
isolation, and clipping parameters are exposed in `find_fwhm`.

At most 50 successful fits are collected, considering up to 500 bright
candidates. Median/MAD 3-sigma clipping is performed within each image, for at
most five iterations. The retained widths are averaged, following the original
mean-summary idea. At least five stars must remain by default; `min_sources`
can be deliberately changed but must be at least three.

## Results, failures, and downstream use

Successful FITS header updates are atomic and preserve image data, WCS, MASK,
CMASK, other HDUs, units, and existing checksums:

| Keyword | Meaning |
| --- | --- |
| FWHM | Clipped mean stellar full width, pixels |
| FWHMSTD | Sample standard deviation between retained stars, pixels |
| FWHMN | Number of retained stars |
| FWHMMET | Moffat |

The comparison plot uses observation fields as categories and colors for
filters. Error bars show between-star scatter, **not uncertainty on the mean**;
labels give retained star counts. Failed frames are not plotted. Do not interpret
it as a time trend or automatically clip genuinely broader images across fields.

`star_measurements` retains x/y positions, fitted gamma/beta, residual RMS,
axis ratio, acceptance, and reasons such as neighbor, saturated, parameter_bound,
poor_residual, elongated, or fwhm_outlier. `background` is the residual local
background after SEP map subtraction, not the full original sky level.

Failures are logged and placed in `failed_paths`, with no replacement FWHM
guessed. Any pre-existing FWHM in a failed image remains unchanged: check failure
records before downstream aperture generation, rather than mistaking an old
header value for a new successful measurement. Successful earlier images are
not rolled back when a later image fails. The returned collection retains both.

Later aperture generation should read each image's FWHM, using the existing
scaling (e.g. aperture **radius** = 3*FWHM). Background-annulus scaling of 5/8 and
photometric aperture correction are separate downstream choices, unchanged here.

Scientific references: [IRAF PSF measurement](https://iraf.readthedocs.io/en/latest/tasks/noao/obsutil/psfmeasure.html),
[Moffat PSF study](https://arxiv.org/abs/astro-ph/0109067), and
[Astropy Moffat2D](https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.Moffat2D.html).
