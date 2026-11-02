# Solve Master Light images locally

Run this after Registration and Light combination. The local Astrometry.net
solver detects sources itself; no FWHM measurement, Sources object, API key,
or upload is needed. Measure each Master Light's FWHM afterward, before creating
FWHM-based photometric apertures.

```python
from photozpy import PlateSolving

# Initialize local WCS solving for the combined working images.
plate_solving: PlateSolving = PlateSolving(
    image_collection=analysis_images,
    overwrite=True,
)

# Solve Master Light images locally and update their WCS headers in analysis/.
analysis_images: ImageCollectionGroup = plate_solving.add_wcs_locally(
    detect_threshold=8,
)

# Check failures before proceeding to photometry.
plate_solving.failed_paths
```

The executable is found automatically: first in PATH, then in the standard
Apple Silicon and Intel Homebrew locations. If it is unavailable, an error is
raised before solving. The optional executable argument is only needed for a
custom installation or selecting a particular version; an explicit choice is
never silently replaced by another installation. An executable and appropriate
local index files must already be installed/configured; this method does not
install or download anything.
Importing photozpy does not require the executable to be present.

## Preserved workflow and updated APIs

The original class and `add_wcs_locally` / `add_wcs_locally_for_file` entry points
remain. Each image is sent to the local `solve-field` program with `--no-plots`
and the original default detection significance of 8. There is no online
fallback or automatic parameter sweep. The online `add_wcs_astrometrynet` method
has been removed. Other modules may still use Astroquery for unrelated services.

The constructor now accepts ImageFileCollection or ImageCollectionGroup and
does not require the legacy Sources/CollectionManager interfaces. Batch calls
return an updated collection, retaining calibration frames and other members.
`save_location=None` writes in place. An explicit directory writes solved copies
there. `overwrite=True` is required to replace any existing destination.

Successful output paths are available in `solved_paths`. Per-image failures are
logged and recorded in `failed_paths`; other images continue, as in the original
loop. Failed originals remain in the returned collection. Do not assume all
images have WCS just because the batch returned. Missing executable, invalid
controls, path collisions, and malformed pointing metadata are rejected before
the batch starts. There is no transaction rolling back earlier successes.

The fixed progress-bar title is `Plate solving`; the current filename appears
at its end. The solver's captured stdout/stderr is included in failure messages.

## Pointing hints and units

Each image supplies its own RA/DEC header hints. With the default `ra_unit="auto"`:

- Numeric RA values and plain numeric strings mean degrees.
- Sexagesimal RA such as `07 13 59.02` or `07:13:59.02` means hours.
- DEC is degrees, including signed sexagesimal strings.

Use `ra_unit="hourangle"` explicitly for decimal-hour RA headers. Search radius
is in degrees (default 10). If either header hint is absent, the solver runs
without a position restriction. Invalid existing hints are not silently ignored.
The per-file method also accepts an explicit RA/DEC pair.

`detect_threshold` is source-detection significance (`--nsigma`), not FWHM,
aperture size, or matching distance. It is forwarded to every image; changing
it changes source detection and should be deliberate.

Optional `extra_args` supports explicit long tuning options, for example:

```python
extra_args=["--downsample", "2"]
# Or, with justified scale bounds for the actual camera/binning:
# extra_args=["--scale-units", "arcsecperpix", "--scale-low", "0.3", "--scale-high", "0.4"]
```

The scale numbers above only illustrate syntax, not recommended SARA values.
Backend/index paths can be supplied with `--config`, `--index-dir`, or
`--index-file`. Paths are resolved before entering the temporary directory.
Output/input redirection, abbreviated options, and duplicates of named controls
are rejected. Unsupported tuning options raise an error rather than being
silently ignored. The full supported list is in `_extra_options`.

## FITS output and failure safety

Every solver call receives a private copy of its input in a fresh temporary
directory. Only that directory is cleaned; neighboring `.wcs`, `.axy`, `.corr`,
or other user files are never swept or removed. The solver does not publish its
own rewritten science FITS. Instead, Photozpy validates its `.solved` success
marker and `.wcs` solution, then atomically replaces the primary WCS header.

Validation requires exit status zero, a successful marker, and a finite,
nonsingular 2D celestial WCS. SIP distortion is retained. Old coordinate matrix,
alternate-WCS, and SIP cards are removed before inserting the new solution.
Observation times and other scientific metadata are preserved. This structural
validation cannot independently prove the sky identification is correct; inspect
catalog overlays/residuals on actual data before scientific photometry.

Science pixels are not shifted, resampled, multiplied, or masked. All original
HDUs, including MASK and CMASK, are retained unchanged. The image-based solver
does not explicitly apply those masks during source extraction. A failed solve
does not publish its destination. An optional `timeout` is a wall-clock limit
in seconds; None adds no Python-side limit.

For the retained single-file API, `clean=True` defaults to in-place WCS update;
`clean=False` defaults to `<stem>solved.fits`. Explicit `output_path` overrides
that naming. Temporary solver products are disposable in both cases. Unlike
the old implementation, in-place updates require explicit `overwrite=True`.

Primary references: [Astrometry.net documentation](https://astrometry.net/doc/readme.html)
and [solver option definitions](https://github.com/dstndstn/astrometry.net/blob/main/solver/augment-xylist.c).
