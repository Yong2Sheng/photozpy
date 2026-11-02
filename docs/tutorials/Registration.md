# Register Light images within each field

First complete bias/dark/flat correction and inspect cosmic-ray removal
(`Artifacts.md`). Registration retains the original astroalign.register
algorithm and groups Light frames by FIELD within each child collection. All
filters in a field share one reference grid. Sources and CollectionManager are
not required. Registration does not solve celestial WCS.

```python
from pathlib import Path

from photozpy import ImageCollectionGroup, Registration
from photozpy.telescope import SARA_CT

# Choose one specific working FITS for every field.
# These filenames illustrate the interface; inspect quality before choosing.
reference_images: dict[str, Path] = {
    "fgl_field": analysis_dir / "4FGL_J0746_SDSS_g_001.fits",
    "SUMSS_field": analysis_dir / "SUMSS_J075625_SDSS_g_001.fits",
    "DLS_field": analysis_dir / "DLS_1052-05_SDSS_g_001.fits",
}

# Align all exposures and filters in each field to its selected reference.
registration: Registration = Registration(
    image_collection=analysis_images,
    telescope=SARA_CT,
    overwrite=True,
)

analysis_images: ImageCollectionGroup = registration.register_images(
    reference_images=reference_images,
    detection_sigma=2,
    min_area=8,
)
```

The dictionary keys must match FIELD headers exactly. Each value selects one
specific Light FITS from the current working collection; it does not select a
whole field or a whole filter. Before-cosmic backups and raw images are not valid
references unless they are themselves in the selected working collection. The
reference's pixel grid stays unchanged; every other image in that field is
transformed to it. Different fields may choose different reference filters.
No quality score or automatic best-image recommendation is implemented.

The detection settings above come from the legacy Data_Analysis notebook, while
the module defaults remain detection_sigma=5 and min_area=5. These parameters
control astroalign star detection, not cosmic-ray detection. Select an observed
reference band explicitly; the first filter in a telescope definition need not
be present in this night's images.

The original `Register_images(filter=...)` entry point remains available.
Without an explicit reference, the first matching candidate in filename order
is chosen separately for each field/directory; a valid existing Reference is
preferred on repeat calls. To choose references yourself, provide
reference_images={"field_name": reference_path}. Each reference must belong to
that field's selected images. This deterministic choice does not rank quality.
When every field has an explicit reference, omit reference_filter. If the
mapping is incomplete, reference_filter is required for the remaining fields;
otherwise validation fails before files are written. Explicit choices override
the filter shortcut. Relative paths are relative to the current process working
directory, so construct them from analysis_dir or use absolute paths.

## Detection settings and progress

Astroalign detects stars separately in each image and reference using SEP. It
estimates a background map and a global background RMS for each image; the
detection threshold passed to SEP is detection_sigma * background.globalrms on
the background-subtracted image. Thus detection_sigma=2 uses each image's own
noise estimate, not the same ADU threshold for every image and not an integrated
2-sigma stellar S/N. Lower values can recover faint common stars but also admit
noise; higher values can leave too few useful stars.

min_area=8 requires at least eight connected pixels in a detection region. It
does not mean eight stars, an aperture radius, or an eight-pixel FWHM. SEP's
default detection filter also participates, so this is not simply eight raw
pixels each above the threshold. Increasing area rejects small features but can
lose compact/faint stars; decreasing it admits smaller candidates and artifacts.
These are source-detection settings, not positional matching tolerances.

The 2/8 example is a starting configuration from the old notebook. Both API
defaults remain 5. The module forwards explicit values unchanged, with no hidden
parameter search after a failure. Astroalign looks for common geometric patterns;
it does not require every star found in a reference to appear in the other image.
By default it uses up to 50 flux-ranked control points per image. Selecting a
reference band is therefore not a rule to always choose the weakest band.

Each FIELD has a fixed progress-bar title. The postfix gives the direction:
`Aligning image.fits to reference.fits`. Only basenames are shown. The reference
is retained, so it is excluded from the bar's image count. Existing compatible
results display `Reusing ...` instead of being resampled again. INFO logging
reports the selected reference filename and final registered/reused counts.

Parameter and algorithm references:
[Astroalign API](https://astroalign.quatrope.org/en/latest/api.html),
[Astroalign examples](https://astroalign.quatrope.org/en/latest/examples.html),
and [SEP extraction](https://sep.readthedocs.io/en/stable/api/sep.extract.html).

## Output behavior and validation

The reference's pixels stay unchanged. Other Light frames are resampled using
astroalign's original find_transform/apply_transform operations when CMASK is
present (otherwise register). CMASK excludes detections only during star finding;
the cleaned science values are resampled normally. The binary CMASK is resampled
separately with nearest-neighbor interpolation for inspection, never merged
into output MASK. The output MASK contains propagated validity flags, uncovered
footprint, and non-finite output pixels. No new flux normalization
or uncertainty model is introduced. The original interpolation is retained;
validate centroid and aperture-flux behavior on actual images before photometry.

ALIGN retains the original Reference/Registered values. REGREF records the
reference basename and REFHASH records a digest of its pixel data, validity mask, and optional CMASK.
Matching registered images can be reused without repeated resampling. A changed
reference, missing legacy reference identity, or another grid requires restoring
pre-registration working copies, rather than deleting markers.

Inputs support a primary 2D image and optional MASK and CMASK, with matching shapes and
units within each field. Uncertainty arrays and unrelated extensions are
rejected. Spatial WCS on targets to be resampled is also rejected: transforming
WCS requires a separate explicit model. An unresampled reference may retain its
WCS. For this pipeline, plate solving can follow Master Light generation.

The constructor writes nothing. In-place updates require overwrite=True;
save_location can instead select a separate output directory. Returned
collections retain calibration/master members. Group/reference checks and
collisions are validated before writes. Each output is atomic, but a later
alignment failure does not roll back earlier completed files. The failed target
is never marked as successfully registered.

Only combine Light frames after registration. Multi-image FIELD/FILTER stacks
must share REGREF/REFHASH; Combine rejects unregistered or mixed-reference
stacks and retains their output mask. Downstream photometry must honor MASK.

Legacy cleaned files storing cosmic detections in MASK without CMASK are rejected
before registration or multi-image Light combination. Restore pre-registration
working copies and separate the detections first; mixed masks in already
resampled files cannot reliably be split by renaming an extension. Do not
resample an already registered image again. Existing observing files are never
automatically migrated.
