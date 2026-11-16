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

For automatic reference selection, replace the dictionary above with:

```python
reference_images = "auto"
```

This selects the first selected Light image of each FIELD in collection order,
separately within each child directory. Bias/dark/flat and Master frames are
excluded. The first frame may be in any observed filter; all filters in that
field use its grid. In auto mode reference_filter is ignored. This is a
convenience selection, not a quality ranking. The selected filenames appear in
the registration log and `registration.reference_paths` after success.
Repeating the call with the same selection/order reuses compatible registered
images. If the first image changes, the existing reference-identity checks still
reject incompatible registration before any writes; restore pre-registration
working copies to select a different grid.

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
parameter search after a failure unless `detection_sigmas` is explicitly supplied.
Astroalign looks for common geometric patterns;
it does not require every star found in a reference to appear in the other image.
By default it uses up to 50 flux-ranked control points per image. Selecting a
reference band is therefore not a rule to always choose the weakest band.

To try several detection thresholds for each image, supply an explicit list:

```python
analysis_images = registration.register_images(
    reference_images=reference_images,
    detection_sigmas=[5, 6, 7, 8, 9, 10, 11, 12],
    min_area=8,
)
```

`detection_sigmas` overrides the single `detection_sigma` value. Each target starts
at the beginning of the list and stops at its first successful alignment. Values
are used in their supplied order; they are not sorted or expanded as endpoints.
Positive finite fractional values are also accepted. The list must be nonempty.
Only star-matching failures (no acceptable transform or too few detected stars)
trigger another attempt. Other errors fail immediately. Attempts always use the
original target pixels, so failed attempts do not accumulate resampling. If all
candidates fail, the error lists the attempted values and chains the last cause;
the failed target remains unwritten. INFO messages report retries and success.
The successful threshold is recorded as `REGSIG` in newly registered FITS headers.
Existing matching registrations are reused even if the candidate list changes;
this option does not select an optimal fit among successful thresholds.

### Optional whole-field smoothing fallback

For split or poorly focused star profiles, enable a second detection recipe:

```python
analysis_images = registration.register_images(
    reference_filter="SDSS_g",
    detection_sigmas=list(range(3, 16)),
    min_area=8,
    fallback_smoothing_sigma=2.0,
    discard_input_wcs=True,
)
```

The default `fallback_smoothing_sigma=None` preserves the existing behavior.
A positive value is the Gaussian standard deviation in **pixels**, independent
of the detection thresholds. For each field, ordinary matching first tries all
supplied thresholds for each target. If one target exhausts them, all staged
results for that field are discarded and **every target** is matched again
using smoothed detection copies and the same threshold list. Earlier ordinary
successes are not reused in the smoothed pass. The reference grid stays fixed.

Only temporary images used to find stars are smoothed. Their MASK and CMASK
pixels are excluded using normalized convolution to avoid spreading masked
artifacts. The transform is applied to the **original unsmoothed science data**
once; the existing validity-mask/CMASK propagation remains unchanged.

A complete field is staged on disk before any of its inputs or final outputs
are replaced. If the smoothed pass also fails, that field remains untouched.
Other fields already completed remain completed. Publication is atomic per
file; an I/O error while publishing can still leave a partial field, which must
be restored from pre-registration copies. Temporary storage must accommodate
one or two field outputs during the attempt; it is cleaned on return/error.

Each output, including the unresampled reference, records `REGSMTH`: `0.0` for
ordinary detection or the Gaussian sigma for smoothed detection. `REGSIG`
continues to record the successful detection threshold on each target.
Repeated calls reuse a **complete, compatible field**, including a field with
only its reference image, without resampling or rewriting it in place. A
reference-only selection preserves its existing mode even if it was previously
part of a larger smoothed field.

Once **any** selected member carries `REGSMTH`, all selected members must be a
complete, consistently marked field. This check also runs when the caller omits
`fallback_smoothing_sigma`; omitting it cannot silently turn a recorded smoothed
field into a mixture. Native exposures added to a recorded field, partial
publication, missing/invalid modes, mixed modes, or changing an existing positive
smoothing width require restoring pre-registration copies. Incremental addition
to an already recorded field is not implemented. `REGSMTH=0` remains compatible
with any positive fallback setting because no smoothing was actually used.

Pure legacy selections without `REGSMTH` retain default reference-based reuse;
opting those already resampled targets into the new fallback requires native
copies. Independently of this option, an absent/empty `ALIGN` with residual
`REGSIG`, `REGSMTH`, `REGREF`, or `REFHASH` is rejected before writing. Never clear
registration markers to bypass this check. If **all** provenance markers are
removed, headers alone cannot establish whether the target pixels were already
resampled. Likewise, manually replacing target pixels while preserving old
markers is outside the supported reuse contract; target pixels have no separate
integrity digest. A before-cosmic-ray backup needs cosmic cleaning again before
registration.

After success, the processor's image_collection and the returned collection
select the output files, so repeating that processor can reuse those completed
outputs. Constructing a new processor from the original native collection with
the same `save_location` instead requires overwrite=True for existing outputs
and recomputes from native pixels; it does not implicitly load destination
images. Select the completed output collection when you intend to reuse it. Matching failures report the field, actual phase,
threshold sequence and smoothing width. Publication errors report the failed
path, number of outputs already published and recovery instructions.

A successful fit is not a quality guarantee. Inspect registered stars and the
master images; this option is a failure fallback, not an automatic quality
ranking or a correction for the original focus/PSF.

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
present or detection smoothing is active (otherwise register). CMASK excludes detections only during star finding;
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
stacks and retains their output mask. Combine does not independently validate
REGSMTH recipe consistency; the registration checks above enforce that contract.
Do not rely on Combine to detect manually assembled mixed-recipe inputs.
Downstream photometry must honor MASK.

Legacy cleaned files storing cosmic detections in MASK without CMASK are rejected
before registration or multi-image Light combination. Restore pre-registration
working copies and separate the detections first; mixed masks in already
resampled files cannot reliably be split by renaming an extension. Do not
resample an already registered image again. Existing observing files are never
automatically migrated.

## Telescope-provided WCS before registration

Some cameras attach an approximate celestial WCS to raw exposures. For a
workflow which solves WCS again after Master stacking, explicitly use:

```python
analysis_images = registration.register_images(
    reference_images="auto",
    detection_sigmas=[5, 6, 7, 8, 9, 10, 11, 12],
    min_area=8,
    discard_input_wcs=True,
)
```

Registration still matches image stars, independently of the input WCS. This
option removes old spatial WCS, alternate solutions and SIP coefficients from
newly resampled output headers and records the removal in HISTORY. It does not
modify the unresampled reference's WCS or WCS on already registered targets
which are reused. It does not solve WCS or propagate the old solution through
the transform. Solve WCS on the resulting Masters before sky-coordinate use.
By default this option is False, preserving the existing pre-write rejection of
WCS-bearing targets. When saving to a separate directory, input files remain
unchanged; in-place registration replaces the analysis working copies.
