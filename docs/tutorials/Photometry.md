# Optical aperture photometry: catalog standards to calibrated magnitudes

Run this after Master Light plate solving and FWHM measurement. The workflow
identifies catalog standards, reviews per-image aperture/background regions in
DS9, measures standard and science sources together, then estimates zero points
and calibrates the targets. The catalog query alone does not measure flux.

For a Chinese notebook-authoring checklist and the validated SARA cell sequence,
see [Photometry_Notebook_Memo.md](Photometry_Notebook_Memo.md).

## Select the catalog query image

The analysis notebook already defines `analysis_dir`, `analysis_images`, and
`DLS_field`. Use the first child collection for this night's images and select
one reference filter from the field's Telescope. Selection uses FIELD, not OBJECT.

```python
from pathlib import Path

from photozpy.standard_stars import StandardStarCatalog

# Select the standard observation field and one reference filter.
standard_field = DLS_field
standard_reference_filter = standard_field.telescope.filters()[-1]

# Get this field/filter's Master Light from the existing analysis collection.
standard_star_path = Path(
    analysis_images[0].files_filtered(
        include_path=True,
        imtype="Master Light",
        field=standard_field.field_name,
        filter=standard_reference_filter,
    )[0]
)
```

This example expects one Master Light for the selected field/filter. Set
standard_reference_filter to standard_field.telescope.resolve_filter_alias("SDSS g")
to choose a particular reference band instead. Region generation still processes
all available filters independently.

## Query the catalog and attach reference magnitudes

Choose the catalog path explicitly. The following uses the Smith Southern Sky
catalog from the previous analysis. The AB label follows the metadata of this
compiled catalog; no AB/Vega conversion or passband transformation is performed.

```python
catalog_dir = Path("/Users/yong/Desktop/Research/Photozpy/Compiled_Standard_Stars")
catalog_path = catalog_dir / "Smith_Southern_Sky.h5"
standard_star = StandardStarCatalog(catalog_path=catalog_path)
standard_star.print_metadata()

standard_sources = standard_star.find_standard_sources(
    fits_path=standard_star_path,
    observation_field=standard_field,
    magnitude_system="AB",
    edge_margin_pix=100.0,
    mag_limit=23.0,
)

print(f"Selected standard-star candidates: {len(standard_sources)}")

# Display all candidates and their catalog measurements in one notebook table.
standard_star.matched_stars
```

`matched_stars` is a pandas DataFrame with one row per selected star. It retains
all catalog columns (including reference magnitudes and uncertainties), initial
pixel coordinates, and the original row index/order. IDs are decoded for display.
Its rows correspond to the Sources returned by the same query, after every cut.
It is a defensive copy: editing the displayed table does not change Sources or
the stored snapshot. Before a query it raises RuntimeError; an empty completed
query returns an empty table with columns intact. A new completed query replaces
the snapshot. A failed query leaves the previous snapshot available. The old
`final_stars` attribute retains its pre-magnitude-limit behavior for compatibility;
use `matched_stars` for the fully filtered overview.

The 100-pixel edge cut and magnitude limit of 23 retain the old notebook's
parameter values. They are not universal quality thresholds. Filter and error
columns are now discovered through the field Telescope's aliases; the primed
SDSS aliases include the compiled catalogs' SDSS_up/gp/rp/ip/zp spellings.
Missing catalog bands remain undefined. PanSTARRS bands are not substituted
for SDSS bands. Ambiguous columns require an explicit override rather than an
arbitrary choice. For custom catalogs, register their filter aliases once in
the FilterSet, or supply magnitude_columns/error_columns as special overrides.

The magnitude cut requires all discovered magnitude columns to satisfy the
limit. SARA defines u/g/r/i/z, so automatic discovery includes u when the catalog
has it, even if the selected image is g. This differs from the old explicit
g/r/i/z-only cut; on the tested DLS image both selections return the same six
stars. Candidate selection uses WCS
footprint/edge geometry, not image-source detection or a measured centroid.
Coverage, saturation, neighbors, and suitability for photometry still need
inspection. A selected candidate is not automatically an accepted calibrator.

Each star has its own SourceInfo and reference FilterMagnitude. Catalog IDs are
decoded from bytes. Instrumental and calibrated slots remain None. The science
`sources` collection is left untouched, and no FITS or region files are written.
Keep `standard_sources` in this notebook; do not recreate it from source_infos
after attaching results, because that creates new undefined measurement slots.

## Inspect candidate positions in interactive DS9

```python
standard_viewer = standard_star.preview_standard_stars(
    standard_star_path,
    standard_sources,
)
```

This opens/reuses the dedicated photozpy-standard-preview DS9 target, displays
the image, and overlays labeled blue crosses at the initial catalog positions.
It replaces that target's current regions. No matplotlib figure, temporary
region file, or Tk control window is created. DS9 remains interactive and open.
These initial reference crosses are locked against editing. This is a visual
candidate check, not the later editable aperture/annulus acceptance workflow;
it does not persist manual decisions. The old plot method remains only for
backward compatibility and is not part of this tutorial's workflow.

## Reuse an already inspected catalog selection

The original DataFrame-returning method remains available. If the notebook has
already produced and inspected `standard_star_table`, use:

```python
standard_sources = standard_star.to_sources(
    standard_star_table,
    observation_field=standard_field,
    magnitude_system="AB",
)
```

Coordinates use catalog RA/Dec units from metadata and default to ICRS, matching
the existing footprint query. `to_sources` also accepts an explicit
`coordinate_frame` for an independently prepared table; it does not infer that
frame from image WCS. The wrapper `find_standard_sources` retains the existing
ICRS query convention.

No catalog-specific missing-value sentinel is translated automatically. NaN
measurements remain NaN; numeric sentinels remain numeric. Error mappings may
omit filters, leaving those errors undefined. Significance is undefined because
this step has not measured an image detection. Mixed-system columns must not be
labeled with a single system unless they have been explicitly transformed.

## Generate regions for all standard-field filters

```python
from photozpy.photometry.regions import PhotometryRegions

# Find this field's Master Lights and centroid the standards on every filter.
# Each image supplies its own FWHM; aperture and annulus files remain independent.
standard_regions = PhotometryRegions.generate_regions(
    image_collection=analysis_images,
    sources=standard_sources,
    root=analysis_dir,
    aperture_scaling=3.0,
    inner_annulus_scaling=5.0,
    outer_annulus_scaling=8.0,
    box_size=19,
)

# Report the number of image/group entries queued for review.
print(f"Standard-field review entries: {len(standard_regions)}")
```

Standards share one pair per field/filter. Different filters are stored in the
same field directory, not separate image directories:

```text
photometry/regions/
    DLS_field_standards/
        SDSS_g'_aperture.reg
        SDSS_g'_background.reg
        SDSS_r'_aperture.reg
        SDSS_r'_background.reg
        ...
```

Each file has comments recording FIELD_NAME, FILTER_NAME, SOURCE_NAMES and
the root-relative IMAGE_PATH. Each circle and
annulus carries `tag={id:<catalog id_name>}` plus its role tag. Order is not used
to associate stars. Names containing DS9 tag delimiters are rejected.

## Generate science-target regions in one operation

```python
# Use the science Sources already defined in the analysis notebook.
# Select their fields' Master Lights automatically, across all available filters.
science_regions = PhotometryRegions.generate_regions(
    image_collection=analysis_images,
    sources=sources,
    root=analysis_dir,
    aperture_scaling=3.0,
    inner_annulus_scaling=5.0,
    outer_annulus_scaling=8.0,
    box_size=25,
)

# Inspect every source/filter record together without iterating over targets.
import pandas as pd

pd.DataFrame(science_regions.records)
```

Science targets are split automatically into per-source directories. Filenames
include FIELD and FILTER, so the same source name observed in different fields
cannot overwrite another field's files. For example:

```text
photometry/regions/
    4FGL_J0746.6-4754/
        fgl_field_SDSS_g'_aperture.reg
        fgl_field_SDSS_g'_background.reg
        fgl_field_SDSS_r'_aperture.reg
        fgl_field_SDSS_r'_background.reg
        ...
```

Both single ImageFileCollection and ImageCollectionGroup inputs work. Explicit
file selections are preserved: directories are not rescanned. Review order is
field order in Sources, source group, then Telescope filter order and image. Each queue
entry edits one image/group pair; multiple science targets sharing a field have
individual file pairs and review entries. Missing fields, duplicate output
names and existing files are rejected before any batch outputs are created.
Computation/I/O failures later in generation can leave earlier entries written.

The single-image generate(image_path, sources, root=...) interface remains
available for one standard group or one science target. The same naming rules
apply to both interfaces. box_size defaults to 19 for standards and 25 for
science in batch mode; an explicit value overrides that choice.

The old CCD_Regions/get_centroids implementation stays in swiftz for existing
users. The new module preserves its numerical main line without importing the
legacy Sources/mImageFileCollection or plotting/network helper module:

1. Project catalog coordinates with image WCS (zero-based pixels).
2. Subtract Background2D using MedianBackground, 20-pixel boxes, 5-box filtering,
   and SigmaClip(sigma=3), as in the old estimate_background helper.
3. Run two centroid_com iterations and a final centroid_quadratic in the local
   window; try centroid_com if the final quadratic fails.
4. Reject warnings, non-finite/invalid centers, or displacement beyond
   max_centroid_shift. Default maximum displacement is box_size // 2 pixels.
   Failure returns the initial position and records the reason. Background
   estimation failure also produces explicit initial-position fallback records.

MASK and non-finite pixels are excluded during background estimation and
centroiding; repaired CMASK pixels are retained. This is the explicit mask
handling added to the legacy algorithm. A finite centroid remains an estimate,
not proof that a faint source was detected. Cross-check it in DS9.

Scaling factors must satisfy 0 < aperture < inner annulus < outer annulus.
These controls set generation defaults; manual review can change the radii.
Partial aperture/annulus coverage is recorded but not automatically rejected or
flux-corrected. Valid-area fractions measure geometry, not the number of stack
exposures contributing to a pixel. No uncertainty or photometry is performed.

```python
import pandas as pd

# Inspect centroid methods, fallback reasons and image-coverage diagnostics.
pd.DataFrame(standard_regions.records)
```

## Edit apertures and annuli in interactive DS9

```python
# Load this image with both editable region files and fixed initial crosses.
# Save & Next validates and saves edits; on the last image it completes review.
standard_regions.review(sync_deletions=True)

# Inspect final coordinates, radii and per-image acceptance after closing controls.
pd.DataFrame(standard_regions.records)
```

Region review preserves the current zoom and image-coordinate view center when
Next or Previous switches filters within the same source/group and field. First
display and a switch to another target/field use zoom-to-fit. Each target's
filters stay together in the batch queue. No extra notebook parameter is needed;
ordinary image preview retains its existing zoom-to-fit default. This display
policy does not change region coordinates or saved scientific data.

The dedicated region-review window uses existing preview navigation and the
macOS isolated Tk worker. Green circles and red annuli are editable. Blue
initial-position crosses are fixed and excluded from measured files. Only the
green aperture displays the source name; red annuli retain hidden source ID tags
for pairing without a second overlapping text label. This also applies when
displaying previously generated background files, without rewriting them.

Centers start concentric but may be moved independently, e.g. to avoid a bright
neighbor in the sky region. Save & Next retains both positions and displays a
red, non-blocking warning for offsets greater than 0.01 pixel (source ID and
offset in pixels). Aperture/annulus overlap also produces a warning. Check that
the edited sky region still represents local background; neither warning
automatically moves regions or stops navigation. Records include background_x,
background_y and center_offset; x/y remain the source-aperture center. Initial
generation coordinates are retained for comparison.

Deleting both shapes rejects that star; deleting just one blocks saving.
With `sync_deletions=True`, Save & Next also excludes that standard star's paired
regions from all unreviewed sessions in the current queue with the same analysis
root, FIELD and source ID. Science targets and already reviewed image/source
decisions remain unchanged. The controls report how many stars and images were
updated, and how many reviewed images still retain those stars. Moves, radius
changes and background adjustments are never copied. Synchronized records remain
unreviewed and retain their own geometry, with `accepted=False` and a
`deletion_sync` record identifying the originating image, filter and manifest.
Both region files and the manifest are saved immediately, so exclusions survive
Quit, interruption and `load_regions()`. Quit and Previous discard unsaved edits
without propagating them. Previously saved exclusions are reapplied idempotently
on Save & Next, including when revisiting a review after a partial write failure;
restoring a star locally does not automatically restore it in other filters.
The default `sync_deletions=False` keeps deletions local to the current image.
Only sessions included in this review queue participate.

Editing a radius is allowed if the radius ordering remains
valid. Preserve the source ID and role tags when editing. Unsupported shapes or
unrecognized IDs block saving, rather than silently corrupting associations.

The review controls include **Reset display**. It restores contrast=1,
bias=0.5, global zscale limits and an asinh stretch. The first displayed image
also initializes these settings so an existing DS9 window cannot carry extreme
contrast/bias into a new review. Later manual display choices are preserved
until explicitly reset. Reset does not reload FITS, change zoom/pan, remove
regions, save/discard edits, or modify science pixels. Use it if the view turns
mostly white/black or strongly emphasizes background noise after colorbar
adjustments. Right-drag controls colorbar contrast/bias, separately from numeric
ADU scale limits.

Save & Next reads DS9 image-coordinate regions, separates circles/annuli,
validates every pair, and publishes both files and the JSON record. It advances
only on success. On the last image, this same button saves before completing.
Previous reloads saved regions and discards current unconfirmed edits; it does
not undo earlier saves. Quit does not save unconfirmed edits. Saved decisions
survive interruption. Files are individually replaced atomically after staging;
this is not a filesystem-wide transaction, and an I/O failure stops navigation.
Keep the generated JSON record together with its two region files.

```python
# Review all science targets and filters with one control window.
science_regions.review()

# Display saved positions, background offsets and acceptance for every entry.
pd.DataFrame(science_regions.records)
```

Each image keeps its own FWHM, centroid and acceptance. The existing review
object can be reopened with review(), without repeating generation. Initial
coordinates and automatic results remain in internal records after manual
editing; SourceInfo coordinates remain unchanged. Internal per-filter JSON
files are managed by the program for record keeping and GUI process transport;
the normal end-to-end notebook does not locate or load them manually.

## Measure reviewed apertures: counts, errors and instrumental magnitudes

If kernel restart removed only the in-memory region objects, do not regenerate
or repeat manual review. Once analysis_images, standard_sources, sources and
analysis_dir are defined again, read the saved associations automatically:

```python
from photozpy.photometry.regions import PhotometryRegions

# Read the saved standard-star regions and review decisions for all filters.
standard_regions = PhotometryRegions.load_regions(
    image_collection=analysis_images,
    sources=standard_sources,
    root=analysis_dir,
)

# Read the saved science-target regions, retaining all manual edits.
science_regions = PhotometryRegions.load_regions(
    image_collection=analysis_images,
    sources=sources,
    root=analysis_dir,
)
```

load_regions shares generation's selection/naming rules but never calls
centroiding, changes radii, writes files or opens DS9. Missing files produce an
explicit error instead of silently generating replacements. No manual JSON path
is needed. Existing source identities must match the saved groups. This restores
region objects, not source metadata/catalog magnitudes or earlier processing.

After completing region review, use one Photometry call for all source roles
and fields. Separate standard/science region files do not require separate
measurements. This step does not estimate zero points or calibrated magnitudes.

```python
from photozpy.photometry.photometry import Photometry

# Use the existing explicit image selection; do not rescan the directory.
photometry = Photometry(image_collection=analysis_images)

# Measure both existing source collections in one operation.
# Shared result records preserve reference values and update the original objects.
photometry_table = photometry.run_photometry(
    sources=(standard_sources, sources),
    regions=(standard_regions, science_regions),
    bkg_clip_sigma=3.0,
    src_detection_sigma=3.0,
    verbose=False,
)

# Display every source/filter measurement in one unit-aware table.
photometry_table
```

```python
# Optionally inspect science-only rows without running photometry again.
photometry_table[photometry_table["source_role"] == "science"]
```

One mixed Sources collection is also accepted. Sources.combine(standard_sources,
sources) creates a view sharing all original measurement slots; unlike rebuilding
Sources from source_infos, it does not reset reference results. A sequence of
region sessions/collections is flattened internally. The same frame may contain
both roles: sessions are grouped by image and data/header/mask are read once per
frame. Each source still uses its own final aperture and background geometry.
There is no standard-frame/science-frame classification in this measurement step.
Independent records claiming the same field/source identity are rejected before
measurement; repeated views of the same shared source record are deduplicated.

```python
# Inspect one source's stored instrumental magnitudes if needed.
# The returned batch table above remains the normal global overview.
sources[0].get_magnitudes(kind="instrumental").qtable
```

### Preserved numerical definitions

The original estimate_source_region_error and counts2mag executable statements
are unchanged. The original run_photometry measurement sequence is retained:

1. Sigma-clipped annulus mean (SigmaClip maxiters=10), with the existing
   mean/median discrepancy warning. No automatic median substitution.
2. Source aperture geometric area and clipped annulus sum_aper_area, keeping the
   old n_ann definition. Source/background regions are read by ID with their
   final independent positions and per-source radii; no FWHM scaling is repeated.
3. Exact aperture sum, subtract mean sky times the source area.
4. Original source+sky/readout/background propagation model, divided by NCOMBINE
   for average stacking, and its separate background-only upper-limit error.
5. Original counts2mag conversion, S/N threshold and background-only flux limit.

Inputs remain ADU, detector gain in electron/ADU, detector read noise in electron,
and dimensionless NCOMBINE. Instrumental magnitudes are -2.5 log10(net ADU counts),
with no exposure-time normalization. No CCDData uncertainty array or new
covariance/calibration-frame noise model is introduced. The inherited model's
assumptions remain unchanged.

### Masks, rejected regions and output tables

Read MASK and non-finite pixels as invalid; do not exclude repaired CMASK pixels.
An incomplete source aperture is reported as invalid_aperture instead of silently
estimating its missing source flux or applying an area correction. Background
invalid pixels are excluded from the original mean estimator and effective area.
An unusable background/noise result has an explicit status and reason. A star
deleted in review is reported as rejected for that image. No automatic replacement
measurement or magnitude is invented for these rows.

Only usable measurements update instrumental slots. Missing/rejected filters in
a measured source remain NaN; a source with no usable measurements has None in
its instrumental slot. Reference/calibrated slots are never changed. Repeating a
batch replaces the instrumental result for sources attempted in that batch, so
pass all desired filters together. Duplicate source/filter measurements fail
instead of overwriting one another. Review is required by default; explicitly
set require_review=False only to measure generated regions without manual review.

Old columns src+bkg, bkg, src, mag_inst, mag_inst_error and src_significance are
retained. Additional columns include source_id, source_role, field/filter, root-relative paths,
source/background centers, radii, n_ap/n_ann, mean_bkg, error_counts,
error_counts_ul, is_upper_limit, status and reason. Counts and count errors have
ADU units; magnitude columns have mag units. Tables keep the original NaN
magnitude error for non-detections. The attached FilterMagnitude alone uses its
established -99 magnitude-error sentinel to identify a valid upper limit.

Photometry.results and phot_table contain the last completed unified batch;
photometry_table is an independent copy for notebook inspection.
Measurements do not write FITS, region or output-table files. Exporting tables
and the zero-point/calibration operations below remain separate from photometry.

## Standard-star zero points and target calibration

Run these steps after the unified photometry call above, using the same
`standard_sources` and `sources` objects. No image paths, region files, filter
mapping or manual transfer of magnitudes is needed at this stage.

```python
# Calculate standard-star zero points separately for each telescope filter.
# Preserve the original single-pass 3-sigma z-score selection and unweighted mean.
zero_point_table = standard_sources.calculate_zero_points(sigma=3.0)

# Inspect the mean, retained scatter and sample counts for every filter.
# Numerical zero_point/scatter values are in magnitudes, not flux units.
zero_point_table
```

```python
# Show every standard/filter offset, including excluded stars and their reasons.
standard_zero_point_table = standard_sources.zero_point_details

# Display the complete pandas DataFrame for a global overview.
standard_zero_point_table
```

```python
# Inspect one panel per measured filter before applying the zero points.
# Red lines are average +/- 0.1 mag visual guides, not uncertainty bounds.
zero_point_figure = standard_sources.plot_zero_point_statistics()
```

```python
# Calibrate all science targets using the same zero points inspected above.
# Preserve instrumental errors/significance and the existing upper-limit flags.
calibrated_table = sources.calibrate_targets(standard_sources=standard_sources)

# Display a compact magnitude/error/significance summary for each source.
sources.show_magnitudes(kind="calibrated")
```

show_magnitudes prints each source name, ICRS RA/Dec, and field once, followed by its magnitude
system and a compact filter/magnitude/error/significance table. It reuses the
existing FilterMagnitude.to_table formatting: missing values are `--` and
upper-limit errors are `UL (-99)`. Undefined results are reported explicitly.
The method only displays existing results; it does not calculate or save them.
Use the same entry point with `kind="instrumental"` or `kind="reference"` to
inspect those slots. If the collection contains both roles, select
`sources.by_role("science")` first to display only science targets.

For per-source copy boxes matching a griz spreadsheet, use
`sources.show_magnitudes(format="google_sheets", filters=("SDSS_g", "SDSS_r", "SDSS_i", "SDSS_z"))`.
The default block contains obs mag, STD correction (calibrated minus instrumental),
corrected mag, error and detection significance. Use
`columns=("magnitude", "error")` to paste directly
into existing magnitude/error columns. See [Sources](Sources.md) for all output
options and missing-value conventions.

calibrated_table remains available as the complete pandas DataFrame for export
or detailed inspection; it does not need to be the notebook's normal display.

The equations follow the original MagnitudeInfo implementation:

```text
ZP_i = reference_magnitude_i - instrumental_magnitude_i
ZP_filter = mean(retained ZP_i)
calibrated_magnitude = instrumental_magnitude + ZP_filter
calibrated_error = instrumental_error
```

Eligibility is evaluated independently per filter: finite reference and
instrumental magnitudes, no instrumental upper-limit flag, and finite detection
significance >= 3 by default. Set `min_significance` explicitly on
calculate_zero_points if a different detection cut is intended. Science-role
sources are ignored when estimating zero points. Standards must share one
Telescope instance and one reference system (AB or Vega); systems are never
converted or mixed implicitly.

Clipping preserves the old `scipy.stats.zscore` calculation: one pass around the
initial mean with population standard deviation, then an unweighted mean of
retained offsets. It is not iterative clipping. With few standards, 3-sigma may
retain a visually separated point. Singleton and exactly constant samples are
retained explicitly, fixing the old undefined-z-score rejection of all values.
If all filters lack eligible/retained stars, calculation raises an error.

The summary DataFrame contains `filter`, `zero_point`, `scatter`,
`n_candidates`, `n_used`, `n_clipped`, and `magnitude_system`. Scatter is the
population standard deviation of retained offsets; it is not an uncertainty
of the mean. Missing filters have NaN zero point/scatter and zero sample counts.
The detail DataFrame includes source ID, field/filter, reference/instrumental
magnitudes, significance, individual zero point, `used` and `reason`. Both
tables are independent copies. `standard_sources.zero_points` retrieves another
copy of the current summary without recomputation.

Plotting preserves the old per-filter scatter/mean layout. Retained offsets are
blue circles, sigma-clipped offsets are gray crosses, other excluded finite
offsets are orange crosses, and the retained mean is green. Purple dashed lines
show the actual initial eligible mean ± sigma × population standard deviation
used by the single-pass z-score rejection. They are stored with the calculation,
not inferred from the retained mean/scatter. Notes show sigma-clipped counts and
unavailable-value reasons. Titles explicitly count retained and rejected offsets, including zero
rejections. Rejected counts include finite offsets excluded by the significance
cut as well as sigma clipping; `n_clipped` in the summary counts only clipping.
Axes include all displayed offsets, even outliers beyond the red guide lines.
Unavailable offsets (for example, a missing reference magnitude or an upper
limit) have no valid zero point to plot and are counted separately when present.
The function returns a Matplotlib Figure and writes no file unless `save_path`
is supplied. It handles one filter or any panel count; filters without finite
offsets are omitted. The +/- 0.1 mag red dotted guides are unrelated to the clipping
threshold, retained scatter or calibration uncertainty.

calibrate_targets reuses the calculated zero points, including custom clipping
thresholds. If no calculation exists, it computes one using defaults. If shared
standard-star measurements have changed (including edits to their QTables), it
requires an explicit calculate_zero_points call before plotting or calibration.
This keeps the inspected and applied zero points consistent.

Calibration operates on science-role sources only, whether standards occupy a
separate field or share a science frame. Target and standard measurements must
use the same Telescope instance. It writes only `calibrated` slots and always
starts from original instrumental values, so repeated calls do not add the
offset twice. Reference/instrumental slots remain available. Missing instrumental
filters remain NaN; sources without measurements have no calibrated result.
A measured target filter without a usable zero point raises before any target
calibrated slot is changed.

Output magnitude/error values in calibrated_table are numerical magnitudes.
Stored results are FilterMagnitude objects with unit-aware QTables. Detection
significance and instrumental errors are copied unchanged, including the -99 mag
upper-limit flag. A limiting magnitude receives the same additive zero point.
Catalog errors and zero-point uncertainty are intentionally not propagated;
adding them is a separate future feature. No exposure-time normalization,
extinction correction, color term or aperture correction is introduced here.


## Optional local-background source masking and checks

Each image keeps its own header FWHM and reviewed source/background geometry.
Enable the following only for the reviewed-region photometry interface:

```python
from photozpy.photometry import LocalBackgroundOptions

photometry_table = photometry.run_photometry(
    sources=(standard_sources, sources),
    regions=(standard_regions, science_regions),
    background_options=LocalBackgroundOptions(
        mask_growth_fwhm=1.5,
        min_area_fraction=0.5,
        min_pixels=50,
        stability_growth_fwhms=(1.0, 1.5, 2.0),
        stability_tolerance=1.0,
    ),
    verbose=False,
)
```

SEP detects source footprints once per image at `detection_sigma=3`, `min_area=5`.
An image-wide background map is used only for this detection. The photometric
sky estimate remains the sigma-clipped **mean inside each local annulus**.
Detected footprints are expanded by the requested multiple of that image's
FWHM and excluded only from annulus statistics. Source aperture sums, saved
regions, image data and calibration are not modified. CMASK is not used as a
validity or neighbor mask. If no sources are detected, only validity masks apply.

Effective annulus area is measured after validity masking, source masking and
sigma clipping. An estimate below either minimum is `invalid_background`, with
no magnitude assigned; `n_ann`, `mean_bkg` and `background_area_fraction` expose
available diagnostics. Growth choices are not used to maximize S/N. The primary
measurement always uses `mask_growth_fwhm`. For each candidate in
`stability_growth_fwhms`, the background is estimated again. The largest absolute
change from the primary mean, multiplied by source aperture area and divided by
the primary **existing model count error**, is `background_stability_sigma`.

`background_stability` is `stable`, `unstable`, `insufficient_area` (a trial has
inadequate sky area), or `not_checked`. `background_trials` contains the candidate
means, areas and validity as JSON; `background_warning` explains flagged results.
Unstable trials do not change primary counts, choose a different mask, or replace
the detection/UL status. Inspect these columns before using a result. The table
metadata stores all options. An empty trial tuple explicitly disables stability
checks; `mask_neighbors=False` enables area checks alone, without requiring FWHM.
Omitting `background_options` preserves the old computation.

These checks do **not** calibrate count errors for interpolation correlations,
residual background structure or unresolved sources. Stable mask growth does not
establish a reliable Gaussian detection significance. Use independent empirical
noise validation before treating model S/N as a calibrated detection confidence.


## Inspect the actual background masks

After a successful `run_photometry(..., background_options=...)`, pass one source
object to inspect **all available filters in one figure**, one filter per row:

```python
source = sources["PMN J0807-7629"]
fig = photometry.plot_background_diagnostics(source)
```

The object supplies source identity, field and telescope filter order. A
single-source `Sources` selection or `SourceInfo` is accepted. You can also pass
an entire `Sources` collection to draw one all-filter figure for each source:

```python
figures = photometry.plot_background_diagnostics(sources)
# For standards: photometry.plot_background_diagnostics(standard_sources)
plt.show()

# A multi-source call returns {(field_name, source_name): Figure}.
fig = figures[("PMN", "PMN J0807-7629")]
```

Source order follows the input collection. Same-name sources in different fields
remain separate. Sources with no snapshots are explicitly reported and skipped;
no snapshots for the whole collection is an error. A single-source call continues
to return one Figure. No source-name, field-name or filter-name
arguments are needed. Only actual diagnostic snapshots from the latest run are
shown: filters without images, rejected sources or invalid source apertures do
not get fabricated masks. A source with no snapshots raises an explicit error.
`photometry.get_background_diagnostics(source)` returns copied snapshots as an
ordered filter-to-diagnostic dictionary for per-filter export. The existing
`plot_background_diagnostic(source_name, filter_name, field_name=...)` remains
available for inspecting a single filter.

The figure contains the original local pixels with reviewed aperture/annulus
outlines, the pixel-center mask used for the sky **mean**, and the mask used for
**exact effective area**. Colors identify retained pixels, invalid pixels,
detected-source exclusions, and sigma-clipped pixels. The exact-area view may
include fractional boundary pixels. Photutils independently clips these two
footprints; do not interpret their pixel counts as the same quantity. The area
sequence shows geometric, on-image, validity-filtered, source-masked, and final
clipped area. Stage differences are disjoint contributions in that order.
`background_geometric_area`, `background_on_image_area`, `background_valid_area`
and `background_source_masked_area` are also recorded in the table; `n_ann`
remains the final exact effective area used by the existing error calculation.

These are snapshots captured from the **same ApertureStats instance** used for
the primary measurement. Plotting does not reread images or rerun detection or
clipping. Snapshot data are cutouts, not full-image copies; source/field/filter
and original zero-based pixel offsets are retained. Callers receive independent
copies. Failed/insufficient primary backgrounds remain inspectable when an
annulus overlaps the image; rejected or invalid source apertures have no snapshot.
A new successful run replaces the snapshots, including clearing them if the run
omits `background_options`. A failed run preserves the previous successful set.

Export a mask snapshot to inspect it after restarting the kernel:

```python
from photozpy.photometry import BackgroundDiagnostic

diagnostic = photometry.get_background_diagnostic(
    "PMN J0807-7629", "SDSS_r'", field_name="PMN",
)
fig = diagnostic.plot()
fig.savefig(analysis_dir / "PMN_r_background_mask.png", dpi=180)
diagnostic.save(analysis_dir / "PMN_r_background_mask.fits")

# No input image, Sources object or region session is required to reload it.
restored = BackgroundDiagnostic.load(analysis_dir / "PMN_r_background_mask.fits")
fig = restored.plot()
```

FITS extensions are `MEANLABEL`, `AREALABEL`, `AREAWEIGHT`, `DETECTED`,
`SOURCEMASK`, and `METADATA`, plus the original local image in PRIMARY. Pixel
labels are 0=outside annulus, 1=used, 2=invalid, 3=source mask, 4=sigma clipped.
Exact area is the sum of `AREAWEIGHT` where `AREALABEL == 1`. Metadata contains
geometry, options, clipping threshold, image association and area stages. Export
refuses to replace an existing file unless `overwrite=True` is explicit. Export
and plotting do not modify scientific image/region files or acceptance decisions.
