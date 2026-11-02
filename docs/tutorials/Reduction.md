# Master frames and bias/dark/flat reduction

`Combine` and `Reduction` retain the original public processing method names and
ccdproc numerical recipes. They now accept the refactored ImageCollectionGroup
or a single ImageFileCollection, and use the collection refresh/filter helpers.
No Sources, legacy Telescope, or CollectionManager is required.

Combine uses Python logging at INFO level to report the input count, wrapped
input filenames, the output master filename, and completion. Full paths are
omitted from this summary. Configure logging in the notebook to display INFO
messages; the library does not install handlers or change global logging levels.
Each group has a FILTER/FIELD title when applicable, a separator, and a blank
line after completion to distinguish it from the next group.
It does not show a master-count progress bar. Set quite=True to suppress messages.
Reduction retains per-image progress bars for the correction stages.
Each bar keeps a fixed stage description and updates its postfix, after the
progress statistics, with the current master and target filenames, for example
`Applying Master_Bias.fits to Dark_001.fits`. It omits full paths and does not
print a separate log line for every image.

```python
import logging

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logging.getLogger("photozpy").setLevel(logging.INFO)
```

## Working images and output

Start with the reviewed, header-corrected copies in `analysis/`, represented by
`analysis_images`. These images must have IMTYPE and BUNIT. Flat and Light also
need FILTER; Light stacking requires FIELD. Nothing below operates on raw data
unless the input collection itself points to raw data.

Keep each returned collection. It explicitly contains newly generated masters
and retains the other working images. Do not assume that refreshing an older
snapshot will discover a master created after that snapshot was selected.

Constructing either object performs no writes. Combine refuses to overwrite
existing master files by default. Reduction also requires an explicit
`overwrite=True` for in-place work. File overwrite permission does not permit a
second scientific subtraction: existing correction markers are checked.

## Notebook sequence

Run and inspect each step separately. This sequence uses your single analysis
working directory, rather than creating directories for each calibration step.

```python
from photozpy import Combine, Reduction
from photozpy.telescope import SARA_CT
```

### 1. Build the master bias

```python
# Sigma-clipped average of Bias frames; write Master_Bias.fits.
analysis_images = Combine(
    analysis_images, telescope=SARA_CT,
).combine_bias_or_dark(image_type="Bias")
```

### 2. Subtract bias

```python
# Subtract Master Bias from Dark, Flat, and Light; mark BIASCORR='Yes'.
analysis_images = Reduction(
    analysis_images, telescope=SARA_CT, overwrite=True,
).apply_bias_correction()
```

### 3. Build the master dark

```python
# Combine bias-corrected Dark frames with one common EXPTIME.
analysis_images = Combine(
    analysis_images, telescope=SARA_CT,
).combine_bias_or_dark(image_type="Dark")
```

### 4. Subtract scaled dark

```python
# Scale Master Dark by target/master EXPTIME; correct Flat and Light.
# Mark DARKCORR='Yes'. Both target and Master Dark must be bias-corrected.
analysis_images = Reduction(
    analysis_images, telescope=SARA_CT, overwrite=True,
).apply_dark_correction()
```

### 5. Build master flats

```python
# Normalize each Flat by its median, reject outliers, and average per FILTER.
analysis_images = Combine(
    analysis_images, telescope=SARA_CT,
).combine_light_or_flat(image_type="Flat")
```

### 6. Apply flat correction

```python
# Divide each Light by its mean-normalized Master Flat; mark FLATCORR='Yes'.
analysis_images = Reduction(
    analysis_images, telescope=SARA_CT, overwrite=True,
).apply_flat_correction()
```

Optional telescope arguments resolve aliases consistently; omit them if all
FILTER cards already match exactly. For calibration images, grouping uses IMTYPE
and FILTER, not OBJECT. Registered Light combinations use FIELD/FILTER and the
original unscaled average; registration must be done separately first.

## Scientific behavior and the flat-rejection update

- Multi-image stacking uses average, sigma clipping, np.ma.median, mad_std,
  default lower/upper thresholds of 5, and float32 calculation.
- Flat stacking retains inverse-median input scaling by default. Single-image
  groups retain the old unscaled copy behavior. Custom flat scaling callables
  supply positive finite scalars; None explicitly disables multi-image scaling.
- For multiple flats, scaling now happens **before sigma clipping**. In the
  previous ccdproc 2.5.1 call using scale=, clipping operated on unscaled data
  and scaling was applied during combination. With varying twilight levels,
  this could hide star pixels among the changing backgrounds. The new order
  is normalization, clipping with the existing median/MAD thresholds, then
  average. This is an intentional scientific behavior update for Flat only.
  Bias, Dark, Light stacking, clipping thresholds, and flat application are
  unchanged. It does not guarantee removal of stars at persistent positions
  or wings that remain below the clipping threshold.
- Bias uses subtract_bias. Dark uses subtract_dark with EXPTIME in seconds and
  scale=True. Flat uses flat_correct without a new normalization or floor.
- Reduction outputs are float32 and deliberately remove mask and uncertainty,
  as the old code did. This is not a new propagated-uncertainty implementation.
- Multi-image Master Flat outputs still drop mask/uncertainty. Master Light
  drops uncertainty but retains its mask so registered footprints stay excluded
  and pixels without valid coverage remain marked. Master Bias/Dark keep
  ccdproc's generated mask/uncertainty.
- Outputs follow CCDData serialization; unrelated input FITS extensions are not
  copied. These calibration operations currently require primary image ext=0.
  Header-only preservation of all HDUs is a different operation.

After flat correction, remove cosmic rays with Artifacts (see `Artifacts.md`),
inspect representative masks, and then run Registration (see `Registration.md`)
before multi-Light stacking. Each
FIELD/FILTER stack must share one recorded REGREF/REFHASH grid; unregistered or
mixed-reference stacks are rejected. A single Light can still be copied without
resampling. Downstream photometry must honor the MASK extension in Master Light.

Normalized flat inputs are held in memory one filter group at a time; raw and
individual calibrated flat files are not rewritten by normalization. Existing
master flats must be explicitly rebuilt to use the new rejection order. If
Light frames were already flat-corrected, restore their pre-flat working state
before applying a newly built master; do not flat-correct them a second time.

## Preview combined images in DS9

Keep the combiner instance that generated the products you want to inspect:

```python
# Initialize stacking for the registered working images.
combiner: Combine = Combine(
    image_collection=analysis_images,
    telescope=SARA_CT,
)

# Sigma-clip and average Light frames separately for each field and filter.
analysis_images: ImageCollectionGroup = combiner.combine_light_or_flat(
    image_type="Light",
)

# Browse the newly combined images in DS9 without moving or modifying files.
combiner.preview_in_ds9()
```

The preview opens or reuses a dedicated DS9 window. Next and Previous browse
only the masters from this object's latest successful combine call, not all
working images. Retry handles loading failures. Quit closes the controls while
leaving DS9 open; the notebook waits until Quit. There is no Reject action.
On macOS the controls run in the existing isolated GUI process. Image values
are displayed directly; this preview does not compare Before/After or display
MASK/CMASK as separate frames. Calling preview again does not repeat stacking.
The same method also works after combining Bias, Dark, or Flat frames.

## Checks and limitations

Selection, master uniqueness, geometry, units, processing flags, relevant
EXPTIME values, and output collisions are checked before writing. Different
night directories use their own masters. A shared save_location must not cause
multiple children to write the same destination basename.

Dark input stacks must have a single common positive EXPTIME. Mixed dark
exposure times are not silently combined or normalized. Scaled dark correction
requires bias-corrected target and master. Flat and Light prior bias/dark states
must agree. No correction stage is automatically added or skipped for you.

Repeated corrections are rejected. To repeat a stage, restore working copies
from an appropriate earlier processing state; do not merely remove its flag.
Numerically non-finite corrected outputs are rejected without adding a new
pixel threshold. Each file is published atomically, but an I/O or numerical
failure in a later file does not roll back earlier completed outputs.

`save_location` on a processing method can direct results to another directory.
Returned selections follow the new corrected paths, retain masters, and put
correction markers in outputs, leaving original source data/headers unchanged.
The default workflow above deliberately updates reviewed working copies in the
same analysis directory.
