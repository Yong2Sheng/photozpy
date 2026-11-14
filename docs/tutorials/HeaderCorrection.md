# HeaderCorrection

Header correction uses the refactored `ObservationField`, `Telescope`,
`ImageCollectionGroup`, and collection refresh helpers. No legacy Sources,
Telescope, or CollectionManager interface is used by this module.

## Example for a night with an independent standard field

The following assumes the notebook has already defined `data_root_dir`,
`FGL_field`, `SUMSS_field`, and `DLS_field`. All fields are supplied directly,
and the standard field can be supplied before any
standard stars have been identified.

```python
from photozpy import HeaderCorrection, ImageCollectionGroup

raw_images = ImageCollectionGroup(directories=data_root_dir)
output_dir = data_root_dir.parent / f"{data_root_dir.name}_headers"

headers = HeaderCorrection(
    image_collection=raw_images,
    observation_fields=(FGL_field, SUMSS_field, DLS_field),
    save_location=output_dir,
    overwrite=False,
)

# Write corrected copies; refuse to overwrite existing output files.
images = headers.correct_headers_by_filename()

# Update FILTER in those output copies, never in raw_images.
images = headers.correct_filter_headers(overwrite=True)
```

Use an output directory that does not already contain outputs from an earlier
run, or explicitly choose `overwrite=True` after reviewing the intended paths.
The constructor writes nothing. Every successful operation updates
`headers.image_collection` and returns the same result. An ImageCollectionGroup
input produces a group output; a single ImageFileCollection produces a single
collection. Output collections represent explicit file selections, not automatic
views of newly created files.

## Metadata rules

- Light filenames match ObservationField.file_pattern. A pattern without glob
  metacharacters is a case-sensitive substring; a pattern containing `*`, `?`,
  or `[` is matched against the complete basename as a glob.
- Supply all fields through `observation_fields`. No SourceInfo or Sources
  objects are required. Repeated fields reuse the same ObservationField instance.
- `FIELD` stores the field name; existing `OBJECT` is preserved. Calibration
  filenames identify Bias, Dark, and Flat; other files must match a field.
- Existing RA/DEC are retained. When both are missing, nominal field pointing
  is written in ICRS sexagesimal notation. `overwrite_pointing=True` explicitly
  replaces both. Nominal pointing is not an image's WCS-derived center.
- DATE-OBS is interpreted in the file's TIMESYS. The resulting start time is
  written to MJD-OBS without an exposure-midpoint adjustment or scale conversion.
  When TIMESYS is absent, supply `time_scale="utc"` (or another known scale)
  explicitly. A fallback that contradicts an existing TIMESYS is rejected.
- GAIN and RDNOISE use the CCD valid on each file's observation date. Gain is
  electron/ADU; read noise is per-pixel RMS in electron/pixel. No image data or
  uncertainty calculation is performed. Existing non-ADU BUNIT is rejected
  rather than silently relabeled.
- Supply `telescope_name="..."` if you want to set TELESCOP. Otherwise it is
  preserved; the code does not infer an instrument name from a Python variable.
- Filter correction uses Telescope.resolve_filter_alias. It requires IMTYPE to
  be populated first. Only Flat and Light filters are changed; other images are
  retained in the output. Raw SDSS_g/r/i/z aliases resolve to their existing
  SDSS primed canonical names. No photometric-system conversion is performed.

## Failure and output handling

Unknown or ambiguous field/filter matches, invalid dates, missing TIMESYS,
invalid card values, existing output files, and duplicate destination paths are
checked before writing the batch. Unknown field files may instead be explicitly
omitted with `unmatched="skip"`; inspect `headers.skipped_files` afterward.

All HDUs and image arrays are preserved by `fits_utilities.update_fits_header`.
Each file is written completely before it replaces its destination. The batch
is not a transaction: a later disk I/O failure can leave earlier outputs written.
A shared output directory for multiple input directories requires unique
basenames; conflicting output paths are rejected rather than overwritten.

Other pipeline stages have not yet been migrated to FIELD. In particular,
legacy selection by OBJECT must be reviewed when those stages are refactored;
header correction does not make those legacy consumers compatible automatically.
