# Cosmic-ray removal before registration

Use `Artifacts` after bias, dark, and flat correction and before Registration.
The implementation retains Yong's original `ccdproc.cosmicray_lacosmic` call:
GAIN and RDNOISE come from each FITS header, and `gain_apply=False` keeps the
pixel values in ADU. CCDProc 2.5.1 still labels its returned CCDData as electrons;
Photozpy restores the ADU unit explicitly, following the original BUNIT policy.
It does not multiply or divide the cleaned pixels by gain.

## Notebook usage

`analysis_images` is the ImageCollectionGroup returned by flat correction.
Only Light frames are cleaned, including independent standard fields;
calibration frames and master frames remain in the returned selection.

```python
from photozpy import Artifacts, ImageCollectionGroup

# Remove cosmic rays from reduced Light working images before registration.
artifacts: Artifacts = Artifacts(
    image_collection=analysis_images,
    overwrite=True,
)

analysis_images: ImageCollectionGroup = artifacts.remove_cosmic_ray(
    plot_mask=True,
    preview=True,
    saturation_adu=65000.0,  # Saturation level in the displayed image's ADU units.
    headers_values={"IMTYPE": "Light"},
)
```

`plot_mask` saves `*_cr_mask.png` beside each FITS output with a fixed grayscale
range: black=0 (unmasked), white=1 (masked). Existing masked pixels, if any, are
included. Display this PNG explicitly if needed; it is not shown in the notebook
automatically. Set plot_mask=False to skip PNG creation; the FITS CMASK is still
saved.

Set preview=True to start an automatic DS9 comparison after the entire cleaning
batch succeeds. It never loads images into DS9 while cleaning is running.
The default preview=False only cleans/saves files. The old plot_before_after
and norm_percent options have been removed; omit them from existing cells.
Processing options are keyword-only to prevent old positional calls from being
silently interpreted as different options.

By default, save_location=None updates the selected working files in place;
overwrite=True is required. To preserve those inputs, pass save_location with
another directory. Raw files are protected only by selecting working copies as
the input collection. Constructors do not perform cleaning or write files.

## Backups and automatic comparison

In-place cleaning first copies all selected Light FITS into a dedicated
before_cosmic_ray/ subdirectory. Working results remain directly in analysis/:

```text
analysis/
    light.fits                  # Cleaned working image, including CMASK.
    light_cr_mask.png           # Optional binary PNG.
    before_cosmic_ray/
        light.fits              # Exact copy immediately before cosmic cleaning.
```

Use backup_dir to override the backup location. Default locations follow each
input directory, so separate nights retain their own backups. Existing backup
files are never overwritten, even with overwrite=True. All copies must succeed
before the first image is cleaned. Backups are excluded from the returned
collection, and are never automatically deleted. After checking, users may
delete the entire backup directory to recover disk space. With save_location
and no explicit backup_dir, unchanged inputs already serve as the Before images;
no redundant backup copy is made. artifacts.backup_paths maps outputs to Before.

preview=True connects to or starts a dedicated DS9 window named
photozpy-cr-preview, requiring DS9 and the XPA tools. A small control window
offers Next, Previous, Sync views, Retry, and Quit; no Good/Reject file moves
occur. The notebook cell waits for Quit. On macOS the controls use the same
isolated process approach as ImageReviewer to keep Tk out of the kernel.

The layout is explicitly three columns in one row: left Before (frame 1),
middle After (frame 2), and right Mask (frame 3). The initial zoom-to-fit is
computed after all three tiles exist, then shared across the frames.
Before/After start with shared
display settings; the mask remains linear with limits 0..1. Adjust either
science frame in DS9 and click Sync views to copy its scale limits, stretch,
contrast, colormap, zoom, and position to its partner. Geometry is also copied
to the mask, while its binary scale remains independent. You can also select
Mask, zoom and pan to a candidate cosmic ray, then click Sync views: its geometry
is copied to Before/After without changing either science frame's intensity
scale, limits, contrast, or colormap. The selected reference frame remains active
after synchronization. Next/Previous read the currently selected frame's live
zoom and image-coordinate center immediately before loading. That geometry is
explicitly restored on each new image, overriding DS9's older per-file view
history. The source frame stays selected. This works even if you changed zoom
or position without clicking Sync views. There is no new zoom-to-fit operation.
Only frame/zoom/pan are queried during navigation; science scale, limits, and
colors are not read or re-synchronized. Use Sync views to match intensity settings.
Loads and geometry changes use individual XPA commands. The attempted -exec
batch route was removed after it left real DS9 frames empty despite command
success. Initial loading now also verifies that each frame contains FITS data
before enabling Next. Intermediate resizing may still be visible; eliminating
that effect has not been established for this DS9 build.
Avoid DS9 global scale locks, which include the mask.

Quit closes the controls, keeps DS9 open, and retains all backups. Reopen the
same completed batch without re-cleaning using:

```python
artifacts.preview_in_ds9()
```

Do this before later in-place pipeline steps alter the After files. No large
image arrays are retained in the notebook; the viewer loads one pair at a time
and DS9 reads CMASK directly from the cleaned FITS extension. No Python-side
science-image read or temporary mask FITS is needed for preview. A DS9
failure does not undo cleaning or remove backups: artifacts.image_collection
still contains the cleaned working images. Fix DS9 and reopen the preview
rather than running cosmic-ray removal again.

## Detection options and inspection

No new detection thresholds are imposed. Extra keyword arguments are forwarded
to CCDProc, for example sigclip or objlim. Omitted arguments use the installed
library defaults; CCDProc 2.5.1 uses sigclip=4.5, sigfrac=0.3, objlim=5.0, niter=4.
GAIN, readnoise, and gain_apply cannot be overridden through these keywords.
Use saturation_adu to specify saturation in the same ADU units as the FITS
image. The example explicitly selects 65000 ADU. Photozpy reads each frame's
GAIN and converts that threshold internally: gain=1 gives 65000 electrons;
gain=2 gives 130000 electrons. Users do not need to perform this conversion.
Output pixels remain ADU; this option does not clip the image at that value.

The new parameter defaults to None to preserve existing calls and the CCDProc
default. Set saturation_adu explicitly to use an ADU-based threshold. Legacy
satlevel in electrons remains accepted only when saturation_adu is None;
providing both raises an error before processing. ADU thresholds must be
positive finite numbers.

Inspect representative masks and compare images in DS9, especially bright or
undersampled stellar cores, before processing an entire night. A detected pixel
is a candidate, not proof of a cosmic ray. To select a subset, supply additional
FIELD or FILTER criteria in headers_values, or provide an explicitly selected
collection. Do not clean that subset twice when later processing the rest.

Outputs save cosmic-ray detections in a separate binary CMASK extension.
The cleaned data participates in registration and Light combination normally.
CMASK is used only for inspection and registration star detection; it is never
merged into the validity MASK used by Combine. A pre-existing input validity
MASK remains separate. FITS cards record CRCORR='Yes', CRMETHOD='LACosmic',
and CRNPIX (newly detected pixel count), also available in
artifacts.cosmic_ray_counts. Registration resamples CMASK separately and adds
its uncovered footprint to the validity MASK.

## Preconditions and failures

Inputs must be finite 2D primary ADU images with positive finite GAIN,
nonnegative finite RDNOISE, and BIASCORR/DARKCORR/FLATCORR='Yes'. Existing ALIGN
markers or CRCORR='Yes' reject cleaning in the wrong order or repeat cleaning.
These checks do not automatically perform any missing calibration step.

Only an optional MASK extension is supported. Uncertainty and unrelated
extensions are rejected, rather than silently discarded. The original output
policy of not creating an uncertainty model is retained.

Header prerequisites, extension names, and output collisions are checked before
the batch starts. A later numerical, plotting, or I/O failure does not roll back
earlier outputs; a plot failure can occur after that frame's cleaned FITS was
written. Each FITS publication is atomic. Restore an appropriate earlier working
copy to repeat cleaning; merely removing CRCORR is not a scientific rollback.
