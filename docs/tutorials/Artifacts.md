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
# Import the public cleaning API and collection type.
from photozpy import Artifacts, ImageCollectionGroup

# Remove cosmic rays from reduced Light working images before registration.
artifacts: Artifacts = Artifacts(
    image_collection=analysis_images,
    overwrite=True,
)

# Back up every selected Light frame before starting two cleaning processes.
analysis_images: ImageCollectionGroup = artifacts.remove_cosmic_ray(
    workers=2,  # Use 1 for serial processing, or try 4 after measuring memory.
    plot_mask=False,  # Skip PNGs; the FITS CMASK is always saved.
    preview=False,  # Open DS9 separately after the batch succeeds.
    saturation_adu=65535.0,  # Saturation threshold in ADU for this detector.
    headers_values={"IMTYPE": "Light"},
)
```

The cell assumes `analysis_images` already contains reduced working images.
Use your detector's established saturation threshold; worker count does not
change any detection parameter. Restart the notebook kernel after updating
Photozpy so the kernel and spawned processes import the same source version.

## Image-level multiprocessing and memory

`workers=1` is the default and preserves the serial path. Set `workers=2` or
`workers=4` explicitly to process independent images concurrently. The local
`spawn` context works with an imported Photozpy API in a Jupyter kernel; no
worker function belongs in a notebook cell, and the global multiprocessing
start method is unchanged. In a standalone Python script, put the batch call
under `if __name__ == "__main__":` to prevent recursive process creation.

The parent checks all selected inputs, output collisions, and backup collisions,
then completes every backup before starting workers. Each worker reads, cleans,
and saves one image at a time using the existing single-image implementation.
Only paths, scalar options, and pixel counts cross process boundaries; the
parent does not preload images or collect cleaned arrays. Nonscalar CCDProc
options (including arrays and lists) require `workers=1`. The process count is
capped at the selected image count and never defaults to all available CPUs.

Each spawned process imports its own scientific libraries and holds its own
working arrays, so memory increases with worker count. Start with 2 and measure
before increasing to 4; spawning workers can outweigh the benefit for very
small batches. Work is dispatched one image at a time, and completed results
update the parent's progress bar. The returned collection, counts, and backup
associations retain input order even when workers finish out of order.

Worker mask PNGs use Matplotlib's noninteractive Agg backend. DS9 starts only
in the parent workflow after the entire batch succeeds. On a worker exception
or notebook interrupt, outstanding processes are terminated and joined before
the error is propagated. Collection and success metadata are published only
after successful completion. Already written outputs and all backups remain;
failure recovery is the same manual restore procedure as for serial processing.

### Measured speed and memory

On an Apple Silicon M5 with 64 GiB RAM and 18 visible CPUs, a 12-image subset
of real pre-cleaning backups (three fields, four filters, 2048 x 2048 float32)
was processed twice per setting with CCDProc 2.5.1 and Astro-SCRAPPY 1.2.0.
The execution order was 1/2/4 then 4/2/1 workers. All runs used temporary input
copies, in-place cleaning with backups, `plot_mask=False`, `preview=False`, and
`saturation_adu=65535.0`. Copying the benchmark inputs was outside the timer;
batch preflight, backups, process startup, cleaning, saving, and collection
updates were included. No internal OpenMP speedup was assumed.

| Workers | Wall time, two runs (s) | Mean CPU time (s) | CPU / wall | Peak aggregate RSS (MiB) | Speedup |
| --- | --- | --- | --- | --- | --- |
| 1 | 30.29 / 30.45 | 30.29 | 1.00 | 686 | 1.00x |
| 2 | 18.00 / 18.02 | 34.57 | 1.92 | 1,983 | 1.69x |
| 4 | 10.11 / 10.04 | 38.43 | 3.82 | 3,224 | 3.02x |

CPU time sums parent and completed worker user/system time. Memory is the
maximum sampled sum of parent and active worker RSS, sampled every 20 ms;
the small multiprocessing resource-tracker process is excluded. Shared pages
can be counted more than once, and short peaks may be missed. The same parent
interpreter was reused between runs, retaining library imports and allocator
memory. Reported peaks are the larger of the two runs and include all imported
libraries and working arrays, not just FITS payloads.

All 72 outputs matched the serial reference exactly for cleaned arrays, masks,
CRNPIX, and detector/processing metadata; each backup matched its source bytes.
For lower memory use, start with 2 workers (about 1.94 GiB here, 41% less wall
time than serial). With sufficient memory available, 4 workers reached about
3.15 GiB and reduced wall time by 67%. These are small-batch measurements,
including startup; the complete 131-image night was not rerun in this benchmark.

`plot_mask` saves `*_cr_mask.png` beside each FITS output with a fixed grayscale
range: black=0 (not detected), white=1 (detected). It shows new cosmic-ray
detections, matching CMASK; the existing validity MASK remains separate.
Display this PNG explicitly if needed; it is not shown in the notebook
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
# Review the completed batch without cleaning the images again.
artifacts.preview_in_ds9()
```

Do this before later in-place pipeline steps alter the After files. No large
image arrays are retained in the notebook; the viewer loads one pair at a time
and DS9 reads CMASK directly from the cleaned FITS extension. No Python-side
science-image read or temporary mask FITS is needed for preview. A DS9
failure does not undo cleaning or remove backups: artifacts.image_collection
still contains the cleaned working images. Fix DS9 and reopen the preview
rather than running cosmic-ray removal again.

## Quality review after cleaning

Use a separate review step to reject entire exposures after inspecting their
calibrated Before/After/Mask views. Run it immediately after cleaning and before
Registration, using the same `artifacts` object:

```python
# Review the completed cleaning batch and apply exposure-quality decisions.
analysis_images = artifacts.review_in_ds9()

# Pass this updated selection to the subsequent registration step.
print(f"Retained cleaned exposures: {len(artifacts.cosmic_ray_counts)}")
```

No `raw_dir` argument is needed. Working images must follow the existing
`root/analysis/image.fits` convention. For each working image, Photozpy infers
the original as `root/image.fits` and the rejection destination as
`root/rejected/image.fits`. Multiple nights infer their roots independently.
Unchanged filenames, DATE-OBS, EXPTIME, and image geometry must match the raw
original, Before, and After. Missing or mismatched originals, other directory
layouts, registered images, and destination collisions stop review rather than
guessing a source or overwriting a previous rejection. Symlinks are rejected.

The review controls provide **Keep**, **Reject**, **Previous**, **Sync views**,
**Retry**, and **Quit**. Reject means the original exposure is unsuitable for
future analyses, not merely that its cosmic-ray mask looks different. Decisions
are staged while the window is open; Previous undoes the preceding decision.
Normal Quit (including closing the controls) applies pending rejects. Quitting
early keeps unvisited exposures. A GUI process failure or interrupt applies no
staged decisions. DS9 stays open. As with preview, macOS keeps Tk in an isolated
process and the notebook waits for the window to close.

For each rejected exposure, Photozpy retains the exact raw original in
`rejected/` and writes an adjacent `image.fits.rejection.json` audit record with
the paths, UTC timestamp, stage, and visual-rejection reason. Existing originals
and audit records in `rejected/` are never overwritten. Only after retention
succeeds is the original removed from root, followed by the known derived
files: the cleaned FITS, its associated Before FITS, and its optional mask PNG.
Regions and other products are never searched for or removed. Accepted files
are unchanged. The returned selection and the object's counts/backup mappings
exclude rejected exposures while preserving other images and their order.

Root-level `ImageCollectionGroup(root)` discovery is nonrecursive, so the
archived raw file is excluded when starting a new pipeline run. Always assign
the review result to `analysis_images`; older collection objects remain snapshots.
If every image is rejected, an empty selection is returned and remains empty on
refresh. Do not proceed to registration when no usable Light images remain.

Rejection file operations are sequential, not a whole-session transaction.
If cleanup fails after raw retention, the original stays in `rejected/`, some
derived files may remain, and the exception is propagated. The object's
`artifacts.image_collection` already excludes exposures whose originals moved;
recover that selection and inspect the recorded paths before continuing.
Earlier completed rejections remain applied. Restore an archived original to
root and regenerate its analysis products if a completed rejection was mistaken.

`preview_in_ds9()` and `preview=True` remain inspection-only and never reject,
move, or delete images. Reopening preview after review shows surviving pairs.

## Detection options and inspection

No new detection thresholds are imposed. Extra keyword arguments are forwarded
to CCDProc, for example sigclip or objlim. Omitted arguments use the installed
library defaults; CCDProc 2.5.1 uses sigclip=4.5, sigfrac=0.3, objlim=5.0, niter=4.
GAIN, readnoise, and gain_apply cannot be overridden through these keywords.
Use saturation_adu to specify saturation in the same ADU units as the FITS
image. The example explicitly selects 65535 ADU. Photozpy reads each frame's
GAIN and converts that threshold internally: gain=1 gives 65535 electrons;
gain=2 gives 131070 electrons. Users do not need to perform this conversion.
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
