# Sources: metadata and optional magnitude results

`SourceInfo` describes one identified source. `Sources` collects those
definitions and associates each one with optional photometric results. No
public `Source` wrapper is needed.

## Create sources before measuring anything

```python
from astropy.coordinates import SkyCoord
from photozpy.sources import ObservationField, SourceInfo, Sources
from photozpy.telescope import SARA_CT

field = ObservationField(
    field_name="example_field",
    telescope=SARA_CT,
    nominal_pointing=SkyCoord(150.0, 2.0, unit="deg"),
    file_pattern="example_field*",
)
science_info = SourceInfo(
    source_name="example_target",
    source_role="science",
    source_coordinate=SkyCoord(150.01, 2.01, unit="deg"),
    observation_field=field,
)
sources = Sources(science_info)
assert sources[0].get_magnitudes(kind="instrumental") is None
```

For multiple definitions, use `Sources(info_a, info_b)` or
`Sources.from_list(existing_infos)`. A dedicated standard field can exist as an
`ObservationField` before its stars are identified; do not invent placeholder
stars to add it to this collection.

## Attach a measurement later

The following values are illustrative, not measurements from real data.

```python
from photozpy.sources import FilterMagnitude

instrumental = FilterMagnitude(
    telescope=SARA_CT,
    filter_names=("SDSS g", "SDSS r"),
    magnitudes=(-10.2, -10.5),
    errors=(0.03, 0.02),
    significance=(36.0, 54.0),
    magnitude_system="instrumental",
)
sources.set_magnitudes(
    source_name="example_target",
    kind="instrumental",
    magnitudes=instrumental,
)
target = sources["example_target"]
print(target.source_info.source_coordinate)
print(target.get_magnitudes(kind="instrumental").to_table())
```

Each source has independent `instrumental`, `reference`, and `calibrated`
slots. `None` means the entire result is undefined; `NaN` inside an existing
`FilterMagnitude` means a particular filter has no value. No zero-filled tables
are created. To clear a slot, pass `magnitudes=None` to `set_magnitudes`.

Measurements must use the same Telescope instance as their observation field.
The instrumental slot accepts the `instrumental` system; the reference and
calibrated slots accept AB or Vega. These labels are validated without any
system conversion, magnitude calculation, or automatic calibration. Assigning
a reference measurement does not calibrate the instrumental measurement.

## Select and iterate

```python
first = sources[0]                       # A single-source Sources
subset = sources[:2]                    # Another Sources
science = sources.by_role("science")
in_field = sources.for_field("example_field")

for single in sources:
    print(single.source_info.source_name)
```

Selections share result records with the parent. Replacing or clearing a result
through `first`, `subset`, or `science` updates the same source in `sources`.
SourceInfo and FilterMagnitude objects are held by reference, not copied.
Separately constructing `Sources(*sources.source_infos)` creates new,
undefined result slots; it is not a copy of the measurements.

## Combine existing collections without losing results

```python
# Create one view while preserving catalog results and sharing future measurements.
all_sources = Sources.combine(standard_sources, science_sources)
```

The combined view preserves input order and shares all three measurement slots
with the originals. Repeated selections of the same record are deduplicated;
independent records with the same field/source identity raise ValueError rather
than silently choosing which results to retain. Standards and science targets
can coexist in the same ObservationField. Photometry accepts this combined view
or a tuple of existing Sources collections directly.

Identity is `(field_name, source_name)`. Duplicate pairs are rejected. Names
can repeat in different fields, but a name-only lookup must be unambiguous:

```python
target = sources.for_field("example_field")["example_target"]
# Equivalent metadata lookup:
info = sources.get_info("example_target", field_name="example_field")
```

Filtering may return an empty collection. Accessing one source's metadata or
results requires exactly one match; missing matches raise KeyError and
ambiguous matches raise ValueError. Slicing follows normal Python rules.

## Display compact magnitude summaries

```python
# Show each source once with its calibrated filter/magnitude/error/significance table.
sources.show_magnitudes(kind="calibrated")

# Inspect instrumental results using the same compact layout.
sources.show_magnitudes(kind="instrumental")

# Select a role before displaying a mixed standard/science collection if needed.
sources.by_role("standard").show_magnitudes(kind="reference")
```

The method prints sources in collection order and returns None. Each block shows
the source name, observation field, magnitude system and the existing
FilterMagnitude.to_table output. Missing filter values appear as `--`, and
upper-limit error flags as `UL (-99)`. Undefined result slots are reported
without creating measurements. Display does not perform photometry/calibration,
alter values, or write files. The default kind is `"calibrated"`.

This is a new API, not a compatibility layer for the legacy Sources class.
Legacy consumers still need their own migration; refactored photometry accepts
these collections and writes instrumental results through the shared slots.

## Calculate zero points and calibrate targets explicitly

```python
# Calculate zero points from standard-role sources with reference/instrumental results.
zero_point_table = standard_sources.calculate_zero_points(sigma=3.0)

# Inspect all star/filter offsets and the reasons for excluding any measurements.
standard_zero_point_table = standard_sources.zero_point_details

# Review the stored per-filter statistics before calibration.
zero_point_figure = standard_sources.plot_zero_point_statistics()

# Apply those zero points to all science-role sources in this collection.
calibrated_table = sources.calibrate_targets(standard_sources=standard_sources)

# Display compact summaries without repeating source names for every filter.
sources.show_magnitudes(kind="calibrated")
```

These are explicit operations; attaching reference or instrumental measurements
does not automatically perform calibration. The numerical recipe remains the
legacy single-pass z-score selection and unweighted mean of reference minus
instrumental magnitudes. Calibration adds this mean to instrumental magnitudes
and copies their errors/significance unchanged. Catalog errors and zero-point
uncertainty are not propagated. Results are written only to calibrated slots.

The tables are pandas DataFrames with numerical magnitude values in mag;
stored magnitudes retain FilterMagnitude's unit-aware representation. Standard
and science sources can share an observation field. All participants must use
the same Telescope instance and the standards one reference magnitude system.
Changed standard measurements require recalculation before reusing their zero
points. See [Photometry.md](Photometry.md#standard-star-zero-points-and-target-calibration)
for selection rules, upper limits, missing filters and diagnostic definitions.
