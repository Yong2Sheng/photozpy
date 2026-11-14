"""Associate source metadata with optional photometric measurements."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import TYPE_CHECKING, Literal, Self, TypeAlias

from .filter_magnitude import FilterMagnitude
from .source_info import SourceInfo, SourceRole, VALID_SOURCE_ROLES

if TYPE_CHECKING:
    from matplotlib.figure import Figure
    from pandas import DataFrame

    from ._zero_points import _ZeroPointStatistics


MagnitudeKind: TypeAlias = Literal["instrumental", "reference", "calibrated"]
MAGNITUDE_KINDS = ("instrumental", "reference", "calibrated")
SourceKey: TypeAlias = tuple[str, str]


class Sources:
    """Collect source definitions and their optional magnitude results.

    Parameters
    ----------
    *source_infos : SourceInfo
        Definitions of the sources to collect. Each (field name, source name)
        pair must be unique. Empty collections are allowed so filtering can
        return an empty result.

    Notes
    -----
    Each source starts with three undefined results (``None``): instrumental,
    reference/catalog, and calibrated magnitudes. A FilterMagnitude may itself
    contain NaN for individual filters; this differs from an absent result.

    Integer indexing, name lookup, slicing, filtering, and iteration all return
    Sources objects, including single-source collections. These selections
    share measurement records with their parent: setting or clearing a result
    through either collection is visible through the other. Metadata and
    FilterMagnitude objects are stored by reference, without copying.

    Photometry is performed separately. Explicit zero-point, diagnostic plotting,
    and target-calibration methods operate on the stored results, without
    converting magnitude systems. Setting a result alone never calibrates it.
    """

    def __init__(self, *source_infos: SourceInfo) -> None:
        self._source_infos = tuple(source_infos)
        self._zero_point_statistics: _ZeroPointStatistics | None = None
        self._magnitudes: dict[
            SourceKey, dict[MagnitudeKind, FilterMagnitude | None]
        ] = {}

        for info in self._source_infos:
            if not isinstance(info, SourceInfo):
                raise TypeError(
                    "Each source must be a SourceInfo object. "
                    "For an iterable, use Sources.from_list(source_infos)."
                )
            key = self._key(info)
            if key in self._magnitudes:
                raise ValueError(
                    f"Duplicate source in observation field: {key!r}."
                )
            self._magnitudes[key] = {
                kind: None for kind in MAGNITUDE_KINDS
            }

    @staticmethod
    def _key(info: SourceInfo) -> SourceKey:
        return (info.observation_field.field_name, info.source_name)

    @classmethod
    def from_list(cls, source_infos: Iterable[SourceInfo]) -> Self:
        """Create a collection from an iterable of source definitions."""
        return cls(*source_infos)

    @classmethod
    def combine(cls, *collections: Sources) -> Self:
        """Combine collections while sharing their existing measurement records.

        Preserve input/source order and all three result slots. Updating the
        combined view updates the original collections, just like a selection.
        Repeated selections of the same shared record are deduplicated. Separate
        records claiming the same (field_name, source_name) raise ValueError
        instead of arbitrarily choosing one owner's results. Source roles may
        coexist within a field; role is not part of the source identity.
        """
        infos: list[SourceInfo] = []
        records: dict[SourceKey, dict[MagnitudeKind, FilterMagnitude | None]] = {}
        for collection in collections:
            if not isinstance(collection, Sources):
                raise TypeError("Each collection must be a Sources object.")
            for info in collection.source_infos:
                key = cls._key(info)
                record = collection._magnitudes[key]
                if key in records:
                    if records[key] is not record:
                        raise ValueError(f"Conflicting independent source records: {key!r}.")
                    continue
                infos.append(info)
                records[key] = record
        combined = cls(*infos)
        combined._magnitudes = records
        return combined

    @property
    def source_infos(self) -> tuple[SourceInfo, ...]:
        """Return source metadata in collection order."""
        return self._source_infos

    @property
    def source_names(self) -> tuple[str, ...]:
        """Return names in order; names may repeat across different fields."""
        return tuple(info.source_name for info in self._source_infos)

    @property
    def nsources(self) -> int:
        """Return the number of sources (also available through len)."""
        return len(self)

    @property
    def source_info(self) -> SourceInfo:
        """Return metadata when this collection contains exactly one source."""
        return self.get_info()

    def __len__(self) -> int:
        return len(self._source_infos)

    def __repr__(self) -> str:
        if len(self) != 1:
            return (
                f"Sources(nsources={len(self)}, "
                f"source_names={self.source_names!r})"
            )

        info = self._source_infos[0]
        coordinates = info.source_coordinate.icrs.to_string(
            "hmsdms", precision=2, pad=True, sep=":"
        )
        results = self._magnitudes[self._key(info)]
        statuses = ", ".join(
            f"{kind}={'undefined' if results[kind] is None else 'defined'}"
            for kind in MAGNITUDE_KINDS
        )
        return (
            "Sources(\n"
            f"  source_name={info.source_name!r},\n"
            f"  coordinates (ICRS): {coordinates},\n"
            f"  magnitudes: {statuses}\n"
            ")"
        )

    def __iter__(self) -> Iterator[Self]:
        """Iterate over single-source collections sharing these results."""
        for index in range(len(self)):
            yield self[index]

    def __getitem__(self, index: int | slice | str) -> Self:
        """Select by position, slice, or an unambiguous source name.

        Missing names raise KeyError. Ambiguous names raise ValueError; use
        ``for_field(field_name)[source_name]`` to disambiguate.
        """
        if isinstance(index, str):
            infos = (self.get_info(index),)
        elif isinstance(index, slice):
            infos = self._source_infos[index]
        elif isinstance(index, int):
            infos = (self._source_infos[index],)
        else:
            raise TypeError("Source index must be an integer, slice, or name.")
        return self._select(infos)

    def _select(self, infos: Iterable[SourceInfo]) -> Self:
        selected = type(self)(*infos)
        selected._magnitudes = {
            self._key(info): self._magnitudes[self._key(info)]
            for info in selected.source_infos
        }
        return selected

    def get_info(
        self,
        source_name: str | None = None,
        *,
        field_name: str | None = None,
    ) -> SourceInfo:
        """Resolve one source, optionally restricting its observation field.

        Omitting source_name is useful for a single-source collection.
        No match raises KeyError; more than one match raises ValueError.
        """
        if source_name is not None and not isinstance(source_name, str):
            raise TypeError("source_name must be a string or None.")
        if field_name is not None and not isinstance(field_name, str):
            raise TypeError("field_name must be a string or None.")
        matches = [
            info for info in self._source_infos
            if (source_name is None or info.source_name == source_name)
            and (field_name is None
                 or info.observation_field.field_name == field_name)
        ]
        if not matches:
            raise KeyError(
                f"No source matches name={source_name!r}, field={field_name!r}."
            )
        if len(matches) != 1:
            raise ValueError(
                "Selection is ambiguous; specify source_name and, "
                "if necessary, field_name to select exactly one source."
            )
        return matches[0]

    def for_field(self, field_name: str) -> Self:
        """Select sources with the given observation-field name."""
        if not isinstance(field_name, str):
            raise TypeError("field_name must be a string.")
        return self._select(
            info for info in self._source_infos
            if info.observation_field.field_name == field_name
        )

    def by_role(self, source_role: SourceRole) -> Self:
        """Select science or standard sources without changing their order."""
        if source_role not in VALID_SOURCE_ROLES:
            raise ValueError("source_role must be 'science' or 'standard'.")
        return self._select(
            info for info in self._source_infos
            if info.source_role == source_role
        )

    @staticmethod
    def _check_kind(kind: MagnitudeKind) -> None:
        if kind not in MAGNITUDE_KINDS:
            raise ValueError(
                "kind must be 'instrumental', 'reference', or 'calibrated'."
            )

    def get_magnitudes(
        self,
        *,
        kind: MagnitudeKind,
        source_name: str | None = None,
        field_name: str | None = None,
    ) -> FilterMagnitude | None:
        """Return a source's selected result, or None when not yet defined."""
        self._check_kind(kind)
        info = self.get_info(source_name, field_name=field_name)
        return self._magnitudes[self._key(info)][kind]

    def set_magnitudes(
        self,
        *,
        kind: MagnitudeKind,
        magnitudes: FilterMagnitude | None,
        source_name: str | None = None,
        field_name: str | None = None,
    ) -> None:
        """Attach or replace a result; pass None to clear it.

        The measurement must use the exact Telescope instance associated with
        the source's observation field. Instrumental results must be labeled
        'instrumental'; reference and calibrated results must be AB or Vega.
        These checks do not convert systems or calculate any values.
        """
        self._check_kind(kind)
        info = self.get_info(source_name, field_name=field_name)
        if magnitudes is not None:
            if not isinstance(magnitudes, FilterMagnitude):
                raise TypeError("magnitudes must be a FilterMagnitude or None.")
            if magnitudes.telescope is not info.observation_field.telescope:
                raise ValueError(
                    "Magnitudes must use the observation field's Telescope "
                    "instance."
                )
            allowed_systems = (
                ("instrumental",) if kind == "instrumental" else ("AB", "Vega")
            )
            if magnitudes.magnitude_system not in allowed_systems:
                raise ValueError(
                    f"Magnitude system {magnitudes.magnitude_system!r} "
                    f"is incompatible with kind={kind!r}."
                )
        self._magnitudes[self._key(info)][kind] = magnitudes

    def show_magnitudes(self, *, kind: MagnitudeKind = "calibrated") -> None:
        """Print one compact magnitude summary per source in collection order.

        Parameters
        ----------
        kind : {"instrumental", "reference", "calibrated"}, optional
            Result slot to display. Default is ``"calibrated"``.

        Notes
        -----
        Show the source name and field once, followed by the magnitude-system
        label and FilterMagnitude.to_table's filter, magnitude, error and
        significance columns. Missing values remain ``--`` and upper-limit
        errors retain the existing ``UL (-99)`` display. An undefined result
        is reported explicitly without creating a placeholder measurement.

        This method only formats existing results. It does not calculate,
        calibrate, modify or save measurements. Use selections such as
        ``by_role("science")`` before calling to restrict the displayed sources.
        """
        self._check_kind(kind)
        summaries: list[str] = []

        for info in self._source_infos:
            magnitudes = self._magnitudes[self._key(info)][kind]
            heading = (
                f"{info.source_name}\n"
                f"Field: {info.observation_field.field_name}"
            )
            table = (
                magnitudes.to_table() if magnitudes is not None
                else f"No {kind} magnitudes available."
            )
            summaries.append(f"{heading}\n{table}")

        print("\n\n".join(summaries) if summaries else "No sources to display.")

    def calculate_zero_points(
        self, sigma: float = 3.0, *, min_significance: float = 3.0
    ) -> DataFrame:
        """Calculate and store per-filter standard-star zero points.

        Use standard-role sources with finite reference/instrumental magnitudes,
        no instrumental upper-limit flag, and detection significance at least
        ``min_significance``. Preserve the legacy method: reference minus
        instrumental magnitude, one pass of scipy.stats.zscore filtering at
        ``sigma``, then an unweighted arithmetic mean. A singleton or constant
        sample is retained rather than rejected by undefined z-scores.

        Return an independent pandas DataFrame with per-filter ``zero_point``
        and population ``scatter`` (both numerical magnitudes), sample counts,
        and magnitude-system labels. Scatter describes retained stars; it is
        not propagated as a zero-point uncertainty. Neither catalog errors nor
        instrumental errors weight the mean. No exposure-time correction,
        color transformation, extinction or aperture correction is applied.

        All standards must share one Telescope instance and reference magnitude
        system. Science-role sources in a mixed collection are ignored. A filter
        without eligible standards has NaN results; no eligible filter raises
        ValueError. Source magnitude slots are never changed by this method.
        """
        from ._zero_points import calculate_zero_points

        result = calculate_zero_points(self, sigma, min_significance)
        self._zero_point_statistics = result
        return result.summary.copy(deep=True)

    def _get_zero_point_statistics(self) -> _ZeroPointStatistics:
        """Require an existing result that still matches the standard inputs."""
        from ._zero_points import standard_signature

        result = self._zero_point_statistics
        if result is None:
            raise RuntimeError("Call calculate_zero_points() first.")
        if result.signature != standard_signature(self):
            raise RuntimeError(
                "Standard-star measurements changed; call calculate_zero_points() again."
            )
        return result

    @property
    def zero_points(self) -> DataFrame:
        """Return a copy of the current per-filter summary; values are in mag."""
        return self._get_zero_point_statistics().summary.copy(deep=True)

    @property
    def zero_point_details(self) -> DataFrame:
        """Return per-star/filter magnitudes, zero points and selection reasons.

        Numeric magnitude columns are in mag. The returned DataFrame is a
        defensive copy and includes excluded measurements for inspection.
        """
        return self._get_zero_point_statistics().details.copy(deep=True)

    def plot_zero_point_statistics(
        self,
        *,
        reference_offset: float = 0.1,
        save_path: str | Path | None = None,
    ) -> Figure:
        """Show one panel per measured filter using the stored zero points.

        Plot all finite individual zero points, with retained points in blue,
        rejected points as orange crosses, and the retained average in green.
        Each panel reports retained/rejected counts and, when present, the
        number of unavailable offsets that have no finite value to plot.
        Red lines at average +/- ``reference_offset`` mag are visual guides,
        not clipping thresholds or uncertainty bounds. Return the Figure for
        further inspection. Save only when ``save_path`` is explicitly supplied.
        Call calculate_zero_points first; modified standard inputs require a
        fresh calculation. Plotting never recalculates or changes the result.
        """
        from ._zero_points import plot_zero_point_statistics

        return plot_zero_point_statistics(
            self._get_zero_point_statistics(), reference_offset, save_path
        )

    def calibrate_targets(self, *, standard_sources: Sources) -> DataFrame:
        """Apply standard-star zero points to all science sources in this view.

        Reuse the supplied collection's inspected zero points. If none have
        been calculated, calculate them with the default thresholds first.
        Changed standard inputs require an explicit recalculation. Standards
        may share the science observation field or occupy a separate field.

        Add each filter's zero point to the original instrumental magnitude
        and write a new FilterMagnitude to the calibrated slot. Preserve the
        instrumental errors and significance exactly, including the -99 mag
        upper-limit error flag. Catalog errors and zero-point uncertainty are
        deliberately not propagated. The output system follows the reference
        magnitudes; there is no AB/Vega conversion.

        Return a pandas DataFrame of all science source/filter results, with
        magnitude/error values in mag. Missing instrumental values remain NaN;
        a source without instrumental results has no calibrated result. A
        finite instrumental value without a usable zero point, or a Telescope
        mismatch, raises before any calibrated slot is updated. Repeated calls
        always use instrumental magnitudes, never already calibrated values.
        """
        from ._zero_points import calibrate_targets

        if not isinstance(standard_sources, Sources):
            raise TypeError("standard_sources must be a Sources object.")
        if standard_sources._zero_point_statistics is None:
            standard_sources.calculate_zero_points()
        result = standard_sources._get_zero_point_statistics()
        return calibrate_targets(self, result)
