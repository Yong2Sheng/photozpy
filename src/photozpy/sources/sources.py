"""Associate source metadata with optional photometric measurements."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import Literal, Self, TypeAlias

from .filter_magnitude import FilterMagnitude
from .source_info import SourceInfo, SourceRole, VALID_SOURCE_ROLES


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

    This class stores results only. It does not perform photometry, magnitude
    system conversion, zero-point estimation, or calibration, and does not
    implement the legacy Sources API.
    """

    def __init__(self, *source_infos: SourceInfo) -> None:
        self._source_infos = tuple(source_infos)
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
