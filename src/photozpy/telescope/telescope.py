from __future__ import annotations

from collections.abc import Sequence
from datetime import date
from itertools import combinations
from typing import Any

from astropy.units import Quantity

from .ccd import CCD
from .filterset import FilterSet, UnknownFilterAliasError


class Telescope:
    """
    Represent a telescope configuration composed of CCD definitions and filter sets.

    A Telescope manages:
        - one or more CCD definitions
        - one or more filter sets
        - optional telescope-level metadata

    During initialization, the input CCDs and filter sets are normalized into
    immutable tuples and validated for internal conflicts:

        - CCD definitions must not have overlapping valid date ranges
        - filter sets must not have overlapping resolvable filter names

    This class provides convenience methods to:
        - print CCD and filter-set summaries
        - query CCD field values by observation date
        - access common CCD parameters such as gain and read noise
        - list all standard filter names available in the telescope
        - resolve a filter alias across all stored filter sets
    """

    def __init__(
        self,
        ccds: CCD | Sequence[CCD],
        filter_sets: FilterSet | Sequence[FilterSet],
        metadata: dict[str, Any] | None = None,
    ) -> None:

        # check if input ccds is valid
        if isinstance(ccds, CCD):
            self._ccds = (ccds,)
        elif isinstance(ccds, Sequence) and not isinstance(ccds, (str, bytes)):
            self._ccds = tuple(ccds)
            if not all(isinstance(ccd, CCD) for ccd in self._ccds):
                raise TypeError("All elements in ccds must be CCD objects.")
        else:
            raise TypeError("ccds must be a CCD or a sequence of CCD objects.")

        # check if ccds have conflicts
        for left, right in combinations(self._ccds, 2):
            if left.conflicts_with(right):
                raise ValueError(
                    f"CCD {left.name!r} and {right.name!r} "
                    f"have an overlap on valid date range."
                )

        # check if input filter_sets is valid
        if isinstance(filter_sets, FilterSet):
            self._filter_sets = (filter_sets,)
        elif isinstance(filter_sets, Sequence) and not isinstance(filter_sets, (str, bytes)):
            self._filter_sets = tuple(filter_sets)
            if not all(isinstance(fs, FilterSet) for fs in self._filter_sets):
                raise TypeError("All elements in filter_sets must be FilterSet objects.")
        else:
            raise TypeError(
                "filter_sets must be a FilterSet or a sequence of FilterSet objects."
            )

        # check if filtersets have conflicts
        for left, right in combinations(self._filter_sets, 2):
            conflicts = left.find_conflicts(right)
            if conflicts:
                raise ValueError(
                    f"Filter sets {left.filterset_name!r} and {right.filterset_name!r} "
                    f"conflict on resolvable names: {conflicts!r}"
                )

        self._metadata = {} if metadata is None else dict(metadata)

    def print_ccd_summary(self) -> None:
        """
        Print summary tables for all CCD definitions in the telescope.
        This is a convenience method for interactive use. It calls
        ``print_summary()`` on each CCD stored in the telescope.
        """

        for ccd in self._ccds:
            ccd.print_summary()

    def print_filterset_summary(self) -> None:
        """
        Print summary tables for all filter sets in the telescope.
        This is a convenience method for interactive use. It calls
        ``print_summary()`` on each filter set stored in the telescope.
        """

        for fs in self._filter_sets:
            fs.print_summary()

    def ccd_value(
        self,
        field_name: str,
        obs_date: date | None = None,
    ) -> object | tuple[object, ...]:
        """
        Return CCD field value(s) from the telescope.

        If ``obs_date`` is not provided, return the requested field from all CCD
        definitions as a tuple. If ``obs_date`` is provided, return the field from
        the unique CCD definition valid on that date.

        Parameters
        ----------
        field_name : str
            Name of the CCD attribute to retrieve, such as ``"gain"`` or
            ``"rdnoise"``.
        obs_date : date | None, optional
            Observation date used to select the valid CCD definition. If ``None``,
            values from all CCD definitions are returned.

        Returns
        -------
        object | tuple[object, ...]
            A single field value if ``obs_date`` is provided, otherwise a tuple of
            values from all CCD definitions.

        Raises
        ------
        AttributeError
            If ``field_name`` is not a valid CCD attribute.
        ValueError
            If no CCD definition is valid for ``obs_date``, or if more than one CCD
            definition is valid for that date.
        """

        if not hasattr(self._ccds[0], field_name):
            raise AttributeError(f"CCD has no attribute {field_name!r}.")

        if obs_date is None:
            return tuple(getattr(ccd, field_name) for ccd in self._ccds)

        matched = tuple(
            ccd for ccd in self._ccds
            if ccd.valid_from <= obs_date <= ccd.valid_to
        )

        if not matched:
            raise ValueError(f"No CCD is valid for obs_date={obs_date!r}.")

        return getattr(matched[0], field_name)

    def gain(
        self,
        obs_date: date | None = None
    ) -> Quantity | tuple[Quantity, ...]:
        """
        Return CCD gain value(s) from the telescope.
        If ``obs_date`` is not provided, return the gain from all CCD definitions
        as a tuple. If ``obs_date`` is provided, return the gain from the CCD
        definition valid on that date.
        """
        return self.ccd_value("gain", obs_date)

    def rdnoise(
        self,
        obs_date: date | None = None
    ) -> Quantity | tuple[Quantity, ...]:
        """
        Return CCD read-noise value(s) from the telescope.

        If ``obs_date`` is not provided, return the read noise from all CCD
        definitions as a tuple. If ``obs_date`` is provided, return the read noise
        from the CCD definition valid on that date.
        """
        return self.ccd_value("rdnoise", obs_date)

    def filters(self) -> tuple[str, ...]:
        """
        Return all standard filter names available in the telescope.

        The returned tuple is flattened across all filter sets stored in the
        telescope.
        """

        return tuple(
            standard_name
            for fs in self._filter_sets
            for standard_name in fs.standard_names
        )

    def resolve_filter_alias(
        self,
        alias: str,
    ) -> str:
        """
        Resolve a filter alias to the telescope-standard filter name.

        The alias is searched across all filter sets stored in the telescope.
        The first matching standard name is returned.

        Raises
        ------
        UnknownFilterAliasError
            If the alias cannot be resolved in any filter set.
        """

        for fs in self._filter_sets:
            try:
                return fs.resolve(alias)
            except UnknownFilterAliasError:
                continue

        raise UnknownFilterAliasError(
            f"Unknown filter alias {alias!r} in telescope."
        )
