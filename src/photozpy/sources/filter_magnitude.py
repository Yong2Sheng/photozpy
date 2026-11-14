"""Store and display filter-aligned magnitude measurements."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Literal, Self, TypeAlias

import astropy.units as u
from astropy.table import QTable
import numpy as np
from tabulate import tabulate

from ..telescope import Telescope


MagnitudeSystem: TypeAlias = Literal["instrumental", "AB", "Vega"]
VALID_SYSTEM = ("instrumental", "AB", "Vega")
# Sentinel required by the LePhare upper-limit convention.
# This is not a physical magnitude uncertainty.
UPPER_LIMIT_ERROR = -99.0 * u.mag


class FilterMagnitude:
    """Represent one magnitude measurement per telescope filter.

    Input filter names are resolved through the associated
    :class:`~photozpy.telescope.Telescope` and stored in the telescope's
    canonical filter order. The internal :class:`astropy.table.QTable` contains
    one row for every telescope filter. Filters without supplied measurements
    remain present with ``NaN`` values.

    Parameters
    ----------
    telescope : Telescope
        Telescope configuration that defines the allowed filter names, aliases,
        and canonical filter order.
    filter_names : Sequence[str]
        Filter names or aliases corresponding to the supplied magnitudes.
    magnitudes : Sequence[float] or astropy.units.Quantity
        One magnitude for each input filter. Unitless values are interpreted as
        magnitudes and stored with unit ``u.mag``.
    magnitude_system : {"instrumental", "AB", "Vega"}, optional
        Photometric system of the supplied magnitudes. This value records
        metadata only; this class does not convert between magnitude systems.
        Default is ``"AB"``.
    errors : Sequence[float], astropy.units.Quantity, or None, optional
        Magnitude uncertainties corresponding to the input filters. Unitless
        values are interpreted as magnitudes. If ``None``, the error column
        remains ``NaN``. The value ``-99 * u.mag`` is reserved for the LePhare
        upper-limit convention and is not a physical uncertainty.
    significance : Sequence[float] or None
        Detection significance corresponding to the input filters. If ``None``,
        the significance column remains ``NaN``.

    Attributes
    ----------
    telescope : Telescope
        Telescope associated with the measurements.
    magnitude_system : MagnitudeSystem
        Magnitude system recorded for the measurements.
    qtable : astropy.table.QTable
        Canonical table containing ``filter``, ``magnitude``, ``error``, and
        ``significance`` columns.

    Notes
    -----
    The QTable is the authoritative data store. The ``filter_names``,
    ``magnitudes``, ``errors``, and ``significance`` properties read directly
    from it.

    This class stores at most one value per telescope filter. Multiple
    per-exposure measurements for variability or light-curve analysis should be
    represented by a separate time-series data model.
    """

    def __init__(
        self,
        *,
        telescope: Telescope,
        filter_names: Sequence[str],
        magnitudes: Sequence[float] | u.Quantity,
        magnitude_system: MagnitudeSystem = "AB",
        errors: Sequence[float] | u.Quantity | None = None,
        significance: Sequence[float] | None,
    ) -> None:

        self.telescope = telescope
        self.magnitude_system = magnitude_system

        # check if input are valid
        if not isinstance(self.telescope, Telescope):
            raise TypeError("'telescope' must be a Telescope object.")

        if self.magnitude_system not in VALID_SYSTEM:
            raise ValueError(
                "'magnitude_system' must be 'instrumental', 'AB', or 'Vega'."
            )

        if isinstance(filter_names, str):
            raise TypeError(
                "'filter_names' must be an iterable of strings, "
                "not a single string."
            )

        input_filter_names = tuple(filter_names)

        if isinstance(magnitudes, u.Quantity):
            magnitude_values = np.atleast_1d(magnitudes)
        else:
            magnitude_values = (
                np.asarray(tuple(magnitudes), dtype=float) * u.mag
            )

        # check if filter_names and magnitudes have the same length
        if len(input_filter_names) != len(magnitude_values):
            raise ValueError(
                "'filter_names' and 'magnitudes' must have the same length."
            )

        if errors is None:
            error_values = None
        elif isinstance(errors, u.Quantity):
            error_values = np.atleast_1d(errors)
        else:
            error_values = (
                np.asarray(tuple(errors), dtype=float) * u.mag
            )

        if significance is None:
            significance_values = None
        else:
            significance_values = np.atleast_1d(significance)

        # standardize filter names
        standard_names = tuple(
            self.telescope.resolve_filter_alias(alias)
            for alias in input_filter_names
        )
        if len(standard_names) != len(set(standard_names)):
            raise ValueError(
                "Multiple filter names resolve to the same standard "
                "filter name! Please check 'filter_names'."
            )

        # generate a null QTable first based in the
        # filters in the Telescope.
        # Note that this filters is from the Telescope,
        # not from the input filter_names.
        filters: tuple[str, ...] = self.telescope.filters()
        n_filters = len(filters)
        self.qtable = QTable(
            {
                "filter": filters,
                "magnitude": np.full(n_filters, np.nan) * u.mag,
                "error": np.full(n_filters, np.nan) * u.mag,
                "significance": np.full(n_filters, np.nan),
            }
        )

        # build canonical row lookup dict for QTable
        filter_indices: dict = {
            filter_name: index
            for index, filter_name in enumerate(filters)
        }

        # fill up the QTable
        for index, standard_name in enumerate(standard_names):

            row_index = filter_indices[standard_name]

            self.qtable["magnitude"][row_index] = magnitude_values[index]

            if error_values is not None:
                self.qtable["error"][row_index] = error_values[index]

            if significance_values is not None:
                self.qtable["significance"][row_index] = significance_values[index]

    @property
    def filter_names(self) -> tuple[str, ...]:
        """Return all canonical filter names in telescope-defined order."""
        return tuple(self.qtable["filter"])

    @property
    def magnitudes(self) -> u.Quantity:
        """Return the magnitude column aligned to ``filter_names``."""
        return self.qtable["magnitude"]

    @property
    def errors(self) -> u.Quantity:
        """Return the magnitude-error column aligned to ``filter_names``."""
        return self.qtable["error"]

    @property
    def significance(self):
        """Return the detection-significance column aligned to ``filter_names``."""
        return self.qtable["significance"]

    def to_table(self) -> str:
        """Format the complete magnitude information as a readable table.

        Missing values are displayed as ``--``. Errors equal to
        :data:`UPPER_LIMIT_ERROR` are displayed as ``UL (-99)`` without changing
        the underlying QTable values.

        Returns
        -------
        str
            Magnitude-system label followed by a ``fancy_grid`` table containing
            the filter, magnitude, error, and significance columns.
        """
        rows = []

        for filter_name, magnitude, error, significance in zip(
            self.filter_names,
            self.magnitudes.to_value(u.mag),
            self.errors.to_value(u.mag),
            self.significance,
        ):
            magnitude_text = (
                "--" if np.isnan(magnitude)
                else f"{magnitude:.4f}"
            )

            if np.isnan(error):
                error_text = "--"
            elif error == UPPER_LIMIT_ERROR.to_value(u.mag):
                error_text = "UL (-99)"
            else:
                error_text = f"{error:.4f}"

            significance_text = (
                "--" if np.isnan(significance)
                else f"{significance:.2f}"
            )

            rows.append(
                [
                    filter_name,
                    magnitude_text,
                    error_text,
                    significance_text,
                ]
            )

        table = tabulate(
            rows,
            headers=[
                "Filter",
                "Magnitude (mag)",
                "Error (mag)",
                "Significance",
            ],
            tablefmt="fancy_grid",
            disable_numparse=True,
        )

        return f"Magnitude system: {self.magnitude_system}\n{table}"

    @classmethod
    def from_qtable(
        cls,
        *,
        telescope: Telescope,
        table: QTable,
        magnitude_system: MagnitudeSystem,
    ) -> Self:
        """Construct filter-aligned magnitudes from a QTable.

        The input rows are passed through the normal constructor, so filter
        aliases are resolved and the resulting table follows the canonical
        filter order defined by ``telescope``.

        Parameters
        ----------
        telescope : Telescope
            Telescope configuration used to validate and order the filters.
        table : astropy.table.QTable
            Input table containing ``filter``, ``magnitude``, ``error``, and
            ``significance`` columns.
        magnitude_system : {"instrumental", "AB", "Vega"}
            Photometric system associated with the table values.

        Returns
        -------
        FilterMagnitude
            A new instance with a canonical, telescope-aligned QTable.

        Raises
        ------
        KeyError
            If any required input column is missing.
        """

        filter_names = table["filter"]
        magnitudes = table["magnitude"]
        errors = table["error"]
        significance = table["significance"]

        return cls(
            telescope=telescope,
            filter_names=filter_names,
            magnitudes=magnitudes,
            errors=errors,
            significance=significance,
            magnitude_system=magnitude_system,
        )
