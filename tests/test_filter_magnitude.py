import astropy.units as u
from astropy.table import QTable
import numpy as np
import pytest

from photozpy.sources.filter_magnitude import (
    FilterMagnitude,
    UPPER_LIMIT_ERROR,
)
from photozpy.telescope import SARA_RM


@pytest.fixture
def magnitude_info() -> FilterMagnitude:
    """Create filter magnitudes with aliases, missing filters, and an upper limit."""
    return FilterMagnitude(
        telescope=SARA_RM,
        filter_names=(
            "SDSS r",
            "Sloan g'2",
            "SDSS i",
        ),
        magnitudes=u.Quantity(
            [19.8421, 20.1342, 21.3500],
            u.mag,
        ),
        errors=u.Quantity(
            [
                0.0285,
                0.0321,
                UPPER_LIMIT_ERROR.to_value(u.mag),
            ],
            u.mag,
        ),
        significance=(35.20, 31.40, 2.10),
        magnitude_system="AB",
    )


def test_init_resolves_aliases_and_uses_telescope_filter_order(
    magnitude_info: FilterMagnitude,
) -> None:
    assert magnitude_info.filter_names == SARA_RM.filters()
    assert magnitude_info.qtable.colnames == [
        "filter",
        "magnitude",
        "error",
        "significance",
    ]

    np.testing.assert_allclose(
        magnitude_info.magnitudes.to_value(u.mag),
        [np.nan, 20.1342, 19.8421, 21.3500, np.nan],
        equal_nan=True,
    )
    np.testing.assert_allclose(
        magnitude_info.errors.to_value(u.mag),
        [np.nan, 0.0321, 0.0285, -99.0, np.nan],
        equal_nan=True,
    )
    np.testing.assert_allclose(
        magnitude_info.significance,
        [np.nan, 31.40, 35.20, 2.10, np.nan],
        equal_nan=True,
    )


def test_init_interprets_unitless_values_as_magnitudes() -> None:
    result = FilterMagnitude(
        telescope=SARA_RM,
        filter_names=("SDSS g",),
        magnitudes=(20.5,),
        errors=(0.1,),
        significance=(12.0,),
    )

    g_index = result.filter_names.index("SDSS_g'")

    assert result.magnitudes.unit == u.mag
    assert result.errors.unit == u.mag
    assert result.magnitudes[g_index] == 20.5 * u.mag
    assert result.errors[g_index] == 0.1 * u.mag
    assert result.magnitude_system == "AB"


def test_init_leaves_optional_measurements_as_nan() -> None:
    result = FilterMagnitude(
        telescope=SARA_RM,
        filter_names=("SDSS g",),
        magnitudes=(20.5,),
        errors=None,
        significance=None,
    )

    assert np.isnan(result.errors.to_value(u.mag)).all()
    assert np.isnan(result.significance).all()


@pytest.mark.parametrize(
    "magnitude_system",
    ["instrumental", "AB", "Vega"],
)
def test_init_accepts_supported_magnitude_systems(
    magnitude_system: str,
) -> None:
    result = FilterMagnitude(
        telescope=SARA_RM,
        filter_names=("SDSS g",),
        magnitudes=(20.5,),
        errors=None,
        significance=None,
        magnitude_system=magnitude_system,  # type: ignore[arg-type]
    )

    assert result.magnitude_system == magnitude_system


def test_init_rejects_invalid_telescope() -> None:
    with pytest.raises(
        TypeError,
        match="must be a Telescope object",
    ):
        FilterMagnitude(
            telescope=object(),  # type: ignore[arg-type]
            filter_names=("SDSS g",),
            magnitudes=(20.5,),
            errors=None,
            significance=None,
        )


def test_init_rejects_unknown_magnitude_system() -> None:
    with pytest.raises(
        ValueError,
        match="must be 'instrumental', 'AB', or 'Vega'",
    ):
        FilterMagnitude(
            telescope=SARA_RM,
            filter_names=("SDSS g",),
            magnitudes=(20.5,),
            errors=None,
            significance=None,
            magnitude_system="ST",  # type: ignore[arg-type]
        )


def test_init_rejects_single_filter_name_string() -> None:
    with pytest.raises(
        TypeError,
        match="not a single string",
    ):
        FilterMagnitude(
            telescope=SARA_RM,
            filter_names="SDSS g",  # type: ignore[arg-type]
            magnitudes=(20.5,),
            errors=None,
            significance=None,
        )


def test_init_rejects_different_filter_and_magnitude_lengths() -> None:
    with pytest.raises(
        ValueError,
        match="must have the same length",
    ):
        FilterMagnitude(
            telescope=SARA_RM,
            filter_names=("SDSS g", "SDSS r"),
            magnitudes=(20.5,),
            errors=None,
            significance=None,
        )


def test_init_rejects_duplicate_canonical_filters() -> None:
    with pytest.raises(
        ValueError,
        match="resolve to the same standard filter name",
    ):
        FilterMagnitude(
            telescope=SARA_RM,
            filter_names=("SDSS_g'", "SDSS g"),
            magnitudes=(20.5, 20.6),
            errors=None,
            significance=None,
        )


def test_init_rejects_unknown_filter_alias() -> None:
    with pytest.raises(
        ValueError,
        match="Unknown filter alias",
    ):
        FilterMagnitude(
            telescope=SARA_RM,
            filter_names=("unknown",),
            magnitudes=(20.5,),
            errors=None,
            significance=None,
        )


def test_properties_read_current_qtable_values(
    magnitude_info: FilterMagnitude,
) -> None:
    g_index = magnitude_info.filter_names.index("SDSS_g'")

    magnitude_info.qtable["magnitude"][g_index] = 18.5 * u.mag
    magnitude_info.qtable["error"][g_index] = 0.02 * u.mag
    magnitude_info.qtable["significance"][g_index] = 50.0

    assert magnitude_info.magnitudes[g_index] == 18.5 * u.mag
    assert magnitude_info.errors[g_index] == 0.02 * u.mag
    assert magnitude_info.significance[g_index] == 50.0


def test_to_table_formats_missing_values_and_upper_limit(
    magnitude_info: FilterMagnitude,
) -> None:
    table = magnitude_info.to_table()

    assert table.startswith("Magnitude system: AB\n")
    assert "Filter" in table
    assert "SDSS_g'" in table
    assert "20.1342" in table
    assert "UL (-99)" in table
    assert "--" in table


def test_from_qtable_preserves_canonical_values(
    magnitude_info: FilterMagnitude,
) -> None:
    restored = FilterMagnitude.from_qtable(
        telescope=SARA_RM,
        table=magnitude_info.qtable,
        magnitude_system="AB",
    )

    assert restored.qtable is not magnitude_info.qtable
    assert restored.filter_names == magnitude_info.filter_names
    np.testing.assert_allclose(
        restored.magnitudes.to_value(u.mag),
        magnitude_info.magnitudes.to_value(u.mag),
        equal_nan=True,
    )
    np.testing.assert_allclose(
        restored.errors.to_value(u.mag),
        magnitude_info.errors.to_value(u.mag),
        equal_nan=True,
    )
    np.testing.assert_allclose(
        restored.significance,
        magnitude_info.significance,
        equal_nan=True,
    )


def test_from_qtable_requires_all_columns(
    magnitude_info: FilterMagnitude,
) -> None:
    table = QTable(magnitude_info.qtable)
    table.remove_column("error")

    with pytest.raises(KeyError, match="error"):
        FilterMagnitude.from_qtable(
            telescope=SARA_RM,
            table=table,
            magnitude_system="AB",
        )
