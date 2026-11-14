from dataclasses import FrozenInstanceError

import astropy.units as u
from astropy.coordinates import SkyCoord
import pytest

from photozpy.sources import ObservationField
from photozpy.telescope import SARA_RM


@pytest.fixture
def nominal_pointing() -> SkyCoord:
    """Return one scalar pointing coordinate."""
    return SkyCoord(
        ra=150.116321 * u.deg,
        dec=2.205830 * u.deg,
        frame="icrs",
    )


@pytest.fixture
def observation_field(
    nominal_pointing: SkyCoord,
) -> ObservationField:
    """Return a valid science-field definition."""
    return ObservationField(
        field_name="SN2025abc",
        telescope=SARA_RM,
        nominal_pointing=nominal_pointing,
        file_pattern="SN2025abc",
    )


def test_init_stores_field_metadata(
    observation_field: ObservationField,
    nominal_pointing: SkyCoord,
) -> None:
    assert observation_field.field_name == "SN2025abc"
    assert observation_field.telescope is SARA_RM
    assert observation_field.nominal_pointing is nominal_pointing
    assert observation_field.file_pattern == "SN2025abc"


def test_init_accepts_missing_file_pattern(
    nominal_pointing: SkyCoord,
) -> None:
    field = ObservationField(
        field_name="standard_field",
        telescope=SARA_RM,
        nominal_pointing=nominal_pointing,
    )

    assert field.file_pattern is None


def test_init_rejects_non_string_field_name(
    nominal_pointing: SkyCoord,
) -> None:
    with pytest.raises(TypeError, match="field_name must be a string"):
        ObservationField(
            field_name=42,  # type: ignore[arg-type]
            telescope=SARA_RM,
            nominal_pointing=nominal_pointing,
        )


@pytest.mark.parametrize("field_name", ["", "   "])
def test_init_rejects_empty_field_name(
    field_name: str,
    nominal_pointing: SkyCoord,
) -> None:
    with pytest.raises(ValueError, match="field_name.*must not be empty"):
        ObservationField(
            field_name=field_name,
            telescope=SARA_RM,
            nominal_pointing=nominal_pointing,
        )


def test_init_rejects_invalid_telescope(
    nominal_pointing: SkyCoord,
) -> None:
    with pytest.raises(TypeError, match="must be a Telescope object"):
        ObservationField(
            field_name="science_field",
            telescope=object(),  # type: ignore[arg-type]
            nominal_pointing=nominal_pointing,
        )


def test_init_rejects_non_skycoord_pointing() -> None:
    with pytest.raises(TypeError, match="must be a SkyCoord object"):
        ObservationField(
            field_name="science_field",
            telescope=SARA_RM,
            nominal_pointing=(150.0, 2.0),  # type: ignore[arg-type]
        )


def test_init_rejects_non_scalar_pointing() -> None:
    pointing = SkyCoord(
        ra=[150.0, 150.1] * u.deg,
        dec=[2.0, 2.1] * u.deg,
        frame="icrs",
    )

    with pytest.raises(ValueError, match="must contain exactly one coordinate"):
        ObservationField(
            field_name="science_field",
            telescope=SARA_RM,
            nominal_pointing=pointing,
        )


def test_init_rejects_non_string_file_pattern(
    nominal_pointing: SkyCoord,
) -> None:
    with pytest.raises(TypeError, match="must be a string or None"):
        ObservationField(
            field_name="science_field",
            telescope=SARA_RM,
            nominal_pointing=nominal_pointing,
            file_pattern=42,  # type: ignore[arg-type]
        )


@pytest.mark.parametrize("file_pattern", ["", "   "])
def test_init_rejects_empty_file_pattern(
    file_pattern: str,
    nominal_pointing: SkyCoord,
) -> None:
    with pytest.raises(ValueError, match="must not be empty when provided"):
        ObservationField(
            field_name="science_field",
            telescope=SARA_RM,
            nominal_pointing=nominal_pointing,
            file_pattern=file_pattern,
        )


def test_instance_is_frozen(
    observation_field: ObservationField,
) -> None:
    with pytest.raises(FrozenInstanceError):
        observation_field.field_name = "replacement"  # type: ignore[misc]
