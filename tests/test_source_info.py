from dataclasses import FrozenInstanceError

import astropy.units as u
from astropy.coordinates import SkyCoord
import pytest

from photozpy.sources import ObservationField, SourceInfo
from photozpy.telescope import SARA_RM


@pytest.fixture
def observation_field() -> ObservationField:
    """Return one observation field shared by multiple sources."""
    return ObservationField(
        field_name="SN2025abc",
        telescope=SARA_RM,
        nominal_pointing=SkyCoord(
            ra=150.116321 * u.deg,
            dec=2.205830 * u.deg,
            frame="icrs",
        ),
        file_pattern="SN2025abc",
    )


@pytest.fixture
def source_coordinate() -> SkyCoord:
    """Return one scalar source coordinate."""
    return SkyCoord(
        ra=150.120000 * u.deg,
        dec=2.210000 * u.deg,
        frame="icrs",
    )


@pytest.fixture
def source_info(
    observation_field: ObservationField,
    source_coordinate: SkyCoord,
) -> SourceInfo:
    """Return a valid science-source definition."""
    return SourceInfo(
        source_name="SN2025abc",
        source_role="science",
        source_coordinate=source_coordinate,
        observation_field=observation_field,
    )


def test_init_stores_source_metadata(
    source_info: SourceInfo,
    observation_field: ObservationField,
    source_coordinate: SkyCoord,
) -> None:
    assert source_info.source_name == "SN2025abc"
    assert source_info.source_role == "science"
    assert source_info.source_coordinate is source_coordinate
    assert source_info.observation_field is observation_field


def test_field_can_contain_science_and_standard_sources(
    observation_field: ObservationField,
    source_coordinate: SkyCoord,
) -> None:
    science_source = SourceInfo(
        source_name="SN2025abc",
        source_role="science",
        source_coordinate=source_coordinate,
        observation_field=observation_field,
    )
    standard_source = SourceInfo(
        source_name="Gaia DR3 123456789",
        source_role="standard",
        source_coordinate=SkyCoord(
            ra=150.130000 * u.deg,
            dec=2.220000 * u.deg,
            frame="icrs",
        ),
        observation_field=observation_field,
    )

    assert science_source.observation_field is standard_source.observation_field
    assert science_source.source_role == "science"
    assert standard_source.source_role == "standard"


def test_init_rejects_non_string_source_name(
    observation_field: ObservationField,
    source_coordinate: SkyCoord,
) -> None:
    with pytest.raises(TypeError, match="source_name must be a string"):
        SourceInfo(
            source_name=42,  # type: ignore[arg-type]
            source_role="science",
            source_coordinate=source_coordinate,
            observation_field=observation_field,
        )


@pytest.mark.parametrize("source_name", ["", "   "])
def test_init_rejects_empty_source_name(
    source_name: str,
    observation_field: ObservationField,
    source_coordinate: SkyCoord,
) -> None:
    with pytest.raises(ValueError, match="source_name.*must not be empty"):
        SourceInfo(
            source_name=source_name,
            source_role="science",
            source_coordinate=source_coordinate,
            observation_field=observation_field,
        )


def test_init_rejects_non_string_source_role(
    observation_field: ObservationField,
    source_coordinate: SkyCoord,
) -> None:
    with pytest.raises(TypeError, match="source_role.*must be a string"):
        SourceInfo(
            source_name="SN2025abc",
            source_role=42,  # type: ignore[arg-type]
            source_coordinate=source_coordinate,
            observation_field=observation_field,
        )


def test_init_rejects_unknown_source_role(
    observation_field: ObservationField,
    source_coordinate: SkyCoord,
) -> None:
    with pytest.raises(ValueError, match="either 'science' or 'standard'"):
        SourceInfo(
            source_name="SN2025abc",
            source_role="comparison",  # type: ignore[arg-type]
            source_coordinate=source_coordinate,
            observation_field=observation_field,
        )


def test_init_rejects_non_skycoord_coordinate(
    observation_field: ObservationField,
) -> None:
    with pytest.raises(TypeError, match="must be a SkyCoord object"):
        SourceInfo(
            source_name="SN2025abc",
            source_role="science",
            source_coordinate=(150.0, 2.0),  # type: ignore[arg-type]
            observation_field=observation_field,
        )


def test_init_rejects_non_scalar_coordinate(
    observation_field: ObservationField,
) -> None:
    coordinate = SkyCoord(
        ra=[150.0, 150.1] * u.deg,
        dec=[2.0, 2.1] * u.deg,
        frame="icrs",
    )

    with pytest.raises(ValueError, match="must contain exactly one coordinate"):
        SourceInfo(
            source_name="SN2025abc",
            source_role="science",
            source_coordinate=coordinate,
            observation_field=observation_field,
        )


def test_init_rejects_invalid_observation_field(
    source_coordinate: SkyCoord,
) -> None:
    with pytest.raises(TypeError, match="must be an ObservationField object"):
        SourceInfo(
            source_name="SN2025abc",
            source_role="science",
            source_coordinate=source_coordinate,
            observation_field=object(),  # type: ignore[arg-type]
        )


def test_instance_is_frozen(source_info: SourceInfo) -> None:
    with pytest.raises(FrozenInstanceError):
        source_info.source_name = "replacement"  # type: ignore[misc]
