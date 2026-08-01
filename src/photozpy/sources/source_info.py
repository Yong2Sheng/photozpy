"""Define immutable metadata for individual astronomical sources."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, TypeAlias

from astropy.coordinates import SkyCoord

from ..telescope import Telescope
from .observation_field import ObservationField

SourceRole: TypeAlias = Literal["science", "standard"]
VALID_SOURCE_ROLES = ("science", "standard")


@dataclass(frozen=True, slots=True, kw_only=True, eq=False)
class SourceInfo:
    """Describe one identified source within an observation field.

    Each instance represents exactly one source with one scalar sky coordinate.
    Source role is independent of the observation field, so the same field may
    contain both science and standard sources. Collections of sources and
    photometric measurements are managed separately.

    Parameters
    ----------
    source_name : str
        Non-empty source identifier, such as a target name or stable catalog
        identifier.
    source_role : {"science", "standard"}
        Role of the source in the analysis workflow.
    source_coordinate : astropy.coordinates.SkyCoord
        One scalar sky coordinate identifying the source position.
    observation_field : ObservationField
        Observation field containing the source.

    Notes
    -----
    Automatically discovered standard stars should be represented by
    ``SourceInfo`` only after the catalog provides a name and coordinate; an
    incomplete placeholder source is not required.

    ``frozen=True`` prevents reassignment of the dataclass attributes, but the
    supplied coordinate and observation field are stored by reference and are
    not defensively copied.
    """

    source_name: str
    source_role: SourceRole
    source_coordinate: SkyCoord
    observation_field: ObservationField

    def __post_init__(self) -> None:

        # a lot of input checks
        if not isinstance(self.source_name, str):
            raise TypeError(
                "source_name must be a string."
            )

        if not self.source_name.strip():
            raise ValueError(
                "'source_name' must not be empty."
            )

        if not isinstance(self.source_role, str):
            raise TypeError(
                "'source_role' must be a string."
            )

        if self.source_role not in VALID_SOURCE_ROLES:
            raise ValueError(
                "'source_role' must be either 'science' or 'standard'."
            )

        if not isinstance(self.source_coordinate, SkyCoord):
            raise TypeError(
                "'source_coordinate' must be a SkyCoord object."
            )

        if not self.source_coordinate.isscalar:
            raise ValueError(
                "'source_coordinate' must contain exactly one coordinate."
            )

        if not isinstance(self.observation_field, ObservationField):
            raise TypeError(
                "'observation_field' must be an ObservationField object."
            )
