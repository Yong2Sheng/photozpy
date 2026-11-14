"""Define immutable metadata for telescope observation fields."""

from __future__ import annotations

from dataclasses import dataclass

from astropy.coordinates import SkyCoord

from ..telescope import Telescope


@dataclass(frozen=True, slots=True, kw_only=True, eq=False)
class ObservationField:
    """Describe one telescope observation field.

    An observation field stores field-level metadata known before image
    processing. It does not contain science targets, standard stars, measured
    source coordinates, or photometric calibration results.

    Parameters
    ----------
    field_name : str
        Non-empty identifier for the observation field. This name may later be
        used to identify the field in FITS metadata or processing outputs.
    telescope : Telescope
        Telescope configuration used to observe the field.
    nominal_pointing : astropy.coordinates.SkyCoord
        One scalar coordinate representing the planned telescope pointing. It
        is not necessarily a science-target coordinate or the WCS-derived center
        of an individual image.
    file_pattern : str or None, optional
        Optional non-empty filename-matching pattern associated with the field.
        Default is ``None``.

    Notes
    -----
    ``frozen=True`` prevents reassignment of the dataclass attributes, but the
    supplied ``Telescope`` and ``SkyCoord`` objects are stored by reference and
    are not defensively copied.
    """

    field_name: str
    telescope: Telescope
    nominal_pointing: SkyCoord
    file_pattern: str | None = None

    def __post_init__(self) -> None:

        # a lot of input checks
        if not isinstance(self.field_name, str):
            raise TypeError(
                "field_name must be a string."
            )

        if not self.field_name.strip():
            raise ValueError(
                "'field_name' must not be empty."
            )

        if not isinstance(self.telescope, Telescope):
            raise TypeError(
                "'telescope' must be a Telescope object."
            )

        if not isinstance(self.nominal_pointing, SkyCoord):
            raise TypeError(
                "'nominal_pointing' must be a SkyCoord object."
            )

        if not self.nominal_pointing.isscalar:
            raise ValueError(
                "'nominal_pointing' must contain exactly one coordinate."
            )

        if self.file_pattern is not None:
            if not isinstance(self.file_pattern, str):
                raise TypeError(
                    "'file_pattern' must be a string or None."
                )

            if not self.file_pattern.strip():
                raise ValueError(
                    "'file_pattern' must not be empty when provided."
                )
