from __future__ import annotations

import logging
from dataclasses import dataclass, field, fields
from datetime import date, datetime
from types import MappingProxyType
from typing import Any, Mapping

import astropy.units as u
from astropy.units import Quantity
from tabulate import tabulate


logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class CCD:
    """
    Data object for a CCD definition.

    This class describes CCD-specific parameters.
    It should not decide whether the current pipeline
    should use it or not.
    """

    name: str
    gain: Quantity
    rdnoise: Quantity  # expected unit: electron / pix

    valid_from: date | None = None
    valid_to: date | None = None

    # Native unbinned angular pixel scale, e.g. arcsec / pixel.
    native_pixel_scale: Quantity | None = None
    nx: int | None = None
    ny: int | None = None

    saturation: Quantity | None = None
    full_well: Quantity | None = None

    # I usually do not use overcan or trim regions,
    # but I want to preserve their position for
    # future application if needed.
    overscan_region: str | None = None
    trim_region: str | None = None
    data_region: str | None = None

    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """
        This will be executed after the __init__ function.
        """

        # freeze the content of the metadata dict,
        # so they can't be modified.
        object.__setattr__(
            self,
            "metadata",
            MappingProxyType(dict(self.metadata)),
        )
        self.validate()

        if self.valid_from is None and self.valid_to is None:
            logger.warning(
                "CCD %r has no valid date range defined; "
                "it will be treated as always valid.",
                self.name,
            )

    @property
    def shape(self) -> tuple[int, int] | None:
        if self.nx is None or self.ny is None:
            return None
        return (self.ny, self.nx)

    @property
    def has_geometry(self) -> bool:
        return self.nx is not None and self.ny is not None

    def is_valid_for(self, obs_date: date | datetime) -> bool:

        if isinstance(obs_date, datetime):
            obs_date = obs_date.date()

        if self.valid_from is not None and obs_date < self.valid_from:
            return False

        if self.valid_to is not None and obs_date > self.valid_to:
            return False

        return True

    def validate(self) -> None:
        """
        Validate if the inputs are appropriate.
        """

        if not self.name:
            raise ValueError("name must not be empty.")

        if not self.gain.unit.is_equivalent(u.electron / u.adu):
            raise ValueError(
                "gain unit must be equivalent to electron/adu."
            )

        # Read noise is stored as per-pixel RMS readout noise.
        if not self.rdnoise.unit.is_equivalent(u.electron / u.pix):
            raise ValueError(
                "rdnoise must be the per-pixel RMS readout noise, "
                "with units equivalent to electron/pixel."
            )

        if self.gain.value <= 0:
            raise ValueError("gain must be positive.")

        if self.rdnoise.value < 0:
            raise ValueError("rdnoise must be non-negative.")

        if self.valid_from is not None and not isinstance(self.valid_from, date):
            raise TypeError("valid_from must be a date or None.")

        if self.valid_to is not None and not isinstance(self.valid_to, date):
            raise TypeError("valid_to must be a date or None.")

        if self.valid_from is not None and self.valid_to is not None:
            if self.valid_from > self.valid_to:
                raise ValueError(
                    "valid_from must be earlier than or equal to valid_to."
                )

        if self.native_pixel_scale is not None:
            if not self.native_pixel_scale.unit.is_equivalent(u.arcsec / u.pix):
                raise ValueError(
                    "native_pixel_scale must have angular-per-pixel units, e.g. arcsec/pix."
                )
            if self.native_pixel_scale.value <= 0:
                raise ValueError("native_pixel_scale must be positive.")

        if self.nx is not None:
            if not isinstance(self.nx, int):
                raise TypeError("nx must be an int or None.")
            if self.nx <= 0:
                raise ValueError("nx must be positive.")

        if self.ny is not None:
            if not isinstance(self.ny, int):
                raise TypeError("ny must be an int or None.")
            if self.ny <= 0:
                raise ValueError("ny must be positive.")

        if self.saturation is not None:
            if not (
                self.saturation.unit.is_equivalent(u.adu)
                or self.saturation.unit.is_equivalent(u.electron)
            ):
                raise ValueError(
                    "saturation must have ADU or electron units."
                )
            if self.saturation.value <= 0:
                raise ValueError("saturation must be positive.")

        if self.full_well is not None:
            if not self.full_well.unit.is_equivalent(u.electron):
                raise ValueError(
                    "full_well must have electron units."
                )
            if self.full_well.value <= 0:
                raise ValueError("full_well must be positive.")

    def to_dict(self) -> dict[str, object]:
        return {
            f.name: getattr(self, f.name)
            for f in fields(self)
        }

    def to_table(self) -> str:
        rows: list[list[object]] = [
            [key, value]
            for key, value in self.to_dict().items()
        ]
        return tabulate(rows, headers=["field", "value"], tablefmt="fancy_grid")

    def print_summary(self) -> None:
        print(self.to_table())
