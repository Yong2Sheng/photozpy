"""Public entry points for the refactored photozpy modules."""

from .calibration import Artifacts, Combine, HeaderCorrection, PlateSolving, Reduction, Registration
from .collections import ImageCollectionGroup
from .interactive_ds9 import ImageReviewer
from .photometry import FWHM
from .sources import FilterMagnitude, ObservationField, SourceInfo, Sources
from .telescope import CCD, FilterSet, Telescope

__version__ = "0.1.0"  # simple, importable version for tooling/tests

__all__ = [
    "CCD",
    "Artifacts",
    "Combine",
    "FilterMagnitude",
    "FilterSet",
    "FWHM",
    "HeaderCorrection",
    "ImageCollectionGroup",
    "ImageReviewer",
    "ObservationField",
    "PlateSolving",
    "Reduction",
    "Registration",
    "SourceInfo",
    "Sources",
    "Telescope",
]
