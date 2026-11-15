"""Public interfaces for refactored image measurements.

Legacy source detection and photometry implementations remain in their modules;
do not import their unrefactored dependencies during FWHM measurement.
"""

from .fwhm import FWHM
from .local_background import LocalBackgroundOptions
from .background_diagnostics import BackgroundDiagnostic
from .regions import PhotometryRegions, PhotometryRegionCollection

__all__ = ["BackgroundDiagnostic", "FWHM", "LocalBackgroundOptions", "PhotometryRegions", "PhotometryRegionCollection"]
