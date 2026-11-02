"""Public interfaces for refactored image measurements.

Legacy source detection and photometry implementations remain in their modules;
do not import their unrefactored dependencies during FWHM measurement.
"""

from .fwhm import FWHM

__all__ = ["FWHM"]
