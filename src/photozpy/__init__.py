from .calibration import (
    Artifacts,
    Combine,
    HeaderCorrection,
    HeaderManipulation,
    PlateSolving,
    Reduction,
    Registration,
)
from .collection_manager import CollectionManager
from .collections import ImageCollectionGroup
from .interactive_ds9 import ImageReviewer
from .mimage_collection import mImageFileCollection
from .photometry import Photometry, SourceDetection, SwiftPhotometry
from .sources import BaseSourceInfo, MagnitudeInfo, Sources
from .standard_stars import StandardStarCatalog
from .swiftz import CCD_Regions, PhotozRegions, SwiftCombine, SwiftDownload
from .telescope import CCD, FilterSet, Telescope

__version__ = "0.1.0"  # simple, importable version for tooling/tests

__all__ = [
    "CCD",
    "Artifacts",
    "BaseSourceInfo",
    "CCD_Regions",
    "CollectionManager",
    "Combine",
    "FilterSet",
    "HeaderCorrection",
    "HeaderManipulation",
    "ImageCollectionGroup",
    "ImageReviewer",
    "MagnitudeInfo",
    "Photometry",
    "PhotozRegions",
    "PlateSolving",
    "Reduction",
    "Registration",
    "SourceDetection",
    "Sources",
    "StandardStarCatalog",
    "SwiftCombine",
    "SwiftDownload",
    "SwiftPhotometry",
    "Telescope",
    "mImageFileCollection",
]
