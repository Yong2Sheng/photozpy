from .swiftz import SwiftDownload, SwiftCombine, CCD_Regions, PhotozRegions
from .telescope import Telescope, CCD, FilterSet
from .collection_manager import CollectionManager
from .collections import ImageCollectionGroup
from .calibration import HeaderCorrection, HeaderManipulation, Combine, Reduction, Registration, PlateSolving, Artifacts
from .photometry import SourceDetection, Photometry, SwiftPhotometry
from .mimage_collection import mImageFileCollection
from .sources import MagnitudeInfo, BaseSourceInfo, Sources
from .standard_stars import StandardStarCatalog
from .interactive_ds9 import ImageReviewer

__version__ = "0.1.0"  # simple, importable version for tooling/tests

__all__ = [
    "SwiftDownload",
    "SwiftCombine",
    "CCD_Regions",
    "Telescope",
    "CollectionManager",
    "HeaderCorrection",
    "HeaderManipulation",
    "Combine",
    "Reduction",
    "Registration",
    "PlateSolving",
    "Artifacts",
    "SourceDetection",
    "Photometry",
    "SwiftPhotometry",
    "mImageFileCollection",
    "MagnitudeInfo",
    "BaseSourceInfo",
    "Sources",
    "PhotozRegions",
    "CCD",
    "StandardStarCatalog",
    "FilterSet",
    "ImageCollectionGroup",
    "ImageReviewer",
]
