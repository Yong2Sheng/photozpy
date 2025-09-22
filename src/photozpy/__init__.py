# Keep package very lean: no top-level re-exports, no star imports.
# Users should import from submodules explicitly, e.g.:
#   from photozpy.photometry import Photometry
#   from photozpy.telescope import Telescope
# This keeps the public API explicit and avoids namespace pollution.

__version__ = "0.1.0"  # simple, importable version for tooling/tests

# Optionally, define __all__ to keep `from photozpy import *` empty on purpose.
__all__ = []  # explicit: top-level exports are intentionally empty


# from .swiftz import *
# from .telescope import Telescope
# from .collection_manager import CollectionManager
# from .calibration import HeaderCorrection, HeaderManipulation, Combine, Reduction, Registration, PlateSolving, Artifacts
# from .photometry import SourceDetection, Photometry, SwiftPhotometry
# from .convenience_functions import *
# from .mimage_collection import *
# from .sources import *
