"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- ASPCORR check
- Combine extensions within one observation file (the same observation ID)
- Combine images with multiple observation files (across multiple observation IDs)
"""

from ..mimage_collection import mImageFileCollection
from ccdproc import ImageFileCollection

class SwiftCombine():

    def __init__(image_collection):

        """
        image_collection : ccdproc.image_collection.ImageFileCollection or photozpy.mimage_collection.mImageFileCollection
            The image collection(s)
        """
            
        # switch image_collection from ImageFileCollection to mImageFileCollection to make it standard for the pipeline
        if isinstance(image_collection, ImageFileCollection):
            self.mcollection = mImageFileCollection(location = image_collection.location, filenames = image_collection.files)

        elif isinstance(image_collection, mImageFileCollection):
            self.mcollection = image_collection

    @staticmethod
    def check_ASPCORR(fits_file):
        pass

