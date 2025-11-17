"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- register a series of images
- only Light type is accecpted
"""

from ccdproc import ImageFileCollection
from .headers import HeaderCorrection, HeaderManipulation
from ..collection_manager import CollectionManager
import numpy as np
from pathlib import Path
from astropy.io import fits
import copy
from astroalign import register
from ccdproc import cosmicray_lacosmic as lacosmic

class Registration():

    def __init__(self, image_collection, telescope):

        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(image_collection, rescan = True)
        self._telescope = telescope


    @staticmethod
    def _Register_images(image_collection, filter, min_area=5, detection_sigma=5):

        """
        Register all the images with same object in the collection.

        Parameters
        ----------
        image_collection

        Returns
        -------
        new_image_collection
        """

        # get the image type and make sure the collection only has Light type images
        image_type = HeaderManipulation.get_header_values(image_collection, header = "imtype")

        if len(image_type) != 1:
            print(f"Only Light image type is supported. You have more than one image type: {image_type}!")
        
        elif image_type[0] != "Light":
            raise TypeError(f"Only Light image type is supported for registration. You have {image_type}")

        # get all the object names, make sure all the object names are the same
        object_name = HeaderManipulation.get_header_values(image_collection, header = "object")
        if len(object_name) != 1:
            raise ValueError("You have more than one object in the image collection!")

        print(f"Aligning {object_name} ......")

        # get the image collection that only contains the reference image
        reference_collection = CollectionManager.filter_collection(image_collection, **{"FILTER": filter})
        reference_image_path = Path(reference_collection.files_filtered(include_path = True)[0])
        refernece_image_name = reference_image_path.name
        
        reference_ext0 = fits.getdata(reference_image_path, ext=0)
        reference_data = reference_ext0.byteswap().newbyteorder()
        #reference_data, mask = lacosmic(reference_data)

        # get the image collection to be registered by removing the one used to register
        register_collection = CollectionManager.delete_images(image_collection, refernece_image_name)
        register_image_paths = register_collection.files_filtered(include_path = True)

        # start image registeration, the files will be over written
        for i in register_image_paths:
            target_file_name = Path(i).name
            print(f"Aligning {target_file_name} using {refernece_image_name}")

            target_ext0 = fits.getdata(i, ext = 0)
            header_target = fits.getheader(i)
            target_data = target_ext0.byteswap().newbyteorder()
            #target_data, mask = lacosmic(target_data)

            # reference_ext0 = fits.getdata(reference_image_path, ext=0)
            # reference_data = reference_ext0.byteswap().newbyteorder()
            # reference_data, mask = lacosmic(reference_data)
            
            img_aligned, footprint = register(target_data, reference_data, min_area=min_area, detection_sigma=detection_sigma)

            header_target["ALIGN"] = "Registered"
            hdu = fits.PrimaryHDU(img_aligned, header_target)
            hdu.writeto(i, overwrite=True)

        # write the alignment keyword for the reference
        with fits.open(reference_image_path, mode = "update") as hdul:
            hdul[0].header["ALIGN"] = "Reference"
            obj_name = hdul[0].header["OBJECT"]
            hdul.flush()

        print(f"Alignment of {obj_name} finished!")
        print("----------------------------------------------------------\n")

        return

    def Register_images(self, filter, min_area=5, detection_sigma=5):

        """
        Register all the images with different objects in the collection.

        Paremeters
        ----------
        filter: string; the filter that the reference image should have. The reason for this is to make sure
                        that the reference image relatively sensitive to have more sources for registration.

        Returns
        -------
        image_collection
        """

        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(self._image_collection)

        # get all the Light type images
        images_collection_to_register = CollectionManager.filter_collection(self._image_collection, **{"IMTYPE": "Light"})

        # get the object names, repeated names are removed.
        object_names = HeaderManipulation.get_header_values(images_collection_to_register, header = "object")

        for object_name in object_names:

            # get the image collection with same object name
            to_register = CollectionManager.filter_collection(images_collection_to_register, **{"OBJECT": object_name})

            # register the images of same object
            Registration._Register_images(to_register, filter = filter, min_area = min_area, detection_sigma=detection_sigma)

        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(self._image_collection)
        

        
        

        

        

        

        