"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- Combine a series of registered images.
- Check image types before combine.
"""

from ccdproc import ImageFileCollection
from ccdproc import combine as ccdpro_combine
from .headers import HeaderCorrection
from .headers import HeaderManipulation
from ..collection_manager import CollectionManager
from astropy.stats import mad_std, sigma_clipped_stats, median_absolute_deviation
import numpy as np
from pathlib import Path

class Combine():
    
    """
    Combine different types of images.
    """
    
    def __init__(self, image_collection, telescope):
        """
        Init the instance.
        
        Parameters
        ----------
        image_collection: ccdproc.ImageFileCollection; the collection of images
        telescope:
        """

        self._image_collection = CollectionManager.refresh_collection(image_collection)
        self._telescope = telescope
        
        
    @staticmethod
    def inv_median(a):
        return 1/np.median(a)

    @staticmethod
    def check_imtype(image_collection, imtype):
        """
        Check if the imtype of the image_collection is the same as the one we input.

        Parameters
        ----------
        image_collection: ccdproc.ImageFileCollection; the collection of images
        imtype: str; the value for IMTYPE header.

        Returns
        -------
        final: boolean; True or False
        """
        
        headers_values = {"IMTYPE": imtype}

        imtype_value = HeaderManipulation.check_headers(image_collection, **headers_values)
        
        return imtype_value

    def combine_bias_or_dark(self, image_type, method = "sigma clip", 
                             sigma_clip_low=5, sigma_clip_high=5, save_location = "", quite = False):
        """
        Combine bias or dark images

        Parameters
        ----------
        image_type; str; Bias or Dark
        method:
        sigma_clip_low
        sigma_clip_high

        Returns
        -------
        """
        
        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(self._image_collection, rescan = True)

        if save_location == "":
            save_location = self._image_collection.location

        # Get the collection to combine and the object names
        collection_to_combine = CollectionManager.filter_collection(self._image_collection, **{"IMTYPE": [image_type]})

        # Get the object names (Bias or Dark)
        object_names = HeaderManipulation.get_header_values(collection_to_combine, header = "object")  # It seems only the lower cases work for the header values
        object_names = [*set(object_names)]  # remove duplicate elements
        # note that object_names inludes all the objects! (targets, bias ,dark and flat)

        for object_name in object_names:  # note that bias and dark don't have filters
            # this is the list of the absolute path to the fits images to be combine
            image_list_to_combine = list(collection_to_combine.files_filtered(object = object_name, include_path = True))
            # This is the list of the file names
            image_filenames = list(collection_to_combine.files_filtered(object = object_name, include_path = False))
            image_number = len(image_filenames)

            if image_number == 0:
                print(f"WARNING: No {object_name} found! Skipping......")

            elif image_number == 1:
                print(f"Only one {object_name} image! Just copying and renaming the file.")
                
                for hdu in collection_temp.hdus(save_with_name = "_master", save_location = save_location, overwrite = True):
                    hdu.header["IMTYPE"] = f"Master {image_type}"
                    hdu.header["NCOMBINE"] = int(image_number)
                
                fname_stem = Path(image_list_to_combine[0]).stem  # get the filename stem: crab_sdss_g.fits --get--> crab_sdss_g
                fname_suffix = Path(image_list_to_combine[0]).suffix
                name_path = Path(save_location) / f"{fname_stem}_master{fname_suffix}"  # the path to the saved file
                new_name_path = name_path.with_name(f"Master_{object_name}{fname_suffix}")  # the new path
                name_path.rename(new_name_path)  # rename the file, notice that you have to use the full path to rename, not just the file name

            elif image_number > 1:
                if quite == False:
                    print(f"Using {image_filenames} to combine the master {object_name} image!")
                if method == "sigma clip":
                    combined = ccdpro_combine(image_list_to_combine,  
                                              method = "average",
                                              sigma_clip=True, sigma_clip_low_thresh=sigma_clip_low, sigma_clip_high_thresh=sigma_clip_high,
                                              sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)
                elif method == "minmax clip":
                    raise ValueError("This function isn't ready yet!")

                else:
                    raise ValueError("Unsupported combining method! What supported now is only sigma clip!")
                
                combined.meta["IMTYPE"] = f"Master {image_type}"
                combined.meta["NCOMBINE"] = int(image_number)
                combined.write(Path(save_location) / f'Master_{object_name}.fits', overwrite=True)
                print(f"{image_type} combined!")
                print("----------------------------------------\n")


        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(self._image_collection, rescan = True)

        return


    def combine_light_or_flat(self, image_type, method = "sigma clip", scale_function = "auto", 
                              sigma_clip_low=5, sigma_clip_high=5, save_location = "", quite = False):
        """
        Combine light or flat images

        Parameters
        ----------
        image_type; str; Flat or Light

        Returns
        -------
        """

        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(self._image_collection, rescan = True)

        if save_location == "":
            save_location = self._image_collection.location

        # flat or light IMTYPE has filters
        filters = self._telescope.filters
            
        # first select image types (Flat or Light)
        collection_to_combine = CollectionManager.filter_collection(self._image_collection, **{"IMTYPE": [image_type]})

        # Get the object names (Flat or various target names)
        object_names = HeaderManipulation.get_header_values(collection_to_combine, header = "object")  # It seems only the lower cases work
        object_names = [*set(object_names)]  # remove duplicate elements
    
        for object_name in object_names:
            for filter in filters:
                # get the image collection with specific filter and object name
                collection_temp = CollectionManager.filter_collection(collection_to_combine, **{"OBJECT": [object_name], "FILTER": [filter]})
                image_list_to_combine = list(collection_temp.files_filtered(include_path = True)) # this is the list of the absolute path to the fits images to be combined
                image_filenames = list(collection_temp.files_filtered()) # this is the list of file names of the images to be combined
                image_number = len(image_filenames)
                
                if image_number == 0:
                    print(f"WARNING: No {object_name} in {filter} filter found! Skipping......")
                    
                elif image_number == 1:
                    print(f"Only one {object_name} image in {filter}! Just copying and renaming the file.")

                    # construct the headers and values to write
                    header_dict = {"IMTYPE": f"Master {image_type}", 
                                   "NCOMBINE": int(image_number)}
                    if image_type == "Light":
                        header_dict["RDNOISE"] = self._telescope.ccd_rdnoise.value
                        header_dict["GAIN"] = self._telescope.ccd_gain.value

                    HeaderManipulation.edit_headers_with_fitsio(collection_temp, 
                                                                save_location = save_location, 
                                                                additional_name = "_master", 
                                                                **header_dict)
                            
                    fname_stem = Path(image_list_to_combine[0]).stem  # get the filename stem: crab_sdss_g.fits --get--> crab_sdss_g, it only contains one fits file
                    fname_suffix = Path(image_list_to_combine[0]).suffix # get the filename suffix: crab_sdss_g.fits --get--> .fits
                    name_path = Path(save_location) / f"{fname_stem}_master{fname_suffix}"  # the path to the copied file
                    new_name_path = name_path.with_name(f"Master_{object_name}_{filter}{fname_suffix}")  # the new path
                    name_path.rename(new_name_path)  # rename the file, notice that you have to use the full path to rename, not just the file name
                    print("----------------------------------------------------------\n")

                elif image_number > 1:
                    print(f"Using {image_filenames} to combine the master {object_name} image in {filter} filter!")
                    if method == "sigma clip":
                        if image_type == "Flat":
                            if scale_function == "auto":
                                scale_function_ = Combine.inv_median
                            
                            combined = ccdpro_combine(image_list_to_combine, 
                                                      method = "average", scale = scale_function_,
                                                      sigma_clip=True, sigma_clip_low_thresh=sigma_clip_low, sigma_clip_high_thresh=sigma_clip_high,
                                                      sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)
                        
                        elif image_type == "Light":
                            combined = ccdpro_combine(image_list_to_combine, 
                                                      method = "average", 
                                                      sigma_clip = True, sigma_clip_low_thresh=sigma_clip_low, sigma_clip_high_thresh=sigma_clip_high, 
                                                      sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)
                    elif method == "minmax clip":
                        raise ValueError("This function isn't ready yet!")

                    else:
                        raise ValueError("Unsupported combining method! What supported now is only sigma clip!")

                    combined.meta["IMTYPE"] = f"Master {image_type}"
                    combined.meta["NCOMBINE"] = int(image_number)
                    if image_type == "Light":
                        rdnoise = combined.meta["RDNOISE"]
                        combined.meta["RDNOISE"] = rdnoise/np.sqrt(image_number)  # update the RDNOISE
                        gain = combined.meta["GAIN"]
                        combined.meta["GAIN"] = gain*image_number  # update the GAIN
                    combined.write(save_location / f"Master_{object_name}_{filter}.fits", overwrite = True)
                    print(f"{object_name} in {filter} filter combined!")
                    print("----------------------------------------------------------\n")

        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(self._image_collection, rescan = True)

        return
