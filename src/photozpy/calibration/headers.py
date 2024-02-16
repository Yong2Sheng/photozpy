"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- correct the headers
- the purpose is to standarize the fits headers since they will be used to identify the images.
"""

from ..telescope import Telescope
from ccdproc import ImageFileCollection
from pathlib import Path
from tqdm import tqdm
from ..collection_manager import CollectionManager
from astropy.time import Time
from pathlib import Path
from astropy.io import fits

class HeaderCorrection():

    def __init__(self, image_collection, telescope, save_location = "", target_dict = None, overwrite = True):

        """
        The input image collection.

        Parameters
        ----------
        image_collection: ccdproc.ImageFileCollection; the collection of images
        telescope:
        save_location: pathlib.Path; the location to save the edited files
        target_dict: dictionary; map the common part of the file names to a target, i.e. all the 

        Returns
        -------
        None
        """

        if not isinstance(image_collection, ImageFileCollection):
            raise TypeError("image_collection shoule be an ImageFileCollection objec!")
        else:
            # refresh the full collection
            self._image_collection = CollectionManager.refresh_collection(image_collection, rescan = True)

        if not isinstance(telescope, Telescope):
            raise TypeError("telescope should be a Telescope oject!")
        else:
            self._telescope = telescope

        self._target_dict = target_dict

        return


    def correct_filter_headers(self, input_collection = None, overwrite = True, save_location = "", from_file_name = False, **filter_dict):
        """
        Correct the filter names to standard ones for future analysis"

        Parameters
        ----------
        save_location: str; the location to save the new files.
        overwrite: bolean; overwrite the files or not if the files already exist.
        filter_dict: customized filter maping dictionary.
        
        Returns
        -------
        new_collection: ccdproc.ImageFileCollection; The new image collection.
        """

        print("Editting filter headers......")
        

        if input_collection == None:
            # refresh the full collection
            self._image_collection = CollectionManager.refresh_collection(self._image_collection, rescan = True)
            input_collection = self._image_collection

        # only include IMTYPE is Flat or Light
        flat_collection = CollectionManager.filter_collection(input_collection, **{"IMTYPE": ["Flat"]})
        light_collection = CollectionManager.filter_collection(input_collection, **{"IMTYPE": ["Light"]})
        fl_collection = ImageFileCollection(input_collection.location, 
                                            filenames = flat_collection.files + light_collection.files)
        # note that although it only works on flat and light image, it will eventually return a image collection
        # that contains all the original files in the input_collection

        for hdu, file_name in tqdm(fl_collection.hdus(overwrite = overwrite, save_location = save_location, return_fname = True)):
            if from_file_name == False:
                old_filter = hdu.header["FILTER"]
                name = hdu.header["OBJECT"]
            else:
                old_filter = file_name  # if the fits file has no value for the FILTER header, it will use the file name as the filter name since the file name contains the filter name as well.
            if old_filter in self._telescope.filters:
                pass
            else:
                new_filter = HeaderManipulation.map_filters(old_filter, **filter_dict)
                hdu.header["FILTER"] = new_filter

        if save_location == "":
            new_collection = ImageFileCollection(location = input_collection.location, filenames = input_collection.files)
        else:
            new_collection = ImageFileCollection(location = save_location, filenames = input_collection.files)
            
        self._image_collection = new_collection

        print("Filter header edition completed!")
        print("----------------------------------------\n")

        return
        

    def correct_headers_by_filename(self, save_location = "", overwrite = True, type = None):
        
        """
        Correct the headers by the filename and/or type (Bias, Dark, Flat, Light). 
        It basically uses the common string in the filename to identify the image type.
        You can also directly declare the image type by the type parameter.

        Parameters
        ----------
        save_location: str; the location to save the new files.
        overwrite: bolean; overwrite the files or not if the files already exist.
        type: str; the type of the image (Bias, Dark, Flat, Light)

        Returns
        -------
        new_collection: ccdproc.ImageFileCollection; The new image collection.
        """
        

        print("Editing headers by file names......")
        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(self._image_collection, rescan = True)

        fits_files = self._image_collection.files
        
        for i in tqdm(fits_files):
            fits_path = self._image_collection.location/i
            if "Bias" in fits_path.stem or "bias" in fits_path.stem or type == "Bias":
                header_dict = {"IMTYPE": ("Bias", None),
                               "GAIN": (self._telescope.ccd_gain.value, self._telescope.ccd_gain.unit.to_string()),
                               "RDNOISE": (self._telescope.ccd_rdnoise.value, self._telescope.ccd_rdnoise.unit.to_string()),
                               "OBJECT":  ("Bias", None),
                               "TELESCOP": self._telescope.telescope, 
                               "BUNIT": "ADU"}
                # create a image collection that only contains this single image
                one_image_collection = ImageFileCollection(location = self._image_collection.location, filenames = i)
                _ = HeaderManipulation.correct_headers(one_image_collection, save_location = save_location, 
                                                       time_transformation = True, **header_dict)

            elif "Dark" in fits_path.stem or "dark" in fits_path.stem or type == "Dark":
                header_dict = {"IMTYPE": ("Dark", None),
                               "GAIN": (self._telescope.ccd_gain.value, self._telescope.ccd_gain.unit.to_string()),
                               "RDNOISE": (self._telescope.ccd_rdnoise.value, self._telescope.ccd_rdnoise.unit.to_string()),
                               "OBJECT": ("Dark", None),
                               "TELESCOP": self._telescope.telescope, 
                               "BUNIT": "ADU"}
                # create a image collection that only contains this single image
                one_image_collection = ImageFileCollection(location = self._image_collection.location, filenames = i)
                _ = HeaderManipulation.correct_headers(one_image_collection, save_location = save_location, 
                                                       time_transformation = True, **header_dict)
                
            elif "Flat" in fits_path.stem or "flat" in fits_path.stem or type == "Flat":
                header_dict = {"IMTYPE": ("Flat", None),
                               "GAIN": (self._telescope.ccd_gain.value, self._telescope.ccd_gain.unit.to_string()),
                               "RDNOISE": (self._telescope.ccd_rdnoise.value, self._telescope.ccd_rdnoise.unit.to_string()),
                               "OBJECT": ("Flat", None),
                               "TELESCOP": self._telescope.telescope, 
                               "BUNIT": "ADU"}
                # create a image collection that only contains this single image
                one_image_collection = ImageFileCollection(location = self._image_collection.location, filenames = i)
                _ = HeaderManipulation.correct_headers(one_image_collection, save_location = save_location, 
                                                       time_transformation = True, **header_dict)
                
            else:
                # For light type, we need to add the target name
                if self._target_dict == None:
                    raise TypeError("Please provide a target dictionary!")
                else:
                    target_name = [value for key,value in self._target_dict.items() if key in fits_path.stem]
                    target_name = [*set(target_name)]
                    if len(target_name) != 1:
                        pass
                        print(f"There is no {fits_path.stem} file or the dictionary doesn't contain this file! Skipping.....")
                    else:
                        target_name = target_name[0]
                    header_dict = {"IMTYPE": ("Light", None),
                                   "GAIN": (self._telescope.ccd_gain.value, self._telescope.ccd_gain.unit.to_string()),
                                   "RDNOISE": (self._telescope.ccd_rdnoise.value, self._telescope.ccd_rdnoise.unit.to_string()),
                                   "OBJECT": target_name,
                                   "TELESCOP": self._telescope.telescope, 
                                   "BUNIT": "ADU"}
                    # create a image collection that only contains this single image
                    one_image_collection = ImageFileCollection(location = self._image_collection.location, filenames = i)
                    _ = HeaderManipulation.correct_headers(one_image_collection, save_location = save_location, 
                                                           time_transformation = True, **header_dict)

        # refresh the image collection
        if save_location == "":
            new_image_collection = ImageFileCollection(location = self._image_collection.location, filenames = fits_files)
        else:
            new_image_collection = ImageFileCollection(location = save_location, filenames = fits_files)
        self._image_collection = new_image_collection

        print("Header edition by file names completed!")
        print("----------------------------------------\n")

        return

    @property
    def image_collection(self):
        return self._image_collection

    @property
    def telescope(self):
        return self._telescope


class HeaderManipulation():

    def __init__():

        pass

    @staticmethod
    def check_headers(image_collection, **headers_values):
        """
        Check if the headers and their values are consistent in the image collection
    
        Parameters
        ----------
        image_collection: ccdproc.ImageFileCollection; the collection of images
        headers_values: dict; the header names and corresponding values. The values should be in tuple.
                              If there is only one value for the header, use (value,)
    
        Returns
        -------
        final : boolean; True or False
        """
        
        results = []
        for header, value in headers_values.items():
            # get the header values from the image collection
            values_in_collection = image_collection.values(header)
            values_in_collection = [*set(values_in_collection)]  # remove the repeated values to speed up
    
            result = []
            for val in value:
                if val in values_in_collection:
                    result += [True]
                else:
                    result += [False]
                    
            results += [all(result)]
    
        final = all(results)
    
        return final

    @staticmethod
    def correct_headers(input_collection, save_location = "", overwrite = True, time_transformation = False, **headers_values):

        """
        Change or add headers and values.

        Parameters
        ----------
        input_collection: ccdproc.ImageFileCollection; the collection of images.
        save_location: str; the location to save the new files.
        overwrite: bolean; overwrite the files or not if the files already exist.
        time_transformation: boolean; transformation the time format from UT1 to MJD and add a MJD-OBS header to meet the fits standard
        headers_values: dict; the dictionary of the headers and their values.

        Returns
        ------
        new_collection: ccdproc.ImageFileCollection; The new image collection.
        """

        if not headers_values:  # python evalute empty object to be False
            raise ValueError("Please provide at least one header and one header values!")

        # crate the save_location if not exist
        if save_location != "":
            Path(save_location).mkdir(parents=False, exist_ok=True)
        
        for hdu in input_collection.hdus(overwrite = overwrite, save_location = save_location):
            for key, value in headers_values.items():
                hdu.header[key] = value
                if time_transformation == True:
                    date_time = Time(hdu.header["DATE-OBS"], format = "isot", scale = "ut1")
                    mjd_time = date_time.mjd
                    hdu.header["MJD-OBS"] = mjd_time

        # Here I want to return a new collection of the modified files because the summary of the input image collection won't be updated by itself
        if save_location == "":
            new_collection = ImageFileCollection(location = input_collection.location, filenames = input_collection.files)
        else:
            new_collection = ImageFileCollection(location = save_location, filenames = input_collection.files)

        return new_collection

    @staticmethod
    def map_filters(name, **filter_dict):
        """
        Map filter names to standard filter names. 

        Paremeters
        ----------
        name: the input filter name to be matched to a standard one.
        filter_dict: customized filter maping dictionary.

        Return
        ------
        mapped_filter: str; the name of the mapped standard filter.
        """

        if filter_dict == {}:
            # use the default mappers
            mappers = {"Sloan g'2": "SDSS_g'", 
                       "Sloan r'2": "SDSS_r'", 
                       "Sloan i'2": "SDSS_i'", 
                       "Sloan z'2": "SDSS_z'", 
                       "SDSS_g": "SDSS_g'", 
                       "SDSS_r": "SDSS_r'", 
                       "SDSS_i": "SDSS_i'", 
                       "SDSS_z": "SDSS_z'",
                       "Bessell I": "Bessell_I"}
        else:
            mappers = filter_dict
            
        mapped_filters = [value for key,value in mappers.items() if key in name]
        mapped_filters = [*set(mapped_filters)] # remove the duplicated filter names
        
        if len(mapped_filters) == 0:
            raise ValueError("Mapping the input filter failed! No mapped filter {name} was found!")
        elif len(mapped_filters) > 1:
            raise ValueError("Mapping the input filter failed! More than two mapped filters are found!")
        else:
            return mapped_filters[0]
    

    @staticmethod
    def get_header_values(image_collection, header):
        """
        Get the header values from an image collection. The repeated values will be removed.

        Paremeters
        ----------
        image_collection
        headers

        Returns
        -------
        header_values: list
        """

        return image_collection.values(header, unique = True)

    @staticmethod
    def edit_headers_with_fitsio(image_collection, save_location = "", additional_name = None, overwrite = True, **headers_values, ):

        """
        This edits the header of multi-extension fits files since the iteration of ImageCollection doesn't work on
        multi-extension fits files.

        Parameters
        ----------
        image_collection
        save_location
        headers_values: dict; the dictionary of the headers you want to add/edit
        additional_name: str; the additional name appending to the original file name

        Returns
        -------
        None
        """

        file_paths = image_collection.files_filtered(include_path = True)

        for file_path in file_paths:

            file_path = Path(file_path)
            if save_location != "":
                save_location = Path(save_location)
                if additional_name == None:
                    save_location = save_location/file_path.name
                else:
                    file_name = file_path.stem + "_master" + file_path.suffix
                    save_location = save_location/file_name
            else:
                if additional_name == None:
                    save_location = file_path
                else:
                    file_name = file_path.stem + "_master" + file_path.suffix
                    save_location = file_path.parent/file_name
    
            hdul = fits.open(file_path)
            for key, value in headers_values.items():
                hdul[0].header[key] = value
            hdul.writeto(save_location, overwrite = overwrite)
            hdul.close()

        
    