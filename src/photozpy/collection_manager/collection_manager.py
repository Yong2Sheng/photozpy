"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- Combine a series of registered images.
- Check image types before combine.
"""

from ccdproc import ImageFileCollection
from pathlib import Path
from itertools import product
import copy
from ..convenience_functions import *
from tqdm import tqdm
from astropy.io import fits

class CollectionManager():

    def __init__(self):

        pass

    @staticmethod
    def refresh_collection(image_collection, rescan = False):
        """
        Refresh the image collection by reading the files. You can also rescan all the fits files to include newly created fits files.

        Paremeters
        ----------
        image_collection: ccdproc.ImageFileCollection; the image collection you want to refresh
        rescan: boolean; set to True if you want to rescan all the fits files in the location of image_collection to include new fits files.

        Returns
        -------
        new_image_collection: ccdproc.ImageFileCollection; the refreshed image collection
        """

        if rescan == False:
            new_image_collection = ImageFileCollection(location = Path(image_collection.location), 
                                                       filenames = image_collection.files)
        elif rescan == True:
            new_image_collection = ImageFileCollection(location = Path(image_collection.location), 
                                                       glob_include = "*.fits")

        return new_image_collection


    @staticmethod
    def unwarp_dictionary(dictionary):
    
        """
        Unwarp a dictionary to produce a dictionary list whose element is a dictionary that contains all the combinations of 
        values for the headers.
    
        Parameters
        ----------
        dictionary: dict; the input dictionary
    
        Returns
        -------
        dict_list: list; the list of all dictionaries that exhausts all the combination of the values.
        """

        # make sure all the values are lists instead of strings or numbers
        for header, value in dictionary.items():
            if not isinstance(value, list):
                dictionary[header] = [value]
        
        # get all the combinations of values from different headers
        values = list(product(*list(dictionary.values())))
    
        headers = list(dictionary.keys())
    
        dict_list = []
        for i in values:
            dict_ = dict(zip(headers, i))  # assemble headers and values. Note that headers are always the same. We just iterate through the combination of the values.
            dict_list += [dict_]
    
        return dict_list
    
    @staticmethod
    def filter_collection(image_collection, **headers_values):
        """
        Filter the image collection based on the headers and values. One header can have multiple corresponding values.
    
        Paremeters
        ----------
        image_collection:
        headers_values: dict; the header names and corresponding values. The values should be in list.
                              {"header": ["value"]} and the list can contain multiple values.
    
        Returns
        -------
        new_image_collection
        """

        # refresh the full collection --> I should NOT refresh the image collection!
        #image_collection = CollectionManager.refresh_collection(image_collection, rescan = True)

        dict_list = CollectionManager.unwarp_dictionary(headers_values)
    
        files = []
        for dict in dict_list:
            files_temp = list(image_collection.files_filtered(**dict))
            files += files_temp
    
        # remove duplicate file names
        files = [*set(files)]
    
        new_image_collection = ImageFileCollection(location = Path(image_collection.location), filenames = files)
    
        return new_image_collection


    @staticmethod
    def delete_images(image_collection, image_list):

        """
        Delete the images in the image_list from the image_collection.

        Parameters
        ----------
        image_collection
        image_list: list; the list of image file names to be deleted

        Returns
        -------
        new_image_collection
        """

        if not isinstance(image_list, list):
            image_list = [image_list]

        file_names = copy.deepcopy(image_collection.files)

        for i in image_list:
            file_names.remove(i)

        new_image_collection = ImageFileCollection(location = Path(image_collection.location), filenames = file_names)

        return new_image_collection
        
    @staticmethod
    def plot_collection(image_collection, save_location = None, norm_percent = 99.9, **headers_values):

        """
        Plot the fits images in an image collection. The image collection can be filtered before plotting.

        Parameters
        ----------
        image_collection : NoneType or ccdproc.ImageFileCollection
            The image collection to be plotted. 
        save_location : NoneType, str or pathlib.Path
            The path to save the plotted fits images.
        headers_values : dict
            The headers and corresponding values to filter the image collection before plotting.

        Returns
        -------
        None
        """

        collection_to_plot = CollectionManager.refresh_collection(image_collection, rescan = True)
        if headers_values is not {}:
            collection_to_plot = CollectionManager.filter_collection(collection_to_plot, **headers_values)
        fits_path_to_plot = list(collection_to_plot.files_filtered(include_path = True))

        for fits_path in tqdm(fits_path_to_plot):
            plot_image(fits_path = fits_path)

        return 

    @staticmethod
    def get_header_values(image_collection, headers, unique = True):

        """
        Get the header valus from the image collection. I wrote this because the ImageFileCollection.values() doesn't work!

        """

        if isinstance(headers, str):
            headers = [headers]

        file_paths = image_collection.files_filtered(include_path = True)

        values_list = []
        for file_path in file_paths:
            header_all = fits.getheader(file_path)
            for header in headers:
                value_ = header_all[header]
                values_list.append(value_)

        if unique:
            values_list = [*set(values_list)]

        return values_list
            

        

        

        

        
        