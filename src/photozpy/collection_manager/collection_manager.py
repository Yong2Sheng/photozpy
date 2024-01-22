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
import matplotlib.pyplot as plt
import matplotlib
from astropy.visualization import astropy_mpl_style, simple_norm
from astropy.io import fits
from matplotlib.colors import LogNorm
from tqdm import tqdm

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
    def _plot_fits_image(fits_path):
        """
        Plot a fits image.

        Parameters
        ----------
        path: str or pathlib.Path; the path to the image.

        Returns
        -------
        None
        """

        plt.style.use(astropy_mpl_style)
        matplotlib.use('Agg')
        
        fits_path = Path(fits_path)  # path of fits file
        image_name = fits_path.stem + ".png"
        image_path = fits_path.parent/image_name
        
        fits_data = fits.getdata(fits_path, ext=0)
        norm = simple_norm(fits_data, 'log', percent= 99.9)
        plt.figure()
        plt.imshow(fits_data, cmap='gray', norm=norm)
        plt.colorbar()
        plt.savefig(image_path, dpi=300)
        plt.clf()
        plt.close("all")

    def plot_collection(image_collection, headers_values = None):

        """
        Plot the fits images in the image_collection. The image collection can be filtered.

        Parameters
        ----------
        image_collection: ccdproc.ImageFileCollection; the image collection you want to plot
        headers_values: dict; the header names and corresponding values. The values should be in list.
                              {"header": ["value"]} and the list can contain multiple values.

        returns
        -------
        collection_to_plot: ccdproc.ImageFileCollection; the image collection you plot
        """

        if headers_values != None:
            collection_to_plot = CollectionManager.filter_collection(image_collection, **headers_values)
        else:
            collection_to_plot = image_collection

        image_path_to_combine = list(collection_to_plot.files_filtered(include_path = True))
        for i in tqdm(image_path_to_combine):
            CollectionManager._plot_fits_image(i)

        return collection_to_plot

        

        
        