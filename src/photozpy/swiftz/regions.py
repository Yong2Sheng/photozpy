"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
Deals with centroid fitting and region generation.
"""

from photutils.centroids import centroid_quadratic, centroid_com, centroid_sources, centroid_2dg
from astropy.nddata import CCDData
import astropy.units as u
import numpy as np
from pathlib import Path
from ..convenience_functions import convert_coords
import pandas as pd

class Regions():

    def __init__(self, image_collection, telescope, source_catalog_path):

        """
        image_collection : ccdproc.image_collection.ImageFileCollection or photozpy.mimage_collection.mImageFileCollection
            The input image collections
        telescope : photozpy.telescope.telescope.Telescope
            The swift telescope spec.
        source_catalog : str or pathlib.Path
            The csv file that contains the source names and source coordinates
        """

        # switch image_collection from ImageFileCollection to mImageFileCollection to make it standard for the pipeline
        if isinstance(image_collection, ImageFileCollection):
            self._mcollection = mImageFileCollection(location = image_collection.location, filenames = image_collection.files)

        elif isinstance(image_collection, mImageFileCollection):
            self._mcollection = image_collection

        self._telescope = telescope

        self._source_catalog_path = source_catalog_path

        self._source_catalog_df = pd.read_csv(_source_catalog_path, sep = ",")

   
    @staticmethod
    def generate_source_regions(region_dir, filter_name, coord, radius = 5.0):

        """
        Generate the circular region files.

        Parameters
        ----------
        region_dir : str or pathlib.Path
            The directory to save the region files
        filter_name : str
            The filter name to be used  to name the region files
        coord : astropy.coordinates.sky_coordinate.SkyCoord
            The sky coordinate of the source
        coord_system : str, optional
            The coordinate system for the region file
        radius : float
            The radius of the circular region in arcsecond
        """

        region_dir = Path(region_dir)
        src_region_fname = f"{filter_name}.reg"
        region_path = region_dir / src_region_fname

        coord_frame = coord.frame.name
        str_coord = coord.to_string().split(" ")
        x = str_coord[0]
        y = str_coord[1]
        
        
        with open(region_path, "w") as f:
            f.write("# Region file format: DS9 \n")
            f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
            f.write(f"{coord_frame}\n")
            f.write(f'circle({x},{y},{radius}")')

        return
            
    
    def get_source_regions(self):

        """
        image_path : str or pathlib.Path
            The path to the fits image file
        coords : astropy.coordinates.sky_coordinate.SkyCoord
            The coordinates of the sources to fit the centroids
        """

        for collection in self._mcollection:
            
            collection_path = collection.location

            for image_path in collection.files_filtered(include_path = True):
                
                # read image data
                ccddata = CCDData.read(image_path)
                array_data = ccddata.data
        
                pixel_coord = convert_coords(image_path = image_path, wcs = None, skycoords = sky_coord, pixelcoords = None, verbose = False)

        
        