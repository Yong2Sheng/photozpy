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
from ..convenience_functions import convert_coords, estimate_background, get_alain_image
import pandas as pd
from ccdproc import ImageFileCollection
from ..mimage_collection import mImageFileCollection
import matplotlib.pyplot as plt
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.coordinates import SkyCoord
import shutil
from regions import PixCoord, CirclePixelRegion, CircleSkyRegion, Regions
import warnings
warnings.simplefilter("ignore", UserWarning)

class PhotozRegions():

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

        self._source_catalog_df = pd.read_csv(self._source_catalog_path, sep = ",")

   
    @staticmethod
    def generate_regions(region_dir, filter_name, coord, radius = 5.0):

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
        str_coord = coord.to_string(precision = 10).split(" ")  
        # this only works for one set of skycoords. There will be errors if you have more than a set of skycoords
        # a set of coord: str_coord.to_string(precision = 6)  ==> '105.130269 -66.179291'  <== this is string
        # two sets of coord: str_coord.to_string(precision = 6)  ==> ['105.130269 -66.179291', '105.169080 -66.158375']  <== this is list
        # since I always do photometry for only one source with Swift, so I am safe for now....
        x = str_coord[0]
        y = str_coord[1]
        
        
        with open(region_path, "w") as f:
            f.write("# Region file format: DS9 \n")
            f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
            f.write(f"{coord_frame}\n")
            f.write(f'circle({x},{y},{radius}")')

        return region_path
    
    def get_initial_bkg_regions(self):
        
        """
        Generate the initial bkg regions based on the catalog coordinate
        
        """
        for collection in self._mcollection:
            
            collection_path = collection.location
            bkg_path = Path(collection_path) / "bkg.reg"
            source_name = collection_path.parts[-1].replace("_", " ")
            if bkg_path.exists():
                print(f"The bkg region file has been generated at {bkg_path}!")
                
            else:
                print(collection_path)
                print(source_name)
                    
                # read the source coordinate from the source catalog
                ra = self._source_catalog_df[self._source_catalog_df["name"] == source_name]["ra"].to_numpy()[0]
                dec = self._source_catalog_df[self._source_catalog_df["name"] == source_name]["dec"].to_numpy()[0]
                sky_coord = SkyCoord(ra = ra, dec = dec, unit = "deg", frame = "fk5")
                print(f"ra = {ra}, dec = {dec}")
                
                PhotozRegions.generate_regions(region_dir = collection_path, filter_name = "bkg", coord = sky_coord, radius = 30)
            
            print("---------------------------------------------------")
            
        return
                
    def get_source_regions(self, box_size = 11, save_image = True, verbose = True, bkg_region_template = False, centroid_method = centroid_quadratic, aladin = True):

        """
        image_path : str or pathlib.Path
            The path to the fits image file
        coords : astropy.coordinates.sky_coordinate.SkyCoord
            The coordinates of the sources to fit the centroids
        """

        for collection in self._mcollection:
            
            collection_path = collection.location
            source_name = collection_path.parts[-1].replace("_", " ")
            print(f"Generating source regions for {source_name}")

            for image_path in collection.files_filtered(include_path = True, **{"SUMTYP" : "FINAL"}):
                
                image_path = Path(image_path)
                # read image data
                ccddata = CCDData.read(image_path, hdu = 1)
                array_data = ccddata.data
                filter_name = ccddata.header["FILTER"]

                # get the bkg substratced data
                bkg, array_data_no_bkg = estimate_background(array_data = array_data)

                # read the source coordinate from the source catalog
                ra = self._source_catalog_df[self._source_catalog_df["name"] == source_name]["ra"].to_numpy()[0]
                dec = self._source_catalog_df[self._source_catalog_df["name"] == source_name]["dec"].to_numpy()[0]
                sky_coord = SkyCoord(ra = ra, dec = dec, unit = "deg", frame = "fk5")

                #convert to the pixel coordinate
                pixel_coord = convert_coords(wcs = ccddata.wcs, skycoords = sky_coord, pixelcoords = None, verbose = False)


                # get the centroid pixelcoords
                x_centroids , y_centroids = centroid_sources(array_data_no_bkg, 
                                                             pixel_coord[0], pixel_coord[1], 
                                                             box_size = box_size, 
                                                             centroid_func = centroid_method)

                if np.isnan(x_centroids.astype(float)).any(): # check is any of the fitted pixel coordinate is NaN
                    
                    print(f"The fitting for {filter_name} using {centroid_method.__name__} failed, switch to {centroid_com.__name__} instead!")
                    
                    x_centroids , y_centroids = centroid_sources(array_data_no_bkg, 
                                                                 pixel_coord[0], pixel_coord[1], 
                                                                 box_size = box_size, 
                                                                 centroid_func = centroid_com)
                    
                if verbose == True:
                    print(f"The centroid pixel for {filter_name} is ({x_centroids},{y_centroids}).")
                
                # convert the centroid pixelcoords to skycoords
                centroid_skycoord = convert_coords(wcs = ccddata.wcs, 
                                                   skycoords = None, 
                                                   pixelcoords = np.array([x_centroids,y_centroids]).T, 
                                                   verbose = False)[0]  # the converted Skycoord is a list, so event if there is only one set skycoord, I have to use [0] to convert it from a list of skycoord to a skycoord


                src_region_path = PhotozRegions.generate_regions(region_dir = image_path.parent, 
                                                                 filter_name = filter_name, 
                                                                 coord = centroid_skycoord, 
                                                                 radius = 5.0)
                # generate the background region for all the filters from the background template
                if bkg_region_template is not None:
                    bkg_template_path = image_path.parent / bkg_region_template  # the template must be stored in the image directory
                    bkg_region_path = image_path.parent  / f"{filter_name}_bkg.reg"
                    shutil.copy(bkg_template_path, bkg_region_path)
                    
                if save_image == True:
                    PhotozRegions.plot_regions(array_data = array_data, wcs = ccddata.wcs, filter_name = filter_name, source_name = source_name,
                                               src_region_path = src_region_path, bkg_region_path = bkg_region_path, 
                                               other_coords = [sky_coord], other_coord_labels = ["original source location"],
                                               save_dir = image_path.parent, save_image = True)
            
            if aladin is True:
                print(f"Querying aladin image for {source_name}")
                sky_region_ = CircleSkyRegion(sky_coord, radius = 5*u.arcsec)  # the original sky region defined by the source catalog
                get_alain_image(wcs = ccddata.wcs, save_dir = image_path.parent, sky_region = sky_region_, 
                                min_cut = 0.5, max_cut = 99.5, logstretch = 10, source_name = source_name, hips = "CDS/P/DSS2/blue")
                
            if verbose == True:
                print("--------------------------------------------------------------------")

        return


    @staticmethod
    def plot_regions(image_path = None, hdu = 1, array_data = None, wcs = None, filter_name = None, source_name = None, 
                     src_region_path = None, bkg_region_path = None, other_coords = None, other_coord_labels = None,
                     image_cutout = ["full_image", 100, 20, 120], log_stretch = 3000, save_dir = None, save_image = False):
        """
        Plot regions
        
        Parameters
        ----------
        other_coords : list
        other_coord_labels : list
        """
    
        tick_fontsize = 12
        subtitle_fontsize = 15
        subplot_labelsize = 12
    
    
        # read image data
        if image_path is not None:
            image_path = Path(image_path)
            ccddata = CCDData.read(image_path, hdu = hdu)
            array_data = ccddata.data
            wcs = ccddata.wcs
            # determine the save location
            if save_dir is None:
                save_dir = image_path.parent
        elif array_data is None or wcs is None:
            raise ValueError("Please provide both array_data and wcs!")
        else:
            if save_dir is None:
                raise ValueError("Please provide save directory of the image.")
            else:
                save_dir = Path(save_dir)
    
        # read regions
        if src_region_path is None and bkg_region_path is None:
            raise TypeError("You must provide at least one region file.")
            
        if src_region_path is not None:
            source_skyregion = Regions.read(src_region_path, format = "ds9")[0]
            source_pixelregion = source_skyregion.to_pixel(wcs)
            
        if bkg_region_path is not None:
            bkg_skyregion = Regions.read(bkg_region_path, format = "ds9")[0]
            bkg_pixelregion = bkg_skyregion.to_pixel(wcs)
    
        # start to make plots
        fig, axs = plt.subplots(2, 2, figsize = (15, 15), subplot_kw=dict(projection = wcs))
        fig.suptitle(f"Region summary for filter {filter_name} of source {source_name}", y=0.92, fontsize = 20)
        norm = ImageNormalize(stretch=LogStretch(log_stretch))
    
        # plot regions to subplots
        if src_region_path is not None:
            # full image, upper left
            axs[0,0].scatter(source_pixelregion.center.x, source_pixelregion.center.y, marker = "+", s = 10, color = "lime", label = "source region center")
            
            # source (and background), upper right
            axs[0,1].scatter(source_pixelregion.center.x, source_pixelregion.center.y, marker = "+", s = 10, color = "lime", label = "source region center")
            source_pixelregion.plot(ax = axs[0,1], color='lime', lw=1.0, label = "source region")
            
            # source zoomed in, lower left
            axs[1,0].scatter(source_pixelregion.center.x, source_pixelregion.center.y, marker = "+", s = 10, color = "lime", label = "source region center")
            source_pixelregion.plot(ax = axs[1,0], color='lime', lw=1.0, label = "source region")
    
        if bkg_region_path is not None:
            # full image, upper left
            bkg_pixelregion.plot(ax = axs[0,0], color='red', lw=1.0, label = "background region")
    
            # source (and background), upper right
            bkg_pixelregion.plot(ax = axs[0,1], color='red', lw=1.0, label = "background region")
    
            # backgroud, lower right
            bkg_pixelregion.plot(ax = axs[1,1], color='red', lw=1.0, label = "background region")
    
    
        # plot images to subplots
        # the full image
        axs[0,0].imshow(array_data, origin='lower', norm=norm, cmap='Greys_r', interpolation='nearest')
        axs[0,0].legend()
        axs[0,0].set_title("Full sky image", fontsize = subtitle_fontsize)
    
        # source (and background), upper right
        axs[0,1].imshow(array_data, origin='lower', norm=norm, cmap='Greys_r', interpolation='nearest')
        axs[0,1].set_xlim(source_pixelregion.center.x - image_cutout[1], source_pixelregion.center.x + image_cutout[1])
        axs[0,1].set_ylim(source_pixelregion.center.y - image_cutout[1], source_pixelregion.center.y+ image_cutout[1])
        axs[0,1].legend()
        axs[0,1].set_title('Source (and background) region)', fontsize = subtitle_fontsize)
    
        # the source regions zoomed in, side wdith = 100
        axs[1,0].imshow(array_data, origin='lower', norm=norm, cmap='Greys_r', interpolation='nearest')
        axs[1,0].set_xlim(source_pixelregion.center.x - image_cutout[2], source_pixelregion.center.x + image_cutout[2])
        axs[1,0].set_ylim(source_pixelregion.center.y - image_cutout[2], source_pixelregion.center.y+ image_cutout[2])
        axs[1,0].legend()
        axs[1,0].set_title('Source region zoomed in', fontsize = subtitle_fontsize)
    
        # the background regions zoomed in side width = 120
        axs[1,1].imshow(array_data, origin='lower', norm=norm, cmap='Greys_r', interpolation='nearest')
        axs[1,1].set_xlim(bkg_pixelregion.center.x - image_cutout[3], bkg_pixelregion.center.x + image_cutout[3])
        axs[1,1].set_ylim(bkg_pixelregion.center.y - image_cutout[3], bkg_pixelregion.center.y+ image_cutout[3])
        axs[1,1].legend()
        axs[1,1].set_title('Background region zoomed in', fontsize = subtitle_fontsize)
        axs[1,1].tick_params(which='major', labelsize=12)
    
        # create RA Dec grids
        for ax in axs.flatten():
            ax.coords.grid(True, color='white', ls='dotted')
            ax.coords[0].set_axislabel('Right Ascension (J2000)', fontsize=tick_fontsize)
            ax.coords[1].set_axislabel('Declination (J2000)', fontsize=tick_fontsize)
            ax.tick_params(which='major', labelsize=tick_fontsize)
            ax.legend(fontsize = subplot_labelsize)
            if other_coords is not None:
                for coord_, label_ in zip(other_coords, other_coord_labels):
                    ax.scatter(coord_.ra.deg, coord_.dec.deg, marker = "+", s = 10, color = "cyan", label = label_)
    
        fig.savefig(save_dir / f"{source_name}_{filter_name}.png", dpi = 300, bbox_inches = "tight")

        plt.close()

        return
            

                
        
        
