"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
Deals with centroid fitting and region generation.
"""

from photutils.centroids import centroid_quadratic, centroid_com, centroid_sources
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
from regions import CircleSkyRegion, Regions
import warnings
warnings.simplefilter("ignore", UserWarning)
import astropy.wcs as wcs

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
                    PhotozRegions._plot_regions(array_data = array_data, wcs = ccddata.wcs, filter_name = filter_name, source_name = source_name,
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
    def _plot_regions(image_path = None, hdu = 1, array_data = None, wcs = None, filter_name = None, source_name = None, 
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


class CCD_Regions():

    def __init__(self, image_collection, sources):

        # switch image_collection from ImageFileCollection to mImageFileCollection to make it standard for the pipeline
        if isinstance(image_collection, ImageFileCollection):
            self._mcollection = mImageFileCollection(image_dir = image_collection.location, filenames = image_collection.files)

        elif isinstance(image_collection, mImageFileCollection):
            self._mcollection = image_collection

        self.sources = sources

    def generate_ccd_regions(self, save_plots = True, hdu = 0, box_size = 51*u.pixel, **kwargs):

        """
        Generate the regions files for CCD analyze. It support circular and annulus regions.

        Parameters
        ----------
        **kwargs
            The keyword arguments are passed to `generate_region_files()`

        Returns
        -------
        """

        for image_collection in self._mcollection:
            for source_name, source_coords in self.sources:
                image_paths = image_collection.files_filtered(include_path = True, **{"OBJECT" : source_name})

                for image_path in image_paths:
                    image_path = Path(image_path)
                    image_parent_path = image_path.parent
                    ccddata = CCDData.read(image_path, hdu = 0)
                    image_filter_name = ccddata.header["FILTER"]
                    image_wcs = ccddata.wcs
                    image_array_data= ccddata.data
                    psf_fwhm = ccddata.header["FWHM"]*u.pixel
                    aperture_radius = psf_fwhm*3
                    inner_radius = psf_fwhm*5
                    outer_radius = psf_fwhm*8

                    # get the centroids of the sky coodinates
                    source_coords_centroids = get_centroids(sky_coords = source_coords, array_data = image_array_data, image_wcs = image_wcs, return_type = "sky_coord", 
                                                            box_size = int(box_size.value), verbose = False, centroid_method = centroid_quadratic)
                    
                    # get the pixel scale 
                    x, y = wcs.utils.proj_plane_pixel_scales(image_wcs)*u.deg/u.pixel  # x and y are pixel scale along two axis. They are usually very close.
                    pixel_scale = u.pixel_scale(x)  # I choose x axis as the pixel scale for conversion.
                    
                    # now the target names are generated automatically
                    target_names = [f"S{i}" for i in np.arange(source_coords_centroids.shape[0])]
                    
                    # generate and save the region files for source aperture and background annulus
                    src_region_path = generate_region_files(region_save_dir = image_parent_path, field_name = source_name, source_name = target_names, filter_name = image_filter_name, 
                                                            sky_coord = source_coords_centroids, region_type = "src", pixel_scale = pixel_scale, inner_radius = aperture_radius)
                    
                    
                    bkg_region_path = generate_region_files(region_save_dir = image_parent_path, field_name = source_name, source_name = target_names, filter_name = image_filter_name, 
                                                            sky_coord = source_coords_centroids, region_type = "bkg", pixel_scale = pixel_scale, inner_radius = inner_radius, outer_radius = outer_radius)
                    
                    # Plot regions
                    if save_plots is True:
                        plot_regions(image_array_data = image_array_data, image_wcs = image_wcs, image_filter_name = image_filter_name, field_name = source_name, 
                                     src_region_path = src_region_path, bkg_region_path = bkg_region_path, aladin_image = None,
                                     aladin_stretch = 3000, data_strentch = 5000, save_dir = image_parent_path, save_image = True)


                    
            

    

def generate_region_files(region_save_dir, field_name, source_name, filter_name, 
                          sky_coord, region_type = "src", pixel_scale = None, inner_radius = 5.0*u.pixel, outer_radius = None):

    """
    Generate the circular region files.

    Parameters
    ----------
    region_save_dir : str or pathlib.Path
        The directory to save the region files
    region_type : str
        The type of the region: src or bkg
    field_name : str
        The field name of the image
    source_name : str or list
        source names for each sky_coord
    filter_name : str
        The filter name to be used  to name the region files
    sky_coord : astropy.coordinates.sky_coordinate.SkyCoord
        The sky coordinate of the source
    pixel_scale : astropy.units.quantity.Quantity
        The pixel scale to convert pixel length to true size in angular units. 
        The pixel scale must have the unit of u.deg/u.pixel.
    inner_radius : float
        The inner radius of the circular region in arcsecond
    outer_radius : float
        The outer radius of the circular region in arcsecond
    """

    region_save_dir = Path(region_save_dir)
    region_fname = f"{field_name}_{filter_name}_{region_type}.reg"
    region_path = region_save_dir / region_fname

    coord_frame = sky_coord.frame.name
    str_coord = sky_coord.to_string(precision = 10)
    ra = [coord.split(" ")[0] for coord in str_coord]
    #print(ra)
    dec = [coord.split(" ")[1] for coord in str_coord]
    #print(dec)

    # check if the length of the source name equals to the length of the ra and dec
    if isinstance(source_name, str):
        source_name = [source_name]
    if len(source_name) != len(ra):
        raise ValueError("The number of ra and dec doesn't match the number of source names.")
    
    with open(region_path, "w") as f:
        f.write("# Region file format: DS9 \n")
        if region_type == "src":
            f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        elif region_type == "bkg":
            f.write('global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=1 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        f.write(f"{coord_frame}\n")
        for _ra, _dec, _name in zip(ra, dec, source_name):

            if inner_radius.unit.name == "pix":
                inner_radius_angular = inner_radius.to(u.arcsec, pixel_scale).value # convert radius from pixel to angular units
            if outer_radius is None:
                f.write(f'circle({_ra},{_dec},{inner_radius_angular}") # text={{{_name}_{region_type}}}\n')
            else:
                if outer_radius.unit.name == "pix":
                    outer_radius_angular = outer_radius.to(u.arcsec, pixel_scale).value
                f.write(f'annulus({_ra},{_dec},{inner_radius_angular}", {outer_radius_angular}") # text={{{_name}}}\n')

    return region_path


def get_centroids(sky_coords=None, pixel_coords=None, image_path=None, hdu=None, array_data=None, image_wcs=None, filter_name=None,
                  return_type = "pix_coord", box_size = 51, verbose = False, centroid_method = centroid_quadratic):

    """
    Find the centroids of the sources.
    
    Parameters
    ----------
    sky_coords : astropy.coordinates.SkyCoord
        The sky coordinates of 

    Returns
    -------
    """

    hdu = 1 if hdu is None else hdu
    
    if image_path is not None:
        
        image_path = Path(image_path)
        # read image data
        ccddata = CCDData.read(image_path, hdu = hdu)
        array_data = ccddata.data
        filter_name = ccddata.header["FILTER"]
        image_wcs = ccddata.wcs

    # get the bkg substratced data
    bkg, array_data_no_bkg = estimate_background(array_data = array_data)

    if sky_coords is not None:
        #convert to the pixel coordinate
        pixel_coords = convert_coords(wcs = image_wcs, skycoords = sky_coords, pixelcoords = None, verbose = False)

    # get the centroid pixelcoords
    x_centroids , y_centroids = centroid_sources(array_data_no_bkg, 
                                                 pixel_coords[:,0], pixel_coords[:,1], 
                                                 box_size = box_size, 
                                                 centroid_func = centroid_method)

    if np.isnan(x_centroids.astype(float)).any(): # check is any of the fitted pixel coordinate is NaN
        
        print(f"The fitting for {filter_name} using {centroid_method.__name__} failed, switch to {centroid_com.__name__} instead!")
        
        x_centroids , y_centroids = centroid_sources(array_data_no_bkg, 
                                                     pixel_coords[0], pixel_coords[1], 
                                                     box_size = box_size, 
                                                     centroid_func = centroid_com)

    if verbose == True:
        print(f"The centroid pixel for {filter_name} is ({x_centroids},{y_centroids}).")
    
    # convert the centroid pixelcoords to skycoords
    centroid_skycoord = convert_coords(wcs = image_wcs, 
                                       skycoords = None, 
                                       pixelcoords = np.array([x_centroids,y_centroids]).T, 
                                       verbose = False)  # the converted Skycoord is a list, so event if there is only one set skycoord, I have to use [0] to convert it from a list of skycoord to a skycoord

    if return_type == "pix_coord":
        return np.array([x_centroids,y_centroids]).T
    elif return_type == "sky_coord":
        return centroid_skycoord

    

def plot_regions(image_path = None, hdu = 1, image_array_data = None, image_wcs = None, image_filter_name = None, field_name = None, 
                 src_region_path = None, bkg_region_path = None, aladin_image = None,
                 image_cutout = ["full_image", 100, 20, 120], aladin_stretch = 3000, data_strentch = 5000, save_dir = None, save_image = False):
    
    if image_path is not None:
        image_path = Path(image_path)
        ccddata = CCDData.read(image_path, hdu = hdu)
        image_filter_name = ccddata.header["FILTER"]
        field_name = ccddata.header["OBJECT"]
        image_wcs = ccddata.wcs
        image_array_data= ccddata.data
        # determine the save location
        if save_dir is None:
            save_dir = image_path.parent
    elif image_array_data is None or image_wcs is None:
        raise ValueError("Please provide both array_data and wcs!")
    else:
        if save_dir is None:
            raise ValueError("Please provide save directory of the image.")
        else:
            save_dir = Path(save_dir)

    # read regions
    if src_region_path is None and bkg_region_path is None:
        raise TypeError("You must provide at least one region file.")
    n_src_regions = None
    n_bkg_regions = None
    # convert src sky regions to pixel regions
    if src_region_path is not None:
        src_sky_regions = Regions.read(src_region_path, format = "ds9")
        src_pixel_regions = [i.to_pixel(image_wcs) for i in src_sky_regions]
        n_src_regions = len(src_pixel_regions)
    # convert bkg sky regions to pixel regions
    if bkg_region_path is not None:
        bkg_sky_regions = Regions.read(bkg_region_path, format = "ds9")
        bkg_pixel_regions = [i.to_pixel(image_wcs) for i in bkg_sky_regions]
        n_bkg_regions = len(bkg_pixel_regions)
    
    # double check if the number of src and bkg regions are equal
    # create None list for the source or background regions that are not defined
    if n_src_regions is None and n_bkg_regions is None:
        if n_src_regions != n_bkg_regions:
            raise ValueError(f"The number of source and background are not equal: source regions: {n_src_regions}; background regions: {n_bkg_regions}.")
        else:
            n_subplots = n_src_regions + 1
            
    elif n_src_regions is None and n_bkg_regions is not None:
        src_pixel_regions = [None]*n_bkg_regions
        n_subplots = n_bkg_regions + 1
        
    elif n_src_regions is not None and n_bkg_regions is None:
        bkg_pixel_regions = [None]*n_src_regions
        n_subplots = n_src_regions + 1
        
    else:
        n_subplots = n_src_regions + 1
    
    # Create figure
    fig, axs = plt.subplots(n_subplots, 2, figsize = (15, 7.5*n_subplots), subplot_kw=dict(projection = image_wcs), constrained_layout=True)
    fig.suptitle(f"Region summary for {image_filter_name} filter of source {field_name}", fontsize = 20)
    text_fontsize = 15
    tick_fontsize = 15
    subplot_labelsize = 15
    subtitle_fontsize = 15
    marker_size = 10
    
    # plot the first raw: aladin image
    if aladin_image is None:
        aladin_result = get_alain_image(wcs = image_wcs, save_image=False)
    axs[0,0].imshow(aladin_result[0].data, origin='lower', norm = ImageNormalize(data = aladin_result[0].data, stretch = LogStretch(aladin_stretch)), cmap='Greys_r', interpolation='nearest')
    axs[0,0].set_title('Aladin Image', fontsize = subtitle_fontsize)
    
    # plot the first raw: data image
    axs[0,1].imshow(image_array_data, origin='lower', norm = ImageNormalize(data = image_array_data, stretch = LogStretch(data_strentch)), cmap='Greys_r', interpolation='nearest')
    if bkg_pixel_regions[0] is not None:  # use bkg region to show and lable the sources as the first option since it has larger radius and easier to see
        for bkg_region in bkg_pixel_regions: 
            bkg_region.plot(ax = axs[0,1], lw=1.0, label = bkg_region.meta["text"], color = np.random.rand(3,))
            axs[0,1].text(bkg_region.center.x, bkg_region.center.y, bkg_region.meta["text"], color = "white", size = text_fontsize)
            axs[0,1].set_title("Data Image", fontsize = subtitle_fontsize)
    else:
        for src_region in src_pixel_regions: # src region is the second option
            src_region.plot(ax = axs[0,1], lw=1.0, label = src_region.meta["text"], color = np.random.rand(3,))
            axs[0,1].text(src_region.center.x, src_region.center.y, src_region.meta["text"], color = "white", size = text_fontsize)
            axs[0,1].set_title("Data Image", fontsize = subtitle_fontsize)
        
    # start plotting the regions
    for src_region, bkg_region, idx in zip(src_pixel_regions, bkg_pixel_regions, np.arange(1, n_subplots)):
        
        if src_region is not None:
            axs[idx,0].imshow(image_array_data, origin='lower', norm = ImageNormalize(data = image_array_data, stretch = LogStretch(data_strentch)), cmap='Greys_r', interpolation='nearest')
            axs[idx,0].scatter(src_region.center.x, src_region.center.y, marker = "+", s = marker_size, color = "lime", label = f"{src_region.meta['text']} source region center")
            src_region.plot(ax = axs[idx,0], color='lime', lw=1.0, label =  f"{src_region.meta['text']} source region")
            axs[idx,0].set_xlim(src_region.center.x - image_array_data.shape[0]*0.05, src_region.center.x + image_array_data.shape[0]*0.05)
            axs[idx,0].set_ylim(src_region.center.y - image_array_data.shape[1]*0.05, src_region.center.y + image_array_data.shape[1]*0.05)
            axs[idx,0].set_title(f"Source Region for {src_region.meta['text']}", fontsize = subtitle_fontsize)
        
        if bkg_region is not None:
            axs[idx,1].imshow(image_array_data, origin='lower', norm = ImageNormalize(data = image_array_data, stretch = LogStretch(data_strentch)), cmap='Greys_r', interpolation='nearest')
            bkg_region.plot(ax = axs[idx,1], color='red', lw=1.0, label =  f"{bkg_region.meta['text']} bakcground annulus")
            axs[idx,1].set_xlim(bkg_region.center.x - image_array_data.shape[0]*0.05, bkg_region.center.x + image_array_data.shape[0]*0.05)
            axs[idx,1].set_ylim(bkg_region.center.y - image_array_data.shape[1]*0.05, bkg_region.center.y + image_array_data.shape[1]*0.05)
            axs[idx,1].set_title(f"Background Region for {bkg_region.meta['text']}", fontsize = subtitle_fontsize)
            
    # create RA Dec grids
    for ax in axs.flatten():
        ax.coords.grid(True, color='white', ls='dotted')
        ax.coords[0].set_axislabel('Right Ascension (J2000)', fontsize=tick_fontsize)
        ax.coords[1].set_axislabel('Declination (J2000)', fontsize=tick_fontsize)
        ax.tick_params(which='major', labelsize=tick_fontsize)
        ax.legend(fontsize = subplot_labelsize)
        ax.legend()
        
    if save_image is True:
        fig.savefig(save_dir / f"{field_name}_{image_filter_name}_regions.png", dpi = 300, bbox_inches = "tight")
        plt.close()
        
    return
        
        
    
                
        
        
