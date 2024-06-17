"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- Perform aperture photometry
"""
from pathlib import Path
from ..collection_manager import CollectionManager
from astropy.nddata import CCDData
from photutils.aperture import CircularAnnulus, CircularAperture, ApertureStats, aperture_photometry
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import SigmaClip
import numpy as np
from ..convenience_functions import *
import re
from ..mimage_collection import mImageFileCollection
from ccdproc import ImageFileCollection
import pandas as pd
import os

class Photometry():

    def __init__(self, image_collection):
        
        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(image_collection, rescan = True)

    # @staticmethod
    # def convert_coords(image_path = None, wcs = None, skycoords = None, pixelcoords = None, verbose = False):

    #     """
    #     Takes a fits file or the astropy.wcs.WCS object as input. 
    #     Convert the sky coordinates to the pixel coordinates or vice versa.

    #     Parameters
    #     ----------
    #     file: str or path.Path object; the path to the object
    #     wcs: the astropy.wcs.WCS object. If both were input, wcs will cover the file.
    #     skycoords: astropy.Skycoord object. The sky coordinate of the object
    #     pixelcoords: 2D numpy array; the pixel coordinate of the object, the first column is the x pixels and the second is the y pixels: [[x pixles],[y pixels]]

    #     Returns
    #     -------
    #     astropy.Skycoords or list
    #     """

    #     # check if the number of input satistifies the calculation

    #     if image_path == None and wcs == None:
    #         raise TypeError("You must give a file path or asrtropy.wcs.WCS obkect as the input!")

    #     if skycoords == None and pixelcoords == None:
    #         raise TypeError("You must give sky coordinates or pixel coordinates as the input!")

    #     elif skycoords != None and pixelcoords != None:
    #         raise TypeError("Please only input the sky coordinates or the pixel coordinates!")

    #     # Read the file and get the wcs object
    #     if wcs != None:
    #         wcs_object = wcs
    #     else:
    #         data = CCDData.read(image_path, unit = "adu")
    #         wcs_object = data.wcs

    #     if skycoords != None and pixelcoords == None:
    #         pixelcoords = data.wcs.world_to_pixel(skycoords)
    #         pixelcoords = np.array((pixelcoords)).T # transfer the array so it's ra/dec in each column
    #         out = pixelcoords # output variable
    #         if verbose == True:
    #             print("Conversion from sky coordiantes to pixel coordinates completed!")

    #     elif skycoord == None and pixelcoord != None:
    #         xpixel_coords = pixelcoord[0]
    #         ypixel_coords = pixelcoord[1]
    #         radec = data.wcs.pixel_to_world(xpixels, ypixel_coords)
    #         out = radec
    #         if verbose == True:
    #             print("Conversion from pixel coordinates to sky coordinates completed!")

    #     return out

    @staticmethod
    def read_fwhm(image_path, keyword = "FWHM"):

        """
        Read the FHWM from the header of the fits file

        Parameters
        ----------
        image_path: string or pathlib.Path; the path to the image
        keyword: string; the keyword name in the header.

        Returns
        -------
        fwhm: float; the read fwhm
        """

        headers = fits.getheader(image_path)
        fwhm = headers[keyword]

        return fwhm
        
    @staticmethod
    def get_aper_centroid(image_path, xy_coords, fwhm):

        """
        Estimate the centroid from aperture.

        Parameters
        Parameters
        ----------
        image_path: string or pathlib.Path; the path to the image
        xycen: numpy.ndarray or list; the physical coordinates of the source: each row: [x, y].

        Returns
        -------
        xycentroid: list; the x, y centroid position
        """
        # get the sigma_clipped median bkg
        data = CCDData.read(image_path)
        
        # fit the centroids first
        aper = CircularAperture(xy_coords, fwhm)
        aperstats = ApertureStats(data, aper)
        xycentroid = np.array([aperstats.xcentroid, aperstats.ycentroid]).T

        return xycentroid
        
    
    @staticmethod
    def get_background(image_path, annulus_aperture, sigma = 3, fwhm = None):

        """
        Estimate the background within the annulus using sigma-clipped median.

        Parameters
        ----------
        image_path: string or pathlib.Path; the path to the image
        xycen: numpy.ndarray or list; the coordinates of the source: [x, y].

        Returns
        -------
        bkg
        """

        if fwhm == None:
            # get the fwhm
            fwhm = Photometry.read_fwhm(image_path)

        # setup sigma clip
        sigclip = SigmaClip(sigma = sigma, maxiters = 10)

        # get the sigma_clipped median bkg
        data = CCDData.read(image_path)
        bkg_stats = ApertureStats(data, annulus_aperture, sigma_clip=sigclip)

        return bkg_stats.median

            
    @staticmethod
    def counts2mag(counts, c_const = 0):
        """
        Calculate the instrumental magnitude.

        Paremeters
        ----------
        counts: float;
        c_const; float; the zero point

        Returns
        -------
        mag: float; magnitude
        """
        mag = -2.5*np.log10(counts.to_value()) + c_const

        return mag
    
    def aperture_photometry(self, sources, bkg_clip_sigma = 3, verbose = True, fhwm_aper_factor = 3.0):

        """
        Does the aperture photometry on the skycoords.

        Parameters
        ----------
        sources: photozpy.photometry.sources.Sources; the source dictionary
        aperture: float; the source aperture
        inner_annulus: float; the inner radius of the annulus
        outer_annulus: float; the outer radius of the annulus

        Returns
        -------
        None
        """

        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(self._image_collection, rescan = True)

        # work on the obejct iteratively
        for object, skycoords in zip(sources.get_objects, sources.get_skycoords):

            # get the image collection to work
            collection_photometry = CollectionManager.filter_collection(self._image_collection, 
                                                                        **{"IMTYPE": "Master Light", "OBJECT": object})
            image_list = collection_photometry.files_filtered(include_path = True)

            for image_path in image_list:
                headers = fits.getheader(image_path)
                filter = headers["FILTER"]
                print(f"Working on photometry of {object} in {filter} ......")

                # get the aperture and annulus aperture
                fwhm = Photometry.read_fwhm(image_path, keyword = "FWHM")
                xy_coords = convert_coords(image_path = image_path, skycoords = skycoords)
                xycentroids = Photometry.get_aper_centroid(image_path = image_path, xy_coords = xy_coords, fwhm = fwhm)
                apertures = CircularAperture(xycentroids, r=fhwm_aper_factor*fwhm)
                annlus_apertures = CircularAnnulus(xycentroids, r_in=5*fwhm, r_out=8*fwhm)

                # get the sigma_clipped background estimation for all the annulus apertures
                bkgs = Photometry.get_background(image_path, annlus_apertures, sigma = bkg_clip_sigma, fwhm = None)

                # perform aperture photometry
                data = CCDData.read(image_path)
                phot_table = aperture_photometry(data, apertures)

                # substract the background from the photometry
                total_bkgs = bkgs * apertures.area
                phot_bkgsub = phot_table['aperture_sum'] - total_bkgs

                # calculate the instrumental magnitude
                m_inst = Photometry.counts2mag(phot_bkgsub)

                # organize the Qtable
                phot_table['total_bkg'] = total_bkgs  # add the column for total background
                phot_table['aperture_sum_bkgsub'] = phot_bkgsub  # add the column for bkg substracted photometry
                phot_table['mag_inst'] = m_inst  # add the column for instrumental magnitude
                
                
                for col in phot_table.colnames:
                    phot_table[col].info.format = '%.8g'  # for consistent table output

                if verbose:
                    print(phot_table)
                    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

                # plot the targets, photometry aperture and background annulus
                _ = plot_image(fits_path = image_path, 
                               save_location = None,
                               skycoords = skycoords, 
                               pixelcoords = None, 
                               circular_apertures = apertures, 
                               annulus_apertures = annlus_apertures,
                               adjust_fov = True,
                               fname_append = "_apertures")


        return

class SwiftPhotometry():

    def __init__(self, image_collection, source_catalog = None):

        if isinstance(image_collection, ImageFileCollection):
            self._mcollection = mImageFileCollection(location = image_collection.location, filenames = image_collection.files)

        elif isinstance(image_collection, mImageFileCollection):
            self._mcollection = image_collection

        self._source_catalog = source_catalog

    @staticmethod
    def extract_mag(out_txt):
    
        """
        Extract the magnitude information from the output txt file
        """
    
        f = open(out_txt, "r")
        lines = f.readlines()
        f.close()
    
        for idx, line in enumerate(lines):
            if "AB system" in line:
                mag_idx = idx + 1
                mag_line = lines[mag_idx]
                if ">" in mag_line:
                    mag = float(re.findall("\d+\.?\d*", mag_line)[0])  # One or more digits (\d+), optional period (\.?), zero or more digits (\d*).
                    mag_error = -99
                else:
                    number_list = re.findall("\d+\.?\d*", mag_line)  # One or more digits (\d+), optional period (\.?), zero or more digits (\d*).
                    mag = float(number_list[0])
                    stat_err = float(number_list[1])
                    sys_err = float(number_list[2])
                    mag_error = round(np.sqrt(stat_err**2+sys_err**2),2)
    
        return mag, mag_error

    def UVOTPhotometry(self):

        """
        Perform Swift photometry on the UVOT filters.
        """

        # init the dataframe to store the results
        df_results = pd.DataFrame(columns = ['name', "UVW2", "UVW2_err", "UVM2", "UVM2_err", "UVW1", "UVW1_err", "U", "U_err", "B", "B_err", "V", "V_err"])
        
        for collection in self._mcollection:

            file_location = Path(collection.location)
            source_name_from_path = file_location.parts[-1].replace("_", " ")
            file_paths = collection.files_filtered(include_path = True, **{"SUMTYP": "FINAL"})  # only use the final summed fits files

            # init the dict to store the photometry for a single source, which will be appended to the main result data frame
            dict_new = {'name' : source_name_from_path,
                        "UVW2" : -99, 
                        "UVW2_err" : -99,
                        "UVM2" : -99,
                        "UVM2_err" : -99,
                        "UVW1" : -99,
                        "UVW1_err" : -99,
                        "U" : -99,
                        "U_err" : -99,
                        "B" : -99,
                        "B_err" : -99,
                        "V" : -99,
                        "V_err" : -99}
            
            for file_path in file_paths:

                # set up the file path, region path, and read the source name from the path
                file_path = Path(file_path)  # file path
                headers = fits.getheader(file_path, ext = 1)
                filter_name = headers["FILTER"]   # filter name
                source_region_path = file_path.parent / f"{filter_name}.reg"  # source region path
                bkg_region_path = file_path.parent / f"{filter_name}_bkg.reg"  # bkg region path
                source_name_from_fits = headers["OBJECT"]

                # check if the two source names match
                if source_name_from_path == source_name_from_fits:
                    pass
                else:
                    raise ValueError(f"The source name from the path ({source_name_from_path}) doesn't match the one from the fits file ({source_name_from_fits})")
                
                print(f"Working on the photometry of source {source_name_from_fits} in {filter_name} filter.")
                
                outfits_path = file_path.parent / f"{filter_name}_result.fits"
                outtxt_path = file_path.parent / f"{filter_name}.result.txt"

                # run the command
                print(f"Running uvotsource image={file_path} srcreg={source_region_path} bkgreg={bkg_region_path} sigma=3 cleanup=y clobber=y outfile={outfits_path} | tee {outtxt_path} >/dev/null")
                os.system(f"uvotsource image={file_path} srcreg={source_region_path} bkgreg={bkg_region_path} sigma=3 cleanup=y clobber=y outfile={outfits_path} | tee {outtxt_path} >/dev/null")

                # extract the magnitudes
                ab_mag, ab_err = SwiftPhotometry.extract_mag(outtxt_path)

                # write the magnitudes to the image file
                with fits.open(file_path, mode = "update") as hdul:
                    hdul[0].header["AB_MAG"] = ab_mag
                    hdul[1].header["AB_MAG_ERR"] = ab_err
                    hdul[0].header["AB_MAG"] = ab_mag,
                    hdul[1].header["AB_MAG_ERR"] = ab_err
                    hdul.flush()

                # write the magnitudes to the dict

                dict_new[filter_name] = ab_mag
                dict_new[f"{filter_name}_err"] = ab_err
                
                print("---------------------------------------------------------------------------")
        
            
            df_new = pd.DataFrame(dict_new, index=[0])
            df_results = pd.concat([df_results, df_new], ignore_index=True)

        if self._source_catalog is not None:
            self._source_catalog = Path(self._source_catalog)
            final_catalog = pd.read_csv(self._source_catalog, sep = ",")
            final_catalog = pd.merge(final_catalog, df_results, on='name')
            final_catalog.to_csv(Path("") / "mag_results.csv", index=False, mode = "w")
        else:
            df_results.to_csv(Path("") / "mag_results.csv", index=False, mode = "w")

        return
                
