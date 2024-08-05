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
from photutils.aperture import CircularAnnulus, CircularAperture, ApertureStats, aperture_photometry, SkyCircularAperture, SkyCircularAnnulus
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
from astropy.coordinates import SkyCoord, concatenate
import astropy.units as u
from regions import PixCoord, CirclePixelRegion, CircleSkyRegion, Regions, CircleAnnulusSkyRegion

class Photometry():

    def __init__(self, image_collection):
        
        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(image_collection, rescan = True)


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
    def get_background(image_array_data, annulus_aperture, sigma = 3):

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

        # setup sigma clip
        sigclip = SigmaClip(sigma = sigma, maxiters = 10)

        # get the sigma_clipped median bkg
        bkg_stats = ApertureStats(image_array_data, annulus_aperture, sigma_clip=sigclip)

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
        mag = -2.5*np.log10(counts) + c_const

        return mag
    
    @staticmethod
    def region2aperture(regions):
    
        if isinstance(regions[0], CircleSkyRegion):
            if len(regions) > 1:
                centers = concatenate([region.center for region in regions])
            else:
                centers = regions[0].center # if the number of regions is one, you can't concatenate the skycoords. Just get the region center
            return SkyCircularAperture(centers, regions[0].radius)
            
        elif isinstance(regions[0], CircleAnnulusSkyRegion):
            if len(regions) > 1:
                centers = concatenate([region.center for region in regions])
            else:
                centers = regions[0].center # if the number of regions is one, you can't concatenate the skycoords. Just get the region cente
            return SkyCircularAnnulus(centers, r_in = regions[0].inner_radius, r_out = regions[0].outer_radius)
                
        
    
    def run_photometry(self, sources, bkg_clip_sigma = 3, hdu = 0, verbose = True):

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
        
        tables = []

        # work on the obejct iteratively
        for source_name, source_coords in sources:

            # get the image collection to work
            collection_photometry = CollectionManager.filter_collection(self._image_collection, 
                                                                        **{"IMTYPE": "Master Light", "OBJECT": source_name})
            image_list = collection_photometry.files_filtered(include_path = True)

            for image_path in image_list:
                image_path = Path(image_path)
                ccddata = CCDData.read(image_path, hdu = hdu)
                image_headers = ccddata.header
                image_filter_name = image_headers["FILTER"]
                image_wcs = ccddata.wcs
                image_array_data = ccddata.data
                print(f"Working on photometry of {source_name} in {image_filter_name} from {image_path.name}")
                print("-------------------------------------------------------------------------------------------------")

                # get the aperture and annulus aperture
            
                src_region_fname = image_path.parent / f"{source_name}_{image_filter_name}_src.reg"
                if not src_region_fname.exists():
                    raise OSError(f"{src_region_fname} not found!")
                else:
                    src_regions = Regions.read(src_region_fname, format='ds9')
                    src_apertures_sky = Photometry.region2aperture(src_regions)
                    src_apertures_pix = src_apertures_sky.to_pixel(image_wcs)
                    
                bkg_region_fname = image_path.parent / f"{source_name}_{image_filter_name}_bkg.reg"
                if not bkg_region_fname.exists():
                    raise OSError(f"{bkg_region_fname} not found!")
                else:
                    bkg_regions = Regions.read(bkg_region_fname, format='ds9')
                    bkg_regions_sky = Photometry.region2aperture(bkg_regions)
                    bkg_annulus_pix = bkg_regions_sky.to_pixel(image_wcs)
                
                # get the sigma_clipped background estimation for all the annulus apertures
                bkgs = Photometry.get_background(image_array_data = image_array_data, annulus_aperture = bkg_annulus_pix, sigma = bkg_clip_sigma)

                # perform aperture photometry
                phot_table = aperture_photometry(ccddata.data, src_apertures_pix)

                # substract the background from the photometry
                total_bkgs = bkgs * src_apertures_pix.area
                phot_table['aperture_sum'].name = "src+bkg"
                phot_bkgsub = phot_table['src+bkg'] - total_bkgs
                phot_bkgsub_error = np.sqrt(phot_table['src+bkg'].value + total_bkgs)

                # calculate the instrumental magnitude
                m_inst = Photometry.counts2mag(phot_bkgsub.value)
                m_inst_error = (2.5/np.log(10))*(phot_bkgsub_error/phot_bkgsub.value)

                # organize the Qtable
                phot_table['bkg'] = total_bkgs  # add the column for total background
                phot_table['src'] = phot_bkgsub  # add the column for bkg substracted photometry
                phot_table['src_error'] = phot_bkgsub_error
                phot_table['mag_inst'] = m_inst  # add the column for instrumental magnitude
                phot_table['mag_inst_error'] = m_inst_error

                phot_table['src+bkg'].unit = u.ct
                phot_table['bkg'].unit = u.ct
                phot_table['src'].unit = u.ct
                phot_table['src_error'].unit = u.ct
                phot_table['mag_inst'].unit = u.mag
                phot_table['mag_inst_error'].unit = u.mag
                phot_table.meta = {"object": source_name,
                                   "filter": image_filter_name}
                
                tables += [phot_table]
                
                for col in phot_table.colnames:
                    phot_table[col].info.format = '%.4f'  # for consistent table output

                if verbose:
                    #nlines = len(phot_table) + 2
                    phot_table.pprint_all()
                    #print(phot_table)
                    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")


        return tables

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
            final_catalog = pd.merge(final_catalog, df_results, on='name', how = "left")
            final_catalog.to_csv(Path("") / "mag_results.csv", index=False, mode = "w")
        else:
            df_results.to_csv(Path("") / "mag_results.csv", index=False, mode = "w")

        return
                
