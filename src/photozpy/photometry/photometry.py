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

class Photometry():

    def __init__(self, image_collection):
        
        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(image_collection, rescan = True)

    @staticmethod
    def convert_coords(image_path = None, wcs = None, skycoords = None, pixelcoords = None):

        """
        Takes a fits file or the astropy.wcs.WCS object as input. 
        Convert the sky coordinates to the pixel coordinates or vice versa.

        Parameters
        ----------
        file: str or path.Path object; the path to the object
        wcs: the astropy.wcs.WCS object. If both were input, wcs will cover the file.
        skycoords: astropy.Skycoord object. The sky coordinate of the object
        pixelcoords: 2D numpy array; the pixel coordinate of the object, the first column is the x pixels and the second is the y pixels: [[x pixles],[y pixels]]

        Returns
        -------
        astropy.Skycoords or list
        """

        # check if the number of input satistifies the calculation

        if image_path == None and wcs == None:
            raise TypeError("You must give a file path or asrtropy.wcs.WCS obkect as the input!")

        if skycoords == None and pixelcoords == None:
            raise TypeError("You must give sky coordinates or pixel coordinates as the input!")

        elif skycoords != None and pixelcoords != None:
            raise TypeError("Please only input the sky coordinates or the pixel coordinates!")

        # Read the file and get the wcs object
        if wcs != None:
            wcs_object = wcs
        else:
            data = CCDData.read(image_path, unit = "adu")
            wcs_object = data.wcs

        if skycoords != None and pixelcoords == None:
            pixelcoords = data.wcs.world_to_pixel(skycoords)
            pixelcoords = np.array((pixelcoords)).T # transfer the array so it's ra/dec in each column
            out = pixelcoords # output variable
            print("Conversion from sky coordiantes to pixel coordinates completed!")

        elif skycoord == None and pixelcoord != None:
            xpixel_coords = pixelcoord[0]
            ypixel_coords = pixelcoord[1]
            radec = data.wcs.pixel_to_world(xpixels, ypixel_coords)
            out = radec
            print("Conversion from pixel coordinates to sky coordinates completed!")

        return out

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
    
    def aperture_photometry(self, sources, bkg_clip_sigma = 3, verbose = True):

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
                xy_coords = Photometry.convert_coords(image_path = image_path, skycoords = skycoords)
                xycentroids = Photometry.get_aper_centroid(image_path = image_path, xy_coords = xy_coords, fwhm = fwhm)
                apertures = CircularAperture(xycentroids, r=fwhm)
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
                    print("----------------------------------------------------------\n")


        return