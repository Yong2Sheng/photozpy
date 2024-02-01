"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- Find local peaks in the image
- Fit the centroids of the local peaks
"""

from ..collection_manager import CollectionManager
from pathlib import Path
from photutils.detection import find_peaks
from photutils.profiles import RadialProfile
from photutils.centroids import centroid_quadratic, centroid_com, centroid_sources, centroid_2dg
from photutils.aperture import CircularAperture
from astropy.nddata import CCDData
from astropy.stats import mad_std, sigma_clipped_stats, median_absolute_deviation
from astropy.visualization import hist, simple_norm
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from astropy.io import fits
from ccdproc import cosmicray_lacosmic as lacosmic


class SourceDetection():

    def __init__(self, image_collection):

        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(image_collection)


    @staticmethod
    def find_local_peaks(image_path, std_clip_threshold = 3, local_peak_detection_threshold = 10, fit_box_size = 41, plot = False, verbose = True):

        """
        Find the local peak positions (pixel coordinates) of an image.

        Parameters
        ----------
        image_path; str; Path, the path to the image to be measured
        plot: boolean, plot the local peaks or not
        verbose: boolean, determines the print of some hints and progress

        Returns
        -------
        local_peak_positions: np.ndarray; the array of the local peak positions.
        """
        
        # read data and headers
        image_path = Path(image_path)
        ccddata = CCDData.read(image_path)
        data = ccddata.data
        #data, mask = lacosmic(data)
        filter = ccddata.header["FILTER"]
        object = ccddata.header["OBJECT"]
        if verbose == False:
            print(f"Determining the local peaks of {object} in {filter} filter.....")
            
        # find the local peak positions
        mean, median, std = sigma_clipped_stats(data, sigma = std_clip_threshold)  # clip the extreme values from data
        threshold = median + (local_peak_detection_threshold * std)  # calculate the threshold used to identify the peaks
        tbl = find_peaks(data, threshold, box_size = fit_box_size)  # Qtable, it's pandas-like object with columns x_peak, y_peak, and peak_value
        local_peak_positions = np.transpose((tbl['x_peak'], tbl['y_peak']))  # the array of pixel coordinates of the local peaks

        return local_peak_positions
        

    def find_centroids(self, method = "median", image_type = "Master Light", 
                       std_clip_threshold = 3, local_peak_detection_threshold = 10, fit_box_size = 41, centroid_function = centroid_com, 
                       show_whole_image = False, show_individual_sources = False, verbose = False):

        
        # refresh the full collection to make sure it includes all the newly generated images (Master Light)
        self._image_collection = CollectionManager.refresh_collection(self._image_collection, rescan = True)

        # Filter the image collection, only keep Master Light type
        image_collection = CollectionManager.filter_collection(self._image_collection, **{"IMTYPE": image_type})

        # get the image list to iterate through
        image_list = image_collection.files_filtered(include_path = True)
        centroids_dict = {"image_path":[],"centroids":[]}
        
        for image_path in image_list:
            # read the image data
            image_path = Path(image_path)
            ccddata = CCDData.read(image_path)
            data = ccddata.data
            filter = ccddata.header["FILTER"]
            object = ccddata.header["OBJECT"]

            # detect local peaks
            local_peak_positions_ = SourceDetection.find_local_peaks(image_path, 
                                                                     std_clip_threshold = std_clip_threshold, 
                                                                     local_peak_detection_threshold = local_peak_detection_threshold, 
                                                                     fit_box_size = fit_box_size, 
                                                                     plot = False, 
                                                                     verbose = verbose)
            # fit and find the centroids using the local peak positions
            x_local_peak_positions = local_peak_positions_[:,0]
            y_local_peak_positions = local_peak_positions_[:,1]
            x_centroids , y_centroids = centroid_sources(data, x_local_peak_positions, y_local_peak_positions, box_size = fit_box_size, centroid_func = centroid_function)
            centroids = np.transpose(np.array([x_centroids, y_centroids]))
            if verbose == False:
                print(f"Determining the centroids of {object} in {filter} filter.....")
                
            centroids_dict["image_path"].append(image_path)
            centroids_dict["centroids"].append(centroids)

            if show_whole_image == True:
                apertures = CircularAperture(centroids, r=15)
                norm = simple_norm(data, 'sqrt', percent=99.9)
                plt.imshow(data, cmap='Greys_r', origin='lower', norm=norm, interpolation='nearest')
                apertures.plot(color='#0547f9', lw=0.5)
                plt.xlim(0, data.shape[1] - 1)
                plt.ylim(0, data.shape[0] - 1)
                save_path = image_path.with_suffix(".png")  # change the suffix of the image to png
                new_stem = f"{save_path.stem}_centroids"
                save_path = save_path.with_stem(new_stem)
                #print(save_path)
                plt.savefig(save_path, dpi = 400, bbox_inches = "tight")
                plt.clf()  # clear canvas

            if show_individual_sources == True:
                for xcen, ycen in zip(x_centroids, y_centroids):
                    centroid = np.tranpose(np.array([xcen, ycen]))
                    apertures = CircularAperture(centroid, r=15)
                    norm = simple_norm(data, 'sqrt', percent=99.9)
                    plt.imshow(data, cmap='Greys_r', origin='lower', norm=norm, interpolation='nearest')
                    apertures.plot(color='#0547f9', lw=1.5)
                    plt.xlim(xcen-30, xcen+30)
                    plt.ylim(ycen-30, ycen+30)
                    save_path = image_path.with_suffix(".png")  # change the suffix of the image to png
                    save_path.with_stem(f"{save_path.stem}_centroid_at_{xcen}_{ycen}")
                    plt.savefig(save_path, dpi = 400, bbox_inches = "tight")
                    plt.clf()  # clear canvas

            print("----------------------------------------\n")

        return centroids_dict


    @staticmethod
    def get_radial_bkg(radial_profile):

        """
        Estimate the background from the radial profile.

        Parameters
        ----------
        radial_profile: photutils.profiles.radial_profile.RadialProfile; the radial profile object for a source

        Returns
        -------
        bkg: float; the estimated background
        """

        std = radial_profile.profile.std()  # standard deviation
        boo = radial_profile.profile < std  # get the boolean list that picks up the radius larger than 1 sigma
        if sum(boo) == 0:
            bkg = None
        else:
            sigma3 = radial_profile.radius[boo][0]*3  # get the radius of 3 sigma
            bkgs = radial_profile.profile[radial_profile.radius>sigma3]
            if bkgs.shape[0] == 0:
                bkg = None
            else:
                bkg = np.mean(radial_profile.profile[radial_profile.radius>sigma3])
        
        return bkg
        
    @staticmethod
    def find_fwhm_from_image(image, xycen, edge_radii_end = 26, edge_radii_step = 0.5):

        """
        Find the FWHM at the xy centroids in a image

        Parameters
        ----------
        image: str, path.Path; the path of the image
        xycen: numpy.ndarray; np.array([xcen, ycen])
        edge_radii_end: int; the radius of the circle used to find the radial profile
        edge_radii_step: float; the step of each radius step

        returns
        -------
        fwhm; float; the measured fwhm
        """

        # read the ccddata
        ccddata = CCDData.read(image)
        data = ccddata.data
        
        edge_radii = np.arange(0, edge_radii_end, edge_radii_step)
        rp = RadialProfile(data, xycen, edge_radii, mask=None)

        # calculate bkg
        bkg = SourceDetection.get_radial_bkg(rp)

        if bkg == None:
            fwhm = 99
        else:
    
            # get half max
            half_max = (rp.profile.max() + bkg)/2
    
            # get FWHM
            interp_f = interpolate.interp1d(rp.radius, rp.profile, kind="linear")
            for i in rp.radius:
                _value = interp_f(i)
                difference = _value - half_max
                if difference <= 0:
                    fwhm = 2*i
                    break
                elif i == rp.radius[-1]:
                    #print(f"The fit of FWHM failed at xcen:{xycen[0]}; ycen:{xycen[1]} for image {image}")
                    fwhm = 99  # tag the failed FWHM fit with a large number 
                    break

        return fwhm

    @staticmethod
    def get_mad_clipped_average(data, mad_sigma = 3):

        """
        Get the mad clipped average of the input data

        Parameters
        ----------
        data: numpy.ndarray
        mad_sigma: int; the sigma value used to reject the extreme values

        Returns
        -------
        data_clipped_average: float
        """

        median = np.median(data.flatten())
        mad = median_absolute_deviation(data.flatten())
        boolean1 = (data - median)/mad < mad_sigma
        boolean2 = (data - median)/mad > -mad_sigma
        boolean = boolean1*boolean2
        data_clipped = data[boolean]
        data_clipped_average = np.mean(data_clipped).round(4)
        
        return data_clipped_average

    def find_fwhm(self, edge_radii_end = 26, edge_radii_step = 0.5, 
                  std_clip_threshold = 3, local_peak_detection_threshold = 10, fit_box_size = 41, centroid_function = centroid_com, 
                  write_fwhm = False):

        """
        Find the averge FWHM out of the sources detected in an image collection.
        Write the FWHM into the fits header

        Parameters
        ----------
        None

        Returns
        -------
        average_fwhm: float; the averaged fwhm found from the image collection
        """

        # refresh the full collection to make sure it includes all the Master Light images
        self._image_collection = CollectionManager.refresh_collection(self._image_collection, rescan = True)

        fwhm_all = []

        # get the centroid dictionary
        centroid_dictionary = self.find_centroids(std_clip_threshold = std_clip_threshold, 
                                                  local_peak_detection_threshold = local_peak_detection_threshold, 
                                                  fit_box_size = fit_box_size, 
                                                  centroid_function = centroid_com)

        for image_path, centroids in zip(centroid_dictionary["image_path"], centroid_dictionary["centroids"]):

            image_name = image_path.name
            print(f"Finding FWHM for {image_name}")

            # read the ccddata
            ccddata = CCDData.read(image_path)
            data = ccddata.data
            
            edge_radii = np.arange(0, edge_radii_end, edge_radii_step)

            for centroid in centroids:

                fwhm = SourceDetection.find_fwhm_from_image(image_path, centroid, 
                                                            edge_radii_end = edge_radii_end, edge_radii_step = edge_radii_step)
                fwhm_all += [fwhm]

        fwhm_all = np.array(fwhm_all)
        fwhm_all = fwhm_all[fwhm_all<99]  # remove the failed FWHM fits
        averaged_fwhm = SourceDetection.get_mad_clipped_average(fwhm_all)
        print(f"The FWHM of the image collection is {averaged_fwhm}")
        print("----------------------------------------\n")

        if write_fwhm == True:

            for image_path in centroid_dictionary["image_path"]:
                with fits.open(image_path, mode = "update") as hdul:
                    hdul[0].header["FWHM"] = averaged_fwhm
                    hdul.flush()

        return

        
        

        
        

        

        

