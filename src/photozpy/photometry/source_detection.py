"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function:
- Find local peaks in the image
- Fit the centroids of the local peaks

Legacy peak/centroid utilities awaiting their own refactoring. FWHM measurement
now lives exclusively in photozpy.photometry.fwhm.FWHM; the superseded radial
FWHM implementation and its private clipping/background helpers were removed.
"""

from ..collection_manager import CollectionManager
from pathlib import Path
from photutils.detection import find_peaks
from photutils.centroids import centroid_com, centroid_sources
from photutils.aperture import CircularAperture
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats
from astropy.visualization import simple_norm
import numpy as np
import matplotlib.pyplot as plt


class SourceDetection():

    def __init__(self, image_collection):

        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(
            image_collection)

    @staticmethod
    def find_local_peaks(image_path, std_clip_threshold=3,
                         local_peak_detection_threshold=10, fit_box_size=41, plot=False, verbose=True):
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
        # data, mask = lacosmic(data)
        filter = ccddata.header["FILTER"]
        object = ccddata.header["OBJECT"]
        if not verbose:
            print(
                f"Determining the local peaks of {object} in {filter} filter.....")

        # find the local peak positions
        mean, median, std = sigma_clipped_stats(
            data, sigma=std_clip_threshold)  # clip the extreme values from data
        # calculate the threshold used to identify the peaks
        threshold = median + (local_peak_detection_threshold * std)
        # Qtable, it's pandas-like object with columns x_peak, y_peak, and
        # peak_value
        tbl = find_peaks(data, threshold, box_size=fit_box_size)
        # the array of pixel coordinates of the local peaks
        local_peak_positions = np.transpose((tbl['x_peak'], tbl['y_peak']))

        return local_peak_positions

    def find_centroids(self, method="median", image_type="Master Light",
                       std_clip_threshold=3, local_peak_detection_threshold=10, fit_box_size=41, centroid_function=centroid_com,
                       show_whole_image=False, show_individual_sources=False, verbose=False):

        # refresh the full collection to make sure it includes all the newly
        # generated images (Master Light)
        self._image_collection = CollectionManager.refresh_collection(
            self._image_collection, rescan=True)

        # Filter the image collection, only keep Master Light type
        image_collection = CollectionManager.filter_collection(
            self._image_collection, **{"IMTYPE": image_type})

        # get the image list to iterate through
        image_list = image_collection.files_filtered(include_path=True)
        centroids_dict = {"image_path": [], "centroids": []}

        for image_path in image_list:
            # read the image data
            image_path = Path(image_path)
            ccddata = CCDData.read(image_path)
            data = ccddata.data
            filter = ccddata.header["FILTER"]
            object = ccddata.header["OBJECT"]

            # detect local peaks
            local_peak_positions_ = SourceDetection.find_local_peaks(image_path,
                                                                     std_clip_threshold=std_clip_threshold,
                                                                     local_peak_detection_threshold=local_peak_detection_threshold,
                                                                     fit_box_size=fit_box_size,
                                                                     plot=False,
                                                                     verbose=verbose)
            # fit and find the centroids using the local peak positions
            x_local_peak_positions = local_peak_positions_[:, 0]
            y_local_peak_positions = local_peak_positions_[:, 1]
            x_centroids, y_centroids = centroid_sources(
                data, x_local_peak_positions, y_local_peak_positions, box_size=fit_box_size, centroid_func=centroid_function)
            centroids = np.transpose(np.array([x_centroids, y_centroids]))
            if not verbose:
                print(
                    f"Determining the centroids of {object} in {filter} filter.....")

            centroids_dict["image_path"].append(image_path)
            centroids_dict["centroids"].append(centroids)

            if show_whole_image:
                apertures = CircularAperture(centroids, r=15)
                norm = simple_norm(data, 'sqrt', percent=99.9)
                plt.imshow(
                    data,
                    cmap='Greys_r',
                    origin='lower',
                    norm=norm,
                    interpolation='nearest')
                apertures.plot(color='#0547f9', lw=0.5)
                plt.xlim(0, data.shape[1] - 1)
                plt.ylim(0, data.shape[0] - 1)
                # change the suffix of the image to png
                save_path = image_path.with_suffix(".png")
                new_stem = f"{save_path.stem}_centroids"
                save_path = save_path.with_stem(new_stem)
                # print(save_path)
                plt.savefig(save_path, dpi=400, bbox_inches="tight")
                plt.clf()  # clear canvas

            if show_individual_sources:
                for xcen, ycen in zip(x_centroids, y_centroids):
                    centroid = np.transpose(np.array([xcen, ycen]))
                    apertures = CircularAperture(centroid, r=15)
                    norm = simple_norm(data, 'sqrt', percent=99.9)
                    plt.imshow(
                        data,
                        cmap='Greys_r',
                        origin='lower',
                        norm=norm,
                        interpolation='nearest')
                    apertures.plot(color='#0547f9', lw=1.5)
                    plt.xlim(xcen - 30, xcen + 30)
                    plt.ylim(ycen - 30, ycen + 30)
                    # change the suffix of the image to png
                    save_path = image_path.with_suffix(".png")
                    save_path.with_stem(
                        f"{save_path.stem}_centroid_at_{xcen}_{ycen}")
                    plt.savefig(save_path, dpi=400, bbox_inches="tight")
                    plt.clf()  # clear canvas

            print("----------------------------------------\n")

        return centroids_dict
