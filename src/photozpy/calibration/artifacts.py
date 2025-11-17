"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- Remove the cosmic ray and other ccd artifacts 
"""
from ..collection_manager import CollectionManager
from pathlib import Path
from astropy.nddata import CCDData
from ccdproc import cosmicray_lacosmic as lacosmic
from ..convenience_functions import plot_image


class Artifacts():

    def __init__(self, image_collection):

        """
        Define the image collection input.

        Parameters
        ----------
        image_collection

        Returns
        -------
        ?
        """

        self._image_collection = CollectionManager.refresh_collection(image_collection, rescan = True)


    @staticmethod
    def _remove_cosmic_ray(fits_path, plot_before_after = False, plot_mask = False, norm_percent = 99.9, save_location = None, **kwargs):

        """
        Remove cosmic rays (and other possible artifacts) from a science image. The science image must bias, dark and flat corrected. This function doesn't work on calibration images (bias, dark and flat)

        Parameters
        ----------
        fits_path : str or pathlib.Path
        plot_before_after : bool, default=False
            Set True to plot the images before and after the cosmic ray removal.
        plot_mask : bool, default=False
            Set True to plot the mask.
        kwargs : dict
            The paramters and arguments passing to `lacosmic`.

        Returns
        -------
        None
        """

        # read ccd data
        fits_path = Path(fits_path)
        ccd = CCDData.read(fits_path)

        if plot_before_after:
            plot_image(fits_path = fits_path, norm_percent = norm_percent, fname_append = "_with_cosmic_ray")

        # define some inputs for lacosmic, which depends on the header values
        gain = ccd.header["GAIN"]  # unit: electrom/adu
        rdnoise = ccd.header["RDNOISE"]  # unit: electron

        newccd = lacosmic(ccd, gain_apply = False, gain = gain, readnoise = rdnoise, **kwargs)
        hdu = newccd.to_hdu(hdu_uncertainty = None)

        if plot_before_after:
            file_name = fits_path.stem + "_without_cosmic_ray" + ".png"
            plot_image(ccd_data = newccd, norm_percent = norm_percent, save_location = fits_path.parent, fname_append = file_name)

        # create hdu object
        hdu = newccd.to_hdu(hdu_uncertainty = None)
        hdu[0].header["BUNIT"] = "adu"

        # save file
        if save_location == None:
            save_path = fits_path
        else:
            save_path = Path(save_location)/fits_path.name
            
        hdu.writeto(save_path, overwrite = True)

        return

    def remove_cosmic_ray(self, plot_before_after = False, plot_mask = False, norm_percent = 99.9, save_location = None, headers_values = None, **kwargs):

        """
        Remove cosmic ray for the images in the collection.

        Parameters
        ----------

        Returns
        -------
        """

        # refresh the image collection to include newly created files
        self._image_collection = CollectionManager.refresh_collection(self._image_collection, rescan = True)

        if headers_values != None:
            image_collection = CollectionManager.filter_collection(self._image_collection, **headers_values)
        else:
            image_collection = self._image_collection

        fits_paths = image_collection.files_filtered(include_path = True)

        for fits_path in fits_paths:
            fits_path = Path(fits_path)
            fname = fits_path.name
            print(f"Removing cosmic ray from {fname}....")
            Artifacts._remove_cosmic_ray(fits_path, plot_before_after = plot_before_after, plot_mask = plot_mask, 
                                        norm_percent = norm_percent, save_location = save_location, **kwargs)

        print("Cosmic removal completed!")
        print("----------------------------------------\n")

        return

            
        

        
        
        