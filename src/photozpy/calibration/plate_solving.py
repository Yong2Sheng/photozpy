"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- apply the calibrations (master bias, master dark and master flat)
"""

from astroquery.astrometry_net import AstrometryNet
from astropy.wcs import WCS
from ..collection_manager import CollectionManager
from pathlib import Path
from astropy.io import fits

class PlateSolving():

    def __init__(self, image_collection, sources):

        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(image_collection, rescan = True)
        self._sources = sources


    def add_wcs(self, api, image_type = "Master Light", fwhm = None, detect_threshold = 5, solve_timeout = 120):

        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(self._image_collection, rescan = True)

        image_collection = CollectionManager.filter_collection(self._image_collection, **{"IMTYPE": image_type})
        all_image_list = image_collection.files_filtered(include_path=True)
        
        # here I want to move the standard sources to the front of the list
        # So once one standard image is finished adding wcs, I can start deciding the standard stars and the standard magnitudes RIGHT AWAY!
        # I don't have to wait until the end of plate solving since it's quiet time consuming.
        std_collection = CollectionManager.filter_collection(self._image_collection, **{"IMTYPE": "Master Light", "OBJECT":  self._sources.standard_stars[0].source_name})
        std_list = std_collection.files_filtered(include_path=True)
        target_list = [i for i in all_image_list if not i in std_list]
        image_list = std_list + target_list
        
        for image in image_list:
            print(f"Adding WCS to {image}")
            if fwhm == None:
                # try to read the fwhm from the header
                with fits.open(image) as hdul:
                    fwhm = hdul[0].header["FWHM"]
            image_file = Path(image).name
            ast = AstrometryNet()
            ast.api_key = api
            wcs_header = ast.solve_from_image(image, fwhm = fwhm, detect_threshold=detect_threshold, solve_timeout = solve_timeout)#, verbose = False)
            plate_solved_wcs_header = WCS(wcs_header).to_header()
            
            with fits.open(image, mode = "update") as hdul:
                for card in plate_solved_wcs_header.cards:
                    hdul[0].header[card[0]] = (card[1], "UPDATED!!")
                hdul.flush()
            print("-----------------------------------------------------------------------------------\n")

        return