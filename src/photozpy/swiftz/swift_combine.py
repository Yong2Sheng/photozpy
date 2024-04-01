"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- ASPCORR check
- Combine extensions within one observation file (the same observation ID)
- Combine images with multiple observation files (across multiple observation IDs)
"""

from ..mimage_collection import mImageFileCollection
from ccdproc import ImageFileCollection
from astropy.io import fits
from pathlib import Path
import os

class SwiftCombine():

    def __init__(image_collection, telescope):

        """
        image_collection : ccdproc.image_collection.ImageFileCollection or photozpy.mimage_collection.mImageFileCollection
            The image collection(s).
        telescope : photozpy.telescope.telescope.Telescope
            The swift telescope spec.
        """
            
        # switch image_collection from ImageFileCollection to mImageFileCollection to make it standard for the pipeline
        if isinstance(image_collection, ImageFileCollection):
            self._mcollection = mImageFileCollection(location = image_collection.location, filenames = image_collection.files)

        elif isinstance(image_collection, mImageFileCollection):
            self._mcollection = image_collection

        self._telescope = telescope



    @staticmethod
    def find_multiext_keywords(file, *keywords):
        
        """
        Find the value of the keywords from a multiextension fits file. It deals with multi-extension fits files
        so you will get a dictionary of all the requested keywords in each extension.
        "ext_No" and "EXTNAME" are mandatory in order to tag the origin of the keyword values.

        Parameters
        ----------
        file: str; the directory to the fits file
        kewords: str(s); the keywords you want to query

        Example run
        -----------
        find_mext_keywords(file_path, "ASPCORR", "HDUCLAS1")
        Out：
        {'ext_No': [1, 2],
        'EXTNAME': ['bb649482961I', 'bb649505811I'],
        'ASPCORR': ['DIRECT', 'DIRECT'],
        'HDUCLAS1': ['IMAGE', 'IMAGE']}
        """

        with fits.open(file) as hdul:
            ext_nums = len(hdul) - 1

            # initialize the dictionary
            dict_ext = {}
            dict_ext["ext_No"] = [i for i in np.arange(1, ext_nums+1)]
            dict_ext["EXTNAME"] = [None] * ext_nums  # [None, None, None, ...]
            for i in keywords:
                dict_ext[i] = [None] * ext_nums  # generate the keys for the keywords like ext_No and EXTNAME

            
            # read and record the keywords in the dictionary
            n_ = 1
            while n_ <= ext_nums:
                dict_ext["EXTNAME"][n_-1] = hdul[n_].header["EXTNAME"]
                for i in keywords:
                    dict_ext[i][n_-1] = hdul[n_].header[i]
                n_ += 1

        return dict_ext

    @staticmethod
    def check_ASPCORR(file):
        """
        Check the ASPCORR keywords.

        Parameters
        ----------
        file: str; the fits file to be exmined.

        Return
        ------
        Boolean
        """
        
        dict_exts = SwiftCombine.find_mext_keywords(file, "ASPCORR")
        asps = dict_exts["ASPCORR"]
        check = ["Yes" for i in asps if "DIRECT" in i]
        if len(check) == len(asps):
            return True
        elif len(check) != len(asps):
            return False

    def sum_extensions(self, fits_path, out_path = None):

        """
        Sum the extensions within the same fits file.

        Paremater
        ---------
        fits_file : str or pathlib.Path
            The path to the fits file
        out_path : str, optional
            The file name of the summed file (the default is `None`, which means it will be named by the filter).

        Return
        ------
        str
            The path of the output file.
        """

        fits_path = Path(fits_path)
        fits_file_name = fits_path.name
        
        # first check if the aspect correction is correct
        if not SwiftCombine.check_ASPCORR(fits_path):
            print(f"The apsect correction is not done for {fits_path}.")
            print(f"Skipping {fits_file_name} since the ASPECT isn't corrected properly! It will be deleted.")
            
        else:
            headers = fits.get_header(file_path)
            filter = self._telescope.map_filters(headers["FILTER"])
            if out_path is None:
                out_path = fits_path.parent / f"{filter}.fits"  # if the out_name is not given, the file will be saved as the filter name.
            
            os.system(f"uvotimsum infile={fits_path} outfile={out_path} | tee -a uvotimsum_log.txt >/dev/null 2>&1")

        os.remove(fits_path) 

        print(f"Extensions in {fits_file_name} has been summed!")

        return out_name

    def sum_fits_files(self, image_collection, out_path = None):

        """
        Sum the fits files. 

        Paremeters
        ----------
        image_collection : ccdproc.ImageFileCollection
            The image collection to sum.
        out_path : str, optional
            The file name of the summed file (the default is `None`, which means it will be named by the filter).

        Return
        ------
        str
            The output file name.
        """


        # check if the filters are the same
        collection_filter = image_collection.values("filters", unique = True)  # collection_filter is list even if there is only one value
        if len(collection_filter) != 1:
            raise ValueError("The fits files to sum have more than one filters!")
        else:
            collection_filter = collection_filter[0]

        if not collection_filter in self._telescope.filters:
            raise ValueError("The filter is not in the telescope spec!")

        file_path = image_collection.files_filtered(include_path = True)
        file_names = image_collection.files_filtered(include_path = False)
        
        # determine the output path
        if out_path is None:
            out_path = Path(files[0]).parent / f"{collection_filter}.fits"

        os.system(f"uvotimsum exclude=NONE infile={file_path} outfile={out_path} | tee -a uvotimsum_log.txt >/dev/null 2>&1")

        print(f"{file_names} has been summed to {collection_filter}.fits")

        os.remove(i) for i in file_path

        return out_path


    def sum_all_files(self):

        """
        Sum all the fits files in the multi image collection.

        Returns
        -------
        photozpy.mcollection.mImageFileCollection
            The new collection that contain the summed fits files.
        """

        for collection in self.mcollection:

            # get the target name
            target_dir = Path(collection.location)
            target_name = target_dir.parts[-1]  # ('data', 'B3 0850+443')[1] is 'B3 0850+443'

            # loop over filters
            for filter in self._telecope:
                
                # First sum all the extensions in the same fits file
                files = collection.files_filtered(include_path = True, filter = filter)

                if len(files) == 0:
                    print(f"No {target_name} in {filter}, skipping ......")

                elif len(files) == 1:
                    self.sum_extensions(fits_path = files[0])

                elif len(files) > 1:
                    num = len(files)

                    # sum the extensions for each file
                    summed_observations = []
                    for i in np.arange(num):
                        _summed = self.sum_extensions(fits_path = files[i], out_path = f"{filter}_{i}")
                        summed_observations += [_summed]  # collect the file path of summed observations

                    # sum the observations
                    _collection = ImageFileCollection(location = target_dir, filenames = summed_observations)

            print(f"{target_name} sum completed!")

        return

        

        
            
        

