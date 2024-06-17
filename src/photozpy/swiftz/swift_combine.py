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
import numpy as np
import shutil

class SwiftCombine():

    def __init__(self, image_collection, telescope):

        """
        image_collection : ccdproc.image_collection.ImageFileCollection or photozpy.mimage_collection.mImageFileCollection
            The image collection(s).
        telescope : photozpy.telescope.telescope.Telescope
            The swift telescope spec.
        """
            
        # switch image_collection from ImageFileCollection to mImageFileCollection to make it standard for the pipeline
        if isinstance(image_collection, ImageFileCollection):
            self._mcollection = mImageFileCollection(location = image_collection.location, filenames = image_collection.files)  # this might contain unwanted fits files

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
        find_multiext_keywords(file_path, "ASPCORR", "HDUCLAS1")
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
        
        dict_exts = SwiftCombine.find_multiext_keywords(file, "ASPCORR")
        asps = dict_exts["ASPCORR"]
        check = ["Yes" for i in asps if "DIRECT" in i]
        if len(check) == len(asps):
            return True
        elif len(check) != len(asps):
            return False

    def sum_extensions(self, fits_file_path, out_path = None, delete_files = False, return_full_path = True):

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

        fits_file_path = Path(fits_file_path)
        fits_file_name = fits_file_path.name
        
        # first check if the aspect correction is correct
        if not SwiftCombine.check_ASPCORR(fits_file_path):
            print(f"The apsect correction is not done for {fits_file_path}.")
            print(f"Skipping {fits_file_name} since the ASPECT isn't corrected properly! It will be deleted.")
            
        else:
            headers = fits.getheader(fits_file_path)
            filter_ = headers["FILTER"]
            if out_path is None:
                out_path = fits_file_path.parent / f"{filter_}.fits"  # if the out_name is not given, the file will be saved as the filter name.
            
            print(f"uvotimsum infile={fits_file_path} outfile={out_path} clobber=yes cleanup=yes | tee -a uvotimsum_log.txt >/dev/null 2>&1")
            os.system(f"uvotimsum infile={fits_file_path} outfile={out_path} clobber=yes cleanup=yes | tee -a uvotimsum_log.txt >/dev/null 2>&1")
            out_path_name = out_path.name

        if delete_files is True:
            os.remove(fits_file_path) 

        print(f"Extensions in {fits_file_name} has been summed to {out_path_name}\n")

        if return_full_path is True:
            return out_path
        else:
            return out_path.name

    def sum_fits_files(self, image_collection, out_path = None, delete_files = False):

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
        collection_filter = mImageFileCollection._get_header_values(image_collection, "FILTER", unique = True)  # collection_filter is list even if there is only one value
        if len(collection_filter) != 1:
            raise ValueError("The fits files to sum have more than one filters!")
        else:
            collection_filter = collection_filter[0]

        if not collection_filter in self._telescope.filters:
            raise ValueError("The filter is not in the telescope spec!")

        fits_file_path = image_collection.files_filtered(include_path = True)
        fits_file_names = image_collection.files_filtered(include_path = False)

        # let't append the observations first before using uvotimsum
        appended_file = Path(fits_file_path[0]).with_name(f"appended_{collection_filter}.fits")
        shutil.copy2(fits_file_path[0], appended_file)
        for i in fits_file_path[1:]:
            os.system(f"fappend {i} {appended_file}")

        # determine the output path
        if out_path is None:
            out_path = Path(fits_file_path[0]).parent / f"{collection_filter}.fits"
        
        os.system(f"uvotimsum exclude=NONE infile={appended_file} outfile={out_path} clobber=yes cleanup=yes | tee -a uvotimsum_log.txt >/dev/null 2>&1")

        print(f"{fits_file_names} has been summed to {collection_filter}.fits\n")

        # remove the files that are no longer needed.
        if delete_files is True:
            for i in fits_file_path + [appended_file]:
                os.remove(i)

        return out_path


    def sum_all_files(self, delete_files = False):

        """
        Sum all the fits files in the multi image collection.

        Returns
        -------
        photozpy.mcollection.mImageFileCollection
            The new collection that contain the summed fits files.
        """

        for collection in self._mcollection:

            # get the target name
            target_dir = Path(collection.location)
            target_name = target_dir.parts[-1]  # ('data', 'B3_0850+443')[-1] is 'B3_0850+443'

            # loop over filters
            for filter_ in self._telescope.filters:
                
                # First sum all the extensions in the same fits file
                add_filter = {"FILTER": filter_,
                              "SUMTYP": "NOTSUM"}
                files = collection.files_filtered(include_path = True, **add_filter)

                if len(files) == 0:
                    print(f"No {target_name} in {filter_}, skipping ......")

                elif len(files) == 1:
                    summed_path = self.sum_extensions(fits_file_path = files[0], return_full_path = True)
                    with fits.open(summed_path, mode = "update") as hdul:
                        hdul[0].header["OBJECT"] = target_name.replace("_", " ")
                        hdul[1].header["OBJECT"] = target_name.replace("_", " ")
                        hdul[0].header["SUMTYP"] = "FINAL"
                        hdul[1].header["SUMTYP"] = "FINAL"
                        hdul.flush()

                elif len(files) > 1:
                    num = len(files)

                    # sum the extensions for each file
                    summed_observations = []
                    for i in np.arange(num):
                        _out_path = Path(files[i]).parent / f"{filter_}_{i}.fits"
                        _summed = self.sum_extensions(fits_file_path = files[i], out_path = _out_path, 
                                                      delete_files = delete_files, return_full_path = False)
                        with fits.open(_out_path, mode = "update") as hdul:
                            hdul[0].header["OBJECT"] = target_name.replace("_", " ")
                            hdul[1].header["OBJECT"] = target_name.replace("_", " ")
                            hdul[0].header["SUMTYP"] = "SUMMED"
                            hdul[1].header["SUMTYP"] = "SUMMED"
                            hdul.flush()
                        summed_observations += [_summed]  # collect the file path of summed observations

                    # sum the observations
                    _collection = ImageFileCollection(location = target_dir, filenames = summed_observations)
                    summed_path = self.sum_fits_files(image_collection = _collection, delete_files = delete_files)
                    with fits.open(summed_path, mode = "update") as hdul:
                        hdul[0].header["OBJECT"] = target_name.replace("_", " ")
                        hdul[1].header["OBJECT"] = target_name.replace("_", " ")
                        hdul[0].header["SUMTYP"] = "FINAL"
                        hdul[1].header["SUMTYP"] = "FINAL"
                        hdul.flush()

            print(f"{target_name} sum completed!")
            print("----------------------------------------------------------------\n")

        return

        

        

        
            
        

