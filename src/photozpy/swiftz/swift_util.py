"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- Download Swift data based on the quires.
- Unzip the data
"""

from pathlib import Path
from swifttools.swift_too import Data, ObsQuery, TOORequests
import pandas as pd
import numpy as np
import copy
from ..convenience_functions import create_folder, del_then_create_folder, ungz_file
import shutil 
from astropy.io import fits

class SwiftDownload():

    def __init__(self, download_dir):

        """
        Initialize the class.
        
        Parameters
        ----------
        root_dir : pathlib.Path or str
            The directory where the data will be stored.
        """

        self.download_dir = Path(download_dir)


    def _standardize_obsquery_input(self, source_dict):
        
        """
        Standardize the obsquery input for this class.
        This function prioritizes the obsquery input in this order: obsid, targetid, and skycoord while the windows (`yyyy-mm-dd`) are always there.

        Parameters
        ----------
        source_dict : dict
            The dictionary contains the source information

        Returns
        -------
        dict
            The dictionary that can be used to query and download data.
        """

        if "obsid" in source_dict.keys():
            return {"name":source_dict["name"],
                    "obsid": source_dict["obsid"],
                    "window_lower": source_dict["window_lower"],
                    "window_upper": source_dict["window_upper"]}
            
        elif "targetid" in source_dict.keys():
            return {"name":source_dict["name"],
                    "targetid": source_dict["targetid"],
                    "window_lower": source_dict["window_lower"],
                    "window_upper": source_dict["window_upper"]}

        elif "ra" and "dec" in source_dict.keys():
            ra_dec = {key: value for key, value in source_dict.items() if key in {"ra", "dec"}}
            skycoord = SkyCoord(**ra_dec, unit = (u.deg, u.deg), frame = "icrs")
            
            return {"name":source_dict["name"],
                    "skycoord": skycoord,
                    "window_lower": source_dict["window_lower"],
                    "window_upper": source_dict["window_upper"]}

    def set_obsquery_info(self, source_dict):

        """
        Set up the obsquery dictionary used to download the data from the input dictionary `source_dict`.
        The source dictionary must include the information of source name, obsid or target id or ra&dec, time window.

        Paremeters
        ----------
        source_dict : dict
            The source dictionary. It must have the keys of `name`, `obsid` or `targetid` or `ra`&`dec`, `window_upper`, and `window_lower`.
            The format for time window is `yyyy-mm-dd`

        Returns
        -------
        dict
            The query dictionary that will be used to query and download the data.
        """

        self.obsquery_info = self._standardize_obsquery_input(source_dict)

        self.master_info = pd.DataFrame.from_dict(self.obsquery_info)

        return self.obsquery_info

    def set_obsquery_info_from_csv(self, csv_file):

        """
        Set up the obsquery dictionary used to download the data from the csv file.
        The csv file must include the information of source name, obsid or target id or ra&dec, time window.

        Parameters
        ----------
        csv_file : pathlib.Path or str
            The csv file (with `,` as separator) must contain the following headers:
            `name`, `obsid` or `targetid` or `ra`&`dec`, `window_upper`, and `window_lower`

        Returns
        -------
        dict
            The query dictionary that will be used to query and download the data.
        """

        self.master_info = pd.read_csv(csv_file, sep = ",")  # this is the DataFrame that contains all the information. New information will be added during the analysis pipeline.
        source_dict = self.master_info.to_dict(orient='list')

        self.obsquery_info = self._standardize_obsquery_input(source_dict)  # this is the dictionary usd to query and download data

        return self.obsquery_info

    @property
    def get_obsquery_info(self):
        
        """
        The dictionary used to query and download data.
        """

        return self.obsquery_info

    @property
    def get_master_info(self):

        """
        The complete master information of the source.
        """

        return self.master_info

    @property
    def get_source_names(self):

        """
        Get the names of the downloaded sources. This is the source name from the meta information used to download the data.
        It's not necessary the source name in the fits header.
        """

        return self.master_info["name"].to_list()
        
    @staticmethod
    def slice_dict(input_dict, slice_index = None, slice_value = None, slice_key = None, remove_keys = None):

        """
        Assume that the values of each keys corresponding to a specific source, and you want to get the information of a specific source. 
        It means that you slice the dictionary vertically respect to the keys.

        For example you want to slice this dictionary:
        
        {'name': ['4FGL J0700.5-6610', 'B3 0850+443', '87GB 213913.0+293303'],
         'targetid': [38456, 14051, 16150],
         'window_lower': ['2024-01-01', '2021-01-01', '2023-07-10'],
         'window_upper': ['2024-03-10', '2021-03-10', '2023-08-10']}

        1. You can use `slice_index = [1, 2]` to get the name, targetid, window_lower and window_upper:
            {'name': ['B3 0850+443', '87GB 213913.0+293303'],
             'targetid': [14051, 16150],
             'window_lower': ['2021-01-01', '2023-07-10'],
             'window_upper': ['2021-03-10', '2023-08-10']}

        2. You can also use `slice_value = ['B3 0850+443', '87GB 213913.0+293303']` and `slice_key = "name"` to get the same output dictionary
            {'name': ['B3 0850+443', '87GB 213913.0+293303'],
             'targetid': [14051, 16150],
             'window_lower': ['2021-01-01', '2023-07-10'],
             'window_upper': ['2021-03-10', '2023-08-10']}
        
        Parameters
        ----------
        input_dict : dict
            The input dictionary. 
        slice_index : iterable object, optional
            The index/indices of the element to be sliced from each key.
        slice_value : iterable object, optional
            The value(s) to be sliced for a specific `slice_key`.
        slice_key : optional
            The key used to find the value(s) to be sliced.
        remove_keys : list or str, optional
            The key(s) to be removed from the new sliced dictionary (the default is `None`, which means no keys will be removed).

        Returns
        -------
        new_dict : dict
            The new sliced dictionary.
        """
        new_dict = {}

        # if the slice_value(s) and slice_key are given, use them to find the index(indices) the slice_value(s)
        if slice_value is not None:
            if not isinstance(slice_value, list):
                slice_value = [slice_value]
            slice_index = [input_dict[slice_key].index(i) for i in slice_value]
        else:
            if isinstance(slice_index, int):
                slice_index = [slice_index]

        # construct the new sliced dictionary by the index/indices
        for key, value in input_dict.items():
            if not isinstance(value, list):
                raise ValueError("The value of each key must be a list to be sliced by index.")   
            else:
                new_dict[key] = [value[i] for i in slice_index]

        if remove_keys is not None:
            if isinstance(remove_keys, str):
                remove_keys = [remove_keys]
            for remove_key in remove_keys:
                del new_dict[remove_key]

        return new_dict

    def download_swift_data(self, radius = 5/60, uvotmode = "0x30ed", organize = False):

        """
        Download the Swift data. Note it will overwrite the source data in the same location you downloaded before. 
        I should make it more flexible ...
        
        Paremeters
        ----------
        radius : float, optional
            The searching radius.
        uvotmode : str, optional
            The uvot mode of the data you need.
        organize : bool, optional
            Organize the downloaded data for further analyais with this pipeline. (the default is `False`, which means do not organize)
        """

        # check how many sources are in the obsquery_info dict
        target_num = len(self.obsquery_info["name"])

        for i in np.arange(target_num):
            target_info = SwiftDownload.slice_dict(input_dict = self.obsquery_info, slice_index = int(i), remove_keys = ["window_lower", "window_upper"])
            target_name = target_info["name"][0]  # use [0] here because the values for each key is always a list event if there is only one element.

            oq = ObsQuery(radius = radius, 
                          begin = self.obsquery_info["window_lower"][i], end = self.obsquery_info["window_upper"][i], 
                          **target_info)
            id_ = "00000000000"  # this variable will be used to avoid downloading the same data file mutiple times by comparing the observation id
            print(f"Downloading data for {target_name} ...")

            target_dir = del_then_create_folder(folder_dir = self.download_dir / target_name.replace(" ", "_"))  # linux doesn't work with spaces in path
            

            for i in np.flip(np.arange(-len(oq), 0)):
                if oq[i].obsid == id_:
                    print(f"{oq[i].obsid} has been downloaded/examined, skipping ...")
                else:
                    if oq[i].uvot_mode == uvotmode:
                        date_ = oq[i].begin.strftime("%Y-%m-%d %H:%M:%S")
                        print(f"{oq[i].obsid} on {date_} is being downloaded, the uvot mode is {oq[i].uvot_mode}")
                        oq[i].download(uvot = True, outdir = target_dir, match = ["*/uvot/image/*sk.img.gz", "*/auxil/*"])
                        id_ = oq[i].obsid
                    else:
                        print(f"{oq[i].obsid} uvod mode {oq[i].uvot_mode} is not the one you requested as {uvotmode}, skipping ...")
                        id_ = oq[i].obsid

        
        if organize is True:
            source_names = self.obsquery_info["name"]
            for source_name in source_names:
                src_dir = self.download_dir / source_name.replace(" ", "_")  # linux doesn't work with spaces in path
                files = list(src_dir.rglob('*sk*'))
                for file in files:
                    fits_path = ungz_file(file_path = file, out_dir = src_dir, return_path = True)
                    
                    # with fits.open(fits_path, mode = "update") as hdul:
                    #     hdul[0].header["OBJECT"] = source_name
                    #     hdul.flush()

        return
                