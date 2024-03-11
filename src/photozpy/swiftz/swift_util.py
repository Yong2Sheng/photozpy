"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- Download Swift data based on the quires.
- Generate metadata of the data
- p
"""

from pathlib import Path
from swifttools.swift_too import Data, ObsQuery, TOORequests
import pandas as pd
import numpy as np
import copy
from ..convenience_functions import create_folder, del_then_create_folder

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


    def _standardize_obsquery_input(self, **kwargs):
        
        """
        Standardize the obsquery input for this class.
        This function prioritizes the obsquery input in this order: obsid, targetid, and skycoord while the windows (`yyyy-mm-dd`) are always there.
        """

        if "obsid" in kwargs.keys():
            return {"name":kwargs["name"],
                    "obsid": kwargs["obsid"],
                    "window_lower": kwargs["window_lower"],
                    "window_upper": kwargs["window_upper"]}
            
        elif "targetid" in kwargs.keys():
            return {"name":kwargs["name"],
                    "targetid": kwargs["targetid"],
                    "window_lower": kwargs["window_lower"],
                    "window_upper": kwargs["window_upper"]}

        elif "ra" and "dec" in kwargs.keys():
            ra_dec = {key: value for key, value in kwargs.items() if key in {"ra", "dec"}}
            skycoord = SkyCoord(**ra_dec, unit = (u.hourangle, u.deg), frame = "icrs")
            
            return {"name":kwargs["name"],
                    "skycoord": skycoord,
                    "window_lower": kwargs["window_lower"],
                    "window_upper": kwargs["window_upper"]}

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
            The query dictionary that will be used to download the data.
        """

        self.obsquery_info = self._standardize_obsquery_input(**source_dict)

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
            The query dictionary that will be used to download the data.
        """

        df = pd.read_csv(csv_file, sep = ",")
        source_dict = df.to_dict(orient='list')

        self.obsquery_info = self._standardize_obsquery_input(**source_dict)

        return self.obsquery_info

    @property
    def get_obsquery_info(self):
        
        """
        The obsquery information.
        """

        return self.obsquery_info
        
    @staticmethod
    def slice_dict_by_index(input_dict, index, remove_keys = None):

        """
        input_dict : dict
            The input dictionary. 
        index : int
            The index of the element to be sliced in the value of each key.
        remove_keys : list or str, optional
            The key(s) to be removed from the new sliced dictionary (the default is `None`, which means no keys will be removed).
        """
        new_dict = {}

        for key, value in input_dict.items():

            if not isinstance(value, list):
                raise ValueError("The value of each key must be a list, so it can be sliced by index.")
            else:
                new_dict[key] = value[index]

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
            target_info = SwiftDownload.slice_dict_by_index(input_dict = self.obsquery_info, index = i, remove_keys = ["window_lower", "window_upper"])
            tagrt_name = target_info["name"]

            oq = ObsQuery(radius = radius, 
                          begin = self.obsquery_info["window_lower"][i], end = self.obsquery_info["window_upper"][i], 
                          **target_info)
            id_ = "00000000000"  # this variable will be used to avoid downloading the same data file mutiple times by comparing the observation id
            print(f"Downloading data for {tagrt_name} ...")

            target_dir = del_then_create_folder(folder_dir = self.download_dir / tagrt_name)
            

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
    
            