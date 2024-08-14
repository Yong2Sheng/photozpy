"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):
"""

import pandas as pd
from pathlib import Path
import numpy as np
from astropy import units as u
from astropy.units.quantity import Quantity
from astropy.table import QTable, Table
import os

class Telescope():
    
    """
    Manages the telecope, ccd and filter information.
    """
    
    def __init__(self, telescope, mode, ccd, ccd_gain = None, ccd_rdnoise = None, filters = None, lib_path = None):
        
        """
        Defines the telescope parameters. It reads the pre-defined telescope library and check if the telescope, ccd and mode combination is in the library.
        It also support the pass the parameters by the keyword argument: ccd_gain, ccd_rdnoise, filters.
        Note that the passed keyword will overwrite the values read from the library!
        
        Parameters
        ----------
        telescope: str; the name of the telescope.
        ccd: str; the name of ccd in the pre-defined telescope library.
        mode: str; mode determines the combination of the filters used.
        ccd_gain: astropy.units.quantity.Quantity; the gain of the ccd. Unit: electron/ADU
        ccd_rdnoise: astropy.units.quantity.Quantity; the read noise of the ccd. Unit: electron/pixel
        lib_path: pathlib.Path or str; the path to the alternative telescope library file.
        filters: list; list of filters.
        
        Returns
        -------
        
        """
        
        if lib_path is not None:
            if isinstance(lib_path,(Path, str)):
                lib_path = Path(lib_path)
            else:
                raise TypeError("Only str or pathlib.Path is supported for lib_path!")
            print("Reading telescope and ccd information from customed library.")
            self.lib_df = pd.read_csv(lib_path, sep = ",", header=0)
        else:
            directory = os.path.dirname(__file__)
            self.lib_df = pd.read_csv(directory+"/telescope_library.csv", sep = ",", header=0)
        
        self._telescope = telescope
        self._mode = mode
        self._ccd = ccd
        self._ccd_gain = ccd_gain
        self._ccd_rdnoise = ccd_rdnoise
        self._filters = filters
                
        return
    
    def get_telescope_parameters(self, para_name):
                 
        """
        Read and return the telescope parameters you want.
        
        Parameters
        ----------
        para_name: str; the parameter name.
        
        Returns
        -------
        para_value: list or float; the value of the parameter
        
        """
        
        if para_name == "ccd_gain":
            para = self._ccd_gain
            
        elif para_name == "ccd_rdnoise":
            para = self._ccd_rdnoise
            
        elif para_name == "filters":
            para = self._filters
                
        else:
            raise ValueError("The parameter name you input doesn't exist! Try ccd_gain, ccd_rdnoise or filters!")
        
        if para is None:
            # if the parammeter is None, it will read from the library.
            filtered = self.lib_df[(self.lib_df["telescope"] == self._telescope) & (self.lib_df["ccd"] == self._ccd) & (self.lib_df["mode"] == self._mode)]
            
            if filtered.shape[0] == 0:
                raise ValueError("The telescope, ccd and mode combination doesn't exist in the library!")
            elif filtered.shape[0] >= 2:
                raise ValueError("More than one set of telescope, ccd and mode combinations in the library!")
            elif filtered.shape[0] == 1:
                if para_name == "ccd_gain" or para_name == "ccd_rdnoise":
                    para_value = float(filtered.loc[:,para_name].to_numpy()[0])
                elif para_name == "filters":
                    para_value = filtered.loc[:,para_name].to_string(header = False, index = False).split("/")
            
        else:
            para_value = para
        
        return para_value
    
    @property
    def telescope(self):
        
        return self._telescope
    
    @property
    def ccd(self):
        
        return self._ccd
    
    @property
    def mode(self):
        
        return self._mode
    
    @property
    def ccd_gain(self):
        
        self._ccd_gain = self.get_telescope_parameters("ccd_gain")
        
        if not isinstance(self._ccd_gain, Quantity):
            self._ccd_gain = self._ccd_gain*u.electron/u.adu
        
        return self._ccd_gain
    
    @property
    def ccd_rdnoise(self):
        
        self._ccd_rdnoise = self.get_telescope_parameters("ccd_rdnoise")
        
        if not isinstance(self._ccd_rdnoise, Quantity):
            self._ccd_rdnoise = self._ccd_rdnoise*u.electron/u.pix
        
        return self._ccd_rdnoise
    
    @property
    def filters(self):
        
        return self.get_telescope_parameters("filters")

    @property
    def telescope_summary(self, save = False, save_path = None):
        
        """
        Produce a summary of the telescope by dictionary.
        """
        
        self.summary = {"telescope": self._telescope,
                        "mode": self._mode,
                        "ccd": self._ccd,
                        "ccd_gain": self.ccd_gain,
                        "ccd_rdnoise": self.ccd_rdnoise,
                        "filters": self.filters}
        
        return self.summary
    
    def save_telescope_summary(self, save_path = None):
        
        """
        Save the telescope_summary as a ecsv file.
        
        Parameters
        ----------
        save_path: str or pathlib.Path; the path to save the ecsv file, including the file name.
        
        Returns
        -------
        None
        """
        
        a = [self.telescope]
        b = [self.mode]
        c = [self.ccd]
        d = [self.ccd_gain]
        e = [self.ccd_rdnoise]
        f = [self.filters]

        qtable = QTable([a, b, c, d, e, f],
                        names=('telescope', 'mode', 'ccd', 'ccd_gain', 'ccd_rdnoise', 'filters'),
                        meta={'name': 'Telescope information'},
                        dtype = [str, str, str, float, float, list])

        qtable.write(save_path,format="ascii.ecsv")

        return

    def map_filters(self, name, **filter_dict):
        """
        Map filter names to standard filter names. 

        Paremeters
        ----------
        name: the input filter name to be matched to a standard one.
        filter_dict: customized filter maping dictionary.

        Return
        ------
        mapped_filter: str; the name of the mapped standard filter.
        """

        if filter_dict == {}:
            # use the default mappers
            mappers = {"Sloan g'2": "SDSS_g'", 
                       "Sloan r'2": "SDSS_r'", 
                       "Sloan i'2": "SDSS_i'", 
                       "Sloan z'2": "SDSS_z'", 
                       "SDSS_g": "SDSS_g'", 
                       "SDSS_r": "SDSS_r'", 
                       "SDSS_i": "SDSS_i'", 
                       "SDSS_z": "SDSS_z'",
                       "Bessell I": "Bessell_I"}#,
                       # "B": "ubb",
                       # "UVM2": "um2",
                       # "U": "uuu",
                       # "V": "uvv",
                       # "UVW1": "uw1",
                       # "UVW2": "uw2"}
        else:
            mappers = filter_dict
            
        mapped_filters = [value for key,value in mappers.items() if key == name]
        mapped_filters = [*set(mapped_filters)] # remove the duplicated filter names
        
        if len(mapped_filters) == 0:
            raise ValueError("Mapping the input filter failed! No mapped filter {name} was found!")
        elif len(mapped_filters) > 1:
            print(mapped_filters)
            raise ValueError("Mapping the input filter failed! More than two mapped filters are found!")
        else:
            return mapped_filters[0]
        
        
    def __eq__(self, other):
        
        if self.telescope_summary == other.telescope_summary:
            return True
        else:
            return False
    
                
        
        