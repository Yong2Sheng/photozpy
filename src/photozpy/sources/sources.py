# from .base_source_info import BaseSourceInfo
# from .magnitude_info import MagnitudeInfo
from itertools import product
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
from .base_source_info import BaseSourceInfo
import logging
logger = logging.getLogger(__name__)

class Sources():
    
    def __init__(self, *sources):

        """
        Parameters
        ----------
        sources : photozpy.new_sources.base_source_info.BaseSourceInfo

        """

        self._sources = sources
        self._source_names = [source.source_name for source in self._sources]
        if not Sources.check_uniqueness(self._source_names):  # check if there are sources with same name
            logger.error("You have at least two sources that have the same source name!")
            
    @classmethod
    def from_list(cls, source_list):
        
        return cls(*source_list)
    
    @classmethod 
    def from_csv(cls, csv_path, telescope, units = (u.hourangle, u.deg)):

        df = pd.read_csv(csv_path, sep = ",")
        
        source_list = []
        
        for index, row in df.iterrows():
        
            if not row["ra"] is np.nan:
                _skycoord = SkyCoord(ra = [row["ra"]], 
                                     dec = [row["dec"]], 
                                     unit = units)
            else:
                _skycoord = None
        
        
            _source_info = BaseSourceInfo(source_name = row["source_name"], 
                                          telescope = telescope, 
                                          source_type = row["source_type"], 
                                          skycoord = _skycoord, 
                                          file_pattern = row["file_pattern"], 
                                          magnitudes = None)
            source_list += [_source_info]
            
        return cls(*source_list)
        

    @staticmethod
    def check_uniqueness(in_list):
    
        length = len(in_list)
    
        idx = 0
    
        while idx <= length-2:
    
            if in_list[idx] in in_list[idx+1:]:
                return False
            else:
                idx += 1
        return True

    @property
    def nsources(self):
        return len(self._sources)
    
    @property 
    def coords(self):
        return [i.skycoord for i in self._sources]

    @property
    def source_names(self):
        return self._source_names

    @property
    def targets(self):
        return self.filter_sources(source_type = "target")
    
    @property 
    def target_names(self):
        return [i.source_name for i in self.targets]
    
    @property 
    def target_coords(self):
        return [i.skycoord for i in self.targets]
    
    @property
    def standard_stars(self):
        return self.filter_sources(source_type = "standard star")
    
    @property 
    def standard_star_names(self):
        return [i.source_name for i in self.standard_stars]
    
    @property 
    def standard_star_coords(self):
        return [i.skycoord for i in self.standard_stars]
    
    @property
    def sources(self):
        return self._sources
    
    @property 
    def telescope(self):
        
        source1 = self._sources[0]
        
        for source in self._sources[1:]:
            if source1.telescope != source.telescope:
                raise ValueError("The telescope is the sources are different!")
            
        return source1.telescope

    @staticmethod
    def unwarp_dictionary(dictionary):
    
        """
        Unwarp a dictionary to produce a dictionary list whose element is a dictionary that contains all the combinations of 
        values for the headers.
    
        Parameters
        ----------
        dictionary: dict; the input dictionary
    
        Returns
        -------
        dict_list: list; the list of all dictionaries that exhausts all the combination of the values.
        """
    
        # make sure all the values are lists instead of strings or numbers
        for header, value in dictionary.items():
            if not isinstance(value, list):
                dictionary[header] = [value]
        
        # get all the combinations of values from different headers
        values = list(product(*list(dictionary.values())))
    
        headers = list(dictionary.keys())
    
        dict_list = []
        for i in values:
            dict_ = dict(zip(headers, i))  # assemble headers and values. Note that headers are always the same. We just iterate through the combination of the values.
            dict_list += [dict_]
    
        return dict_list

    @staticmethod
    def _filter_dict(dictionary, **headers_values):
    
        """
        Check if the headers_values match the ones in dictionary. 
        Return False if there is more than one that do not match.
        Return True only if all the key and values match.
    
        dictionary = {'source_type': 'target', 'source_name': 'PKS'}
        headers_values = {"source_type": "target", 'source_name': 'PKS'}
        Will return True
        """
    
        for key, value in headers_values.items():
    
            if not key in dictionary.keys():  # check if the keys in headers_values are also in the dictionary
                logger.warning("The keys of the dictionary do not match!")
                return False
    
            else:
                if dictionary[key] != value:
                    return False
                
        return True
            
    def filter_sources(self, **headers_values):

        """
        Find the sources that match the keys and values in headers_values.
        """

        dict_filters = Sources.unwarp_dictionary(headers_values)

        filtered_sources = []
        
        for source in self._sources:
            
            for dict_filter in dict_filters:

                if Sources._filter_dict(source.source_dict, **dict_filter):
                    filtered_sources += [source]
                    
        return filtered_sources


    @property
    def zpoints(self):
        try:
            self._zpoints
        except AttributeError:
            return self.calculate_zpoints()
        else:
            return self._zpoints


    def calculate_zpoints(self, clip_sigma = 3):

        """
        It only works with the case that one standard star image in all the sources
        """

        standard_star_magnitudes = self.standard_stars[0].magnitudes
        self._zpoints = standard_star_magnitudes.calculate_zero_points(update = True, sigma = clip_sigma)

        return self._zpoints  # This is a QTable

    def calibrate_targets(self, zpoints = None):

        targets = self.targets

        for target in targets:

            mag = target.magnitudes
            mag.zero_points = self.zpoints
            mag.calibrate_inst_mags()

        return
    
    def list_calibrated_targets(self):
        
        telescope = self.telescope
        
        for target in self.targets:
            print(target.source_name)
            # print magnitudes
            print(f"The magnitudes in {telescope.filters} are:")
            for filter_names in telescope.filters:
                for i in np.around(target.magnitudes.cab_mags[filter_names].value, 4):
                    print(i)
            
            # print magnitude errors
            print(f"The magnitude errors in {telescope.filters} are:")
            for filter_names in telescope.filters:
                for i in np.around(target.magnitudes.cab_mag_errors[filter_names].value, 4):
                    print(i)
                    
            # print detection significance
            print(f"The detection significance in {telescope.filters} are:")
            for filter_names in telescope.filters:
                for i in np.around(target.magnitudes.detection_significance[filter_names].value, 4):
                    print(i)
                
            print("======================================================================")
            print("======================================================================")
    
    
    def __getitem__(self, value):
        
        if isinstance(value, str): 
            
            return self.filter_sources(source_name = value)
        
        elif isinstance(value, int):
            return self._sources[value]
        
        elif isinstance(value, slice):
            
            # set up the start of the slicing
            if value.start is None:
                start = 0
            elif value.start < 0:
                start = value.start + self.nsources
            else:
                start = value.start
            
            # set up the stop of the slicin
            if value.stop is None:
                stop = -1 + self.nsources
            elif value.stop < 0:
                stop = value.stop + self.nsources
            else:
                stop = value.stop
                
            indices = np.arange(start, stop, value.step)
            return [self._sources[idx] for idx in indices]
    