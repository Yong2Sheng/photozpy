import astropy.units as u
from astropy.table import QTable
from astropy.units.quantity import Quantity
import copy
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import logging
logger = logging.getLogger(__name__)

class MagnitudeInfo():

    def __init__(self, filters, mag_type, inst_mags = None, inst_mag_errors = None, zero_points = None, cab_mags = None, cab_mag_errors = None, calibrated = False):

        """

        Parameters
        ----------
        filters : list
        mag_type : str
            standard star or target
        inst_mags : None or list: [[g_mags], [r_mags], [i_mags],[z_mags]]
        zero_points : astropy.table.table.QTable
        cab_mags : astropy.units.quantity.Quantity
        """

        self._filters = filters
        self._mag_type = mag_type
        self._inst_mags = MagnitudeInfo._standarize_mags(inst_mags, self._filters)
        self._inst_mag_errors = MagnitudeInfo._standarize_mags(inst_mag_errors, self._filters)
        self._zero_points = MagnitudeInfo._standarize_mags(zero_points, self._filters)
        self._cab_mags = MagnitudeInfo._standarize_mags(cab_mags, self._filters)
        self._cab_mag_errors = MagnitudeInfo._standarize_mags(cab_mag_errors, self._filters)
        self._calibrated= calibrated
        self._diff_dict = {}
        

    @staticmethod
    def _standarize_mags(mags, filters):
        
        """
        Standarize the magnitudes input.
        
        Parameters
        ----------
        mags : list or QTable
        filters : list
        """
        
        if mags is None:
            return QTable([[-99]*u.mag, [-99]*u.mag, [-99]*u.mag, [-99]*u.mag], # the default is -99
                          names = filters)
  
        elif isinstance(mags, list):
            
            if isinstance(mags[0], (Quantity, list)):  # [[-99, 23]*u.mag, [-99, 56]*u.mag, [-99, 46]*u.mag, [-99, 233]*u.mag], [[-99, 23], [-99, 56], [-99, 46]]
                
                rows = [len(i) for i in mags]
                
                if rows.count(rows[0]) != len(rows):  # check if the number of magnitudes in different filters are equal
                    raise ValueError("The number of magnitudes in filters are not equal!")
                    
                else:
                    if isinstance(mags[0], Quantity):  # [[-99, 23]*u.mag, [-99, 56]*u.mag, [-99, 46]*u.mag, [-99, 233]*u.mag]
                        return QTable(mags, 
                                      names = filters)
                    
                    elif isinstance(mags[0], list):  # [[-99, 23], [-99, 56], [-99, 46]]
                        return QTable([i*u.mag for i in mags], 
                                      names = filters)
            
            elif isinstance(mags[0], float):  #[-99, 233]] or [23, 56, 46, 233]
                return QTable(mags*u.mag, 
                              names = filters)
            
        elif isinstance(mags, QTable):
            return mags
        
        else:
            raise TypeError("Only QTable and list are supported for magnitude inputs!")
        

    @property
    def inst_mags(self):
        return self._inst_mags
    
    @inst_mags.setter
    def inst_mags(self, new_inst_mags):
        self._inst_mags = MagnitudeInfo._standarize_mags(new_inst_mags, self._filters)
        
    @property
    def inst_mag_errors(self):
        return self._inst_mag_errors
    
    @inst_mag_errors.setter
    def inst_mag_errors(self, new_inst_mag_errors):
        self._inst_mag_errors = MagnitudeInfo._standarize_mags(new_inst_mag_errors, self._filters)

    @property
    def cab_mags(self):
        return self._cab_mags
    
    @cab_mags.setter
    def cab_mags(self, new_cab_mags):
        self._cab_mags = MagnitudeInfo._standarize_mags(new_cab_mags, self._filters)

    @property
    def cab_mag_errors(self):
        return self._cab_mag_errors
    
    @cab_mag_errors.setter
    def cab_mag_errors(self, new_cab_mag_errors):
        self._cab_mag_errors = MagnitudeInfo._standarize_mags(new_cab_mag_errors, self._filters)
        
    @property 
    def diff_zero_points(self):
        return self._diff_dict

    @property
    def zero_points(self):
        return self._zero_points
    
    @zero_points.setter
    def zero_points(self, new_zero_points):
        self._zero_points = MagnitudeInfo._standarize_mags(new_zero_points, self._filters)

    @staticmethod
    def single_mag_operator(mag1, mag2, operation):

        """
        mag1 : mag type
        mag2 : mag type
        """

        if mag1 == -99*u.mag or mag2 == -99*u.mag:
            return -99*u.mag
        else:
            if operation == "+":
                return mag1 + mag2
            elif operation == "-":
                return mag1 - mag2
            else:
                raise TypeError(f"{operation} not supported!")
                
    @staticmethod
    def mag_operator(mag_list1, mag_list2, operation):
        
        list_out = []
        
        if len(mag_list1) == len(mag_list2):
            pass
        
        elif len(mag_list1) == 1:
            mag_list1 = np.repeat(mag_list1.value[0], len(mag_list2))*u.mag
            
        elif len(mag_list2) == 1:
            mag_list2 = np.repeat(mag_list2.value[0], len(mag_list1))*u.mag
            
        else:
            raise ValueError("Can't broadcase magnitude quantity with more than one value!")
            
        for mag1, mag2 in zip(mag_list1, mag_list2):
            
            _mag = MagnitudeInfo.single_mag_operator(mag1, mag2, operation)
            
            list_out += [_mag.value] # first store calculated mags by float
        
        return list_out*u.mag  # return mag quantity
                
                
    @staticmethod
    def remove_outlier(array, sigma = 3, verbose = False):
    
        zscores = stats.zscore(array)
    
        new_array = []
    
        for element, zscore in zip(array, zscores):
            if np.abs(zscore) <= sigma:
                new_array += [element]
            else:
                if verbose is True:
                    logger.info(f"{element} is an outlier with zscore of {zscore}!")
                else:
                    pass
    
        return np.array(new_array)
    
    def calculate_zero_points(self, sigma = 3, update = True):
        
        if len(self._cab_mags) <= 8:          
            logger.warning(f"You only have {len(self._cab_mags)} standard stars, which is less than 8.")
            logger.warning("Please be careful that the calculated zero points might be biased.")
        
        if self._mag_type == "target":
            
            raise TypeError("You should NOT use target magnitudes to find zero points! Please use standard stars instead.")
            
        
        calculated_zero_points = copy.deepcopy(self._zero_points)
        
        for filter_name in calculated_zero_points.colnames:
            
            
            zero = MagnitudeInfo.mag_operator(mag_list1 = self._cab_mags[filter_name], # .value return a list even with one number: [-18]
                                              mag_list2 = self._inst_mags[filter_name], 
                                              operation = "-")
            
            self._diff_dict[filter_name] = zero.value
            
            zero = MagnitudeInfo.remove_outlier(zero.value, sigma = sigma)  # remove outlier and give it an unit
            
            zero = np.mean(zero)*u.mag
            
            calculated_zero_points[filter_name] = zero  # fill up the QTable 
            
            
            
            
        if update is True:
            self._zero_points = calculated_zero_points
        
        return calculated_zero_points

    
    def calibrate_inst_mags(self, update = True):
        
        calibrated_inst_mags = {}
        calibrated_inst_mag_erros = {}
        
        
        for filter_name in self._inst_mags.colnames:
            
            if self._zero_points[filter_name].value[0] == -99:
                
                logger.warning(f"The zero point in {filter_name} is -99!")
            
            calibrated = MagnitudeInfo.mag_operator(mag_list1 = self._zero_points[filter_name],  # .value return a list even with one number: [-18]
                                                    mag_list2 = self._inst_mags[filter_name], 
                                                    operation = "+")
            
            calibrated_inst_mags[filter_name] = calibrated
            calibrated_inst_mag_erros = copy.deepcopy(self._inst_mag_errors)

        if update is True:
            self._cab_mags = QTable(calibrated_inst_mags)
            self._cab_mag_errors = QTable(calibrated_inst_mag_erros)
            
        if self._calibrated is True:
            
            logger.warning("The targets are already calibrated!")
            
        self._calibrated = True

        return
    
    def plot_zero_point_statistics(self, save = True):
        
        nfilters = len(self._filters)
        nrows = int(np.ceil(nfilters/2))
        
        fig, axes = plt.subplots(nrows, 2, figsize=(12,10), sharex=False) ##sharex=True
        
        for idx, filter_name in enumerate(self._filters):
            
            axes[int(np.floor(idx/2)), idx % 2].scatter(x = [filter_name]*len(self._diff_dict[filter_name]), 
                                                   y = self._diff_dict[filter_name], 
                                                   label = "offsets", color = "skyblue")
            
            axes[int(np.floor(idx/2)), idx % 2].scatter(x = [filter_name], 
                                                  y = self._zero_points[filter_name], 
                                                  label = r"$\sigma$-clipped average", color = "lime")
            
            axes[int(np.floor(idx/2)), idx % 2].hlines(y=self._zero_points[filter_name].value[0]-0.1, color='r', xmin = 0, xmax = 1, label = "sigma-clipped - 0.1")
            axes[int(np.floor(idx/2)), idx % 2].hlines(y=self._zero_points[filter_name].value[0]+0.1, color='r', xmin = 0, xmax = 1, label = "sigma-clipped + 0.1")
            
            axes[int(np.floor(idx/2)), idx % 2].set_title(filter_name)
            
            
            axes[int(np.floor(idx/2)), idx % 2].legend()
            
            
        plt.show()
            
        if save is True:
            fig.savefig("Zero_points_statistics", dpi=300)
            
        #plt.close()
        
        return
        
        
            