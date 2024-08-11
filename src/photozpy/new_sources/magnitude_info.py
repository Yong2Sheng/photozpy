import astropy.units as u
from astropy.table import QTable
import logging
import copy
logger = logging.getLogger(__name__)

class MagnitudeInfo():

    def __init__(self, filters, mag_type, inst_mags = None, inst_mag_errors = None, zero_points = None, cab_mags = None, cab_mag_errors = None):

        """

        Parameters
        ----------
        filters : list
        mag_type : str
            standard star or source
        inst_mags : astropy.units.quantity.Quantity
        zero_points : astropy.table.table.QTable
        cab_mags : astropy.units.quantity.Quantity
        """

        self._filters = filters
        self._mag_type = None

        ################################################################################################################ 
        
        if inst_mags is None:
            self._inst_mags = QTable(data = [-99]*len(self._filters)*u.mag, # the default is -99
                                     names = self._filters)
            
        elif isinstance(inst_mags, list):
            
            if len(self._filters) != len(inst_mags):
                logger.error("The number of instrumental magnitudes is not equal the number of filters!")
            else:
                self._inst_mags = QTable(data = inst_mags*u.mag, 
                                         names = self._filters)

        elif isinstance(inst_mags, QTable):
            self._inst_mags = inst_mags

        else:
            logger.error("Only QTable and list are supported for instrumental magnitudes!")
            
        ################################################################################################################    
        
        if inst_mag_errors is None:
            self._inst_mag_errors = QTable(data = [-99]*len(self._filters)*u.mag, # the default is -99
                                         names = self._filters)
            
        elif isinstance(inst_mag_errors, list):
            
            if len(self._filters) != len(inst_mag_errors):
                logger.error("The number of instrumental magnitude errors is not equal the number of filters!")
            else:
                self._inst_mag_errors = QTable(data = inst_mag_errors*u.mag,
                                               names = self._filters)
                
        elif isinstance(inst_mag_errors, QTable):
            self._inst_mag_errors = inst_mag_errors
        
        else:
            logger.error("Only QTable and list are supported for instrumental magnitude errors!")

        ################################################################################################################ 

        if zero_points is None:
            self._zero_points = QTable(data = [-99, -99, -99, -99]*u.mag, 
                                       names = self._filters)
        elif isinstance(zero_points, list):
            
            if len(self._filters) != len(zero_points):
                logger.error("The number of zero points is not equal the number of filters!")
            else:
                self._zero_points = QTable(data = zero_points*u.mag,
                                               names = self._filters)
                
        elif isinstance(zero_points, QTable):
            self._zero_points = zero_points
        
        else:
            logger.error("Only QTable and list are supported for zero points!")

        ################################################################################################################ 

        if cab_mags is None:
            self._cab_mags = QTable(data = [-99]*len(self._filters)*u.mag, # the default is -99
                                     names = self._filters)
        elif isinstance(cab_mags, list):

            if len(self._filters) != len(cab_mags):
                logger.error("The number of calibrated magnitudes is not equal the number of filters!")
            else:
                self._cab_mags = QTable(data = cab_mags*u.mag, 
                                         names = self._filters)
        elif isinstance(cab_mags, QTable):
            self._cab_mags = cab_mags
        else:
            logger.error("Only QTable and list are supported for calibrated magnitudes!")

        ################################################################################################################ 
    
        if cab_mag_errors is None:
            self._cab_mag_errors = QTable(data = [-99]*len(self._filters)*u.mag, # the default is -99
                                     names = self._filters)
            
        elif isinstance(cab_mag_errors, list):
            
            if len(self._filters) != len(cab_mag_errors):
                logger.error("The number of calibrated magnitudes is not equal the number of filters!")
            else:
                self._cab_mag_errors = QTable(data = cab_mag_errors*u.mag, 
                                         names = self._filters)
                
        elif isinstance(cab_mag_errors, QTable):
            self._cab_mag_errors = cab_mag_errors
        
        else:
            logger.error("Only QTable and list are supported for calibrated magnitude errors!")

        ################################################################################################################ 
    

    @property
    def inst_mags(self):
        return self._inst_mags

    @property
    def inst_mag_errors(self):
        return self._inst_mag_errors

    @property
    def cab_mags(self):
        return self._cab_mags

    @property
    def cab_mag_errors(self):
        return self._cab_mag_errors

    @property
    def zero_points(self):
        return self._zero_points

    @staticmethod
    def mag_operator(mag1, mag2, operation):

        """
        mag1 : float
        mag2 : float
        """

        if mag1 == -99 or mag2 == -99:
            return -99
        else:
            if operation == "+":
                return round(mag1 + mag2, 4)
            elif operation == "-":
                return round(mag1 - mag2, 4)
            else:
                logger.error(f"{operation} not supported!")
        

    def calculate_zero_points(self):


        for filter_name in self._inst_mags.colnames:
            
            zero = MagnitudeInfo.mag_operator(mag1 = self._cab_mags[filter_name].value[0], # .value return an array: [-18]
                                                    mag2 = self._inst_mags[filter_name].value[0], 
                                                    operation = "-")
            
            self._zero_points[filter_name] = zero*u.mag
            
            if self._mag_type == "source":
                logger.warning("You are using source magnitudes to find zero points!")
                
        return self._zero_points
    
    
    def calibrate_inst_mags(self):
        
        for filter_name in self._inst_mags.colnames:
            
            calibrated = MagnitudeInfo.mag_operator(mag1 = self._zero_points[filter_name].value[0], # .value return an array: [-18]
                                                    mag2 = self._inst_mags[filter_name].value[0], 
                                                    operation = "+")
            
            self._cab_mags[filter_name] = calibrated*u.mag
            self._cab_mag_errors = copy.deepcopy(self._inst_mag_errors)

        return
            