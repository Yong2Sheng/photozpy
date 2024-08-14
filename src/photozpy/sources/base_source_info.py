from .magnitude_info import MagnitudeInfo

class BaseSourceInfo():

    def __init__(self, source_name, telescope, source_type, skycoord = None, file_pattern = None, magnitudes = None):

        """
        Define the source information dictionary

        Parameters
        source_name : str
            The source name
        telescope : photozpy.telescope.Telescope
            The telescope used to make the observation
        source_type : str
            The type of the source: standard star or target
        skycoord : astropy.coordinates.SkyCoord
            The coordinates of the sources
        file_pattern : str
            The unix-style file name pattern
        magnitudes : photozpy.new_sources.magnitude_info.MagnitudeInfo
        ----------
        """

        self._source_name = source_name
        
        self._source_type = source_type

        self._telescope = telescope
        
        self._skycoord = skycoord
        
        if self._skycoord is not None:
            self._skycoord.info.name = self._source_name
        
        self._file_pattern = file_pattern

        self._magnitudes = MagnitudeInfo(filters = self._telescope.filters,
                                         mag_type = self._source_type, 
                                         inst_mags = None, 
                                         inst_mag_errors = None, 
                                         zero_points = None, 
                                         cab_mags = None, 
                                         cab_mag_errors = None)

        self._source_dict = {"source_name": self._source_name, 
                             "source_type": self._source_type, 
                             "telescope": self._telescope, 
                             "skycoord": self._skycoord,
                             "file_pattern": self._file_pattern, 
                             "magnitudes": self._magnitudes}

    @property
    def source_dict(self):
        return self._source_dict
    
    @property
    def source_name(self):
        return self._source_name

    @property
    def telescope(self):
        return self._telescope

    @property
    def source_type(self):
        return self._source_type

    @property
    def skycoord(self):
        return self._skycoord
    
    @skycoord.setter
    def skycoord(self, new_skycoord):
        self._skycoord = new_skycoord
        self._skycoord.info.name = self._source_name

    @property
    def file_pattern(self):
        return self._file_pattern

    @property
    def magnitudes(self):
        return self._magnitudes