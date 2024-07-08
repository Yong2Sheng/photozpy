"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- handles the source information
"""

from astropy.coordinates import SkyCoord

class Sources():

    def __init__(self, object_name = None, skycoords = None):

        """
        Define the initial source dictionary.

        Parameters
        ----------
        object: string or list; the object name
        skycoords: astropy.coordinates.SkyCoord; the skycoords of the sources in the object

        Returns
        -------
        source_dict: dictionary; the source dictionary
        """

        # initialize the source dictionary
        ## The object name is the value of the OBJECT header of the fits file.
        ## The object name is used to identify the image files.
        ## One object can correspond to multiple SkyCoords, for example, standard stars.
        if object_name == None and skycoords == None:
            self.source_dict = {"Object":[], 
                                "SkyCoords":[]}
        elif object_name != None and skycoords != None:
            if not isinstance(object_name, list):
                object_name = [object_name]
            if not isinstance(skycoords, list):
                skycoords = [skycoords]
            self.source_dict = {"Object":object_name, 
                                "SkyCoords":skycoords}
        else:
            raise ValueError("Both object and skycoords should be None or not None!")

    def add_object(self, object_name, skycoords):

        """
        Add a new obejct and skycoords to the source dictionary

        Parameters
        ----------
        object: string or list; the object name
        skycoords: astropy.coordinates.SkyCoord; the skycoords of the sources in the object

        Returns
        -------
        source_dict
        """

        if not isinstance(object_name, list):
            object_name = [object_name]

        if not isinstance(skycoords, list):
            skycoords = [skycoords]

        self.source_dict["Object"] = self.source_dict["Object"] + object_name
        self.source_dict["SkyCoords"] = self.source_dict["SkyCoords"] + skycoords

        return self.source_dict

    @property
    def get_objects(self):
        return self.source_dict["Object"]

    @property
    def get_skycoords(self):
        return self.source_dict["SkyCoords"]
            