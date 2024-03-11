"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
This is the extension of astropy.ccdproc.ImageFileCollection to fit the needs for photozpy
"""

from ccdproc import ImageFileCollection
from pathlib import Path

class mImageFileCollection(ImageFileCollection):
    
    def __init__(self, parent_location = None, keywords = None,
                 find_fits_by_reading = None,
                 filenames = None, glob_include = None, glob_exclude = None, ext = 0,
                 sub_dir = None):
        
        """
        sub_dir: list of pathlib.Path or str; the sub directories contained in the parent location that holds the images
        """
        
        self.parent_location = Path(parent_location)

        
        if sub_dir == None:
            # if the sub directory is not define or the list length is one, there is no tree strcutures for the images.
            # it's essentially the same as using ccdproc.ImageFileCollection
            self.mlocation = self.parent_location  # use self.mlocation to distinguish from the self.location of the ImageFileCollection class
            self.mcollection = super().__init__(location = self.mlocation, keywords = keywords, 
                                                find_fits_by_reading = find_fits_by_reading, 
                                                filenames = filenames, glob_include = glob_include, glob_exclude = glob_exclude, ext = 0)

        elif isinstance(sub_dir, str) or isinstance(sub_dir, Path):
            self.sub_dir = Path(sub_dir)
            self.mlocation = self.parent_location / self.sub_dir
            self.mcollection = super().__init__(location = self.mlocation, keywords = keywords, 
                                                find_fits_by_reading = find_fits_by_reading, 
                                                filenames = filenames, glob_include = glob_include, glob_exclude = glob_exclude, ext = 0)
            
        elif isinstance(sub_dir, list):
            if len(sub_dir) == 1:
                self.sub_dir = Path(sub_dir[0])
                self.mlocation = self.parent_location / self.sub_dir
                self.mcollection = super().__init__(location = self.mlocation, keywords = keywords, 
                                                    find_fits_by_reading = find_fits_by_reading, 
                                                    filenames = filenames, glob_include = glob_include, glob_exclude = glob_exclude, ext = 0)

            elif len(sub_dir) > 1:
                self.sub_dir = []
                self.mlocation = []
                self.mcollection = []
                for i in sub_dir:
                    _sub_dir = Path(i)
                    self.sub_dir += [_sub_dir]

                    _m_location = self.parent_location / _sub_dir
                    print(_m_location)
                    self.mlocation += [_m_location]

                    _mcollection = super().__init__(location = _m_location, keywords = keywords,
                                                    find_fits_by_reading = find_fits_by_reading, 
                                                    filenames = filenames, glob_include = glob_include, glob_exclude = glob_exclude, ext = 0)
                    
                    self.mcollection += [_mcollection]

    def __getitem__(self, index):

        if isinstance(self.mcollection, ImageFileCollection):
            if index != 0:
                raise IndexError("There is only one image collection!")
            else:
                return self.mcollection

        elif isintance(self.mcollection, list):
            return self.mcollection[index]
                    

        
            