"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
This is the extension of astropy.ccdproc.ImageFileCollection to fit the needs for photozpy
"""

from ccdproc import ImageFileCollection
from pathlib import Path

class mImageFileCollection():
    
    def __init__(self, image_dir = None, keywords = None,
                 find_fits_by_reading = False,
                 filenames = None, glob_include = None, glob_exclude = None, ext = 0):
        
        """
        image_dir : list, str or pathlib.Path, optional
            The directories that have the image file. 
        keywords : list of str, "*", optional
            Keywords that should be used as column headings in the summary table. If the value is or includes ‘*’ then all keywords that appear in any of the FITS headers of the files in the collection become table columns. Default value is ‘*’ unless info_file is specified. Default is `None`.
        find_fits_by_reading : bool, optional
            If `True`, read each file in location to check whether the file is a FITS file and include it in the collection based on that, rather than by file name. Compressed files, e.g. image.fits.gz, will NOT be properly detected. Will be ignored if `filenames` is not `None`.
        filenames : str, list of str, optional
            List of the names of FITS files which will be added to the collection. The filenames may either be in location or the name can be a relative or absolute path to the file. Default is `None`.
        glob_include: str, optional
            Unix-style filename pattern to select filenames to include in the file collection. Can be used in conjunction with `glob_exclude` to easily select subsets of files in the target directory. Default is `None`.
        glob_exclude : str, optional
            Unix-style filename pattern to select filenames to exclude from the file collection. Can be used in conjunction with `glob_include` to easily select subsets of files in the target directory. Default is `None`.
        ext: str or int, optional
            The extension from which the header and data will be read in all files.Default is `0`.
        """
        
        # standarize the directory list
        if isinstance(image_dir, str):
            self.image_dir = [Path(image_dir)]
            self.ncollection = 1

        elif isinstance(image_dir, Path):
            self.image_dir = [image_dir]
            self.ncollection = 1

        elif isinstance(image_dir, list):
            self.image_dir = [Path(i) for i in image_dir]
            self.ncollection = len(self.image_dir)

        
        # create image collection(s)
        if self.ncollection == 1:
            self.mcollection = ImageFileCollection(location = self.image_dir[0], keywords = keywords, 
                                                   find_fits_by_reading = find_fits_by_reading, 
                                                   filenames = filenames, glob_include = glob_include, glob_exclude = glob_exclude, ext = 0)
            self.mcollection = [self.mcollection]
        
        elif self.ncollection > 1:
            
            self.mcollection = []
            
            for i in self.image_dir:

                _mcollection = ImageFileCollection(location = i, keywords = keywords,
                                                   find_fits_by_reading = find_fits_by_reading, 
                                                   filenames = filenames, glob_include = glob_include, glob_exclude = glob_exclude, ext = 0)
                
                self.mcollection += [_mcollection]

    def __getitem__(self, index):

        if self.ncollection == 1:
            if index != 0:
                raise IndexError("There is only one image collection!")
            else:
                return self.mcollection[0]

        else:
            return self.mcollection[index]
                    

        
            