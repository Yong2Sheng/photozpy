"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
This is the extension of astropy.ccdproc.ImageFileCollection to fit the needs for photozpy
"""

from ccdproc import ImageFileCollection
from pathlib import Path
from astropy.io import fits

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


    def filter_master_collection(self, header, value):

        """
        Filter the master collection based on the headers and values to get one image collection.

        Parameters
        ----------
        header : str
            The header name.
        value : str or list
            The corresponding value(s).

        Returns
        -------
        ccdproc.ImageFileCollection
            The filtered image collection.
        """

        pass

    
    def get_header_values(self, header_name, index = None):
        
        """
        Get the header values from an image collection. The repeated values will be removed.

        Paremeters
        ----------
        header_name
        headers

        Returns
        -------
        list
            The values of the header
        """

        return image_collection.values(header, unique = True)

    def refresh_collections(self, index = None):

        """
        Refresh the collections to reflect the changes of the fits files.

        Paremeters
        ----------
        index : int, optional
            The index of the collection to be refreshed. (the default is `None`, which means all the collections will be refreshed).

        """
        pass


    @staticmethod
    def _get_header_values(image_collection, headers, unique = True):

        """
        Get the header valus from the image collection. I wrote this because the ImageFileCollection.values() doesn't work!

        """

        if isinstance(headers, str):
            headers = [headers]

        file_paths = image_collection.files_filtered(include_path = True)

        values_list = []
        for file_path in file_paths:
            header_all = fits.getheader(file_path)
            for header in headers:
                value_ = header_all[header]
                values_list.append(value_)

        if unique:
            values_list = [*set(values_list)]

        return values_list

    def get_collection_header_values(self, headers, index = None, unique = True):

        if isinstance(headers, str):
            headers = [headers]

        file_paths = image_collection.files_filtered(include_path = True)

        if index is not None:
            _collection = self.mcollection[index]
            header_values = mImageFileCollection._get_header_values(_collection, headers = headers)

        else:
            header_values = []
            for i in np.arange(self.ncollection):
                _collection = self.mcollection[i]
                _header_values = mImageFileCollection._get_header_values(_collection, headers = headers)
                header_values += _header_values
                
        if unique:
            header_values = [*set(header_values)]

        return header_values
        
        
                    

        
            