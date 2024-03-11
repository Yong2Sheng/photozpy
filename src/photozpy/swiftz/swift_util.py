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

class Download():

    def __init__(self, download_dir):

        """
        Initialize the class.
        
        Parameters
        ----------
        root_dir: pathlib.Path or str; the directory where the data will be stored.

        Returns
        -------
        None
        """

        