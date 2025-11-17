"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function:
- apply the calibrations (master bias, master dark and master flat)
"""

from astroquery.astrometry_net import AstrometryNet
from astropy.wcs import WCS
from ..collection_manager import CollectionManager
from pathlib import Path
from astropy.io import fits
import subprocess


class PlateSolving():

    def __init__(self, image_collection, sources):

        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(
            image_collection, rescan=True)
        self._sources = sources

    def add_wcs_astrometrynet(self, api, image_type="Master Light", fwhm=None, detect_threshold=5,
                              solve_timeout=120, ra_header_name="RA", dec_header_name="DEC", ra_dec_units=("hour", "deg")):

        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(
            self._image_collection, rescan=True)

        image_collection = CollectionManager.filter_collection(
            self._image_collection, **{"IMTYPE": image_type})
        all_image_list = image_collection.files_filtered(include_path=True)

        # here I want to move the standard sources to the front of the list
        # So once one standard image is finished adding wcs, I can start deciding the standard stars and the standard magnitudes RIGHT AWAY!
        # I don't have to wait until the end of plate solving since it's quiet
        # time consuming.
        std_collection = CollectionManager.filter_collection(
            self._image_collection, **{
                "IMTYPE": "Master Light", "OBJECT": self._sources.standard_stars[0].source_name})
        std_list = std_collection.files_filtered(include_path=True)
        target_list = [i for i in all_image_list if i not in std_list]
        image_list = std_list + target_list

        for image in image_list:
            print(f"Adding WCS to {image}")
            if fwhm is None:
                # try to read the fwhm from the header
                with fits.open(image) as hdul:
                    fwhm = hdul[0].header["FWHM"]
            ast = AstrometryNet()
            ast.api_key = api
            wcs_header = ast.solve_from_image(
                image,
                fwhm=fwhm,
                detect_threshold=detect_threshold,
                solve_timeout=solve_timeout)  # , verbose = False)
            plate_solved_wcs_header = WCS(wcs_header).to_header()

            with fits.open(image, mode="update") as hdul:
                for card in plate_solved_wcs_header.cards:
                    hdul[0].header[card[0]] = (card[1], "UPDATED!!")
                hdul.flush()
            print(
                "-----------------------------------------------------------------------------------\n")

        return

    @staticmethod
    def add_wcs_locally_for_file(file_path, ra=None, dec=None, radius=10, detect_threshold=8,
                                 extra_args=None, clean=True):
        """
        Solve WCS headers using a locally installed Astrometry.net engine.

        Notes
        -----
        - The conda-forge package of Astrometry.net does not support ARM Macs,
          so Homebrew must be used on those systems.
        - You must download the appropriate index files and place them into the
          data directory of the local installation, e.g.,
          `/opt/homebrew/Cellar/astrometry-net/0.97/data`.
        - Example index download script (for a SARA telescope; adjust based on FOV):
            for j in 5203 5204 5205; do
              for ((i=0; i<48; i++)); do
                I=$(printf %02i $i)
                wget https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-$j-$I.fits
              done
            done

        Parameters
        ----------
        file_path : str or pathlib.Path
            Path to the FITS file to be solved.
        ra : str, optional
            Right Ascension (e.g., 'hh:mm:ss' or degrees) of the field centre.
            If None, the RA will be extracted from the FITS header.
        dec : str, optional
            Declination (e.g., '+dd:mm:ss' or degrees) of the field centre.
            If None, the Dec will be extracted from the FITS header.
        radius : float, optional
            Radius (in degrees) around (ra, dec) within which to search indexes.
            Default is 10°.
        detect_threshold : float, optional
            Number of sigma for source detection. Defaults to the internal default if None.
        extra_args : string, optional
            Extra args passed to solve_field. Defaults to the internal default if None.
        clean : bool, optional
            If True, delete the original input file and rename/replace it
            with the solved output (effectively overwriting the original).
            If False, keep both original and solved files. Default is True.


        Returns
        -------
        subprocess.CompletedProcess
            The result of the `solve-field` command including return code,
            stdout and stderr where available.

        Raises
        ------
        FileNotFoundError
            If the input file does not exist.
        RuntimeError
            If the WCS solving process fails (non-zero return code or missing WCS keywords).

        """

        file_path = Path(file_path)

        if not file_path.exists():
            raise FileNotFoundError(
                f"Input file does not exist: '{file_path}'")

        if ra is None:
            try:
                headers = fits.getheader(file_path, ext=0)
                ra = headers["RA"]
            except KeyError:
                pass

        if dec is None:
            try:
                headers = fits.getheader(file_path, ext=0)
                dec = headers["DEC"]
            except KeyError:
                pass

        # basic commands
        output_file = file_path.parent / (file_path.stem + "solved.fits")
        cmd = ["solve-field",
               str(file_path),
               "--overwrite",
               "--no-plots",
               "--new-fits", str(output_file),
               "--nsigma", str(detect_threshold)]

        if (ra is not None) and (dec is not None):
            cmd += ["--ra", str(ra), "--dec", str(dec),
                    "--radius", str(radius)]

        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True)

        # if 0, then solving succeeds; if not 0, the solving fails
        if result.returncode != 0:
            raise RuntimeError(f"solve-field failed for '{file_path}'. "
                               f"stdout={result.stdout!r}, stderr={result.stderr!r}")

        # clean log and intermediate files
        # delete the origina file and rename the solved file to the origian
        # file name
        if clean is True:
            patterns = [
                "*.xyls",
                "*.corr",
                "*.match",
                "*.rdls",
                "*.wcs",
                "*.solved",
                "*.axy"]
            for pattern in patterns:
                for f in file_path.parent.glob(pattern):
                    try:
                        f.unlink()
                    except Exception as e:
                        print(f"Failed to delete {f}: {e}")

            try:
                output_file.rename(file_path)
            except FileExistsError:
                output_file.replace(file_path)

        return result

    def add_wcs_locally(self, image_type="Master Light", detect_threshold=8):
        """
        Solve all the fits images in a ImageFileCollection filtered by the
        `image_type` (IMTYPE).

        Parameters
        ----------
        image_type : string, optional
            To filter the fits files with that value of IMTYPE header.
        detect_threshold : float, optional
            Number of sigma for source detection. Defaults to the internal default if None.

        Returns
        -------
        None

        Notes
        -----
        This function loops through the files in the internal ImageFileCollection,
        selects only those whose IMTYPE matches `image_type`, and for each file
        invokes the local Astrometry.net solver (via add_wcs_locally_for_file).
        Files that fail solving are logged, and processing continues.
        """

        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(
            self._image_collection, rescan=True)

        image_collection = CollectionManager.filter_collection(
            self._image_collection, **{"IMTYPE": image_type})
        all_image_list = image_collection.files_filtered(include_path=True)

        # here I want to move the standard sources to the front of the list
        # So once one standard image is finished adding wcs, I can start deciding the standard stars and the standard magnitudes RIGHT AWAY!
        # I don't have to wait until the end of plate solving since it's quiet
        # time consuming.
        std_collection = CollectionManager.filter_collection(
            self._image_collection, **{
                "IMTYPE": "Master Light", "OBJECT": self._sources.standard_stars[0].source_name})
        std_list = std_collection.files_filtered(include_path=True)
        target_list = [i for i in all_image_list if i not in std_list]
        image_list = std_list + target_list
        total_number = len(image_list)

        failure_list = []
        for count, image in enumerate(image_list, start=1):
            image_path = Path(image)
            try:
                _ = PlateSolving.add_wcs_locally_for_file(image_path)
            except FileNotFoundError as e:
                print(f"[{count}/{total_number}]:  {image_path.name} not found.")
                failure_list.append((image_path, e))
                continue
            except RuntimeError as e:
                print(
                    f"[{count}/{total_number}]: {image_path.name} WCS solving failed.")
                failure_list.append((image_path, e))
                continue
            else:
                print(
                    f"[{count}/{total_number}]: {image_path.name} WCS solved successfully.")

        # print failures if present
        if failure_list:
            print(f"{len(failure_list)} files failed:")
            for fp, err in failure_list:
                print(f"  {fp} → {err}")
        return
