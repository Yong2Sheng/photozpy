"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function:
- Perform aperture photometry
"""
from pathlib import Path
from ..collection_manager import CollectionManager
from astropy.nddata import CCDData
from photutils.aperture import CircularAperture, ApertureStats, aperture_photometry, SkyCircularAperture, SkyCircularAnnulus
from astropy.io import fits
from astropy.units import UnitsError
from astropy.stats import SigmaClip
import numpy as np
import re
from ..mimage_collection import mImageFileCollection
from ccdproc import ImageFileCollection
import pandas as pd
import os
from astropy.coordinates import concatenate
import astropy.units as u
from regions import CircleSkyRegion, Regions, CircleAnnulusSkyRegion
from astropy.table import QTable
import logging
import warnings
logger = logging.getLogger(__name__)


class Photometry():

    def __init__(self, image_collection):

        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(
            image_collection, rescan=True)

    @staticmethod
    def read_fwhm(image_path, keyword="FWHM"):
        """
        Read the FHWM from the header of the fits file

        Parameters
        ----------
        image_path: string or pathlib.Path; the path to the image
        keyword: string; the keyword name in the header.

        Returns
        -------
        fwhm: float; the read fwhm
        """

        headers = fits.getheader(image_path)
        fwhm = headers[keyword]

        return fwhm

    @staticmethod
    def get_aper_centroid(image_path, xy_coords, fwhm):
        """
        Estimate the centroid from aperture.

        Parameters
        Parameters
        ----------
        image_path: string or pathlib.Path; the path to the image
        xycen: numpy.ndarray or list; the physical coordinates of the source: each row: [x, y].

        Returns
        -------
        xycentroid: list; the x, y centroid position
        """
        # get the sigma_clipped median bkg
        data = CCDData.read(image_path)

        # fit the centroids first
        aper = CircularAperture(xy_coords, fwhm)
        aperstats = ApertureStats(data, aper)
        xycentroid = np.array([aperstats.xcentroid, aperstats.ycentroid]).T

        return xycentroid

    @staticmethod
    def estimate_mean_background(
            image_array_data, annulus_aperture, clip_sigma=3, tolerance=0.1):
        """
        Estimate the mean background level within a background annulus,
        using sigma clipping and a mean–median consistency check.

        Parameters
        ----------
        image_array_data : CCDData.data
            The 2D pixel array from a CCD image on which aperture
            statistics will be calculated.

        annulus_aperture : Aperture
            The annulus aperture defining the background region
            (e.g., a `CircularAnnulus`). This is passed to
            `ApertureStats` for background estimation.

        clip_sigma : float, optional
            Sigma level to use for iterative sigma clipping
            (default = 3). Passed to `SigmaClip`.

        tolerance : float, optional
            Consistency threshold between the mean and median
            background values, expressed as a fraction of the
            sigma-clipped per-pixel RMS of the annulus. If
            |mean − median| > tolerance × sigma_ann, a warning
            is raised (default = 0.1).

        Returns
        -------
        mean_bkg : float
            The sigma-clipped mean background per pixel within the annulus.

        bkg_stats : ApertureStats
            The full `ApertureStats` object containing all background
            statistics (mean, median, std, etc.). This can be used for
            further testing, debugging, or switching to the median later.

        Notes
        -----
        - This function does **not** automatically switch to the median
          background if the mean and median differ. It only raises a
          warning. Users may decide how to handle such cases.
        - The tolerance is expressed in units of sigma_ann (the clipped
          per-pixel RMS). Using a sigma-scaled tolerance makes this
          criterion robust across different sky levels and instruments.
        """

        sigclip = SigmaClip(sigma=clip_sigma, maxiters=10)
        bkg_stats = ApertureStats(
            image_array_data,
            annulus_aperture,
            sigma_clip=sigclip,
            sum_method="exact")

        mean_bkg = bkg_stats.mean
        median_bkg = bkg_stats.median
        sigma_ann = bkg_stats.std

        diff = np.abs(mean_bkg - median_bkg)
        bad = diff > tolerance * sigma_ann

        if np.any(bad):
            # get idx of the bkg annlus when differnece larger than tolerance
            bad_idx = np.where(bad)[0] if bad.ndim > 0 else np.array([0], dtype=int)

            lines = []
            lines.append(
                f"Background mean and median differ beyond tolerance for {bad_idx.size} annuli/annulus. "
                f"indices={bad_idx.tolist()}"
            )

            # write the warning text
            for i in bad_idx:
                i = int(i)
                lines.append(
                    "  "
                    f"annulus index={i}  "
                    f"mean={float(mean_bkg[i]):.4f}  "
                    f"median={float(median_bkg[i]):.4f}  "
                    f"diff={float(diff[i]):.4f}  "
                    f"sigma={float(sigma_ann[i]):.4f}  "
                    f"tol*sigma={float(tolerance * sigma_ann[i]):.4f}"
                )
    
            logger.warning("\n".join(lines))

        return mean_bkg * (u.adu / u.pix), bkg_stats

    @staticmethod
    def estimate_source_region_error(gain,                # e-/ADU
                                     readout_noise,       # e- (RMS per pixel)
                                     n_ap,                # effective aperture pixels (A_eff), float or Quantity[pix]
                                     n_ann,               # effective annulus pixels (A_eff), float or Quantity[pix]
                                     mean_bkg,            # mean background in ADU
                                     aperture_counts,     # total counts (source+sky) in aperture
                                     n_com,               # the number of exposures to stack the image, float
                                     ):
        r"""
        Estimate the error of the source aperture sum in the average-stacked image.

        It has three terms:
        1. Shot-noise variance from all electrons that actually landed inside the aperture (source + sky).
        2. Readout noise from reading the source aperture pixels.
        3. Error propagation from the background subtraction from the background annulus.

        The variance in the average-stacked image (in electrons) is：

        \sigma_{e,\mathrm{ap,sum}}^2 =
        \frac{1}{N_{\rm com}}\left[
            c_{e,\rm ap}
          + N_{\mathrm{ap}}\sigma_{\mathrm{readout}}^2
          + \frac{N_{\mathrm{ap}}^2}{N_{\mathrm{ann}}}
            \big(\hat{b}_{e,\mathrm{ann}} + \sigma_{\mathrm{readout}}^2\big)
        \right].


        Final error equation for backgroud-only case:
        \sigma^2_{e, \rm UL,avg} =
        \frac{1}{N_{\rm com}}\left[
            N_{\rm ap}\,\hat b_{e,\rm ann}  %shotnoise term
          + N_{\mathrm{ap}}\sigma_{\mathrm{readout}}^2
          + \frac{N_{\mathrm{ap}}^2}{N_{\mathrm{ann}}}
            \big(\hat{b}_{e,\mathrm{ann}} + \sigma_{\mathrm{readout}}^2\big)
        \right].


        For derivation, see Aperture Photometry Error Estimation and Upper Limits in SiYuan note.

        Parameters
        ----------
        gain : astropy.units.Quantity
            CCD gain in units of electron/ADU.
        readout_noise : astropy.units.Quantity
            CCD readout noise per pixel, in electrons.
        n_ap : astropy.units.Quantity
            Effective number of pixels in the source aperture (A_eff).
        n_ann : astropy.units.Quantity
            Effective number of pixels in the background annulus (A_eff).
        mean_bkg : astropy.units.Quantity
            Mean background in ADU/pix.
        aperture_counts : astropy.units.Quantity
            Total source+sky counts inside the aperture in ADU.

        """

        # units check is essential!
        if not gain.unit.is_equivalent(u.electron / u.adu):
            raise UnitsError(
                f"Gain should be in electron/ADU, got {gain.unit}")

        if not readout_noise.unit.is_equivalent(u.electron):
            raise UnitsError(
                f"Readout noise must be electrons (RMS per pixel), got {readout_noise.unit}")

        if not n_ap.unit.is_equivalent(u.pix):
            raise UnitsError(
                f"Source aperture pixel number should be in pix, got {n_ap.unit}")

        if not n_ann.unit.is_equivalent(u.pix):
            raise UnitsError(
                f"Background Annulus pixel number should be in pix, got {n_ann.unit}")

        if not mean_bkg.unit.is_equivalent(u.adu / u.pix):
            raise UnitsError(
                f"Background mean should be in adu/pix, got {mean_bkg.unit}")

        if not aperture_counts.unit.is_equivalent(u.adu):
            raise UnitsError(
                f"Total source+sky counts should be in adu, got {aperture_counts.unit}")

        if hasattr(n_com, "unit"):
            raise UnitsError("n_com must be dimensionless (number of exposures).")

        # Shot-noise variance
        # treat C_ap,e as a variance-like term (Poisson), attach e^-2 for
        # consistency when summing with read-noise variances.
        shot_noise_variance = (aperture_counts * gain).value * (u.electron**2)
        # shotnoise term for non-detection case
        shot_noise_variance_ul = (n_ap * mean_bkg * gain).value * (u.electron**2)

        # readout noise of reading the source aperture pixels
        # note here n_ap is just counting the number of readout noise to add.
        readout_noise_aperture = n_ap.value * readout_noise**2

        # background error propagation variance
        bkg_propagation = (n_ap.value**2 / n_ann.value * (mean_bkg.value *
                           gain.value + readout_noise.value**2)) * u.electron**2

        # final variance
        # variance for detection
        var_stacked = (shot_noise_variance + readout_noise_aperture + bkg_propagation) / n_com
        # variance for non-detection. The first turn changes to background only counts
        var_stacked_ul = (shot_noise_variance_ul + readout_noise_aperture + bkg_propagation) / n_com

        # error
        # after finishing the calculationes, we turn the error back to adu
        error_adu = np.sqrt(var_stacked) / gain
        error_ul_adu = np.sqrt(var_stacked_ul) / gain

        # double check if the unit conversion goes through automatically by itself.
        if not error_adu.unit.is_equivalent(u.adu):
            raise UnitsError(
                f"The final estimated error should be in adu, got {error_adu.unit}")

        # return shot_noise_variance, readout_noise, bkg_annulus_pix for future
        # debugging purpose
        return np.atleast_1d(error_adu), np.atleast_1d(error_ul_adu)

    @staticmethod
    def counts2mag(src_counts, error_counts, error_counts_ul,
                   detection_significance=3):
        r"""
        Convert source counts to instrumental magnitudes and propagate errors.

        The instrumental magnitude is defined as
        $$
        m_{\rm inst} = -2.5 \,\log_{10}(s) + ZP .
        $$
        In this function we compute *instrumental* magnitudes, so we set
        $ZP = 0$ by definition. Absolute zero points should be derived later
        using standard stars.

        Error propagation follows:
        $$
        \sigma_m^2
        &= \left(\frac{\partial m}{\partial s}\right)^2 \sigma_s^2\\

        &= \left(\frac{2.5}{\ln 10}\right)^2 \left(\frac{\sigma_{s}}{s}\right)^2,
        $$

        where:
          - $\sigma_s$ is the uncertainty on $s$ (in the same units as $s$);
          - $\sigma_m$ is the magnitude uncertainty.

        Note
        ----
        The ratio $\sigma_s / s$ is dimensionless, so the choice of units
        for `src_counts` (ADU vs electrons) does not affect the final
        magnitude error, as long as `src_counts`, `error_counts`, and
        `error_counts_ul` are expressed in the **same** units.

        This function also supports non-detections. If the measured
        significance falls below ``detection_significance`` (e.g.
        ``src_counts / error_counts < detection_significance`` or
        ``src_counts <= 0``), a magnitude is not reported and a
        background-only upper limit based on ``error_counts_ul`` should
        be used instead.

        Parameters
        ----------
        src_counts : astropy.units.Quantity
            Net source counts (ADU or electrons).
        error_counts : astropy.units.Quantity
            Uncertainty on the source counts (same units as `src_counts`).
        error_counts_ul:astropy.units.Quantity
            The background-only Uncertainty (same units as `src_counts`).
            It must be strictly positive for all non-detections.
        detection_significance : int or float
            The significance threshold (in sigma) to determine a detection.
            If `src_counts / error_counts < detection_significance`,
            the source is considered non-detected, and `error_counts_ul`
            will be used for upper limit calculation.

        Returns
        -------
        mag : Quantity
            Instrumental magnitude. If the source is not significantly detected,
            this may represent an upper limit magnitude calculated as
            `-2.5 * log10(k * error_counts_ul)`, where k is `detection_significance`.
        error_mag : Quantity
            Magnitude uncertainty. For detected sources, this is propagated from
            `error_counts`. For non-detections, this is typically set to NaN,
            as the magnitude value itself represents an upper limit.
        snr : numpy.ndarray
            Signal-to-noise ratio ``src_counts / error_counts``.
        """

        # make sure all the input are array-like instead of scalars
        src_counts = np.atleast_1d(src_counts)
        error_counts = np.atleast_1d(error_counts)
        error_counts_ul = np.atleast_1d(error_counts_ul)

        # make sure input are Quantity objects
        if not all(isinstance(x, u.Quantity) for x in (src_counts, error_counts, error_counts_ul)):
            raise TypeError("src_counts, error_counts, and error_counts_ul "
                            "must be astropy Quantity objects.")

        # We require *identical* units (not just equivalent values), because we are
        # working with raw counts (ADU or electrons).
        if not (src_counts.unit == error_counts.unit == error_counts_ul.unit):
            raise UnitsError("src_counts, error_counts, and error_counts_ul "
                             "must all have the same units."
                             f"(got src={src_counts.unit}, err={error_counts.unit}, ul={error_counts_ul.unit})")

        # make sure the input shapes are the same
        # if there is a global error_counts_ul or error_counts, it can be broadcasted in the future
        if not (src_counts.shape == error_counts.shape == error_counts_ul.shape):
            raise ValueError(f"Shapes must match: src={src_counts.shape}, "
                             f"err={error_counts.shape}, error_counts_ul={error_counts_ul.shape}")

        # # signal-to-noise ratio in counts-space
        snr = src_counts.value / error_counts.value
        nd_mask = (snr < detection_significance) | (src_counts.value <= 0)  # mask for the non-detection case

        # ignore warnings when src_counts <= 0
        with np.errstate(divide='ignore', invalid='ignore'):
            mag = -2.5 * np.log10(src_counts.value)
            error_mag = (2.5 / np.log(10)) * (error_counts.value / src_counts.value)

        mag[nd_mask] = -2.5 * np.log10(detection_significance * error_counts_ul.value[nd_mask])
        error_mag[nd_mask] = np.nan

        return np.atleast_1d(mag * u.mag), np.atleast_1d(error_mag * u.mag), np.atleast_1d(snr)

    @staticmethod
    def region2aperture(regions):

        if isinstance(regions[0], CircleSkyRegion):
            if len(regions) > 1:
                centers = concatenate([region.center for region in regions])
            else:
                # if the number of regions is one, you can't concatenate the
                # skycoords. Just get the region center
                centers = regions[0].center
            return SkyCircularAperture(centers, regions[0].radius)

        elif isinstance(regions[0], CircleAnnulusSkyRegion):
            if len(regions) > 1:
                centers = concatenate([region.center for region in regions])
            else:
                # if the number of regions is one, you can't concatenate the
                # skycoords. Just get the region cente
                centers = regions[0].center
            return SkyCircularAnnulus(
                centers, r_in=regions[0].inner_radius, r_out=regions[0].outer_radius)

    def run_photometry(self, sources, bkg_clip_sigma=3,
                       src_detection_sigma=3, hdu=0, verbose=True):
        """
        Does the aperture photometry on the skycoords.

        Parameters
        ----------
        sources: include standard stars and targets
        aperture: float; the source aperture
        inner_annulus: float; the inner radius of the annulus
        outer_annulus: float; the outer radius of the annulus

        Returns
        -------
        None
        """

        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(
            self._image_collection, rescan=True)

        # work on the object iteratively
        for source in sources:

            source_name = source.source_name
            telescope = source.telescope

            # get the image collection to work
            collection_photometry = CollectionManager.filter_collection(self._image_collection,
                                                                        **{"IMTYPE": "Master Light", "OBJECT": source_name})
            image_list = collection_photometry.files_filtered(
                include_path=True)
            # print(image_list)

            # initialize mag and error dict
            # we use [np.nan] avoid case when a source only has some of the filters
            # len(np.nan) returns an error, so I put it into a list
            mag_dict = {filter_name: [np.nan]
                        for filter_name in telescope.filters}
            mag_err_dict = {filter_name: [np.nan]
                            for filter_name in telescope.filters}
            significance_dict = {filter_name: [np.nan]
                                 for filter_name in telescope.filters}

            for image_path in image_list:
                image_path = Path(image_path)
                ccddata = CCDData.read(image_path, hdu=hdu)
                image_headers = ccddata.header
                image_filter_name = image_headers["FILTER"]
                gain = image_headers["GAIN"] * u.electron / u.adu
                readout_noise = image_headers["RDNOISE"] * u.electron
                n_com = image_headers["NCOMBINE"]

                print(
                    f"Working on photometry of {source_name} in {image_filter_name} from {image_path.name}"
                )

                if image_filter_name not in telescope.filters:
                    raise ValueError(
                        "The image filter is not in the filters of the telescope defined in sources!")
                image_wcs = ccddata.wcs
                image_array_data = ccddata.data
                # print("-------------------------------------------------------------------------------------------------")

                # get the aperture and annulus aperture

                src_region_fname = image_path.parent / f"{source_name}_{image_filter_name}_src.reg"
                if not src_region_fname.exists():
                    raise OSError(f"{src_region_fname} not found!")
                else:
                    src_regions = Regions.read(src_region_fname, format='ds9')
                    src_apertures_sky = Photometry.region2aperture(src_regions)
                    src_apertures_pix = src_apertures_sky.to_pixel(image_wcs)
                    # src_apertures_pix.area is dimensionless so here I added
                    # unit
                    n_ap = src_apertures_pix.area * u.pix

                bkg_region_fname = image_path.parent / f"{source_name}_{image_filter_name}_bkg.reg"
                if not bkg_region_fname.exists():
                    raise OSError(f"{bkg_region_fname} not found!")
                else:
                    bkg_regions = Regions.read(bkg_region_fname, format='ds9')
                    bkg_regions_sky = Photometry.region2aperture(bkg_regions)
                    bkg_annulus_pix = bkg_regions_sky.to_pixel(image_wcs)

                # get the sigma_clipped background estimation for all the annulus apertures
                # Important! If the annulus aperture contains multiple annulus (standard star case), the returned bkgs will be an array
                # If the annulus aperture contains only one annulus (target
                # case), the returned bkgs will be a float
                mean_bkg, bkg_stats = Photometry.estimate_mean_background(
                    image_array_data=image_array_data, annulus_aperture=bkg_annulus_pix, clip_sigma=bkg_clip_sigma)
                # get the effective number of pixels in the annulus after sigma
                # clip; bkg_stats.sum_aper_area has unit pix^2
                n_ann_eff = bkg_stats.sum_aper_area / u.pix
                total_bkg = mean_bkg * src_apertures_pix.area * u.pix

                # perform aperture photometry
                phot_table = aperture_photometry(
                    ccddata.data, src_apertures_pix, method="exact")
                # rename it for clarity
                phot_table['aperture_sum'].name = "src+bkg"
                phot_table['src+bkg'].unit = u.adu

                # remove background from the aperture region
                phot_src = phot_table['src+bkg'] - total_bkg

                # estimate the error in adu
                error_adu, error_ul_adu = Photometry.estimate_source_region_error(gain=gain,
                                                                                  readout_noise=readout_noise,
                                                                                  n_ap=n_ap,
                                                                                  n_ann=n_ann_eff,
                                                                                  aperture_counts=phot_table['src+bkg'],
                                                                                  mean_bkg=mean_bkg,
                                                                                  n_com=n_com)

                # calculate the instrumental magnitude
                m_inst, m_inst_error, significance = Photometry.counts2mag(src_counts=phot_src,
                                                                           error_counts=error_adu,
                                                                           error_counts_ul=error_ul_adu,
                                                                           detection_significance=src_detection_sigma)

                # organize the Qtable
                # add the column for total background
                phot_table['bkg'] = total_bkg
                # add the column for bkg subtracted photometry
                phot_table['src'] = phot_src
                # phot_table['error'] = m_inst_error
                # add the column for instrumental magnitude
                phot_table['mag_inst'] = m_inst
                phot_table['mag_inst_error'] = m_inst_error
                phot_table["src_significance"] = significance

                phot_table.meta = {"object": source_name,
                                   "filter": image_filter_name}

                for colname in ["xcenter", "ycenter"]:
                    phot_table[colname].info.format = "%9.4f"

                for colname in ["src+bkg", "bkg", "src"]:
                    phot_table[colname].info.format = "%8d"

                # phot_table["error"].info.format = "%9.4f"

                phot_table["mag_inst"].info.format = "%4f"

                phot_table["mag_inst_error"].info.format = "%4f"

                mag_dict[image_filter_name] = phot_table["mag_inst"]
                mag_err_dict[image_filter_name] = phot_table["mag_inst_error"]
                significance_dict[image_filter_name] = phot_table["src_significance"]

                phot_table.pprint_all()
                self.phot_table = phot_table
            print(
                "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
            print("\n")

            # # replace np.nan with -99
            # for key, value in mag_dict.items():
            #     for idx, v in enumerate(value):
            #         if np.isnan(v.value):
            #             mag_dict[key][idx] = -99*u.mag

            #  # replace np.nan with -99
            # for key, value in mag_err_dict.items():
            #     for idx, v in enumerate(value):
            #         if np.isnan(v.value):
            #             mag_err_dict[key][idx] = -99*u.mag

            source.magnitudes.inst_mags = QTable(mag_dict)
            source.magnitudes.inst_mag_errors = QTable(mag_err_dict)
            source.magnitudes.detection_significance = QTable(
                significance_dict)

        return


class SwiftPhotometry():

    def __init__(self, image_collection, source_catalog=None):

        if isinstance(image_collection, ImageFileCollection):
            self._mcollection = mImageFileCollection(
                location=image_collection.location,
                filenames=image_collection.files)

        elif isinstance(image_collection, mImageFileCollection):
            self._mcollection = image_collection

        self._source_catalog = source_catalog

    @staticmethod
    def extract_mag(out_txt):
        """
        Extract the magnitude information from the output txt file
        """

        f = open(out_txt, "r")
        lines = f.readlines()
        f.close()

        for idx, line in enumerate(lines):
            if "AB system" in line:
                mag_idx = idx + 1
                mag_line = lines[mag_idx]
                if ">" in mag_line:
                    # One or more digits (\d+), optional period (\.?), zero or
                    # more digits (\d*).
                    mag = float(re.findall("\\d+\\.?\\d*", mag_line)[0])
                    mag_error = -99
                else:
                    # One or more digits (\d+), optional period (\.?), zero or
                    # more digits (\d*).
                    number_list = re.findall("\\d+\\.?\\d*", mag_line)
                    mag = float(number_list[0])
                    stat_err = float(number_list[1])
                    sys_err = float(number_list[2])
                    mag_error = round(np.sqrt(stat_err**2 + sys_err**2), 2)

        return mag, mag_error

    def UVOTPhotometry(self):
        """
        Perform Swift photometry on the UVOT filters.
        """

        # init the dataframe to store the results
        df_results = pd.DataFrame(
            columns=[
                'name',
                "UVW2",
                "UVW2_err",
                "UVM2",
                "UVM2_err",
                "UVW1",
                "UVW1_err",
                "U",
                "U_err",
                "B",
                "B_err",
                "V",
                "V_err"])

        for collection in self._mcollection:

            file_location = Path(collection.location)
            source_name_from_path = file_location.parts[-1].replace("_", " ")
            # only use the final summed fits files
            file_paths = collection.files_filtered(
                include_path=True, **{"SUMTYP": "FINAL"})

            # init the dict to store the photometry for a single source, which
            # will be appended to the main result data frame
            dict_new = {'name': source_name_from_path,
                        "UVW2": -99,
                        "UVW2_err": -99,
                        "UVM2": -99,
                        "UVM2_err": -99,
                        "UVW1": -99,
                        "UVW1_err": -99,
                        "U": -99,
                        "U_err": -99,
                        "B": -99,
                        "B_err": -99,
                        "V": -99,
                        "V_err": -99}

            for file_path in file_paths:

                # set up the file path, region path, and read the source name
                # from the path
                file_path = Path(file_path)  # file path
                headers = fits.getheader(file_path, ext=1)
                filter_name = headers["FILTER"]   # filter name
                source_region_path = file_path.parent / \
                    f"{filter_name}.reg"  # source region path
                bkg_region_path = file_path.parent / \
                    f"{filter_name}_bkg.reg"  # bkg region path
                source_name_from_fits = headers["OBJECT"]

                # check if the two source names match
                if source_name_from_path == source_name_from_fits:
                    pass
                else:
                    raise ValueError(
                        f"The source name from the path ({source_name_from_path}) doesn't match the one from the fits file ({source_name_from_fits})")

                print(
                    f"Working on the photometry of source {source_name_from_fits} in {filter_name} filter.")

                outfits_path = file_path.parent / f"{filter_name}_result.fits"
                outtxt_path = file_path.parent / f"{filter_name}.result.txt"

                # run the command
                print(
                    f"Running uvotsource image={file_path} srcreg={source_region_path} bkgreg={bkg_region_path} sigma=3 cleanup=y clobber=y outfile={outfits_path} | tee {outtxt_path} >/dev/null")
                os.system(
                    f"uvotsource image={file_path} srcreg={source_region_path} bkgreg={bkg_region_path} sigma=3 cleanup=y clobber=y outfile={outfits_path} | tee {outtxt_path} >/dev/null")

                # extract the magnitudes
                ab_mag, ab_err = SwiftPhotometry.extract_mag(outtxt_path)

                # write the magnitudes to the image file
                with fits.open(file_path, mode="update") as hdul:
                    hdul[0].header["AB_MAG"] = ab_mag
                    hdul[1].header["AB_MAG_ERR"] = ab_err
                    hdul[0].header["AB_MAG"] = ab_mag,
                    hdul[1].header["AB_MAG_ERR"] = ab_err
                    hdul.flush()

                # write the magnitudes to the dict

                dict_new[filter_name] = ab_mag
                dict_new[f"{filter_name}_err"] = ab_err

                print(
                    "---------------------------------------------------------------------------")

            df_new = pd.DataFrame(dict_new, index=[0])
            df_results = pd.concat([df_results, df_new], ignore_index=True)

        if self._source_catalog is not None:
            self._source_catalog = Path(self._source_catalog)
            final_catalog = pd.read_csv(self._source_catalog, sep=",")
            final_catalog = pd.merge(
                final_catalog, df_results, on='name', how="left")
            final_catalog.to_csv(
                Path("") /
                "mag_results.csv",
                index=False,
                mode="w")
        else:
            df_results.to_csv(
                Path("") /
                "mag_results.csv",
                index=False,
                mode="w")

        return
