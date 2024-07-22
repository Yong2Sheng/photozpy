from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
from astropy.visualization import astropy_mpl_style, simple_norm
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from astropy.io import fits
from astropy.nddata import CCDData
from astropy import units as u
import os
from pathlib import Path
import shutil
import gzip
import numpy as np
from regions import PixCoord, CirclePixelRegion, CircleSkyRegion, Regions
from astroquery.hips2fits import hips2fits
from matplotlib.colors import Colormap
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import LogStretch

def get_alain_image(image_path = None, hdu = 1, wcs = None, save_dir = None, sky_region = None, min_cut = 0.5, max_cut = 99.5, logstretch = 10, source_name = None, hips = "CDS/P/DSS2/blue", save_image = True):

    """
    Get the corresponding aladin image from wcs or the image file that contains wcs.
    """

    if wcs is None:
        wcs = CCDData.read(image_path, hdu = hdu).wcs

    
    result = hips2fits.query_with_wcs(hips=hips,
                                      wcs=wcs,
                                      get_query_payload=False,
                                      format='fits',
                                      min_cut=0.5,
                                      max_cut=99.5,
                                      cmap=Colormap('viridis'))

    if save_image is True:
        if save_dir is None:
            save_dir = Path(image_path).parent

        x_npix = result[0].data.shape[0]
        y_npix = result[0].data.shape[1]
    
        tick_fontsize = 15
        subplot_labelsize = 15
        subtitle_fontsize = 15
        
        fig = plt.figure(figsize = (15, 15))
        if source_name is not None:
            fig.suptitle(f"Region summary for {source_name}", y=0.92, fontsize = 20)
        norm = ImageNormalize(stretch=LogStretch(logstretch))
    
        # plot the full image
        ax0 = fig.add_subplot(211, projection=wcs)
        ax0.imshow(result[0].data, origin='lower', norm=norm, cmap='Greys_r', interpolation='nearest')
    
        ax0.coords.grid(True, color='white', ls='dotted')
        ax0.coords[0].set_axislabel('Right Ascension (J2000)', fontsize=tick_fontsize)
        ax0.coords[1].set_axislabel('Declination (J2000)', fontsize=tick_fontsize)
        ax0.tick_params(which='major', labelsize=tick_fontsize)
        ax0.legend(fontsize = subplot_labelsize)
    
        if sky_region is None:
            fig.savefig(save_dir / "Aladin_image.png", dpi = 300, bbox_inches = "tight")
        else:
            # if sky region is given, it will plot the sky region in the full and zoomed image.
    
            # plot the region in the full image
            pixereg = sky_region.to_pixel(wcs)
            pixereg.plot(ax = ax0, color='red', lw=0.5, label = "Source region")
    
            # plotted the zoomed image
            ax1 = fig.add_subplot(212, projection=wcs)
            ax1.imshow(result[0].data, origin='lower', norm=norm, cmap='Greys_r', interpolation='nearest')
            pixereg.plot(ax = ax1, color='red', lw=0.5, label = "Zoomed source region")
            ax1.set_xlim(pixereg.center.x - x_npix*0.1, pixereg.center.x + x_npix*0.1)
            ax1.set_ylim(pixereg.center.y - y_npix*0.1, pixereg.center.y + y_npix*0.1)
    
            for ax in [ax0, ax1]:
                ax.coords.grid(True, color='white', ls='dotted')
                ax.coords[0].set_axislabel('Right Ascension (J2000)', fontsize=tick_fontsize)
                ax.coords[1].set_axislabel('Declination (J2000)', fontsize=tick_fontsize)
                ax.tick_params(which='major', labelsize=tick_fontsize)
                ax.legend(fontsize = subplot_labelsize)
    
            fig.savefig(save_dir / "Aladin_image.png", dpi = 300, bbox_inches = "tight")
            plt.close()

    return result

def estimate_background(image_path = None, array_data = None, hdu = 0, unit = None, sigma = 3, box_size = (20, 20), filter_size = (5, 5), **kwargs):

    """
    Estimate the 2D background of the fits image or array data by `photutils.background.MedianBackgroun` with sigma clip on the data.

    Parameters
    ----------
    image_path : str or pathlib.Path, optional
        The path to the fits file. It will overwrite `array_data` if both are provided.
    array_data : numpy.ndarray,optional
        The array of the image data.
    hdu : int, optional
        The hdu index of the image data in the fits file.
    unit : str, optional
        The unit of the image array data.
    sigma : float, optional
        The number of standard deviations to use for both the lower
        and upper clipping limit. These limits are overridden by
        ``sigma_lower`` and ``sigma_upper``, if input. The default is 3.
    box_size : int or array_like (int), optional
        The box size along each axis. If ``box_size`` is a scalar then
        a square box of size ``box_size`` will be used. If ``box_size``
        has two elements, they must be in ``(ny, nx)`` order. For best
        results, the box shape should be chosen such that the ``data``
        are covered by an integer number of boxes in both dimensions.
        When this is not the case, see the ``edge_method`` keyword for
        more options.
    **kwargs
        The parameters passed to `photutils.background.Background2D`.

    Returns
    -------
    estimated_bkg : photutils.background.background_2d.Background2D
        The estimated background array.
    array_data_no_bkg : numpy.ndarray
    """

    if image_path is not None:
        ccddata = CCDData.read(image_path, hdu = hdu, unit = unit)
        array_data = ccddata.data
    elif array_data is None:
        raise TypeError("You must provide fits path or array data!")

    sigma_clip = SigmaClip(sigma = sigma)
    bkg_estimator = MedianBackground()
    estimated_bkg = Background2D(array_data, box_size = box_size, filter_size = filter_size, 
                                 sigma_clip=sigma_clip, bkg_estimator=bkg_estimator, **kwargs)
    
    array_data_no_bkg = array_data - estimated_bkg.background

    return estimated_bkg, array_data_no_bkg

def ungz_file(file_path, out_dir = None, return_path = False, overwrite = True):

    """
    Ungz file and change the extension from img to fits.

    Paremeters
    ----------
    file_path : str or pathlib.Path
        The dir of the file to be 
    """

    gz_path = Path(file_path)
    ungzed_path = gz_path.with_suffix("").with_suffix(".fits")

    if out_dir is not None:
        ungzed_path = out_dir / ungzed_path.name
    
    # check if the ungzed file exists
    if ungzed_path.exists():
        if overwrite is True:
            os.remove(ungzed_path)
        else:
            raise OSError(f"File alreay exists: {ungzed_path}.")

    gz_file = gzip.GzipFile(gz_path)
    open(ungzed_path, "wb+").write(gz_file.read())
    gz_file.close()

    if return_path == False:
        return
    else:
        return ungzed_path


def create_folder(folder_dir):
    """
    Check if the folder exists and create the folder if not.

    Parameters
    ----------
    folder_path : str
        The path of the folder you want to create.

    Returns
    -------
    path_exist : boolean
        False if the folder_path doesn't exist;
        True if the folder_path already exist.
    """

    if os.path.isdir(folder_dir) != True:
        os.mkdir(folder_dir)
        path_exist = False
    else:
        path_exist = True
        
    return path_exist

def del_then_create_folder(folder_dir):
    """
    Delete folder if exists and then create the folder.

    Parameters
    ----------
    folder_dir : pathlib.path or str
        The directory of the folder
    """

    folder_dir = Path(folder_dir)
    
    if folder_dir.exists():
        shutil.rmtree(folder_dir)
    folder_dir.mkdir()
    
    return folder_dir


def convert_coords(image_path = None, wcs = None, skycoords = None, pixelcoords = None, verbose = False):

    """
    Takes a fits file or the astropy.wcs.WCS object as input. 
    Convert the sky coordinates to the pixel coordinates or vice versa.

    Parameters
    ----------
    file: str or path.Path object; the path to the object
    wcs: the astropy.wcs.WCS object. If both were input, wcs will cover the file.
    skycoords: astropy.Skycoord object. The sky coordinate of the object
    pixelcoords: 2D numpy array; the pixel coordinate of the object, the first column is the x pixels and the second is the y pixels: [[x pixles],[y pixels]]

    Returns
    -------
    astropy.Skycoords or list
    """

    # check if the number of input satistifies the calculation

    if image_path == None and wcs == None:
        raise TypeError("You must give a file path or asrtropy.wcs.WCS obkect as the input!")

    if skycoords == None and pixelcoords is None:
        raise TypeError("You must give sky coordinates or pixel coordinates as the input!")

    elif skycoords != None and pixelcoords is not None:
        raise TypeError("Please only input the sky coordinates or the pixel coordinates!")

    # Read the file and get the wcs object
    if wcs != None:
        wcs_object = wcs
    else:
        data = CCDData.read(image_path)#, unit = "adu")
        wcs_object = data.wcs

    if skycoords != None and pixelcoords is None:
        pixelcoords = wcs_object.world_to_pixel(skycoords)
        pixelcoords = np.array((pixelcoords)).T # transfer the array so it's ra/dec in each column
        out = pixelcoords # output variable
        if verbose == True:
            print("Conversion from sky coordiantes to pixel coordinates completed!")

    elif skycoords == None and pixelcoords is not None:
        xpixel_coords = pixelcoords[:,0]
        ypixel_coords = pixelcoords[:,1]
        radec = wcs_object.pixel_to_world(xpixel_coords, ypixel_coords)
        out = radec
        if verbose == True:
            print("Conversion from pixel coordinates to sky coordinates completed!")

    return out


def plot_image(fits_path = None,
               ccd_data = None, 
               save_location = None,
               skycoords = None, 
               pixelcoords = None, 
               circular_apertures = None, 
               annulus_apertures = None,
               adjust_fov = True,
               norm_percent = 99.9,
               fname_append = False):

    """
    Plot a image with skycoord marks and apertures.

    Parameters
    ----------
    image_path : str or pathlib.Path
        The path to 
    skycoords : NoneType or astropy.coordinates.SkyCoord, optional
        The sky coordinates to the plotted (the default is None, which implies no sky coordinates will be plotted).
    pixelcoords : NoneType or np.ndarray, optional
        The pixel coordinates to the plotted (the default is None, which implies no pixel coordinates will be plotted).
    circular_apertures : NoneType or photutils.aperture.CircularAnnulus, optional
        The circular apertures to be plotted (the default is None, which implies no circular apertures will be plotted).
    annulus_apertures : NoneType or photutils.aperture.CircularAperture, optional
        The annulus apertures to be plotted (the default is None, which implies no annulus apertures will be plotted)
    adjust_fov : bool, default=True
        When full_fov is True, the fov will be adjusted to cut off the unrelated regions to emphasize the sky coordinates and apertures.

    Returns
    -------
    None
    """

    plt.style.use(astropy_mpl_style)
    matplotlib.use('Agg')

    # initiate the array of centers
    centers_all = np.array([[0 ,0]])  # [0,0] is just used to initiate a 2d array
    
    if fits_path != None:
        # Get file names
        fits_path = Path(fits_path)  # path of image file
        if fname_append:
            image_name = fits_path.stem + f"{fname_append}" + ".png"
        else:
            image_name = fits_path.stem + ".png"
        if save_location == None:
            image_path = fits_path.parent/image_name
        else:
            image_path = Path(save_location)/image_name
    
        # Read fits data
        fits_data = fits.getdata(fits_path, ext=0)

    elif ccd_data != None:
        fits_data = ccd_data.data
        if fname_append:
            image_name = ccd_data.header["OBJECT"] + "_" + ccd_data.header["FILTER"] + fname_append + ".png"
        else:
            image_name = ccd_data.header["OBJECT"] + "_" + ccd_data.header["FILTER"] + fname_append + ".png"
        if save_location == None:
            image_path = Path(os.getcwd())/image_name
        else:
            image_path = Path(save_location)/image_name
    else:
        raise ValueError("You must provide either fits file path or CCDData!")

    norm = simple_norm(fits_data, 'log', percent= norm_percent) # define stretched norm
    
    # Plot the image
    plt.figure()
    plt.imshow(fits_data, norm=norm, cmap='gray')
    plt.colorbar()

    # Plot the pixel coordinates converted from sky coordinates
    if skycoords != None:
        # convert sky coordinates to pixel coordinates
        pixelcoords_from_skycoords = convert_coords(fits_path, skycoords = skycoords)
        centers_all = np.concatenate((centers_all, pixelcoords_from_skycoords), axis=0)
        
        plt.scatter(pixelcoords_from_skycoords[:,0], pixelcoords_from_skycoords[:,1], marker = "*", s=5, color = "lime", label = "True Coordinates")
        
        for row_idx in np.arange(skycoords.shape[0]):
            x_ = pixelcoords_from_skycoords[row_idx][0]
            y_ = pixelcoords_from_skycoords[row_idx][1]
            ra = skycoords[row_idx].ra.to_string(unit=u.hourangle, sep=':')
            dec = skycoords[row_idx].dec.to_string(unit=u.degree, sep=':')
            plt.annotate(f"RA={ra}, Dec={dec}", (x_, y_), color = "lime", size = 10)
            
    if isinstance(pixelcoords, np.ndarray):
        centers_all = np.concatenate((centers_all, pixelcoords), axis=0)
        plt.scatter(pixelcoords[:,0], pixelcoords[:,1], marker = "*", s=5, color = "lime", label = "True Coordinates")
        for row_idx in np.arange(pixelcoords.shape[0]):
            x_ = pixelcoords[row_idx][0]
            y_ = pixelcoords[row_idx][1]
            plt.annotate(f"x={x_}, y={y_}", (x_, y_), color = "lime", size = 10)

    if circular_apertures != None:
        centers_all = np.concatenate((centers_all, circular_apertures.positions), axis=0)
        ap_patches = circular_apertures.plot(color='orange', lw=1, label='Photometry aperture')
        circular_centers = circular_apertures.positions
        # plot the centers
        for row_idx in np.arange(circular_centers.shape[0]):
            x_ = circular_centers[row_idx][0]
            y_ = circular_centers[row_idx][1]
            plt.scatter(x_, y_, marker = "x", s=5, color = "orange")
            #plt.annoate(f"x={x_}, y={y_}", (x_, y_))

    if annulus_apertures != None:
        ann_patches = annulus_apertures.plot(color='red', lw=1, label='Background annulus')

    # Determine the final size of image (xlim and ylim)
    if skycoords == None and pixelcoords == None and circular_apertures == None and annulus_apertures == None:
        plt.savefig(image_path, dpi=300)
        plt.clf()
        plt.close("all")
    elif adjust_fov == False:
        plt.savefig(image_path, dpi=300)
        plt.clf()
        plt.close("all")
    else:
        data_xmax = fits_data.shape[0]
        data_ymax = fits_data.shape[1]

        # get plotting range of x
        xmin = centers_all[1:,0].min() - 100  # need to exlude the first row of [0,0]
        if xmin < 0:
            xmin = 0
        xmax = centers_all[1:,0].max() + 100
        if xmax > data_xmax:
            xmax = data_xmax

        # get plotting range of y
        ymin = centers_all[1:,1].min() - 100
        if ymin < 0:
            ymin = 0
        ymax = centers_all[1:,1].max() + 100
        if ymax > data_ymax:
            ymax = data_ymax

        plt.xlim(xmin, xmax)
        plt.ylim(ymax, ymin)
        plt.savefig(image_path, dpi=300)
        plt.clf()
        plt.close("all")

    return
