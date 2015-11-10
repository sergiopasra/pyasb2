import logging
from astropy.io import fits
import numpy as np

from pyasb.astrometry import xy2horiz

def create_synthetic_dark(dark_data, bias_data, t_sci, t_dark):
    mdarkdata = dark_data - bias_data
    mdarkdata *= (t_sci / t_dark)
    mdarkdata += bias_data
    return mdarkdata


def basic_calibration(fname, configs):
    used_filter = 'Johnson_V'
    section = 'calibrations_%s' % used_filter

    do_dark = True
    do_flat = True
    substract_corners_background = True
    # Better a empty string ""?
    darkname = ''
    flatname = ''
    biasname = ''

    if configs.has_section(section):
        darkname = configs.get(section, 'darkframe')
        flatname = configs.get(section, 'flatfield')
        biasname = configs.get(section, 'biasname')
        if darkname == '':
            do_dark = False

        if flatname == '':
            do_flat = False
    else:
        # Not really interesting
        do_dark = False
        do_flat = False

    if not (do_dark or do_dark):
        # we stop here
        raise TypeError('no flat, no dark')

    class ImageInfo(object):
        resolution = [2500, 2500] # There are NAXIS1, NAXIS2 from images
        delta_x = delta_x
        delta_y = delta_y
        radial_factor = radial_factor
        azimuth_zeropoint = azimuth_zeropoint


    image_info = ImageInfo()

    msdata, msheader = fits.getdata(fname, header=True)
    t_exp_science = float(msheader['EXPOSURE'])
    # Horrible non standard format "%Y%m%d_%H%M%S"
    fits_date = msheader['DATE']
    import datetime

    date_array = datetime.datetime.strptime(fits_date, "%Y%m%d_%H%M%S")
    date_string = date_array.strftime("%Y/%m/%d %H:%M:%S")

    print(fits_date)
    print(date_array)
    print(date_string)

    fsdata = msdata[:]

    if do_dark:
        mdarkdata, mdarkheader = fits.getdata(darkname, header=True)
        t_exp_dark = float(mdarkheader['EXPOSURE'])
        # TODO: better comparation here
        if t_exp_dark != t_exp_science:
            print('trying to create a syn dark')
            print('we need a bias')
            if biasname:
                bdata, bheader = fits.getdata(biasname, header=True)
                mdarkdata = create_synthetic_dark(mdarkdata, bdata, t_exp_science, t_exp_dark)
            else:
                print('bias name is invalid')
                raise TypeError('bias invalid')

        print('corecting bias')
        fsdata = dark_correction(fsdata, mdarkdata)

    if substract_corners_background:
        print('subtract background in the corner')
        im_coor = image_coordinates(image_info)
        altitude_map = im_coor[1]
        fsdata, med, err = do_substract_corners_background(fsdata, altitude_map, limit=-20)
        print("Removed: {:.2f} +/- {:.2f} counts from measured background".format(med, err))
    if do_flat:
        logging.info('loading master flat')
        mflatdata, mflatheader = fits.getdata(flatname, header=True)

        print('corecting flat field')
        fsdata = flat_field_correction(fsdata, mflatdata)

    return fsdata


def dark_correction(data, masterdark_data):
    return data - masterdark_data


def flat_field_correction(data, masterflat_data):
    return data / masterflat_data


def do_substract_corners_background(array, altitude_map, limit=-20):
    data_corners = array[altitude_map < limit]
    bias_image_median = np.median(data_corners)
    bias_image_std = np.std(data_corners)
    bias_image_err = bias_image_std / np.sqrt(np.size(data_corners))
    newarray = array - bias_image_median
    return newarray, bias_image_median, bias_image_err

def image_coordinates(image_info):
    """Reimplementation with numpy arrays (fast on large arrays).
     We need it as we will use a very large array"""

    # Image coordinates
    x = np.arange(image_info.resolution[0])
    y = np.arange(image_info.resolution[1])
    xx, yy = np.meshgrid(x, y)

    # Unreal/projected altitude and azimuth
    az, alt = xy2horiz(xx, yy, image_info, derotate=False)
    azimuth_map = np.array(az, dtype='float32')
    altitude_map = np.array(alt, dtype='float32')
    return azimuth_map, altitude_map
