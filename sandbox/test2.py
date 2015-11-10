__author__ = 'spr'

import logging

import ConfigParser as configparser

import numpy as np
import numpy
from astropy.io import fits
from astropy.table import Table

from pyasb.user import main
from pyasb.astrometry import xy2horiz
from pyasb.star_calibration import Star

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
        delta_x = -18.63912476
        delta_y = -31.06643504
        radial_factor = 14.19766968
        azimuth_zeropoint = 88.64589921


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

def build_star_record():
    pass


def create_stars_list(catfile):
    catalogrec = numpy.recfromcsv(catfile, delimiter=';', skiprows=0)

    table = Table(catalogrec)

    # max_star_number is a parameter
    stars_tot = process_catalog_general(catalogrec, max_star_number=300)

    print(" - Total stars: %d" % len(stars_tot))



def process_catalog_general(catalogrec, max_star_number=-1):
    """
    Returns the processed catalog with
    all the starts that should be visible.
    """

    if max_star_number < 0:
        max_star_number = len(catalogrec)
        print 'Using default for ', max_star_number

    stars_tot = []

    class OtherImageInfo(object):
        used_filter = 'Johnson_V'
        date_string = "2013/09/12 01:17:09"
        #
        latitude = 40.450941
        longitude = -3.726065


    image_info = OtherImageInfo()

    for each_rec in catalogrec[:max_star_number]:
        star = Star(each_rec, image_info)
        stars_tot.append(star)

    return stars_tot




defaults = {
    'biasname': '',
    'flatfield': '',
    'darkframe': ''
}

if __name__ == '__main__':

    input_options = main(['-i', 'Johnson_V20130912_011709.fit.gz'])

    configs = configparser.SafeConfigParser(defaults=defaults)
    configs.read('config.ini')


    for fname in input_options.fits_filename_list:
        print('basic calibration')

        final = basic_calibration(fname, configs)

        fits.writeto('test.fits', final, clobber=True)

        catfile = 'catalog.txt'
        print('load catalog from {}'.format(catfile))

        calib_stars = create_stars_list(catfile)