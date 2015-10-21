import sys

sys.path.append('../')


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

def build_star_record():
    pass


def create_stars_list(catfile):
    catalogrec = numpy.recfromcsv(catfile, delimiter=';', skiprows=0)

    table = Table(catalogrec)

    print(table)

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
def horiz2eq(az, alt, image_info=None, sidtime=None, lat=None, lon=None):
    '''
    Calculate equatorial coordinates for the given observation site
    and the given point in the sky
    The coordinates must be given in degrees or hours
    '''

    if lat is None:
        lat = image_info.latitude
    if lon is None:
        lon = image_info.longitude
    if sidtime is None:
        sidtime = image_info.sidereal_time

    # Sidereal Time to Local Sidereal Time
    sidtime = sidtime + lon / 15.

    lat = lat * np.pi / 180.
    lon = lon * np.pi / 180.
    sidtime = sidtime * np.pi / 12.
    az = az * np.pi / 180.
    alt = alt * np.pi / 180.

    _sindec = np.sin(alt) * np.sin(lat) + np.cos(alt) * \
        np.cos(lat) * np.cos(az)
    dec = np.arcsin(_sindec)
    _cosdec = np.cos(dec)
    _sinH = -np.sin(az) * np.cos(alt) / _cosdec
    _cosH = (np.sin(alt) - _sindec * np.sin(lat)) / (_cosdec * np.cos(lat))

    H = np.arctan2(_sinH, _cosH)
    ra = sidtime - H

    ra = (ra * 12. / np.pi) % 24
    dec = dec * 180. / np.pi

    return(ra, dec)

def eq2horiz(ra, dec, image_info=None, sidtime=None, lat=None, lon=None):
    '''
    Calculate horizontal coordinates for the given observation site
    and the given point in the sky.
    The coordinates must be given in degrees or hours
    '''

    if lat is None:
        lat = image_info.latitude
    if lon is None:
        lon = image_info.longitude
    if sidtime is None:
        sidtime = image_info.sidereal_time

    # Sidereal Time to Local Sidereal Time
    sidtime = sidtime + lon / 15.

    lat = lat * np.pi / 180.
    lon = lon * np.pi / 180.
    sidtime = sidtime * np.pi / 12.
    ra = ra * np.pi / 12.
    dec = dec * np.pi / 180.

    H = sidtime - ra

    _sina = np.sin(dec) * np.sin(lat) + np.cos(dec) * np.cos(lat) * np.cos(H)
    alt = np.arcsin(_sina)
    _cosa = np.cos(alt)
    _sinA = -np.sin(H) * np.cos(dec) / _cosa
    _cosA = (np.sin(dec) - np.sin(lat) * _sina) / (_cosa * np.cos(lat))
    az = np.arctan2(_sinA, _cosA)

    az = (az * 180. / np.pi) % 360
    alt = alt * 180. / np.pi

    return az, alt

def horiz2xy(azimuth, altitude, image_info, derotate=True):
    '''
    Return X,Y position in the image from azimuth/altitude horizontal coord.
    azimuth and altitude must be in degrees.
    '''

    if derotate == True and (image_info.latitude_offset != 0 or image_info.longitude_offset != 0):
        # We have a real azimuth and altitude coordinates. If the camera is not
        # pointing to the zenith, we need to derotate the image.
        print 'derotate'
        print image_info.latitude_offset
        ra_appa, dec_appa = horiz2eq(
            azimuth, altitude,
            image_info,
            lat=image_info.latitude,
            lon=image_info.longitude)

        azimuth, altitude = eq2horiz(
            ra_appa, dec_appa,
            image_info,
            lat=image_info.latitude - image_info.latitude_offset,
            lon=image_info.longitude - image_info.longitude_offset)

    Rfactor = image_info.radial_factor * \
        (180.0 / np.pi) * np.sqrt(2 * (1 - np.sin(altitude * np.pi / 180.0)))
    X = image_info.resolution[0] / 2 + image_info.delta_x -\
        Rfactor * \
        np.cos(azimuth * np.pi / 180.0 -
               image_info.azimuth_zeropoint * np.pi / 180.0)
    Y = image_info.resolution[1] / 2 + image_info.delta_y +\
        Rfactor * \
        np.sin(azimuth * np.pi / 180.0 -
               image_info.azimuth_zeropoint * np.pi / 180.0)
    return X, Y

if __name__ == '__main__':

    from astropy import units as u
    from astropy.coordinates import SkyCoord

    from astropy.coordinates import SkyCoord, EarthLocation, AltAz
    from astropy.coordinates import FK5
    from astropy.time import Time
    from astropy.table import Column
    import numpy
    from astropy.io import fits
    #from wcsaxes import WCS
    from astropy.wcs import WCS

    input_options = main(['-i', 'Johnson_V20130912_011709.fit.gz'])

    configs = configparser.SafeConfigParser(defaults=defaults)
    configs.read('config.ini')

    magfield = 'vamg'
    colorfields = ['uv', 'bv', 'rv', 'iv']

    min_altitude = 15
    max_magnitude = 5 # In V band, for the momment

    latitude = 40.450941
    longitude = -3.726065
    height = 667

    delta_x = -18.63912476
    delta_y = -31.06643504
    radial_factor = 14.19766968
    azimuth_zeropoint = 88.64589921

    for fname in input_options.fits_filename_list:
        print('basic calibration')

        #final = basic_calibration(fname, configs)

        #fits.writeto('test.fits', final, clobber=True)

        catfile = 'catalog.txt'
        print('load catalog from {}'.format(catfile))

        catalogrec = numpy.recfromcsv(catfile, delimiter=';', skiprows=0)

        table = Table(catalogrec)
        print table.colnames

        # Only the first 300
        table = table[:300]

        # Double stars
        double_mask = np.logical_not(np.char.isspace(table['doub']))

        # Variable stars
        variable_mask = np.logical_not(np.char.isspace(table['var']))

        # Bad Photoemetry
        bad_phot_mask = np.logical_not(np.char.isspace(table['bad_phot_']))

        # Incomplete photometry
        phot_nan_mask = np.zeros_like(table['vmag'], dtype='bool')
        for colname in ['uv', 'bv', 'rv', 'iv']:
            phot_nan_mask = np.logical_or(phot_nan_mask, np.isnan(table[colname]))

        # photometric stars are
        photo_mask = ~double_mask & ~variable_mask & ~bad_phot_mask & ~phot_nan_mask



        ntable = table[photo_mask]
        extreme_colors_mask = (ntable['bv'] < -1) & (ntable['bv'] > 2)
        bright_stars_mask = ntable['vmag'] <= max_magnitude

        table_photometric = ntable[~extreme_colors_mask & bright_stars_mask]

        print 'photometric', len(table_photometric)
        print 'full', len(table)


        mm = SkyCoord(ra=table['raj1950'], dec=table['dej1950'], unit=(u.hourangle, u.deg), frame=FK5(equinox='J1950'))

        mm2 = SkyCoord(ra=table['_raj2000'], dec=table['_dej2000'], unit=(u.hourangle, u.deg), frame=FK5(equinox='J2000'))


        mm1 = mm.transform_to(FK5(equinox='J2000'))
        # Correct

        # mm1 ~ mm2

        location = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=height*u.m)

        time = Time("2013-09-12 01:17:09")
        mm_altz = mm.transform_to(AltAz(obstime=time,location=location))
        mm2_altz = mm2.transform_to(AltAz(obstime=time,location=location))

        print mm_altz[4].alt.degree

        table.add_column(Column(mm_altz.alt.degree, name='alt'))
        table.add_column(Column(mm_altz.az.degree, name='az'))
        table.add_column(Column(mm.dec.degree, name='dec'))
        table.add_column(Column(mm.ra.degree, name='ra'))

        visibility_mask = mm_altz.alt.degree > min_altitude

        visible_catalog = table[visibility_mask]
        print visible_catalog['recno', 'commonname', 'dec', 'ra', 'alt', 'az'][:5]
        print len(visible_catalog)
        #print mm_altz.secz
        used_filter = 'Johnson_V'

        class ImageInfo3(object):
            latitude = latitude
            longitude = longitude

            latitude_offset = -0.64102456
            longitude_offset = 0.67447422
            # Scale and rotation
            radial_factor = radial_factor
            azimuth_zeropoint = azimuth_zeropoint
            sidereal_time = 5 # Seems irrelevant?
            resolution = [2500, 2500] # There are NAXIS1, NAXIS2 from images
            delta_x = delta_x
            delta_y = delta_y

        imageinfo3 = ImageInfo3()

        print horiz2xy(300.29175212, 25.5948991205, imageinfo3, derotate=False)


        cdelt = 1.0 / radial_factor

        refx = 2500.0 / 2 + delta_x
        refy = 2500.0 / 2 + delta_y

        w = WCS(naxis=2)
        w.wcs.crpix = [refx, refy]
        w.wcs.cdelt = [cdelt, cdelt]
        w.wcs.crval = [90, 0]
        #w.wcs.crval = [85, 0]
        w.wcs.ctype = ["DEC--ZEA", "RA---ZEA"]
        w.wcs.lonpole = 180-azimuth_zeropoint
        w.wcs.latpole = 90.0

        print w.all_world2pix([[25.5948991205, 300.29175212]], 1)
