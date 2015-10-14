#
# PyASB astrometry functions.
#
# Convert from one coordinate system to another.
#____________________________
#
# This module is part of the PyASB project,
# created and maintained by Miguel Nievas [UCM].
#____________________________
#

"""Astrometry with Pyephem"""

import numpy as np

import ephem

# Setup Pyephem Observatory


def pyephem_setup_common(image_info):
    obs_pyephem = ephem.Observer()
    obs_pyephem.pressure = 0  # Dont consider atmospheric effects
    obs_pyephem.date = image_info.date_string
    return obs_pyephem


def pyephem_setup_image(image_info):
    obs_pyephem = pyephem_setup_common(image_info)
    obs_pyephem.lat = (
        image_info.latitude - image_info.latitude_offset) * np.pi / 180
    obs_pyephem.lon = (
        image_info.longitude - image_info.longitude_offset) * np.pi / 180
    return obs_pyephem


def pyephem_setup_real(image_info):
    obs_pyephem = pyephem_setup_common(image_info)
    obs_pyephem.lat = (image_info.latitude) * np.pi / 180
    obs_pyephem.lon = (image_info.longitude) * np.pi / 180
    return obs_pyephem

# Standalone functions.
# To be used on single points


def horiz2xy_old(azimuth, altitude, image_info):
    '''
    Return X,Y position in the image from azimuth/altitude horizontal coord.
    azimuth and altitude must be in degrees.
    '''

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
    return(X, Y)


def horiz2xy(azimuth, altitude, image_info, derotate=True):
    '''
    Return X,Y position in the image from azimuth/altitude horizontal coord.
    azimuth and altitude must be in degrees.
    '''

    if derotate == True and (image_info.latitude_offset != 0 or image_info.longitude_offset != 0):
        # We have a real azimuth and altitude coordinates. If the camera is not
        # pointing to the zenith, we need to derotate the image.
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
    return(X, Y)


def xy2horiz(X, Y, image_info, derotate=True):
    '''
    Return horizontal coordinates from X,Y position in the image.
    azimuth and altitude are in degrees.
    '''

    X = X - image_info.resolution[0] / 2. - image_info.delta_x
    Y = Y - image_info.resolution[1] / 2. - image_info.delta_y
    rfactor = np.sqrt(X ** 2 + Y ** 2) / image_info.radial_factor

    if np.size(rfactor) > 1:
        rfactor[rfactor > 360. / np.pi] = 360. / np.pi

    alt_factor = np.array(1 - 0.5 * (np.pi * rfactor / 180.0) ** 2)
    altitude = (180.0 / np.pi) * np.arcsin(alt_factor)
    azimuth = (360 + image_info.azimuth_zeropoint +
               180.0 * np.arctan2(Y, -X) / np.pi) % 360

    if derotate and (image_info.latitude_offset != 0 or image_info.longitude_offset != 0):
        # We have a real azimuth and altitude coordinates. If the camera is not
        # pointing to the zenith, we need to rotate the image.

        ra_real, dec_real = horiz2eq(
            azimuth, altitude,
            image_info,
            lat=image_info.latitude - image_info.latitude_offset,
            lon=image_info.longitude - image_info.longitude_offset)

        azimuth, altitude = eq2horiz(
            ra_real, dec_real,
            image_info,
            lat=image_info.latitude,
            lon=image_info.longitude)

    return(azimuth, altitude)


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

    return(az, alt)


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


def zenith_position(image_info):
    # Return X,Y position of zenith in the image.
    return horiz2xy(0, 90, image_info)


def optical_axis(image_info):
    # Return horizontal coordinates of the optical axis
    return xy2horiz(image_info.resolution[0] / 2, image_info.resolution[1] / 2, image_info)


def atmospheric_refraction(altitude, mode):
    # Return apparent (non-corrected from refraction) or
    # real (corrected from refraction) altitude.
    # Garfinkel (1967), http://en.wikipedia.org/wiki/Atmospheric_refraction
    def cot(x):
        # Return cotangent of the given value
        return np.cos(x) / np.sin(x)

    if mode == 'dir':
        # Return apparent altitude from the real one.
        return altitude + (1.02 / 60) * cot(altitude + 10.3 / (altitude + 5.11))
    elif mode == 'inv':
        # Return real altitude from the apparent one
        return altitude - (1.00 / 60) * cot(altitude + 7.31 / (altitude + 4.4))
    else:
        print 'Unknow mode ' + mode + '. Cannot correct from atmospheric refraction.'
        return altitude


def calculate_airmass(altitude):
    # Estimate airmass from apparent altitude using Pickering (2002) model.
    # zdist in degrees
    return 1 / np.sin((altitude + 244. / (165 + 47 * altitude ** 1.1)) * np.pi / 180.)

# Vectorial functions.
# Generate a class that contains a map of coordinates
# that match the Image pixels


class ImageCoordinates(object):

    def __init__(self, image_info):
        self.calculate_altaz(image_info)

    def calculate_altaz(self, image_info):
        ''' Reimplementation with numpy arrays (fast on large arrays). 
            We need it as we will use a very large array'''

        # Image coordinates
        x = np.arange(image_info.resolution[0])
        y = np.arange(image_info.resolution[1])
        X, Y = np.meshgrid(x, y)

        # Unreal/projected altitude and azimuth
        az, alt = xy2horiz(X, Y, image_info, derotate=False)
        self.azimuth_map = np.array(az, dtype='float16')
        self.altitude_map = np.array(alt, dtype='float16')
