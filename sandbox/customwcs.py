import numpy as np

class CustomWCS(object):

    def __init__(self):
        self.latitude = 40.450941
        self.longitude = -3.726065

        self.latitude_offset = -0.64102456
        self.longitude_offset = 0.67447422
        # Scale and rotation
        self.delta_x = -18.63912476
        self.delta_y = -31.06643504
        self.radial_factor = 14.19766968
        self.azimuth_zeropoint = 88.64589921
        self.sidereal_time = 0.700320860922
        self.resolution = [2500, 2500] # These are NAXIS1, NAXIS2 from images
        self.radial_factor = 14.19766968
        self.cdelt = 1.0 / self.radial_factor

        self.refx = self.resolution[1] / 2 + self.delta_x
        self.refy = self.resolution[0] / 2 + self.delta_y

        self.lat1, self.lon1 = self.latitude, self.longitude
        self.lat2 = self.latitude - self.latitude_offset
        self.lon2 = self.longitude - self.longitude_offset

        self.lat1_d = self.lat1 * np.pi / 180.0
        self.lat2_d = self.lat2 * np.pi / 180.0
        self.lon1_d = self.lon1 * np.pi / 180.0
        self.lon2_d = self.lon2 * np.pi / 180.0
        self.difflon = self.lon2_d - self.lon1_d

    def horiz2xy(self, azimuth, altitude, derotate=True):
        '''
        Return X,Y position in the image from azimuth/altitude horizontal coord.
        azimuth and altitude must be in degrees.
        '''

        if derotate:
            # We have a real azimuth and altitude coordinates. If the camera is not
            # pointing to the zenith, we need to derotate the image.

            azimuth, altitude = self.derot_1s(azimuth, altitude)


        Rfactor = self.radial_factor * (180.0 / np.pi) * np.sqrt(2 * (1 - np.sin(altitude * np.pi / 180.0)))
        angle = (azimuth  - self.azimuth_zeropoint) * np.pi / 180.0
        xx = self.refx - Rfactor * np.cos(angle)
        yy = self.refy + Rfactor * np.sin(angle)
        return xx, yy

    def derot_1s(self, az, alt):

        az_d = az * np.pi / 180.
        alt_d = alt * np.pi / 180.

        az2, alt2 = self.derot_1s_rad(az_d, alt_d, self.lat1_d, self.lat2_d, self.difflon)

        az = (az2 * 180. / np.pi) % 360
        alt = alt2 * 180. / np.pi

        return az, alt

    def derot_1s_rad(self, az, alt, lat1, lat2, angle2):

        sindec = np.sin(alt) * np.sin(lat1) + np.cos(alt) * np.cos(lat1) * np.cos(az)
        cosdec = np.sqrt(1-sindec**2)
        sinH1 = -np.sin(az) * np.cos(alt) / cosdec
        cosH1 = (np.sin(alt) - sindec * np.sin(lat1)) / (cosdec * np.cos(lat1))

        cosH2 = cosH1 * np.cos(angle2) - sinH1 * np.sin(angle2)
        sinH2 = sinH1 * np.cos(angle2) + cosH1 * np.sin(angle2)

        _sina = sindec * np.sin(lat2) + cosdec * np.cos(lat2) * cosH2
        alt2 = np.arcsin(_sina)
        _cosa = np.cos(alt2)
        _sinA = -sinH2 * cosdec / _cosa
        _cosA = (sindec - np.sin(lat2) * _sina) / (_cosa * np.cos(lat2))
        az2 = np.arctan2(_sinA, _cosA)

        return az2, alt2

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

    if derotate:
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
    return X, Y
