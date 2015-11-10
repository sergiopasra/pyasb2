import numpy as np

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

    return ra, dec

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

def horiz2eq_(az, alt, lat, lon):
    '''
    Calculate equatorial coordinates for the given observation site
    and the given point in the sky
    The coordinates must be given in degrees or hours
    '''

    lat = lat * np.pi / 180.
    local = lon * np.pi / 180.0
    az = az * np.pi / 180.
    alt = alt * np.pi / 180.

    _sindec = np.sin(alt) * np.sin(lat) + np.cos(alt) * np.cos(lat) * np.cos(az)
    dec = np.arcsin(_sindec)
    _cosdec = np.cos(dec)
    _sinH = -np.sin(az) * np.cos(alt) / _cosdec
    _cosH = (np.sin(alt) - _sindec * np.sin(lat)) / (_cosdec * np.cos(lat))

    H = np.arctan2(_sinH, _cosH)
    ra = local - H

    ra = (ra * 180. / np.pi) % 360
    dec = dec * 180. / np.pi

    return ra, dec


def eq2horiz_(ra, dec, lat, lon):
    '''
    Calculate horizontal coordinates for the given observation site
    and the given point in the sky.
    The coordinates must be given in degrees or hours
    '''

    # Sidereal Time to Local Sidereal Time
    sidtime = lon / 15.

    lat = lat * np.pi / 180.
    local = lon * np.pi / 180
    ra = ra * np.pi / 180.0
    dec = dec * np.pi / 180.

    H = local - ra

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
    angle = azimuth * np.pi / 180.0 - image_info.azimuth_zeropoint * np.pi / 180.0
    X = image_info.resolution[0] / 2 + image_info.delta_x - Rfactor * np.cos(angle)
    Y = image_info.resolution[1] / 2 + image_info.delta_y + Rfactor * np.sin(angle)
    return X, Y


def derot_1s(az, alt, image_info):

    lat1, lon1 = image_info.latitude, image_info.longitude
    lat2 = image_info.latitude - image_info.latitude_offset
    lon2 = image_info.longitude - image_info.longitude_offset

    lat1_d = lat1 * np.pi / 180.0
    lat2_d = lat2 * np.pi / 180.0
    lon1_d = lon1 * np.pi / 180.0
    lon2_d = lon2 * np.pi / 180.0
    az_d = az * np.pi / 180.
    alt_d = alt * np.pi / 180.

    az2, alt2 = derot_1s_rad(az_d, alt_d, lat1_d, lat2_d, lon2_d-lon1_d)

    az = (az2 * 180. / np.pi) % 360
    alt = alt2 * 180. / np.pi

    return az, alt

def derot_1s_deg(az, alt, lat1, lat2, difflon):

    lat1 *= np.pi / 180.0
    lat2 *= np.pi / 180.0
    difflon *= np.pi / 180.0
    az *= np.pi / 180.
    alt *= np.pi / 180.

    az2, alt2 = derot_1s_rad(az, alt, lat1, lat2, difflon)

    az = (az2 * 180. / np.pi) % 360
    alt = alt2 * 180. / np.pi

    return az, alt


def derot_1s_rad(az, alt, lat1, lat2, difflon):

    angle2 = difflon

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


def derot_1s_rad_zen(lat1, lat2, angle2):

    _sina = np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(angle2)

    alt2 = np.arcsin(_sina)
    _cosa = np.cos(alt2)
    _sinA = -np.sin(angle2) * np.cos(lat1) / _cosa
    _cosA = (np.sin(lat1) - np.sin(lat2) * _sina) / (_cosa * np.cos(lat2))
    az2 = np.arctan2(_sinA, _cosA)

    return az2, alt2

def horiz2xy_(azimuth, altitude, image_info, derotate=True):
    '''
    Return X,Y position in the image from azimuth/altitude horizontal coord.
    azimuth and altitude must be in degrees.
    '''

    if derotate:
        # We have a real azimuth and altitude coordinates. If the camera is not
        # pointing to the zenith, we need to derotate the image.

        azimuth, altitude = derot_1s(azimuth, altitude, image_info)


    Rfactor = image_info.radial_factor * \
        (180.0 / np.pi) * np.sqrt(2 * (1 - np.sin(altitude * np.pi / 180.0)))
    angle = azimuth * np.pi / 180.0 - image_info.azimuth_zeropoint * np.pi / 180.0
    X = image_info.resolution[0] / 2 + image_info.delta_x - Rfactor * np.cos(angle)
    Y = image_info.resolution[1] / 2 + image_info.delta_y + Rfactor * np.sin(angle)
    return X, Y


class ImageInfo3(object):
    latitude = 40.450941
    longitude = -3.726065

    latitude_offset = -0.64102456
    longitude_offset = 0.67447422
    # Scale and rotation
    delta_x = -18.63912476
    delta_y = -31.06643504
    radial_factor = 14.19766968
    azimuth_zeropoint = 88.64589921
    sidereal_time = 0.700320860922
    resolution = [2500, 2500] # There are NAXIS1, NAXIS2 from images

if __name__ == '__main__':

    delta_x = -18.63912476
    delta_y = -31.06643504
    radial_factor = 14.19766968


    from astropy.wcs import WCS
    import matplotlib.pyplot as plt

    image_info = ImageInfo3()
    azimuth_zeropoint = image_info.azimuth_zeropoint

    cdelt = 1.0 / radial_factor

    refx = 2500.0 / 2 + delta_x
    refy = 2500.0 / 2 + delta_y

    w0 = WCS(naxis=2)
    w0.wcs.crpix = [refx, refy]
    w0.wcs.cdelt = [cdelt, cdelt]
    w0.wcs.crval = [0, 90]
    w0.wcs.ctype = ["pLON-ZEA", "pLAT-ZEA"]

    w1 = WCS(naxis=2)
    w1.wcs.crpix = [refx, refy]
    w1.wcs.cdelt = [cdelt, cdelt]
    w1.wcs.crval = [90, 0]
    w1.wcs.ctype = ["pLAT-ZEA", "pLON-ZEA"]
    w1.wcs.lonpole = 180 - azimuth_zeropoint

    ref_az, ref_alt = derot_1s(0, 90, image_info)

    import itertools

    for m in itertools.product([0,90,180,270], [1,0,-1], [1,0,-1], [0,90,180,270], [1,0,-1], [1,0,-1]):
         #ref_az, ref_alt = 321.67057465, 89.180358574
        #print ref_az, ref_alt
        a1,a2,a3,a4,a5,a6 = m

        w2 = WCS(naxis=2)
        w2.wcs.crpix = [refx, refy]
        w2.wcs.cdelt = [cdelt, cdelt]
        w2.wcs.crval = [ref_alt, a1 + a2*ref_az+a3*azimuth_zeropoint]
        w2.wcs.ctype = ["pLAT-ZEA", "pLON-ZEA"]
        w2.wcs.lonpole =  a4 + a5*ref_az+a6*azimuth_zeropoint

        nn = np.zeros((100, 2))
        nn[:,1] = np.arange(100)
        for i in [89]:
        #for i in [49]:
            nn[:,0] = i


            x0, y0 = horiz2xy_(nn[:,1], nn[:,0], image_info, derotate=False)

            x1, y1 = horiz2xy_(nn[:,1], nn[:,0], image_info, derotate=True)

            xy2 = w1.all_world2pix(nn, 1)
            xy3 = w2.all_world2pix(nn, 1)

            plt.plot(x0, y0, color='black')
            plt.plot(xy2[:,0], xy2[:,1], color='green')

            plt.plot(x1, y1, color='red')
            plt.plot(xy3[:,0], xy3[:,1], color='blue')

        print m
        plt.savefig("file_%s.png" % (m,))
        plt.close()
