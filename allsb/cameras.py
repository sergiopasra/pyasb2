
import logging

from astropy.coordinates import EarthLocation
from astropy.time import Time

_logger = logging.getLogger(__name__)


class ExpoInfo:
    def __init__(self):
        self.obstime = None
        self.location = None
        self.exposure_time = None
        self.focal = None
        self.aperture = None


def override_location(parsed_args, expoinfo):
    llat = parsed_args.location_latitude
    llon = parsed_args.location_longitude
    lhei = parsed_args.location_height

    if llat is not None and llon is not None and lhei is not None:
        # lat and lon in deg
        # lhei in m
        location = EarthLocation(
            lat=llat,
            lon=llon,
            height=lhei
        )
        expoinfo.location = location

    return expoinfo


def override_time(parsed_args, expoinfo):
    ltime = parsed_args.location_timestamp

    if ltime is not None:
        obstime = Time(ltime)
        expoinfo.obstime = obstime

    return expoinfo


class CameraInfo:
    def __init__(self, name, saturation):
        self.name = name
        self.saturation = saturation


def read_obsgeo(header):

    _logger.debug('try to read OBSGEO-')
    geokeys = ['OBSGEO-X', 'OBSGEO-Y', 'OBSGEO-Z']

    geovals = []
    for key in geokeys:
        if key in header:
            geovals.append(header[key])
        else:
            raise KeyError('key {} not present'.format(key))
    location = EarthLocation.from_geocentric(geovals[0], geovals[1], geovals[2], unit='m')
    return location



def reader_canon(header):
    _logger.debug('using function %s', 'reader_canon')

    if header['OBSTIME']:
        # Time in UTC
        obstime = Time(header['OBSTIME'])
    else:
        obstime = None

    exptime = float(header['EXPTIME']) # Seconds
    focal = float(header['FOCAL'])  # mm
    aperture = float(header['APERTUR'])  # ratio f/APERTUR

    try:
        location = read_obsgeo(header)
        _logger.debug('location from OBSGEO: %s', location)
    except KeyError:
        location = None

    if location is None:
        _logger.debug('location from header undefined')

    expo = ExpoInfo()
    expo.obstime = obstime
    expo.exposure_time = exptime
    expo.focal = focal
    expo.aperture = aperture
    expo.location = location
    return expo


CAMERAS = {}

camera_info_canon_5d = CameraInfo("Canon EOS 5D Mark II", saturation=65535)
camera_info_canon_6d = CameraInfo("Canon EOS 6D", saturation=65535)

CAMERAS['Canon EOS 5D Mark II'] = (reader_canon, camera_info_canon_5d)
CAMERAS['Canon EOS 6D'] = (reader_canon, camera_info_canon_6d)