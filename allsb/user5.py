
import argparse
import logging
import os.path

import matplotlib.pyplot as plt
import astropy.io.fits as fits

import allsb.utils
import allsb.cameras as cams
from allsb.testastro import wcs_calibrate_astrometry_net

_logger = logging.getLogger(__name__)


def basic_process_files(cameras, filename):

    with fits.open(filename) as hdulist:
        # get metadata
        # get data
        data = hdulist[0].data
        thisdatafile = dict(wcs=False, shape=data.shape, filename=filename)

        # Get CAMERA?
        camera_model = hdulist[0].header['CAMERA']

        header_reader, camera_info = cameras[camera_model]
        expoinfo = header_reader(hdulist[0].header)
        thisdatafile["expoinfo"] = expoinfo
        thisdatafile["camerainfo"] = camera_info
        thisdatafile["saturation_mask"] = data >= camera_info.saturation

        wcs_catalog = allsb.utils.insert_adj_filename(filename, 'corrfile')
        if os.path.isfile(wcs_catalog):
            thisdatafile['wcs'] = True
            thisdatafile['wcs_type'] = 'catalog'
            thisdatafile['wcs_catalog'] = wcs_catalog
            thisdatafile['wcs_catalog_type'] = 'AN'

        return thisdatafile


def calc_aaframe(location, obstime):
    from astropy.coordinates import AltAz
    import astropy.units as u

    aaframe = AltAz(
        obstime=obstime,
        location=location,
        temperature=10 * u.deg_C,  # Values to get refraction
        pressure=101325 * u.Pa,
        obswl=0.5 * u.micron # This value can come from the filter...
    )
    return aaframe


def compute_initial_wcs(image_info, center_ij=None):
    import math

    rad0 = 1100

    res = (1282.4836, 1911.9343, math.sqrt(2) / rad0)

    if center_ij is None:
        shape = image_info['shape']
        center_ij = (int(shape[0] // 2), int(shape[1] // 2))

    _logger.debug('centroid, result=%s', res)
    # image_info['res'] = res
    if not image_info['wcs']:
        _logger.debug('calling wcs_calibrate_astrometry_net')
        wcs_calibrate_astrometry_net(image_info['filename'], center_ij)
        # update image_info
        image_info['wcs'] = True
        image_info['wcs_type'] = 'catalog'
        filename = image_info['filename']
        wcs_catalog = allsb.utils.insert_adj_filename(filename, 'corrfile')
        image_info['wcs_catalog'] = wcs_catalog
        image_info['wcs_catalog_type'] = 'AN'


def main(args=None):

    catfile = '/home/spr/devel/github/pyasb2/catalog2.txt'

    parser = argparse.ArgumentParser()
    loc_group = parser.add_argument_group('location')
    loc_group.add_argument('--location-latitude', type=float)
    loc_group.add_argument('--location-longitude', type=float)
    loc_group.add_argument('--location-height', type=float)
    loc_group.add_argument('--location-timestamp')
    parser.add_argument('filename')

    parsed_args = parser.parse_args(args)

    logging.basicConfig(level=logging.DEBUG)

    _logger.debug('filename is: %s', parsed_args.filename)

    # Basic characterize
    image_info = basic_process_files(cams.CAMERAS, parsed_args.filename)

    # Override location
    cams.override_location(parsed_args, image_info['expoinfo'])
    cams.override_time(parsed_args, image_info['expoinfo'])

    _logger.debug('location %s', image_info['expoinfo'].location)
    _logger.debug('obstime %s', image_info['expoinfo'].obstime)

    altaz_frame = calc_aaframe(image_info['expoinfo'].location, image_info['expoinfo'].obstime)
    _logger.debug('altaz frame %s', altaz_frame)

    compute_initial_wcs(image_info)

    # print(image_info)

    # plt.matshow(result['saturation_mask'])
    # plt.show()


if __name__ == '__main__':

    main()
