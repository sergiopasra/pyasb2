
import argparse
import logging
import os.path

import astropy.io.fits as fits

# FIXME: UHMMM
from allsb.testastro import wcs_calibrate_astrometry_net

import allsb.coords_bann
import allsb.utils as U
import allsb.maxcircle as M
import allsb.catalog


_logger = logging.getLogger(__name__)


def main(args=None):

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

    datafiles = basic_process_files(parsed_args.filename)

    for datafile in datafiles:
        # res = calc_maxcircle(datafile, do_plots=True)
        res = (2855.3553458859865, 1852.4162698338164, 3512.439362575322)
        _logger.debug('centroid, result=%s', res)
        datafile['res'] = res

        if not datafile['wcs']:
            wcs_calibrate_astrometry_net(datafile)
            # update datafile
            datafile['wcs'] = True
            datafile['wcs_type'] = 'catalog'
            filename = datafile['filename']
            wcs_catalog = U.insert_adj_filename(filename, 'corrfile')
            datafile['wcs_catalog'] = wcs_catalog
            datafile['wcs_catalog_type'] = 'AN'
        #
        table_obs = calc_obstable(datafile)

        val = solve_plate_zen_eqa(table_obs, res, 6e-4)
        for par in ['x0', 'y0', 'scale', 'a0', 'E', 'eps']:
            print(par, val.params[par].value)

        print('---')
        val = solve_plate_zen_eqa(table_obs, res, 1e-5)
        for par in ['x0', 'y0', 'scale', 'a0', 'E', 'eps']:
            print(par, val.params[par].value)
        print(val.success)


def solve_plate_zen_eqa(table_obs, res, scale):
    import allsb.fitting
    return allsb.fitting.calc_projection_zen_eqa(table_obs, x0=res[0], y0=res[1], scale=scale)


def solve_plate_zen_pol(table_obs, res):
    import allsb.fitting
    val = allsb.fitting.calc_projection_zen_pol(table_obs, x0=res[0], y0=res[1])
    return val

def calc_obstable(datafile):
    from astropy.coordinates import EarthLocation, AltAz
    import astropy.units as u
    from astropy.time import Time

    # Time and location
    with fits.open(datafile['filename']) as hdulist:
        dateobs = hdulist[0].header['DATE-OBS']
        latitude = hdulist[0].header['LAT']
        longitude = hdulist[0].header['LON']
        height = hdulist[0].header['height']

    # Location Madrid
    latitude = 40.450941
    longitude = -3.726065
    height = 667

    # Physics, UCM
    latitude = 40.4509261
    longitude = -3.7262174

    location = EarthLocation(
        lat=latitude * u.deg,
        lon=longitude * u.deg,
        height=height * u.m
    )
    # print(EarthLocation.of_site('gbt'))
    location_timestamp = Time(dateobs)

    _logger.debug('location %s', location)
    _logger.debug('obs date %s', location_timestamp)

    aaframe = AltAz(
        obstime=location_timestamp,
        location=location,
        temperature=10 * u.deg_C,  # Values to get refraction
        pressure=101325 * u.Pa,
        obswl=0.5 * u.micron
    )

    # read catalog
    filename = datafile['wcs_catalog']
    table_obs = allsb.catalog.read_an_axfile(filename, aaframe)

    #print(table_obs.columns)

    if False:
        import matplotlib.pyplot as plt
        f, ax = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
        allsb.catalog.plt_catalog(ax, table_obs)
        plt.show()
    return table_obs


def calc_maxcircle(datafile, do_plots=False):
    with fits.open(datafile['filename']) as hdulist:
        return M.maxcircle(hdulist, do_plots)


def basic_process_files(filename):

    datafiles = []

    with fits.open(filename) as hdulist:
        # get metadata
        # get data
        data = hdulist[0].data
        thisdatafile = {
            'wcs': False,
            'shape': data.shape, 'filename': filename
        }

        wcs_catalog = U.insert_adj_filename(filename, 'corrfile')
        if os.path.isfile((wcs_catalog)):
            thisdatafile['wcs'] = True
            thisdatafile['wcs_type'] = 'catalog'
            thisdatafile['wcs_catalog'] = wcs_catalog
            thisdatafile['wcs_catalog_type'] = 'AN'
        datafiles.append(thisdatafile)

    return datafiles


def plot_func(data):
    import matplotlib.pyplot as plt

    plt.imshow(data)
    plt.show()


def wcs_calibrate(datafile):
    _logger.debug('shape is {shape}'.format(**datafile))
    plot_func(datafile['data'])


if __name__ == '__main__':

    main()
