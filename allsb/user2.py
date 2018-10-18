
import argparse
import logging
import os.path

import astropy.io.fits as fits
import matplotlib.pyplot as plt

# FIXME: UHMMM
from allsb.testastro import wcs_calibrate_astrometry_net
from allsb.photometry import filter_phot_catalogue, prepare_phot_catalogue
from allsb.projection import proj_zen_eqa_inv

import allsb.coords_bann
import allsb.utils as U
import allsb.maxcircle as M
import allsb.catalog


_logger = logging.getLogger(__name__)


def test_direct3(x0, y0, scale, a0, E, eps):
    import numpy as np
    import math

    x00 = 2855.3553458859865
    x01 = 1852.4162698338164
    rad = 3512.439362575322 / 2.0

    npoints = 200
    ang = np.linspace(0, 3 * np.pi / 2 - 0.5, npoints)
    print('rad is', rad, 1 / rad)

    az = np.linspace(0, 2 * math.pi)

    for d in [0, 20, 40, 60, 80]:
        alt = np.deg2rad(d * np.ones_like(az))
        x, y = proj_zen_eqa_inv(alt, az, x0, y0, scale, a0, E, eps)
        plt.plot(x, y, '-.')

    x1 = x00 + rad * np.cos(ang)
    y1 = x01 + rad * np.sin(ang)

    plt.plot(x1, y1, 'r-')

    plt.xlim([0, 5633])
    plt.ylim([0, 3753])
    plt.show()



def solve_plate_zen_eqa(table_obs, res):
    import allsb.fitting
    import lmfit
    import math
    scale = 0.00083 # 2.0 / res[2]
    x0 = res[0]
    y0 = res[1]
    max0 = 1.5 * scale
    #print('range for scale', 0.5 * scale , scale, max0)
    print('range for scale', 0.0007, 0.00083, 0.0009)
    params = lmfit.Parameters()
    params.add('scale', value=scale, min=0.5 * scale, max=max0)
    params.add('x0', value=x0, min=x0-100, max=x0+100)
    params.add('y0', value=y0, min=y0-100, max=y0+100)
    params.add('a0', value=0, vary=True, min=-math.pi, max=math.pi)
    params.add('E', value=0, vary=True, min=-math.pi, max=math.pi)
    params.add('eps', value=0, vary=True, min=-math.pi / 2.0, max=math.pi / 2.0)

    # FIXME: hardcoded
    x = table_obs['field_x'] + 2355
    y = table_obs['field_y'] + 1352
    alt = table_obs['alt']
    az = table_obs['az']

    # print('fit projection')
    res = allsb.fitting.calc_projection_zen_eqa(params, x, y, alt, az)
    return res


def solve_plate_zen_pol(table_obs, res):
    import allsb.fitting
    val = allsb.fitting.calc_projection_zen_pol(table_obs, x0=res[0], y0=res[1])
    return val


def calc_aaframe(datafile):
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
    return aaframe


def calc_obstable(wcs_catalog_filename, aaframe, do_plot=False):

    # read catalog
    table_obs = allsb.catalog.read_an_axfile(wcs_catalog_filename, aaframe)

    if do_plot:
        f, ax = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
        allsb.catalog.plt_catalog(ax, table_obs)
        plt.show()
    return table_obs


def calc_maxcircle(datafile, do_plots=False):
    with fits.open(datafile['filename']) as hdulist:
        return M.maxcircle(hdulist, do_plots)


def basic_process_files(filename):

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

        return thisdatafile


def plot_func(data):
    import matplotlib.pyplot as plt

    plt.imshow(data)
    plt.show()


def wcs_calibrate(datafile):
    _logger.debug('shape is {shape}'.format(**datafile))
    plot_func(datafile['data'])


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

    datafile = basic_process_files(parsed_args.filename)

    #res = calc_maxcircle(datafile, do_plots=True)
    res = (2855.3553458859865, 1852.4162698338164, 3512.439362575322)
    _logger.debug('centroid, result=%s', res)
    datafile['res'] = res
    if not datafile['wcs']:
        _logger.debug('calling wcs_calibrate_astrometry_net')
        wcs_calibrate_astrometry_net(datafile)
        # update datafile
        datafile['wcs'] = True
        datafile['wcs_type'] = 'catalog'
        filename = datafile['filename']
        wcs_catalog = U.insert_adj_filename(filename, 'corrfile')
        datafile['wcs_catalog'] = wcs_catalog
        datafile['wcs_catalog_type'] = 'AN'
    #
    # AltAz reference frame
    _logger.debug('compute AltAz reference frame')
    aaframe = calc_aaframe(datafile)

    _logger.debug('calling calc_obstable')

    table_obs = calc_obstable(datafile['wcs_catalog'], aaframe)
    _logger.debug('calling solve_plate_zen_eqa')
    val = solve_plate_zen_eqa(table_obs, res)

    for par in ['x0', 'y0', 'scale', 'a0', 'E', 'eps']:
        print(par, val.params[par].value)

    # Star catalogue and more
    table = filter_phot_catalogue(catfile)
    table = prepare_phot_catalogue(table, aaframe, min_altitude=25)
    _logger.debug('we have %s photo stars', len(table))

    # Compute expected X/Y coordinates of photo stars
    import numpy as np

    itl = np.array([table['alt'], table['az']]).T

    print(itl.shape)
    x0 = val.params['x0'].value
    y0 = val.params['y0'].value
    scale = val.params['scale'].value
    a0 = val.params['a0'].value
    E = val.params['E'].value
    eps = val.params['eps'].value

    # test_direct3(x0, y0, scale, a0, E, eps)

    xm, ym = proj_zen_eqa_inv(table_obs['alt'], table_obs['az'], x0, y0, scale, a0, E, eps)
    plt.scatter(xm, ym)
    xm, ym = proj_zen_eqa_inv(table['alt'], table['az'], x0, y0, scale, a0, E, eps)
    plt.scatter(xm, ym)

    for a, b in zip(xm, ym):
        print(a, b)
    #plt.scatter(table_obs['field_x'] + 2355, table_obs['field_y'] + 1352)
    plt.show()



if __name__ == '__main__':

    main()
