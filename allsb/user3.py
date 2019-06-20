
import argparse
import logging
import os.path
import math

import astropy.io.fits as fits
import matplotlib.pyplot as plt

# FIXME: UHMMM
from allsb.testastro import wcs_calibrate_astrometry_net
from allsb.photometry import filter_phot_catalogue, filter_catalogue, prepare_phot_catalogue, prepare_astrometry_catalogue
from allsb.projection import proj_zen_eqa_inv, create_wcs, create_wcs_radec

import allsb.coords_bann
import allsb.utils as U
import allsb.maxcircle as M
import allsb.catalog


_logger = logging.getLogger(__name__)


def test_direct3(x0, y0, scale, a0, E, eps):
    import numpy as np

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
    scale = 2.0 / res[2]
    x0 = res[0]
    y0 = res[1]
    max0 = 1.5 * scale
    print('range for scale', 0.5 * scale , scale, max0)
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
    # print(datafile)

    #
    # AltAz reference frame
    _logger.debug('compute AltAz reference frame')
    aaframe = calc_aaframe(datafile)

    min_magnitude = 6.0
    min_altitude = 25

    npoints = 40

    table = filter_catalogue(catfile, min_magnitude)
    table = prepare_astrometry_catalogue(table, aaframe, min_altitude)
    _logger.debug('we have %s stars for astrometry', len(table))

    table_n = table[:npoints]

    x0 = 2825.1933894004724
    y0 = 1894.1714870521146
    scale = 0.0007703000466395117
    #scale = 0.00075
    a0 = 1.4315396267055513
    E = -0.0014445185050160347
    eps = -0.0700522925367364

    #x0 = 2855.3553458859865
    #y0 = 1852.4162698338164
    #E = 0
    #eps = 0
    #a0 = math.pi / 2
    xm, ym = proj_zen_eqa_inv(table_n['alt'], table_n['az'], x0, y0, scale, a0, E, eps)
    import numpy as np
    import sep

    #_logger.debug('calling calc_obstable')
    data = fits.getdata(datafile['filename'])
    data = data.astype('float32')
    y,x = np.ogrid[0:data.shape[0], 0:data.shape[1]]
    mask = (y - 1852.4162698338164)**2 + (x - 2855.3553458859865)**2 >= (2800 / 2.0)**2

    m, s = np.mean(data), np.std(data)
    #plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-s, vmax=m+s, origin='lower')
    # plt.imshow(mask)
    #plt.colorbar()
    #plt.show()

    bkg = sep.Background(data, mask=mask)
    data_sub = data - bkg
    objects = sep.extract(data_sub, thresh=3.0, err=bkg.globalrms, mask=mask, minarea=10.0)

    # how many objects were detected
    _logger.debug('detected %d objects', len(objects))
    idx_by_flux = np.flip(np.argsort(objects['flux']))
    objects_by_flux = objects[idx_by_flux]

    from matplotlib.patches import Ellipse

    # plot background-subtracted image
    plot_p = False
    N = npoints
    objects_coords = []
    for obj in objects_by_flux[:N]:
        objects_coords.append((obj['x'], obj['y']))

    if plot_p:
        fig, ax = plt.subplots()
        m, s = np.mean(data_sub), np.std(data_sub)
        im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower', alpha=0.3)

        # plot an ellipse for each object
        for obj in objects_by_flux[:N]:
            e = Ellipse(xy=(obj['x'], obj['y']),
                        width=6*obj['a'],
                        height=6*obj['b'],
                        angle=obj['theta'] * 180. / np.pi)
            e.set_facecolor('none')
            e.set_edgecolor('red')
            ax.add_artist(e)

        ax.scatter(xm.data, ym.data, marker='.')
        plt.show()

    from scipy.spatial import KDTree

    computed_objs = np.column_stack((xm.data, ym.data))
    kdtree = KDTree(computed_objs)

    observed_objs = np.array(objects_coords)
    dis, index = kdtree.query(observed_objs, distance_upper_bound=50.0)

    object_p = False
    x2 = []
    y2 = []
    alt2 = []
    az2 = []
    for idx, (iidx, d) in enumerate(zip(index, dis)):
        if d < np.inf:
            obj = objects_by_flux[idx]
            # print(obj['x'], obj['y'], iidx, d, computed_objs[iidx])
            x2.append(obj['x'])
            y2.append(obj['y'])
            alt2.append(table_n['alt'][iidx])
            az2.append(table_n['az'][iidx])
            # plot an ellipse for each object
            if object_p:
                fig, ax = plt.subplots()
                im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
                               vmin=m - s, vmax=m + s, origin='lower', alpha=0.3)

                e = Ellipse(xy=(obj['x'], obj['y']),
                            width=6 * obj['a'],
                            height=6 * obj['b'],
                            angle=obj['theta'] * 180. / np.pi)
                e.set_facecolor('none')
                e.set_edgecolor('red')
                ax.add_artist(e)
                ax.scatter(xm.data[iidx], ym.data[iidx], marker='.')
                plt.xlim(int(obj['x'])- 50, int(obj['x']) + 50)
                plt.ylim(int(obj['y']) - 50, int(obj['y']) + 50)
                plt.show()

    _logger.debug('calling solve_plate_zen_eqa (iter2)')

    import allsb.fitting
    import lmfit
    import math

    params = lmfit.Parameters()
    params.add('scale', value=scale, min=0.0007, max=0.0009)
    params.add('x0', value=x0, min=x0-100, max=x0+100)
    params.add('y0', value=y0, min=y0-100, max=y0+100)
    params.add('a0', value=a0, vary=True, min=-math.pi, max=math.pi)
    params.add('E', value=E, vary=True, min=-math.pi, max=math.pi)
    # params.add('eps', value=eps, vary=True, min=-math.pi / 2.0, max=math.pi / 2.0)
    params.add('eps', value=abs(eps), vary=True, min=0, max=math.pi)

    _logger.debug('fit projection iter2')
    val2 = allsb.fitting.calc_projection_zen_eqa(params, x2, y2, alt2, az2)

    # x0 2828.0618227093737
    # y0 1875.4853799821208
    # scale 0.0007862965093399066
    # a0 1.4352118701541068
    # E -3.0405862093935454
    # eps 0.06124805298368204

    # x0 2828.0776745053095
    # y0 1875.2830144745894
    # scale 0.0007862979768017623
    # a0 1.4352134570547
    # E 0.1018521605640128
    # eps -0.0610795008432643

    for par in ['x0', 'y0', 'scale', 'a0', 'E', 'eps']:
        _logger.debug("%s %f", par, val2.params[par].value)

    x0_2 = val2.params['x0'].value
    y0_2 = val2.params['y0'].value
    scale_2 = val2.params['scale'].value
    a0_2 = val2.params['a0'].value
    E_2 = val2.params['E'].value
    eps_2 = val2.params['eps'].value

    # create wcs
    wcs = create_wcs(x0_2, y0_2, scale_2, a0_2, E_2, eps_2)
    header = wcs.to_header()
    header.tofile('headerM.txt', overwrite=True)

    _logger.debug('write astrometry in WCS format')
    with fits.open(datafile['filename']) as hdul:
        hdul[0].header.extend(header)
        hdul.writeto('alt.fits', overwrite=True)

    plot_p2 = False
    if plot_p2:
        table2 = filter_catalogue(catfile, min_magnitude=3.5)
        min_altitude = 20
        table2 = prepare_phot_catalogue(table2, aaframe, min_altitude)
        _logger.debug('we have %s photo stars', len(table))
        table_n2 = table2[:]
        xm2, ym2 = proj_zen_eqa_inv(table_n2['alt'], table_n2['az'], x0_2, y0_2, scale_2, a0_2, E_2, eps_2)

        # pole = SkyCoord(0, 90, unit='deg')
        # pole_coords = (, )
        xz_2, yz_2 = proj_zen_eqa_inv(np.array([math.pi / 2, np.deg2rad(40.43915023)]), np.array([0.0, 0]), x0_2, y0_2, scale_2, a0_2, E_2, eps_2)

        fig = plt.figure()
        ax = fig.add_subplot(1 ,1 ,1, projection=wcs)
        m, s = np.mean(data_sub), np.std(data_sub)
        ax.imshow(data_sub, interpolation='nearest', cmap='gray',
                  vmin=m-s, vmax=m+s, origin='lower', alpha=1)
        ax.grid(color='green', ls='solid')
        # plot an ellipse for each object
        for obj in objects_by_flux:
            e = Ellipse(xy=(obj['x'], obj['y']),
                        width=6*obj['a'],
                        height=6*obj['b'],
                        angle=obj['theta'] * 180. / np.pi)
            e.set_facecolor('none')
            e.set_edgecolor('red')
            ax.add_artist(e)

        ax.scatter(xm2.data, ym2.data, marker='.')
        ax.scatter([x0_2], [y0_2], marker='x')
        ax.scatter([xz_2], [yz_2], marker='o')
        plt.show()

    table3 = filter_phot_catalogue(catfile, min_magnitude=3.5)
    table3 = prepare_phot_catalogue(table3, aaframe, min_altitude)

    plot_p3 = True
    table_n3 = table3[:]

    xm3, ym3 = proj_zen_eqa_inv(table_n3['alt'], table_n3['az'], x0_2, y0_2, scale_2, a0_2, E_2, eps_2)

    if plot_p3:
        _logger.debug('we have %s photo stars', len(table3))
        fig = plt.figure()
        # ax = fig.add_subplot(1, 1, 1, projection=wcs)
        ax = fig.add_subplot(1, 1, 1)
        m, s = np.mean(data_sub), np.std(data_sub)
        ax.imshow(data_sub, interpolation='nearest', cmap='gray',
                  vmin=m-s, vmax=m+s, origin='lower', alpha=1)
        # ax.grid(color='green', ls='solid')
        # plot an ellipse for each object
        for obj in objects_by_flux:
            e = Ellipse(xy=(obj['x'], obj['y']),
                        width=6*obj['a'],
                        height=6*obj['b'],
                        angle=obj['theta'] * 180. / np.pi)
            e.set_facecolor('none')
            e.set_edgecolor('red')
            ax.add_artist(e)

        ax.scatter(xm3.data, ym3.data, color='green', marker='.')
        plt.show()


    from photutils import CircularAperture, aperture_photometry, CircularAnnulus

    positions = np.array([xm3.data, ym3.data]).T
    apertures = CircularAperture(positions, r=8.0)
    annulus_apertures = CircularAnnulus(positions, r_in=8.0, r_out=10.0)
    apers = [apertures, annulus_apertures]
    phot_table = aperture_photometry(data_sub, apers)
    print(phot_table.colnames)
    bkg_mean = phot_table['aperture_sum_1'] / annulus_apertures.area()
    bkg_sum = bkg_mean * apertures.area()
    final_sum = phot_table['aperture_sum_0'] - bkg_sum

    # find excess
    m1a = bkg_sum > 25000

    # phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
    # print(phot_table)

    # yb2 = table_n3['vmag']  + 2.5 * np.log10(phot_table['aperture_sum_0'])
    exptime = 16.0
    yb = table_n3['vmag'] + 2.5 * np.log10(final_sum / exptime)
    xb = 1.0 / np.cos(math.pi / 2 - table_n3['alt'])

    # Use only airmass upto 1.8
    airmass_mask = xb < 2.5
    t_mask = airmass_mask & ~m1a

    res = np.polyfit(xb[t_mask], yb[t_mask], deg=1)
    print(res)
    xl = np.linspace(1.0, 2.5)
    yl = np.polyval(res, xl)
    plt.xlabel('airmass')
    plt.ylabel('m0 + 2.5 log F')
    plt.scatter(xb[t_mask], yb[t_mask], marker='.', color='red')
    # plt.scatter(xb[t_mask], yb2[t_mask], marker='x', color='blue')
    plt.text(1.6 - 0.2, 15.5, "m0 + 2.5 log F = {0[1]:5.2f} {0[0]:5.2f} sec z ".format(res))
    plt.plot(xl, yl)
    plt.show()

#    sigma_clip = SigmaClip(sigma=3., iters=10)
#    bkg_estimator = MedianBackground()
#    bkg = Background2D(data_sub, (50, 50), filter_size=(3, 3),
#                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
#                       mask=mask)

#    print(bkg.background_median)
#    print(bkg.background_rms_median)

#    plt.imshow(bkg.background, origin='lower', cmap='Greys_r')
#    plt.show()


def miniplot(data):
    import numpy as np
    m, s = np.mean(data), np.std(data)
    plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-s, vmax=m+s, origin='lower')
    plt.colorbar()
    plt.show()


if __name__ == '__main__':

    main()