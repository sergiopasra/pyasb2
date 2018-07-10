

"""Command line interface"""

import argparse
import configparser
import logging
import math

import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
from astropy.wcs import WCS

# import datetime

from allsb.reduction import reduction


def nwarp(angles):
    # return angles
    ang = np.fmod(angles, 2 * math.pi)
    neg = ang < 0
    ang[neg] += 2 * math.pi
    return ang


def main():

    config = configparser.ConfigParser()

    parser = argparse.ArgumentParser(description='Process ASTMON images.')
    parser.add_argument('-c', '--config',
                        help='Configuration file',
                        default='config.ini')
    parser.add_argument('images', nargs='+', help='ASTMON image')

    args = parser.parse_args()
    config.read(args.config)

    for filename in args.images:
        print('filename:', filename)

        with fits.open(filename) as hdul:
            header = hdul[0].header
            # date_str = header['DATE']
            # image_dt = datetime.datetime.strptime(date_str, "%Y%m%d_%H%M%S")
            image_filter = header['FILTER']
            # calibrations

        calib_d = config['calibrations_{}'.format(image_filter)]
        fc = reduction(filename, calib_d['darkframe'], calib_d['flatfield'])
        fc.write('test.fits', overwrite=True)



    latitude = 40.450941
    refx, refy = 1231.36087524, 1218.93356496
    cdeltx, cdelty = 0.0704340939421, 0.0704340939421
    ref_alt, ref_az = 89.180358574, 141.230103929
    az_zero0 = 90 - 88.64589921

    w2 = WCS(naxis=2)
    w2.wcs.crpix = [refx, refy]
    w2.wcs.cdelt = [cdeltx, cdelty]
    w2.wcs.crval = [ref_alt, ref_az]
    w2.wcs.ctype = ["pLAT-ZEA", "pLON-ZEA"]
    w2.wcs.lonpole = ref_az + 90 + az_zero0

    load(w2)

    # catalogs

def apply_refraction(real_alt):
    # https://en.wikipedia.org/wiki/Atmospheric_refraction
    # This model is valid at 1 Atm and 10 C
    ang = (real_alt + 10.3 / (real_alt + 5.11))
    cot_ang = (90 - ang) / 180.0 * np.pi
    rr = 1.02 * np.tan(cot_ang) # In arc minutes
    app_alt = real_alt + rr / 60.0
    return app_alt


def correct_refraction(app_alt):
    ang = app_alt + 7.31 / (app_alt + 4.4)
    cot_ang = (90 - ang) / 180.0 * np.pi
    rr = np.tan(cot_ang) # Arc minutes
    real_alt = app_alt - rr / 60.0
    return real_alt


def create_wcs1(customwcs):
    cdelt = customwcs.cdelt

    refx = customwcs.refx
    refy = customwcs.refy

    w1 = WCS(naxis=2)
    w1.wcs.crpix = [refx, refy]
    w1.wcs.cdelt = [cdelt, cdelt]
    w1.wcs.crval = [90, 0]
    w1.wcs.ctype = ["pLAT-ZEA", "pLON-ZEA"]
    w1.wcs.lonpole = 180 - customwcs.azimuth_zeropoint

    return w1

def create_wcs2(customwcs):
    ref_az, ref_alt = customwcs.derot_1s(0, 90)
    cdelt = customwcs.cdelt

    refx = customwcs.refx
    refy = customwcs.refy

    w2 = WCS(naxis=2)
    w2.wcs.crpix = [refx, refy]
    w2.wcs.cdelt = [cdelt, cdelt]
    w2.wcs.crval = [ref_alt, ref_az]
    w2.wcs.ctype = ["pLAT-ZEA", "pLON-ZEA"]
    w2.wcs.lonpole = 180-customwcs.azimuth_zeropoint + ref_az
    return w2


def filter_phot(table):

    table.add_column(Column(np.zeros_like(table['vmag'], dtype='bool'), name='photo'))

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
    table['photo'] = photo_mask

#    good_phot = table.group_by('photo')
#
#    ntable = good_phot.groups[1]
#
#    extreme_colors_mask = (ntable['bv'] < -1) & (ntable['bv'] > 2)
#    bright_stars_mask = ntable['vmag'] <= max_magnitude

    # Check this can be done
#    ntable['photo'][extreme_colors_mask] = False
#    ntable['photo'][~bright_stars_mask] = True

#    table_photometric = table
    return table


def load(wcs):

    from astropy import units as u
    from astropy.coordinates import SkyCoord, EarthLocation, AltAz
    from astropy.coordinates import FK5
    from astropy.time import Time
    from astropy.table import Column
    from astropy.table import Table
    import numpy
    from astropy.io import fits

    magfield = 'vmag'
    colorfields = ['uv', 'bv', 'rv', 'iv']

    min_altitude = 15
    max_magnitude = 5 # In V band, for the moment

    # Location
    latitude = 40.450941
    longitude = -3.726065
    height = 667

    location = EarthLocation(
        lat=latitude * u.deg,
        lon=longitude * u.deg,
        height=height * u.m
    )

    time = Time("2013-09-12 01:17:09")

    #customwcs = CustomWCS()

    catfile = 'catalog2.txt'
    table = Table.read(catfile, format='ascii.csv')

    if True:

        table = table[:300]

        # full            | visible
        # |- Photometric  | visible

        mm = SkyCoord(ra=table['raj1950'], dec=table['dej1950'], unit=(u.hourangle, u.deg), frame=FK5(equinox='J1950'))

        mm_altz = mm.transform_to(AltAz(obstime=time,location=location))

        table.add_column(Column(mm_altz.alt.degree, name='alt'))
        table.add_column(Column(mm_altz.az.degree, name='az'))
        table.add_column(Column(mm.dec.degree, name='dec'))
        table.add_column(Column(mm.ra.degree, name='ra'))

        visibility_mask = mm_altz.alt.degree > min_altitude

        visible_catalog = table[visibility_mask]

        # And photometry

        # Handle refraction
        # Astropy does ref correction using ERFA complex models.
        # Perhaps we should use that
        # It requires Pressure, Temp, Humidity, etc
        app_altz = apply_refraction(visible_catalog['alt'])

        visible_catalog.add_column(Column(app_altz, name='alt_app'))

        itl = np.array([visible_catalog['alt_app'], visible_catalog['az']]).T

        res = wcs.all_world2pix(itl, 1)

        #vx, vy = customwcs.horiz2xy(visible_catalog['az'], visible_catalog['alt_app'], derotate=True)

        visible_catalog.add_column(Column(res[:,0], name='x'))
        visible_catalog.add_column(Column(res[:,1], name='y'))

        by_phot = visible_catalog.group_by('photo')
        visible_catalog_phot = by_phot.groups[1]

        print('filtering catalog done')
        print('photometric', len(visible_catalog_phot))
        print('full', len(visible_catalog))

        fig = plt.figure()
        ax1 = fig.add_axes([0.0, 0.0, 1.0, 1.0], projection=wcs)

        ax1.coords.grid(color='blue', alpha=1, linestyle='solid')
        overlay1 = ax1.get_coords_overlay(wcs)
        overlay1.grid(color='white', linestyle='solid', alpha=1)

        data_sub = fits.getdata('test.fits')
        im = ax1.imshow(data_sub, interpolation='nearest', cmap='gray',
                        origin='lower')

        ax1.scatter(visible_catalog_phot['x'], visible_catalog_phot['y'])
        plt.show()

        visible_catalog.write('cat1.csv', format='ascii.csv', overwrite=True)
        numpy.savetxt('stars.txt', res)

        with open('dum.reg', 'w') as fd:
            fd.write("# Region file format: DS9 version 4.1\n")
            fd.write('global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
            fd.write('physical\n')
            for m in visible_catalog_phot:
                print(m['name'])
                fd.write('circle({},{},8.0) # text={{ {} }}\n'.format(m['x'], m['y'], m['name']))


if __name__ == '__main__':

    main()