

"""Command line interface"""

import argparse
import configparser
import logging

import matplotlib.pyplot as plt
import numpy as np
import math

import pyregion
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import FK5
from astropy.time import Time
from astropy.table import Column
from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import fits

# import datetime
_logger = logging.getLogger(__name__)


def nwarp(angles):
    # return angles
    ang = np.fmod(angles, 2 * math.pi)
    neg = ang < 0
    ang[neg] += 2 * math.pi
    return ang


def filter_phot_catalogue(catfile, min_magnitude=6.0):
    _logger.debug('read photo catalog')
    table = Table.read(catfile, format='ascii.csv')
    _logger.debug('filter stars by mag/photomety min_mag=%s', min_magnitude)

    by_phot = table.group_by('photo')

    catalog_phot = by_phot.groups[1]
    # visible_catalog_phot = visible_catalog
    current_mag = catalog_phot['vmag']
    magnitude_mask = current_mag < min_magnitude
    cpl = catalog_phot[magnitude_mask]
    return cpl


def filter_catalogue(catfile, min_magnitude=6.0):
    _logger.debug('read catalog')
    table = Table.read(catfile, format='ascii.csv')
    _logger.debug('filter stars by mag min_mag=%s', min_magnitude)

    current_mag = table['vmag']
    magnitude_mask = current_mag < min_magnitude
    cpl = table[magnitude_mask]
    return cpl


def prepare_phot_catalogue(table, aaframe, min_altitude=25):

    # load astrometry

    # Load star catalog
    _logger.debug('convert coordinates from J1950')
    #print(table)
    coords_sky = SkyCoord(
        ra=table['raj1950'],
        dec=table['dej1950'],
        unit=(u.hourangle, u.deg),
        frame=FK5(equinox='J1950')
    )

    # star postition predictions
    _logger.debug('add RADec columns')
    table.add_column(Column(coords_sky.dec.degree, name='dec'))
    table.add_column(Column(coords_sky.ra.degree, name='ra'))

    coords_altz = coords_sky.transform_to(aaframe)

    _logger.debug('add AltAz columns')
    table.add_column(Column(coords_altz.alt.radian, name='alt'))
    table.add_column(Column(coords_altz.az.radian, name='az'))
    table.add_column(Column(coords_altz.alt.degree, name='alt_deg'))
    table.add_column(Column(coords_altz.az.degree, name='az_deg'))

    _logger.debug('filter columns for visibility, alt>%s', min_altitude)
    # Use color if needed
    visibility_mask = coords_altz.alt.degree > min_altitude
    visible_catalog = table[visibility_mask]
    return visible_catalog


def prepare_astrometry_catalogue(table, aaframe, min_altitude=25):

    # load astrometry

    # Load star catalog
    _logger.debug('convert coordinates from J1950')
    #print(table)
    coords_sky = SkyCoord(
        ra=table['raj1950'],
        dec=table['dej1950'],
        unit=(u.hourangle, u.deg),
        frame=FK5(equinox='J1950')
    )

    # star postition predictions
    _logger.debug('add RADec columns')
    table.add_column(Column(coords_sky.dec.degree, name='dec'))
    table.add_column(Column(coords_sky.ra.degree, name='ra'))

    coords_altz = coords_sky.transform_to(aaframe)

    _logger.debug('add AltAz columns')
    table.add_column(Column(coords_altz.alt.radian, name='alt'))
    table.add_column(Column(coords_altz.az.radian, name='az'))
    table.add_column(Column(coords_altz.alt.degree, name='alt_deg'))
    table.add_column(Column(coords_altz.az.degree, name='az_deg'))

    _logger.debug('filter columns for visibility, alt>%s', min_altitude)
    # Use color if needed
    visibility_mask = coords_altz.alt.degree > min_altitude
    visible_catalog = table[visibility_mask]
    return visible_catalog



def main():

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
    aaframe = AltAz(obstime=time, location=location,
                    temperature=10 * u.deg_C,
                    pressure=101325 * u.Pa,
                    obswl=0.5 * u.micron
                    )

    min_magnitude = 6.0
    min_altitude = 25
    catfile = '/home/spr/devel/github/pyasb2/catalog2.txt'

    table = filter_phot_catalogue(catfile, min_magnitude)
    table = prepare_phot_catalogue(table, aaframe, min_altitude)
    _logger.debug('we have %s photo stars', len(table))

    # load astrometry
    #header = fits.Header.fromfile('header2.txt')
    #wcs = WCS(header)

    print('catalog', len(table))
    #print(vcpl)
    #print(vcpl['alt_deg'])
    #itl = np.array([table['alt_deg'], table['az_deg']]).T

    #print(itl.shape)
    #print(wcs)

    #res = wcs.all_world2pix(itl, 1)

    #table.add_column(Column(res[:, 0], name='x'))
    #table.add_column(Column(res[:, 1], name='y'))

    nom_theta = table['alt']
    nom_phi = table['az']

    if False:
        f, axs = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
        axs.set_theta_offset(-math.pi / 2)
        axs.set_theta_direction(-1)
        axs.scatter(nom_phi, 90 - nom_theta / math.pi * 180)
        axs.set_rmax(90)
        plt.show()

    with open('photometry.reg', 'w') as fd:
        fd.write("# Region file format: DS9 version 4.1\n")
        fd.write('global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        fd.write('physical\n')
        for m in vcpl:
            # print(m['name'])
            fd.write('circle({},{},8.0) # text={{ {} }}\n'.format(m['x'], m['y'], m['name']))

    data, header = fits.getdata('test.fits', header=True)

    p_table = photometry_aper(data, vcpl, wcs)
    print(p_table)

    mag_col = vcpl['vmag']
    aper_sum = p_table['aperture_sum']
    mask1 = aper_sum > 0
    print('Mask sum', sum(mask1))
    term3 = 2.5 * np.log10(aper_sum[mask1])
    term2 = term3 + mag_col[mask1]
    term1 = 1 / np.sin(vcpl['alt'][mask1])

    fitt = np.polyfit(term1, term2, 1)
    print(fitt)

    ux = np.linspace(1, 2.5)
    uy = np.polyval(fitt, ux)
    #plt.xlim([1, 2])
    #plt.ylim([0, 20000])
    plt.scatter(term1, term2)
    plt.plot(ux, uy, 'r')
    plt.show()


def create_wcs(pr):
    r0 = pr[0]
    r1 = pr[1]
    ref_delta_p = math.pi / 2 - pr[2]
    ref_alpha_p = pr[3]
    ref_phi_p = pr[4]
    delta_p = np.rad2deg(ref_delta_p)
    alpha_p = np.rad2deg(ref_alpha_p)
    phi_p = np.rad2deg(ref_phi_p)
    #m = total(nom_theta, nom_phi, coords_x, coords_y, r0, r1, math.pi / 2, 0, ref_phi_p)
    #print('after', m)

    cdeltx, cdelty = 0.0704340939421, 0.0704340939421

    w1 = WCS(naxis=2)
    w1.wcs.crpix = [r0, r1]
    w1.wcs.cdelt = [cdeltx, cdelty]
    w1.wcs.crval = [delta_p, alpha_p]
    w1.wcs.ctype = ["pLAT-ZEA", "pLON-ZEA"]
    w1.wcs.lonpole = phi_p - 90
    return w1


def photometry_aper(data, vcpl, wcs):
    from astropy.stats import sigma_clipped_stats
    from astropy.stats import SigmaClip
    from photutils import CircularAperture, aperture_photometry
    from photutils import CircularAnnulus
    from photutils import Background2D, MedianBackground

    if False:
        fig = plt.figure()
        ax1 = fig.add_axes([0.0, 0.0, 1.0, 1.0], projection=wcs)

        ax1.coords.grid(color='blue', alpha=1, linestyle='solid')
        overlay1 = ax1.get_coords_overlay(wcs)
        overlay1.grid(color='white', linestyle='solid', alpha=1)

        m, s = np.mean(data), np.std(data)
        mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5)
        print(mean, median, std)
        im = ax1.imshow(data, interpolation='nearest', cmap='gray',
                        vmin=m - s, vmax=m + s,
                        origin='lower')

        ax1.scatter(vcpl['x'], vcpl['y'])
        plt.show()

    print('background estimation')
    sigma_clip = SigmaClip(sigma=3., iters=10)
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (50, 50), filter_size=(3, 3),
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    print('done')
    data2 = data - bkg.background
    if False:
        fig = plt.figure()
        ax1 = fig.add_axes([0.0, 0.0, 1.0, 1.0], projection=wcs)

        ax1.coords.grid(color='blue', alpha=1, linestyle='solid')
        overlay1 = ax1.get_coords_overlay(wcs)
        overlay1.grid(color='white', linestyle='solid', alpha=1)

        m, s = np.mean(data2), np.std(data2)
        mean, median, std = sigma_clipped_stats(data2, sigma=3.0, iters=5)
        im = ax1.imshow(data, interpolation='nearest', cmap='gray',
                        vmin=m - s, vmax=m + s,
                        origin='lower')

        ax1.scatter(vcpl['x'], vcpl['y'])
        plt.show()

    conv = []
    for m in vcpl['x', 'y']:
        conv.append((m[0], m[1]))

    apertures = CircularAperture(conv, r=5.0)
    #annulus_apertures = CircularAnnulus(conv, r_in=6., r_out=8.)
    #apers = [apertures, annulus_apertures]
    #p_table = aperture_photometry(data2, apers, method='subpixel', subpixels=5)
    p_table = aperture_photometry(data2, apertures, method='subpixel', subpixels=5)
    #bkg_mean = p_table['aperture_sum_1'] / annulus_apertures.area()
    #bkg_sum = bkg_mean * apertures.area()
    #final_sum = p_table['aperture_sum_0'] - bkg_sum
    #p_table['residual_aperture_sum'] = final_sum
    return p_table


if __name__ == '__main__':

    main()