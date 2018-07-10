

"""Command line interface"""

import sys
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
import astropy.table
from astropy.wcs import WCS
from astropy.io import fits

# import datetime


def table_reg(pyreg):
    names = []
    coords_x = []
    coords_y = []
    for rr in pyreg:
        names.append(rr.attr[1]['text'].strip())
        coords_x.append(rr.coord_list[0])
        coords_y.append(rr.coord_list[1])

    print('create table from inputs')
    table_name = Table()
    table_name['name'] = names
    table_name['x'] = coords_x
    table_name['y'] = coords_y
    return table_name

def filter_catalog(table_obs, location, time, min_altitude=25.0, min_magnitude=6.0):
    print('convert coordinates')

    coords_sky = SkyCoord(
        ra=table_obs['raj1950'],
        dec=table_obs['dej1950'],
        unit=(u.hourangle, u.deg),
        frame=FK5(equinox='J1950')
    )
    #
    aaframe = AltAz(obstime=time, location=location,
                    temperature=10 * u.deg_C,
                    pressure=101325 * u.Pa,
                    obswl=0.5 * u.micron
                    )

    # star position predictions
    coords_altz = coords_sky.transform_to(aaframe)

    table_obs.add_column(Column(coords_altz.alt.radian, name='alt'))
    table_obs.add_column(Column(coords_altz.az.radian, name='az'))
    table_obs.add_column(Column(coords_altz.alt.degree, name='alt_deg'))
    table_obs.add_column(Column(coords_altz.az.degree, name='az_deg'))
    table_obs.add_column(Column(coords_sky.dec.degree, name='dec'))
    table_obs.add_column(Column(coords_sky.ra.degree, name='ra'))

    visibility_mask = table_obs['alt_deg'] > min_altitude
    visible_catalog = table_obs[visibility_mask]

#    by_phot = table.group_by('photo')

#    catalog_phot = by_phot.groups[1]
    visible_catalog_phot = visible_catalog

    current_mag = visible_catalog_phot['vmag']
    magnitude_mask = current_mag < min_magnitude
    cpl = visible_catalog_phot[magnitude_mask]
    return cpl


def table_xy(table, wcs):
    coord2 = np.array([table['alt_deg'], table['az_deg']]).T
    x2 = wcs.wcs_world2pix(np.asarray(coord2), 1)
    result = Table()
    result['name'] = table['name']
    result['x'] = x2[:,0]
    result['y'] = x2[:,1]
    return result

def table_to_region(table, filename, color='green'):
    pre = '# Region file format: DS9 version 4.1 global color=green \
    dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 \
    dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 physical\n'
    with open(filename, 'w') as fd:
        fd.write(pre)
        for n,x,y in table:
            fd.write('circle({},{},8) # color={} text={{ {} }}\n'.format(x,y, color, n))



def nwarp(angles):
    # return angles
    ang = np.fmod(angles, 2 * math.pi)
    neg = ang < 0
    ang[neg] += 2 * math.pi
    return ang


def main(argv=None):

    parser = argparse.ArgumentParser()
    loc_group = parser.add_argument_group('location')
    loc_group.add_argument('--location-latitude', type=float)
    loc_group.add_argument('--location-longitude', type=float)
    loc_group.add_argument('--location-height', type=float)
    loc_group.add_argument('--location-timestamp')

    parser.add_argument('--catalog', default='catalog2.txt')
    parser.add_argument('region_file')
    parser.add_argument('--fits', default='test.fits')

    args = parser.parse_args(args=argv)
    print(args)

    if (args.location_latitude is None
            or args.location_longitude is None
            or args.location_height is None):
        # location is not set from CLI, error
        print("location is not set from CLI, error")
        sys.exit(1)

    # Create time
    if args.location_timestamp is not None:
        # time = Time("2013-09-12 01:17:09")
        time = Time(args.location_timestamp)
    else:
        # get it from header, not now
        print("get time from header, not now")
        sys.exit(1)

    # Location Madrid
    # latitude = 40.450941
    # longitude = -3.726065
    # height = 667

    # Location VdD
    # latitude = 41.00222222222222
    # longitude = -2.4833333333333334
    # height = 1149

    latitude = args.location_latitude
    longitude = args.location_longitude
    height = args.location_height

    location = EarthLocation(
        lat=latitude * u.deg,
        lon=longitude * u.deg,
        height=height * u.m
    )

    print('location', location)
    print('time', time)

    # Loading stars in image from region file
    region_name = args.region_file

    print('load region file')

    pyreg = pyregion.open(region_name)
    table_name = table_reg(pyreg)

    #
    # Load star catalog
    catfile = args.catalog
    print('read catalog')
    table = Table.read(catfile, format='ascii.csv')

    print('join on names')
    table_obs = astropy.table.join(table, table_name)
    print('done')

    print('convert coordinates')
    #print(table)
    coords_sky = SkyCoord(
        ra=table_obs['raj1950'],
        dec=table_obs['dej1950'],
        unit=(u.hourangle, u.deg),
        frame=FK5(equinox='J1950')
    )
    #
    aaframe = AltAz(obstime=time, location=location,
                    temperature=10 * u.deg_C,
                    pressure=101325 * u.Pa,
                    obswl=0.5 * u.micron
                    )

    # star position predictions
    coords_altz = coords_sky.transform_to(aaframe)
    print('add columns')
    table_obs.add_column(Column(coords_altz.alt.radian, name='alt'))
    table_obs.add_column(Column(coords_altz.az.radian, name='az'))
    table_obs.add_column(Column(coords_altz.alt.degree, name='alt_deg'))
    table_obs.add_column(Column(coords_altz.az.degree, name='az_deg'))
    table_obs.add_column(Column(coords_sky.dec.degree, name='dec'))
    table_obs.add_column(Column(coords_sky.ra.degree, name='ra'))

    print('filter columns')

    sub2 = table_obs['name',  'alt',   'az', 'alt_deg', 'az_deg', 'dej1950', 'x', 'y']
    coords_x = sub2['x']
    coords_y = sub2['y']

    # Initial astrometry
    rad_pix = 1228.18
    rad_ang = np.rad2deg(math.sqrt(2))
    rad_scale = rad_ang / rad_pix

    xcenter = 1934
    ycenter = 1262

    # Offset angles
    # For plotting
    plt_n_off_deg = -180 - 28.0
    plt_n_off = np.deg2rad(plt_n_off_deg)

    # For WCS
    az_pole = np.deg2rad(270 + 28.0)
    # WCS lonpole is az_pole - 90
    # LONPOLE is 180 + 28.0
    # plt_n_off_def is -LONPOLE
    #
    # For "compute" function
    az2_pole_deg = -90 - 28.0
    az2_pole = np.deg2rad(az2_pole_deg)

    pbase = [xcenter, ycenter, rad_scale, 0.0, 0.0, az_pole]
    # Initial WCS
    wcs0 = create_wcs2(pbase)

    coord2 = np.array([sub2['alt_deg'], sub2['az_deg']]).T
    # np.savetxt('coords2.txt', x2)
    hdr0 = wcs0.to_header()
    # nheader2.extend(h2)
    hdr0.tofile('header0.txt', overwrite=True)

    if False:
        x2 = wcs0.wcs_world2pix(np.asarray(coord2), 1)
        circle1 = plt.Circle((xcenter, ycenter), rad_pix, fill=False, color='r')
        #circle2 = plt.Circle((0.5, 0.5), 0.2, color='blue')
        #circle3 = plt.Circle((1, 1), 0.2, color='g', clip_on=False)

        fig, ax = plt.subplots()  # note we must use plt.subplots, not plt.subplot
        # (or if you have an existing figure)
        # fig = plt.gcf()
        # ax = fig.gca()

        ax.add_artist(circle1)
        #ax.add_artist(circle2)
        #ax.add_artist(circle3)
        plt.axis('equal')
        plt.xlim([0, 3866])
        plt.ylim([0, 2574])
        ax.scatter(coords_x, coords_y)
        ax.scatter(x2[:,0], x2[:,1], color='r')
        plt.show()

    nom_theta = sub2['alt']
    nom_phi = sub2['az']

    if False:
        f, axs = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
        axs.set_theta_offset(plt_n_off)
        axs.set_theta_direction(-1)
        axs.scatter(nom_phi, 90 - nom_theta / math.pi * 180)
        axs.set_rmax(90)
        plt.show()

    print('Scale is', rad_scale, 'deg / pix')
    print('center', xcenter, ycenter)
    print('AZ', np.rad2deg(az_pole))

    # angles of stars from astrometry
    p0 = [xcenter, ycenter, az_pole]

    r0 = p0[0]
    r1 = p0[1]
    #az2_pole = - math.pi / 2 - 0.45
    #az2_pole_deg = np.rad2deg(az2_pole)
    print('AZs', az2_pole_deg)
    # radians
    nom_theta = sub2['alt']
    nom_phi = sub2['az']

    ang_theta, ang_phi = compute(coords_x, coords_y, r0, r1, rad_scale, az=az2_pole)
    # ang_theta_w, ang_phi_w = compute(coords_x, coords_y, r0, r1, rad_scale)

    if False:
        f, axs = plt.subplots(1, 2, subplot_kw=dict(projection='polar'))
        axs[0].set_theta_offset(plt_n_off)
        axs[0].set_theta_direction(-1)
        axs[0].scatter(nom_phi, 90 - nom_theta / math.pi * 180)
        axs[0].set_rmax(90)
        axs[1].set_theta_offset(plt_n_off)
        axs[1].set_theta_direction(-1)
        axs[1].scatter(ang_phi, 90 - (ang_theta / math.pi * 180))
        axs[1].set_rmax(90)
        plt.show()

    join = np.array([coords_x, coords_y]).T
    calc0 = wcs0.wcs_pix2world(join, 1)

    if False:
        f, axs = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
        axs.set_theta_offset(plt_n_off)
        axs.set_theta_direction(-1)
        axs.scatter(ang_phi[:1], 90 - (ang_theta[:1] / math.pi * 180))
        axs.scatter(nom_phi[:1], 90 - (nom_theta[:1] / math.pi * 180), color='g')
        #axs.scatter(np.deg2rad(calc0[:,1]), 90 - calc0[:,0], color='r')
        axs.set_rmax(90)
        plt.show()

    dist = distance(nom_theta, nom_phi, ang_theta, ang_phi)
    d = np.array(dist)
    #print(d)
    print(np.array(nom_theta))
    print(ang_theta)
    print('dist', np.sqrt(dist.dot(dist)))

    # rotation
    # ang_delta, ang_alpha = rotation(ang_theta, ang_phi, delta_p=math.pi / 2, alpha_p=0.0, phi_p=math.pi)
    # dis = distance(ang_delta, ang_alpha, ang_theta, ang_phi)
    # for a1,b1,a2,b2,d in zip(ang_delta, ang_alpha, ang_theta, ang_phi, dis):
    #    print('_A_', a1, a2, '_B_', b1, b2, d)

    #d = total(nom_theta, nom_phi, coords_x, coords_y, r0, r1,
    #          delta_p=math.pi/2, alpha_p=0.0, phi_p=math.pi)
    #print('before', d)

    from scipy.optimize import minimize

    # OLD
    # p0 = [1229, 1226, scale, 0.0, 0.0, 3.16027342]
    p0 = pbase #[xcenter, ycenter, scale, 0.0, 0.0, 3.16027342]
    pbase2 = [xcenter, ycenter, rad_scale, 0.1, 0.1, az_pole]

    def mini_func(pr, xdata):

        r0, r1, scale, z_delta_p, alpha_p, phi_p = pr
        print(z_delta_p)
        return mini_func_base(xdata, r0, r1, scale, z_delta_p, alpha_p, phi_p)

    M = len(nom_theta)
    xdata = np.empty((4, M))
    xdata[0] = nom_theta
    xdata[1] = nom_phi
    xdata[2] = coords_x
    xdata[3] = coords_y

    pbase2 = [1900.0, 1250, 0.05, 0.1, 0.1, 0.1]

    res = minimize(mini_func, pbase2, args=(xdata,),
                   bounds=[
                       (1800, 2000),
                       (1100, 1400),
                       (0, 0.1),
                       (0, 0.1),
                       (-math.pi, math.pi),
                       (0, 2 * math.pi),
                    ]
                   )
    print(res)
    print('RES', res.success)
    print('RES', res.x)
    print('x0, y0', res.x[0], res.x[1])
    print('scale', res.x[2])
    print('angles', np.rad2deg(res.x[3:6]))
    print('-----------------------')
    print('x0, y0', p0[0], p0[1])
    print('scale', p0[2])
    print('angles', np.rad2deg(p0[3:6]))

    wcs1 = create_wcs2(res.x)
    data, header = fits.getdata(args.fits, header=True)

    nheader1 = fits.Header(header)
    h1 = wcs1.to_header()
    nheader1.extend(h1)
    h1.tofile('header1.txt', overwrite=True)
    fits.writeto('test-w1.fits', data, header=nheader1, overwrite=True)

    # return region file
    table_ex = filter_catalog(table, location, time, min_magnitude=3.0, min_altitude=25)
    #print(table)
    xy_vals = table_xy(table_ex, wcs1)
    table_to_region(xy_vals, "dum1.reg", color='blue')

    xy_vals = table_xy(table_ex, wcs0)
    table_to_region(xy_vals, "dum0.reg", color='green')


    return

    join = np.array([coords_x, coords_y]).T
    coord2 = np.array([sub2['alt_deg'], sub2['az_deg']]).T

    # Dont fit
    # w2 = create_wcs2(res.x)
    w2 = create_wcs2(pbase)

    data, header = fits.getdata(args.fits, header=True)
    print(w2)
    nheader2 = fits.Header(header)
    h2 = w2.to_header()
    nheader2.extend(h2)
    h2.tofile('header1.txt', overwrite=True)
    fits.writeto('test-w1.fits', data, header=nheader2, overwrite=True)

    x2 = w2.wcs_world2pix(np.asarray(coord2), 1)
    np.savetxt('coords1.txt', x2)
    return

    m2 = w2.wcs_pix2world(join, 1)

    plt.scatter(sub2['az_deg'],sub2['alt_deg']- m2[:,0])
    plt.show()

    plt.scatter(sub2['alt_deg'],sub2['alt_deg']- m2[:,0])
    plt.show()

    if True:
        fig = plt.figure()
        ax1 = fig.add_axes([0.0, 0.0, 1.0, 1.0], projection=w1)
        ax1.coords.grid(color='blue', alpha=1, linestyle='solid')
        overlay1 = ax1.get_coords_overlay(w1)
        overlay1.grid(color='white', linestyle='solid', alpha=1)

        m, s = np.mean(data_sub), np.std(data_sub)
        im = ax1.imshow(data_sub, interpolation='nearest', cmap='gray',
                    vmin=m - s, vmax=m + s, origin='lower')
        plt.show()

    if True:
        fig = plt.figure()
        ax1 = fig.add_axes([0.0, 0.0, 1.0, 1.0], projection=w2)
        ax1.coords.grid(color='blue', alpha=1, linestyle='solid')
        overlay1 = ax1.get_coords_overlay(w2)
        overlay1.grid(color='white', linestyle='solid', alpha=1)

        m, s = np.mean(data_sub), np.std(data_sub)
        im = ax1.imshow(data_sub, interpolation='nearest', cmap='gray',
                    vmin=m - s, vmax=m + s, origin='lower')
        plt.show()


def mini_func_base(xdata, r0, r1, scale, z_delta_p, alpha_p, phi_p):
    nom_theta = xdata[0]
    nom_phi = xdata[1]
    coords_x = xdata[2]
    coords_y = xdata[3]
    # returns radians
    delta_p = math.pi / 2 - z_delta_p

    #print('S', r0, r1, np.rad2deg([delta_p, alpha_p, phi_p, alpha_p + phi_p]))

    ang_theta, ang_phi = compute(coords_x, coords_y, r0, r1, scale)
    ang_delta, ang_alpha = rotation(ang_theta, ang_phi, delta_p, alpha_p, phi_p)
    d = distance(nom_theta, nom_phi, ang_delta, ang_alpha)
    return np.sum(d.dot(d))


def fit_func(xdata, r0, r1, z_delta_p, alpha_p, phi_p):
    nom_theta = xdata[0]
    nom_phi = xdata[1]
    coords_x = xdata[2]
    coords_y = xdata[3]
    # returns radians
    delta_p = math.pi / 2 - z_delta_p

    #print('R', r0, r1, np.rad2deg([delta_p, alpha_p, phi_p, alpha_p + phi_p]))

    ang_theta, ang_phi = compute(coords_x, coords_y, r0, r1)

    ang_delta, ang_alpha = rotation(ang_theta, ang_phi, delta_p, alpha_p, phi_p)

    d = distance(nom_theta, nom_phi, ang_delta, ang_alpha)
    return np.sum(d.dot(d))


def rotation(ang_theta, ang_phi, delta_p=math.pi / 2, alpha_p=0.0, phi_p=math.pi):
    # rotation
    ang_phi = np.asarray(ang_phi)
    t11 = np.sin(ang_theta) * np.cos(delta_p)
    t12 = np.cos(ang_theta) * np.sin(delta_p) * np.cos(ang_phi - phi_p)
    t3 = -np.cos(ang_theta) * np.sin(ang_phi - phi_p)

    m1 = np.sin(ang_theta) * np.sin(delta_p)
    m2 = np.cos(ang_theta) * np.cos(delta_p) * np.cos(ang_phi - phi_p)

    ang_alpha = nwarp(alpha_p + np.arctan2(t3, t11 - t12))
    # ang_alpha = alpha_p + np.arctan2(t3, t11 - t12)
    ang_delta = np.arcsin(m1 + m2)
    return ang_delta, ang_alpha


def total(nom_theta, nom_phi, coords_x, coords_y, r0, r1, delta_p, alpha_p, phi_p):
    ang_theta, ang_phi = compute(coords_x, coords_y, r0, r1)
    # rotation

    ang_delta, ang_alpha = rotation(ang_theta, ang_phi, delta_p, alpha_p, phi_p)
    #ang_delta, ang_alpha = ang_theta, ang_phi

    d = distance(nom_theta, nom_phi, ang_delta, ang_alpha)
    return np.sqrt(d.dot(d))


def compute(x, y, r0, r1, rad, az=0.0):
    # Returns radians
    # radial_factor = 1 / 14.19766968
    radial_factor = rad
    pr0 = np.asarray(x) - r0
    pr1 = np.asarray(y) - r1

    # IN equations, X axis is inverted (grows to the left)
    x2 = -radial_factor * pr0
    y2 =  radial_factor * pr1

    rad_theta = np.hypot(x2, y2) # degrees
    ang_phi = np.arctan2(x2, -y2) + az

    factor = 1 - (rad_theta / 180 * math.pi) ** 2 / 2.0
    ang_theta = np.arcsin(factor)

    return ang_theta, ang_phi


def distance(nom_theta_rad, nom_phi_rad, ang_theta_rad, ang_phi_rad):
    # radians
    t1 = np.sin(nom_theta_rad) * np.sin(ang_theta_rad)
    t2 = np.cos(nom_theta_rad) * np.cos(ang_theta_rad)
    ct = t1 + t2 * np.cos(np.abs(nom_phi_rad - ang_phi_rad))
    dist = np.arccos(ct)
    return dist


def create_wcs2(pr):
    r0 = pr[0]
    r1 = pr[1]
    ref_delta_p = math.pi / 2 - pr[3] # field 3 is zenith distance
    ref_alpha_p = pr[4]
    ref_phi_p = pr[5]
    delta_p = np.rad2deg(ref_delta_p)
    alpha_p = np.rad2deg(ref_alpha_p)
    phi_p = np.rad2deg(ref_phi_p)
    #m = total(nom_theta, nom_phi, coords_x, coords_y, r0, r1, math.pi / 2, 0, ref_phi_p)
    #print('after', m)

    cdeltx, cdelty = pr[2], pr[2]

    w1 = WCS(naxis=2)
    w1.wcs.crpix = [r0, r1]
    w1.wcs.cdelt = [cdeltx, cdelty]
    w1.wcs.crval = [delta_p, alpha_p]
    w1.wcs.ctype = ["pLAT-ZEA", "pLON-ZEA"]
    w1.wcs.lonpole = phi_p - 90
    return w1


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


if __name__ == '__main__':

    main()