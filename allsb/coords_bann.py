

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

# import datetime

from allsb.reduction import reduction

import numpy
import numpy.linalg as linalg

from allsb.projection import compute


def fit_offset_and_rotation(coords0, coords1):
    """Fit a rotation and a traslation between two sets points.

    Fit a rotation matrix and a traslation bewtween two matched sets
    consisting of M N-dimensional points

    Parameters
    ----------
    coords0 : (M, N) array_like
    coords1 : (M, N) array_lke

    Returns
    -------
    offset : (N, ) array_like
    rotation : (N, N) array_like

    Notes
    ------
    Fit offset and rotation using Kabsch's algorithm[1]_ [2]_

    .. [1] Kabsch algorithm: https://en.wikipedia.org/wiki/Kabsch_algorithm

    .. [2] Also here: http://nghiaho.com/?page_id=671

    """

    cP = coords0.mean(axis=0)
    cQ = coords1.mean(axis=0)

    P0 = coords0 - cP
    Q0 = coords1 - cQ

    crossvar = numpy.dot(P0.T, Q0)

    u, s, vt = linalg.svd(crossvar)

    d = linalg.det(u) * linalg.det(vt)

    if d < 0:
        s[-1] = -s[-1]
        vt[:, -1] = -vt[:, -1]

    rot = numpy.dot(vt, u)
    off = -numpy.dot(rot, cP) + cQ

    return off, rot


def fit_rotation(coords0, coords1):
    """Fit a rotation and a traslation between two sets points.

    Fit a rotation matrix and a traslation bewtween two matched sets
    consisting of M N-dimensional points

    Parameters
    ----------
    coords0 : (M, N) array_like
    coords1 : (M, N) array_lke

    Returns
    -------
    offset : (N, ) array_like
    rotation : (N, N) array_like

    Notes
    ------
    Fit offset and rotation using Kabsch's algorithm[1]_ [2]_

    .. [1] Kabsch algorithm: https://en.wikipedia.org/wiki/Kabsch_algorithm

    .. [2] Also here: http://nghiaho.com/?page_id=671

    """

    P0 = coords0
    Q0 = coords1

    crossvar = numpy.dot(P0.T, Q0)

    u, s, vt = linalg.svd(crossvar)

    d = linalg.det(u) * linalg.det(vt)

    if d < 0:
        s[-1] = -s[-1]
        vt[:, -1] = -vt[:, -1]

    rot = numpy.dot(vt, u)

    return rot


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

    # Loading stars in image from region file
    region_name = "stars_r.reg"
    r = pyregion.open(region_name)
    names = []
    coords_x = []
    coords_y = []
    for rr in r:
        #print(rr.coord_list[:2], rr.attr[1]['text'])
        names.append(rr.attr[1]['text'].strip())
        coords_x.append(rr.coord_list[0])
        coords_y.append(rr.coord_list[1])

    # Load star catalog
    catfile = 'catalog2.txt'
    print('read table')
    table = Table.read(catfile, format='ascii.csv')
    print('filter names')
    # names = ['HD  29139']
    mask_names = np.isin(table['name'], names)
    # print(mask_names, names)
    table = table[mask_names]
    print('convert coordinates')

    #print(table)
    coords_sky = SkyCoord(
        ra=table['raj1950'],
        dec=table['dej1950'],
        unit=(u.hourangle, u.deg),
        frame=FK5(equinox='J1950')
    )
    #
    aaframe = AltAz(obstime=time, location=location,
                    temperature=10 * u.deg_C,
                    pressure=101325 * u.Pa,
                    obswl=0.5 * u.micron
                    )
    # print(aaframe)
    # star postition predictions
    coords_altz = coords_sky.transform_to(aaframe)
    print('add columns')
    table.add_column(Column(coords_altz.alt.radian, name='alt'))
    table.add_column(Column(coords_altz.az.radian, name='az'))
    table.add_column(Column(coords_altz.alt.degree, name='alt_deg'))
    table.add_column(Column(coords_altz.az.degree, name='az_deg'))
    table.add_column(Column(coords_sky.dec.degree, name='dec'))
    table.add_column(Column(coords_sky.ra.degree, name='ra'))

    print('filter columns')

    sub2 = table['name',  'alt',   'az', 'alt_deg', 'az_deg']
    #print(sub2)
    x0 = 1228.75538203
    y0 = 1227.05500371
    pol = [0, 0.0704340939421 / 180 * math.pi, 0]
    a0 = math.pi / 2 # rotaion of the plane
    E = 0 * math.pi / 2 # az of pryection
    eps = 0.001 # zenith distance of projection
    # ang_theta, ang_phi = calcm(coords_x, coords_y, x0, y0, pol, a0, E, eps)

    #print(ang_phi)
    print('--------')
    ang_theta, ang_phi = compute(coords_x, coords_y, x0, y0)
    #print(ang_phi)
    print(ang_phi, ang_theta)
    print('--------')
    nom_theta = sub2['alt']
    nom_phi = sub2['az']
    if True:
        f, axs = plt.subplots(1, 2, subplot_kw=dict(projection='polar'))
        axs[0].set_theta_offset(-math.pi / 2)
        axs[0].set_theta_direction(-1)
        axs[0].scatter(nom_phi, 90 - nom_theta / math.pi * 180)
        axs[0].set_rmax(90)
        axs[1].set_theta_offset(-math.pi / 2)
        axs[1].set_theta_direction(-1)
        axs[1].scatter(ang_phi, 90 - (ang_theta / math.pi * 180))
        axs[1].set_rmax(90)
        plt.show()

    return


def calcm(x, y, x0, y0, pol, a0, E, eps):
    x = np.asarray(x)
    y = np.asarray(y)
    xx = x - x0
    yy = y - y0
    r = np.hypot(xx, yy)

    # radial_factor = 1 / 14.19766968
    # rad_theta = radial_factor * r
    # ang_theta = np.arcsin(1 - (rad_theta / 180 * math.pi) ** 2 / 2)
    #u = math.pi / 2 - ang_theta
    u = pol[2] + r * (pol[1] * r + pol[0])

    ang = np.arctan2(yy, -xx)  # clockwise angle
    b = a0 - E + ang

    cos_z = np.cos(u) * math.cos(eps) - np.sin(u) * math.sin(eps) * np.cos(b)

    z_dis = np.arccos(cos_z)
    sin_z = np.sin(z_dis)
    sin_az_E = np.sin(b) * np.sin(u) / sin_z
    az_E = np.arcsin(sin_az_E)

    # az_Es = (np.cos(u) - math.cos(eps) * cos_z) /(math.sin(eps) * sin_z)
    az_Es = (np.cos(u) - math.cos(eps) * cos_z)
    mask = az_Es < 0
    # print('len', len(mask))
    az_E[mask] = math.pi - az_E[mask]
    az = az_E + E
    alt = math.pi / 2 - z_dis

    return alt, az


if __name__ == '__main__':

    main()