
import math

import numpy as np

from ..projection import proj_zen_eqa, proj_zen_eqa_inv, create_wcs


def gen_circle(npoints, x0, y0, rad):
    ang = np.linspace(0, 3 * np.pi / 2 -0.5 , npoints)
    x = x0 + rad * np.cos(ang)
    y = y0 + rad * np.sin(ang)
    return x, y


def test_wcs1():
    x0 = 1534.0
    y0 = 1523.0
    scale = 0.0006
    a0 = np.deg2rad(30.0)  # South West to North
    E = np.deg2rad(20)
    eps = np.deg2rad(10)

    az = np.linspace(0, 3*math.pi / 2-1, 100)
    alt = math.pi / 4 * np.ones_like(az)

    x, y = proj_zen_eqa_inv(alt, az, x0, y0, scale, a0, E, eps)

    wcs = create_wcs(x0, y0, scale, a0, E, eps)
    u, v = wcs.wcs_world2pix(np.rad2deg(alt), np.rad2deg(az), 1)

    assert np.allclose(x, u)
    assert np.allclose(y, v)


def test_wcs_inv():
    x0 = 1534.0
    y0 = 1523.0
    scale = 0.0006
    a0 = np.deg2rad(30.0)  # South West to North
    E = np.deg2rad(20)
    eps = np.deg2rad(10)
    wcs = create_wcs(x0, y0, scale, a0, E, eps)

    ang = np.linspace(0, 3*math.pi / 2-1, 100)
    x = 1000 + 200 * np.cos(ang)
    y = 1000 + 200 * np.sin(ang)

    alt, az = proj_zen_eqa(x, y, x0, y0, scale, a0, E, eps)
    altw, azw = wcs.wcs_pix2world(x, y, 1)

    altw_r = np.deg2rad(altw)
    azw_r = np.deg2rad(azw)

    # Put in [0, 2 PI] range
    az[az < 0 ] = 2 * math.pi + az[az < 0 ]

    assert np.allclose(alt, altw_r)
    assert np.allclose(az, azw_r)


def test_wcs_eq():
    from astropy.wcs import WCS

    import astropy.units as u
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation, AltAz

    m33 = SkyCoord.from_name('M33')
    bear_mountain = EarthLocation(lat=41.3 * u.deg, lon=-74 * u.deg, height=390 * u.m)
    utcoffset = -4 * u.hour  # Eastern Daylight Time
    time = Time('2012-7-12 23:00:00') - utcoffset
    frame = AltAz(obstime=time, location=bear_mountain)
    m33altaz = m33.transform_to(frame)


    print(m33)
    print(m33altaz)

    print(u)
    x0 = 2828.0618227093737
    y0 = 1875.4853799821208
    scale = 0.0007862965093399066
    a0 = 1.4352118701541068
    E = -3.0405862093935454
    eps = 0.06124805298368204

    alt_p = np.rad2deg(math.pi / 2 - eps)
    az_p = np.rad2deg(E)

    mum = SkyCoord([az_p], [alt_p], frame=frame, unit='deg')
    print(mum)
    mol = mum.transform_to('fk5')
    print(mol)
    phi_p = np.rad2deg(a0)
    s = np.rad2deg(scale)

    # Compute delta, alpha of alt_p, az_p
    #
    delta_p = 45
    alpha_p = 0

    w1 = WCS(naxis=2)
    w1.wcs.crpix = [x0, y0]
    w1.wcs.cdelt = [s, s]
    w1.wcs.crval = [delta_p, alpha_p]
    w1.wcs.ctype = ["DEC--ZEA", "RA---ZEA"]
    w1.wcs.lonpole = 180 + alpha_p - phi_p

    assert False
