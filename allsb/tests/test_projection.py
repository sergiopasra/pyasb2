

import numpy as np
import math
# import matplotlib.pyplot as plt

from allsb.projection import proj_zen_eqa, proj_zen_eqa_inv, distance
from allsb.projection import proj_zen_eqa0, proj_zen_eqa_inv0

alt_sample = [
1.4518938212810586,
1.3515498208376238,
 1.266300210016938,
1.3386561552058769,
 1.347948332611059,
1.4421935076565948,
 1.376632875163342,
1.3110975488849679,
1.5079876968458417,
1.3751831860325696,
1.5329246084887183,
1.1468600516042478,
1.3830787132565663,
 1.147913476394376,
1.3375205585544754,
1.3801123549078351,
1.3593522971036953,
1.4173296542991722,
 1.505529186600587,
1.3418664430420326,
1.2580656655680658,
1.3461720974133566,
             ]

az_sample = [
 0.5715867188911156,
 1.1085214304443107,
 0.5640249504875221,
 2.6198741383810686,
 1.7372265644489646,
 3.1145632519774553,
  4.760556001189695,
  5.333345075284374,
0.31420767476995665,
  5.490412608116122,
 1.1889307567288459,
 3.3485536010436947,
   1.61419039379402,
 3.4920902719370024,
  4.979293185962311,
  5.602746564129527,
  4.405296155390457,
  5.529425614141132,
  4.958366605335738,
  5.643369464514795,
 0.9285414102909622,
  1.037658812496056,
]

x_sample = [
2714.7506103515625,
  2549.21337890625,
2562.7451934814453,
2700.5631103515625,
 2541.572723388672,
 2828.535369873047,
3051.9091186523438,
3052.2183227539062,
2776.9178771972656,
2962.1609497070312,
 2765.009735107422,
 3006.855224609375,
 2578.243377685547,
 3080.108154296875,
3083.4904174804688,
2936.1072998046875,
  3074.94775390625,
2922.8573608398438,
2886.3308715820312,
2952.0812377929688,
 2461.952102661133,
2549.8931274414062,
    ]

y_sample = [
1687.4722290039062,
1714.0656127929688,
1506.7544250488281,
2073.8910522460938,
1886.1613159179688,
 1960.042724609375,
 1756.553466796875,
1574.9295043945312,
1729.2329711914062,
 1605.482192993164,
1788.0726928710938,
  2339.07373046875,
1843.7195434570312,
2306.1134033203125,
1686.0067138671875,
1595.2093048095703,
1847.6611633300781,
1643.3155212402344,
1769.2496948242188,
1547.3372039794922,
 1611.599853515625,
1693.1781921386719,
]

def plt_catalog(ax, alt, az):
    nom_theta = alt
    nom_phi = az
    ax.set_theta_offset(-math.pi / 2)
    ax.set_theta_direction(-1)
    ax.scatter(nom_phi, 90 - np.rad2deg(nom_theta), s=8)
    ax.set_rmax(90)
    return ax


def gen_random(npoints):
    x = np.random.uniform(200, 800, npoints)
    y = np.random.uniform(200, 800, npoints)
    return x, y


def gen_circle(npoints, x0, y0):
    ang = np.linspace(0, 3 * np.pi / 2 -0.5 , npoints)
    rad = 500.0
    x = x0 + rad * np.cos(ang)
    y = y0 + rad * np.sin(ang)
    return x, y


def test_invert():

    x0 = 534.0
    y0 = 523.0
    a0 = np.deg2rad(20.0) # South West to North
    scale = 1 / 1000.0
    E = np.deg2rad(90)
    eps = np.deg2rad(30)

    npoints = 200
    ang = np.linspace(0, 3 * np.pi / 2 -0.5 , npoints)
    rad = 500.0
    x = x0 + rad * np.cos(ang)
    y = y0 + rad * np.sin(ang)

    alt, az = proj_zen_eqa(x, y, x0, y0, scale, a0, E, eps)
    x2, y2 = proj_zen_eqa_inv(alt, az, x0, y0, scale, a0, E, eps)

    assert np.allclose(x, x2)
    assert np.allclose(y, y2)


def test_center():

    x0 = 534.0
    y0 = 523.0
    a0 = np.deg2rad(20.0) # South West to North
    scale = 1 / 1000.0
    E = np.deg2rad(90)
    eps = np.deg2rad(30)


    x = np.array([x0])
    y = np.array([y0])

    alt, az = proj_zen_eqa(x, y, x0, y0, scale, a0, E, eps)
    alt0, az0 = proj_zen_eqa0(a0, E, eps)

    assert np.allclose(alt, alt0)
    assert np.allclose(az, az0)
    assert np.allclose(math.pi / 2 - eps, alt0)
    assert np.allclose(E, az0)


def test_center_inv():

    x0 = 534.0
    y0 = 523.0
    a0 = np.deg2rad(20.0) # South West to North
    scale = 1 / 1000.0
    E = np.deg2rad(90)
    eps = np.deg2rad(30)


    alt = np.array([math.pi / 2])
    az = np.array([0.0])

    x, y = proj_zen_eqa_inv(alt, az, x0, y0, scale, a0, E, eps)
    xc, yc = proj_zen_eqa_inv0(x0, y0, scale, a0, E, eps)

    assert np.allclose(x, xc)
    assert np.allclose(y, yc)


def test_invert_2():
    x0 = 2806.7309947662493
    y0 = 1825.2414491531129
    scale = 0.0004270536357103684
    a0 = 1.4319622080500167
    E = -0.00927220010107499
    eps = -4.545596554672571e-05

    npoints = 200
    ang = np.linspace(0, 3 * np.pi / 2 -0.5 , npoints)
    rad = 500.0
    x = x0 + rad * np.cos(ang)
    y = y0 + rad * np.sin(ang)

    alt, az = proj_zen_eqa(x, y, x0, y0, scale, a0, E, eps)
    x2, y2 = proj_zen_eqa_inv(alt, az, x0, y0, scale, a0, E, eps)

    assert np.allclose(x, x2)
    assert np.allclose(y, y2)


def test_direct1():
    x0 = 2806.7309947662493
    y0 = 1825.2414491531129
    scale = 0.0004270536357103684
    a0 = 1.4319622080500167
    E = -0.00927220010107499
    eps = -4.545596554672571e-05

    alt = np.array(alt_sample)
    az = np.array(az_sample)

    x, y = proj_zen_eqa_inv(alt, az, x0, y0, scale, a0, E, eps)
    alt2, az2 = proj_zen_eqa(x, y, x0, y0, scale, a0, E, eps)

    neg_az2 = az2 < 0
    az2[neg_az2] = 2 * math.pi + az2[neg_az2]

    assert np.allclose(alt, alt2)
    assert np.allclose(az, az2)


def test_direct2():
    x0 = 2806.7309947662493
    y0 = 1825.2414491531129
    scale = 0.0004270536357103684
    a0 = 1.4319622080500167
    E = -0.00927220010107499
    eps = -4.545596554672571e-05

    a0 = 0
    E = 0
    eps = 0

    az = np.linspace(0, 2 * math.pi)
    alt = np.deg2rad(0 * np.ones_like(az))

    x, y = proj_zen_eqa_inv(alt, az, x0, y0, scale, a0, E, eps)

    assert True


def test_direct3():

    x00 = 2855.3553458859865
    x01 = 1852.4162698338164
    rad = 3512.439362575322 / 2.0

    npoints = 200
    ang = np.linspace(0, 3 * np.pi / 2 - 0.5, npoints)
    # print('rad is', rad, 1 / rad)
    x0 = 2855.3
    y0 = 1852.0
    scale = 0.00065
    scale = 1 / rad
    a0 = 0.0
    E = 0.0
    eps = 0.0

    x0 = 2806.7309947662493
    y0 = 1825.2414491531129
    #scale = 0.0004270536357103684
    a0 = 1.4319622080500167
    E = -0.00927220010107499
    eps = -4.545596554672571e-05

    az = np.linspace(0, 2 * math.pi)
    alt = np.deg2rad(20 * np.ones_like(az))

    x, y = proj_zen_eqa_inv(alt, az, x0, y0, scale, a0, E, eps)

    #plt.plot(x, y, '.')

    x1 = x00 + rad * np.cos(ang)
    y1 = x01 + rad * np.sin(ang)

    #plt.plot(x1, y1, '.')

    #plt.plot(x_sample, y_sample, 'r+')
    #plt.xlim([0, 5633])
    #plt.ylim([0, 3753])
    #plt.show()

    assert True


def test_simul_invert():

    x0 = 534.0
    y0 = 523.0
    a0 = np.deg2rad(20.0) # South West to North
    scale = 1 / 1000.0
    E = np.deg2rad(90)
    eps = np.deg2rad(30)

    npoints = 200
    ang = np.linspace(0, 3 * np.pi / 2 -0.5 , npoints)
    rad = 500.0
    x = x0 + rad * np.cos(ang)
    y = y0 + rad * np.sin(ang)

    alt, az = proj_zen_eqa(x, y, x0, y0, scale, a0, E, eps)
    x2, y2 = proj_zen_eqa_inv(alt, az, x0, y0, scale, a0, E, eps)

    assert np.allclose(x, x2)
    assert np.allclose(y, y2)


def test_simul3():
    np.random.seed(10299900)

    x0 = 1034.0
    y0 = 1101.0
    a0 = np.deg2rad(23.0) # South West to North
    scale = 1 / 700.0
    E = np.deg2rad(121)
    eps = np.deg2rad(8)

    npoints = 60
    az = np.random.uniform(0, 2 * math.pi, npoints)
    alt = np.arccos(np.random.uniform(0, 1, npoints))

    x, y = proj_zen_eqa_inv(alt, az, x0, y0, scale, a0, E, eps)
    d_alt, d_az = proj_zen_eqa(x, y, x0, y0, scale, a0, E, eps)

    k1 = distance(az, alt, d_az, d_alt)
    assert np.allclose(k1, 0, atol=1e-5)

