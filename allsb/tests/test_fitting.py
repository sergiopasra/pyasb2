
import lmfit
import numpy as np
import math
import matplotlib.pyplot as plt

from allsb.fitting import calc_projection_zen_eqa
from allsb.projection import proj_zen_eqa, proj_zen_eqa_inv, distance


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


def main():
    x0 = 1034.0
    y0 = 1023.0
    a0 = np.deg2rad(0.0) # South West to North
    scale = 1 / 700.0
    E = np.deg2rad(0)
    eps = np.deg2rad(0)

    x_cop = x0
    y_cop = y0
    alt_cop, az_cop = proj_zen_eqa(x_cop, y_cop, x0, y0, scale, a0, E, eps)

    npoints = 200
    #ang = np.linspace(0, 3 * np.pi / 2 -0.5 , 200)
    #rad = 500.0
    #x = x0 + rad * np.cos(ang)
    #y = y0 + rad * np.sin(ang)
    x = np.random.uniform(0, 2000, npoints)
    y = np.random.uniform(0, 2000, npoints)

    alt, az = proj_zen_eqa(x, y, x0, y0, scale, a0, E, eps)

    # test_invert()

    # rad_theta = np.linspace(0, 2, 100)
    # sin_u = 1 - rad_theta ** 2 / 2
    # u = np.arcsin(sin_u)
    # plt.plot(rad_theta, u)
    # plt.plot(rad_theta, -rad_theta + math.pi / 2)
    # plt.show()
    # return
    #f, ax = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
    f = plt.figure(figsize=(10, 5))
    ax = f.add_subplot(121, projection='polar')
    plt_catalog(ax, alt, az)
    #plt_catalog(ax, alt_cop, az_cop)
    #plt.show()
    ax2 = f.add_subplot(122)
    ax2.plot(x, y, 'r.')
    plt.show()
    # u1 = np.sin(ang)
    # u2 = np.cos(ang)
    # ang0 = np.arctan2(u1, u2)
    # nm = ang0 < 0
    # ang0[nm] = 2 * math.pi + ang0[nm]
    # plt.scatter(ang, ang0)
    # plt.show()


def simul2():
    np.random.seed(10299900)

    x0 = 1034.0
    y0 = 1101.0
    a0 = np.deg2rad(23.0) # South West to North
    scale = 1 / 700.0
    E = np.deg2rad(121)
    eps = np.deg2rad(8)

    npoints = 6
    az = np.random.uniform(0, 2 * math.pi, npoints)
    alt = np.arccos(np.random.uniform(0, 1, npoints))

    print('compute projection')
    x, y = proj_zen_eqa_inv(alt, az, x0, y0, scale, a0, E, eps)

    params = lmfit.Parameters()

    params.add('scale', value=0.5*scale, min=0, max=4*scale)
    params.add('x0', value=1000, min=0.0, max=2000.0)
    params.add('y0', value=1000, min=0, max=2000)
    params.add('a0', value=0, vary=True, min=-math.pi, max=math.pi)
    params.add('E', value=0, vary=True, min=-math.pi, max=math.pi)
    params.add('eps', value=0, vary=True, min=-math.pi / 2.0, max=math.pi / 2.0)

    print('fit projection')
    res2 = calc_projection_zen_eqa(params, x, y, alt, az)
    print('results')
    print('x0', res2.params['x0'].value, x0)
    print('y0', res2.params['y0'].value, y0)
    print('scale', res2.params['scale'].value, scale)

    print('a0', np.rad2deg(res2.params['a0'].value), np.rad2deg(a0))
    print('eps', np.rad2deg(res2.params['eps'].value), np.rad2deg(eps))
    print('E', np.rad2deg(res2.params['E'].value), np.rad2deg(E))


    c_x0 = res2.params['x0'].value
    c_y0 = res2.params['y0'].value
    c_scale = res2.params['scale'].value
    c_a0 = res2.params['a0'].value
    c_E = res2.params['E'].value
    c_eps = res2.params['eps'].value
    c_alt, c_az = proj_zen_eqa(x, y, c_x0, c_y0, c_scale, c_a0, c_E, c_eps)
    d_alt, d_az = proj_zen_eqa(x, y, x0, y0, scale, a0, E, eps)

    k0 = distance(az, alt, c_az, c_alt)
    k1 = distance(az, alt, d_az, d_alt)

    if False:
        plt.scatter(az, k0)
        plt.scatter(az, k1)
        plt.show()

        plt.scatter(alt, k0)
        plt.scatter(alt, k1)
        plt.show()


    #print('distance0', (k0 * k0).sum() / len(k0))
    #print('distance1', (k1 * k1).sum() / len(k1))

    f = plt.figure(figsize=(10, 5))
    ax = f.add_subplot(111, projection='polar')
    plt_catalog(ax, alt, az)
    ax.scatter(1.00 * c_az, 90 - np.rad2deg(c_alt), s=8)
    #ax.scatter(1.00 * d_az, 90 - np.rad2deg(d_alt), s=8)
    #plt_catalog(ax, c_alt, c_az)
    plt.show()
    return


def simul3():
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

    print('compute projection')
    x, y = proj_zen_eqa_inv(alt, az, x0, y0, scale, a0, E, eps)
    d_alt, d_az = proj_zen_eqa(x, y, x0, y0, scale, a0, E, eps)

    k1 = distance(az, alt, d_az, d_alt)

    plt.scatter(az, k1)
    plt.show()

    plt.scatter(alt, k1)
    plt.show()

    #print('distance1', (k1 * k1).sum() / len(k1))

    #f = plt.figure(figsize=(10, 5))
    #ax = f.add_subplot(111, projection='polar')
    #plt_catalog(ax, alt, az)
    # ax.scatter(1.01*c_az, 90 - np.rad2deg(c_alt), s=8)
    #ax.scatter(1.001 * d_az, 90 - np.rad2deg(d_alt), s=8)
    #plt.show()
    return


# main()

# simul2()
