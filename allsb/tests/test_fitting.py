
import numpy as np
import math
import matplotlib.pyplot as plt

from allsb.fitting import proj_rad


def proj_zen_eqa(x, y, x0, y0, scale, a0, E, eps):

    sin_eps = math.sin(eps)
    cos_eps = math.cos(eps)

    rad_theta, ang = proj_rad(x, y, scale, x0, y0)
    # print('A rad_theta', rad_theta)
    # print('A ang', ang)
    # 0 <= rad_theta <= 2
    #rarg = 1 - rad_theta**2 / 2

    #sin_u = np.clip(rarg, -1, 1)
    #cos_u = np.sqrt(1 - sin_u**2)
    # Alternative
    # 0 <= rad_theta <= 2

    # u = math.pi / 2 - 2 * np.arcsin(rad_theta / 2.0)
    u = 2 * np.arcsin(rad_theta / 2.0)
    # print('A u', u)
    sin_u = np.sin(u)
    cos_u = np.cos(u)

    b = a0 - E + ang
    # print('A b', b)
    sin_b = np.sin(b)
    cos_b = np.cos(b)

    # cosine law
    # z, angZ = math.pi + b
    cos_z = cos_u * cos_eps - sin_u * sin_eps * cos_b
    # print('A cos_z', cos_z)
    # sin(a - E) * sin(z) = sin(b) * sin(u)
    # cos(a - E) * sin(z) * sin(eps) = cos(u) - cos(eps)cos(z)
    # Developing, to avoid eps = 0
    # cos(a - E) * sin(z) = cos(u) * sin(eps) + sin(u) cos(eps) cos(b)

    sin_aE_sin_z = sin_b * sin_u
    cos_aE_sin_z = cos_u * sin_eps + sin_u * cos_eps * cos_b
    azE = np.arctan2(sin_aE_sin_z, cos_aE_sin_z)
    # print('A az - E', azE)
    az = azE + E
    # print('A az', az)
    z_dis = np.arccos(cos_z)
    # print('A z', z_dis)
    alt = math.pi / 2 - z_dis
    # print('A alt', alt)
    return alt, az


def proj_zen_eqa_inv(alt, az, x0, y0, scale, a0, E, eps):

    sin_eps = math.sin(eps)
    cos_eps = math.cos(eps)

    z_dis = math.pi / 2 - alt
    # print('B z', z_dis)
    cos_z = np.cos(z_dis)
    sin_z = np.sin(z_dis)
    azE = az - E
    # print('B az - E', azE)
    cos_azE = np.cos(azE)
    sin_azE = np.sin(azE)

    cos_u = cos_eps * cos_z + sin_z * sin_eps * cos_azE
    u = np.arccos(cos_u)
    # print('B u', u)
    # sin_u = np.sin(u)

    sin_b_sinu = sin_azE * sin_z
    cos_b_sinu = (cos_azE * sin_z - cos_u * sin_eps) / cos_eps
    b = np.arctan2(sin_b_sinu, cos_b_sinu)
    # print('B b', b)
    ang = b - a0 + E
    # print('B ang', ang)
    rad_theta = 2 * np.sin(u / 2)
    # print('B rad_theta', rad_theta)
    # There is a - (minus) sign here
    x = x0 - rad_theta / scale * np.cos(ang)
    y = y0 + rad_theta / scale * np.sin(ang)
    return x, y


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


def main():
    x0 = 534.0
    y0 = 523.0
    a0 = np.deg2rad(00.0) # South West to North
    scale = 1 / 1000.0
    E = np.deg2rad(90)
    eps = np.deg2rad(30)

    x_cop = x0
    y_cop = y0
    alt_cop, az_cop = proj_zen_eqa(x_cop, y_cop, x0, y0, scale, a0, E, eps)

    npoints = 200
    #ang = np.linspace(0, 3 * np.pi / 2 -0.5 , 200)
    #rad = 500.0
    #x = x0 + rad * np.cos(ang)
    #y = y0 + rad * np.sin(ang)
    #x = np.random.uniform(200, 800, npoints)
    #y = np.random.uniform(200, 800, npoints)
    x = [105, 505, 1005]
    y = [100, 500, 100]

    alt, az = proj_zen_eqa(x, y, x0, y0, scale, a0, E, eps)
    x2, y2 = proj_zen_eqa_inv(alt, az, x0, y0, scale, a0, E, eps)
    print(x2, y2)

    test_invert()
    return
    #f, ax = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
    f = plt.figure(figsize=(10, 5))
    ax = f.add_subplot(121, projection='polar')
    plt_catalog(ax, alt, az)
    plt_catalog(ax, alt_cop, az_cop)
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

main()