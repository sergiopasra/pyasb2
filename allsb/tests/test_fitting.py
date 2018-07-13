
import numpy as np
import math
import matplotlib.pyplot as plt

from allsb.fitting import proj_rad


def proj_zen_eqa(x, y, x0, y0, scale, a0, E, eps):

    sin_eps = math.sin(eps)
    cos_eps = math.cos(eps)

    rad_theta, ang = proj_rad(x, y, scale, x0, y0)
    # 0 <= rad_theta <= 2
    #rarg = 1 - rad_theta**2 / 2

    #sin_u = np.clip(rarg, -1, 1)
    #cos_u = np.sqrt(1 - sin_u**2)
    # Alternative
    # 0 <= rad_theta <= 2

    # u = math.pi / 2 - 2 * np.arcsin(rad_theta / 2.0)
    u = 2 * np.arcsin(rad_theta / 2.0)

    sin_u = np.sin(u)
    cos_u = np.cos(u)

    b = a0 - E + ang
    sin_b = np.sin(b)
    cos_b = np.cos(b)

    cos_z = cos_u * cos_eps - sin_u * sin_eps * cos_b

    # sin(a - E) * sin(z) = sin(b) * sin(u)
    # cos(a - E) * sin(z) * sin(eps) = cos(u) - cos(eps)cos(z)
    # Developing, to avoid eps = 0
    # cos(a - E) * sin(z) = cos(u) * sin(eps) + sin(u) cos(eps) cos(b)

    sin_aE_sin_z = sin_b * sin_u
    cos_aE_sin_z = cos_u * sin_eps + sin_u * cos_eps * cos_b
    azE = np.arctan2(sin_aE_sin_z, cos_aE_sin_z)
    az = azE + E
    z_dis = np.arccos(cos_z)
    alt = math.pi / 2 - z_dis
    return alt, az


def plt_catalog(ax, alt, az):
    nom_theta = alt
    nom_phi = az
    ax.set_theta_offset(-math.pi / 2)
    ax.set_theta_direction(-1)
    ax.scatter(nom_phi, 90 - np.rad2deg(nom_theta), s=8)
    ax.set_rmax(90)
    return ax


def main():
    x0 = 534.0
    y0 = 523.0
    a0 = np.deg2rad(00.0) # South West to North
    scale = 1 / 1000.0
    E = np.deg2rad(0)
    eps = np.deg2rad(0)

    x_cop = x0
    y_cop = y0
    alt_cop, az_cop = proj_zen_eqa(x_cop, y_cop, x0, y0, scale, a0, E, eps)

    ang = np.linspace(0, 3 * np.pi / 2 -0.5 , 200)
    rad = 500.0
    x = x0 + rad * np.cos(ang)
    y = y0 + rad * np.sin(ang)

    alt, az = proj_zen_eqa(x, y, x0, y0, scale, a0, E, eps)
    
    f, ax = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
    plt_catalog(ax, alt, az)
    plt_catalog(ax, alt_cop, az_cop)
    plt.show()

    # u1 = np.sin(ang)
    # u2 = np.cos(ang)
    # ang0 = np.arctan2(u1, u2)
    # nm = ang0 < 0
    # ang0[nm] = 2 * math.pi + ang0[nm]
    # plt.scatter(ang, ang0)
    # plt.show()

main()