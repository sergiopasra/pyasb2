
import lmfit
import numpy as np
import math


from .projection import proj_zen_eqa, distance, distance_hav


def calc_projection_zen_eqa(params, x, y, alt, az):

    res = lmfit.minimize(residual_zen_eqa, params, args=(x, y, alt, az),
                         iter_cb=iter_cb,
                         maxfev=1000)
    return res


def iter_cb(params, iter, resid, *args, **kws):
    pass # print(iter, resid)


def residual_zen_eqa(params, x, y, alt, az):

    x0 = params['x0']
    y0 = params['y0']
    scale = params['scale']
    a0 = params['a0']
    E = params['E']
    eps = params['eps']

    c_alt, c_az = proj_zen_eqa(x, y, x0, y0, scale, a0, E, eps)

    # dist = np.rad2deg(distance(alt, az, c_alt, c_az))
    dist = distance_hav(alt, az, c_alt, c_az)
    return dist


if __name__ == '__main__':
    from astropy.wcs import WCS
    import matplotlib.pyplot as plt
    from .projection import total_transform, total_cost

    crpix_i = [485.042114258, 645.945454915]
    cd_ij = np.multiply([[1.0,  0.0], [0.0, 1.0]], 0.08)
    delta_p = np.deg2rad(36.89199116+20)
    alpha_p = np.deg2rad(74.33805674-20)
    phi_p = 0.5*math.pi

    wcs = WCS('/home/spr/devel/github/pyasb2/data/2/IMG_1206-G_newfits.fits')

    ang = np.linspace(0, 2*math.pi, 1000)
    rad = 450.0
    xy0 = 500.0
    x = xy0 + rad * np.cos(ang)
    y = xy0 + rad * np.sin(ang)
    xy = np.column_stack((x, y))

    res = wcs.all_pix2world(x,y,1)

    p3 = total_transform(xy, crpix_i, cd_ij, delta_p, alpha_p, phi_p)

    c = total_cost(xy, np.deg2rad(res), crpix_i, cd_ij, delta_p, alpha_p, phi_p)
    print(c)
    plt.subplot(projection=wcs)
    plt.scatter(res[0], res[1])
    plt.scatter(np.rad2deg(p3[0]), np.rad2deg(p3[1]))
    # plt.grid(color='white', ls='solid')
    plt.show()
