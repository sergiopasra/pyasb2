
import lmfit
import numpy as np
import math


def proj_rad(x, y, scale, x0, y0):
    x = np.asarray(x)
    y = np.asarray(y)
    xx = x - x0
    yy = y - y0
    r = scale * np.hypot(xx, yy)
    # clockwise angle
    # axis inverted
    ang = np.arctan2(yy, -xx)  # clockwise angle
    return r, ang


def proj_zen_pol(x, y, x0, y0, pol0, pol1, pol2, a0, E, eps):

    r, ang = proj_rad(x, y, 1.0, x0, y0)

    # radial_factor = 1 / 14.19766968
    # rad_theta = radial_factor * r
    # ang_theta = np.arcsin(1 - (rad_theta / 180 * math.pi) ** 2 / 2)
    #u = math.pi / 2 - ang_theta
    u = pol2 + r * (pol1 * r + pol0)

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


def proj_zen_eqa(x, y, x0, y0, scale, a0, E, eps):

    # This angle is measured clock-wise...
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

    cos_z = cos_u * math.cos(eps) - sin_u * math.sin(eps) * np.cos(b)

    # sin(a - E) * sin(z) = sin(b) * sin(u)
    # cos(a - E) * sin(z) * sin(eps) = cos(u) - cos(eps)cos(z)
    # Developing, to avoid eps = 0
    # cos(a - E) * sin(z) = cos(u) * sin(eps) + sin(u) cos(eps) cos(b)

    sin_aE_sin_z = np.sin(b) * sin_u
    cos_aE_sin_z = cos_u * math.sin(eps) + sin_u * math.cos(eps) * np.cos(b)
    az_ba = np.arctan2(sin_aE_sin_z, cos_aE_sin_z)
    az = az_ba + E
    z_dis = np.arccos(cos_z)
    alt = math.pi / 2 - z_dis
    return alt, az


def residual_zen_pol(params, x, y, alt, az):

    x0 = params['x0']
    y0 = params['y0']
    pol0 = params['pol0']
    pol1 = params['pol1']
    pol2 = params['pol2']
    a0 = params['a0']
    E = params['E']
    eps = params['eps']

    c_alt, c_az = proj_zen_pol(x, y, x0, y0, pol0, pol1, pol2, a0, E, eps)

    dist = distance(alt, az, c_alt, c_az)
    return dist


def residual_zen_eqa(params, x, y, alt, az):

    x0 = params['x0']
    y0 = params['y0']
    scale = params['scale']
    a0 = params['a0']
    E = params['E']
    eps = params['eps']

    c_alt, c_az = proj_zen_eqa(x, y, x0, y0, scale, a0, E, eps)

    dist = distance(alt, az, c_alt, c_az)
    return dist


def calc_projection_zen_pol(table_obs, x0=2000, y0=2000):

    params = lmfit.Parameters()

    params.add('pol0', value=1)
    params.add('pol1', value=0)
    params.add('pol2', value=0)

    params.add('x0', value=x0, min=x0 - 100, max=x0 + 100)
    params.add('y0', value=y0, min=y0 - 100, max=y0 + 100)
    params.add('a0', value=0.0, min=-math.pi, max=math.pi)
    params.add('E', value=0.0, min=-math.pi, max=math.pi)
    params.add('eps', value=0.0, min=-math.pi, max=math.pi)

    x = table_obs['field_x']
    y = table_obs['field_y']
    alt = table_obs['field_y']
    az = table_obs['field_y']

    lmfit.minimize(residual_zen_pol, params, args=(x, y, alt, az), iter_cb=icb)

    res = lmfit.minimize(residual_zen_pol, params, args=(x, y, alt, az), iter_cb=icb)
    return res


def calc_projection_zen_eqa(table_obs, x0=2000, y0=2000, scale=1e-4):

    x = table_obs['field_x'] #/ 1000.0
    y = table_obs['field_y'] #/ 1000.0
    alt = table_obs['field_y']
    az = table_obs['field_y']

    # FIXME: hardcoded
    rmax = 4000.0
    params = lmfit.Parameters()
    params.add('x0', value=x0 ,min=x0-100, max=x0+100)
    params.add('y0', value=y0 ,min=y0-100, max=y0+100)
    params.add('scale', value=scale, min=0, max=2 / rmax)
    # Rotation of camera
    params.add('a0', value=0.0, min=-math.pi, max=math.pi)
    # Parameters of center of projecttion
    # E azimut
    # eps zenith distance
    params.add('E', value=0.0, min=-math.pi, max=math.pi)
    params.add('eps', value=0.0, min=-math.pi, max=math.pi)

    res = lmfit.minimize(residual_zen_eqa, params, args=(x, y, alt, az), iter_cb=icb)
    return res

def icb(params, iter, resid, *args, **kws):
    # print(params)
    pass


def nwarp(angles):
    # return angles
    ang = np.fmod(angles, 2 * math.pi)
    neg = ang < 0
    ang[neg] += 2 * math.pi
    return ang


def distance(nom_theta_rad, nom_phi_rad, ang_theta_rad, ang_phi_rad):
    # radians
    t1 = np.sin(nom_theta_rad) * np.sin(ang_theta_rad)
    t2 = np.cos(nom_theta_rad) * np.cos(ang_theta_rad)
    ct = t1 + t2 * np.cos(np.abs(nom_phi_rad - ang_phi_rad))
    dist = np.arccos(ct)
    return dist


def rotation_axis(matrix):
    a1 = matrix[2,1] - matrix[1,2]
    a2 = matrix[0,2] - matrix[2,0]
    a3 = matrix[1,0] - matrix[0,1]
    return np.array([a1, a2, a3])


def rotation_angle(matrix):
    cost = 0.5 * (np.trace(matrix) - 1)
    return math.acos(cost)


def matrix_vector(ang_theta, ang_phi):
    cos_theta = np.cos(ang_theta)
    sin_theta = np.sin(ang_theta)
    cos_phi = np.cos(ang_phi)
    sin_phi = np.sin(ang_phi)
    l = cos_theta * cos_phi
    m = cos_theta * sin_phi
    n = sin_theta
    out = np.array([l, m, n])
    return out


def part1(xy, crpix_i, cdelt_i, pc_ij):

    arr = np.array(xy)
    xx = arr - crpix_i
    xx = np.dot(pc_ij, xx.T).T
#    xx = np.dot(xx, np.transpose(pc_ij))
    xx *= cdelt_i
    return xx


def part1b(xy, crpix_i, cd_ij):

    arr = np.array(xy)
    xx = arr - crpix_i
    xx = np.dot(cd_ij, xx.T).T
#    xx = np.dot(xx, np.transpose(cd_ij))
    return xx


def part2(xy):
    x2 = xy[:, 0]
    y2 = xy[:, 1]
    rad_theta = np.hypot(x2, y2) # degrees
    ang_phi = np.arctan2(x2, -y2)
    factor = np.deg2rad(rad_theta) / 2
    ang_theta = math.pi / 2 - np.arcsin(factor)
    return ang_theta, ang_phi


def part3(ang_theta, ang_phi, delta_p=math.pi / 2, alpha_p=0.0, phi_p=math.pi):
    return rotation1(ang_theta, ang_phi, delta_p, alpha_p, phi_p)


def total_transform(xy, crpix_i, cd_ij, delta_p=math.pi / 2, alpha_p=0.0, phi_p=math.pi):
    p1 = part1b(xy, crpix_i, cd_ij)
    p2 = part2(p1)
    p3 = part3(p2[0], p2[1], delta_p, alpha_p, phi_p)
    return p3


def distance_cost(nom_theta, nom_phi, coords_x, coords_y, r0, r1, delta_p, alpha_p, phi_p):
    ang_theta, ang_phi = compute(coords_x, coords_y, r0, r1)
    # rotation

    ang_delta, ang_alpha = rotation1(ang_theta, ang_phi, delta_p, alpha_p, phi_p)
    #ang_delta, ang_alpha = ang_theta, ang_phi

    d = distance(nom_theta, nom_phi, ang_delta, ang_alpha)
    return np.sqrt(d.dot(d))


def total_cost(xy, expected, crpix_i, cd_ij, delta_p=math.pi / 2, alpha_p=0.0, phi_p=math.pi):
    p1 = part1b(xy, crpix_i, cd_ij)
    p2 = part2(p1)
    p3 = part3(p2[0], p2[1], delta_p, alpha_p, phi_p)

    nom_theta = expected[0]
    nom_phi = expected[1]
    ang_delta = p3[0]
    ang_alpha = p3[1]

    d = distance(nom_theta, nom_phi, ang_delta, ang_alpha)
    return np.sqrt(d.dot(d))


def rotation1(ang_theta, ang_phi, delta_p=math.pi / 2, alpha_p=0.0, phi_p=math.pi):
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


def rotation2(ang_delta, ang_alpha, delta_p=math.pi / 2, alpha_p=0.0, phi_p=math.pi):
    # rotation
    ang_alpha = np.asarray(ang_alpha)
    t11 = np.sin(ang_delta) * np.cos(delta_p)
    t12 = np.cos(ang_delta) * np.sin(delta_p) * np.cos(ang_alpha - alpha_p)
    t3 = -np.cos(ang_delta) * np.sin(ang_alpha - alpha_p)

    m1 = np.sin(ang_delta) * np.sin(delta_p)
    m2 = np.cos(ang_delta) * np.cos(delta_p) * np.cos(ang_alpha - alpha_p)

    ang_phi = nwarp(phi_p + np.arctan2(t3, t11 - t12))
    ang_phi = nwarp(phi_p + np.arctan2(t3, t11 - t12))
    # ang_phi = phi_p + np.arctan2(t3, t11 - t12)
    ang_theta = np.arcsin(m1 + m2)
    return ang_theta, ang_phi


def compute(x, y, r0, r1):
    # Returns radians
    radial_factor = 1 / 14.19766968
    pr0 = np.asarray(x) - r0
    pr1 = np.asarray(y) - r1

    # IN equations, X axis is inverted (grows to the left)
    x2 = -radial_factor * pr0
    y2 =  radial_factor * pr1

    rad_theta = np.hypot(x2, y2) # degrees
    ang_phi = np.arctan2(x2, -y2)

    factor = np.deg2rad(rad_theta) / 2
    ang_theta = math.pi / 2 - np.arcsin(factor)

    return ang_theta, ang_phi


def compute2(x, y, r0, r1, rad, az=0.0):
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


if __name__ == '__main__':
    from astropy.wcs import WCS
    import matplotlib.pyplot as plt

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
    #plt.subplot(projection=wcs)
    #plt.scatter(res[0], res[1])
    #plt.scatter(np.rad2deg(p3[0]), np.rad2deg(p3[1]))
    # plt.grid(color='white', ls='solid')
    #plt.show()