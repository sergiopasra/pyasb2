

import math
import numpy as np
import matplotlib.pyplot as plt
from allsb.projection import proj_zen_eqa, proj_zen_eqa_inv, distance


def generate_stars(n=15):
    az = np.random.uniform(0, 2 * math.pi, n)
    alt = np.arccos(np.random.uniform(0, 1, n))
    return alt, az


def project_points(alt, az):

    x0 = 2806.7309947662493
    y0 = 1825.2414491531129
    scale = 0.00085
    a0 = 1.4319622080500167
    E = -0.00927220010107499
    eps = -4.545596554672571e-05


    x0 = 0
    y0 = 0
    scale = 1.0 #0.00085
    a0 = 0
    E = 0
    eps = 0
    x, y = proj_zen_eqa_inv(alt, az, x0, y0, scale, a0, E, eps)
    return x, y


def project_points_alt(alt, az):

    z_dis = math.pi / 2 - alt
    rad_theta = 2 * np.sin(z_dis / 2.0)

    # The - sign here is from some weird angle convention
    x = -rad_theta * np.cos(az)
    y = rad_theta * np.sin(az)
    return x, y


def comb_index(n, k):
    from itertools import combinations, chain

    index = np.fromiter(chain.from_iterable(combinations(range(n), k)),
                        int, count=-1)
    return index.reshape(-1, k)


def projected_areas(x, y):

    # Triangles
    superindex1 = comb_index(len(x), 3)

    # method 1
    vec = np.column_stack((x, y))
    vec_re = vec[superindex1]
    # Triplets of pairs
    # Remove first vector
    ab_ac = vec_re - vec_re[:, 0::3]
    # Split AB and AC vectors
    ab = np.squeeze(ab_ac[:, 1::3])
    ac = np.squeeze(ab_ac[:, 2::3])
    # AREA = 0.5 * abs (AB x AC)
    cc = 0.5 * np.abs(np.cross(ab, ac))

    # method 2
    other_method = False
    if other_method:
        x_coords = x[superindex1]
        y_coords = y[superindex1]
        x_coords = x_coords - x_coords[:,0:1]
        y_coords = y_coords - y_coords[:,0:1]

        ab = np.column_stack((x_coords[:, 1], y_coords[:, 1]))
        ac = np.column_stack((x_coords[:, 2], y_coords[:, 2]))

        cc = 0.5 * np.abs(np.cross(ab,ac))

    return cc


def sky_areas(alt, az):
    from numpy.core.umath_tests import inner1d

    # Unit vector
    z = np.sin(alt)
    zc = np.cos(alt)
    y = zc * np.cos(az)
    x = zc * np.sin(az)
    vec = np.column_stack((x, y, z))
    #print('sky_areas', vec)
    # Triangles
    superindex1 = comb_index(len(alt), 3)
    vec_re = vec[superindex1]

    # angles between vectors

    b0 = np.squeeze(vec_re[:, 0::3])
    b1 = np.squeeze(vec_re[:, 1::3])
    b2 = np.squeeze(vec_re[:, 2::3])

    # In this case
    # inner1d(b0, b1)  == np.sum(b0 * b1, axis=1)
    angle_a = np.arctan2(np.linalg.norm(np.cross(b0, b1)), inner1d(b0, b1))
    angle_b = np.arctan2(np.linalg.norm(np.cross(b1, b2)), inner1d(b1, b2))
    angle_c = np.arctan2(np.linalg.norm(np.cross(b2, b0)), inner1d(b2, b0))

    # print('a', angle_a)
    # print('b', angle_b)
    # print('c', angle_c, b2, b0)
    s = 0.5 * (angle_a + angle_b + angle_c)

    f1 = np.tan(s / 2) * np.tan(0.5 * (s - angle_a)) * np.tan(0.5 * (s - angle_b)) * np.tan(0.5 * (s - angle_c))
    area = 4 * np.arctan(np.sqrt(f1))
    return area


def test_sky_area(alt1, az1, alt2, az2, alt3, az3):
    arc_a = np.arccos(np.sin(alt1) * np.sin(alt2) + np.cos(alt1) * np.cos(alt2) * np.cos(az1 - az2))
    arc_b = np.arccos(np.sin(alt2) * np.sin(alt3) + np.cos(alt2) * np.cos(alt3) * np.cos(az2 - az3))
    arc_c = np.arccos(np.sin(alt3) * np.sin(alt1) + np.cos(alt3) * np.cos(alt1) * np.cos(az3 - az1))
    s = 0.5 * (arc_a + arc_b +  arc_c)
    fa = np.sin(s - arc_a) / np.sin(arc_a)
    fb = np.sin(s - arc_b) / np.sin(arc_b)
    fc = np.sin(s - arc_c) / np.sin(arc_c)

    A = 2 * np.arcsin(np.sqrt(fb * fc))
    B = 2 * np.arcsin(np.sqrt(fa * fc))
    C = 2 * np.arcsin(np.sqrt(fb * fa))
    area = A + B + C - math.pi
    return area


def main():
    #np.random.seed(213827)
    N = 12
    alt, az = generate_stars(n=N)
    #alt = np.array([0, math.pi / 2, 0])
    #az = np.array([0, 0, math.pi / 2])
    #plt.plot(alt, az, '.')
    #plt.show()
    do_plot = False

    x, y = project_points(alt, az)
    x_a, y_a = project_points_alt(alt, az)
    # limit
    if do_plot:
        az_lim = np.linspace(0, 2 * math.pi)
        alt_lim = np.zeros_like(az_lim)
        alt_lim2 = -math.pi / 2 + np.zeros_like(az_lim)
        x_lim, y_lim = project_points(alt_lim, az_lim)
        x_lim2, y_lim2 = project_points(alt_lim2, az_lim)

        #plt.xlim([0, 5200])
        #plt.ylim([0, 3600])
        plt.plot(x, y, 'b.')
        plt.plot(x_a, y_a, 'g*')
        plt.plot(x_lim, y_lim, 'r-')
        plt.plot(x_lim2, y_lim2, 'r-')
        plt.show()


    cc = projected_areas(x, y)
    #plt.plot(cc, 'r.')
    #plt.show()

    #alt = [0, math.pi / 2, 0]
    #az = [0, 0, math.pi / 2]
    dd = sky_areas(alt, az)
    #plt.plot(dd, cc, 'r.')
    #plt.show()

    #print(alt)
    #print(az)
    #print('test', test_sky_area(alt[0], az[0], alt[1], az[1], alt[2], az[2]))
    #print('dd', dd)
    print(cc, dd, cc / dd, (cc / dd)**2)
    #area = 0.5 * |ax * (by-cy) + bx * (cy-ay) + cx * (ay-by)|


def one_triangle(alt, az):
    do_plot = False

    x, y = project_points(alt, az)
    x_a, y_a = project_points_alt(alt, az)
    # limit
    if do_plot:
        az_lim = np.linspace(0, 2 * math.pi)
        alt_lim = np.zeros_like(az_lim)
        alt_lim2 = -math.pi / 2 + np.zeros_like(az_lim)
        x_lim, y_lim = project_points(alt_lim, az_lim)
        x_lim2, y_lim2 = project_points(alt_lim2, az_lim)

        #plt.xlim([0, 5200])
        #plt.ylim([0, 3600])
        plt.plot(x, y, 'b.')
        plt.plot(x_a, y_a, 'g*')
        plt.plot(x_lim, y_lim, 'r-')
        plt.plot(x_lim2, y_lim2, 'r-')
        plt.show()


    cc = projected_areas(x, y)

    dd = sky_areas(alt, az)

    print(cc, dd, cc / dd, (cc / dd)**2)



#az = np.array([0, 0, math.pi / 2])

#one_triangle(np.array([0, math.pi / 2, 0]), az)
#one_triangle(np.array([0, 1.1, 0]), az)
#one_triangle(np.array([0, 0.8, 0]), az)


main()
