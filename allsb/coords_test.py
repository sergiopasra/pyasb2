

"""Command line interface"""


import numpy as np
import math

def arg(x, y):
    return math.atan2(y, x)

def warp(angle):
    ang = math.fmod(angle, 2 * math.pi)
    if ang < 0:
        return ang + 2 * math.pi
    else:
        return ang

def nwarp(angles):
    # return angles
    ang = np.fmod(angles, 2 * math.pi)
    neg = ang < 0
    ang[neg] += 2 * math.pi
    return ang


def main():

    phi = 0.3
    phi_p = math.pi
    print(phi, phi_p)
    print(warp(phi - phi_p - math.pi))

    ang_theta = [1.1]
    ang_phi = [phi]
    print(ang_phi)
    ang_delt, ang_alpha = rotation(ang_theta,
                                   ang_phi,
                                   delta_p=math.pi / 2,
                                   alpha_p=0.0,
                                   phi_p=phi_p
                                   )
    print(ang_alpha)

def rotation(ang_theta, ang_phi, delta_p=math.pi / 2, alpha_p=0.0, phi_p=math.pi):
    # rotation
    ang_phi = np.asarray(ang_phi)
    t11 = np.sin(ang_theta) * np.cos(delta_p)
    t12 = np.cos(ang_theta) * np.sin(delta_p) * np.cos(ang_phi - phi_p)

    t3 = -np.cos(ang_theta) * np.sin(ang_phi - phi_p)
    print('A', t11-t12)
    print('B', t3)

    m1 = np.sin(ang_theta) * np.sin(delta_p)
    m2 = np.cos(ang_theta) * np.cos(delta_p) * np.cos(ang_phi - phi_p)

    ang_alpha = nwarp(alpha_p + np.arctan2(t3, t11 - t12))
    ang_delta = np.arcsin(m1 + m2)
    return ang_delta, ang_alpha

if __name__ == '__main__':

    main()