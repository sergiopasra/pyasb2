import numpy as np

if __name__ == '__main__':

    delta_x = -18.63912476
    delta_y = -31.06643504
    radial_factor = 14.19766968

    from astropy.wcs import WCS
    import matplotlib.pyplot as plt
    import customwcs


    customwcs = customwcs.CustomWCS()

    azimuth_zeropoint = customwcs.azimuth_zeropoint

    cdelt = customwcs.cdelt

    refx = customwcs.refx
    refy = customwcs.refy

    w0 = WCS(naxis=2)
    w0.wcs.crpix = [refx, refy]
    w0.wcs.cdelt = [cdelt, cdelt]
    w0.wcs.crval = [0, 90]
    w0.wcs.ctype = ["pLON-ZEA", "pLAT-ZEA"]

    w1 = WCS(naxis=2)
    w1.wcs.crpix = [refx, refy]
    w1.wcs.cdelt = [cdelt, cdelt]
    w1.wcs.crval = [90, 0]
    w1.wcs.ctype = ["pLAT-ZEA", "pLON-ZEA"]
    w1.wcs.lonpole = 180 - azimuth_zeropoint

    ref_az, ref_alt = customwcs.derot_1s(0, 90)
    #ref_az, ref_alt = 321.67057465, 89.180358574
    #print ref_az, ref_alt
    w2 = WCS(naxis=2)
    w2.wcs.crpix = [refx, refy]
    w2.wcs.cdelt = [cdelt, cdelt]
    w2.wcs.crval = [ref_alt, 270 + ref_az -azimuth_zeropoint]
    w2.wcs.ctype = ["pLAT-ZEA", "pLON-ZEA"]
    w2.wcs.lonpole =  -azimuth_zeropoint + ref_az
    print '---------------------'
    print refx, refy
    print cdelt, cdelt
    print ref_alt, ref_az
    print azimuth_zeropoint
    print '---------------------'

    xy2 = w1.all_world2pix([[90, 0]], 1)
    print 'C', customwcs.horiz2xy(0, 90, derotate=True)
    xy3 = w2.all_world2pix([[90, 0]], 1)
    print 'W', xy3
    xy4 = customwcs.horiz2xy(0, 90, derotate=True)
    print 'W2', xy4

    nn = np.zeros((360, 2))
    nn[:,1] = np.arange(360)
    for i in [89, 60, 40, 0]:
    #for i in [49]:
        nn[:,0] = i

        x0, y0 = customwcs.horiz2xy(nn[:,1], nn[:,0], derotate=False)

        x1, y1 = customwcs.horiz2xy(nn[:,1], nn[:,0], derotate=True)

        xy2 = w1.all_world2pix(nn, 1)
        xy3 = w2.all_world2pix(nn, 1)

        plt.plot(x0, y0, color='black')
        plt.plot(xy2[:,0], xy2[:,1], color='green')

        plt.plot(x1, y1, color='red')
        plt.plot(xy3[:,0], xy3[:,1], color='blue')

    plt.show()

    x, y = customwcs.horiz2xy(296.180713807, 49.2133806433, derotate=False)
    print x,y
    xy3 = w1.all_world2pix([[49.2133806433, 296.180713807]], 1)
    print xy3
    x, y = customwcs.horiz2xy(296.180713807, 49.2133806433, derotate=True)
    print x,y
    xy3 = w2.all_world2pix([[49.2133806433, 296.180713807]], 1)
    print xy3