import sys

sys.path.append('../')

import ConfigParser as configparser

import numpy as np

from astropy.table import Table
import astropy

from pyasb.user import main


defaults = {
    'biasname': '',
    'flatfield': '',
    'darkframe': ''
}

def apply_refraction(real_alt):
    # https://en.wikipedia.org/wiki/Atmospheric_refraction
    # This model is valid at 1 Atm and 10 C
    ang = (real_alt + 10.3 / (real_alt + 5.11))
    cot_ang = (90 - ang) / 180.0 * np.pi
    rr = 1.02 * np.tan(cot_ang) # In arc minutes
    app_alt = real_alt + rr / 60.0
    return app_alt


def correct_refraction(app_alt):
    ang = app_alt + 7.31 / (app_alt + 4.4)
    cot_ang = (90 - ang) / 180.0 * np.pi
    rr = np.tan(cot_ang) # Arc minutes
    real_alt = app_alt - rr / 60.0
    return real_alt


def create_wcs1(customwcs):
    cdelt = customwcs.cdelt

    refx = customwcs.refx
    refy = customwcs.refy

    w1 = WCS(naxis=2)
    w1.wcs.crpix = [refx, refy]
    w1.wcs.cdelt = [cdelt, cdelt]
    w1.wcs.crval = [90, 0]
    w1.wcs.ctype = ["pLAT-ZEA", "pLON-ZEA"]
    w1.wcs.lonpole = 180 - customwcs.azimuth_zeropoint

    return w1

def create_wcs2(customwcs):
    ref_az, ref_alt = customwcs.derot_1s(0, 90)
    cdelt = customwcs.cdelt

    refx = customwcs.refx
    refy = customwcs.refy

    w2 = WCS(naxis=2)
    w2.wcs.crpix = [refx, refy]
    w2.wcs.cdelt = [cdelt, cdelt]
    w2.wcs.crval = [ref_alt, ref_az]
    w2.wcs.ctype = ["pLAT-ZEA", "pLON-ZEA"]
    w2.wcs.lonpole = 180-customwcs.azimuth_zeropoint + ref_az
    return w2


def photometric_radius(magnitude, airmass, exposure, latitude, dec, base_radius=0.8):
    ''' Needs astrometry properties, photometric filter properties and image_info
        Returns R1,R2 and R3 '''

        # Returns R1,R2,R3. Needs photometric properties and astrometry.
    MF_magn = 10 ** (-0.4 * magnitude)
    MF_reso = 0.5 #0.5 * (min(resolution) / 2500)
    MF_airm = 0.7 * airmass
    if latitude >= 0:
        MF_decl = 0.2 * exposure * np.abs(1. - np.divide(dec, 90.0))
    else:
        MF_decl = 0.2 * exposure * np.abs(1. + np.divide(dec, 90.0))

    MF_totl = 1 + MF_magn + MF_reso + MF_decl + MF_airm

    R1 = base_radius * MF_totl
    R2 = R1 * 1.5 + 1
    R3 = R1 * 3.0 + 3
    return R1, R2, R3


def filter_phot(table):

    table.add_column(Column(np.zeros_like(table['vmag'], dtype='bool'), name='photo'))

    double_mask = np.logical_not(np.char.isspace(table['doub']))

    # Variable stars
    variable_mask = np.logical_not(np.char.isspace(table['var']))

    # Bad Photoemetry
    bad_phot_mask = np.logical_not(np.char.isspace(table['bad_phot_']))

    # Incomplete photometry
    phot_nan_mask = np.zeros_like(table['vmag'], dtype='bool')
    for colname in ['uv', 'bv', 'rv', 'iv']:
        phot_nan_mask = np.logical_or(phot_nan_mask, np.isnan(table[colname]))

    # photometric stars are
    photo_mask = ~double_mask & ~variable_mask & ~bad_phot_mask & ~phot_nan_mask
    table['photo'] = photo_mask

#    good_phot = table.group_by('photo')
#
#    ntable = good_phot.groups[1]
#
#    extreme_colors_mask = (ntable['bv'] < -1) & (ntable['bv'] > 2)
#    bright_stars_mask = ntable['vmag'] <= max_magnitude

    # Check this can be done
#    ntable['photo'][extreme_colors_mask] = False
#    ntable['photo'][~bright_stars_mask] = True

#    table_photometric = table
    return table


if __name__ == '__main__':

    import matplotlib.pyplot as plt
    import matplotlib.patches as ptc

    from astropy import units as u
    from astropy.coordinates import SkyCoord

    from astropy.coordinates import SkyCoord, EarthLocation, AltAz
    from astropy.coordinates import FK5
    from astropy.time import Time
    from astropy.table import Column
    import numpy
    from astropy.io import fits
    from wcsaxes import WCS
    #from astropy.wcs import WCS
    from customwcs import CustomWCS

    input_options = main(['-i', 'Johnson_V20130912_011709.fit.gz'])

    #configs = configparser.SafeConfigParser(defaults=defaults)
    #configs.read('config.ini')

    magfield = 'vmag'
    colorfields = ['uv', 'bv', 'rv', 'iv']

    min_altitude = 15
    max_magnitude = 5 # In V band, for the momment

    latitude = 40.450941
    longitude = -3.726065
    height = 667

    delta_x = -18.63912476
    delta_y = -31.06643504
    radial_factor = 14.19766968
    azimuth_zeropoint = 88.64589921
    azimuth_zeropoint = 90

    latitude = 40.450941
    longitude = -3.726065
    height = 667

    delta_x = -18.63912476
    delta_y = -31.06643504
    radial_factor = 14.19766968
    azimuth_zeropoint = 88.64589921

    customwcs = CustomWCS()

    if False:
        catfile = 'catalog.txt'
        print('load catalog from {}'.format(catfile))

        catalogrec = numpy.recfromcsv(catfile, delimiter=';', skiprows=0)

        table = Table(catalogrec)

        table2 = filter_phot(table)

        table2.write('catalog2.txt', format='ascii.csv')
    else:
        catfile = 'catalog2.txt'
        catalogrec = numpy.recfromcsv(catfile, delimiter=',', skiprows=0)
        table = Table(catalogrec)

    for fname in input_options.fits_filename_list:
        print('basic calibration')

        #final = basic_calibration(fname, configs)

        #fits.writeto('corrected.fits', final, clobber=True)

        #print table.colnames

        # Only the first 300
        table = table[:300]

        # full            | visible
        # |- Photometric  | visible

        mm = SkyCoord(ra=table['raj1950'], dec=table['dej1950'], unit=(u.hourangle, u.deg), frame=FK5(equinox='J1950'))
        location = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=height*u.m)

        time = Time("2013-09-12 01:17:09")
        mm_altz = mm.transform_to(AltAz(obstime=time,location=location))

        table.add_column(Column(mm_altz.alt.degree, name='alt'))
        table.add_column(Column(mm_altz.az.degree, name='az'))
        table.add_column(Column(mm.dec.degree, name='dec'))
        table.add_column(Column(mm.ra.degree, name='ra'))

        visibility_mask = mm_altz.alt.degree > min_altitude

        visible_catalog = table[visibility_mask]

        # And photometry

        # Handle refraction
        # Astropy does ref correction using ERFA complex models.
        # Perhaps we should use that
        # It requires Pressure, Temp, Humidity, etc
        app_altz = apply_refraction(visible_catalog['alt'])

        visible_catalog.add_column(Column(app_altz, name='alt_app'))

        vx, vy = customwcs.horiz2xy(visible_catalog['az'], visible_catalog['alt_app'], derotate=True)

        visible_catalog.add_column(Column(vx, name='x'))
        visible_catalog.add_column(Column(vy, name='y'))

        by_phot = visible_catalog.group_by('photo')
        visible_catalog_phot = by_phot.groups[1]
        print 'filtering catalog done'
        print 'photometric', len(visible_catalog_phot)
        print 'full', len(visible_catalog)

        #print visible_catalog['name','photo'][:6]

        #sys.exit(0)
        used_filter = 'Johnson_V'

        print 'recenter stars'
        from numina.array.recenter import centering_centroid, wc_to_pix_1d
        data = fits.getdata('corrected.fits')

        recentered_x = []
        recentered_y = []

        show_recenter = False

        nstars = 40

        for x, y in visible_catalog['x', 'y'][:nstars]:

            nx,ny, sb, status,_ = centering_centroid(data, x, y, box=(5, 5))

            if status == 1:
                recentered_x.append(nx)
                recentered_y.append(ny)
            else:
                recentered_x.append(x)
                recentered_y.append(y)
            #continue
            s = 10
            ixp = wc_to_pix_1d(x)
            iyp = wc_to_pix_1d(y)

            if show_recenter:
                fig = plt.figure()
                ax = fig.add_subplot(111)

                cut0 = data[iyp-s:iyp+s+1,ixp-s:ixp+s+1]

                ax.imshow(cut0, origin='lower', interpolation='nearest',
                          extent=[ixp-s-0.5, ixp+s+0.5, iyp-s-0.5, iyp+s+0.5])

                ax.add_patch(ptc.Circle((x,y), 2,
                    facecolor='none', edgecolor=(0, 0, 0.8),
                    linewidth=1, fill=False, alpha=0.5,
                    label='Origin'))
                ax.add_patch(ptc.Circle((nx, ny), 2,
                    facecolor='none', edgecolor=(0.8, 0, 0.8),
                    linewidth=1, fill=False, alpha=0.5,
                    label='Detected'))
                ax.autoscale(False)
                plt.show()

        print 'done recentering'

        data0 = fits.getdata('Johnson_V20130912_011709.fit.gz')
        fig = plt.figure()
        #ax1 = fig.add_axes([0.0, 0.0, 1.0, 1.0], projection=w1)
        ax1 = fig.add_axes([0.0, 0.0, 1.0, 1.0])
        #ax1.coords.grid(color='blue', alpha=0.5, linestyle='solid')
        ax1.imshow(data0, origin='lower', interpolation='nearest')
        ax1.autoscale(False)
        #
        ax1.plot(visible_catalog['x'][:nstars], visible_catalog['y'][:nstars], 'r*')
        ax1.plot(recentered_x, recentered_y, 'b*')

        magnitude = visible_catalog_phot['vmag']
        airmass = 1.0 / np.cos(visible_catalog_phot['alt_app'])
        exposure = 1.0
        lat = 40
        dec_vis_deg = mm[visibility_mask].dec.degree

        #R1, R2, R3 = photometric_radius(magnitude, airmass, exposure, latitude, dec_vis_deg, base_radius=0.8)

        # Looping seems the only way
        for x, y, phot in zip(recentered_x, recentered_y, visible_catalog['photo'][:nstars]):
            if phot:
                ax1.add_patch(ptc.Circle((x,y), 5,
                        facecolor='none', edgecolor='red',
                        linewidth=1, fill=False, alpha=0.5,
                        label='Origin'))
            else:
                ax1.add_patch(ptc.Circle((x,y), 5,
                        facecolor='none', edgecolor='blue',
                        linewidth=1, fill=False, alpha=0.5,
                        label='Origin'))
            if phot:
                ax1.add_patch(ptc.Circle((x,y), 9,
                        facecolor='none', edgecolor="red",
                        linewidth=1, fill=False, alpha=0.5,
                        label='Origin'))
                ax1.add_patch(ptc.Circle((x,y), 11,
                        facecolor='none', edgecolor="green",
                        linewidth=1, fill=False, alpha=0.5,
                        label='Origin'))

        plt.show()

        fig = plt.figure()
        #ax1 = fig.add_axes([0.0, 0.0, 1.0, 1.0], projection=w1)
        #ax2 = fig.add_axes([0.0, 0.0, 1.0, 1.0])

        plt.plot(recentered_x, recentered_y, 'b*')
        plt.show()

        fig = plt.figure()
        #ax1 = fig.add_axes([0.0, 0.0, 1.0, 1.0], projection=w1)
        ax1 = fig.add_subplot(111)
        #ax1.coords.grid(color='blue', alpha=0.5, linestyle='solid')


        print 'coordinate difference'
        #ax1.plot(vx, vy, 'bo')
        xx = visible_catalog['x'][:nstars]
        yy = visible_catalog['y'][:nstars]
        difx = recentered_x - xx
        dify = recentered_y - yy
        ax1.quiver(xx, yy, difx, dify, angles="xy", scale_units="xy", scale=1.0 / 50,
                   width=5e-3,)
        plt.show()


        print 'photometry'

        from photutils import CircularAperture, aperture_photometry, CircularAnnulus

        r1 = 5.0
        r2 = 9.0
        r3 = 11.0
        x = recentered_x
        y = recentered_y
        positions = zip(x,y)
        aper1 = CircularAperture(positions, r=r1)

        aper3 = CircularAnnulus(positions, r_in=r2, r_out=r3)

        rawflux_table = aperture_photometry(data, aper1)
        bkgflux_table = aperture_photometry(data, aper3)
        phot_table = astropy.table.hstack([rawflux_table, bkgflux_table], table_names=['raw', 'bkg'])

        aperture_area = np.pi * r1 ** 2
        annulus_area = np.pi * (r3 ** 2 - r2 ** 2)
        bkg_sum = phot_table['aperture_sum_bkg'] * aperture_area / annulus_area
        final_sum = phot_table['aperture_sum_raw'] - bkg_sum
        phot_table['residual_aperture_sum'] = final_sum
        print(phot_table['residual_aperture_sum'])
