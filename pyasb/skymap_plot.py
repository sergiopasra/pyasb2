
# SkyMap module
#
# Auxiliary functions to plot the SkyMap
# ____________________________
#
# This module is part of the PyASB project,
# created and maintained by Miguel Nievas [UCM].
# ____________________________


import math

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpp
import matplotlib as mpl

from .astrometry import horiz2xy, zenith_position


class SkyMap(object):

    ''' SkyMap class '''

    def __init__(self, image_info, fits_image):
        # Set image_info as local sub-object, we will use
        # it a lot.
        self.image_info = image_info
        self.star_catalog = None # Empy constructor better?
        self.stretch_data(
            fits_image.fits_data_notcalibrated,
            image_info.perc_low,
            image_info.perc_high)
        # self.setup_skymap()
        # self.draw_catalog_stars()
        # self.draw_detected_stars()
        # self.astrometry_solver()
        # self.draw_polar_axes()
        # self.show_figure()

    def setup_skymap(self):
        """
        To be executed at the beginning (no stars)
        """

        self.define_skymap()
        self.draw_skymap_data()
        self.skyfigure.canvas.draw()
        self.skyfigure.canvas.flush_events()

        plt.show(block=False)

    def complete_skymap(self):
        """
        To be executed when an astrometric solution is found
        """
        self.draw_catalog_stars()
        self.draw_detected_stars()
        self.draw_polar_axes()
        self.skyfigure.canvas.draw()
        self.skyfigure.canvas.flush_events()
        self.show_figure()

    def set_starcatalog(self, star_catalog):
        self.star_catalog = star_catalog

    def draw_catalog_stars(self):
        for star in self.star_catalog.stars_tot:
            self.draw_annotate_star(star, type=0)

    def draw_detected_stars(self):
        for star in self.star_catalog.StarList_Det:
            self.draw_annotate_star(star, type=1)
        for star in self.star_catalog.StarList_Phot:
            self.draw_annotate_star(star, type=2)

    def stretch_data(self, fits_data, pmin, pmax):
        #log_fits_data = np.log(fits_data-np.min(fits_data)+1,dtype="float32")
        log_fits_data = np.arcsinh(
            fits_data - np.min(fits_data) + 1, dtype="float32")
        valuemin = np.percentile(log_fits_data, pmin)
        valuemax = np.percentile(log_fits_data, pmax)
        self.stretched_fits_data = log_fits_data.clip(valuemin, valuemax)

    def define_skymap(self):
        """Create figure and self.skyimage subplot."""
        self.skyfigure = plt.figure(figsize=(8, 8))
        self.skyimage = self.skyfigure.add_subplot(111)
        self.skyfigure.canvas.draw()  # (block=False)

    def mouse_press_callback(self, event):
        """Coordinate input """
        if event.button == 3:
            ix, iy = event.xdata, event.ydata
            print('x = %d, y = %d' % (ix, iy))
            self.identified_stars.append(
                [self.name, self.azim, self.alti, ix, iy])
            self.star_index += 1
            self.astrometry_optimizer(full=(self.star_index > 3))
            self.scatter_stars.append(
                self.skyimage.scatter(ix, iy, marker='o', c='red', alpha=0.2))
            self.label_stars.append(
                self.skyimage.annotate(
                    self.name, xy=(ix, iy),
                    xycoords='data', xytext=(0, 3),
                    textcoords='offset points', fontsize=8, alpha=0.8))

            self.name = self.star_catalog.stars_tot[self.star_index].name
            self.azim = self.star_catalog.stars_tot[self.star_index].azimuth
            self.alti = self.star_catalog.stars_tot[
                self.star_index].altit_real
            px, py = horiz2xy(
                self.azim, self.alti, self.image_info, derotate=True)

            try:
                self.preliminary_star.remove()
            except:
                pass

            self.preliminary_star = \
                self.skyimage.scatter(
                    px, py, marker='o', c='yellow', alpha=0.5)
            print('Name: %s, Az: %s, Alt: %s' %
                  (self.name, self.azim, self.alti))
            self.skyfigure.canvas.draw()
            self.skyfigure.canvas.flush_events()

        return(None)

    def key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes:
            return(None)

        if event.key == 'n':
            print('Next star')
            self.star_index += 1
            self.name = self.star_catalog.stars_tot[self.star_index].name
            self.azim = self.star_catalog.stars_tot[self.star_index].azimuth
            self.alti = self.star_catalog.stars_tot[
                self.star_index].altit_real
            px, py = horiz2xy(
                self.azim, self.alti, self.image_info, derotate=True)
            self.preliminary_star.remove()
            self.preliminary_star = \
                self.skyimage.scatter(
                    px, py, marker='o', c='yellow', alpha=0.5)
            print('Name: %s, Az: %s, Alt: %s' %
                  (self.name, self.azim, self.alti))
        elif event.key == 'p':
            print('Previous star')
            self.preliminary_star.remove()
            self.scatter_stars[-1].remove()
            self.label_stars[-1].remove()
            self.scatter_stars.pop()
            self.label_stars.pop()
            self.identified_stars.pop()
            self.star_index -= 1
            self.name = self.star_catalog.stars_tot[self.star_index].name
            self.azim = self.star_catalog.stars_tot[self.star_index].azimuth
            self.alti = self.star_catalog.stars_tot[
                self.star_index].altit_real
            px, py = horiz2xy(
                self.azim, self.alti, self.image_info, derotate=True)
            self.preliminary_star = \
                self.skyimage.scatter(
                    px, py, marker='o', c='yellow', alpha=0.5)
            self.skyfigure.canvas.draw()
            self.skyfigure.canvas.flush_events()
            print('Name: %s, Az: %s, Alt: %s' %
                  (self.name, self.azim, self.alti))

        if event.key == 'q':
            print('End')
            self.skyfigure.canvas.mpl_disconnect(self.cid_mouse)
            self.skyfigure.canvas.mpl_disconnect(self.cid_keyboard)
            print(self.identified_stars)
            plt.close()

        self.astrometry_optimizer(full=(self.star_index > 3))

        return(None)

    def astrometry_optimizer(self, full=True):
        from scipy.optimize import minimize
        from astrometry import horiz2xy

        def horiz2xy_chi2(sol, az, alt, x, y):
            self.image_info.radial_factor = sol[0]
            self.image_info.azimuth_zeropoint = sol[1]
            if (full == True):
                self.image_info.delta_x = sol[2]
                self.image_info.delta_y = sol[3]
                self.image_info.latitude_offset = sol[4]
                self.image_info.longitude_offset = sol[5]
            else:
                self.image_info.delta_x = 0
                self.image_info.delta_y = 0
                self.image_info.latitude_offset = 0
                self.image_info.longitude_offset = 0

            xf, yf = horiz2xy(az, alt, self.image_info, derotate=True)
            return(np.sum((xf - x) ** 2 + (yf - y) ** 2))

        coords = np.array(self.identified_stars)[:, 1:]  # Remove star name
        coords = np.array(coords, dtype=float)          # Convert to float
        [_az, _alt, _x, _y] = np.transpose(coords)        # Transpose and split
        print('Solving equation system')

        if (full == True):
            initial = [10, 0, 0, 0, 0, 0]
        else:
            initial = [0, 0]

        res = minimize(
            horiz2xy_chi2, initial, args=(_az, _alt, _x, _y), tol=1e-3)

        # Fix negative radial factor
        if (res.x[0] < 0):
            res.x[0] = -res.x[0]
            res.x[1] = 180 - res.x[1]

        print(
            "Parameters (radial_factor, azimuth_zeropoint, delta_x, delta_y, lat_offset, lon_offset): ")
        print(res.x)
        print("Score [sum(dev^2)] = %.3f" %
              horiz2xy_chi2(res.x, _az, _alt, _x, _y))
        print("Success: %s" % res.success)

    def astrometry_solver(self):
        print(
            '*** star select tool. Press right-click to begin. *** \n' +
            'Right-click: assign star coords. \n' +
            'n:           next star (skip current). \n' +
            'p:           previous star (remove last entry). \n' +
            'q:           quit star select tool. \n')

        self.identified_stars = []
        self.scatter_stars = []
        self.label_stars = []
        self.star_index = 0
        self.completed = 0

        # For the northern hemisphere, put Polaris as the first star
        if self.image_info.latitude > 0:
            polaris_index = [
                star.HDcode for star in self.star_catalog.stars_tot].index("HD8890")
            AuxStar = self.star_catalog.stars_tot[polaris_index]
            self.star_catalog.stars_tot[
                polaris_index] = self.star_catalog.stars_tot[0]
            self.star_catalog.stars_tot[0] = AuxStar

        self.name = self.star_catalog.stars_tot[0].name
        self.azim = self.star_catalog.stars_tot[0].azimuth
        self.alti = self.star_catalog.stars_tot[0].altit_real
        print('Name: %s, Az: %s, Alt: %s' % (self.name, self.azim, self.alti))

        self.cid_mouse = self.skyfigure.canvas.mpl_connect(
            'button_press_event', self.mouse_press_callback)
        self.cid_keyboard = self.skyfigure.canvas.mpl_connect(
            'key_press_event', self.key_press_callback)
        plt.show(block=True)

    def draw_skymap_data(self):
        ''' Draw image '''
        self.skyimage.imshow(self.stretched_fits_data, cmap=mpl.cm.gray)

        self.skyimage.axis(
            [0, self.image_info.resolution[0], 0, self.image_info.resolution[1]])
        information = str(self.image_info.date_string) + " UTC\n" + str(self.image_info.latitude) + 5 * " " +\
            str(self.image_info.longitude) + "\n" + self.image_info.used_filter

        self.skyimage.text(0.010, 0.010, information, fontsize='small', color='white',
                           transform=self.skyimage.transAxes, backgroundcolor=(0, 0, 0, 0.75))

        plt.draw()

    def draw_polar_axes(self):
        ''' Draws meridian and altitude isolines. '''

        zenith_xy = zenith_position(self.image_info)

        for each_altitude in np.arange(0, 90, 15):
            coord_altitude_0 = horiz2xy(0, each_altitude, self.image_info)
            radius = math.sqrt(
                (coord_altitude_0[0] - zenith_xy[0]) ** 2 +
                (coord_altitude_0[1] - zenith_xy[1]) ** 2)
            self.skyimage.add_patch(
                mpp.Circle((zenith_xy[0], zenith_xy[1]), radius,
                           facecolor='k', fill=False, alpha=0.2, label='_nolegend_'))
            self.skyimage.annotate(
                str(each_altitude),
                xy=(radius + zenith_xy[0], zenith_xy[1]),
                alpha=0.2,
                fontsize=10)

        key_azimuths = {0: "N", 90: "E", 180: "S", 270: "W"}

        for each_azimuth in np.arange(0, 360, 30):
            coord_azimuth_0 = horiz2xy(each_azimuth, 0, self.image_info)
            self.skyimage.plot(
                [zenith_xy[0], coord_azimuth_0[0]],
                [zenith_xy[1], coord_azimuth_0[1]],
                color='k',
                alpha=0.2,)

            if each_azimuth in key_azimuths:
                azimuth_label = str(key_azimuths[each_azimuth])
            else:
                azimuth_label = str(each_azimuth)
            self.skyimage.annotate(
                azimuth_label,
                xy=horiz2xy(
                    each_azimuth, self.image_info.min_altitude, self.image_info),
                color='k',
                alpha=0.2,
                fontsize=10)

    def draw_annotate_star(self, Star, type=0):
        # Draw identified stars and measuring circles.
        # Annotate HD catalog code and Magnitude for each star.

        if(type == 0):
            self.skyimage.scatter(Star.Xcoord, Star.Ycoord,
                                  marker='+', c='red', alpha=0.2, label='Catalog')
            self.skyimage.annotate(
                Star.name, xy=(Star.Xcoord, Star.Ycoord),
                xycoords='data', xytext=(0, 3),
                textcoords='offset points', fontsize=8, alpha=0.8)
        elif(type == 1):
            self.skyimage.add_patch(mpp.Circle(
                (Star.Xcoord, Star.Ycoord), Star.R1,
                facecolor='none', edgecolor=(0, 0, 0.8),
                linewidth=1, fill=False, alpha=0.5,
                label='Detected'))
        elif(type == 2):
            self.skyimage.add_patch(mpp.Circle(
                (Star.Xcoord, Star.Ycoord), Star.R2,
                facecolor='none', edgecolor=(0, 0.8, 0),
                linewidth=1, fill=False, alpha=0.5,
                label='Photometric'))
            self.skyimage.add_patch(mpp.Circle(
                (Star.Xcoord, Star.Ycoord), Star.R3,
                facecolor='none', edgecolor=(0.8, 0, 0),
                linewidth=1, fill=False, alpha=0.5,
                label='_nolegend_'))
            self.skyimage.annotate(Star.FilterMag, xy=(Star.Xcoord, Star.Ycoord),
                                   xycoords='data', xytext=(0, -10),
                                   textcoords='offset points', fontsize=8)

    def show_figure(self):
        #self.skyimage.legend(('Catalog','Detected','Photometric'),loc='upper right')
        if True:#self.image_info.skymap_path == "screen":
            plt.show()
            # self.skyfigure.canvas.draw()
            # self.skyfigure.canvas.flush_events()
        else:
            skymap_filename = str("%s/SkyMap_%s_%s_%s.png" % (
                self.image_info.skymap_path, self.image_info.obs_name,
                self.image_info.fits_date, self.image_info.used_filter))

            plt.tight_layout(pad=0)
            plt.savefig(skymap_filename)

        # plt.clf()
        # plt.close('all')
