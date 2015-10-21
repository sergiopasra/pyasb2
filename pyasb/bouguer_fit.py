

# Bouguer fitting module
#
# Fit fluxes and star data to an extinction law to obtain
# extinction and instrument zeropoint.
# ____________________________
#
#This module is part of the PyASB project,
#created and maintained by Miguel Nievas [UCM].
# ____________________________


DEBUG = False


import inspect
import math

import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np

import pyasb.astrometry


class BouguerFit(object):

    def __init__(self, image_info, PhotometricCatalog):
        print('Calculating Instrument zeropoint and extinction ...')


        # Get Zero Point from image_info and Stars from the Catalog
        self.bouguer_fixedy(image_info)
        self.bouguer_data(PhotometricCatalog)
        # Set the default values
        self.bouguer_setdefaults(image_info)
        # Try to fit the data

        self.bouguer_fit(image_info)
        self.bouguer_plot(image_info)

        print("Bouguer extinction fit results: \n" +
              " -> C=%.3f+/-%.3f, K=%.3f+/-%.3f, r=%.3f"
              % (self.regression.mean_zeropoint, self.regression.error_zeropoint,
                 self.regression.extinction, self.regression.error_extinction,
                 self.regression.kendall_tau))

    def bouguer_data(self, StarCatalog):
        ''' Get Star data from the catalog '''
        self.xdata = np.array([Star.airmass
                               for Star in StarCatalog.StarList_Phot])
        self.ydata = np.array([Star.m25logF
                               for Star in StarCatalog.StarList_Phot])
        self.yerr = np.array([Star.m25logF_unc
                              for Star in StarCatalog.StarList_Phot])

    def bouguer_fixedy(self, image_info):
        ''' Try to get the fixed Y (zero point)'''
        #try:
        self.fixed_y = image_info.used_zero_point[0]
        self.fixed_y_unc = image_info.used_zero_point[1]
        self.yisfixed = True
        #except:
        #    print(' dont fix the Zero Point')
        #    self.yisfixed = False

    def bouguer_setdefaults(self, image_info):
        ''' Set default values (just in case that the bouguer fit fails '''
        if self.yisfixed == True:
            class Regression(object):
                mean_zeropoint = self.fixed_y
                error_zeropoint = self.fixed_y_unc
                mean_slope = 10.0
                error_slope = 10.0
                extinction = 10.0
                error_extinction = 10.0
                kendall_tau = 0.0
                Nstars_initial = 0
                Nstars_final = 0
                Nstars_rel = 0

            self.regression = Regression()
            self.can_continue = True

    def bouguer_fit(self, image_info):
        '''
        Fit measured fluxes to an extinction model
        Return regression parameters (ZeroPoint, Extinction)
        '''

        if self.yisfixed:
            self.regression = TheilSenRegression(
                Xpoints=self.xdata,
                Ypoints=self.ydata,
                image_info=image_info,
                y0=self.fixed_y,
                y0err=self.fixed_y_unc)
        else:

            self.regression = TheilSenRegression(
                Xpoints=self.xdata,
                Ypoints=self.ydata,
                image_info=image_info)


        # Apply bad point filter to data
        self.xdata = self.xdata[self.regression.badfilter]
        self.ydata = self.ydata[self.regression.badfilter]
        self.yerr = self.yerr[self.regression.badfilter]

    def bouguer_plot(self, image_info):
        if image_info.bouguerfit_path == False:
            # Don't draw anything
            print('Skipping BouguerFit Graph')
            return(None)

        # Plot photometric data from the bouguer fit

        xfit = np.linspace(
            1, pyasb.astrometry.calculate_airmass(image_info.min_altitude), 10)
        yfit = np.polyval(
            [self.regression.mean_slope, self.regression.mean_zeropoint], xfit)

        bouguerfigure = plt.figure(figsize=(8, 6))
        bouguerplot = bouguerfigure.add_subplot(111)
        bouguerplot.set_title('Bouguer extinction law fit\n', size="xx-large")
        bouguerplot.set_xlabel('Airmass')
        bouguerplot.set_ylabel(r'$m_0+2.5\log_{10}(F)$', size="large")
        bouguerplot.errorbar(
            self.xdata, self.ydata, yerr=self.yerr, fmt='*', ecolor='g')
        bouguerplot.plot(xfit, yfit, 'r-')


        plot_infotext = \
            image_info.date_string + "\n" + str(image_info.latitude) + 5 * " " + str(image_info.longitude) + "\n" +\
            image_info.used_filter + 4 * " " + "Rcorr=" + str("%.3f" % float(self.regression.kendall_tau)) + "\n" +\
            "C=" + str("%.3f" % float(self.regression.mean_zeropoint)) +\
            "+/-" + str("%.3f" % float(self.regression.error_zeropoint)) + "\n" +\
            "K=" + str("%.3f" % float(self.regression.extinction)) + "+/-"\
            + str("%.3f" % float(self.regression.error_slope)) + "\n" +\
            str("%.0f" % (self.regression.Nstars_rel)) + "% of " +\
            str(self.regression.Nstars_initial) + \
            " photometric measures shown"
        bouguerplot.text(
            0.05, 0.05, plot_infotext, fontsize='x-small', transform=bouguerplot.transAxes)


        # Show or save the bouguer plot
        if image_info.bouguerfit_path == "screen":
            plt.show()
        else:
            bouguer_filename = str("%s/BouguerFit_%s_%s_%s.png" % (
                image_info.bouguerfit_path, image_info.obs_name,
                image_info.fits_date, image_info.used_filter))
            plt.tight_layout(pad=0)
            plt.savefig(bouguer_filename, bbox_inches='tight')

        # plt.clf()
        # plt.close('all')


class TheilSenRegression(object):
    # Robust Theil Sen estimator, instead of the classic least-squares.

    def __init__(self, Xpoints, Ypoints, image_info, y0=None, y0err=None, x0=None, x0err=None):
        if (len(Xpoints) != len(Ypoints) or len(Ypoints) <= 2):
            raise ValueError("len(Xpoints) != len(Ypoints) or len(Ypoints) > 2)")
        self.Xpoints = np.array(Xpoints)
        self.Ypoints = np.array(Ypoints)
        if y0 != None:
            self.fixed_zp = True
            self.y0 = y0
            if y0err != None:
                self.y0err = y0err
            else:
                self.y0err = 0.0

            if x0 != None:
                self.x0 = x0
                if x0err != None:
                    self.x0err = 0.0
            else:
                self.x0 = 0.0
                self.x0err = 0.0
        else:
            self.fixed_zp = False
        self.Nstars_initial = len(self.Ypoints)
        self.Nstars_final = self.Nstars_initial
        self.pair_blacklist = []
        # Perform the regression
        self.perform_regression()
        # Delete bad points
        self.delete_bad_points(image_info)
        # Try to improve the regression with filtered data
        self.perform_regression()

        self.Nstars_final = sum(self.badfilter)
        self.Nstars_rel = 100. * self.Nstars_final / self.Nstars_initial

    def perform_regression(self):
        # Prepare data for regression
        self.build_matrix_values()
        self.build_complementary_matrix()
        self.build_slopes_matrix()
        self.upper_diagonal_slope_matrix_values()
        # Slope
        self.calculate_mean_slope()
        # Zero point
        self.build_zeropoint_array()
        self.calculate_mean_zeropoint()
        # Errors and fit quality
        self.calculate_residuals()
        self.calculate_kendall_tau()
        self.calculate_errors()
        if self.fixed_zp == True:
            self.mean_zeropoint = self.y0
            self.error_zeropoint = self.y0err

        self.Nstars_final = len(self.Ypoints)

    def build_matrix_values(self):
        self.X_matrix_values = \
            np.array([[column for column in self.Xpoints]
                      for line in self.Xpoints])
        self.Y_matrix_values = \
            np.array([[line for line in self.Ypoints]
                      for line in self.Ypoints])

    def build_complementary_matrix(self):
        if self.fixed_zp == False:
            self.X_complementary_values = self.X_matrix_values.transpose()
            self.Y_complementary_values = self.Y_matrix_values.transpose()
        if self.fixed_zp == True:
            self.X_complementary_values = np.array([[self.x0
                                                     for column in self.Xpoints] for line in self.Xpoints])
            self.Y_complementary_values = np.array([[self.y0
                                                     for column in self.Ypoints] for line in self.Ypoints])

    def build_slopes_matrix(self):
        self.slopes_matrix = \
            ((self.Y_matrix_values - self.Y_complementary_values + 1e-20) /
             (self.X_matrix_values - self.X_complementary_values + 1e-20))
        # +1e-20 lets us hide Numpy warning with 0/0

    def upper_diagonal_slope_matrix_values(self):
        self.upper_diag_slopes = \
            np.array([self.slopes_matrix[l][c]
                      for l in xrange(len(self.slopes_matrix))
                      for c in xrange(len(self.slopes_matrix[0])) if c > l])

    def calculate_mean_slope(self):
        self.mean_slope = np.median(self.upper_diag_slopes)
        self.extinction = -self.mean_slope

    def build_zeropoint_array(self):
        self.zeropoint_array = self.Ypoints - self.Xpoints * self.mean_slope

    def calculate_mean_zeropoint(self):
        self.mean_zeropoint = np.median(self.zeropoint_array)

    def calculate_residuals(self):
        self.residuals = self.zeropoint_array - self.mean_zeropoint

    def delete_bad_points(self, image_info):
        # 3*std_residuals threshold
        std_residual = np.std(self.residuals)
        self.badfilter = np.abs(self.residuals) < np.abs(
            image_info.lim_Kendall_tau * std_residual)
        self.Xpoints = self.Xpoints[self.badfilter]
        self.Ypoints = self.Ypoints[self.badfilter]

    def calculate_errors(self):
        xmedcuad = np.median(self.Xpoints) ** 2
        xcuaddif = self.Xpoints ** 2 - xmedcuad
        xdensity = np.sum(xcuaddif)
        sigma2_res = (1. / (self.Nstars_final - 2)) * abs(sum(self.residuals))
        sigma2_slope = sigma2_res / abs(xdensity)
        sigma2_int = sigma2_res * \
            (1. / self.Nstars_final + 1. * xmedcuad / abs(xdensity))

        self.error_slope = stats.t.ppf(
            0.975, self.Nstars_final - 2) * math.sqrt(sigma2_slope)
        self.error_zeropoint = stats.t.ppf(
            0.975, self.Nstars_final - 2) * math.sqrt(sigma2_int)
        self.error_extinction = self.error_slope

    def calculate_kendall_tau(self):
        self.kendall_tau = \
            (1. * np.sum(self.upper_diag_slopes > 0) - 1. * np.sum(self.upper_diag_slopes < 0))\
            / (1. * np.size(self.upper_diag_slopes))
