
# SkyMap module
#
# Auxiliary functions to plot the SkyMap
# ____________________________
#
# This module is part of the PyASB project,
# created and maintained by Miguel Nievas [UCM].
# ____________________________


from __future__ import print_function


# NOTE: The following 2 functions should be moved to separate file or at least to a new class
# NOTE: Maybe should be rewrite as follows?:
# 1.) Create the file with the header
# 2.) Iterative add lines

class Summary(object):

    def __init__(self, image, input_options, image_analysis, instrument_calibration,
                 image_sky_brightness, cloud_coverage):
        self.summarize_results(input_options, image, image_analysis,
                               instrument_calibration, image_sky_brightness,
                               cloud_coverage)
        self.save_summary_to_file(image.image_info)

    def summarize_results(self, InputOptions, image, ImageAnalysis,
                          InstrumentCalibration, ImageSkyBrightness, CloudCoverage):

        sum_date = str(image.image_info.fits_date)
        sum_filter = str(image.image_info.used_filter)
        sum_stars = str(
            InstrumentCalibration.bouguer_fit.Regression.Nstars_initial)
        sum_gstars = str(
            "%.1f" % float(InstrumentCalibration.bouguer_fit.Regression.Nstars_rel))
        sum_zpoint = \
            str("%.3f" % float(InstrumentCalibration.bouguer_fit.Regression.mean_zeropoint)) + ' +/- ' +\
            str("%.3f" % float(
                InstrumentCalibration.bouguer_fit.Regression.error_zeropoint))
        sum_extinction = \
            str("%.3f" % float(InstrumentCalibration.bouguer_fit.Regression.extinction)) + ' +/- ' +\
            str("%.3f" % float(
                InstrumentCalibration.bouguer_fit.Regression.error_extinction))
        sum_skybrightness = \
            str("%.3f" % float(ImageSkyBrightness.SBzenith)) + ' +/- ' +\
            str("%.3f" % float(ImageSkyBrightness.SBzenith_err))
        sum_cloudcoverage = \
            str("%.3f" % float(CloudCoverage.mean_cloudcover)) + ' +/- ' +\
            str("%.3f" % float(CloudCoverage.error_cloudcover))

        self.summary_content = \
            [sum_date, sum_filter, sum_stars, sum_gstars,
             sum_zpoint, sum_extinction, sum_skybrightness, sum_cloudcoverage]

    def save_summary_to_file(self, image_info):

        spath = getattr(image_info, 'summary_path', None)

        if spath is None or spath is False:
            print('Skipping write summary to file')
        else:
            print('Write summary to file')

            content = [
                '#Date, Filter, Stars, % Good Stars, ZeroPoint, Extinction, SkyBrightness, CloudCoverage\n']
            for line in self.summary_content:
                content_line = ""
                for element in line:
                    content_line += element
                content.append(content_line + ", ")

            if spath == "screen":
                print(content)
            else:
                summary_filename = str("%s/Summary_%s_%s_%s.txt" % (
                    image_info.summary_path, image_info.obs_name,
                    image_info.fits_date, image_info.used_filter))
                summary_filename1 = "%s/Summary_%s_%s_%s.txt" % (
                    image_info.summary_path, image_info.obs_name,
                    image_info.fits_date, image_info.used_filter)

                assert summary_filename == summary_filename1

                summaryfile = open(summary_filename, 'w+')
                summaryfile.writelines(content)
                summaryfile.close()
