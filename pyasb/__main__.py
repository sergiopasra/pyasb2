
# PyASB launcher module
#
# Concatenate processes
# ____________________________
#
# This module is part of the PyASB project,
# created and maintained by Miguel Nievas [UCM].

import sys
import inspect

from .input_options import ReadOptions
from .image_info import ImageInfo
from .help import PlatformHelp

from .load_fitsimage import FitsImage
from .astrometry import ImageCoordinates
from .bouguer_fit import BouguerFit
from .sky_brightness import SkyBrightness, SkyBrightnessGraph
from .skymap_plot import SkyMap
from .cloud_coverage import CloudCoverage, StarCatalog
from .write_summary import Summary
from .read_config import ConfigOptions
from .user import main

config_file_default = 'config.cfg'


#@profile
class LoadImage(object):

    def __init__(self, input_options, image_info,
                 config_options, configs,
                 input_file=None):
        # Load Image file list

        if input_file == None:
            input_file = input_options.fits_filename_list[0]

        # Load fits image '''
        self.fits_image = FitsImage(input_file)
        # Local copy of image_info. We will process it further.
        self.image_info = image_info
        self.image_info.read_header(self.fits_image.fits_Header)
        self.image_info.config_processing_specificfilter(config_options, configs)

        try:
            self.fits_image.subtract_corners_background = True
            self.fits_image.reduce_science_frame(
                self.image_info.darkframe,
                self.image_info.sel_flatfield,
                master_bias=None,
                image_info=self.image_info)
        except:
            print('Cannot reduce science frame')
            raise

        # Flip image if needed
        self.fits_image.flip_image_if_needed(self.image_info)

        self.fits_image.__clear__()
        self.output_paths(input_options)

    def output_paths(self, InputOptions):
        # Output file paths (NOTE: should be moved to another file or at least separated function)
        # Photometric table

        path_list = [
            "photometry_table_path", "skymap_path", "bouguerfit_path",
            "skybrightness_map_path", "skybrightness_table_path",
            "cloudmap_path", "clouddata_path", "summary_path"]

        for path in path_list:
            try:
                setattr(self.image_info, path, getattr(InputOptions, path))
            except:
                try:
                    getattr(InputOptions, path)
                except:
                    setattr(self.image_info, path, False)

#@profile


class ImageAnalysis(object):

    def __init__(self, image):
        ''' Analize image and perform star astrometry & photometry. 
            Returns image_info and StarCatalog'''
        self.StarCatalog = StarCatalog(image.image_info)

        if (image.image_info.calibrate_astrometry == True):
            image.image_info.skymap_path = "screen"
            sky_map = SkyMap(image.image_info, image.FitsImage)
            sky_map.setup_skymap()
            sky_map.set_starcatalog(self.StarCatalog)
            sky_map.astrometry_solver()

        self.StarCatalog.process_catalog_specific(
            image.FitsImage, image.image_info)
        self.StarCatalog.save_to_file(image.image_info)
        sky_map = SkyMap(image.image_info, image.FitsImage)
        sky_map.setup_skymap()
        sky_map.set_starcatalog(self.StarCatalog)
        sky_map.complete_skymap()


class InstrumentCalibration(object):

    def __init__(self, image_info, star_catalog):

        try:
            self.bouguer_fit = BouguerFit(image_info, star_catalog)
        except Exception:
            print('Cannot perform the Bouguer Fit. Error is: ')
            raise


#@profile
class MeasureSkyBrightness(object):

    def __init__(self, FitsImage, image_info, BouguerFit):
        ImageCoordinates_ = ImageCoordinates(image_info)
        TheSkyBrightness = SkyBrightness(
            FitsImage, image_info, ImageCoordinates_, BouguerFit)
        TheSkyBrightnessGraph = SkyBrightnessGraph(
            TheSkyBrightness, image_info, BouguerFit)

        '''
        TheSkyBrightness = SkyBrightness(image_info)
        TheSkyBrightness.load_mask(altitude_cut=10)
        TheSkyBrightness.load_sky_image(FitsImage)
        #TheSkyBrightness.calibrate_image(FitsImage,image_info,BouguerFit)
        TheSkyBrightness.zernike_decomposition(BouguerFit,npoints=5000,order=10)
        '''

        self.SBzenith = TheSkyBrightness.SBzenith
        self.SBzenith_err = TheSkyBrightness.SBzenith_err

#@profile


def perform_complete_analysis(input_options, image_info_common,
                              config_options, configs, input_file):


    print 'LOAD IMAGE'
    image = LoadImage(
        input_options, image_info_common, config_options, configs, input_file)

    # Look for stars that appears in the catalog, measure their fluxes. Generate starmap.

    import sys
    sys.exit(1)
    print 'IMAGE ANALYSIS'
    ImageAnalysis_ = ImageAnalysis(image)

    print('Image date: ', image.image_info.date_string,
          ', Image filter: ', image.image_info.used_filter)

    # Create the needed classes for the summary write
    class InstrumentCalibration_(object):

        class BouguerFit(object):

            class Regression(object):
                mean_zeropoint = -1
                error_zeropoint = -1
                extinction = -1
                error_extinction = -1
                Nstars_rel = -1
                Nstars_initial = -1

    try:
        # Calibrate instrument with image. Generate fit plot.
        # Clean (no leaks)
        instrumentcalibration_ = InstrumentCalibration(
            image.image_info,
            ImageAnalysis_.StarCatalog)
    except:
        class ImageSkyBrightness:
            SBzenith = '-1'
            SBzenith_err = '-1'

    else:
        # Measure sky brightness / background. Generate map.
        ImageSkyBrightness = MeasureSkyBrightness(
            image.FitsImage,
            image.image_info,
            instrumentcalibration_.bouguer_fit)

    #
    #    Even if calibration fails,
    #    we will try to determine cloud coverage
    #    and write the summary
    #

    # Detect clouds on image
    ImageCloudCoverage = CloudCoverage(
        image,
        ImageAnalysis_,
        instrumentcalibration_.bouguer_fit)

    print 'SUMMARY'
    Summary_ = Summary(image, input_options, ImageAnalysis_,
                       instrumentcalibration_, ImageSkyBrightness, ImageCloudCoverage)


if __name__ == '__main__':

    import ConfigParser as configparser

    input_options = main()
    print input_options
    #PlatformHelp_ = PlatformHelp()
    #input_options = ReadOptions(sys.argv)

    configs = configparser.SafeConfigParser(defaults={'a': '100'})
    configs.read('config.ini')

    config_options = ConfigOptions(input_options.configfile)

    image_info_common = ImageInfo()
    image_info_common.config_processing_common(config_options, input_options)

    for input_file in input_options.fits_filename_list:
        perform_complete_analysis(
            input_options, image_info_common, config_options, configs, input_file)
