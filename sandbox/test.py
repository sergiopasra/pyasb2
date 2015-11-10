
# PyASB launcher module
#
# Concatenate processes
# ____________________________
#
# This module is part of the PyASB project,
# created and maintained by Miguel Nievas [UCM].

import sys
import inspect

from pyasb.input_options import ReadOptions
from pyasb.image_info import ImageInfo
from pyasb.help import PlatformHelp

from pyasb.load_fitsimage import FitsImage
from pyasb.astrometry import ImageCoordinates
from pyasb.bouguer_fit import BouguerFit
from pyasb.sky_brightness import SkyBrightness, SkyBrightnessGraph
from pyasb.skymap_plot import SkyMap
from pyasb.cloud_coverage import CloudCoverage, StarCatalog
from pyasb.write_summary import write_summary
from pyasb.read_config import ConfigOptions
from pyasb.user import main

config_file_default = 'config.cfg'


#@profile
class LoadImage(object):

    def __init__(self, input_options, image_info,
                 config_options, configs,
                 input_file=None):
        # Load Image file list

        if input_file is None:
            input_file = input_options.fits_filename_list[0]

        # Load fits image '''
        self.fits_image = FitsImage(input_file)
        # Local copy of image_info. We will process it further.
        self.image_info = image_info
        self.image_info.read_header(self.fits_image.fits_Header)
        self.image_info.config_processing_specificfilter(config_options, configs)

        self.fits_image.subtract_corners_background = True
        print ("Reducing image")
        self.fits_image.reduce_science_frame(
            self.image_info.darkframe,
            self.image_info.sel_flatfield,
            master_bias=None,
            image_info=self.image_info)
        print ("Done")

        # Flip image if needed

        value = configs.getboolean('pyasb', 'flip_image')
        print 'flip_image value', value
        if self.image_info.flip_image:
            print ("Flip Image")
            self.fits_image.flip_image()

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


def image_analysis(image, starcatalog):

    if image.image_info.calibrate_astrometry:
        image.image_info.skymap_path = "screen"

        sky_map = SkyMap(image.image_info, image.fits_image)

        sky_map.setup_skymap()
        sky_map.set_starcatalog(starcatalog)
        sky_map.astrometry_solver()

    starcatalog.process_catalog_specific(
        image.fits_image, image.image_info)

    starcatalog.save_to_file(image.image_info)

    print('Star Map plot ...')
    sky_map = SkyMap(image.image_info, image.fits_image)
    sky_map.setup_skymap()
    sky_map.set_starcatalog(starcatalog)
    sky_map.complete_skymap()

    return starcatalog


def calibrate_instrument(image_info, star_catalog):
    return BouguerFit(image_info, star_catalog)


class MeasureSkyBrightness(object):

    def __init__(self, FitsImage, image_info, BouguerFit):
        ImageCoordinates_ = ImageCoordinates(image_info)
        TheSkyBrightness = SkyBrightness(
            FitsImage, image_info, ImageCoordinates_, BouguerFit)
        TheSkyBrightnessGraph = SkyBrightnessGraph(
            TheSkyBrightness, image_info, BouguerFit)

        self.SBzenith = TheSkyBrightness.SBzenith
        self.SBzenith_err = TheSkyBrightness.SBzenith_err


def perform_complete_analysis(input_options, image_info_common,
                              config_options, configs, input_file):


    print 'LOAD IMAGE'
    image = LoadImage(
        input_options, image_info_common, config_options, configs, input_file)

    # Look for stars that appears in the catalog, measure their fluxes. Generate starmap.
    print 'image:', image

    print 'IMAGE ANALYSIS'

    print('Creating Star Catalog ...')
    starcatalog = StarCatalog(image.image_info)


    print('Image analysys')
    starcatalog = image_analysis(image, starcatalog)
    sys.exit()
    print('Image date: ', image.image_info.date_string,
          ', Image filter: ', image.image_info.used_filter)
    print 'HARL'
    # Calibrate instrument with image. Generate fit plot.
    # Clean (no leaks)
    instrument_calibration = calibrate_instrument(
        image.image_info,
        starcatalog)
    print 'HARL'
    # Measure sky brightness / background. Generate map.
    imageskybrightness = MeasureSkyBrightness(
        image.fits_image,
        image.image_info,
        instrument_calibration)

    #
    #    Even if calibration fails,
    #    we will try to determine cloud coverage
    #    and write the summary
    #

    # Detect clouds on image
    ImageCloudCoverage = CloudCoverage(
        image,
        starcatalog,
        instrument_calibration)

    print 'SUMMARY'
    write_summary(image, input_options, starcatalog,
                  instrument_calibration, imageskybrightness, ImageCloudCoverage)


if __name__ == '__main__':

    import ConfigParser as configparser

    input_options = main()

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
