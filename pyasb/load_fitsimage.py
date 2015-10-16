
# Load FITS image and header
#
# This module loads the AllSky FITS image and returns both
# the Image binary data and the plain-text header.
# ____________________________
#
# This module is part of the PyASB project,
# created and maintained by Miguel Nievas [UCM].
# ____________________________


import os
import numpy as np
import warnings

import astropy.io.fits as fits
from .astrometry import ImageCoordinates


class ImageTest(object):

    '''Perform some test on the image header and extract information'''

    @staticmethod
    def correct_exposure(file_header):

        return file_header['EXPOSURE']

    @staticmethod
    def correct_date(file_header):
        # Date and time

        return file_header['DATE']


    @staticmethod
    def correct_resolution(file_header):
        resolution = [
            int(file_header['NAXIS1']), int(file_header['NAXIS2'])]
        return resolution


    @staticmethod
    def correct_filter(file_header):
        # Test if there's a known filter

        used_filter = file_header['FILTER']
        return used_filter


class FitsImage(ImageTest):

    def __init__(self, input_file):
        self.load_science(input_file)
        # Backup original data
        print('Backup original (non-calibrated) data')
        self.fits_data_notcalibrated = self.fits_data[:]

    def load_science(self, input_file):
        print('Loading ScienceFrame [{}]...'.format(input_file))

        file_opened = fits.open(input_file)
        self.fits_data = file_opened[0].data
        self.fits_Header = file_opened[0].header
        self.fits_Texp = float(
            ImageTest.correct_exposure(self.fits_Header))

        print('OK')

    def load_dark(self, MasterDark):
        print('Loading MasterDark ...'),

        MasterDark_HDU = fits.open(MasterDark)
        self.MasterDark_Data = MasterDark_HDU[0].data
        self.MasterDark_Header = MasterDark_HDU[0].header
        self.MasterDark_Texp = MasterDark_HDU[0].header['EXPOSURE']
        return MasterDark_HDU
        print('OK')

    def load_flat(self, masterflat_filename):

        with fits.open(masterflat_filename) as hdul:
            mflatdata = hdul[0].data
            mflatdata_norm = mflatdata / np.mean(mflatdata)
            return mflatdata_norm

    def load_bias(self, MasterBias):
        print('Loading MasterBias ...'),

        MasterBias_HDU = fits.open(MasterBias)
        self.MasterBias_Data = MasterBias_HDU[0].data
        self.MasterBias_Header = MasterBias_HDU[0].header
        self.MasterBias_Texp = float(
            ImageTest.correct_exposure(self.MasterBias_Header))
        return MasterBias_HDU

    def reduce_science_frame(self, master_dark=None, master_flat=None, master_bias=None, image_info=None):
        '''
        Load MasterDark and MasterFlat. MasterBias is neccesary only if working
        with different exposures between Dark and Science frames
        '''

        skip_dark = False
        skip_flat = True
        mdarkdata = 0.0
        mflatdata = 1.0

        print image_info
        # Load FLAT Field

        if master_flat is not None:
            try:
                print('Loading MasterFlat ...'),
                mflatdata = self.load_flat(master_flat)
                skip_flat = False
                print('OK')
            except StandardError:
                warnings.warn('MasterFlat cannot be loaded, SKIP the flat calibration', RuntimeWarning)
                skip_flat = True

            # Load DARK Frame

        if master_dark is None:
            skip_dark = True

        try:
            mdark_hdul = self.load_dark(master_dark)
            mdarkdata = mdark_hdul[0].data
        except:
            # Try to use MasterDark as a fixed offset value
            try:
                mdarkdata = float(master_dark)
            except ValueError:
                warnings.warn('MasterDark cannot be loaded, SKIP the dark calibration', RuntimeWarning)
                skip_dark = True
            else:
                msg = 'MasterDark used as a fixed offset value.\n'\
                     'Its *STRONGLY* recommended to use a proper MasterDark'
                warnings.warn(msg, RuntimeWarning)
                skip_dark = False
        else:
            texpdark = float(mdark_hdul[0].header['EXPOSURE'])
            if texpdark == self.fits_Texp:
                mdarkdata = mdark_hdul[0].data
            else:
                if master_bias is None:
                    msg = "Science and Dark don't have the same exposure ({} vs {})!".format(texpdark, self.fits_Texp)
                    warnings.warn(msg, RuntimeWarning)
                    mdarkdata = mdark_hdul[0].data
                else:
                    print('Creating synthetic Dark ...'),
                    mbias_hdul = self.load_bias(master_bias)
                    mdarkdata = mdark_hdul[0].data - mbias_hdul[0].data
                    mdarkdata *= (self.fits_Texp / texpdark)
                    mdarkdata += mbias_hdul[0].data
                    print('OK')

            skip_dark = False

        print('Calibrating image with MasterFlat and MasterDark ...'),

        # Subtract dark frame
        if not skip_dark:
            self.fits_data = self.dark_correction(self.fits_data, mdarkdata)

        # Subtract background / bias (measure it in the non-illuminated corners
        # of the image).
        if self.subtract_corners_background and image_info is not None:
            pass
        else:
            print ('Other correction')
            image_coordinates = ImageCoordinates(image_info)
            data_corners = self.fits_data[image_coordinates.altitude_map < -20]
            self.bias_image_median = np.median(data_corners)
            self.bias_image_std = np.std(data_corners)
            self.bias_image_err = self.bias_image_std / \
                np.sqrt(np.size(data_corners))
            self.fits_data = self.fits_data - self.bias_image_median
            print("Removed: %.2f +/- %.2f counts from measured background"
                  % (self.bias_image_median, self.bias_image_err))

            if hasattr(image_info, 'summary_path') and \
                            image_info.summary_path not in [False, "False", "false", "F", "screen"]:
                if not os.path.exists(image_info.summary_path):
                    os.makedirs(image_info.summary_path)
                measured_bias_log = open(
                    image_info.summary_path + '/measured_image_bias.txt', 'a+')
                text_to_log = str(image_info.date_string) + ',' + str(image_info.used_filter) + ',' +\
                    str(self.bias_image_median) + ',' + \
                    str(self.bias_image_err) + '\r\n'
                measured_bias_log.write(text_to_log)
                measured_bias_log.close()

        # Flat field correction
        if not skip_flat:
            self.fits_data = self.flat_field_correction(self.fits_data, mflatdata)

        print('OK')

    def dark_correction(self, data, masterdark_data):
        return data - masterdark_data

    def flat_field_correction(self, data, masterflat_data):
        return data / masterflat_data

    def flip_image(self):

        try:
            self.fits_data_notcalibrated = np.fliplr(
                self.fits_data_notcalibrated)
        except StandardError:
            print('Warning. Cannot flip raw image as requested')

        try:
            self.fits_data = np.fliplr(self.fits_data)
        except StandardError:
            print('Warning. Cannot flip calibrated image as requested')

    def __clear__(self):
        backup_attributes = [
            "fits_data", "fits_Header", "fits_data_notcalibrated"]

        for atribute in list(self.__dict__):
            # if atribute[0]!="_" and atribute not in backup_attributes:
            if atribute not in backup_attributes:
                del vars(self)[atribute]
