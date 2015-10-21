
# FITS operations
#
# Calculate the sum, difference, multiply or divide fits images.
# The header is taken from the first image, with added comment.
# ____________________________
#
# This module is part of the PyASB project,
# created and maintained by Miguel Nievas [UCM].
# ____________________________


import inspect

import numpy as np
import astropy.io.fits as pyfits
import datetime


def verbose(function, *args):
    '''
    Run a function in verbose mode
    '''
    try:
        out = function(*args)
    except StandardError:
        # Something happened while runing function
        raise
    else:
        return(out)


def load_image(image):
    print('Loading Data and Header for the given Image ...'),
    Image_HDU = pyfits.open(image)
    Image_Data = Image_HDU[0].data
    Image_Header = Image_HDU[0].header
    # Convert to string
    #Image_Header_text = Image_Header.tostring()
    #Image_Header_text = encode_utf8_to_iso88591(Image_Header_text)
    # Image_Header.fromstring(Image_Header_text)
    # print(Image_Header)
    print('OK')
    return(Image_Data, Image_Header)


class FitsOperation(object):

    def __init__(self, fits1, fits2):
        self.loadfits(fits1, fits2)

    @staticmethod
    def get_datetime_filename(self):
        return(str(datetime.datetime.now()).replace(" ", "_").replace(":", "-").split(".")[0])

    def loadfits(self, fits1, fits2):
        self.Data1, self.Header1 = \
            verbose(load_image, fits1)
        self.Data2, self.Header2 = \
            verbose(load_image, fits2)

    def sumfits(self):
        self.DataResult = self.Data1 + self.Data2
        self.HeaderResult = self.Header1

    def subtractfits(self):
        self.DataResult = self.Data1 - self.Data2
        self.HeaderResult = self.Header1

    def multiplyfits(self):
        self.DataResult = self.Data1 * self.Data2
        self.HeaderResult = self.Header1

    def dividefits(self):
        self.DataResult = self.Data1 * 1. / self.Data2
        self.HeaderResult = self.Header1

    def normalizefits(self):
        self.DataResult = self.DataResult * 1. / np.median(self.DataResult)

    def addnote(self, note="Fits edited with PyASB.fits_operator.py"):
        self.HeaderResult.add_comment(note)

    def savefits(self, filename=None):
        if filename == None:
            filename = self.get_datetime_filename
        pyfits.writeto(
            filename, self.DataResult, self.HeaderResult, clobber=True)
