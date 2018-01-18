

"""Basic image reduction"""


import ccdproc
from astropy import units as u


def reduction(raw_name, dark_name, flat_name):
    """Perform basic reduction"""
    raw_image = ccdproc.CCDData.read(raw_name, unit='adu')
    dark_image = ccdproc.CCDData.read(dark_name, unit='adu')
    flat_image = ccdproc.CCDData.read(flat_name, unit='adu')

    dark_sub = ccdproc.subtract_dark(raw_image, dark_image,
                                     exposure_time='EXPOSURE',
                                     exposure_unit=u.second)

    flat_corr = ccdproc.flat_correct(dark_sub, flat_image)
    return flat_corr
