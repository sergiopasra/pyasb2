
import argparse
import logging
import astropy.io.fits as fits

# FIXME: UHMMM
from allsb.testastro import wcs_calibrate_astrometry_net

_logger = logging.getLogger(__name__)


ASTROMETRY_API_KEY = 'gmbziuhfdkfnpysg'


def main(args=None):

    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parsed_args = parser.parse_args(args)

    logging.basicConfig(level=logging.DEBUG)

    _logger.debug('filename is: %s', parsed_args.filename)

    datafiles = basic_process_files(parsed_args.filename)

    for datafile in datafiles:
        if not datafile['wcs']:
            wcs_calibrate_astrometry_net(datafile)


def basic_process_files(filename):

    with fits.open(filename) as hdulist:
        # get metadata
        # get data
        data = hdulist[0].data
        return [{'data': data, 'wcs': False,
                 'shape': data.shape, 'filename': filename}]


def plot_func(data):
    import matplotlib.pyplot as plt

    plt.imshow(data)
    plt.show()


def wcs_calibrate(datafile):
    _logger.debug('shape is {shape}'.format(**datafile))
    plot_func(datafile['data'])


if __name__ == '__main__':

    main()
