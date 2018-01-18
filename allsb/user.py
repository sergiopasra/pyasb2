

"""Command line interface"""

import argparse
import configparser
import logging

import astropy.io.fits as fits

# import datetime

from allsb.reduction import reduction

def main():

    config = configparser.ConfigParser()

    parser = argparse.ArgumentParser(description='Process ASTMON images.')
    parser.add_argument('-c', '--config',
                        help='Configuration file',
                        default='config.ini')
    parser.add_argument('images', nargs='+', help='ASTMON image')

    args = parser.parse_args()
    config.read(args.config)

    for filename in args.images:
        print('filename:', filename)

        with fits.open(filename) as hdul:
            header = hdul[0].header
            # date_str = header['DATE']
            # image_dt = datetime.datetime.strptime(date_str, "%Y%m%d_%H%M%S")
            image_filter = header['FILTER']
            # calibrations

        calib_d = config['calibrations_{}'.format(image_filter)]
        fc = reduction(filename, calib_d['darkframe'], calib_d['flatfield'])
        fc.write('test.fits', overwrite=True)


if __name__ == '__main__':

    main()