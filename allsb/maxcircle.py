
import logging

import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches


import skimage.morphology as M
import skimage.filters  as F
from skimage.measure import label, regionprops


_logger = logging.getLogger(__name__)


def maxcircle(hdulist, do_plots=False):
    return maxcircle_img(hdulist[0].data, do_plots)


def maxcircle_img(image, do_plots=False):
    _logger.debug('maxcircle with %s', image.dtype)

    # sobel works because the border is bright

    markers = np.zeros_like(image)
    markers[0:10,0:10] = 1
    markers[2000:2010, 2000:2010] = 2
    #print('T mim', F.threshold_minimum(image))

    th_o = F.threshold_otsu(image)
    _logger.debug('reference threshold, Otsu method %s', th_o)

    #markers[image < th_o *0.9] = 1
    #markers[image > th_o *1.1] = 2

    #plt.hist(image.ravel()[::100], bins='auto')
    #plt.show()

    emap = F.sobel(image)
    _logger.debug('sobel done')
    if do_plots:
        plt.imshow(emap)
        plt.show()
    _logger.debug('segmentation')
    segmentation = M.watershed(emap, markers)
    _logger.debug('segmentation done')
    if do_plots:
        plt.imshow(segmentation)
        plt.show()

    _logger.debug('label')
    segmentation2 = M.dilation(segmentation)
    if do_plots:
        plt.imshow(segmentation2)
        plt.show()
    label_img, nlabels = label(segmentation2, background=1, return_num=True)

    _logger.debug('labeling done')
    _logger.debug('properties')
    regions = regionprops(label_img, coordinates='xy', cache=True)
    _logger.debug('properties done')
    _logger.debug('nregions %d', len(regions))
    fig, ax = plt.subplots()
    ax.imshow(image, cmap=plt.cm.gray)

    x0 = 0
    y0 = 0
    equivalent_diameter = 0
    for props in regions:
        y0, x0 = props.centroid
        equivalent_diameter = props.equivalent_diameter

        e1 = patches.Ellipse((x0, y0), props.equivalent_diameter, props.equivalent_diameter,
                             angle=0, edgecolor='b', linewidth=2, fill=False, zorder=2)

        ax.add_patch(e1)

    if do_plots:
        plt.show()

    return x0, y0, equivalent_diameter


def main(args=None):
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parsed_args = parser.parse_args(args)

    with fits.open(parsed_args.filename) as hdulist:
        maxcircle(hdulist)


if __name__ == '__main__':

    main()
