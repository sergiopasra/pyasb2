
import math

import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches


import skimage.morphology as M
import skimage.filters  as F
from skimage.measure import label, regionprops



def maxcircle(hdulist):
    maxcircle_img(hdulist[0].data)


def maxcircle_img(image):
    print('maxcircle', image.dtype)

    # sobel works because the border is bright

    markers = np.zeros_like(image)
    markers[0:10,0:10] = 1
    markers[2000:2010, 2000:2010] = 2
    #print('T mim', F.threshold_minimum(image))

    th_o = F.threshold_otsu(image)
    print('T ot', th_o)

    #markers[image < th_o *0.9] = 1
    #markers[image > th_o *1.1] = 2

    #plt.hist(image.ravel()[::100], bins='auto')
    #plt.show()

    emap = F.sobel(image)
    print('sobel done')
    plt.imshow(emap)
    plt.show()
    print('segmentation')
    segmentation = M.watershed(emap, markers)
    print('segmentation done')
    plt.imshow(segmentation)
    plt.show()
    print(segmentation.max())
    print('label')
    segmentation2 = M.dilation(segmentation)
    plt.imshow(segmentation2)
    plt.show()
    label_img, nlabels = label(segmentation2, background=1, return_num=True)
    print(nlabels)

    print('labeling done')
    print('properties')
    regions = regionprops(label_img, coordinates='xy', cache=True)
    print('properties done')
    print('nregions', len(regions))
    fig, ax = plt.subplots()
    ax.imshow(image, cmap=plt.cm.gray)

    for props in regions:
        y0, x0 = props.centroid
        print(y0, x0)
        print(props.equivalent_diameter)

        e1 = patches.Ellipse((x0, y0), props.equivalent_diameter, props.equivalent_diameter,
                             angle=0, edgecolor='b', linewidth=2, fill=False, zorder=2)

        ax.add_patch(e1)

    plt.show()


    return


def main(args=None):
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parsed_args = parser.parse_args(args)

    with fits.open(parsed_args.filename) as hdulist:
        maxcircle(hdulist)


if __name__ == '__main__':

    main()



