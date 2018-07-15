
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

#from wcsaxes import WCS

latitude = 40.450941
refx, refy = 1231.36087524, 1218.93356496
cdeltx, cdelty = 0.0704340939421, 0.0704340939421
ref_alt, ref_az = 89.180358574, 141.230103929

az_zero0 = 90 - 88.64589921

header = fits.Header()

w1 = WCS(header=header, naxis=2)
w1.wcs.crpix = [refx, refy]
w1.wcs.cdelt = [cdeltx, cdelty]
w1.wcs.crval = [90, 0]
w1.wcs.ctype = ["pLAT-ZEA", "pLON-ZEA"]
w1.wcs.lonpole = 90

print(w1)

w1.all_world2pix([[0,0]], 1)
header1 = w1.to_header()

w2 = WCS(naxis=2)
w2.wcs.crpix = [refx, refy]
w2.wcs.cdelt = [cdeltx, cdelty]
w2.wcs.crval = [ref_alt, ref_az]
w2.wcs.ctype = ["pLAT-ZEA", "pLON-ZEA"]
w2.wcs.lonpole = ref_az + 90 + az_zero0

header2 = w2.to_header(key='A')

w3 = WCS(naxis=2)
w3.wcs.crpix = [refx, refy]
w3.wcs.cdelt = [cdeltx, cdelty]
w3.wcs.crval = [latitude, 0]
w3.wcs.ctype = ["DEC--ZEA", "RA---ZEA"]
w3.wcs.lonpole = az_zero0 - 90

header3 = w3.to_header(key='B')


print ref_alt, ref_az

header.update(header1)
header.update(header2)
header.update(header3)

data0 = np.zeros((2500, 2500))

data = fits.getdata('Johnson_V20130912_011709.fit.gz')

hdu = fits.PrimaryHDU(data, header=header)
hdul = fits.HDUList([hdu])

hdul.writeto('test.fits', clobber=True)

fig = plt.figure()

ax1 = fig.add_axes([0.0, 0.0, 1.0, 1.0], projection=w1)

ax1.coords.grid(color='blue', alpha=1, linestyle='solid')
overlay1 = ax1.get_coords_overlay(w1)
overlay1.grid(color='white', linestyle='solid', alpha=1)
#overlay3 = ax1.get_coords_overlay(w3)
#overlay3.grid(color='green', linestyle='solid', alpha=1)
ax1.imshow(data, origin='lower')

plt.show()
