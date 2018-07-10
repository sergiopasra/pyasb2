
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np
import sep
from astropy.wcs import WCS

latitude = 40.450941
refx, refy = 1231.36087524, 1218.93356496
cdeltx, cdelty = 0.0704340939421, 0.0704340939421
ref_alt, ref_az = 89.180358574, 141.230103929
az_zero0 = 90 - 88.64589921


# data = fits.getdata('test.fits')[600:1900,600:1900]
data = fits.getdata('test.fits')

ii, jj = np.mgrid[0:2500, 0:2500]

ss = np.hypot(ii - 1220, jj - 1230)
rr1 = ss > 1142.84
rr2 = ss > (1142.84 - 200)
rr3 = ss > (1142.84 - 300)

fits.writeto('mask1.fits', rr1.astype('int'), clobber=True)
fits.writeto('mask2.fits', rr2.astype('int'), clobber=True)

# m, s = np.mean(data), np.std(data)
# plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-s, vmax=m+s, origin='lower')
# plt.colorbar()
# plt.show()

# background
data_c = data.copy()
data_c = data_c.byteswap().newbyteorder()
bkg = sep.Background(data_c, mask=rr2)
bkg_image = bkg.back()

# plt.imshow(bkg_image, interpolation='nearest', cmap='gray', origin='lower')
# plt.colorbar()
# plt.show()

data_sub = data_c - bkg

objects = sep.extract(data_sub, 1.5, err=bkg.globalrms, mask=rr3)

print(len(objects))

# plot background-subtracted image
w2 = WCS(naxis=2)
w2.wcs.crpix = [refx, refy]
w2.wcs.cdelt = [cdeltx, cdelty]
w2.wcs.crval = [ref_alt, ref_az]
w2.wcs.ctype = ["pLAT-ZEA", "pLON-ZEA"]
w2.wcs.lonpole = ref_az + 90 + az_zero0


fig = plt.figure()
ax1 = fig.add_axes([0.0, 0.0, 1.0, 1.0], projection=w2)

ax1.coords.grid(color='blue', alpha=1, linestyle='solid')
overlay1 = ax1.get_coords_overlay(w2)
overlay1.grid(color='white', linestyle='solid', alpha=1)

m, s = np.mean(data_sub), np.std(data_sub)
im = ax1.imshow(data_sub, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

from matplotlib.patches import Ellipse

# plot an ellipse for each object
for i in range(len(objects)):
    e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                width=6*objects['a'][i],
                height=6*objects['b'][i],
                angle=objects['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax1.add_artist(e)

plt.show()
