from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

header = fits.getheader('..\FITS_file\CygX_E_OTFMAP.fits')
w = WCS(header)
px = np.linspace(200., 300., 10)
py = np.linspace(200., 300., 10)
wx, wy = w.wcs_pix2world(px, py, 1)
print(wx.shape)