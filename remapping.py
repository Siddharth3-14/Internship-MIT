from os import sep
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from scipy.optimize import curve_fit
# from PolarSi import *


def generate_RA_DEC_mesh(hdr):
    """generate_RA_DEC_mesh

    Generates the RA and DEC grid for the intensity map

    """
    if 'CDELT1' in hdr.header:
        RA_delt = hdr.header['CDELT1']
        DEC_delt = hdr.header['CDELT2']

    if 'CD1_1' in hdr.header:
        RA_delt = hdr.header['CD1_1']
        DEC_delt = hdr.header['CD2_2']

    RA_ref = (hdr.header['CRPIX1'])
    DEC_ref = (hdr.header['CRPIX2'])
    RA_ref_value = hdr.header['CRVAL1']
    DEC_ref_value = hdr.header['CRVAL2']
    RA_axis_len = hdr.header['NAXIS1']
    DEC_axis_len = hdr.header['NAXIS2']

    RA_axis = np.arange(1,RA_axis_len+1)
    DEC_axis = np.arange(1,DEC_axis_len+1)
    DEC_axis_modified = np.arange(1,RA_axis_len+1)
    
    DEC_array = (DEC_axis - DEC_axis_len/2)*DEC_delt + DEC_ref_value
    DEC_array_modified = (DEC_axis_modified - RA_axis_len/2)*DEC_delt + DEC_ref_value
    RA_array = RA_ref_value-(RA_axis - RA_axis_len/2)*(RA_delt*(-1)/np.cos(DEC_array_modified*0.01745))

    # #making a meshgrid from the arrays
    DEC_grid,RA_grid = np.meshgrid(DEC_array,RA_array , sparse=False, indexing='ij')
    return DEC_grid,RA_grid

FITS1 = '../FITS_file/CygX_E_OTFMAP.fits'
hdul = fits.open(FITS1)
# print(hdul.info())
MapPolAngle = hdul[11]
BlankedMapPolAngle = MapPolAngle.copy()
############## generating the RA and DEC mesh
DEC_grid,RA_grid = generate_RA_DEC_mesh(hdul[0])
PA_grid = MapPolAngle.copy()
df = pd.read_csv('DR21_CSO.dat',delimiter='\s+',header = None,low_memory=False)
# df.head()

RA_array = df[6]
DEC_array = df[7]
PA_array = df[17]+360
ePA_array = df[18]
def mapping_func(grid_xx,grid_yy,x,y,f_x_y,tol):
    """mapping_func

    This function maps f(x,y) onto (x,y) grid from a data array of x , y and f(x,y)
    
    Args:
        grid_xx (2D array): x coordinates mesh 
        grid-yy (2D array): y coordinates mesh 
        x (1D array): x coordinate value at which f is evaluated
        y (1D array): y coordinate value at which f is evaluated
        f_x_y (1D array): value of f evaluated at x and y
    Returns:
        mappped_func (2D array): f(x,y) mapped into (x,y)
    """
    mappped_func = np.zeros((grid_xx.shape[0],grid_xx.shape[1]))
    for i in range(f_x_y.shape[0]):
        temp = (abs(grid_xx - x[i])<tol/2)*(abs(grid_yy - y[i])<tol/2)*f_x_y[i]
        mappped_func += temp*(mappped_func==0)
        # plt.imshow(temp)
        # plt.show()
    return mappped_func

tol =  DEC_grid[1,0] -  DEC_grid[0,0] 
print(tol)
optical_mapped_pol= mapping_func(RA_grid,DEC_grid,RA_array,DEC_array,PA_array,tol)
PA_grid.data = optical_mapped_pol
selector = (optical_mapped_pol==0) 
PA_grid.data[selector] = np.nan
plt.figure()
plt.imshow(PA_grid.data)
plt.show() 
PA_grid.data = PA_grid.data - 360
plt.figure()
plt.imshow(PA_grid.data)
plt.show()
PA_grid.writeto('remapped_pa.fits')