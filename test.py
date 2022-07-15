Internship-MIT
from os import sep
import numpy as np
import matplotlib.pyplot as plt
# import aplpy
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord


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

    RA_ref = int(hdr.header['CRPIX1'])
    DEC_ref = int(hdr.header['CRPIX2'])
    RA_ref_value = hdr.header['CRVAL1']
    DEC_ref_value = hdr.header['CRVAL2']
    RA_axis_len = hdr.header['NAXIS1']
    DEC_axis_len = hdr.header['NAXIS2']

    RA_axis = np.arange(0,RA_axis_len,1)
    DEC_axis = np.arange(0,DEC_axis_len,1)
    RA_axis = RA_ref_value - RA_delt*(RA_ref - RA_axis)
    DEC_axis = DEC_ref_value + DEC_delt*(DEC_axis - DEC_ref)
    DEC_grid,RA_grid = np.meshgrid(DEC_axis,RA_axis , sparse=False, indexing='ij')


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

def Calc_l(ra1,dec1,ra2,dec2):

    c1 = SkyCoord(ra1,dec1,unit = 'deg')
    c2 = SkyCoord(ra2,dec2,unit = 'deg')
    sep = c1.separation(c2)
    return sep.arcminute

FITS1 = 'FITS_file/OMC_BandE.fits'
hdul = fits.open(FITS1)
# print(hdul.info())
Pol_angle = hdul[11].data

# # evt_data = Table(hdul[18].data)
# # print(evt_data.colnames)
DEC_grid,RA_grid = generate_RA_DEC_mesh(hdul[0])
# # print(1.5*0.5*(DEC_grid[1,0] - DEC_grid[0,0])*60)

# # plt.figure()
# # ax1 = plt.subplot(121)
# # ax1.imshow(DEC_grid)
# # ax2 = plt.subplot(122)
# # ax2.imshow(RA_grid)
# # plt.show()

s_map = np.zeros_like(Pol_angle)
seperation = Calc_l(RA_grid[50,50],DEC_grid[50,50],RA_grid,DEC_grid)
seperation = seperation*(seperation<1.5*0.5)*(0.5*0.5<seperation)
plt.imshow(seperation)
plt.show()
Pol_mask = Pol_angle*(seperation>0)
plt.imshow(Pol_mask)
plt.show()
Pol_diff =  ((Pol_angle-Pol_angle[50,50])**2)*(seperation>0)
plt.imshow(Pol_diff)
plt.show()
S = np.sqrt(np.nansum(Pol_diff)/np.count_nonzero(Pol_diff))
print(S)
s_map[50,50] = S


s_map = np.zeros_like(Pol_angle)
for i in range(RA_grid.shape[0]):
    for j in range(RA_grid.shape[1]):
        seperation = Calc_l(RA_grid[i,j],DEC_grid[i,j],RA_grid,DEC_grid)
        seperation = seperation*(seperation<1.5*0.5)*(0.5*0.5<seperation)
        # plt.imshow(seperation)
        # plt.show()
        Pol_mask = Pol_angle*(seperation>0)
        # plt.imshow(Pol_mask)
        # plt.show()
        Pol_diff =  ((Pol_angle- Pol_angle[i,j])**2)*(seperation>0)
        # plt.imshow(Pol_diff)
        # plt.show()
        S = np.sqrt(np.nansum(Pol_diff)/np.count_nonzero(Pol_diff))
        s_map[i,j] = S



################### 1st attempt at plotting  S vs P ###################################
# Pol = hdul[7].data
# s_map = s_map + Pol*(Pol == np.nan)
# Pol_mod = Pol*(Pol<20)
# plt.figure()
# ax1 = plt.subplot(121)
# ax1.set_title('Pol [%]')
# ax1.set_xlabel('RA')
# ax1.set_ylabel('DEC')
# ax1.imshow(Pol_mod,origin='lower')
# ax2 = plt.subplot(122)
# ax2.imshow(s_map,origin='lower')
# ax2.set_title('S [deg]')
# ax2.set_xlabel('RA')
# ax2.set_ylabel('DEC')
# plt.show()

# s_map_array = s_map.flatten()
# Pol_array = Pol_mod.flatten()
# new_s_map_array = [x for x in s_map_array if np.isnan(x) == False]
# new_Pol_array = [x for x in Pol_array if np.isnan(x) == False]

# x_min = np.min(new_Pol_array)
# x_max = np.max(new_Pol_array)
  
# y_min = np.min(new_s_map_array)
# y_max = np.max(new_s_map_array)
  
# x_bins = np.linspace(x_min, x_max, 50)
# y_bins = np.linspace(y_min, y_max, 50)

# fig = plt.subplots(figsize =(10, 10))
# # Creating plot
# plt.hist2d(new_Pol_array, new_s_map_array,bins =[x_bins, y_bins])
# plt.title("S X P 2D histogram")
# plt.ylabel('S [deg]')
# plt.xlabel('P [%]')
  
# # show plot
# plt.show()