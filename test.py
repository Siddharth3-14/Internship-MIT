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
    """Calc_l

    Find the seperation between two points on the celestial sphere located at points (ra1,dec1) and (ra2,dec2)
    returns the value in armin

    """
    c1 = SkyCoord(ra1,dec1,unit = 'deg')
    c2 = SkyCoord(ra2,dec2,unit = 'deg')
    sep = c1.separation(c2)
    return sep.arcminute

########## importing and testing the file
FITS1 = '../FITS_file/OMC_BandE.fits'
hdul = fits.open(FITS1)
# print(hdul.info())
MapStokesI = hdul[0]
MapStokesQ = hdul[2]
MapStokesU = hdul[4]
MapDebPol = hdul[8]
MapPolAngle = hdul[11]
MapPolFlux = hdul[13]
MapPolFluxError = hdul[14]


MapPolSNR = MapPolFlux.copy()
BlankedMapPol = MapDebPol.copy()
BlankedMapPolAngle = MapPolAngle.copy()
BlankedMapStokesI = MapStokesI.copy()
BlankedMapStokesQ = MapStokesQ.copy()
BlankedMapStokesU = MapStokesU.copy()

######## taking points only with singal to noise ratio more than 2
MapPolSNR.data[:] = np.nan
MapPolSNR.data = MapPolFlux.data/MapPolFluxError.data
Selector = (MapPolSNR.data < 2)

BlankedMapPol.data[Selector] = np.nan
BlankedMapPolAngle.data[Selector] = np.nan
BlankedMapStokesI.data[Selector] = np.nan
BlankedMapStokesQ.data[Selector] = np.nan
BlankedMapStokesU.data[Selector] = np.nan

############## generating the RA and DEC mesh
DEC_grid,RA_grid = generate_RA_DEC_mesh(hdul[0])
seperation = MapPolAngle.copy()


############## Testing the algorithm at point x_index,y_index

x_index = 50
y_index = 50

############## making the filter for selecting points withing the ring
seperation.data = Calc_l(RA_grid[x_index,y_index],DEC_grid[x_index,y_index],RA_grid,DEC_grid)
seperation_selector = (seperation.data<0.5*0.5)
seperation.data[seperation_selector] = np.nan
seperation_selector = (seperation.data>1.5*0.5)
seperation.data[seperation_selector] = np.nan


############## first version
AngleDiff = BlankedMapPolAngle.data - BlankedMapPolAngle.data[x_index,y_index]
Angle_selector =AngleDiff>90
AngleDiff[Angle_selector] = AngleDiff[Angle_selector] - 180
Angle_selector = AngleDiff<-90
AngleDiff[Angle_selector] = AngleDiff[Angle_selector] + 180
seperation_selector = (seperation.data >0)
S = np.nanmean(AngleDiff[seperation_selector]**2)**0.5
print(S)


############ second version
tempa = BlankedMapStokesQ.data*BlankedMapStokesU.data[x_index,y_index] - BlankedMapStokesQ.data[x_index,y_index]*BlankedMapStokesU.data
tempb = BlankedMapStokesQ.data*BlankedMapStokesQ.data[x_index,y_index] + BlankedMapStokesU.data*BlankedMapStokesU.data[x_index,y_index]

AngleDiff_v2 = (180/np.pi)*0.5*np.arctan(tempa/tempb)
Angle_selector_v2 =AngleDiff_v2>90
AngleDiff_v2[Angle_selector_v2] = AngleDiff_v2[Angle_selector_v2] - 180
Angle_selector_v2 = AngleDiff_v2<-90
AngleDiff_v2[Angle_selector_v2] = AngleDiff_v2[Angle_selector_v2] + 180
S2 = np.nanmean(AngleDiff_v2[seperation_selector]**2)**0.5
print(S2)


########## Running the algorithm on each cell
set_delta = 0.5   # in arcminute
S_map = BlankedMapPolAngle.copy()
S_map_v2 = BlankedMapPolAngle.copy()
for i in range(RA_grid.shape[0]):
    for j in range(RA_grid.shape[1]):

        ##### seperation filter
        seperation.data = Calc_l(RA_grid[i,j],DEC_grid[i,j],RA_grid,DEC_grid)
        seperation_selector = (seperation.data<0.5*set_delta)
        seperation.data[seperation_selector] = np.nan
        seperation_selector = (seperation.data>1.5*set_delta)
        seperation.data[seperation_selector] = np.nan
        seperation_selector = (seperation.data >0)

        
        ##### first version
        AngleDiff = BlankedMapPolAngle.data - BlankedMapPolAngle.data[i,j]
        Angle_selector =AngleDiff>90
        AngleDiff[Angle_selector] = AngleDiff[Angle_selector] - 180
        Angle_selector = AngleDiff<-90
        AngleDiff[Angle_selector] = AngleDiff[Angle_selector] + 180
        
        ## once more to take care of > 180 degree angle difference
        Angle_selector =AngleDiff>90
        AngleDiff[Angle_selector] = AngleDiff[Angle_selector] - 180
        Angle_selector = AngleDiff<-90
        AngleDiff[Angle_selector] = AngleDiff[Angle_selector] + 180
        
        S = np.nanmean(AngleDiff[seperation_selector]**2)**0.5
        S_map.data[i,j] = S

        
        ##### second version
        tempa = BlankedMapStokesQ.data*BlankedMapStokesU.data[i,j] - BlankedMapStokesQ.data[i,j]*BlankedMapStokesU.data
        tempb = BlankedMapStokesQ.data*BlankedMapStokesQ.data[i,j] + BlankedMapStokesU.data*BlankedMapStokesU.data[i,j]
        AngleDiff_v2 = 0.5 * (180/np.pi)*np.arctan2(tempa,tempb)
        #Angle_selector_v2 =AngleDiff_v2>90
        #AngleDiff_v2[Angle_selector_v2] = AngleDiff_v2[Angle_selector_v2] - 180
        #Angle_selector_v2 = AngleDiff_v2<-90
        #AngleDiff_v2[Angle_selector_v2] = AngleDiff_v2[Angle_selector_v2] + 180
        S_v2 = np.nanmean(AngleDiff_v2[seperation_selector]**2)**0.5
        S_map_v2.data[i,j] = S_v2

plt.figure()
ax1 = plt.subplot(131)
ax1.imshow(S_map.data,origin='lower')
ax2 = plt.subplot(132)
ax2.imshow(S_map_v2.data,origin='lower')
ax3 = plt.subplot(133)
ax3.imshow(S_map.data-S_map_v2.data,origin='lower')
ax1.set_title('S map old')
ax2.set_title('S map new')
ax3.set_title('S map old - S map new')
plt.show()

adf_v1 = S_map.data.flatten()
adf_v2 = S_map_v2.data.flatten()

plt.figure()
ax1 = plt.subplot(121)
ax1.scatter(adf_v1,adf_v2)
ax2 = plt.subplot(122)
ax2.plot(adf_v1-adf_v2)
ax1.set_title('S map old vs S map new')
ax2.set_title('S map old - S map new')
plt.show()
# s_map_array = S_map.data.flatten()
# Pol_array = BlankedMapPol.data.flatten()

# y_min = np.nanmin(Pol_array)
# y_max = 20
  
# x_min = np.nanmin(s_map_array)
# x_max = np.nanmax(s_map_array)
  
# x_bins = np.arange(x_min, x_max, 0.5)
# y_bins = np.arange(y_min, y_max, 0.5)

# fig = plt.subplots(figsize =(10, 10))
# # Creating plot
# plt.hist2d(s_map_array,Pol_array,bins =[x_bins, y_bins])
# plt.title("P X S 2D histogram")
# plt.ylabel('P [%]')
# plt.xlabel('S [deg]')
# plt.show()


# I_map_array = BlankedMapStokesI.data.flatten()

# y_min = np.nanmin(Pol_array)
# y_max = 20
  
# x_min = np.nanmin(I_map_array)
# x_max = np.nanmax(I_map_array)
  
# x_bins = np.arange(x_min, x_max, 1)
# y_bins = np.arange(y_min, y_max, 0.25)

# fig = plt.subplots(figsize =(10, 10))
# # Creating plot
# plt.hist2d(I_map_array,Pol_array,bins =[x_bins, y_bins])
# plt.title("P X I 2D histogram")
# plt.ylabel('P [%]')
# plt.xlabel('I ')
# plt.show()



# ################### 1st attempt at plotting  S vs P ###################################
# # Pol = hdul[7].data
# # s_map = s_map + Pol*(Pol == np.nan)
# # Pol_mod = Pol*(Pol<20)
# # plt.figure()
# # ax1 = plt.subplot(121)
# # ax1.set_title('Pol [%]')
# # ax1.set_xlabel('RA')
# # ax1.set_ylabel('DEC')
# # ax1.imshow(Pol_mod,origin='lower')
# # ax2 = plt.subplot(122)
# # ax2.imshow(s_map,origin='lower')
# # ax2.set_title('S [deg]')
# # ax2.set_xlabel('RA')
# # ax2.set_ylabel('DEC')
# # plt.show()

# # s_map_array = s_map.flatten()
# # Pol_array = Pol_mod.flatten()
# # new_s_map_array = [x for x in s_map_array if np.isnan(x) == False]
# # new_Pol_array = [x for x in Pol_array if np.isnan(x) == False]

# # x_min = np.min(new_Pol_array)
# # x_max = np.max(new_Pol_array)
  
# # y_min = np.min(new_s_map_array)
# # y_max = np.max(new_s_map_array)
  
# # x_bins = np.linspace(x_min, x_max, 50)
# # y_bins = np.linspace(y_min, y_max, 50)

# # fig = plt.subplots(figsize =(10, 10))
# # # Creating plot
# # plt.hist2d(new_Pol_array, new_s_map_array,bins =[x_bins, y_bins])
# # plt.title("S X P 2D histogram")
# # plt.ylabel('S [deg]')
# # plt.xlabel('P [%]')
  
# # # show plot
# # plt.show()