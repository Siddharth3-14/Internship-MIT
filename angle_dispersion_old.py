from os import sep
import numpy as np
import matplotlib.pyplot as plt
import aplpy
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from scipy.optimize import curve_fit


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

    c1 = SkyCoord(ra1,dec1,unit = 'deg')
    c2 = SkyCoord(ra2,dec2,unit = 'deg')
    sep = c1.separation(c2)
    return sep.arcminute
    
def wrapper(Angle_grid):
    while ((np.nanmax(Angle_grid)>90) or (np.nanmin(Angle_grid)<-90)):
        Angle_selector =Angle_grid>90
        Angle_grid[Angle_selector] = Angle_grid[Angle_selector] - 180
        Angle_selector = Angle_grid<-90
        Angle_grid[Angle_selector] = Angle_grid[Angle_selector] + 180
    return Angle_grid

########## importing and testing the file
########## importing and testing the file
FITS1 = '../FITS_file/CygX_E_OTFMAP.fits'
hdul = fits.open(FITS1)
# print(hdul.info())
MapStokesI = hdul[0]
MapStokesQ = hdul[2]
MapStokesU = hdul[4]
MapDebPol = hdul[8]
MapDebPolError = hdul[9]
MapPolAngle = hdul[11]
MapPolFlux = hdul[13]
MapPolFluxError = hdul[14]


MapPolSNR = MapDebPol.copy()
BlankedMapPol = MapDebPol.copy()
BlankedMapPolAngle = MapPolAngle.copy()
BlankedMapStokesI = MapStokesI.copy()
BlankedMapStokesQ = MapStokesQ.copy()
BlankedMapStokesU = MapStokesU.copy()

######## taking points only with singal to noise ratio more than 2
MapPolSNR.data[:] = np.nan
MapPolSNR.data = MapDebPol.data/MapDebPolError.data
Selector = (MapPolSNR.data < 3)

BlankedMapPol.data[Selector] = np.nan
BlankedMapPolAngle.data[Selector] = np.nan
BlankedMapStokesI.data[Selector] = np.nan
BlankedMapStokesQ.data[Selector] = np.nan
BlankedMapStokesU.data[Selector] = np.nan

############## generating the RA and DEC mesh
DEC_grid,RA_grid = generate_RA_DEC_mesh(hdul[0])
seperation = MapPolAngle.copy()

# plt.figure()
# ax1 = plt.subplot(121)
# ax2 = plt.subplot(122)
# ax2.imshow(BlankedMapStokesI.data)
# ax1.imshow(MapStokesI.data)
# plt.show()

# ############## Testing the algorithm at point x_index,y_index

# x_index = 50
# y_index = 50

# ############## making the filter for selecting points withing the ring
# seperation.data = Calc_l(RA_grid[x_index,y_index],DEC_grid[x_index,y_index],RA_grid,DEC_grid)
# seperation_selector = (seperation.data<0.5*0.5)
# seperation.data[seperation_selector] = np.nan
# seperation_selector = (seperation.data>1.5*0.5)
# seperation.data[seperation_selector] = np.nan


# ############## first version
# AngleDiff = BlankedMapPolAngle.data - BlankedMapPolAngle.data[x_index,y_index]
# Angle_selector =AngleDiff>90
# AngleDiff[Angle_selector] = AngleDiff[Angle_selector] - 180
# Angle_selector = AngleDiff<-90
# AngleDiff[Angle_selector] = AngleDiff[Angle_selector] + 180
# seperation_selector = (seperation.data >0)
# S = np.nanmean(AngleDiff[seperation_selector]**2)**0.5
# print(S)


# ############ second version
# tempa = BlankedMapStokesQ.data*BlankedMapStokesU.data[x_index,y_index] - BlankedMapStokesQ.data[x_index,y_index]*BlankedMapStokesU.data
# tempb = BlankedMapStokesQ.data*BlankedMapStokesQ.data[x_index,y_index] + BlankedMapStokesU.data*BlankedMapStokesU.data[x_index,y_index]

# AngleDiff_v2 = (180/np.pi)*0.5*np.arctan(tempa/tempb)
# Angle_selector_v2 =AngleDiff_v2>90
# AngleDiff_v2[Angle_selector_v2] = AngleDiff_v2[Angle_selector_v2] - 180
# Angle_selector_v2 = AngleDiff_v2<-90
# AngleDiff_v2[Angle_selector_v2] = AngleDiff_v2[Angle_selector_v2] + 180
# S2 = np.nanmean(AngleDiff_v2[seperation_selector]**2)**0.5
# print(S2)


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
        # AngleDiff = BlankedMapPolAngle.data - BlankedMapPolAngle.data[i,j]
        # AngleDiff = wrapper(AngleDiff)
        # # Angle_selector =AngleDiff>90
        # # AngleDiff[Angle_selector] = AngleDiff[Angle_selector] - 180
        # # Angle_selector = AngleDiff<-90
        # # AngleDiff[Angle_selector] = AngleDiff[Angle_selector] + 180
        
        # # ## once more to take care of > 180 degree angle difference
        # # Angle_selector =AngleDiff>90
        # # AngleDiff[Angle_selector] = AngleDiff[Angle_selector] - 180
        # # Angle_selector = AngleDiff<-90
        # # AngleDiff[Angle_selector] = AngleDiff[Angle_selector] + 180

        # S = np.nanmean(AngleDiff[seperation_selector]**2)**0.5
        # S_map.data[i,j] = S

        
        # ##### second version
        tempa = BlankedMapStokesQ.data*BlankedMapStokesU.data[i,j] - BlankedMapStokesQ.data[i,j]*BlankedMapStokesU.data
        tempb = BlankedMapStokesQ.data*BlankedMapStokesQ.data[i,j] + BlankedMapStokesU.data*BlankedMapStokesU.data[i,j]
        AngleDiff_v2 = 0.5 * (180/np.pi)*np.arctan2(tempa,tempb)
        #Angle_selector_v2 =AngleDiff_v2>90
        #AngleDiff_v2[Angle_selector_v2] = AngleDiff_v2[Angle_selector_v2] - 180
        #Angle_selector_v2 = AngleDiff_v2<-90
        #AngleDiff_v2[Angle_selector_v2] = AngleDiff_v2[Angle_selector_v2] + 180
        S_v2 = np.nanmean(AngleDiff_v2[seperation_selector]**2)**0.5
        S_map_v2.data[i,j] = S_v2

# S_map_v2.writeto('dispersion.fits')
# BlankedMapStokesI.writeto('intensity.fits')
# BlankedMapPol.writeto('polarization_frac.fits')

# DEC_array = DEC_grid.flatten()
# RA_array = RA_grid.flatten()
# StokesI_array = BlankedMapStokesI.data.flatten()
# S_array = S_map_v2.data.flatten()
# P_array = BlankedMapPol.data.flatten()

def lin_fit(x, a, b):
    return a + b*x

s_array = S_map_v2.data.flatten()
p_array = BlankedMapPol.data.flatten()
I_array = BlankedMapStokesI.data.flatten()

log_s = np.log(s_array)
log_p = np.log(p_array)
log_I = np.log(I_array)


p_min = np.nanmin(log_p)
p_max = np.log(50)
s_min = np.nanmin(log_s)
s_max = np.nanmax(log_s)
I_min = np.nanmin(log_I)
I_max = np.nanmax(log_I)
  

p_bins = np.arange(p_min, p_max, 0.075)
s_bins = np.arange(s_min, s_max, 0.075)
I_bins = np.arange(I_min, I_max, 0.075)

df_log = pd.DataFrame({'logp': log_p,'logs':log_s,'logI':log_I})
# df = df.dropna()
df_log = df_log.dropna()

# PS_param, PS_param_cov = curve_fit(lin_fit, df_log['logs'], df_log['logp'])
# PS_FitFunc = lin_fit(s_bins,PS_param[0],PS_param[1])
# print(PS_param[0],PS_param[1])


# # Plotting log p vs log s 
# fig = plt.subplots(figsize =(15, 10))
# plt.hist2d(log_s,log_p,bins =[s_bins, p_bins])
# label_temp = r'log(p) = C + $\alpha_s$log(S){linebreak} $\alpha_s$: {alpha_s:.4f} C: {C:.04f}'.format(alpha_s = PS_param[1],C = PS_param[0],linebreak='\n')
# plt.plot(s_bins,PS_FitFunc,'r',linewidth=3,label = label_temp)
# plt.title("log p X log S 2D histogram")
# plt.ylabel('log p ')
# plt.xlabel('log S ')
# plt.legend()
# plt.tight_layout()
# plt.show()

PI_param, PI_param_cov = curve_fit(lin_fit, df_log['logI'], df_log['logp'])
PI_FitFunc = lin_fit(I_bins,PI_param[0],PI_param[1])
# print(PI_param[0],PI_param[1])

# Plotting log p vs log I 
fig = plt.subplots(figsize =(10, 10))
plt.hist2d(log_s,log_p,bins =[s_bins, p_bins])
# label_temp = r'log(p) = C + $\alpha_I$log(I){linebreak} $\alpha_I$: {alpha_I:.4f} C: {C:.4f}'.format(alpha_I = PI_param[1],C = PI_param[0],linebreak='\n')
# plt.plot(I_bins,PI_FitFunc,'r',linewidth=3,label = label_temp)
plt.title("log p X log I 2D histogram")
plt.ylabel('log p')
plt.xlabel('log I')
# plt.legend()
plt.tight_layout()
plt.show()

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







# # Launch APLpy figure of 2D cube
# InputHDU = fits.open('scupollegacy_dr21_cube.fits')
# img = aplpy.FITSFigure(InputHDU[0],convention='wells',slices=[0])
# img.show_colorscale(vmin=0, vmax=0.02, cmap='gray') #, stretch='sqrt')


# # Modify the tick labels for precision and format
# img.axis_labels.set_font(size='xx-large')
# img.tick_labels.set_font(size='xx-large')

# #Revert to the default font size for the axis labels and for the ticks labels:
# img.axis_labels.set_font(size=None)
# img.tick_labels.set_font(size=None)