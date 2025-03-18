#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import aplpy
from aplpy import FITSFigure  
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from scipy.optimize import curve_fit
from astropy.stats import bootstrap
import matplotlib.ticker as ticker
import astropy.units as u 
from scipy.stats import bootstrap as bootstrapscipy
plt.rcParams.update({'font.size': 18})
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
import Functions
import math

# %matplotlib inline
#%%
FITS1 = '../2025/STOKES I.fits'
FITS2 = '../2025/ERROR I.fits'
FITS3 = '../2025/DEBIASED PERCENT POL.fits'
FITS4 = '../2025/ERROR PERCENT POL.fits'
FITS5 = '../2025/ROTATED POL ANGLE.fits'
FITS6 = '../2025/ERROR POL ANGLE.fits'
FITS7 = '../2025/DR21_full_NH2_Repr.fits'
FITS8 = '../2025/DR21_full_Tdust_Repr.fits'
FITS9 = '../2025/DR21_full_IRAC4_Repr.fits'
FITS10 = '../2025/DR21_full_Fil_Mask.fits'
FITS11 = '../2025/DR21_full_Fil_OF_Mask.fits'
FITS12 = '../2025/STOKES Q.fits'
FITS13 = '../2025/STOKES U.fits'



hdul1 = fits.open(FITS1)
hdul2 = fits.open(FITS2)
hdul3= fits.open(FITS3)
hdul4 = fits.open(FITS4)
hdul5 = fits.open(FITS5)
hdul6 = fits.open(FITS6)
hdul7 = fits.open(FITS7)
hdul8 = fits.open(FITS8)
hdul9 = fits.open(FITS9)
hdul10 = fits.open(FITS10)
hdul11 = fits.open(FITS11)
hdul12 =fits.open(FITS12)
hdul13 =fits.open(FITS13)






# print(hdul1.info())
# print(hdul2.info())
# print(hdul3.info())
# print(hdul4.info())
# print(hdul5.info())
# print(hdul6.info())
# print(hdul7.info())
# print(hdul8.info())
# print(hdul9.info())
# print(hdul10.info())
# print(hdul11.info())
# print(hdul12.info())
# print(hdul13.info())



MapStokesI_2025 = hdul1[0]
MapStokesIError_2025 = hdul2[1]
MapDebPol_2025 = hdul3[1]
MapDebPolError_2025 = hdul4[1]
MapPolAngle_2025 = hdul5[1]
MapPolAngleError_2025 = hdul6[1]
MapColumndensity_2025 = hdul7[0]
MapTemperature_2025 = hdul8[0]
Map8Micron_2025 = hdul9[0]
Mask_2025 = hdul10[0]
MaskOF_2025 = hdul11[0]
MapStokesQ_2025 = hdul12[1]
MapStokesU_2025 = hdul13[1]



# plt.figure(figsize=(6,8))
# # plt.imshow(CheckMapTemperature.data,origin='lower')
# plt.imshow(MapPolAngle_2025.data ,origin='lower')
# plt.tight_layout()
# plt.show()





MapPolSNR_2025 = MapDebPol_2025.copy()
BlankedMapPol_2025 = MapDebPol_2025.copy()
BlankedMapDebPolError_2025 = MapDebPolError_2025.copy()
BlankedMapPolAngle_2025 = MapPolAngle_2025.copy()
BlankedMapPolAngleError_2025 = MapPolAngleError_2025.copy()
BlankedMapStokesI_2025 = MapStokesI_2025.copy()
BlankedMapStokesIError_2025 = MapStokesIError_2025.copy()
BlankedMapStokesQ_2025 = MapStokesQ_2025.copy()
BlankedMapStokesU_2025 = MapStokesU_2025.copy()
BlankedMapColumnDensity_2025 = MapColumndensity_2025.copy()
BlankedMapTemperature_2025 = MapTemperature_2025.copy()
BlankedMap8Mircon_2025 = Map8Micron_2025.copy()


MapPolSNR_2025.data[:] = np.nan
MapPolSNR_2025.data = MapDebPol_2025.data/MapDebPolError_2025.data


Selector = (MapPolSNR_2025.data < 3)
# BlankedMapStokesI_2025.data[Selector] = np.nan
# BlankedMapStokesIError_2025.data[Selector] = np.nan
# BlankedMapDebPol_2025.data[Selector] = np.nan
# BlankedMapDebPolError_2025.data[Selector] = np.nan
# BlankedMapPolAngle_2025.data[Selector] = np.nan
# BlankedMapPolAngleError_2025.data[Selector] = np.nan
BlankedMapPol_2025.data[Selector] = np.nan
BlankedMapDebPolError_2025.data[Selector] = np.nan
BlankedMapPolAngle_2025.data[Selector] = np.nan
BlankedMapPolAngleError_2025.data[Selector] = np.nan
BlankedMapStokesI_2025.data[Selector] = np.nan
BlankedMapStokesIError_2025.data[Selector] = np.nan
BlankedMapStokesQ_2025.data[Selector] = np.nan
BlankedMapStokesU_2025.data[Selector] = np.nan
BlankedMapColumnDensity_2025.data[Selector] = np.nan
BlankedMapTemperature_2025.data[Selector] = np.nan
BlankedMap8Mircon_2025.data[Selector] = np.nan


############## removing any points with pfrac above 50
Selector = (BlankedMapPol_2025.data>50)
BlankedMapPol_2025.data[Selector] = np.nan
BlankedMapDebPolError_2025.data[Selector] = np.nan
BlankedMapPolAngle_2025.data[Selector] = np.nan
BlankedMapPolAngleError_2025.data[Selector] = np.nan
BlankedMapStokesI_2025.data[Selector] = np.nan
BlankedMapStokesIError_2025.data[Selector] = np.nan
BlankedMapStokesQ_2025.data[Selector] = np.nan
BlankedMapStokesU_2025.data[Selector] = np.nan
BlankedMapColumnDensity_2025.data[Selector] = np.nan
BlankedMapTemperature_2025.data[Selector] = np.nan
BlankedMap8Mircon_2025.data[Selector] = np.nan



############ removing any data points with I/I_error < 100
Selector = MapStokesI_2025.data/MapStokesIError_2025.data < 100
BlankedMapPol_2025.data[Selector] = np.nan
BlankedMapDebPolError_2025.data[Selector] = np.nan
BlankedMapPolAngle_2025.data[Selector] = np.nan
BlankedMapPolAngleError_2025.data[Selector] = np.nan
BlankedMapStokesI_2025.data[Selector] = np.nan
BlankedMapStokesIError_2025.data[Selector] = np.nan
BlankedMapStokesQ_2025.data[Selector] = np.nan
BlankedMapStokesU_2025.data[Selector] = np.nan
BlankedMapColumnDensity_2025.data[Selector] = np.nan
BlankedMapTemperature_2025.data[Selector] = np.nan
BlankedMap8Mircon_2025.data[Selector] = np.nan


# plt.figure(figsize=(6,8))
# # plt.imshow(CheckMapTemperature.data,origin='lower')
# plt.imshow(MapStokesI_2025.data ,origin='lower')
# plt.tight_layout()
# plt.show()



# plt.figure(figsize=(6,8))
# # plt.imshow(CheckMapTemperature.data,origin='lower')
# plt.imshow(BlankedMapStokesI_2025.data ,origin='lower')
# plt.tight_layout()
# plt.show()


DEC_grid,RA_grid = Functions.generate_RA_DEC_mesh(BlankedMapStokesI_2025)
seperation = MapPolAngle_2025.copy()

# print(RA_grid2025.shape)
# print(DEC_grid2025.shape)

# CheckMapTemperature = BlankedMapTemperature.copy()
# selector = BlankedMapTemperature.data>40
# CheckMapTemperature.data[selector] = 5
# CheckMapTemperature = BlankedMapTemperature.copy()
# selector = BlankedMapTemperature.data<30
# CheckMapTemperature.data[selector] = 5

# print(1)
# plt.figure(figsize=(6,8))
# # plt.imshow(CheckMapTemperature.data,origin='lower')
# plt.imshow(BlankedMapColumnDensity_2025.data ,origin='lower',vmin=80)
# plt.tight_layout()
# plt.show()


# plt.figure(figsize=(6,8))
# # plt.imshow(CheckMapTemperature.data,origin='lower')
# im = plt.imshow(BlankedMapPol_2025.data ,origin='lower')
# plt.colorbar(im)
# plt.title("Debiased P% 2025")
# # plt.tight_layout()
# plt.show()


# plt.figure(figsize=(6,8))
# # plt.imshow(CheckMapTemperature.data,origin='lower')
# im = plt.imshow(BlankedMapStokesI_2025.data ,origin='lower')
# plt.colorbar(im)
# plt.title("Stokes I 2025")
# # plt.tight_layout()
# plt.show()


# plt.figure(figsize=(6,8))
# # plt.imshow(CheckMapTemperature.data,origin='lower')
# im = plt.imshow(BlankedMapStokesQ_2025.data ,origin='lower',cmap = 'RdBu')
# plt.colorbar(im)
# plt.title("Stokes Q 2025")
# plt.tight_layout()
# plt.show()

# plt.figure(figsize=(6,8))
# # plt.imshow(CheckMapTemperature.data,origin='lower')
# im = plt.imshow(BlankedMapStokesU_2025.data ,origin='lower',cmap = 'RdBu')
# plt.colorbar(im)
# plt.title("Stokes U 2025")
# plt.tight_layout()
# plt.show()

# plt.figure()
# plt.hist(BlankedMapStokesQ_2025.data.flatten(),50)
# plt.title('Distribution of Stokes Q 2025')
# plt.xlabel('Stokes Q')
# plt.ylabel('Number of data points')
# plt.show()

# plt.figure()
# plt.hist(BlankedMapStokesU_2025.data.flatten(),50)
# plt.title('Distribution of Stokes U 2025')
# plt.xlabel('Stokes U')
# plt.ylabel('Number of data points')
# plt.show()

#%%
# dense_cores = pd.read_csv('../data/dense_cores - Copy.tsv',delimiter='\t',header=None)
# ra_core = np.array(dense_cores[7])
# dec_core = np.array(dense_cores[8])
# FWHMa = np.array(dense_cores[2])
# FWHMb = np.array(dense_cores[3])
# #DR21 clump N46
# #DR21OH clump N44
# #DR21OH-W clump N38
# #DR21OH-S clump N47



# angle_FWHMa = np.degrees(np.arctan(FWHMa/(2*1700)))
# angle_3FWHMa = np.degrees(np.arctan(3*FWHMa/(2*1700)))

# angle_FWHMb = (FWHMb/1700)
# print(angle_FWHMa)
# plt.hist(angle_FWHMa,50)

# def Calc_l(ra1,dec1,ra2,dec2):

#     c1 = SkyCoord(ra1,dec1,unit = 'deg')
#     c2 = SkyCoord(ra2,dec2,unit = 'deg')
#     sep = c1.separation(c2)
#     return sep.deg

# mask_core = np.zeros_like(RA_grid)
# for i in range(ra_core.shape[0]):
# 	seperation = Calc_l(ra_core[i],dec_core[i],RA_grid,DEC_grid)
# 	mask_core[seperation<angle_3FWHMa[i]] = 1
# 	mask_core = mask_core.astype(bool)



# plt.figure()
# plt.imshow(mask_core,origin='lower')
# plt.colorbar()
# plt.show()
	
#%%

def get_data(infile1,infile2,infile3,infile4,infile5):
	#get data
	hawc1 = fits.open(infile1)
	hawc2 = fits.open(infile2)
	hawc3 = fits.open(infile3)
	hawc4 = fits.open(infile4)
	hawc5 = fits.open(infile5)
	p    = hawc1[1]
	perr = hawc2[1]
	#pa   = hawc['ROTATED POL ANGLE']
	pa   = hawc3[1]
	stkI = hawc4[0]
	stkIerr = hawc5[1]
	# pi   = hawc[15]
	# pierr   = hawc[14]

	#Jy/px to Jy/sqarcsec
	pxscale = stkI.header['CDELT2']*3600
	stkI.data /= pxscale**2
	# pi.data /= pxscale**2 
	stkIerr.data /= pxscale**2
	return p,perr,pa,stkI,stkIerr,pxscale

def quality_cuts(stkI,stkIerr,p,perr,SNRp_cut,p_cut,SNRi_cut):
	#snr in P
	SNRp = p.data/perr.data
	mask_snrp = np.where(SNRp < SNRp_cut)
	p.data[mask_snrp] = np.nan
	#p_cut
	maskp = np.where(p.data > p_cut)
	p.data[maskp] = np.nan
	#snr in P
	SNRi = stkI.data/stkIerr.data
	mask_snri = np.where(SNRi < SNRi_cut)
	p.data[mask_snri] = np.nan
	return p


# filename1='../2025/STOKES I.fits'
# filename2='../2025/DEBIASED PERCENT POL.fits'
# filename3='../2025/ROTATED POL ANGLE.fits'
# filename4='../2025/STOKES I.fits'



# hawcp1 = fits.open(filename1)
# hawcp2 = fits.open(filename2)
# hawcp3 = fits.open(filename3)
# hawcp4 = fits.open(filename4)

# Herschel = fits.open(filename5)
# MapHer250 = Herschel[0]

# # x = hawcp['STOKES I'].data
# # y = hawcp['DEBIASED PERCENT POL'].data
# # y1=hawcp['ROTATED POL ANGLE'].data
# # y2 = hawcp['DEBIASED POL FLUX'].data


filename= '../2025/DR21_full_NH2_Repr.fits'
Herschel = fits.open(filename)
MapHer250 = Herschel[0]

MapHer250log = MapHer250.copy()
MapHer250log.data = np.log10(MapHer250.data)
title = 'SIMPLIFI'

# plt.figure()
# plt.hist(MapHer250log.data)
# plt.show()

SNRp_cut = 3.0
p_cut = 50
width  = 50
height = 50
cmap = 'plasma'
title_size = 16
tick_labels = 15
label_plot = 15
label_colorbar = 15
tick_colorbar = 15
label_fontsize = 20
SNRi_cut = 100
scalevec = 1.5 #1px = scalevec * 1% pol 
vec_legend = 5.0

filename1 = '../2025/DEBIASED PERCENT POL.fits'
filename2 = '../2025/ERROR PERCENT POL.fits'
filename3 = '../2025/ROTATED POL ANGLE.fits'
filename4 = '../2025/STOKES I.fits'
filename5 = '../2025/ERROR I.fits'
p,perr,pa,stkI,stkIerr,pxscale = get_data(filename1,filename2,filename3,filename4,filename5)
p = quality_cuts(stkI,stkIerr,p,perr,SNRp_cut,p_cut,SNRi_cut)

# RA = (stkI.header['OBSRA']*u.hourangle).to(u.deg)
# DEC = stkI.header['OBSDEC']*u.deg


#### SCRIPT
fig = plt.figure(figsize=(13,10))
gc = FITSFigure(MapHer250log,figure=fig,vmin = 22.1,vmax = 24.5)
gc.show_colorscale(cmap='default')
gc.add_colorbar(location='right', width=0.2, pad=0.25, ticks=None,axis_label_text= 'log Column Density')
gc.show_contour(colors = 'white',levels = 7)
# gc.show_vectors(p,pa,scale=scalevec,step=4,color='black',linewidth=3.5)
# gc.show_vectors(p,pa,scale=scalevec,step=4,color='yellow',linewidth=2.0)
vecscale = scalevec * pxscale/3600
gc.add_scalebar(vec_legend*vecscale,r'$p\%$ ='+str(vec_legend),corner='bottom right',frame=True,color='black')
plt.savefig('../plots/2025/19/Column_density.png')
plt.close()
# plt.show()


# fig = plt.figure(figsize=(13,10))
# gc = FITSFigure(MapHer250log,figure=fig,vmin = 22.1,vmax = 24.5)
# gc.show_colorscale(cmap='default')
# gc.add_colorbar(location='right', width=0.2, pad=0.25, ticks=None,axis_label_text= 'log Column Density')
# gc.show_contour(colors = 'white',levels = 7)
# # gc.show_circles(ra_core,dec_core,angle_FWHMa)
# # gc.show_vectors(p,pa,scale=scalevec,step=4,color='black',linewidth=3.5)
# # gc.show_vectors(p,pa,scale=scalevec,step=4,color='yellow',linewidth=2.0)
# vecscale = scalevec * pxscale/3600
# gc.add_scalebar(vec_legend*vecscale,r'$p\%$ ='+str(vec_legend),corner='bottom right',frame=True,color='black')
# plt.show()

# %%

###### making the S map
print('started calculating angle dispersion')
set_delta = (19.2)/60   # in arcminute
S_map = BlankedMapPolAngle_2025.copy()
sigma_S_map = BlankedMapPolAngleError_2025.copy()

for i in range(RA_grid.shape[0]):
    for j in range(RA_grid.shape[1]):
# for i in range(10):
#     for j in range(10):

        ##### seperation filter
        seperation.data = Functions.Calc_l(RA_grid[i,j],DEC_grid[i,j],RA_grid,DEC_grid)
        seperation_selector = (seperation.data<0.5*set_delta)
        seperation.data[seperation_selector] = np.nan
        seperation_selector = (seperation.data>1.5*set_delta)
        seperation.data[seperation_selector] = np.nan
        seperation_selector = (seperation.data >0)

        ##### making the dispersion map
        tempa = BlankedMapStokesQ_2025.data*BlankedMapStokesU_2025.data[i,j] - BlankedMapStokesQ_2025.data[i,j]*BlankedMapStokesU_2025.data
        tempb = BlankedMapStokesQ_2025.data*BlankedMapStokesQ_2025.data[i,j] + BlankedMapStokesU_2025.data*BlankedMapStokesU_2025.data[i,j]
        AngleDiff_v2 = 0.5 * (180/np.pi)*np.arctan2(tempa,tempb)
        S = np.nanmean(AngleDiff_v2[seperation_selector]**2)**0.5
        S_map.data[i,j] = S

        ##### making the dispersion error map
        # sigma_S = np.nanmean(BlankedMapPolAngleError.data[seperation_selector]**2)**0.5
        # sigma_S_map.data[i,j] = sigma_S
        N = np.nansum(~np.isnan(seperation.data))
        # print(N)
        temp1 = (( 1 /(S_map.data[i,j]*N) )**2)
        # temp2 = np.nansum((AngleDiff_v2[seperation_selector])**2)*(BlankedMapPolAngleError_2025.data[i,j]**2)
        # temp3 =  np.nansum((AngleDiff_v2[seperation_selector]*BlankedMapPolAngleError_2025.data[seperation_selector])**2)

        temp2 = (BlankedMapPolAngleError_2025.data[i,j]**2)*(np.nansum(AngleDiff_v2[seperation_selector])**2)
        temp3 =  np.nansum(((BlankedMapPolAngleError_2025.data[seperation_selector])**2)*(AngleDiff_v2[seperation_selector]**2))

        
        # temp1 = ( ( BlankedMapPolAngleError.data[i,j]/( S_map.data[i,j]*N ) )**2)*(np.nansum(AngleDiff_v2[seperation_selector])**2)
        # temp2 = (BlankedMapPolAngleError.data*AngleDiff_v2)
        # temp3 = (( 1 /(S_map.data[i,j]*N) )**2)*( np.nansum( (temp2[seperation_selector])**2 ) )

        sigma_S =  temp1*(temp2 + temp3)
        sigma_S_map.data[i,j] = np.sqrt(sigma_S)

S_map_deb = S_map.copy()

S_map_deb.data = np.sqrt(S_map.data**2 - sigma_S_map.data**2)

# S_map_deb.writeto('../FITS_file/new_fits/S_map_deb_new_2025v2.fits')
# sigma_S_map.writeto('../FITS_file/new_fits/sigma_S_new_2025v2.fits')
print('Finished calculating angle dispersion')
#%%
# S_map_deb =  fits.open('../FITS_file/new_fits/S_map_deb_new_2025v2.fits')[1]
# sigma_S_map = fits.open('../FITS_file/new_fits/sigma_S_new_2025v2.fits')[1]

Selector = S_map_deb.data/sigma_S_map.data < 3
S_map_deb.data[Selector] = np.nan
sigma_S_map.data[Selector] = np.nan


# Selector = S_map_deb1.data/sigma_S_map1.data < 5
# S_map_deb1.data[Selector] = np.nan
# sigma_S_map1.data[Selector] = np.nan


fig = plt.figure(figsize=(13,10))
gc1 = FITSFigure(S_map_deb,figure=fig)
gc1.show_colorscale(cmap='default')
# gc1.show_circles(ra_core,dec_core,angle_FWHMa)
gc1.add_colorbar(location='top', width=0.2, pad=0.15, ticks=None,axis_label_text= 'Angle Dispersion')
gc1.show_contour(colors = 'white',levels = 5)
plt.tight_layout()
plt.savefig('../plots/2025/19/Angle_dispersion.png')
plt.close()
# plt.show()


#%%

plt.figure()
plt.hist(S_map_deb.data.flatten(),50)
plt.title('Distribution of debiased angle dispersion 2025 data')
plt.xlabel('debiased S')
plt.ylabel('Number')
plt.savefig('../plots/2025/19/hist_debiased_S.png')
plt.close()
# plt.show()

plt.figure()
plt.hist(sigma_S_map.data.flatten(),50)
plt.title('Distribution of debiased angle dispersion 2025 data')
plt.xlabel('debiased S')
plt.ylabel('Number')
plt.savefig('../plots/2025/19/hist_sigma_S.png')
plt.close()
# plt.show()

# fig = plt.figure(figsize=(13,10))
# gc1 = FITSFigure(S_map_deb1,figure=fig)
# gc1.show_colorscale(cmap='default')
# gc1.add_colorbar(location='top', width=0.2, pad=0.15, ticks=None,axis_label_text= 'Angle Dispersion')
# gc1.show_contour(colors = 'white',levels = 5)
# plt.tight_layout()
# plt.show()

# %%

######### Case 1

############### changing data to 1D array

s_array = S_map_deb.data.flatten()[::4]
es_array = sigma_S_map.data.flatten()[::4]
p_array = BlankedMapPol_2025.data.flatten()[::4]
ep_array = BlankedMapDebPolError_2025.data.flatten()[::4]
I_array = BlankedMapStokesI_2025.data.flatten()[::4]
eI_array = BlankedMapStokesIError_2025.data.flatten()[::4]
nh2_array = BlankedMapColumnDensity_2025.data.flatten()[::4]
micron8_array = BlankedMap8Mircon_2025.data.flatten()[::4]
temp_array = BlankedMapTemperature_2025.data.flatten()[::4]

valid = np.isnan(s_array)|np.isnan(es_array)|np.isnan(p_array)|np.isnan(ep_array)|np.isnan(I_array)|np.isnan(eI_array)|np.isnan(nh2_array)|np.isnan(micron8_array)|np.isnan(temp_array)
s_array = s_array[~valid]	
es_array = es_array[~valid]
p_array = p_array[~valid]
ep_array = ep_array[~valid]
I_array = I_array[~valid]
eI_array = eI_array[~valid]
nh2_array = nh2_array[~valid]
micron8_array = micron8_array[~valid]
temp_array = temp_array[~valid]
#%%

# from multiprocessing import Pool

from joblib import Parallel, delayed
import numpy as np
from scipy.optimize import curve_fit

def get_wstd(vals, weights):
	vals_n_weights = [(vals[i], weights[i]) for i in range(0, len(weights))]
	def get_wmean(vals, weights):
		weighted_vals = []
		for tup in vals_n_weights:
			weighted_vals.append(round(tup[0]*tup[1]/sum(weights), 3))    
		answer = sum(weighted_vals)
		return answer 

	numerator = []
	for i in range(0, len(weights)):
		numerator.append(weights[i]*(vals[i]-get_wmean(vals, weights))**2)
	var = sum(numerator)/(((len(weights)-1)*sum(weights))/len(weights))
	wstdev = math.sqrt(var)
	return wstdev

def binning_equal_width(array1,array2,error_array1,error_array2,Nbins,min_n = 1):
	Nbins += 1 
	log_filtered1 = np.log10(array1)
	log_filtered2 = np.log10(array2)

	bins = np.linspace(np.nanmin(log_filtered1),np.nanmax(log_filtered1),Nbins)

	bin_centres = []
	binned_data = []
	error_bar = []

	for i in range(0,(bins.shape[0]-1)):
		temp_array1 = array1.copy()
		temp_array2 = array2.copy()
		temp_error1 = error_array1.copy()
		temp_error2 = error_array2.copy()

		Selector = (log_filtered1 > bins[i])&(log_filtered1 < bins[i+1])
		temp_array1_fil = temp_array1[Selector]
		temp_array2_fil = temp_array2[Selector] 
		temp_error1_fil = temp_error1[Selector]
		temp_error2_fil = temp_error2[Selector] 


		valid = ~(np.isnan(temp_array1_fil)|np.isnan(temp_array2_fil)|np.isnan(temp_error1_fil)|np.isnan(temp_error2_fil))
		temp_array1_fil = temp_array1_fil[valid]
		temp_array2_fil = temp_array2_fil[valid]
		temp_error1_fil = temp_error1_fil[valid]
		temp_error2_fil = temp_error2_fil[valid]

		# print(np.sum(valid))
		if np.sum(valid) == 0:
			bin_centres.append(np.nan)
			binned_data.append(np.nan)
			error_bar.append(np.nan)

		elif (np.sum(valid) < min_n):
			bin_centres.append(np.nan)
			binned_data.append(np.nan)
			error_bar.append(np.nan)

		elif (np.sum(valid)  ==  1):
			bin_centres.append(temp_array1_fil[0])
			binned_data.append(temp_array2_fil[0])
			error_bar.append(temp_error2_fil[0])

		else:
			bin_centre_temp = np.average(temp_array1_fil,weights= 1.0/temp_error1_fil**2)
			bin_data_temp = np.average(temp_array2_fil,weights=1.0/temp_error2_fil**2)
			binned_error = get_wstd(temp_array2_fil,weights=1.0/temp_error2_fil**2)
			
			bin_centres.append(bin_centre_temp)
			binned_data.append(bin_data_temp)
			error_bar.append(binned_error)

	bin_centres = np.array(bin_centres)
	binned_data = np.array(binned_data)
	error_bar = np.array(error_bar)
	return bin_centres,binned_data,error_bar	

def plotting_equal_width_bins(array1,array2,error_array1,error_array2,Nbins,min_n = 1,color_bins='#06402B',color_error='#8C564B'):
	bin_centres,binned_data,error_bar = binning_equal_width(array1,array2,error_array1,error_array2,Nbins,min_n = 1)
	plt.errorbar(bin_centres,binned_data,yerr=error_bar,c=color_error)
	plt.scatter(bin_centres,binned_data,s = 75,c=color_bins)
	plt.plot(bin_centres,binned_data,c=color_bins)

def fitting_binned(array1,array2,error_array1,error_array2,Nbins,min_n = 1):

	def lin_fit(x, a, b):
		return a + b*x
	bin_centres,binned_data,error_bar = binning_equal_width(array1,array2,error_array1,error_array2,Nbins,min_n = 1)
	bin_centres = np.log10(np.array(bin_centres))
	binned_data = np.log10(np.array(binned_data))
	valid = ~(np.isnan(bin_centres)|np.isnan(binned_data))

	param, PS_param_cov = curve_fit(lin_fit, bin_centres[valid], binned_data[valid])
	level_bins = np.linspace(np.nanmin(bin_centres),np.nanmax(bin_centres),11)
	PS_FitFunc = lin_fit(level_bins,param[0],param[1])
	return param,10**level_bins,10**PS_FitFunc

def fitting_binned_bootstrap(array1,array2,error_array1,error_array2,Nbins,min_n = 1):

	def lin_fit(x, a, b):
		return a + b*x
	bin_centres,binned_data,error_bar = binning_equal_width(array1,array2,error_array1,error_array2,Nbins,min_n = 1)
	bin_centres = np.log10(np.array(bin_centres))
	binned_data = np.log10(np.array(binned_data))
	valid = ~(np.isnan(bin_centres)|np.isnan(binned_data))

	param, PS_param_cov = curve_fit(lin_fit, bin_centres[valid], binned_data[valid])
	level_bins = np.linspace(np.nanmin(bin_centres),np.nanmax(bin_centres),11)
	PS_FitFunc = lin_fit(level_bins,param[0],param[1])
	return param[1]


def bootstrap_fitting(array1, array2, error_array1, error_array2, Nbins, min_n=1, n_bootstrap=1000):
	def lin_fit(x, a, b):
		return a + b*x

	n_data = len(array1)
	bootstrap_indices = np.random.choice(n_data, (n_bootstrap, n_data), replace=True)

	def bootstrap_iteration(indices):
		"""Performs a single bootstrap iteration and returns the slope parameter."""
		resampled_array1 = array1[indices]
		resampled_array2 = array2[indices]
		resampled_error1 = error_array1[indices]
		resampled_error2 = error_array2[indices]

		try:
			param = fitting_binned_bootstrap(resampled_array1, resampled_array2, resampled_error1, resampled_error2, Nbins, min_n)
			return param  # Return only the slope
		except:
			return None  # Return None for failed fits

	# Run bootstrap iterations in parallel
	params_bootstrap = Parallel(n_jobs=-1, backend="loky")(
		delayed(bootstrap_iteration)(indices) for indices in bootstrap_indices
	)

	# Remove failed fits (None values)
	params_bootstrap = np.array([p for p in params_bootstrap if p is not None])

	if params_bootstrap.shape[0] == 0:
		raise ValueError("Bootstrap failed: No successful fits.")

	param_mean = np.mean(params_bootstrap, axis=0)
	param_std = np.std(params_bootstrap, axis=0)
	# plt.figure()
	# plt.hist(params_bootstrap,50)
	# plt.show()
	return param_mean, param_std


#%%
#### single paramter fitting
print('plotting I_vs_N')
array1 = nh2_array
array2 = I_array
error_array1 = np.ones_like(nh2_array)
error_array2 = eI_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2)
plotting_equal_width_bins(array1, array2,error_array1,error_array2,10)
param,x,y  = fitting_binned(array1,array2,error_array1,error_array2,10)
data = (array1,array2,error_array1,error_array2)
param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
label_temp = r'log(I) = C + $\alpha$log($N$){linebreak}  $\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(frac = r'$\%$',I = r'$_{I}$',alpha = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label =label_temp)
plt.title(r'$I$ vs $N$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel(r'$I$')
plt.tight_layout()
plt.legend()
plt.savefig('../plots/2025/19/I_vs_N.png')
plt.close()
# plt.show()


#%%
print('plotting p_vs_I')

array1 = I_array
array2 = p_array
error_array1 = eI_array
error_array2 = ep_array



plt.figure(figsize=(12,8))
plt.scatter(array1,array2)
plotting_equal_width_bins(array1, array2,error_array1,error_array2,10)
param,x,y  = fitting_binned(array1,array2,error_array1,error_array2,10)
data = (array1,array2,error_array1,error_array2)
param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
label_temp = r'log(p{frac}) = C + $\alpha$log($I$){linebreak}  $\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(frac = r'$\%$',I = r'$_{I}$',alpha = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
plt.title(r'$p\%$ vs $I$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('I')
plt.ylabel(r'$p\%$')
plt.tight_layout()
plt.legend()
plt.savefig('../plots/2025/19/p_vs_I.png')
plt.close()
# plt.show()


# plt.figure(figsize=(12,8))
# plt.scatter(I_array,p_array)
# plotting_equal_width_bins(I_array,p_array,eI_array,ep_array,10)
# param,x,y  = fitting_binned(I_array,p_array,eI_array,ep_array,10)
# print(param[1])
# data = (I_array,p_array,eI_array,ep_array)
# param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
# label_temp = r'log(p{frac}) = C + $\alpha$log($I$){linebreak}  $\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(frac = r'$\%$',I = r'$_{I}$',alpha = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = 'alpha = {alpha:.3f}'.format(alpha = param[1]))
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
# plt.title(r'$p\%$ vs $I$')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('I')
# plt.ylabel(r'$p\%$')
# plt.tight_layout()
# plt.legend()
# plt.show()

#%%
print('plotting p_vs_N')
array1 = nh2_array
array2 = p_array
error_array1 = np.ones_like(nh2_array)
error_array2 = ep_array



plt.figure(figsize=(12,8))
plt.scatter(array1,array2)
plotting_equal_width_bins(array1, array2,error_array1,error_array2,10)
param,x,y  = fitting_binned(array1,array2,error_array1,error_array2,10)
data = (array1,array2,error_array1,error_array2)
param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
label_temp = r'log(p{frac}) = C + $\alpha$log($N$){linebreak}  $\alpha$: {alpha_nh2:.3f}{pm}{error:.3f}'.format(frac = r'$\%$',nh2 = r'$_{N}$',alpha_nh2= param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
plt.title(r'$p\%$ vs $N$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel(r'$p\%$')
plt.tight_layout()
plt.legend()
plt.savefig('../plots/2025/19/p_vs_N.png')
plt.close()
# plt.show()



# plt.figure(figsize=(12,8))
# plt.scatter(nh2_array,p_array,c = "#1F77B4")
# plotting_equal_width_bins(nh2_array,p_array,np.ones_like(nh2_array),ep_array,10)
# param,x,y  = fitting_binned(nh2_array,p_array,np.ones_like(nh2_array),ep_array,10)
# print(param[1])

# # data = (nh2_array,p_array,np.ones_like(nh2_array),ep_array)
# # param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
# # label_temp = r'log(p{frac}) = C + $\alpha$log($N$){linebreak}  $\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(frac = r'$\%$',I = r'$_{I}$',alpha = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = 'alpha = {alpha:.3f}'.format(alpha = param[1]))
# plt.title(r'$p\%$ vs $N$')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('N')
# plt.ylabel(r'$p\%$')
# plt.tight_layout()
# plt.legend()
# plt.show()


#%%

print('plotting S_vs_P')

array1 = p_array
array2 = s_array
error_array1 = ep_array
error_array2 = es_array



plt.figure(figsize=(12,8))
plt.scatter(array1,array2)
plotting_equal_width_bins(array1, array2,error_array1,error_array2,10)
param,x,y  = fitting_binned(array1,array2,error_array1,error_array2,10)
data = (array1,array2,error_array1,error_array2)
param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
label_temp = r'log(S) = C + $\alpha$log(p{frac}){linebreak} $\alpha$: {alpha_s:.3f}{pm}{error:.3f}'.format(frac = r'$\%$',alpha_s = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
plt.title(r'$S$ vs $p\%$')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('S')
plt.xlabel(r'$p\%$')
plt.tight_layout()
plt.legend()
plt.savefig('../plots/2025/19/S_vs_p.png')
plt.close()
# plt.show()

# s_bins = np.linspace(np.nanmin(np.log10(s_array)),np.nanmax(np.log10(s_array)),25)
# p_bins = np.linspace(np.nanmin(np.log10(p_array)),np.nanmax(np.log10(p_array)),25)

# plt.figure(figsize=(12,8))
# plt.scatter(p_array,s_array)
# plotting_equal_width_bins(p_array,s_array,ep_array,es_array,10)
# param,x,y  = fitting_binned(p_array,s_array,ep_array,es_array,10)
# print(param[1])
# data = (p_array,s_array,ep_array,es_array)
# param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
# label_temp = r'log(S) = C + $\alpha$log(p{frac}){linebreak} $\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(frac = r'$\%$',I = r'$_{I}$',alpha = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = 'alpha = {alpha:.3f}'.format(alpha = param[1]))
# plt.title(r'$S$ vs $p\%$')
# plt.xscale('log')
# plt.yscale('log')
# plt.ylabel('S')
# plt.xlabel(r'$p\%$')
# plt.tight_layout()
# plt.legend()
# plt.show()


#%%
print('plotting S_vs_N')

array1 = nh2_array
array2 = s_array
error_array1 = np.ones_like(nh2_array)
error_array2 = es_array



plt.figure(figsize=(12,8))
plt.scatter(array1,array2)
plotting_equal_width_bins(array1, array2,error_array1,error_array2,10)
param,x,y  = fitting_binned(array1,array2,error_array1,error_array2,10)
data = (array1,array2,error_array1,error_array2)
param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
label_temp = r'log(S) = C + $\alpha$log($N$){linebreak}$\alpha$: {alpha_nh2:.3f}{pm}{error:.3f}'.format(frac = r'$\%$',alpha_nh2 = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
plt.title(r'$S$ vs $N$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel(r'$S$')
plt.tight_layout()
plt.legend()
plt.savefig('../plots/2025/19/S_vs_N.png')
plt.close()
# plt.show()


# plt.figure(figsize=(16,8))
# plt.scatter(nh2_array,s_array)

# plotting_equal_width_bins((nh2_array),(s_array),np.ones_like(nh2_array),es_array,10)
# param,x,y  = fitting_binned((nh2_array),(s_array),np.ones_like(nh2_array),es_array,10)
# print(param[1])

# # data = ((nh2_array),(s_array),np.ones_like(nh2_array),es_array)
# # param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
# # label_temp = r'log(S) = C + $\alpha$log($N$){linebreak}$\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(frac = r'$\%$',I = r'$_{I}$',alpha = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
# # plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = 'alpha = {alpha:.3f}'.format(alpha = param[1]))
# plt.title(r'$S$ vs $N$')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('N')
# plt.ylabel(r'$S$')
# plt.tight_layout()
# plt.legend()
# plt.show()


#%%

print('plotting S_vs_I')

array1 = I_array
array2 = s_array
error_array1 = eI_array
error_array2 = es_array



plt.figure(figsize=(12,8))
plt.scatter(array1,array2)
plotting_equal_width_bins(array1, array2,error_array1,error_array2,10)
param,x,y  = fitting_binned(array1,array2,error_array1,error_array2,10)
data = (array1,array2,error_array1,error_array2)
param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
label_temp = r'log(S) = C + $\alpha$log($I$){linebreak}$\alpha$: {alpha_nh2:.3f}{pm}{error:.3f}'.format(frac = r'$\%$',alpha_nh2 = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
plt.title(r'$S$ vs $I$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('I')
plt.ylabel(r'$S$')
plt.tight_layout()
plt.legend()
plt.savefig('../plots/2025/19/S_vs_I.png')
plt.close()
# plt.show()


# plt.figure(figsize=(16,8))
# plt.scatter(I_array,s_array)


# plotting_equal_width_bins((I_array),(s_array),eI_array,es_array,10)
# param,x,y  = fitting_binned((I_array),(s_array),eI_array,es_array,10)
# print(param[1])
# # data = ((I_array),(s_array),eI_array,es_array)
# # param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
# # label_temp = r'log(S) = C + $\alpha$log($N$){linebreak}$\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(frac = r'$\%$',I = r'$_{I}$',alpha = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')

# # label_temp = r'log(S) = C + $\alpha$log($I$){linebreak}$\alpha$: {alpha_nh2:.3f}{pm}{error:.3f}'.format(frac = r'$\%$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# # plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = 'alpha = {alpha:.3f}'.format(alpha = param[1]))
# plt.title(r'$S$ vs $I$')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('I')
# plt.ylabel(r'$S$')
# plt.tight_layout()
# plt.legend()
# plt.show()

#%%
print('plotting pS_vs_N')
array1 = nh2_array
array2 = p_array*s_array
error_array1 = np.ones_like(nh2_array)
error_array2 =  np.sqrt((p_array*es_array)**2+ (ep_array*s_array)**2)



plt.figure(figsize=(12,8))
plt.scatter(array1,array2)
plotting_equal_width_bins(array1, array2,error_array1,error_array2,10)
param,x,y  = fitting_binned(array1,array2,error_array1,error_array2,10)
data = (array1,array2,error_array1,error_array2)
param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
label_temp = r'log(p{frac}*S) = C + $\alpha$log($N$){linebreak}$\alpha$: {alpha_nh2:.3f}{pm}{error:.3f}'.format(frac = r'$\%$',alpha_nh2 = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
plt.title(r'$p\%*S$ vs $N$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel(r'$p\%*S$')
plt.tight_layout()
plt.legend()
plt.savefig('../plots/2025/19/pS_vs_N.png')
plt.close()
# plt.show()



# p_array_copy = p_array.copy()
# s_array_copy = s_array.copy()
# error = np.sqrt((p_array*es_array)**2+ (ep_array*s_array)**2)
# # Selector = (np.log10(nh2_array)<22)|(np.log10(nh2_array)>23)
# # p_array_copy[Selector] = np.nan
# # s_array_copy[Selector] = np.nan
# plt.figure(figsize=(16,8))
# plt.scatter(nh2_array,p_array*s_array)
# plotting_equal_width_bins((nh2_array),(p_array*s_array),np.ones_like(nh2_array),error,10)
# param,x,y  = fitting_binned((nh2_array),(p_array*s_array),np.ones_like(nh2_array),error,10)
# print(param[1])
# # param,x,y = Functions.binning_equal_width((nh2_array),(p_array*s_array),np.ones_like(nh2_array),(p_array*es_array + ep_array*s_array),15)
# # param,x,y = Functions.binning_equal_width((nh2_array),(p_array*s_array),np.ones_like(nh2_array),error,10)
# # valid = ~(np.isnan(nh2_array)|np.isnan(p_array*s_array)|np.isnan(np.ones_like(nh2_array))|np.isnan(p_array*es_array + ep_array*s_array))
# # valid = ~(np.isnan(nh2_array)|np.isnan(p_array*s_array)|np.isnan(np.ones_like(nh2_array))|np.isnan(error))
# # data = ((nh2_array)[valid],(p_array*s_array)[valid],np.ones_like(nh2_array)[valid],(p_array*es_array + ep_array*s_array)[valid])
# # data = ((nh2_array)[valid],(p_array*s_array)[valid],np.ones_like(nh2_array)[valid],(error)[valid])
# # label_temp = r'log(p{frac}*S) = C + $\alpha$log($N$){linebreak}$\alpha$: {alpha_nh2:.3f}{pm}{error:.3f}'.format(frac = r'$\%$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = 'alpha = {alpha:.3f}'.format(alpha = param[1]))
# # plt.plot(x,np.ones_like(x)*np.nanmean(p_array_copy*s_array_copy),'tab:red',linewidth = 2.5,label = "mean : {mean:.3f}".format(mean = np.nanmean(p_array_copy*s_array_copy)))
# plt.title(r'$p\%*S$ vs $N$')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('N')
# plt.ylabel(r'$p\%*S$')
# plt.tight_layout()
# plt.legend()
# plt.show()

#%%
print('plotting pS_vs_I')
array1 = I_array
array2 = p_array*s_array
error_array1 = eI_array
error_array2 =  np.sqrt((p_array*es_array)**2+ (ep_array*s_array)**2)



plt.figure(figsize=(12,8))
plt.scatter(array1,array2)
plotting_equal_width_bins(array1, array2,error_array1,error_array2,10)
param,x,y  = fitting_binned(array1,array2,error_array1,error_array2,10)
data = (array1,array2,error_array1,error_array2)
param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
label_temp = r'log(p{frac}*S) = C + $\alpha$log($I$){linebreak}$\alpha$: {alpha_nh2:.3f}{pm}{error:.3f}'.format(frac = r'$\%$',alpha_nh2 = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
plt.title(r'$p\%*S$ vs $N$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel(r'$p\%*S$')
plt.tight_layout()
plt.legend()
plt.savefig('../plots/2025/19/pS_vs_I.png')
plt.close()
# plt.show()


# I_array_copy = I_array.copy()

# error = np.sqrt((p_array*es_array)**2+ (ep_array*s_array)**2)

# # p_array_copy = p_array.copy()
# # s_array_copy = s_array.copy()
# # Selector = (np.log10(nh2_array)<22)|(np.log10(nh2_array)>23)
# # p_array_copy[Selector] = np.nan
# # s_array_copy[Selector] = np.nan


# plt.figure(figsize=(16,8))
# plt.scatter(I_array_copy,p_array*s_array)

# plotting_equal_width_bins(I_array_copy,(p_array*s_array),eI_array,error,10)
# param,x,y  = fitting_binned(I_array_copy,(p_array*s_array),eI_array,error,10)
# print(param[1])

# # param,x,y = Functions.binning_equal_width((nh2_array),(p_array*s_array),np.ones_like(nh2_array),(p_array*es_array + ep_array*s_array),15)
# # param,x,y = Functions.binning_equal_width(I_array_copy,(p_array*s_array),eI_array,error,10)

# # valid = ~(np.isnan(nh2_array)|np.isnan(p_array*s_array)|np.isnan(np.ones_like(nh2_array))|np.isnan(p_array*es_array + ep_array*s_array))
# # valid = ~(np.isnan(I_array_copy)|np.isnan(p_array*s_array)|np.isnan(eI_array)|np.isnan(error))

# # data = ((nh2_array)[valid],(p_array*s_array)[valid],np.ones_like(nh2_array)[valid],(p_array*es_array + ep_array*s_array)[valid])
# # data = (I_array_copy[valid],(p_array*s_array)[valid],eI_array[valid],(error)[valid])

# # label_temp = r'log(p{frac}*S) = C + $\alpha$log($I$){linebreak}$\alpha$: {alpha_nh2:.3f}{pm}{error:.3f}'.format(frac = r'$\%$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = 'alpha = {alpha:.3f}'.format(alpha = param[1]))
# # plt.plot(x,np.ones_like(x)*np.nanmean(p_array_copy*s_array_copy),'tab:red',linewidth = 2.5,label = "mean : {mean:.3f}".format(mean = np.nanmean(p_array_copy*s_array_copy)))
# plt.title(r'$p\%*S$ vs $I$')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('I')
# plt.ylabel(r'$p\%*S$')
# plt.tight_layout()
# plt.legend()



#%%


plt.figure(figsize=(6,8))
# plt.imshow(CheckMapTemperature.data,origin='lower')
im =plt.imshow(MaskOF_2025.data ,origin='lower',cmap = 'RdBu')
plt.colorbar(im)
plt.title("Mask")
plt.tight_layout()
plt.show()
#%%

s_C3_OF = S_map_deb.copy()
es_C3_OF = sigma_S_map.copy()
p_C3_OF = BlankedMapPol_2025.copy()
ep_C3_OF = BlankedMapDebPolError_2025.copy()
I_C3_OF = BlankedMapStokesI_2025.copy()
eI_C3_OF = BlankedMapStokesIError_2025.copy()
nh2_C3_OF = BlankedMapColumnDensity_2025.copy()
micron8_C3_OF = BlankedMap8Mircon_2025.copy()
temp_C3_OF = BlankedMapTemperature_2025.copy()

s_C3_rest = S_map_deb.copy()
es_C3_rest = sigma_S_map.copy()
p_C3_rest = BlankedMapPol_2025.copy()
ep_C3_rest = BlankedMapDebPolError_2025.copy()
I_C3_rest = BlankedMapStokesI_2025.copy()
eI_C3_rest = BlankedMapStokesIError_2025.copy()
nh2_C3_rest = BlankedMapColumnDensity_2025.copy()
micron8_C3_rest = BlankedMap8Mircon_2025.copy()
temp_C3_rest = BlankedMapTemperature_2025.copy()


Selector = (MaskOF_2025.data != 7)

s_C3_OF.data[Selector] = np.nan
es_C3_OF.data[Selector] = np.nan
p_C3_OF.data[Selector] = np.nan
ep_C3_OF.data[Selector] = np.nan
I_C3_OF.data[Selector] = np.nan
eI_C3_OF.data[Selector] = np.nan
nh2_C3_OF.data[Selector] = np.nan
micron8_C3_OF.data[Selector] = np.nan
temp_C3_OF.data[Selector] = np.nan

Selector = ~Selector
s_C3_rest.data[Selector] = np.nan
es_C3_rest.data[Selector] = np.nan
p_C3_rest.data[Selector] = np.nan
ep_C3_rest.data[Selector] = np.nan
I_C3_rest.data[Selector] = np.nan
eI_C3_rest.data[Selector] = np.nan
nh2_C3_rest.data[Selector] = np.nan
micron8_C3_rest.data[Selector] = np.nan
temp_C3_rest.data[Selector] = np.nan


s_C3_OF_array = s_C3_OF.data.flatten()[::4]
es_C3_OF_array = es_C3_OF.data.flatten()[::4]
p_C3_OF_array = p_C3_OF.data.flatten()[::4]
ep_C3_OF_array = ep_C3_OF.data.flatten()[::4]
I_C3_OF_array = I_C3_OF.data.flatten()[::4]
eI_C3_OF_array = eI_C3_OF.data.flatten()[::4]
nh2_C3_OF_array = nh2_C3_OF.data.flatten()[::4]
micron8_C3_OF_array = micron8_C3_OF.data.flatten()[::4]
temp_C3_OF_array = temp_C3_OF.data.flatten()[::4]

valid = np.isnan(s_C3_OF_array)|np.isnan(es_C3_OF_array)|np.isnan(p_C3_OF_array)|np.isnan(ep_C3_OF_array)|np.isnan(I_C3_OF_array)|np.isnan(eI_C3_OF_array)|np.isnan(nh2_C3_OF_array)|np.isnan(micron8_C3_OF_array)|np.isnan(temp_C3_OF_array)			
s_C3_OF_array = s_C3_OF_array[~valid]
es_C3_OF_array = es_C3_OF_array[~valid]		
p_C3_OF_array = p_C3_OF_array[~valid]
ep_C3_OF_array = ep_C3_OF_array[~valid]
I_C3_OF_array = I_C3_OF_array[~valid]
eI_C3_OF_array = eI_C3_OF_array[~valid]
nh2_C3_OF_array = nh2_C3_OF_array[~valid]
micron8_C3_OF_array = micron8_C3_OF_array[~valid]
temp_C3_OF_array = temp_C3_OF_array[~valid]


s_C3_rest_array = s_C3_rest.data.flatten()[::4]
es_C3_rest_array = es_C3_rest.data.flatten()[::4]
p_C3_rest_array = p_C3_rest.data.flatten()[::4]
ep_C3_rest_array = ep_C3_rest.data.flatten()[::4]
I_C3_rest_array = I_C3_rest.data.flatten()[::4]
eI_C3_rest_array = eI_C3_rest.data.flatten()[::4]
nh2_C3_rest_array = nh2_C3_rest.data.flatten()[::4]
micron8_C3_rest_array = micron8_C3_rest.data.flatten()[::4]
temp_C3_rest_array = temp_C3_rest.data.flatten()[::4]

valid = np.isnan(s_C3_rest_array)|np.isnan(es_C3_rest_array)|np.isnan(p_C3_rest_array)|np.isnan(ep_C3_rest_array)|np.isnan(I_C3_rest_array)|np.isnan(eI_C3_rest_array)|np.isnan(nh2_C3_rest_array)|np.isnan(micron8_C3_rest_array)|np.isnan(temp_C3_rest_array)
s_C3_rest_array = s_C3_rest_array[~valid]
es_C3_rest_array = es_C3_rest_array[~valid]
p_C3_rest_array = p_C3_rest_array[~valid]
ep_C3_rest_array = ep_C3_rest_array[~valid]
I_C3_rest_array = I_C3_rest_array[~valid]
eI_C3_rest_array = eI_C3_rest_array[~valid]
nh2_C3_rest_array = nh2_C3_rest_array[~valid]
micron8_C3_rest_array = micron8_C3_rest_array[~valid]
temp_C3_rest_array = temp_C3_rest_array[~valid]

# fig, ax = plt.subplots(3, 1, figsize=(12, 12))

# # First plot
# ax[0].imshow(I_C3_R1.data,origin = 'lower')
# ax[0].set_title('Ridge')
# # ax[0].legend()
# # ax[0].grid()

# # Second plot
# ax[1].imshow(I_C3_OF.data,origin = 'lower')
# ax[1].set_title('outflow')
# # ax[1].legend()
# # ax[1].grid()

# ax[2].imshow(I_C3_rest.data,origin = 'lower')
# ax[2].set_title('rest')

# # Adjust layout
# plt.tight_layout()

# # Show the plots
# plt.show()


fig, ax = plt.subplots(1, 2, figsize=(12, 12))
ax[0].imshow(es_C3_OF.data,origin = 'lower')
ax[0].set_title('outflow')
ax[0].grid()
ax[1].imshow(es_C3_rest.data,origin = 'lower')
ax[1].set_title('rest')
ax[1].grid()
plt.tight_layout()
plt.close()
# plt.show()
#%%
print('plotting I_vs_N_notOF')
array1 = nh2_C3_rest_array
array2 = I_C3_rest_array
error_array1 = np.ones_like(nh2_C3_rest_array)
error_array2 = eI_C3_rest_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,label = r'Not_OF{linebreak}log(I) = C + $\alpha$log($N$)'.format(linebreak='\n',nh2 = r'$_{N}$'))
plotting_equal_width_bins(array1, array2,error_array1,error_array2,10)
param,x,y  = fitting_binned(array1,array2,error_array1,error_array2,10)
data = (array1,array2,error_array1,error_array2)
param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
label_temp = r'$\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(alpha = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label =label_temp)
plt.title(r'I vs N')
plt.xscale('log')
plt.yscale('log')
# plt.yticks([1e0,1e1])
plt.xticks(np.array([1e22,1e23,1e24]))
plt.xlabel('N')
plt.ylabel(r'I')
plt.legend()
plt.tight_layout()
plt.savefig('../plots/2025/19/I_vs_N_notOF.png')
plt.close()
# plt.show()



# array3 =nh2_C3_rest_array
# weights3 = np.ones_like(nh2_C3_rest_array)

# array4 =I_C3_rest_array
# weights4 = eI_C3_rest_array

# array5 = nh2_C3_OF_array
# weights5 = np.ones_like(nh2_C3_OF_array)

# array6 = I_C3_OF_array
# weights6 = eI_C3_OF_array


# plt.figure(figsize=(12,8))
# plt.scatter(array3,array4,label = r'Not_OF{linebreak}log(I) = C + $\alpha$log($N$)'.format(linebreak='\n',nh2 = r'$_{N}$'))
# # plt.scatter(array5,array6,c = '#7fbc41',label = r'OF{linebreak}log(p{frac}) = C + $\alpha$log($N$)'.format(frac = r'$\%$',linebreak='\n',nh2 = r'$_{N}$'))

# plotting_equal_width_bins(array3,array4,weights3,weights4,10)
# param,x,y  = fitting_binned(array3,array4,weights3,weights4,10)
# print(param[1])
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = 'alpha = {alpha:.3f}'.format(alpha = param[1]))



# param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,15)
# valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
# data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
# res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=1000,paired = True,vectorized=False, random_state=0)
# label_temp = r'$\alpha$: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R = r'$_{Not_OF}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'#abdda4',linewidth = 3.5)

# param,x,y = Functions.binning_equal_width(array5,array6,weights5,weights6,10)
# valid = ~(np.isnan(array5)|np.isnan(array6)|np.isnan(weights5)|np.isnan(weights6))
# data = (array5[valid],array6[valid],weights5[valid],weights6[valid])
# res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=1000,paired = True,vectorized=False, random_state=0)
# label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R = r'$_{OF}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'#ffff99',linewidth = 3.5,label = label_temp)

# plt.title(r'I vs N')
# plt.xscale('log')
# plt.yscale('log')
# # plt.yticks([1e0,1e1])
# plt.xticks(np.array([1e22,1e23,1e24]))
# plt.xlabel('N')
# plt.ylabel(r'I')
# plt.legend()
# plt.tight_layout()
# plt.show()

#%%
print('plotting p_vs_I_notOF')
array1 = I_C3_rest_array
array2 = p_C3_rest_array
error_array1 = eI_C3_rest_array
error_array2 =ep_C3_rest_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,label = r'Not_OF{linebreak}log(p{frac}) = C + $\alpha$log(I)'.format(frac = r'$\%$',linebreak='\n'))
plotting_equal_width_bins(array1, array2,error_array1,error_array2,10)
param,x,y  = fitting_binned(array1,array2,error_array1,error_array2,10)
data = (array1,array2,error_array1,error_array2)
param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
label_temp = r'$\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(alpha = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label =label_temp)
plt.title(r'$p\%$ vs I')
plt.xscale('log')
plt.yscale('log')
plt.yticks([1e0,1e1])
plt.xlabel('I')
plt.ylabel(r'$p\%$')
plt.legend()
plt.tight_layout()
plt.savefig('../plots/2025/19/p_vs_I_notOF.png')
plt.close()
# plt.show()



# array3 =I_C3_rest_array
# weights3 = eI_C3_rest_array

# array4 = p_C3_rest_array
# weights4 = ep_C3_rest_array

# array5 = I_C3_OF_array
# weights5 = eI_C3_OF_array

# array6 = p_C3_OF_array
# weights6 = ep_C3_OF_array

# plt.figure(figsize=(12,8))
# plt.scatter(array3,array4,label = r'Not_OF{linebreak}log(p{frac}) = C + $\alpha$log(I)'.format(frac = r'$\%$',linebreak='\n'))

# plotting_equal_width_bins(array3,array4,weights3,weights4,10)
# param,x,y  = fitting_binned(array3,array4,weights3,weights4,10)
# print(param[1])
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = 'alpha = {alpha:.3f}'.format(alpha = param[1]))


# # plt.scatter(array5,array6,c = '#7fbc41',label = r'OF{linebreak}log(p{frac}) = C + $\alpha$log(I)'.format(frac = r'$\%$',linebreak='\n'))


# # param,x,y = binning_equal_width(array3,array4,weights3,weights4,15)
# # valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
# # data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
# # res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=1000,paired = True,vectorized=False, random_state=0)
# # label_temp = r'$\alpha$: {alpha_I:.3f}{pm}{error:.3f}'.format(R = r'$_{Not_OF}$',alpha_I = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# # plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

# # param,x,y = binning_equal_width(array5,array6,weights5,weights6,10)
# # valid = ~(np.isnan(array5)|np.isnan(array6)|np.isnan(weights5)|np.isnan(weights6))
# # data = (array5[valid],array6[valid],weights5[valid],weights6[valid])
# # res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=1000,paired = True,vectorized=False, random_state=0)
# # label_temp = r'$\alpha${R}: {alpha_I:.3f}{pm}{error:.3f}'.format(R = r'$_{OF}$',alpha_I = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# # plt.plot(x,y,'#ffff99',linewidth = 3.5,label = label_temp)

# plt.title(r'$p\%$ vs I')
# plt.xscale('log')
# plt.yscale('log')
# plt.yticks([1e0,1e1])
# plt.xlabel('I')
# plt.ylabel(r'$p\%$')
# plt.legend()
# plt.tight_layout()
# plt.show()


# %%


print('plotting p_vs_N_notOF')
array1 = nh2_C3_rest_array
array2 = p_C3_rest_array
error_array1 = np.ones_like(nh2_C3_rest_array)
error_array2 =ep_C3_rest_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,label = r'Not_OF{linebreak}log(p{frac}) = C + $\alpha$log($N$)'.format(frac = r'$\%$',linebreak='\n',nh2 = r'$_{N}$'))
plotting_equal_width_bins(array1, array2,error_array1,error_array2,10)
param,x,y  = fitting_binned(array1,array2,error_array1,error_array2,10)
data = (array1,array2,error_array1,error_array2)
param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
label_temp = r'$\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(alpha = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label =label_temp)
plt.title(r'$p\%$ vs N')
plt.xscale('log')
plt.yscale('log')
plt.yticks([1e0,1e1])
plt.xticks(np.array([1e22,1e23,1e24]))
plt.xlabel('N')
plt.ylabel(r'$p\%$')
plt.legend()
plt.tight_layout()
plt.savefig('../plots/2025/19/p_vs_N_notOF.png')
plt.close()
# plt.show()



# array3 =nh2_C3_rest_array
# weights3 = np.ones_like(nh2_C3_rest_array)

# array4 = p_C3_rest_array
# weights4 = ep_C3_rest_array

# array5 = nh2_C3_OF_array
# weights5 = np.ones_like(nh2_C3_OF_array)

# array6 = p_C3_OF_array
# weights6 = ep_C3_OF_array

# plt.figure(figsize=(12,8))
# plt.scatter(array3,array4,label = r'Not_OF{linebreak}log(p{frac}) = C + $\alpha$log($N$)'.format(frac = r'$\%$',linebreak='\n',nh2 = r'$_{N}$'))

# plotting_equal_width_bins(array3,array4,weights3,weights4,10)
# param,x,y  = fitting_binned(array3,array4,weights3,weights4,10)
# print(param[1])
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = 'alpha = {alpha:.3f}'.format(alpha = param[1]))


# plt.scatter(array5,array6,c = '#7fbc41',label = r'OF{linebreak}log(p{frac}) = C + $\alpha$log($N$)'.format(frac = r'$\%$',linebreak='\n',nh2 = r'$_{N}$'))

# param,x,y = binning_equal_width(array3,array4,weights3,weights4,15)
# valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
# data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
# res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=1000,paired = True,vectorized=False, random_state=0)
# label_temp = r'$\alpha$: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R = r'$_{Not_OF}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

# param,x,y = binning_equal_width(array5,array6,weights5,weights6,10)
# valid = ~(np.isnan(array5)|np.isnan(array6)|np.isnan(weights5)|np.isnan(weights6))
# data = (array5[valid],array6[valid],weights5[valid],weights6[valid])
# res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=1000,paired = True,vectorized=False, random_state=0)
# label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R = r'$_{OF}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'#ffff99',linewidth = 3.5,label = label_temp)

# plt.title(r'$p\%$ vs N')
# plt.xscale('log')
# plt.yscale('log')
# plt.yticks([1e0,1e1])
# plt.xticks(np.array([1e22,1e23,1e24]))
# plt.xlabel('N')
# plt.ylabel(r'$p\%$')
# plt.legend()
# plt.tight_layout()
# plt.show()


#%%
print('plotting S_vs_N_notOF')

array1 = nh2_C3_rest_array
array2 = s_C3_rest_array
error_array1 = np.ones_like(nh2_C3_rest_array)
error_array2 = es_C3_rest_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,label = r'NotOF{linebreak}log(S) = C + $\alpha$log($N$)'.format(frac = r'$\%$',linebreak='\n',nh2 = r'$_{N}$'))
plotting_equal_width_bins(array1, array2,error_array1,error_array2,10)
param,x,y  = fitting_binned(array1,array2,error_array1,error_array2,10)
data = (array1,array2,error_array1,error_array2)
param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
label_temp = r'$\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(alpha = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label =label_temp)
plt.title(r'S vs N')
plt.xscale('log')
plt.yscale('log')
plt.yticks([1e0,1e1])
plt.xticks(np.array([1e22,1e23,1e24]))
plt.xlabel('N')
plt.ylabel(r'S')
plt.legend()
plt.tight_layout()
plt.savefig('../plots/2025/19/S_vs_N_notOF.png')
plt.close()
# plt.show()


# array3 =nh2_C3_rest_array
# weights3 = np.ones_like(nh2_C3_rest_array)

# array4 = s_C3_rest_array
# weights4 = es_C3_rest_array

# array5 = nh2_C3_OF_array
# weights5 = np.ones_like(nh2_C3_OF_array)

# array6 = s_C3_OF_array
# weights6 = es_C3_OF_array

# plt.figure(figsize=(12,8))
# plt.scatter(array3,array4,label = r'NotOF{linebreak}log(S) = C + $\alpha$log($N$)'.format(frac = r'$\%$',linebreak='\n',nh2 = r'$_{N}$'))
# # plt.scatter(array5,array6,c = '#7fbc41',label = r'OF{linebreak}log(S) = C + $\alpha$log($N$)'.format(frac = r'$\%$',linebreak='\n',nh2 = r'$_{N}$'))

# plotting_equal_width_bins(array3,array4,weights3,weights4,10)
# param,x,y  = fitting_binned(array3,array4,weights3,weights4,10)
# print(param[1])
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = 'alpha = {alpha:.3f}'.format(alpha = param[1]))




# param,x,y = binning_equal_width(array3,array4,weights3,weights4,15)
# valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
# data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
# res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=1000,paired = True,vectorized=False, random_state=0)
# label_temp = r'$\alpha$: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R = r'$_{Not_OF}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

# param,x,y = binning_equal_width(array5,array6,weights5,weights6,10)
# valid = ~(np.isnan(array5)|np.isnan(array6)|np.isnan(weights5)|np.isnan(weights6))
# data = (array5[valid],array6[valid],weights5[valid],weights6[valid])
# res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=1000,paired = True,vectorized=False, random_state=0)
# label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R = r'$_{OF}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'#ffff99',linewidth = 3.5,label = label_temp)

# plt.title(r'S vs N')
# plt.xscale('log')
# plt.yscale('log')
# plt.yticks([1e0,1e1])
# plt.xticks(np.array([1e22,1e23,1e24]))
# plt.xlabel('N')
# plt.ylabel(r'S')
# plt.legend()
# plt.tight_layout()
# plt.show()

#%%
print('plotting S_vs_I_notOF')
array1 = I_C3_rest_array
array2 = s_C3_rest_array
error_array1 = eI_C3_rest_array
error_array2 = es_C3_rest_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,label = r'Not_OF{linebreak}log(S) = C + $\alpha$log($I$)'.format(frac = r'$\%$',linebreak='\n',nh2 = r'$_{N}$'))
plotting_equal_width_bins(array1, array2,error_array1,error_array2,10)
param,x,y  = fitting_binned(array1,array2,error_array1,error_array2,10)
data = (array1,array2,error_array1,error_array2)
param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
label_temp = r'$\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(alpha = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label =label_temp)
plt.title(r'S vs I')
plt.xscale('log')
plt.yscale('log')
plt.yticks([1e0,1e1])
plt.xlabel('I')
plt.ylabel(r'S')
plt.legend()
plt.tight_layout()
plt.savefig('../plots/2025/19/S_vs_I_notOF.png')
plt.close()
# plt.show()



# array3 =I_C3_rest_array
# weights3 = eI_C3_rest_array

# array4 = s_C3_rest_array
# weights4 = es_C3_rest_array

# array5 = I_C3_OF_array
# weights5 = eI_C3_OF_array

# array6 = s_C3_OF_array
# weights6 = es_C3_OF_array

# plt.figure(figsize=(12,8))
# plt.scatter(array3,array4,label = r'Not_OF{linebreak}log(S) = C + $\alpha$log($I$)'.format(frac = r'$\%$',linebreak='\n',nh2 = r'$_{N}$'))
# # plt.scatter(array5,array6,c = '#7fbc41',label = r'OF{linebreak}log(S) = C + $\alpha$log($N$)'.format(frac = r'$\%$',linebreak='\n',nh2 = r'$_{N}$'))

# plotting_equal_width_bins(array3,array4,weights3,weights4,10)
# param,x,y  = fitting_binned(array3,array4,weights3,weights4,10)
# print(param[1])
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = 'alpha = {alpha:.3f}'.format(alpha = param[1]))




# param,x,y = binning_equal_width(array3,array4,weights3,weights4,15)
# valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
# data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
# res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=1000,paired = True,vectorized=False, random_state=0)
# label_temp = r'$\alpha$: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R = r'$_{Not_OF}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

# param,x,y = binning_equal_width(array5,array6,weights5,weights6,10)
# valid = ~(np.isnan(array5)|np.isnan(array6)|np.isnan(weights5)|np.isnan(weights6))
# data = (array5[valid],array6[valid],weights5[valid],weights6[valid])
# res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=1000,paired = True,vectorized=False, random_state=0)
# label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R = r'$_{OF}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'#ffff99',linewidth = 3.5,label = label_temp)

# plt.title(r'S vs I')
# plt.xscale('log')
# plt.yscale('log')
# plt.yticks([1e0,1e1])
# # plt.xticks(np.array([1e22,1e23,1e24]))
# plt.xlabel('I')
# plt.ylabel(r'S')
# plt.legend()
# plt.tight_layout()
# plt.show()



#%%
print('plotting S_vs_p_notOF')

array1 =  p_C3_rest_array
array2 = s_C3_rest_array
error_array1 = ep_C3_rest_array
error_array2 = es_C3_rest_array
plt.figure(figsize=(12,8))
plt.scatter(array1,array2,label = r'Not_OF{linebreak}log(S) = C + $\alpha$log(p{frac})'.format(frac = r'$\%$',linebreak='\n'))
plotting_equal_width_bins(array1, array2,error_array1,error_array2,10)
param,x,y  = fitting_binned(array1,array2,error_array1,error_array2,10)
data = (array1,array2,error_array1,error_array2)
param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
label_temp = r'$\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(alpha = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label =label_temp)
plt.title(r'S vs $p\%$ ')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('S')
plt.xlabel(r'$p\%$')
plt.tight_layout()
plt.legend()
plt.savefig('../plots/2025/19/S_vs_p_notOF.png')
plt.close()
# plt.show()


# array3 = p_C3_rest_array
# weights3 = ep_C3_rest_array


# array4 =s_C3_rest_array
# weights4 = es_C3_rest_array

# array5 = p_C3_OF_array
# weights5 = ep_C3_OF_array

# array6 = s_C3_OF_array
# weights6 = es_C3_OF_array



# plt.figure(figsize=(12,8))
# plt.scatter(array3,array4,label = r'Not_OF{linebreak}log(S) = C + $\alpha$log(p{frac})'.format(frac = r'$\%$',linebreak='\n'))

# plotting_equal_width_bins(array3,array4,weights3,weights4,10)
# param,x,y  = fitting_binned(array3,array4,weights3,weights4,10)
# print(param[1])
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = 'alpha = {alpha:.3f}'.format(alpha = param[1]))



# plt.scatter(array5,array6,c = '#7fbc41',label = r'OF{linebreak}log(S) = C + $\alpha$log(p{frac})'.format(frac = r'$\%$',linebreak='\n'))

# param,x,y = binning_equal_width(array3,array4,weights3,weights4,10)
# x,y,z = binning_equal_width_latest(array3,array4,weights3,weights4,10)
# print(fitting_binned(x,y,z))
# valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
# data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
# res = bootstrapscipy(data, bootstrap_binned,n_resamples=10,paired = True,vectorized=True, random_state=0)
# label_temp = r'$\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(R = r'$_{Not_OF}}$',alpha = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)



# param,x,y = binning_equal_width(array5,array6,weights5,weights6,10)
# x,y,z = binning_equal_width_latest(array5,array6,weights5,weights6,10)
# # print(fitting_binned(x,y,z))

# valid = ~(np.isnan(array5)|np.isnan(array6)|np.isnan(weights5)|np.isnan(weights6))
# data = (array5[valid],array6[valid],weights5[valid],weights6[valid])
# res = bootstrapscipy(data, binning_equal_width_vbootstrap,n_resamples=10,paired = True,vectorized=False, random_state=0)
# label_temp = r'$\alpha${R}: {alpha:.3f}{pm}{error:.3f}'.format(R = r'$_{OF}$',alpha = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'#ffff99',linewidth = 3.5)

# plt.title(r'S vs $p\%$ ')
# plt.xscale('log')
# plt.yscale('log')
# plt.ylabel('S')
# plt.xlabel(r'$p\%$')
# plt.tight_layout()
# plt.legend()
# plt.show()



# %%
print('plotting pS_vs_I_notOF')
array1 =  I_C3_rest_array
array2 =  p_C3_rest_array*s_C3_rest_array
error_array1 = eI_C3_rest_array
error_array2 = np.sqrt((ep_C3_rest_array*s_C3_rest_array)**2 + (es_C3_rest_array*p_C3_rest_array)**2)
plt.figure(figsize=(12,8))
plt.scatter(array1,array2,label = r'Not_OF{linebreak}log(S) = C + $\alpha$log(p{frac})'.format(frac = r'$\%$',linebreak='\n'))
plotting_equal_width_bins(array1, array2,error_array1,error_array2,10)
param,x,y  = fitting_binned(array1,array2,error_array1,error_array2,10)
data = (array1,array2,error_array1,error_array2)
param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
label_temp = r'$\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(alpha = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label =label_temp)
plt.title(r'$p\%*S$ vs $I$')
plt.xlabel('I')
plt.ylabel(r'$p\%*S$')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.legend()
plt.savefig('../plots/2025/19/pS_vs_I_notOF.png')
plt.close()
# plt.show()


# array3 =I_C3_rest_array
# weights3 = eI_C3_rest_array
# array4 = p_C3_rest_array*s_C3_rest_array
# weights4 = np.sqrt((ep_C3_rest_array*s_C3_rest_array)**2 + (es_C3_rest_array*p_C3_rest_array)**2)
# array5 = I_C3_OF_array
# weights5 = eI_C3_OF_array
# array6 = p_C3_OF_array*s_C3_OF_array
# weights6 =  np.sqrt((ep_C3_OF_array*s_C3_OF_array)**2 + (es_C3_OF_array*p_C3_OF_array)**2)

# plt.figure(figsize=(12,8))
# plt.scatter(array3,array4,label = r'Not_OF{linebreak}log(p{frac}*S) = C + $\alpha$log($I$)'.format(frac = r'$\%$',linebreak='\n'))
# # plt.scatter(array5,array6,c = '#7fbc41',label = r'OF{linebreak}log(p{frac}*S) = C + $\alpha$log($N$)'.format(frac = r'$\%$',linebreak='\n'))

# plotting_equal_width_bins(array3,array4,weights3,weights4,10)
# param,x,y  = fitting_binned(array3,array4,weights3,weights4,10)
# print(param[1])
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = 'alpha = {alpha:.3f}'.format(alpha = param[1]))




# param,x,y = binning_equal_width(array3,array4,weights3,weights4,15)
# valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
# data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
# res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=1000,paired = True,vectorized=False, random_state=0)
# label_temp = r'$\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(R = r'$_{rest}$',alpha = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

# param,x,y = binning_equal_width(array5,array6,weights5,weights6,10)
# valid = ~(np.isnan(array5)|np.isnan(array6)|np.isnan(weights5)|np.isnan(weights6))
# data = (array5[valid],array6[valid],weights5[valid],weights6[valid])
# res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=1000,paired = True,vectorized=False, random_state=0)
# label_temp = r'$\alpha${R}: {alpha:.3f}{pm}{error:.3f}'.format(R = r'$_{OF}$',alpha = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'#ffff99',linewidth = 3.5,label = label_temp)
# plt.title(r'$p\%*S$ vs $I$')
# plt.xlabel('I')
# plt.ylabel(r'$p\%*S$')
# plt.xscale('log')
# plt.yscale('log')
# plt.tight_layout()
# plt.legend()
# plt.show()
#%%
print('plotting pS_vs_N_notOF')

array1 =  nh2_C3_rest_array
array2 =  p_C3_rest_array*s_C3_rest_array
error_array1 = np.ones_like(nh2_C3_rest_array)
error_array2 = np.sqrt((ep_C3_rest_array*s_C3_rest_array)**2 + (es_C3_rest_array*p_C3_rest_array)**2)
plt.figure(figsize=(12,8))
plt.scatter(array1,array2,label = r'Not_OF{linebreak}log(p{frac}*S) = C + $\alpha$log($N$)'.format(frac = r'$\%$',linebreak='\n'))
plotting_equal_width_bins(array1, array2,error_array1,error_array2,10)
param,x,y  = fitting_binned(array1,array2,error_array1,error_array2,10)
data = (array1,array2,error_array1,error_array2)
param_mean,param_std  = bootstrap_fitting(data[0],data[1],data[2],data[3],10)
label_temp = r'$\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(alpha = param_mean,pm = r'$\pm$',error=param_std,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label =label_temp)
plt.title(r'$p\%*S$ vs $N$')
plt.xlabel('N')
plt.ylabel(r'$p\%*S$')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.legend()
plt.savefig('../plots/2025/19/pS_vs_N_notOF.png')
plt.close()
# plt.show()



# array3 =nh2_C3_rest_array
# weights3 = np.ones_like(nh2_C3_rest_array)


# array4 = p_C3_rest_array*s_C3_rest_array
# weights4 = np.sqrt((ep_C3_rest_array*s_C3_rest_array)**2 + (es_C3_rest_array*p_C3_rest_array)**2)

# array5 = nh2_C3_OF_array
# weights5 = np.ones_like(nh2_C3_OF_array)

# array6 = p_C3_OF_array*s_C3_OF_array
# weights6 =  np.sqrt((ep_C3_OF_array*s_C3_OF_array)**2 + (es_C3_OF_array*p_C3_OF_array)**2)

# plt.figure(figsize=(12,8))
# plt.scatter(array3,array4,label = r'Not_OF{linebreak}log(p{frac}*S) = C + $\alpha$log($N$)'.format(frac = r'$\%$',linebreak='\n'))
# # plt.scatter(array5,array6,c = '#7fbc41',label = r'OF{linebreak}log(p{frac}*S) = C + $\alpha$log($N$)'.format(frac = r'$\%$',linebreak='\n'))
# plotting_equal_width_bins(array3,array4,weights3,weights4,10)
# param,x,y  = fitting_binned(array3,array4,weights3,weights4,10)
# print(param[1])
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = 'alpha = {alpha:.3f}'.format(alpha = param[1]))



# param,x,y = binning_equal_width(array3,array4,weights3,weights4,15)
# valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
# data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
# res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=1000,paired = True,vectorized=False, random_state=0)
# label_temp = r'$\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(R = r'$_{rest}$',alpha = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

# param,x,y = binning_equal_width(array5,array6,weights5,weights6,10)
# valid = ~(np.isnan(array5)|np.isnan(array6)|np.isnan(weights5)|np.isnan(weights6))
# data = (array5[valid],array6[valid],weights5[valid],weights6[valid])
# res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=1000,paired = True,vectorized=False, random_state=0)
# label_temp = r'$\alpha${R}: {alpha:.3f}{pm}{error:.3f}'.format(R = r'$_{OF}$',alpha = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'#ffff99',linewidth = 3.5,label = label_temp)
# plt.title(r'$p\%*S$ vs $N$')
# plt.xlabel('N')
# plt.ylabel(r'$p\%*S$')
# plt.xscale('log')
# plt.yscale('log')
# plt.tight_layout()
# plt.legend()
# plt.show()
# %%
