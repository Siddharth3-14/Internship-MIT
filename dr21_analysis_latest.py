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
# %matplotlib inline

#%%
########## importing and testing the file
FITS1 = '../FITS_file/new_fits/DR21_OTF_full_pipeline.fits'
FITS2 = '../FITS_file/new_fits/DR21_full_NH2_Repr.fits'
FITS3 = '../FITS_file/new_fits/DR21_full_Tdust_Repr.fits'
FITS4 = '../FITS_file/new_fits/DR21_full_IRAC4_Repr.fits'
FITS5 = '../FITS_file/new_fits/DR21_full_Her250_Repr.fits'
FITS6 = '../FITS_file/new_fits/DR21_full_Fil_Mask.fits'

hdul = fits.open(FITS1)
hdul2 = fits.open(FITS2)
hdul3 = fits.open(FITS3)
hdul4 = fits.open(FITS4)
hdul5 = fits.open(FITS5)
hdul6 = fits.open(FITS6)

print(hdul.info())

MapStokesI = hdul[0]
MapStokesIError = hdul[1]
MapStokesQ = hdul[2]
MapStokesU = hdul[4]
MapDebPol = hdul[8]
MapDebPolError = hdul[9]
MapPolAngleNonRotated = hdul[10]
MapPolAngle = hdul[11]
MapPolAngleError = hdul[12]
MapPolFlux = hdul[13]
MapPolFluxError = hdul[14]
MapColumndensity = hdul2[0]
MapTemperature = hdul3[0]
Map8Micron = hdul4[0]
MapHer250 = hdul5[0]
Mask = hdul6[0]

MapPolSNR = MapDebPol.copy()
BlankedMapPol = MapDebPol.copy()
BlankedMapPolAngle = MapPolAngle.copy()
BlankedMapPolAngleError = MapPolAngleError.copy()
BlankedMapPolAngleNonRotated = MapPolAngleNonRotated.copy() 
BlankedMapStokesI = MapStokesI.copy()
BlankedMapStokesIError = MapStokesIError.copy()
BlankedMapStokesQ = MapStokesQ.copy()
BlankedMapStokesU = MapStokesU.copy()
BlankedMapColumnDensity = MapColumndensity.copy()
BlankedMapTemperature = MapTemperature.copy()
BlankedMap8Mircon = Map8Micron.copy()
BlankedMapHer250 = MapHer250.copy()
BlankedMapDebPolError = MapDebPolError.copy()


######## taking points only with singal to noise ratio more than 3
MapPolSNR.data[:] = np.nan
MapPolSNR.data = MapDebPol.data/MapDebPolError.data


Selector = (MapPolSNR.data < 3)
BlankedMapPol.data[Selector] = np.nan
BlankedMapPolAngle.data[Selector] = np.nan
BlankedMapPolAngleError.data[Selector] = np.nan
BlankedMapStokesI.data[Selector] = np.nan
BlankedMapStokesIError.data[Selector] = np.nan
BlankedMapStokesQ.data[Selector] = np.nan
BlankedMapStokesU.data[Selector] = np.nan
BlankedMapPolAngleNonRotated.data[Selector] = np.nan
BlankedMapColumnDensity.data[Selector] = np.nan
BlankedMapTemperature.data[Selector] = np.nan
BlankedMap8Mircon.data[Selector] = np.nan
BlankedMapHer250.data[Selector] = np.nan
BlankedMapDebPolError.data[Selector] = np.nan

############## removing any points with pfrac above 50
Selector = (BlankedMapPol.data>50)
BlankedMapPol.data[Selector] = np.nan
BlankedMapPolAngle.data[Selector] = np.nan
BlankedMapPolAngleError.data[Selector] = np.nan
BlankedMapStokesI.data[Selector] = np.nan
BlankedMapStokesIError.data[Selector] = np.nan
BlankedMapStokesQ.data[Selector] = np.nan
BlankedMapStokesU.data[Selector] = np.nan
BlankedMapPolAngleNonRotated.data[Selector] = np.nan
BlankedMapColumnDensity.data[Selector] = np.nan
BlankedMapTemperature.data[Selector] = np.nan
BlankedMap8Mircon.data[Selector] = np.nan
BlankedMapHer250.data[Selector] = np.nan
BlankedMapDebPolError.data[Selector] = np.nan


############ removing any data points with I/I_error < 100
Selector = MapStokesI.data/MapStokesIError.data < 100
BlankedMapPol.data[Selector] = np.nan
BlankedMapPolAngle.data[Selector] = np.nan
BlankedMapPolAngleError.data[Selector] = np.nan
BlankedMapStokesI.data[Selector] = np.nan
BlankedMapStokesIError.data[Selector] = np.nan
BlankedMapStokesQ.data[Selector] = np.nan
BlankedMapStokesU.data[Selector] = np.nan
BlankedMapPolAngleNonRotated.data[Selector] = np.nan
BlankedMapColumnDensity.data[Selector] = np.nan
BlankedMapTemperature.data[Selector] = np.nan
BlankedMap8Mircon.data[Selector] = np.nan
BlankedMapHer250.data[Selector] = np.nan
BlankedMapDebPolError.data[Selector] = np.nan

# plt.figure(figsize=(6,6))
# plt.imshow(np.log10(MapHer250.data),origin='lower',vmin = 0 , vmax = 3)

# Selector = BlankedMap8Mircon.data < 80
# BlankedMap8Mircon.data[Selector] = np.nan


# BlankedMapColumnDensity.data = BlankedMapColumnDensity.data*(BlankedMapStokesI.data/BlankedMapStokesI.data)
# BlankedMapTemperature.data = BlankedMapTemperature.data*(BlankedMapStokesI.data/BlankedMapStokesI.data)
# BlankedMap8Mircon.data = BlankedMap8Mircon.data*(BlankedMapStokesI.data/BlankedMapStokesI.data)
# BlankedMapPolAngleError.data = BlankedMapPolAngleError.data*(BlankedMapStokesI.data/BlankedMapStokesI.data)
# BlankedMapHer250.data = BlankedMapHer250.data*(BlankedMapStokesI.data/BlankedMapStokesI.data)

############## generating the RA and DEC mesh
DEC_grid,RA_grid = Functions.generate_RA_DEC_mesh(hdul[0])
seperation = MapPolAngle.copy()

# CheckMapTemperature = BlankedMapTemperature.copy()
# selector = BlankedMapTemperature.data>40
# CheckMapTemperature.data[selector] = 5
# CheckMapTemperature = BlankedMapTemperature.copy()
# selector = BlankedMapTemperature.data<30
# CheckMapTemperature.data[selector] = 5

print(1)
plt.figure(figsize=(6,8))
# plt.imshow(CheckMapTemperature.data,origin='lower')
plt.imshow(BlankedMap8Mircon.data ,origin='lower',vmin=80)
plt.tight_layout()
plt.show()


#%%

def get_data(infile):
	#get data
	hawc = fits.open(infile)
	p    = hawc[8]
	perr = hawc[9]
	#pa   = hawc['ROTATED POL ANGLE']
	pa   = hawc[11]
	stkI = hawc[0]
	stkIerr = hawc[1]
	pi   = hawc[15]
	pierr   = hawc[14]

	#Jy/px to Jy/sqarcsec
	pxscale = stkI.header['CDELT2']*3600
	stkI.data /= pxscale**2
	pi.data /= pxscale**2 
	stkIerr.data /= pxscale**2
	return p,perr,pa,stkI,stkIerr,pi,pierr,pxscale

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


filename1='../FITS_file/new_fits/DR21_OTF_full_pipeline.fits'
filename2= '../FITS_file/new_fits/DR21_full_log_NH2_Repr.fits'

hawcp = fits.open(filename1)
Herschel = fits.open(filename2)
MapHer250 = Herschel[0]

x = hawcp['STOKES I'].data
y = hawcp['DEBIASED PERCENT POL'].data
y1=hawcp['ROTATED POL ANGLE'].data
y2 = hawcp['DEBIASED POL FLUX'].data



Herschel = fits.open(filename2)
MapHer250 = Herschel[0]

title = 'SIMPLIFI'

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


p,perr,pa,stkI,stkIerr,pi,pierr,pxscale = get_data(filename1)
p = quality_cuts(stkI,stkIerr,p,perr,SNRp_cut,p_cut,SNRi_cut)

RA = (stkI.header['OBSRA']*u.hourangle).to(u.deg)
DEC = stkI.header['OBSDEC']*u.deg


#### SCRIPT
fig = plt.figure(figsize=(13,10))
gc = FITSFigure(MapHer250,figure=fig)
gc.show_colorscale(cmap='default',vmin = 22.1,vmax = 24.5)
gc.add_colorbar(location='right', width=0.2, pad=0.25, ticks=None,axis_label_text= 'log Column Density')
gc.show_contour(colors = 'white',levels = 7)
# gc.show_vectors(p,pa,scale=scalevec,step=4,color='black',linewidth=3.5)
# gc.show_vectors(p,pa,scale=scalevec,step=4,color='yellow',linewidth=2.0)
vecscale = scalevec * pxscale/3600
# gc.add_scalebar(vec_legend*vecscale,r'$p_{frac}$ ='+str(vec_legend),corner='bottom right',frame=True,color='black')
plt.show()

# %%

#### testing

# filename2= '../FITS_file/new_fits/DR21_full_log_NH2_Repr.fits'

# Herschel = fits.open(filename2)
# MapHer250 = Herschel[0]

# title = 'SIMPLIFI'

# # figure
# width  = 50
# height = 50
# cmap = 'plasma'

# title_size = 16
# tick_labels = 15
# label_plot = 15
# label_colorbar = 15
# tick_colorbar = 15
# label_fontsize = 20

# SNRi_cut = 100
# scalevec = 0.5 #1px = scalevec * 1% pol 
# vec_legend = 5.0

# #### SCRIPT
# fig = plt.figure(figsize=(13,10))
# gc = FITSFigure(BlankedMapPol,figure=fig)
# gc.show_colorscale(cmap='default',vmin = 0,vmax = 18)
# gc.add_colorbar(location='right', width=0.2, pad=0.15, ticks=None,axis_label_text= 'Polarization Fraction')
# # gc.show_contour(colors = 'white',levels = 7)
# gc.show_circles(RA_array,DEC_array,0.002)
# plt.show()

# %%

###### making the S map

set_delta = 9/60   # in arcminute
S_map = BlankedMapPolAngle.copy()
sigma_S_map = BlankedMapPolAngleError.copy()

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
        tempa = BlankedMapStokesQ.data*BlankedMapStokesU.data[i,j] - BlankedMapStokesQ.data[i,j]*BlankedMapStokesU.data
        tempb = BlankedMapStokesQ.data*BlankedMapStokesQ.data[i,j] + BlankedMapStokesU.data*BlankedMapStokesU.data[i,j]
        AngleDiff_v2 = 0.5 * (180/np.pi)*np.arctan2(tempa,tempb)
        S = np.nanmean(AngleDiff_v2[seperation_selector]**2)**0.5
        S_map.data[i,j] = S

        ##### making the dispersion error map
        # sigma_S = np.nanmean(BlankedMapPolAngleError.data[seperation_selector]**2)**0.5
        # sigma_S_map.data[i,j] = sigma_S
        N = np.nansum(~np.isnan(seperation.data))
        print(N)
        temp1 = (( 1 /(S_map.data[i,j]*N) )**2)
        temp2 = np.nansum((AngleDiff_v2[seperation_selector])**2)*(BlankedMapPolAngleError.data[i,j]**2)
        temp3 =  np.nansum((AngleDiff_v2[seperation_selector]*BlankedMapPolAngleError.data[seperation_selector])**2)

        
        # temp1 = ( ( BlankedMapPolAngleError.data[i,j]/( S_map.data[i,j]*N ) )**2)*(np.nansum(AngleDiff_v2[seperation_selector])**2)
        # temp2 = (BlankedMapPolAngleError.data*AngleDiff_v2)
        # temp3 = (( 1 /(S_map.data[i,j]*N) )**2)*( np.nansum( (temp2[seperation_selector])**2 ) )

        sigma_S =  temp1*(temp2 + temp3)
        sigma_S_map.data[i,j] = np.sqrt(sigma_S)

S_map_deb = S_map.copy()

S_map_deb.data = np.sqrt(S_map.data**2 - sigma_S_map.data**2)

#%%
S_map_deb =  fits.open('../FITS_file/new_fits/S_map_deb_latest.fits')[1]
sigma_S_map = fits.open('../FITS_file/new_fits/sigma_S_latest.fits')[1]

Selector = S_map_deb.data/sigma_S_map.data < 3
S_map_deb.data[Selector] = np.nan
sigma_S_map.data[Selector] = np.nan


# Selector = S_map_deb1.data/sigma_S_map1.data < 5
# S_map_deb1.data[Selector] = np.nan
# sigma_S_map1.data[Selector] = np.nan


fig = plt.figure(figsize=(13,10))
gc1 = FITSFigure(S_map_deb,figure=fig)
gc1.show_colorscale(cmap='default')
gc1.add_colorbar(location='top', width=0.2, pad=0.15, ticks=None,axis_label_text= 'Angle Dispersion')
gc1.show_contour(colors = 'white',levels = 5)
plt.tight_layout()
plt.show()



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
p_array = BlankedMapPol.data.flatten()[::4]
ep_array = BlankedMapDebPolError.data.flatten()[::4]
I_array = BlankedMapStokesI.data.flatten()[::4]
eI_array = BlankedMapStokesIError.data.flatten()[::4]
nh2_array = BlankedMapColumnDensity.data.flatten()[::4]
micron8_array = BlankedMap8Mircon.data.flatten()[::4]
temp_array = BlankedMapTemperature.data.flatten()[::4]

#%%

#### single paramter fitting

plt.figure(figsize=(12,8))
plt.scatter(I_array,p_array)
param,x,y = Functions.binning_equal_width(I_array,p_array,eI_array,ep_array,15)
valid = ~(np.isnan(I_array)|np.isnan(p_array)|np.isnan(eI_array)|np.isnan(ep_array))
data = (I_array[valid],p_array[valid],eI_array[valid],ep_array[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'log(p{frac}) = C + $\alpha$log($I$){linebreak}  $\alpha$: {alpha_I:.3f}{pm}{error:.3f}'.format(frac = r'$_{frac}$',I = r'$_{I}$',alpha_I = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
plt.plot(x,y,'tab:red',linewidth = 2.5,label = label_temp)
plt.title(r'$p_{frac}$ vs $I$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('I')
plt.ylabel(r'$p_{frac}$')
plt.tight_layout()
plt.legend()
plt.show()
#%%

# nh2_bins = np.linspace(np.nanmin(np.log10(nh2_array)),np.nanmax(np.log10(nh2_array)),25)
# p_bins = np.linspace(np.nanmin(np.log10(p_array)),np.nanmax(np.log10(p_array)),25)

plt.figure(figsize=(12,8))
plt.scatter(nh2_array,p_array,c = "#1F77B4")
# plt.hist2d(nh2_array,p_array,bins = [10**nh2_bins,10**p_bins])
NH2_param,x,y = Functions.binning_equal_width(nh2_array,p_array,np.ones_like(nh2_array),ep_array,15)
valid = ~(np.isnan(nh2_array)|np.isnan(p_array)|np.isnan(np.ones_like(nh2_array))|np.isnan(ep_array))
data = (nh2_array[valid],p_array[valid],np.ones_like(nh2_array)[valid],ep_array[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'log(p{frac}) = C + $\alpha$log($N$){linebreak}  $\alpha$: {alpha_nh2:.3f}{pm}{error:.3f}'.format(frac = r'$_{frac}$',nh2 = r'$_{N}$',alpha_nh2 = NH2_param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
plt.plot(x,y,'#D62728',linewidth = 3.5,label = label_temp)
plt.title(r'$p_{frac}$ vs $N$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel(r'$p_{frac}$')
plt.tight_layout()
plt.legend()
plt.show()

#%%
p_array1 = p_array.copy()
s_array1 = s_array.copy()
ep_array1 = ep_array.copy()
es_array1 = es_array.copy()

Selector = p_array < 0.5
p_array1[Selector] = np.nan
s_array1[Selector] = np.nan
ep_array1[Selector] = np.nan
es_array1[Selector] = np.nan

# s_bins = np.linspace(np.nanmin(np.log10(s_array)),np.nanmax(np.log10(s_array)),25)
# p_bins = np.linspace(np.nanmin(np.log10(p_array)),np.nanmax(np.log10(p_array)),25)

plt.figure(figsize=(12,8))
plt.scatter(s_array1,p_array1)
# plt.hist2d(s_array,p_array,bins=[10**s_bins,10**p_bins])
param,x,y = Functions.binning_equal_width(s_array1,p_array1,es_array1,ep_array1,15)
valid = ~(np.isnan(s_array1)|np.isnan(p_array1)|np.isnan(es_array1)|np.isnan(ep_array1))
data = (s_array1[valid],p_array1[valid],es_array1[valid],ep_array1[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'log(p{frac}) = C + $\alpha$log(S){linebreak} $\alpha$: {alpha_s:.3f}{pm}{error:.3f}'.format(frac = r'$_{frac}$',alpha_s = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
plt.plot(x,y,'tab:red',linewidth = 2.5,label = label_temp)
plt.title(r'$p_{frac}$ vs $S$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('S')
plt.ylabel(r'$p_{frac}$')
plt.tight_layout()
plt.legend()
plt.show()

#%%


# s_bins = np.linspace(np.nanmin(np.log10(s_array)),np.nanmax(np.log10(s_array)),25)
# p_bins = np.linspace(np.nanmin(np.log10(p_array)),np.nanmax(np.log10(p_array)),25)

plt.figure(figsize=(12,8))
plt.scatter(p_array,s_array)
# plt.hist2d(p_array,s_array,bins=[10**p_bins,10**s_bins])
param,x,y = Functions.binning_equal_width(p_array,s_array,ep_array,es_array,15)
valid = ~(np.isnan(p_array)|np.isnan(s_array)|np.isnan(ep_array)|np.isnan(es_array))
data = (p_array[valid],s_array[valid],ep_array[valid],es_array[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'log(S) = C + $\alpha$log(p{frac}){linebreak} $\alpha$: {alpha_s:.3f}{pm}{error:.3f}'.format(frac = r'$_{frac}$',alpha_s = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
plt.plot(x,y,'tab:red',linewidth = 2.5,label = label_temp)
plt.title(r'$S$ vs $p_{frac}$')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('S')
plt.xlabel(r'$p_{frac}$')
plt.tight_layout()
plt.legend()
plt.show()


#%%
plt.figure(figsize=(16,8))
plt.scatter(nh2_array,s_array1)
param,x,y = Functions.binning_equal_width((nh2_array),(s_array1),np.ones_like(nh2_array),es_array1,15)
valid = ~(np.isnan(nh2_array)|np.isnan(s_array1)|np.isnan(np.ones_like(nh2_array))|np.isnan(es_array1))
data = ((nh2_array)[valid],(s_array1)[valid],np.ones_like(nh2_array)[valid],(es_array1)[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'log(S) = C + $\alpha$log($N$){linebreak}$\alpha$: {alpha_nh2:.3f}{pm}{error:.3f}'.format(frac = r'$_{frac}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
plt.plot(x,y,'tab:red',linewidth = 2.5,label = label_temp)
plt.title(r'$S$ vs $N$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel(r'$S$')
plt.tight_layout()
plt.legend()

#%%

p_array_copy = p_array.copy()
s_array_copy = s_array.copy()


Selector = (np.log10(nh2_array)<22)|(np.log10(nh2_array)>23)
p_array_copy[Selector] = np.nan
s_array_copy[Selector] = np.nan


plt.figure(figsize=(16,8))
plt.scatter(nh2_array,p_array*s_array)
param,x,y = Functions.binning_equal_width((nh2_array),(p_array*s_array),np.ones_like(nh2_array),(p_array*es_array + ep_array*s_array),15)
# valid = ~(np.isnan(nh2_array)|np.isnan(p_array*s_array)|np.isnan(np.ones_like(nh2_array))|np.isnan(p_array*es_array + ep_array*s_array))
# data = ((nh2_array)[valid],(p_array*s_array)[valid],np.ones_like(nh2_array)[valid],(p_array*es_array + ep_array*s_array)[valid])
# res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
# label_temp = r'log(p{frac}*S) = C + $\alpha$log($N$){linebreak}$\alpha$: {alpha_nh2:.3f}{pm}{error:.3f}'.format(frac = r'$_{frac}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
plt.plot(x,np.ones_like(x)*np.nanmean(p_array_copy*s_array_copy),'tab:red',linewidth = 2.5,label = "mean : {mean:.3f}".format(mean = np.nanmean(p_array_copy*s_array_copy)))
plt.title(r'$p_{frac}*S$ vs $N$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel(r'$p_{frac}*S$')
plt.tight_layout()
plt.legend()
#%%

plt.figure(figsize=(16,8))
plt.scatter(nh2_array,temp_array)
param,x,y = Functions.binning_equal_width((nh2_array),(temp_array),np.ones_like(nh2_array),np.ones_like(temp_array),15)
valid = ~(np.isnan(nh2_array)|np.isnan(p_array)|np.isnan(np.ones_like(nh2_array))|np.isnan(np.ones_like(temp_array)))
data = (nh2_array[valid],temp_array[valid],np.ones_like(nh2_array)[valid],np.ones_like(temp_array)[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'log(T) = C + $\alpha$log($N$){linebreak}$\alpha$: {alpha_nh2:.3f}{pm}{error:.3f}'.format(frac = r'$_{frac}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
plt.plot(x,y,'tab:red',linewidth = 2.5,label = label_temp)
plt.title(r'$T$ vs $N$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel(r'$T$')
plt.tight_layout()
plt.legend()
plt.show()

#%%
plt.figure(figsize=(10,6))
plt.scatter(micron8_array,p_array)
param,x,y = Functions.binning_equal_width(micron8_array,p_array,np.ones_like(micron8_array),ep_array,15)
valid = ~(np.isnan(micron8_array)|np.isnan(p_array)|np.isnan(np.ones_like(micron8_array))|np.isnan(ep_array))
data = (micron8_array[valid],p_array[valid],np.ones_like(micron8_array)[valid],ep_array[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'log(p{frac}) = C + $\alpha$log(8micron){linebreak}$\alpha$: {alpha:.3f}{pm}{error:.3f}'.format(frac = r'$_{frac}$',alpha = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
plt.title('whole region binned Pfrac vs 8micron')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('8micron')
plt.ylabel('pfrac')
plt.legend()
plt.tight_layout()
plt.show()

#%%
plt.figure(figsize=(12,8))
plt.scatter(temp_array,p_array)
param,x,y = Functions.binning_equal_width(temp_array,p_array,np.ones_like(temp_array),ep_array,15)
valid = ~(np.isnan(temp_array)|np.isnan(p_array)|np.isnan(np.ones_like(temp_array))|np.isnan(ep_array))
data = (temp_array[valid],p_array[valid],np.ones_like(temp_array)[valid],ep_array[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'log(p{frac}) = C + $\alpha$log(T){linebreak} $\alpha$: {alpha_T:.3f}{pm}{error:.3f}'.format(frac = r'$_{frac}$',alpha_T = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
plt.title(r'$p_{frac}$ vs T')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('T')
plt.ylabel(r'$p_{frac}$')
plt.tight_layout()
plt.legend()
plt.show()

P_n = p_array*((np.nanmedian(nh2_array)/nh2_array)**(NH2_param[1]))

plt.figure(figsize=(12,8))
plt.scatter(temp_array,P_n)
param,x,y = Functions.binning_equal_width(temp_array,P_n,np.ones_like(temp_array),ep_array,15)
valid = ~(np.isnan(temp_array)|np.isnan(P_n)|np.isnan(np.ones_like(temp_array))|np.isnan(ep_array))
data = (temp_array[valid],P_n[valid],np.ones_like(temp_array)[valid],ep_array[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'log(p{frac}) = C + $\alpha$log(T){linebreak} $\alpha$: {alpha_T:.3f}{pm}{error:.3f}'.format(frac = r'$_{frac}^{N_{dec}}$',alpha_T = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'tab:green',linewidth = 2.5,label = label_temp)
plt.title(r'$p_{frac}^{N_{dec}}$ vs T')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('T')
plt.ylabel(r'$p_{frac}^{N_{dec}}$')
plt.tight_layout()
plt.legend()
plt.show()
plt.show()

#%%
########## Double parameter fitting

# Functions.binning_equal_width_2D(s_array,nh2_array,p_array,es_array,np.ones_like(nh2_array),ep_array,50,50)

valid = ~(np.isnan(s_array)|np.isnan(I_array)|np.isnan(p_array)|np.isnan(es_array)|np.isnan(eI_array)|np.isnan(ep_array))
PSI_param, PSI_param_cov = Functions.binning_equal_width_2D(s_array[valid],I_array[valid],p_array[valid],es_array[valid],eI_array[valid],ep_array[valid],50,50)

print('Double parameter fitting')
print('p vs IS','\n','C :',PSI_param[0])
print('s index :',PSI_param[1])
print('I index :',PSI_param[2])
print('s index error :', np.sqrt(PSI_param_cov[1,1]))
print('I index error :', np.sqrt(PSI_param_cov[2,2]),'\n')

#%%
valid = ~(np.isnan(s_array)|np.isnan(nh2_array)|np.isnan(p_array)|np.isnan(es_array)|np.isnan(np.ones_like(nh2_array))|np.isnan(ep_array))
PSNH2_param, PSNH2_param_cov = Functions.binning_equal_width_2D(s_array[valid],nh2_array[valid],p_array[valid],es_array[valid],np.ones_like(nh2_array)[valid],ep_array[valid],50,50)
print('p vs NS','\n','C :',PSNH2_param[0])
print('s index :',PSNH2_param[1])
print('nh2 index :',PSNH2_param[2])
print('s index error :', np.sqrt(PSNH2_param_cov[1,1]))
print('nh2 index error :', np.sqrt(PSNH2_param_cov[2,2]))

#%%

################ Case 2

s_C2_R1 = S_map_deb.copy()
es_C2_R1 = sigma_S_map.copy()
p_C2_R1 = BlankedMapPol.copy()
ep_C2_R1 = BlankedMapDebPolError.copy()
I_C2_R1 = BlankedMapStokesI.copy()
eI_C2_R1 = BlankedMapStokesIError.copy()
nh2_C2_R1 = BlankedMapColumnDensity.copy()
micron8_C2_R1 = BlankedMap8Mircon.copy()
temp_C2_R1 = BlankedMapTemperature.copy()

s_C2_rest = S_map_deb.copy()
es_C2_rest = sigma_S_map.copy()
p_C2_rest = BlankedMapPol.copy()
ep_C2_rest = BlankedMapDebPolError.copy()
I_C2_rest = BlankedMapStokesI.copy()
eI_C2_rest = BlankedMapStokesIError.copy()
nh2_C2_rest = BlankedMapColumnDensity.copy()
micron8_C2_rest = BlankedMap8Mircon.copy()
temp_C2_rest = BlankedMapTemperature.copy()

Selector = (Mask.data != 1)

s_C2_R1.data[Selector] = np.nan
es_C2_R1.data[Selector] = np.nan
p_C2_R1.data[Selector] = np.nan
ep_C2_R1.data[Selector] = np.nan
I_C2_R1.data[Selector] = np.nan
eI_C2_R1.data[Selector] = np.nan
nh2_C2_R1.data[Selector] = np.nan
micron8_C2_R1.data[Selector] = np.nan
temp_C2_R1.data[Selector] = np.nan


Selector = ~Selector
s_C2_rest.data[Selector] = np.nan
es_C2_rest.data[Selector] = np.nan
p_C2_rest.data[Selector] = np.nan
ep_C2_rest.data[Selector] = np.nan
I_C2_rest.data[Selector] = np.nan
eI_C2_rest.data[Selector] = np.nan
nh2_C2_rest.data[Selector] = np.nan
micron8_C2_rest.data[Selector] = np.nan
temp_C2_rest.data[Selector] = np.nan

Selector = p_C2_R1.data < 0.5
s_C2_R1.data[Selector] = np.nan
es_C2_R1.data[Selector] = np.nan



Selector = p_C2_rest.data < 0.5
s_C2_rest.data[Selector] = np.nan
es_C2_rest.data[Selector] = np.nan



s_C2_R1_array = s_C2_R1.data.flatten()[::4]
es_C2_R1_array = es_C2_R1.data.flatten()[::4]
p_C2_R1_array = p_C2_R1.data.flatten()[::4]
ep_C2_R1_array = ep_C2_R1.data.flatten()[::4]
I_C2_R1_array = I_C2_R1.data.flatten()[::4]
eI_C2_R1_array = eI_C2_R1.data.flatten()[::4]
nh2_C2_R1_array = nh2_C2_R1.data.flatten()[::4]
micron8_C2_R1_array = micron8_C2_R1.data.flatten()[::4]
temp_C2_R1_array = temp_C2_R1.data.flatten()[::4]

s_C2_rest_array = s_C2_rest.data.flatten()[::4]
es_C2_rest_array = es_C2_rest.data.flatten()[::4]
p_C2_rest_array = p_C2_rest.data.flatten()[::4]
ep_C2_rest_array = ep_C2_rest.data.flatten()[::4]
I_C2_rest_array = I_C2_rest.data.flatten()[::4]
eI_C2_rest_array = eI_C2_rest.data.flatten()[::4]
nh2_C2_rest_array = nh2_C2_rest.data.flatten()[::4]
micron8_C2_rest_array = micron8_C2_rest.data.flatten()[::4]
temp_C2_rest_array = temp_C2_rest.data.flatten()[::4]

# %%
# array1 = I_C2_R1_array
# weights1 = eI_C2_R1_array

# array2 = p_C2_R1_array
# weights2 = ep_C2_R1_array

# array3 =I_C2_rest_array
# weights3 = eI_C2_rest_array

# array4 = p_C2_rest_array
# weights4 = ep_C2_rest_array

# plt.figure(figsize=(12,8))
# plt.scatter(array1,array2,c='#4575b4',label = 'R1')
# plt.scatter(array3,array4, c='#d73027',label = r'R2,R3,R4,R5,R6{linebreak}log(p{frac}) = C + $\alpha${I}log($I$)'.format(frac = r'$_{frac}$',linebreak='\n',I = r'$_{I}$'))

# param,x,y = Functions.binning_equal_width(array1,array2,weights1,weights2,10)
# valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(weights1)|np.isnan(weights2))
# data = (array1[valid],array2[valid],weights1[valid],weights2[valid])
# res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
# label_temp = r'$\alpha${R}: {alpha_I:.3f}{pm}{error:.3f}'.format(R = r'$_{R1}$',alpha_I = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

# param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,10)
# valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
# data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
# res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
# label_temp = r'$\alpha${R}: {alpha_I:.3f}{pm}{error:.3f}'.format(R = r'$_{rest}$',alpha_I = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
# plt.plot(x,y,'#fdae61',linewidth = 3.5,label = label_temp)
# plt.title(r'$p_{frac}$ vs I')
# plt.xscale('log')
# plt.yscale('log')
# plt.yticks([1e0,1e1])
# # plt.xticks(np.array([1e22,1e23,1e24]))
# plt.xlabel('I')
# plt.ylabel(r'$p_{frac}$')
# plt.legend()
# plt.tight_layout()
# plt.show()

#%%

array1 = nh2_C2_R1_array
weights1 = np.ones_like(nh2_C2_R1_array)

array2 = p_C2_R1_array
weights2 = ep_C2_R1_array

array3 =nh2_C2_rest_array
weights3 = np.ones_like(nh2_C2_rest_array)

array4 = p_C2_rest_array
weights4 = ep_C2_rest_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,c='#4575b4',label = 'R1')
plt.scatter(array3,array4, c='#BCBD22',label = r'R2,R3,R4,R5,R6{linebreak}log(p{frac}) = C + $\alpha${nh2}log($N$)'.format(frac = r'$_{frac}$',linebreak='\n',nh2 = r'$_{N}$'))

param,x,y = Functions.binning_equal_width(array1,array2,weights1,weights2,15)
valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(weights1)|np.isnan(weights2))
data = (array1[valid],array2[valid],weights1[valid],weights2[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R= r'$_{R1}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'tab:red',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,15,"#E377C2","#7F7F7F")
valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R= r'$_{rest}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#17BECF',linewidth = 3.5,label = label_temp)
plt.title(r'$p_{frac}$ vs N')
plt.xscale('log')
plt.yscale('log')
plt.yticks([1e0,1e1])
plt.xticks(np.array([1e22,1e23,1e24]))
plt.xlabel('N')
plt.ylabel(r'$p_{frac}$')
plt.legend()
plt.tight_layout()
plt.show()

#%%

array1 = nh2_C2_R1_array
weights1 = np.ones_like(nh2_C2_R1_array)

array2 = s_C2_R1_array
weights2 = es_C2_R1_array

array3 =nh2_C2_rest_array
weights3 = np.ones_like(nh2_C2_rest_array)

array4 = s_C2_rest_array
weights4 = es_C2_rest_array


plt.figure(figsize=(12,8))
plt.scatter(array1,array2,c='#4575b4',label = 'R1')
plt.scatter(array3,array4, c='#BCBD22',label = r'R2,R3,R4,R5,R6{linebreak}log(S) = C + $\alpha$log($N$)'.format(linebreak='\n'))


param,x,y = Functions.binning_equal_width(array1,array2,weights1,weights2,15)
valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(weights1)|np.isnan(weights2))
data = (array1[valid],array2[valid],weights1[valid],weights2[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R= r'$_{R1}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'tab:red',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,15,"#E377C2","#7F7F7F")
valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R= r'$_{rest}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#17BECF',linewidth = 3.5,label = label_temp)
plt.title(r'$S$ vs N')
plt.xscale('log')
plt.yscale('log')
plt.yticks([1e0,1e1])
plt.xticks(np.array([1e22,1e23,1e24]))
plt.xlabel('N')
plt.ylabel(r'$S$')
plt.legend()
plt.tight_layout()
plt.show()


#%%

array1 = s_C2_R1_array
weights1 = es_C2_R1_array

array2 = p_C2_R1_array
weights2 = ep_C2_R1_array

array3 =s_C2_rest_array
weights3 = es_C2_rest_array

array4 = p_C2_rest_array
weights4 = ep_C2_rest_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,c='#4575b4',label = 'R1')
plt.scatter(array3,array4, c='#d73027',label = r'R2,R3,R4,R5,R6{linebreak}log(p{frac}) = C + $\alpha_s$log(S)'.format(frac = r'$_{frac}$',linebreak='\n'))

param,x,y = Functions.binning_equal_width(array1,array2,weights1,weights2,10)
valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(weights1)|np.isnan(weights2))
data = (array1[valid],array2[valid],weights1[valid],weights2[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_s:.3f}{pm}{error:.3f}'.format(R= r'$_{R1}$',alpha_s = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,10)
valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_s:.3f}{pm}{error:.3f}'.format(R= r'$_{rest}$',alpha_s = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#fdae61',linewidth = 3.5,label = label_temp)

plt.title(r'$p_{frac}$ vs S')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('S')
plt.ylabel(r'$p_{frac}$')
plt.tight_layout()
plt.legend()
plt.show()

#%%
array2 = s_C2_R1_array.copy()
weights2 = es_C2_R1_array.copy()

array1 = p_C2_R1_array.copy()
weights1 = ep_C2_R1_array.copy()

array4 =s_C2_rest_array.copy()
weights4 = es_C2_rest_array.copy()

array3 = p_C2_rest_array.copy()
weights3 = ep_C2_rest_array.copy()

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,c='#4575b4',label = 'R1')
plt.scatter(array3,array4, c='#BCBD22',label = r'R2,R3,R4,R5,R6{linebreak}log(S) = C + $\alpha$log(p{frac})'.format(frac = r'$_{frac}$',linebreak='\n'))


param,x,y = Functions.binning_equal_width(array1,array2,weights1,weights2,15)
valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(weights1)|np.isnan(weights2))
data = (array1[valid],array2[valid],weights1[valid],weights2[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R= r'$_{R1}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'tab:red',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,15,"#E377C2","#7F7F7F")
valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R= r'$_{rest}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#17BECF',linewidth = 3.5,label = label_temp)
# plt.xlim(1.3,max(np.nanmax(p_C2_R1_array),np.nanmax(p_C2_rest_array)))
plt.xscale('log')
plt.yscale('log')
plt.ylabel('S')
plt.xlabel(r'$p_{frac}$')
plt.tight_layout()
plt.legend()
plt.show()

#%%

array1 = nh2_C2_R1_array
weights1 = np.ones_like(nh2_C2_R1_array)

array2 = p_C2_R1_array*s_C2_R1_array
weights2 = ep_C2_R1_array*s_C2_R1_array + es_C2_R1_array*p_C2_R1_array

array3 =nh2_C2_rest_array
weights3 = np.ones_like(nh2_C2_rest_array)

array4 = p_C2_rest_array*s_C2_rest_array
weights4 = ep_C2_rest_array*s_C2_rest_array + es_C2_rest_array*p_C2_rest_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,c='#4575b4',label = 'R1')
plt.scatter(array3,array4, c='#BCBD22',label = r'R2,R3,R4,R5,R6{linebreak}log(p{frac}*S) = C + $\alpha$log($N$)'.format(frac = r'$_{frac}$',linebreak='\n'))


param,x,y = Functions.binning_equal_width(array1,array2,weights1,weights2,15)
valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(weights1)|np.isnan(weights2))
data = (array1[valid],array2[valid],weights1[valid],weights2[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R= r'$_{R1}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'tab:red',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,15,"#E377C2","#7F7F7F")
valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R= r'$_{rest}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#17BECF',linewidth = 3.5,label = label_temp)

plt.title(r'$p_{frac}*S$ vs N')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel(r'$p_{frac}*S$')
plt.tight_layout()
plt.legend()
plt.show()

#%%

array1 = temp_C2_R1_array
weights1 = np.ones_like(temp_C2_R1_array)

array2 = p_C2_R1_array
weights2 = ep_C2_R1_array

array3 =temp_C2_rest_array
weights3 = np.ones_like(temp_C2_rest_array)

array4 = p_C2_rest_array
weights4 = ep_C2_rest_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,c='#4575b4',label = 'R1')
plt.scatter(array3,array4, c='#d73027',label = r'R2,R3,R4,R5,R6{linebreak}log(p{frac}) = C + $\alpha$log($T$)'.format(frac = r'$_{frac}$',linebreak='\n'))

param,x,y = Functions.binning_equal_width(array1,array2,weights1,weights2,10)
valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(weights1)|np.isnan(weights2))
data = (array1[valid],array2[valid],weights1[valid],weights2[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_T:.3f}{pm}{error:.3f}'.format(R= r'$_{R1}$',alpha_T = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,10)
valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_T:.3f}{pm}{error:.3f}'.format(R= r'$_{rest}$',alpha_T = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#fdae61',linewidth = 3.5,label = label_temp)

plt.title(r'$p_{frac}$ vs T')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('T')
plt.ylabel(r'$p_{frac}$')
plt.tight_layout()
plt.legend()
plt.show()

#%%

array2 = temp_C2_R1_array
weights2 = np.ones_like(temp_C2_R1_array)

array1 = nh2_C2_R1_array
weights1 = np.ones_like(nh2_C2_R1_array)

array4 =temp_C2_rest_array
weights4 = np.ones_like(temp_C2_rest_array)

array3 = nh2_C2_rest_array
weights3 = np.ones_like(nh2_C2_rest_array)

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,c='#4575b4',label = 'R1')
plt.scatter(array3,array4, c='#d73027',label = r'R2,R3,R4,R5,R6{linebreak}log(T) = C + $\alpha$log($N$)'.format(frac = r'$_{frac}$',linebreak='\n'))

param,x,y = Functions.binning_equal_width(array1,array2,weights1,weights2,10)
valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(weights1)|np.isnan(weights2))
data = (array1[valid],array2[valid],weights1[valid],weights2[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_T:.3f}{pm}{error:.3f}'.format(R= r'$_{R1}$',alpha_T = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,10)
valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_T:.3f}{pm}{error:.3f}'.format(R= r'$_{rest}$',alpha_T = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#fdae61',linewidth = 3.5,label = label_temp)

plt.title(r'$T$ vs N')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel(r'$T$')
plt.tight_layout()
plt.legend()
plt.show()


#%%

array1 = temp_C2_R1_array
weights1 = np.ones_like(temp_C2_R1_array)

array2 = p_C2_R1_array
weights2 = ep_C2_R1_array

array3 =temp_C2_rest_array
weights3 = np.ones_like(temp_C2_rest_array)

array4 = p_C2_rest_array
weights4 = ep_C2_rest_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,c='#4575b4',label = 'R1')
plt.scatter(array3,array4, c='#d73027',label = r'R2,R3,R4,R5,R6{linebreak}log(P{frac}) = C + $\alpha$log($T$)'.format(frac = r'$_{frac}$',linebreak='\n'))

param,x,y = Functions.binning_equal_width(array1,array2,weights1,weights2,10)
valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(weights1)|np.isnan(weights2))
data = (array1[valid],array2[valid],weights1[valid],weights2[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_T:.3f}{pm}{error:.3f}'.format(R= r'$_{R1}$',alpha_T = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,10)
valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_T:.3f}{pm}{error:.3f}'.format(R= r'$_{rest}$',alpha_T = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#fdae61',linewidth = 3.5,label = label_temp)

plt.title(r'$P_{frac}$ vs T')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('T')
plt.xlabel(r'$P_{frac}$')
plt.tight_layout()
plt.legend()
plt.show()

#%%

array1 = temp_C2_R1_array
weights1 = np.ones_like(temp_C2_R1_array)

NH2_param ,_,_ = Functions.binning_equal_width(nh2_C2_R1_array,p_C2_R1_array,np.ones_like(nh2_C2_R1_array),ep_C2_R1_array,10)
array2 = p_C2_R1_array*((np.nanmedian(nh2_C2_R1_array)/nh2_C2_R1_array)**(NH2_param[1]))
weights2 = ep_C2_R1_array

array3 =temp_C2_rest_array
weights3 = np.ones_like(temp_C2_rest_array)

NH2_param ,_,_ = Functions.binning_equal_width(nh2_C2_rest_array,p_C2_rest_array,np.ones_like(nh2_C2_rest_array),ep_C2_rest_array,10)
array4 = p_C2_rest_array*((np.nanmedian(nh2_C2_rest_array)/nh2_C2_rest_array)**(NH2_param[1]))
weights4 = ep_C2_rest_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,c='#4575b4',label = 'R1')
plt.scatter(array3,array4, c='#d73027',label = r'R2,R3,R4,R5,R6{linebreak}log(p{frac}) = C + $\alpha$log($T$)'.format(frac = r'$_{frac}^{N_{dec}}$',linebreak='\n'))

param,x,y = Functions.binning_equal_width(array1,array2,weights1,weights2,10)
valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(weights1)|np.isnan(weights2))
data = (array1[valid],array2[valid],weights1[valid],weights2[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_T:.3f}{pm}{error:.3f}'.format(R= r'$_{R1}$',alpha_T = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,10)
valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_T:.3f}{pm}{error:.3f}'.format(R= r'$_{rest}$',alpha_T = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#fdae61',linewidth = 3.5,label = label_temp)

plt.title(r'$p_{frac}^{N_{dec}}$ vs T')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('T')
plt.ylabel(r'$p_{frac}^{N_{dec}}$')
plt.tight_layout()
plt.legend()
plt.show()
plt.show()



#%%

########## double parameter fitting

array1 =s_C2_R1_array
weights1 = es_C2_R1_array

array2 = I_C2_R1_array
weights2 = eI_C2_R1_array

array3 = p_C2_R1_array
weights3 = ep_C2_R1_array


valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(array3)|np.isnan(weights1)|np.isnan(weights2)|np.isnan(weights3))
PSI_param, PSI_param_cov = Functions.binning_equal_width_2D(array1[valid],array2[valid],array3[valid],weights1[valid],weights2[valid],weights3[valid],50,50)
print('Ridge')
print('Double parameter fitting')
print('p vs IS','\n','C :',PSI_param[0])
print('s index :',PSI_param[1])
print('I index :',PSI_param[2])
print('s index error :', np.sqrt(PSI_param_cov[1,1]))
print('I index error :', np.sqrt(PSI_param_cov[2,2]),'\n')


array1 =s_C2_rest_array
weights1 = es_C2_rest_array

array2 = I_C2_rest_array
weights2 = eI_C2_rest_array

array3 = p_C2_rest_array
weights3 = ep_C2_rest_array


valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(array3)|np.isnan(weights1)|np.isnan(weights2)|np.isnan(weights3))
PSI_param, PSI_param_cov = Functions.binning_equal_width_2D(array1[valid],array2[valid],array3[valid],weights1[valid],weights2[valid],weights3[valid],50,50)
print("Rest")
print('Double parameter fitting')
print('p vs IS','\n','C :',PSI_param[0])
print('s index :',PSI_param[1])
print('I index :',PSI_param[2])
print('s index error :', np.sqrt(PSI_param_cov[1,1]))
print('I index error :', np.sqrt(PSI_param_cov[2,2]),'\n')


#%%
array1 =s_C2_R1_array
weights1 = es_C2_R1_array

array2 = nh2_C2_R1_array
weights2 = np.ones_like(nh2_C2_R1_array)

array3 = p_C2_R1_array
weights3 = ep_C2_R1_array



valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(array3)|np.isnan(weights1)|np.isnan(weights2)|np.isnan(weights3))
PSNH2_param, PSNH2_param_cov = Functions.binning_equal_width_2D(array1[valid],array2[valid],array3[valid],weights1[valid],weights2[valid],weights3[valid],50,50)
print('Ridge')
print('p vs NS','\n','C :',PSNH2_param[0])
print('s index :',PSNH2_param[1])
print('nh2 index :',PSNH2_param[2])
print('s index error :', np.sqrt(PSNH2_param_cov[1,1]))
print('nh2 index error :', np.sqrt(PSNH2_param_cov[2,2]))


array1 =s_C2_rest_array
weights1 = es_C2_rest_array

array2 = nh2_C2_rest_array
weights2 = np.ones_like(nh2_C2_rest_array)

array3 = p_C2_rest_array
weights3 = ep_C2_rest_array



valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(array3)|np.isnan(weights1)|np.isnan(weights2)|np.isnan(weights3))
PSNH2_param, PSNH2_param_cov = Functions.binning_equal_width_2D(array1[valid],array2[valid],array3[valid],weights1[valid],weights2[valid],weights3[valid],50,50)
print("Rest")
print('p vs NS','\n','C :',PSNH2_param[0])
print('s index :',PSNH2_param[1])
print('nh2 index :',PSNH2_param[2])
print('s index error :', np.sqrt(PSNH2_param_cov[1,1]))
print('nh2 index error :', np.sqrt(PSNH2_param_cov[2,2]))



# %%


########### case 3


s_C3_R1 = S_map_deb.copy()
es_C3_R1 = sigma_S_map.copy()
p_C3_R1 = BlankedMapPol.copy()
ep_C3_R1 = BlankedMapDebPolError.copy()
I_C3_R1 = BlankedMapStokesI.copy()
eI_C3_R1 = BlankedMapStokesIError.copy()
nh2_C3_R1 = BlankedMapColumnDensity.copy()
micron8_C3_R1 = BlankedMap8Mircon.copy()
temp_C3_R1 = BlankedMapTemperature.copy()

s_C3_rest = S_map_deb.copy()
es_C3_rest = sigma_S_map.copy()
p_C3_rest = BlankedMapPol.copy()
ep_C3_rest = BlankedMapDebPolError.copy()
I_C3_rest = BlankedMapStokesI.copy()
eI_C3_rest = BlankedMapStokesIError.copy()
nh2_C3_rest = BlankedMapColumnDensity.copy()
micron8_C3_rest = BlankedMap8Mircon.copy()
temp_C3_rest = BlankedMapTemperature.copy()

s_C3_R6 = S_map_deb.copy()
es_C3_R6 = sigma_S_map.copy()
p_C3_R6 = BlankedMapPol.copy()
ep_C3_R6 = BlankedMapDebPolError.copy()
I_C3_R6 = BlankedMapStokesI.copy()
eI_C3_R6 = BlankedMapStokesIError.copy()
nh2_C3_R6 = BlankedMapColumnDensity.copy()
micron8_C3_R6 = BlankedMap8Mircon.copy()
temp_C3_R6 = BlankedMapTemperature.copy()


Selector = (Mask.data != 1)
s_C3_R1.data[Selector] = np.nan
es_C3_R1.data[Selector] = np.nan
p_C3_R1.data[Selector] = np.nan
ep_C3_R1.data[Selector] = np.nan
I_C3_R1.data[Selector] = np.nan
eI_C3_R1.data[Selector] = np.nan
nh2_C3_R1.data[Selector] = np.nan
micron8_C3_R1.data[Selector] = np.nan
temp_C3_R1.data[Selector] = np.nan


Selector = (Mask.data == 1) + (Mask.data == 6)
s_C3_rest.data[Selector] = np.nan
es_C3_rest.data[Selector] = np.nan
p_C3_rest.data[Selector] = np.nan
ep_C3_rest.data[Selector] = np.nan
I_C3_rest.data[Selector] = np.nan
eI_C3_rest.data[Selector] = np.nan
nh2_C3_rest.data[Selector] = np.nan
micron8_C3_rest.data[Selector] = np.nan
temp_C3_rest.data[Selector] = np.nan

Selector = Mask.data != 6
s_C3_R6.data[Selector] = np.nan
es_C3_R6.data[Selector] = np.nan
p_C3_R6.data[Selector] = np.nan
ep_C3_R6.data[Selector] = np.nan
I_C3_R6.data[Selector] = np.nan
eI_C3_R6.data[Selector] = np.nan
nh2_C3_R6.data[Selector] = np.nan
micron8_C3_R6.data[Selector] = np.nan
temp_C3_R6.data[Selector] = np.nan
#%%

s_C3_R1_array = s_C3_R1.data.flatten()[::4]
es_C3_R1_array = es_C3_R1.data.flatten()[::4]
p_C3_R1_array = p_C3_R1.data.flatten()[::4]
ep_C3_R1_array = ep_C3_R1.data.flatten()[::4]
I_C3_R1_array = I_C3_R1.data.flatten()[::4]
eI_C3_R1_array = eI_C3_R1.data.flatten()[::4]
nh2_C3_R1_array = nh2_C3_R1.data.flatten()[::4]
micron8_C3_R1_array = micron8_C3_R1.data.flatten()[::4]
temp_C3_R1_array = temp_C3_R1.data.flatten()[::4]

s_C3_R6_array = s_C3_R6.data.flatten()[::4]
es_C3_R6_array = es_C3_R6.data.flatten()[::4]
p_C3_R6_array = p_C3_R6.data.flatten()[::4]
ep_C3_R6_array = ep_C3_R6.data.flatten()[::4]
I_C3_R6_array = I_C3_R6.data.flatten()[::4]
eI_C3_R6_array = eI_C3_R6.data.flatten()[::4]
nh2_C3_R6_array = nh2_C3_R6.data.flatten()[::4]
micron8_C3_R6_array = micron8_C3_R6.data.flatten()[::4]
temp_C3_R6_array = temp_C3_R6.data.flatten()[::4]

s_C3_rest_array = s_C3_rest.data.flatten()[::4]
es_C3_rest_array = es_C3_rest.data.flatten()[::4]
p_C3_rest_array = p_C3_rest.data.flatten()[::4]
ep_C3_rest_array = ep_C3_rest.data.flatten()[::4]
I_C3_rest_array = I_C3_rest.data.flatten()[::4]
eI_C3_rest_array = eI_C3_rest.data.flatten()[::4]
nh2_C3_rest_array = nh2_C3_rest.data.flatten()[::4]
micron8_C3_rest_array = micron8_C3_rest.data.flatten()[::4]
temp_C3_rest_array = temp_C3_rest.data.flatten()[::4]


#%%

array1 = I_C3_R1_array
weights1 = eI_C3_R1_array

array2 = p_C3_R1_array
weights2 = ep_C3_R1_array

array3 =I_C3_rest_array
weights3 = eI_C3_rest_array

array4 = p_C3_rest_array
weights4 = ep_C3_rest_array

array5 = I_C3_R6_array
weights5 = eI_C3_R6_array

array6 = p_C3_R6_array
weights6 = ep_C3_R6_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,c='#4575b4',label = 'R1')
plt.scatter(array3,array4,c='#d73027',label = 'R2,R3,R4,R5')
plt.scatter(array5,array6,c = '#7fbc41',label = r'R6{linebreak}log(p{frac}) = C + $\alpha$log(I)'.format(frac = r'$_{frac}$',linebreak='\n'))

param,x,y = Functions.binning_equal_width(array1,array2,weights1,weights2,10)
valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(weights1)|np.isnan(weights2))
data = (array1[valid],array2[valid],weights1[valid],weights2[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_I:.3f}{pm}{error:.3f}'.format(R = r'$_{R1}$',alpha_I = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#fdae61',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,10)
valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_I:.3f}{pm}{error:.3f}'.format(R = r'$_{rest}$',alpha_I = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array5,array6,weights5,weights6,10)
valid = ~(np.isnan(array5)|np.isnan(array6)|np.isnan(weights5)|np.isnan(weights6))
data = (array5[valid],array6[valid],weights5[valid],weights6[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_I:.3f}{pm}{error:.3f}'.format(R = r'$_{R6}$',alpha_I = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#ffff99',linewidth = 3.5,label = label_temp)

plt.title(r'$p_{frac}$ vs I')
plt.xscale('log')
plt.yscale('log')
plt.yticks([1e0,1e1])
plt.xlabel('I')
plt.ylabel(r'$p_{frac}$')
plt.legend()
plt.tight_layout()
plt.show()


# %%


array1 = nh2_C3_R1_array
weights1 = np.ones_like(nh2_C3_R1_array)

array2 = p_C3_R1_array
weights2 = ep_C3_R1_array

array3 =nh2_C3_rest_array
weights3 = np.ones_like(nh2_C3_rest_array)

array4 = p_C3_rest_array
weights4 = ep_C3_rest_array

array5 = nh2_C3_R6_array
weights5 = np.ones_like(nh2_C3_R6_array)

array6 = p_C3_R6_array
weights6 = ep_C3_R6_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,c='#4575b4',label = 'R1')
plt.scatter(array3,array4,c='#d73027',label = 'R2,R3,R4,R5')
plt.scatter(array5,array6,c = '#7fbc41',label = r'R6{linebreak}log(p{frac}) = C + $\alpha$log($N$)'.format(frac = r'$_{frac}$',linebreak='\n',nh2 = r'$_{N}$'))

param,x,y = Functions.binning_equal_width(array1,array2,weights1,weights2,10)
valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(weights1)|np.isnan(weights2))
data = (array1[valid],array2[valid],weights1[valid],weights2[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R = r'$_{R1}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#fdae61',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,10)
valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R = r'$_{rest}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array5,array6,weights5,weights6,10)
valid = ~(np.isnan(array5)|np.isnan(array6)|np.isnan(weights5)|np.isnan(weights6))
data = (array5[valid],array6[valid],weights5[valid],weights6[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R = r'$_{R6}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#ffff99',linewidth = 3.5,label = label_temp)

plt.title(r'$p_{frac}$ vs N')
plt.xscale('log')
plt.yscale('log')
plt.yticks([1e0,1e1])
plt.xticks(np.array([1e22,1e23,1e24]))
plt.xlabel('N')
plt.ylabel(r'$p_{frac}$')
plt.legend()
plt.tight_layout()
plt.show()


#%%

array1 = nh2_C3_R1_array
weights1 = np.ones_like(nh2_C3_R1_array)

array2 = s_C3_R1_array
weights2 = es_C3_R1_array

array3 =nh2_C3_rest_array
weights3 = np.ones_like(nh2_C3_rest_array)

array4 = s_C3_rest_array
weights4 = es_C3_rest_array

array5 = nh2_C3_R6_array
weights5 = np.ones_like(nh2_C3_R6_array)

array6 = s_C3_R6_array
weights6 = es_C3_R6_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,c='#4575b4',label = 'R1')
plt.scatter(array3,array4,c='#d73027',label = 'R2,R3,R4,R5')
plt.scatter(array5,array6,c = '#7fbc41',label = r'R6{linebreak}log(p{frac}) = C + $\alpha$log($N$)'.format(frac = r'$_{frac}$',linebreak='\n',nh2 = r'$_{N}$'))

param,x,y = Functions.binning_equal_width(array1,array2,weights1,weights2,10)
valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(weights1)|np.isnan(weights2))
data = (array1[valid],array2[valid],weights1[valid],weights2[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R = r'$_{R1}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#fdae61',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,10)
valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R = r'$_{rest}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array5,array6,weights5,weights6,10)
valid = ~(np.isnan(array5)|np.isnan(array6)|np.isnan(weights5)|np.isnan(weights6))
data = (array5[valid],array6[valid],weights5[valid],weights6[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_nh2:.3f}{pm}{error:.3f}'.format(R = r'$_{R6}$',alpha_nh2 = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#ffff99',linewidth = 3.5,label = label_temp)

plt.title(r'$p_{frac}$ vs N')
plt.xscale('log')
plt.yscale('log')
plt.yticks([1e0,1e1])
plt.xticks(np.array([1e22,1e23,1e24]))
plt.xlabel('N')
plt.ylabel(r'$p_{frac}$')
plt.legend()
plt.tight_layout()
plt.show()


# %%

array1 = s_C3_R1_array
weights1 = es_C3_R1_array

array2 = p_C3_R1_array
weights2 = ep_C3_R1_array

array3 =s_C3_rest_array
weights3 = es_C3_rest_array

array4 = p_C3_rest_array
weights4 = ep_C3_rest_array

array5 = s_C3_R6_array
weights5 = es_C3_R6_array

array6 = p_C3_R6_array
weights6 = ep_C3_R6_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,c='#4575b4',label = 'R1')
plt.scatter(array3,array4,c='#d73027',label = 'R2,R3,R4,R5')
plt.scatter(array5,array6,c = '#7fbc41',label = r'R6{linebreak}log(p{frac}) = C + $\alpha$log(S)'.format(frac = r'$_{frac}$',linebreak='\n'))

param,x,y = Functions.binning_equal_width(array1,array2,weights1,weights2,10)
valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(weights1)|np.isnan(weights2))
data = (array1[valid],array2[valid],weights1[valid],weights2[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_s:.3f}{pm}{error:.3f}'.format(R = r'$_{R1}$',alpha_s = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#fdae61',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,10)
valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_s:.3f}{pm}{error:.3f}'.format(R = r'$_{rest}$',alpha_s = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array5,array6,weights5,weights6,10)
valid = ~(np.isnan(array5)|np.isnan(array6)|np.isnan(weights5)|np.isnan(weights6))
data = (array5[valid],array6[valid],weights5[valid],weights6[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_s:.3f}{pm}{error:.3f}'.format(R = r'$_{R6}$',alpha_s = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#ffff99',linewidth = 3.5,label = label_temp)
plt.title(r'$p_{frac}$ vs S')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('S')
plt.ylabel(r'$p_{frac}$')
plt.tight_layout()
plt.legend()
plt.show()

# %%

array1 = nh2_C3_R1_array
weights1 = np.ones_like(nh2_C3_R1_array)

array2 = p_C3_R1_array*s_C3_R1_array
weights2 = ep_C3_R1_array*s_C3_R1_array + es_C2_R1_array*p_C3_R1_array

array3 =nh2_C3_rest_array
weights3 = np.ones_like(nh2_C3_rest_array)

array4 = p_C3_rest_array*s_C3_rest_array
weights4 = ep_C3_rest_array*s_C3_rest_array + es_C3_rest_array*p_C3_rest_array

array5 = nh2_C3_R6_array
weights5 = np.ones_like(nh2_C3_R6_array)

array6 = p_C3_R6_array*s_C3_R6_array
weights6 = ep_C3_R6_array*s_C3_R6_array + es_C3_R6_array*p_C3_R6_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,c='#4575b4',label = 'R1')
plt.scatter(array3,array4,c='#d73027',label = 'R2,R3,R4,R5')
plt.scatter(array5,array6,c = '#7fbc41',label = r'R6{linebreak}log(p{frac}*S) = C + $\alpha$log($N$)'.format(frac = r'$_{frac}$',linebreak='\n'))

param,x,y = Functions.binning_equal_width(array1,array2,weights1,weights2,10)
valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(weights1)|np.isnan(weights2))
data = (array1[valid],array2[valid],weights1[valid],weights2[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha:.3f}{pm}{error:.3f}'.format(R = r'$_{R1}$',alpha = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#fdae61',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,10)
valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha:.3f}{pm}{error:.3f}'.format(R = r'$_{rest}$',alpha = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array5,array6,weights5,weights6,10)
valid = ~(np.isnan(array5)|np.isnan(array6)|np.isnan(weights5)|np.isnan(weights6))
data = (array5[valid],array6[valid],weights5[valid],weights6[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha:.3f}{pm}{error:.3f}'.format(R = r'$_{R6}$',alpha = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#ffff99',linewidth = 3.5,label = label_temp)
plt.title(r'$p_{frac}*S$ vs $N$')
plt.xlabel('N')
plt.ylabel(r'$p_{frac}*S$')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.legend()
plt.show()


# %%


array1 = temp_C3_R1_array
weights1 = np.ones_like(temp_C3_R1_array)

array2 = p_C3_R1_array
weights2 = ep_C3_R1_array

array3 =temp_C3_rest_array
weights3 = np.ones_like(temp_C3_rest_array)

array4 = p_C3_rest_array
weights4 = ep_C3_rest_array

array5 = temp_C3_R6_array
weights5 = np.ones_like(temp_C3_R6_array)

array6 = p_C3_R6_array
weights6 = ep_C3_R6_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,c='#4575b4',label = 'R1')
plt.scatter(array3,array4,c='#d73027',label = 'R2,R3,R4,R5')
plt.scatter(array5,array6,c = '#7fbc41',label = r'R6{linebreak}log(p{frac}) = C + $\alpha$log(T)'.format(frac = r'$_{frac}$',linebreak='\n'))

param,x,y = Functions.binning_equal_width(array1,array2,weights1,weights2,10)
valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(weights1)|np.isnan(weights2))
data = (array1[valid],array2[valid],weights1[valid],weights2[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_T:.3f}{pm}{error:.3f}'.format(R = r'$_{R1}$',alpha_T = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#fdae61',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,10)
valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_T:.3f}{pm}{error:.3f}'.format(R = r'$_{rest}$',alpha_T = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array5,array6,weights5,weights6,10)
valid = ~(np.isnan(array5)|np.isnan(array6)|np.isnan(weights5)|np.isnan(weights6))
data = (array5[valid],array6[valid],weights5[valid],weights6[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_T:.3f}{pm}{error:.3f}'.format(R = r'$_{R6}$',alpha_T = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#ffff99',linewidth = 3.5,label = label_temp)
plt.title(r'$p_{frac}$ vs T')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('T')
plt.ylabel(r'$p_{frac}$')
plt.tight_layout()
plt.legend()
plt.show()

# %%


array1 = temp_C3_R1_array
weights1 = np.ones_like(temp_C3_R1_array)

NH2_param ,_,_ = Functions.binning_equal_width(nh2_C2_R1_array,p_C2_R1_array,np.ones_like(nh2_C2_R1_array),ep_C2_R1_array,10)
array2 = p_C3_R1_array*((np.nanmedian(nh2_C2_R1_array)/nh2_C2_R1_array)**(NH2_param[1]))
weights2 = ep_C3_R1_array

array3 =temp_C3_rest_array
weights3 = np.ones_like(temp_C3_rest_array)

NH2_param ,_,_ = Functions.binning_equal_width(nh2_C3_rest_array,p_C3_rest_array,np.ones_like(nh2_C3_rest_array),ep_C3_rest_array,10)
array4 = p_C3_rest_array*((np.nanmedian(nh2_C3_rest_array)/nh2_C3_rest_array)**(NH2_param[1]))
weights4 = ep_C3_rest_array

array5 = temp_C3_R6_array
weights5 = np.ones_like(temp_C3_R6_array)

NH2_param ,_,_ = Functions.binning_equal_width(nh2_C3_R6_array,p_C3_R6_array,np.ones_like(nh2_C3_R6_array),ep_C3_R6_array,10)
array6 = p_C3_R6_array*((np.nanmedian(nh2_C3_R6_array)/nh2_C3_R6_array)**(NH2_param[1]))
weights6 = ep_C3_R6_array

plt.figure(figsize=(12,8))
plt.scatter(array1,array2,c='#4575b4',label = 'R1')
plt.scatter(array3,array4,c='#d73027',label = 'R2,R3,R4,R5')
plt.scatter(array5,array6,c = '#7fbc41',label = r'R6{linebreak}log(p{frac}) = C + $\alpha$log(T)'.format(frac = r'$_{frac}^{N_{dec}}$',linebreak='\n'))

param,x,y = Functions.binning_equal_width(array1,array2,weights1,weights2,10)
valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(weights1)|np.isnan(weights2))
data = (array1[valid],array2[valid],weights1[valid],weights2[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_T:.3f}{pm}{error:.3f}'.format(R = r'$_{R1}$',alpha_T = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#fdae61',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array3,array4,weights3,weights4,10)
valid = ~(np.isnan(array3)|np.isnan(array4)|np.isnan(weights3)|np.isnan(weights4))
data = (array3[valid],array4[valid],weights3[valid],weights4[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_T:.3f}{pm}{error:.3f}'.format(R = r'$_{rest}$',alpha_T = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#abdda4',linewidth = 3.5,label = label_temp)

param,x,y = Functions.binning_equal_width(array5,array6,weights5,weights6,10)
valid = ~(np.isnan(array5)|np.isnan(array6)|np.isnan(weights5)|np.isnan(weights6))
data = (array5[valid],array6[valid],weights5[valid],weights6[valid])
res = bootstrapscipy(data, Functions.binning_equal_width_vbootstrap,n_resamples=10000,paired = True,vectorized=False, random_state=0)
label_temp = r'$\alpha${R}: {alpha_T:.3f}{pm}{error:.3f}'.format(R = r'$_{R6}$',alpha_T = param[1],pm = r'$\pm$',error=res.standard_error,linebreak='\n')
plt.plot(x,y,'#ffff99',linewidth = 3.5,label = label_temp)
plt.title(r'$p_{frac}^{N_{dec}}$ vs T')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('T')
plt.ylabel(r'$p_{frac}^{N_{dec}}$')
plt.tight_layout()
plt.legend()
plt.show()

# %%


array1 =s_C3_R1_array
weights1 = es_C3_R1_array

array2 = I_C3_R1_array
weights2 = eI_C3_R1_array

array3 = p_C3_R1_array
weights3 = ep_C3_R1_array


valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(array3)|np.isnan(weights1)|np.isnan(weights2)|np.isnan(weights3))
PSI_param, PSI_param_cov = Functions.binning_equal_width_2D(array1[valid],array2[valid],array3[valid],weights1[valid],weights2[valid],weights3[valid],50,50)
print('Ridge')
print('Double parameter fitting')
print('p vs IS','\n','C :',PSI_param[0])
print('s index :',PSI_param[1])
print('I index :',PSI_param[2])
print('s index error :', np.sqrt(PSI_param_cov[1,1]))
print('I index error :', np.sqrt(PSI_param_cov[2,2]),'\n')


array1 =s_C3_rest_array
weights1 = es_C3_rest_array

array2 = I_C3_rest_array
weights2 = eI_C3_rest_array

array3 = p_C3_rest_array
weights3 = ep_C3_rest_array


valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(array3)|np.isnan(weights1)|np.isnan(weights2)|np.isnan(weights3))
PSI_param, PSI_param_cov = Functions.binning_equal_width_2D(array1[valid],array2[valid],array3[valid],weights1[valid],weights2[valid],weights3[valid],50,50)
print("Rest")
print('Double parameter fitting')
print('p vs IS','\n','C :',PSI_param[0])
print('s index :',PSI_param[1])
print('I index :',PSI_param[2])
print('s index error :', np.sqrt(PSI_param_cov[1,1]))
print('I index error :', np.sqrt(PSI_param_cov[2,2]),'\n')


array1 =s_C3_R6_array
weights1 = es_C3_R6_array

array2 = I_C3_R6_array
weights2 = eI_C3_R6_array

array3 = p_C3_R6_array
weights3 = ep_C3_R6_array


valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(array3)|np.isnan(weights1)|np.isnan(weights2)|np.isnan(weights3))
PSI_param, PSI_param_cov = Functions.binning_equal_width_2D(array1[valid],array2[valid],array3[valid],weights1[valid],weights2[valid],weights3[valid],50,50)
print("R6")
print('Double parameter fitting')
print('p vs IS','\n','C :',PSI_param[0])
print('s index :',PSI_param[1])
print('I index :',PSI_param[2])
print('s index error :', np.sqrt(PSI_param_cov[1,1]))
print('I index error :', np.sqrt(PSI_param_cov[2,2]),'\n')

# %%


array1 =s_C3_R1_array
weights1 = es_C3_R1_array

array2 = nh2_C3_R1_array
weights2 = np.ones_like(nh2_C3_R1_array)

array3 = p_C3_R1_array
weights3 = ep_C3_R1_array

valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(array3)|np.isnan(weights1)|np.isnan(weights2)|np.isnan(weights3))
PSNH2_param, PSNH2_param_cov = Functions.binning_equal_width_2D(array1[valid],array2[valid],array3[valid],weights1[valid],weights2[valid],weights3[valid],50,50)
print('Ridge')
print('p vs NS','\n','C :',PSNH2_param[0])
print('s index :',PSNH2_param[1])
print('nh2 index :',PSNH2_param[2])
print('s index error :', np.sqrt(PSNH2_param_cov[1,1]))
print('nh2 index error :', np.sqrt(PSNH2_param_cov[2,2]))


array1 =s_C3_rest_array
weights1 = es_C3_rest_array

array2 = nh2_C3_rest_array
weights2 = np.ones_like(nh2_C3_rest_array)

array3 = p_C3_rest_array
weights3 = ep_C3_rest_array

valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(array3)|np.isnan(weights1)|np.isnan(weights2)|np.isnan(weights3))
PSNH2_param, PSNH2_param_cov = Functions.binning_equal_width_2D(array1[valid],array2[valid],array3[valid],weights1[valid],weights2[valid],weights3[valid],50,50)
print("Rest")
print('p vs NS','\n','C :',PSNH2_param[0])
print('s index :',PSNH2_param[1])
print('nh2 index :',PSNH2_param[2])
print('s index error :', np.sqrt(PSNH2_param_cov[1,1]))
print('nh2 index error :', np.sqrt(PSNH2_param_cov[2,2]))

array1 =s_C3_R6_array
weights1 = es_C3_R6_array

array2 = nh2_C3_R6_array
weights2 = np.ones_like(nh2_C3_R6_array)

array3 = p_C3_R6_array
weights3 = ep_C3_R6_array

valid = ~(np.isnan(array1)|np.isnan(array2)|np.isnan(array3)|np.isnan(weights1)|np.isnan(weights2)|np.isnan(weights3))
PSNH2_param, PSNH2_param_cov = Functions.binning_equal_width_2D(array1[valid],array2[valid],array3[valid],weights1[valid],weights2[valid],weights3[valid],50,50)
print("R6")
print('p vs NS','\n','C :',PSNH2_param[0])
print('s index :',PSNH2_param[1])
print('nh2 index :',PSNH2_param[2])
print('s index error :', np.sqrt(PSNH2_param_cov[1,1]))
print('nh2 index error :', np.sqrt(PSNH2_param_cov[2,2]))

# %%



# s_C3_R1 = S_map_deb.copy()
# es_C3_R1 = sigma_S_map.copy()
p_C3_R1 = BlankedMapPol.copy()
# ep_C3_R1 = BlankedMapDebPolError.copy()
I_C3_R1 = BlankedMapStokesI.copy()
# eI_C3_R1 = BlankedMapStokesIError.copy()
nh2_C3_R1 = BlankedMapColumnDensity.copy()
# micron8_C3_R1 = BlankedMap8Mircon.copy()
# temp_C3_R1 = BlankedMapTemperature.copy()

# s_C3_R2 = S_map_deb.copy()
# es_C3_R2 = sigma_S_map.copy()
p_C3_R2 = BlankedMapPol.copy()
# ep_C3_R2 = BlankedMapDebPolError.copy()
I_C3_R2 = BlankedMapStokesI.copy()
# eI_C3_R2 = BlankedMapStokesIError.copy()
nh2_C3_R2 = BlankedMapColumnDensity.copy()
# micron8_C3_R2 = BlankedMap8Mircon.copy()
# temp_C3_R2 = BlankedMapTemperature.copy()


# s_C3_R3 = S_map_deb.copy()
# es_C3_R3 = sigma_S_map.copy()
p_C3_R3 = BlankedMapPol.copy()
# ep_C3_R3 = BlankedMapDebPolError.copy()
I_C3_R3 = BlankedMapStokesI.copy()
# eI_C3_R3 = BlankedMapStokesIError.copy()
nh2_C3_R3 = BlankedMapColumnDensity.copy()
# micron8_C3_R3 = BlankedMap8Mircon.copy()
# temp_C3_R3 = BlankedMapTemperature.copy()


# s_C3_R4 = S_map_deb.copy()
# es_C3_R4 = sigma_S_map.copy()
p_C3_R4 = BlankedMapPol.copy()
# ep_C3_R4 = BlankedMapDebPolError.copy()
I_C3_R4 = BlankedMapStokesI.copy()
# eI_C3_R4 = BlankedMapStokesIError.copy()
nh2_C3_R4 = BlankedMapColumnDensity.copy()
# micron8_C3_R4 = BlankedMap8Mircon.copy()
# temp_C3_R4 = BlankedMapTemperature.copy()

# s_C3_R5 = S_map_deb.copy()
# es_C3_R5 = sigma_S_map.copy()
p_C3_R5 = BlankedMapPol.copy()
# ep_C3_R5 = BlankedMapDebPolError.copy()
I_C3_R5 = BlankedMapStokesI.copy()
# eI_C3_R5 = BlankedMapStokesIError.copy()
nh2_C3_R5 = BlankedMapColumnDensity.copy()
# micron8_C3_R5 = BlankedMap8Mircon.copy()
# temp_C3_R5 = BlankedMapTemperature.copy()

# s_C3_R6 = S_map_deb.copy()
# es_C3_R6 = sigma_S_map.copy()
p_C3_R6 = BlankedMapPol.copy()
# ep_C3_R6 = BlankedMapDebPolError.copy()
I_C3_R6 = BlankedMapStokesI.copy()
# eI_C3_R6 = BlankedMapStokesIError.copy()
nh2_C3_R6 = BlankedMapColumnDensity.copy()
# micron8_C3_R6 = BlankedMap8Mircon.copy()
# temp_C3_R6 = BlankedMapTemperature.copy()


Selector = (Mask.data != 1)
# s_C3_R1.data[Selector] = np.nan
# es_C3_R1.data[Selector] = np.nan
p_C3_R1.data[Selector] = np.nan
# ep_C3_R1.data[Selector] = np.nan
I_C3_R1.data[Selector] = np.nan
# eI_C3_R1.data[Selector] = np.nan
nh2_C3_R1.data[Selector] = np.nan
# micron8_C3_R1.data[Selector] = np.nan
# temp_C3_R1.data[Selector] = np.nan


# s_C3_R1_array = s_C3_R1.data.flatten()[::4]
# es_C3_R1_array = es_C3_R1.data.flatten()[::4]
p_C3_R1_array = p_C3_R1.data.flatten()[::4]
# ep_C3_R1.data.flatten()[::4]
I_C3_R1_array = I_C3_R1.data.flatten()[::4]
# eI_C3_R1.data.flatten()[::4]
nh2_C3_R1_array = nh2_C3_R1.data.flatten()[::4]
# micron8_C3_R1.data.flatten()[::4]
# temp_C3_R1.data.flatten()[::4]



Selector = (Mask.data != 2)
# s_C3_R2.data[Selector] = np.nan
# es_C3_R2.data[Selector] = np.nan
p_C3_R2.data[Selector] = np.nan
# ep_C3_R2.data[Selector] = np.nan
I_C3_R2.data[Selector] = np.nan
# eI_C3_R2.data[Selector] = np.nan
nh2_C3_R2.data[Selector] = np.nan
# micron8_C3_R2.data[Selector] = np.nan
# temp_C3_R2.data[Selector] = np.nan

# s_C3_R2.data.flatten()[::4]
# es_C3_R2.data.flatten()[::4]
p_C3_R2_array = p_C3_R2.data.flatten()[::4]
# ep_C3_R2.data.flatten()[::4]
I_C3_R2_array = I_C3_R2.data.flatten()[::4]
# eI_C3_R2.data.flatten()[::4]
nh2_C3_R2_array = nh2_C3_R2.data.flatten()[::4]
# micron8_C3_R2.data.flatten()[::4]
# temp_C3_R2.data.flatten()[::4]


Selector = (Mask.data != 3)
# s_C3_R3.data[Selector] = np.nan
# es_C3_R3.data[Selector] = np.nan
p_C3_R3.data[Selector] = np.nan
# ep_C3_R3.data[Selector] = np.nan
I_C3_R3.data[Selector] = np.nan
# eI_C3_R3.data[Selector] = np.nan
nh2_C3_R3.data[Selector] = np.nan
# micron8_C3_R3.data[Selector] = np.nan
# temp_C3_R3.data[Selector] = np.nan


# s_C3_R3.data.flatten()[::4]
# es_C3_R3.data.flatten()[::4]
p_C3_R3_array = p_C3_R3.data.flatten()[::4]
# ep_C3_R3.data.flatten()[::4]
I_C3_R3_array = I_C3_R3.data.flatten()[::4]
# eI_C3_R3.data.flatten()[::4]
nh2_C3_R3_array = nh2_C3_R3.data.flatten()[::4]
# micron8_C3_R3.data.flatten()[::4]
# temp_C3_R3.data.flatten()[::4]


Selector = (Mask.data != 4)
# s_C3_R4.data[Selector] = np.nan
# es_C3_R4.data[Selector] = np.nan
p_C3_R4.data[Selector] = np.nan
# ep_C3_R4.data[Selector] = np.nan
I_C3_R4.data[Selector] = np.nan
# eI_C3_R4.data[Selector] = np.nan
nh2_C3_R4.data[Selector] = np.nan
# micron8_C3_R4.data[Selector] = np.nan
# temp_C3_R4.data[Selector] = np.nan


# s_C3_R4.data.flatten()[::4]
# es_C3_R4.data.flatten()[::4]
p_C3_R4_array = p_C3_R4.data.flatten()[::4]
# ep_C3_R4.data.flatten()[::4]
I_C3_R4_array = I_C3_R4.data.flatten()[::4]
# eI_C3_R4.data.flatten()[::4]
nh2_C3_R4_array = nh2_C3_R4.data.flatten()[::4]
# micron8_C3_R4.data.flatten()[::4]
# temp_C3_R4.data.flatten()[::4]


Selector = (Mask.data != 5)
# s_C3_R5.data[Selector] = np.nan
# es_C3_R5.data[Selector] = np.nan
p_C3_R5.data[Selector] = np.nan
# ep_C3_R5.data[Selector] = np.nan
I_C3_R5.data[Selector] = np.nan
# eI_C3_R5.data[Selector] = np.nan
nh2_C3_R5.data[Selector] = np.nan
# micron8_C3_R5.data[Selector] = np.nan
# temp_C3_R5.data[Selector] = np.nan


# s_C3_R5.data.flatten()[::4]
# es_C3_R5.data.flatten()[::4]
p_C3_R5_array = p_C3_R5.data.flatten()[::4]
# ep_C3_R5.data.flatten()[::4]
I_C3_R5_array = I_C3_R5.data.flatten()[::4]
# eI_C3_R5.data.flatten()[::4]
nh2_C3_R5_array = nh2_C3_R5.data.flatten()[::4]
# micron8_C3_R5.data.flatten()[::4]
# temp_C3_R5.data.flatten()[::4]

Selector = (Mask.data != 6)
# s_C3_R6.data[Selector] = np.nan
# es_C3_R6.data[Selector] = np.nan
p_C3_R6.data[Selector] = np.nan
# ep_C3_R6.data[Selector] = np.nan
I_C3_R6.data[Selector] = np.nan
# eI_C3_R6.data[Selector] = np.nan
nh2_C3_R6.data[Selector] = np.nan
# micron8_C3_R6.data[Selector] = np.nan
# temp_C3_R6.data[Selector] = np.nan


# s_C3_R6.data.flatten()[::4]
# es_C3_R6.data.flatten()[::4]
p_C3_R6_array = p_C3_R6.data.flatten()[::4]
# ep_C3_R6.data.flatten()[::4]
I_C3_R6_array = I_C3_R6.data.flatten()[::4]
# eI_C3_R6.data.flatten()[::4]
nh2_C3_R6_array = nh2_C3_R6.data.flatten()[::4]
# micron8_C3_R6.data.flatten()[::4]
# temp_C3_R6.data.flatten()[::4]

plt.figure(figsize=(12,8))
plt.scatter(nh2_C3_R1_array,p_C3_R1_array,c='#4575b4',label = 'R1',s = 50)
plt.scatter(nh2_C3_R2_array,p_C3_R2_array,c='#d73027',label = 'R2',s = 45)
plt.scatter(nh2_C3_R3_array,p_C3_R3_array,c = 'y',label = 'R3',s = 40)
plt.scatter(nh2_C3_R4_array,p_C3_R4_array,c = 'orange',label = 'R4',s = 35)
plt.scatter(nh2_C3_R5_array,p_C3_R5_array,label = 'R5',s = 35)
plt.scatter(nh2_C3_R6_array,p_C3_R6_array,c = '#7fbc41',label = 'R6',s = 35)
plt.title(r'$p_{frac}$ vs N')
plt.xscale('log')
plt.yscale('log')
# plt.yticks([1e0,1e1])
plt.xlabel('N')
plt.ylabel(r'$p_{frac}$')
plt.legend()
plt.tight_layout()
plt.show()

# %%
