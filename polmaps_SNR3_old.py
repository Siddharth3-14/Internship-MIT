from aplpy import FITSFigure  
import matplotlib.pyplot as plt 
import astropy.io.fits as fits
import astropy.units as u 
import numpy as np 
import os as os
import six

def get_data(infile):
	#get data
	hawc = fits.open(infile)
	p    = hawc['DEBIASED PERCENT POL']
	perr = hawc['ERROR PERCENT POL']
	#pa   = hawc['ROTATED POL ANGLE']
	pa   = hawc['POL ANGLE']
	stkI = hawc['STOKES I']
	stkIerr = hawc['ERROR I']
	pi   = hawc['DEBIASED POL FLUX']
	pierr   = hawc['ERROR POL FLUX']

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



filename1='../FITS_file/CygX_E_OTFMAP.fits'

hawcp = fits.open(filename1)

x = hawcp['STOKES I'].data
y = hawcp['DEBIASED PERCENT POL'].data
y1=hawcp['ROTATED POL ANGLE'].data
y2 = hawcp['DEBIASED POL FLUX'].data

#plt.scatter(x,y2,color='k')


title = 'SIMPLIFI'


SNRp_cut = 3.0
p_cut = 50

# figure
width  = 50
height = 50
#vmax = 1.5
#vmin = 0.3
cmap = 'plasma'

title_size = 16
tick_labels = 15
label_plot = 15
label_colorbar = 15
tick_colorbar = 15
label_fontsize = 20

SNRi_cut = 100
scalevec = 0.5 #1px = scalevec * 1% pol 
vec_legend = 5.0



#### SCRIPT
fig = plt.figure(figsize=(13,10))

p,perr,pa,stkI,stkIerr,pi,pierr,pxscale = get_data(filename1)

#central coordinates
RA = (stkI.header['OBSRA']*u.hourangle).to(u.deg)
DEC = stkI.header['OBSDEC']*u.deg

#print(RA,DEC)

###Figure
gc = FITSFigure(stkI,figure=fig)

##STOKES I
#colorscale
#gc.show_colorscale(cmap=cmap,vmin=vmin,vmax=vmax,smooth=1,kernel='gauss')
gc.show_colorscale(smooth=1,kernel='gauss')

#colorbae
gc.add_colorbar()
gc.colorbar.set_axis_label_text('I (Jy/sqarcsec)')
gc.colorbar.set_axis_label_font(size=label_colorbar)
gc.colorbar.set_font(size=tick_colorbar)
#recenter

# plt.show()

# #contours
# sigmaI = np.nanstd(stkIerr.data)
# levelsI = sigmaI * 2**(np.arange(2,12,1))
# #levelsI = np.arange(3,155,3)*sigmaI
# gc.show_contour(colors='white',levels=levelsI,linewidth=0.1,\
# 				filled=False,smooth=1,kernel='box',alpha=0.4)


# #gc.show_contour(levels=levelsI,\
# #				filled=True,smooth=1,kernel='box',cmap='viridis')
# #beam
gc.add_beam(color='white',edgecolor='black',linestyle='solid', linewidth=3)
# #title
# gc.set_title(title,fontsize=title_size)



##Pol map
p = quality_cuts(stkI,stkIerr,p,perr,SNRp_cut,p_cut,SNRi_cut)
# p.data = p.data




#quaterbeam
gc.show_vectors(p,pa,scale=scalevec,step=2,color='black',linewidth=4.0)
gc.show_vectors(p,pa,scale=scalevec,step=2,color='white',linewidth=2.0)


# print(p,pa)


# #show label
# #gc.add_label(0.2,0.95,'Total Flux',fontsize=label_fontsize,\
# #	color='white',weight='bold',relative=True)

# #legend vector
vecscale = scalevec * pxscale/3600
gc.add_scalebar(vec_legend*vecscale,r'$p_{frac}$ ='+np.str(vec_legend),corner='bottom right',frame=True,color='black',facecolor='blue')

# #Figure parameters
gc.tick_labels.set_font(size=tick_labels)
gc.axis_labels.set_font(size=label_plot)



# #fig.savefig(figname,dpi=300)
# #os.system('open '+figname)


plt.show()

