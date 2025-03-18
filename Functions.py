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
from scipy.stats import bootstrap as bootstrapscipy
import math
plt.rcParams.update({'font.size': 18})
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16


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
    
def lin_fit(x, a, b):
    return a + b*x

# def single_curve_fitting(x):
#     param, param_cov = curve_fit(lin_fit,(x[:,2],x[:,1]),x[:,0])
#     return param

def DoubleParamFunc(X, a, b, c):
    x,y = X
    return a + b*x + c*y

def curve_fitting(x):
    param, param_cov = curve_fit(DoubleParamFunc,(x[:,2],x[:,1]),x[:,0])
    return param


def lin_fit(x, a, b):
    return a + b*x

def remove_nan(array1,array2):
    selector = ~np.isnan(array1)
    array1_fil = array1[selector]
    array2_fil = array2[selector]

    selector = ~np.isnan(array2_fil)
    array1_fil = array1_fil[selector]
    array2_fil = array2_fil[selector]
    return array1_fil,array2_fil




##########
def binning_equal_width(array1,array2,error_array1,error_array2,Nbins,color_bins='#06402B',color_error='#8C564B'):
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

        Selector = log_filtered1 < bins[i]
        temp_array1[Selector] = np.nan
        temp_array2[Selector] = np.nan
        temp_error1[Selector] = np.nan
        temp_error2[Selector] =  np.nan

        Selector =  log_filtered1 > bins[i+1]
        temp_array1[Selector] = np.nan
        temp_array2[Selector] = np.nan
        temp_error1[Selector] = np.nan
        temp_error2[Selector] =  np.nan

        ma_array = np.ma.MaskedArray(temp_array1, mask=np.isnan(temp_array1))
        ma_errors = np.ma.MaskedArray(temp_error1, mask=np.isnan(temp_error1))
        average_temp = np.ma.average(ma_array, weights=ma_errors)
        bin_centres.append(average_temp)

        ma_array = np.ma.MaskedArray(temp_array2, mask=np.isnan(temp_array2))
        ma_errors = np.ma.MaskedArray(temp_error2, mask=np.isnan(temp_error2))
        average_temp = np.ma.average(ma_array, weights=ma_errors)
        binned_data.append(average_temp)

        error_bar.append(np.nanstd(temp_array2))
        print(error_bar)
    bin_centres = np.array(bin_centres)
    binned_data = np.array(binned_data)
    error_bar = np.array(error_bar)
    
    plt.errorbar(bin_centres,binned_data,yerr=error_bar,c=color_error)
    plt.scatter(bin_centres,binned_data,s = 75,c=color_bins)
    plt.plot(bin_centres,binned_data,c=color_bins)

    bin_centres = np.log10(np.array(bin_centres))
    binned_data = np.log10(np.array(binned_data))
    valid = ~(np.isnan(bin_centres)|np.isnan(binned_data))

    # bin_centres,binned_data = remove_nan(bin_centres,binned_data)
    level_bins = np.linspace(np.amin(bin_centres),np.amax(bin_centres),11)
    param, PS_param_cov = curve_fit(lin_fit, bin_centres[valid], binned_data[valid])
    PS_FitFunc = lin_fit(level_bins,param[0],param[1])
    return param,10**level_bins,10**PS_FitFunc

def binning_equal_width_vbootstrap(array1,array2,weights_array1,weights_array2,Nbins=15):
    log_filtered1 = np.log10(array1)
    log_filtered2 = np.log10(array2)

    bins = np.linspace(np.nanmin(log_filtered1),np.nanmax(log_filtered1),Nbins)


    bin_centres = []
    binned_data = []
    error_bar = []

    for i in range(0,(bins.shape[0]-1)):
        temp_array1 = array1.copy()
        temp_array2 = array2.copy()
        temp_weights1 = weights_array1.copy()
        temp_weights2 = weights_array2.copy()

        Selector = log_filtered1 < bins[i]
        temp_array1[Selector] = np.nan
        temp_array2[Selector] = np.nan
        temp_weights1[Selector] = np.nan
        temp_weights2[Selector] =  np.nan

        Selector =  log_filtered1 > bins[i+1]
        temp_array1[Selector] = np.nan
        temp_array2[Selector] = np.nan
        temp_weights1[Selector] = np.nan
        temp_weights2[Selector] =  np.nan

        ma_array = np.ma.MaskedArray(temp_array1, mask=np.isnan(temp_array1))
        ma_weights = np.ma.MaskedArray(temp_weights1, mask=np.isnan(temp_weights1))
        average_temp = np.ma.average(ma_array, weights=ma_weights)
        bin_centres.append(average_temp)

        ma_array = np.ma.MaskedArray(temp_array2, mask=np.isnan(temp_array2))
        ma_weights = np.ma.MaskedArray(temp_weights2, mask=np.isnan(temp_weights2))
        average_temp = np.ma.average(ma_array, weights=ma_weights)
        # average_temp = np.ma.average(ma_array, weights=ma_weights)

        binned_data.append(average_temp)

        error_bar.append(np.nanstd(temp_array2))

    bin_centres = np.array(bin_centres)
    binned_data = np.array(binned_data)
    error_bar = np.array(error_bar)

    bin_centres = np.log10(np.array(bin_centres))
    binned_data = np.log10(np.array(binned_data))
    valid = ~(np.isnan(bin_centres)|np.isnan(binned_data))

    # bin_centres,binned_data = remove_nan(bin_centres,binned_data)
    param, PS_param_cov = curve_fit(lin_fit, bin_centres[valid], binned_data[valid])
    return param[1]





###########


def weightedmean(vals, valerrs, method = 'best', Niter = 200):
    '''returns weighted mean and standard error on weighted mean
        the method used to calculate the standard error on the mean differs:
        ML: Not recommended! maximum likelihood method from http://ned.ipac.caltech.edu/level5/Leo/Stats4_5.html
            this can severely underestimate the standard error on the mean
        best: standard error by Cochran 1977 as used here: https://www.cs.tufts.edu/~nr/cs257/archive/donald-gatz/weighted-standard-error.pdf
              fast, as accurate as bootstrap for most cases
        bootstrap: standard deviation of distribution of weighted means produced by bootstrap resampling, 
                    needs input Number of iterations to use Niter (prefer N > 200)
                    slowest but most accurate
    '''
    #vals, valerrs = np.array(vals), np.array(valerrs)
    weights = 1./valerrs**2
    wmerr2 = 1./np.nansum(weights)
    # weighted mean
    wm = wmerr2*np.nansum(vals*weights)
    
    if method == 'ML':
        wmerr = np.sqrt(wmerr2)
        
    elif method == 'best':
        n = len(vals)
        meanweight = np.nanmean(weights)
        termA = np.sum( np.power(( weights*vals -  meanweight* wm ),2) )
        termB = - 2*wm* np.nansum( (weights - meanweight) * (weights*vals - meanweight*wm))
        termC = wm**2 * np.nansum( np.power((weights - meanweight),2) )
        wmerr = np.sqrt( n/(n-1) * wmerr2**2 * (termA + termB + termC) )
        
    elif method == 'bootstrap':
        xmeans = np.zeros(Niter)
        sxm = np.zeros(Niter)
        for ii in range(Niter):
            # resample the measurements 
            resample_inds = bootstrap_resample(vals)

            # calculate weighted mean
            a, b = vals[resample_inds],valerrs[resample_inds]
            xmeans[ii],sxm[ii] = weightedmean(a,b)
            
        #standard deviation of xmeans distr (best estimate of the error on the mean by bootstrapping)
        wmerr = np.std(xmeans)
        
    return wm, wmerr



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



# def binning_equal_width(array1,array2,error_array1,error_array2,Nbins,min_n = 1,color_bins='#06402B',color_error='#8C564B'):
#     log_filtered1 = np.log10(array1)
#     log_filtered2 = np.log10(array2)

#     bins = np.linspace(np.nanmin(log_filtered1),np.nanmax(log_filtered1),Nbins)


#     bin_centres = []
#     binned_data = []
#     error_bar = []

#     for i in range(0,(bins.shape[0]-1)):
#         temp_array1 = array1.copy()
#         temp_array2 = array2.copy()
#         temp_error1 = error_array1.copy()
#         temp_error2 = error_array2.copy()
        
#         Selector = (log_filtered1 < bins[i])|(log_filtered1 > bins[i+1])
#         temp_array1_fil = temp_array1[Selector]
#         temp_array2_fil = temp_array2[Selector] 
#         temp_error1_fil = temp_error1[Selector]
#         temp_error2_fil = temp_error2[Selector] 

#         valid = ~(np.isnan(temp_array1_fil)|np.isnan(temp_array2_fil)|np.isnan(temp_error1_fil)|np.isnan(temp_error2_fil))
#         print(np.sum(valid))
#         temp_array1_fil = temp_array1_fil[valid]
#         temp_array2_fil = temp_array2_fil[valid]
#         temp_error1_fil = temp_error1_fil[valid]
#         temp_error2_fil = temp_error2_fil[valid]

#         if np.sum(valid) == 0:
#             bin_centres.append(np.nan)
#             binned_data.append(np.nan)
#             error_bar.append(np.nan)

#         elif (np.sum(valid) < min_n):
#             bin_centres.append(np.nan)
#             binned_data.append(np.nan)
#             error_bar.append(np.nan)

#         elif (np.sum(valid)  ==  1):

#             bin_centres.append(temp_array1_fil[0])
#             binned_data.append(temp_array2_fil[0])
#             error_bar.append(temp_error2_fil[0])

#             # bin_centres.append(temp_array1_fil[0])
#             # binned_data.append(temp_array2_fil[0])
#             # error_bar.append(temp_weights2_fil[0])
#         else:
#             bin_centre_temp,_ = np.average(temp_array1_fil,1.0/temp_error1_fil**2)
#             # bin_data_temp,_ = weightedmean(temp_array2_fil,temp_error2_fil)
#             bin_data_temp,_ = np.average(temp_array2_fil,1.0/temp_error2_fil**2)
#             # bin_data_temp,_ = weightedmean(temp_array2_fil,temp_weights2_fil)
#             binned_error = get_wstd(temp_array2_fil,1.0/temp_error2_fil**2)
            
#             bin_centres.append(bin_centre_temp)
#             binned_data.append(bin_data_temp)
#             error_bar.append(binned_error)

#     bin_centres = np.array(bin_centres)
#     binned_data = np.array(binned_data)
#     error_bar = np.array(error_bar)
    
#     plt.errorbar(bin_centres,binned_data,yerr=error_bar,c=color_error)
#     plt.scatter(bin_centres,binned_data,s = 75,c=color_bins)
#     plt.plot(bin_centres,binned_data,c=color_bins)

#     bin_centres = np.log10(np.array(bin_centres))
#     binned_data = np.log10(np.array(binned_data))
#     valid = ~(np.isnan(bin_centres)|np.isnan(binned_data))

#     # bin_centres,binned_data = remove_nan(bin_centres,binned_data)
#     level_bins = np.linspace(np.nanmin(bin_centres),np.nanmax(bin_centres),11)
#     param, PS_param_cov = curve_fit(lin_fit, bin_centres[valid], binned_data[valid])
#     PS_FitFunc = lin_fit(level_bins,param[0],param[1])
#     return param,10**level_bins,10**PS_FitFunc




# def binning_equal_width_vbootstrap(array1,array2,weights_array1,weights_array2,Nbins=15):
#     log_filtered1 = np.log10(array1)
#     log_filtered2 = np.log10(array2)

#     bins = np.linspace(np.nanmin(log_filtered1),np.nanmax(log_filtered1),Nbins)


#     bin_centres = []
#     binned_data = []
#     error_bar = []

#     for i in range(0,(bins.shape[0]-1)):
#         temp_array1 = array1.copy()
#         temp_array2 = array2.copy()
#         temp_error1 = error_array1.copy()
#         temp_error2 = error_array2.copy()
        
#         Selector = (log_filtered1 < bins[i])|(log_filtered1 > bins[i+1])
#         temp_array1_fil = temp_array1[Selector]
#         temp_array2_fil = temp_array2[Selector] 
#         temp_error1_fil = temp_error1[Selector]
#         temp_error2_fil = temp_error2[Selector] 

#         valid = ~(np.isnan(temp_array1_fil)|np.isnan(temp_array2_fil)|np.isnan(temp_error1_fil)|np.isnan(temp_error2_fil))
#         print(np.sum(valid))
#         temp_array1_fil = temp_array1_fil[valid]
#         temp_array2_fil = temp_array2_fil[valid]
#         temp_error1_fil = temp_error1_fil[valid]
#         temp_error2_fil = temp_error2_fil[valid]

#         if np.sum(valid) == 0:
#             bin_centres.append(np.nan)
#             binned_data.append(np.nan)
#             error_bar.append(np.nan)

#         elif (np.sum(valid) < min_n):
#             bin_centres.append(np.nan)
#             binned_data.append(np.nan)
#             error_bar.append(np.nan)

#         elif (np.sum(valid)  ==  1):

#             bin_centres.append(temp_array1_fil[0])
#             binned_data.append(temp_array2_fil[0])
#             error_bar.append(temp_error2_fil[0])

#         else:
#             bin_centre_temp,_ = np.average(temp_array1_fil,1.0/temp_error1_fil**2)
#             bin_data_temp,_ = np.average(temp_array2_fil,1.0/temp_error2_fil**2)
#             binned_error = get_wstd(temp_array2_fil,1.0/temp_error2_fil**2)
            
#             bin_centres.append(bin_centre_temp)
#             binned_data.append(bin_data_temp)
#             error_bar.append(binned_error)

#     bin_centres = np.array(bin_centres)
#     binned_data = np.array(binned_data)
#     error_bar = np.array(error_bar)

#     bin_centres = np.log10(np.array(bin_centres))
#     binned_data = np.log10(np.array(binned_data))
#     valid = ~(np.isnan(bin_centres)|np.isnan(binned_data))

#     # bin_centres,binned_data = remove_nan(bin_centres,binned_data)
#     param, PS_param_cov = curve_fit(lin_fit, bin_centres[valid], binned_data[valid])
#     return param[1]



############# 2D binning
def binning_equal_width_2D(array1,array2,array3,weights_array1,weights_array2,weights_array3,Nbins1,Nbins2):
    log_filtered1 = np.log10(array1)
    log_filtered2 = np.log10(array2)
    log_filtered3 = np.log10(array3)


    bins1 = np.linspace(np.nanmin(log_filtered1),np.nanmax(log_filtered1),Nbins1)
    bins2 = np.linspace(np.nanmin(log_filtered2),np.nanmax(log_filtered2),Nbins2)


    bin_centres1 = []
    bin_centres2 = []

    binned_data = []
    error_bar = []

    for i in range(0,(bins1.shape[0]-1)):
        for j in range(0,(bins2.shape[0]-1)):

            temp_array1 = array1.copy()
            temp_array2 = array2.copy()
            temp_array3 = array3.copy()

            temp_weights1 = weights_array1.copy()
            temp_weights2 = weights_array2.copy()
            temp_weights3 = weights_array3.copy()


            Selector = log_filtered1 < bins1[i]
            temp_array1[Selector] = np.nan
            temp_array2[Selector] = np.nan
            temp_array3[Selector] = np.nan

            temp_weights1[Selector] = np.nan
            temp_weights2[Selector] =  np.nan
            temp_weights3[Selector] =  np.nan


            Selector =  log_filtered1 > bins1[i+1]
            temp_array1[Selector] = np.nan
            temp_array2[Selector] = np.nan
            temp_array3[Selector] = np.nan

            temp_weights1[Selector] = np.nan
            temp_weights2[Selector] =  np.nan
            temp_weights3[Selector] =  np.nan


            Selector = log_filtered2 < bins2[j]
            temp_array1[Selector] = np.nan
            temp_array2[Selector] = np.nan
            temp_array3[Selector] = np.nan

            temp_weights1[Selector] = np.nan
            temp_weights2[Selector] =  np.nan
            temp_weights3[Selector] =  np.nan


            Selector =  log_filtered2 > bins2[j+1]
            temp_array1[Selector] = np.nan
            temp_array2[Selector] = np.nan
            temp_array3[Selector] = np.nan

            temp_weights1[Selector] = np.nan
            temp_weights2[Selector] =  np.nan
            temp_weights3[Selector] =  np.nan


            ma_array = np.ma.MaskedArray(temp_array1, mask=np.isnan(temp_array1))
            ma_weights = np.ma.MaskedArray(temp_weights1, mask=np.isnan(temp_weights1))
            average_temp = np.ma.average(ma_array, weights=ma_weights)
            bin_centres1.append(average_temp)

            ma_array = np.ma.MaskedArray(temp_array2, mask=np.isnan(temp_array2))
            ma_weights = np.ma.MaskedArray(temp_weights2, mask=np.isnan(temp_weights2))
            average_temp = np.ma.average(ma_array, weights=ma_weights)
            bin_centres2.append(average_temp)

            ma_array = np.ma.MaskedArray(temp_array3, mask=np.isnan(temp_array3))
            ma_weights = np.ma.MaskedArray(temp_weights3, mask=np.isnan(temp_weights3))
            average_temp = np.ma.average(ma_array, weights=ma_weights)
            binned_data.append(average_temp)

            error_bar.append(np.nanstd(temp_array3))

    bin_centres1 = np.array(bin_centres1)
    bin_centres2 = np.array(bin_centres2)
    binned_data = np.array(binned_data)
    error_bar = np.array(error_bar)

    plt.figure()
    plt.scatter(array1, array2,c = array3)
    plt.xscale('log')
    plt.yscale('log')
    plt.colorbar()
    plt.show()

    bin_centres1 = np.log10(np.array(bin_centres1))
    bin_centres2 = np.log10(np.array(bin_centres2))
    binned_data = np.log10(np.array(binned_data))
    valid = ~(np.isnan(bin_centres1) | np.isnan(bin_centres2) | np.isnan(binned_data))
    param,param_cov = curve_fit(DoubleParamFunc,(bin_centres1[valid],bin_centres2[valid]),binned_data[valid])
    modelled_data = DoubleParamFunc((bin_centres1[valid],bin_centres2[valid]), param[0], param[1], param[2])

    # plt.figure()
    # plt.scatter(bin_centres1[valid],bin_centres2[valid],c = binned_data[valid] - modelled_data)
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.colorbar()
    # plt.show()

    plt.figure()
    plt.hist(binned_data[valid] - modelled_data,bins = 20)
    plt.show()
    return param,param_cov


def binning_equal_width_2D_vbootstrap(array1,array2,array3,weights_array1,weights_array2,weights_array3,Nbins1=15,Nbins2=15):
    log_filtered1 = np.log10(array1)
    log_filtered2 = np.log10(array2)
    log_filtered3 = np.log10(array3)

    bins1 = np.linspace(np.nanmin(log_filtered1),np.nanmax(log_filtered1),Nbins1)
    bins2 = np.linspace(np.nanmin(log_filtered2),np.nanmax(log_filtered2),Nbins2)

    bin_centres1 = []
    bin_centres2 = []

    binned_data = []
    error_bar = []

    for i in range(0,(bins1.shape[0]-1)):
        for j in range(0,(bins2.shape[0]-1)):

            temp_array1 = array1.copy()
            temp_array2 = array2.copy()
            temp_array3 = array3.copy()

            temp_weights1 = weights_array1.copy()
            temp_weights2 = weights_array2.copy()
            temp_weights3 = weights_array3.copy()


            Selector = log_filtered1 < bins1[i]
            temp_array1[Selector] = np.nan
            temp_array2[Selector] = np.nan
            temp_array3[Selector] = np.nan

            temp_weights1[Selector] = np.nan
            temp_weights2[Selector] =  np.nan
            temp_weights3[Selector] =  np.nan


            Selector =  log_filtered1 > bins1[i+1]
            temp_array1[Selector] = np.nan
            temp_array2[Selector] = np.nan
            temp_array3[Selector] = np.nan

            temp_weights1[Selector] = np.nan
            temp_weights2[Selector] =  np.nan
            temp_weights3[Selector] =  np.nan


            Selector = log_filtered2 < bins2[j]
            temp_array1[Selector] = np.nan
            temp_array2[Selector] = np.nan
            temp_array3[Selector] = np.nan

            temp_weights1[Selector] = np.nan
            temp_weights2[Selector] =  np.nan
            temp_weights3[Selector] =  np.nan


            Selector =  log_filtered2 > bins2[j+1]
            temp_array1[Selector] = np.nan
            temp_array2[Selector] = np.nan
            temp_array3[Selector] = np.nan

            temp_weights1[Selector] = np.nan
            temp_weights2[Selector] =  np.nan
            temp_weights3[Selector] =  np.nan


            ma_array = np.ma.MaskedArray(temp_array1, mask=np.isnan(temp_array1))
            ma_weights = np.ma.MaskedArray(temp_weights1, mask=np.isnan(temp_weights1))
            average_temp = np.ma.average(ma_array, weights=ma_weights)
            bin_centres1.append(average_temp)

            ma_array = np.ma.MaskedArray(temp_array2, mask=np.isnan(temp_array2))
            ma_weights = np.ma.MaskedArray(temp_weights2, mask=np.isnan(temp_weights2))
            average_temp = np.ma.average(ma_array, weights=ma_weights)
            bin_centres2.append(average_temp)

            ma_array = np.ma.MaskedArray(temp_array3, mask=np.isnan(temp_array3))
            ma_weights = np.ma.MaskedArray(temp_weights3, mask=np.isnan(temp_weights3))
            average_temp = np.ma.average(ma_array, weights=ma_weights)
            binned_data.append(average_temp)

            error_bar.append(np.nanstd(temp_array3))

    bin_centres1 = np.array(bin_centres1)
    bin_centres2 = np.array(bin_centres2)
    binned_data = np.array(binned_data)
    error_bar = np.array(error_bar)



    bin_centres1 = np.log10(np.array(bin_centres1))
    bin_centres2 = np.log10(np.array(bin_centres2))
    binned_data = np.log10(np.array(binned_data))
    valid = ~(np.isnan(bin_centres1) | np.isnan(bin_centres2) | np.isnan(binned_data))
    param,param_cov = curve_fit(Functions.DoubleParamFunc,(bin_centres1[valid],bin_centres2[valid]),binned_data[valid])

    return param[1]


# def binning_datav1(array1,array2,delt_bin):
#     levels = np.linspace(np.nanmin(array1),np.nanmax(array1),delt_bin)
#     levels_centres = (levels[:-1] + levels[1:])/2
#     binned_data = []
#     error_bar = []
#     for i in range(levels.shape[0]-1):
#         temp_array2 = array2.copy()
#         Selector = array1 < levels[i]
#         temp_array2[Selector] = np.nan
#         Selector =  array1 > levels[i+1]
#         temp_array2[Selector] =np.nan
#         binned_data.append(np.nanmean(temp_array2))
#         error_bar.append(np.nanstd(10**temp_array2))
#     binned_data = np.array(binned_data)
#     error_bar = np.array(error_bar)
#     plt.figure(figsize=(12,8))
#     plt.errorbar(10**levels_centres,10**binned_data,yerr=error_bar,c='grey')
#     plt.scatter(10**levels_centres,10**binned_data,s = 40,c='r')
#     plt.plot(10**levels_centres,10**binned_data,c='k')

# def binning_datav2(array1,array2,delt_bin):
#     levels = np.linspace(np.nanmin(array1),np.nanmax(array1),delt_bin)
#     levels_centres = (levels[:-1] + levels[1:])/2
#     binned_data = []
#     error_bar = []
#     for i in range(levels.shape[0]-1):
#         temp_array2 = array2.copy()
#         Selector = array1 < levels[i]
#         temp_array2[Selector] = np.nan
#         Selector =  array1 > levels[i+1]
#         temp_array2[Selector] =np.nan
#         binned_data.append(np.nanmean(temp_array2))
#         error_bar.append(np.nanstd(10**temp_array2))
#     binned_data = np.array(binned_data)
#     error_bar = np.array(error_bar)
#     # plt.figure(figsize=(10,6))
#     plt.errorbar(10**levels_centres,10**binned_data,yerr=error_bar,c='grey')
#     plt.scatter(10**levels_centres,10**binned_data,s = 40,c='r')
#     plt.plot(10**levels_centres,10**binned_data,c='k')

# def binning_datav3(array1,array2,delt_bin):
#     levels = np.linspace(np.nanmin(array1),np.nanmax(array1),delt_bin)
#     levels_centres = (levels[:-1] + levels[1:])/2
#     binned_data = []
#     error_bar = []
#     for i in range(levels.shape[0]-1):
#         temp_array2 = array2.copy()
#         Selector = array1 < levels[i]
#         temp_array2[Selector] = np.nan
#         Selector =  array1 > levels[i+1]
#         temp_array2[Selector] =np.nan
#         binned_data.append(np.nanmean(temp_array2))
#         error_bar.append(np.nanstd(10**temp_array2))
#     level_bins = np.linspace(np.amin(levels_centres),np.amax(levels_centres),10)
#     binned_data = np.array(binned_data)
#     error_bar = np.array(error_bar)
#     levels_centres = np.array(levels_centres)
#     param, PS_param_cov = curve_fit(lin_fit, levels_centres, binned_data,sigma=error_bar)
#     PS_FitFunc = lin_fit(level_bins,param[0],param[1])
#     return param,level_bins,PS_FitFunc



# def binning_equal_data(array1,array2,bins):
#     filtered1,filtered2 = remove_nan(array1,array2)
#     def sortlinked(array1,array2):
#         sort = np.argsort(array1)
#         return array1[sort],array2[sort]
#     bin_centres = []
#     data_binned = []
#     error_bar = []
#     sorted1,sorted2 = sortlinked(filtered1,filtered2)
#     bin_size = int((sorted2.shape[0])/bins)
    
#     for i in range(0,bins+1):
#         bin_centres.append(np.nanmean(sorted1[i*bin_size:i*bin_size+bin_size]))
#         data_binned.append(np.nanmean(sorted2[i*bin_size:i*bin_size+bin_size]))
#         error_bar.append(np.nanstd(10**sorted2))
#     bin_centres = np.array(bin_centres)
#     data_binned = np.array(data_binned)
#     error_bar = np.array(error_bar)
#     plt.errorbar(10**bin_centres,10**data_binned,yerr=error_bar,c='grey')
#     plt.scatter(10**bin_centres,10**data_binned,s = 40,c='r')
#     plt.plot(10**bin_centres,10**data_binned,c='k')

# def binning_equal_data_fits(array1,array2,bins):
#     filtered1,filtered2 = remove_nan(array1,array2)
#     def sortlinked(array1,array2):
#         sort = np.argsort(array1)
#         return array1[sort],array2[sort]
#     if np.nan in array1:
#         print(1)
#     if np.nan in array2:
#         print(2)
#     bin_centres = []
#     data_binned = []
#     error_bar = []


#     sorted1,sorted2 = sortlinked(filtered1,filtered2)
#     bin_size = int((sorted2.shape[0])/bins)
    
#     for i in range(0,bins+1):
#         bin_centres.append(np.nanmean(sorted1[i*bin_size:i*bin_size+bin_size]))
#         data_binned.append(np.nanmean(sorted2[i*bin_size:i*bin_size+bin_size]))
#         error_bar.append(np.nanstd(10**sorted2))
#     bin_centres = np.array(bin_centres)
#     data_binned = np.array(data_binned)
#     error_bar = np.array(error_bar)
#     level_bins = np.linspace(np.amin(bin_centres),np.amax(bin_centres),10)
#     param, PS_param_cov = curve_fit(lin_fit, bin_centres, data_binned,sigma=error_bar)
#     PS_FitFunc = lin_fit(level_bins,param[0],param[1])
#     return param,10**level_bins,10**PS_FitFunc


# def binning_equal_width_fits(array1,array2,Nbins = 10):
#     filtered1,filtered2 = remove_nan(array1,array2)
#     log_filtered1 = np.log10(filtered1)
#     log_filtered2 = np.log10(filtered2)

#     bins = np.linspace(np.nanmin(log_filtered1),np.nanmax(log_filtered1),Nbins)
#     # bin_centres = (bins[:-1] + bins[1:])/2

#     bin_centres = []
#     binned_data = []
#     error_bar = []

#     for i in range(0,(bins.shape[0]-1)):
#         temp_array1 = filtered1.copy()
#         temp_array2 = filtered2.copy()

#         Selector = log_filtered1 < bins[i]
#         temp_array1[Selector] = np.nan
#         temp_array2[Selector] = np.nan
#         Selector =  log_filtered1 > bins[i+1]
#         temp_array1[Selector] = np.nan
#         temp_array2[Selector] = np.nan
#         bin_centres.append(np.nanmean(temp_array1))
#         binned_data.append(np.nanmean(temp_array2))
#         error_bar.append(np.nanstd(temp_array2))


#     bin_centres = np.log10(np.array(bin_centres))
#     binned_data = np.log10(np.array(binned_data))
#     bin_centres,binned_data = remove_nan(bin_centres,binned_data)
#     level_bins = np.linspace(np.amin(bin_centres),np.amax(bin_centres),11)
#     param, PS_param_cov = curve_fit(lin_fit, bin_centres, binned_data)
#     PS_FitFunc = lin_fit(level_bins,param[0],param[1])
#     return param,10**level_bins,10**PS_FitFunc
