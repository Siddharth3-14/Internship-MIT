from astropy.io import fits
from astropy.convolution import Ring2DKernel
from scipy.ndimage import generic_filter
import numpy as np
import matplotlib.pyplot as plt

def AngleDifferenceRMS(PolarizationAngleArray):
    # assure correct data format
    LengthInputArray = len(PolarizationAngleArray)
    print(LengthInputArray)
    # assert ((LengthInputArray % 2) > 0), 'input array should contain an uneven number of pixels'
    
    # separate central pixel from rest of data
    PolarizationAngleCentralPixel = PolarizationAngleArray[int(LengthInputArray/2.)]
    PolarizationAngleArray[int(LengthInputArray/2.)] = np.nan
     
    # calculate angle differences and remove wraps
    AngleDifferenceArray = np.array(PolarizationAngleArray) - PolarizationAngleCentralPixel
    Selector = (AngleDifferenceArray > 90.)
    AngleDifferenceArray[Selector] = np.array(AngleDifferenceArray[Selector]) - 180.
    Selector = (AngleDifferenceArray < -90.)
    AngleDifferenceArray[Selector] = np.array(AngleDifferenceArray[Selector]) + 180.
   
    # calculate RMS of angle differences
    RMSAngleDifferences = ( np.nanmean(AngleDifferenceArray**2.) )**0.5
    # return result
    return(RMSAngleDifferences)

def test_func(PolarizationAngleArray):
    LengthInputArray = len(PolarizationAngleArray)
    print(LengthInputArray)
    PolarizationAngleCentralPixel = PolarizationAngleArray[int(LengthInputArray/2.)]
    PolarizationAngleArray[int(LengthInputArray/2.)] = np.nan
    # calculate angle differences and remove wraps
    AngleDifferenceArray = np.array(PolarizationAngleArray) - PolarizationAngleCentralPixel
    Selector = (AngleDifferenceArray > 90.)
    AngleDifferenceArray[Selector] = np.array(AngleDifferenceArray[Selector]) - 180.
    Selector = (AngleDifferenceArray < -90.)
    AngleDifferenceArray[Selector] = np.array(AngleDifferenceArray[Selector]) + 180.
   
    # calculate RMS of angle differences
    RMSAngleDifferences = ( np.nanmean(AngleDifferenceArray**2.) )**0.5
    return (RMSAngleDifferences)

FileIn = fits.open('../FITS_file/OMC_BandE.fits')
MapStokesI = FileIn[0]
MapPolAngle = FileIn[11]
MapPolFlux = FileIn[13]
MapPolFluxError = FileIn[14]

MapPolSNR = MapPolFlux.copy()
MapPolSNR.data[:] = np.nan
MapPolSNR.data = MapPolFlux.data / MapPolFluxError.data
# FileIn.info()




# plt.imshow(MapPolAngle.data, origin='lower')
# plt.colorbar()
# plt.show()



# ## Implement Calculation

# ### Actual Calculation

# This sets up the calculation and produces output.
ApertureMask = Ring2DKernel(2, 5)
plt.imshow(ApertureMask)
plt.show()


BlankedMapPolAngle = MapPolAngle.copy()


# flag data at low SNR
Selector = (MapPolSNR.data < 2.)
BlankedMapPolAngle.data[Selector] = np.nan
# plt.imshow(BlankedMapPolAngle.data)
# plt.show()

# assure that data are within +/- 90 deg
Selector = (BlankedMapPolAngle.data > 90.)
BlankedMapPolAngle.data[Selector] = BlankedMapPolAngle.data[Selector] - 180.

Selector = (BlankedMapPolAngle.data < -90.)
BlankedMapPolAngle.data[Selector] = BlankedMapPolAngle.data[Selector] + 180.



# calculate angle differences
RMSAngleDifferences = generic_filter(BlankedMapPolAngle.data,
                                     AngleDifferenceRMS,
                                     footprint=ApertureMask.array,
                                     mode='constant', cval=np.nan)


# visualize results
plt.imshow(RMSAngleDifferences, origin='lower')
plt.colorbar()
plt.show()

# # output to FITS
# MapRMSAngleDifference = MapPolAngle.copy()
# MapRMSAngleDifference.data = RMSAngleDifferences
# MapRMSAngleDifference.writeto('Dummy.fits')


# At this point the variable `RMSAngleDifferences` contains an array that holds the RMS angle differences between one pixel and its neighbors. It can be written out into a FITS file or similar.

