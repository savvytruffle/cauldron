'''This program plots an HR diagram with isochrones for a single binary star system, in 
   this case KIC 6778289. It prints out a calculated temperature ratio and teff for the 
   secondary star. It also plots the isochrone read in with logg vs. logTeff axes.'''

import numpy as np
from scipy.interpolate import interp1d
from sys import argv
import matplotlib.pyplot as plt

####################################### 6778289 ###########################################

starId = 6778289                                        # Kepler Input Catalog #
aspcapTeff = 6572                                       # Teff of starId from ASPCAP
teff1err = 162  
logTeff1err = [2.209515]                                # error on the above
#logTeff2err = [1.886255]                                # log(err) on Teff2 which is calculated below
fluxRatio = [0.4621566737]                              # flux ratios from BF areas
fluxRatioErr = [0.02931]                                # error on the above
logg1 = [4.291286]                                      # log(g) of star 1 from BF flux ratio
logg1_err = [0.049574]                                  # error on the above
logg2 = [4.254075]                                      # log(g) of star 2 from BF flux ratio
logg2_err = [0.055558]                                  # error on the above 

####################################### 6778289 ###########################################


print(' ')
print('Running analysis for star', starId)

#### Put filenames and corresponding labels for the isochrones you want to plot here ####

isofiles = ['fehm05afep0_age1.txt', 'fehm00afem2_age1.txt', 'fehm07afep0_age1.txt', 'fehm07afep0_age2.txt']
labels = ['1 Gyr $Z=0.56$','1 Gyr $Z=0.6$', '1 Gyr $Z=0.7$', '2 Gyr $Z=0.7$']
isochroneLogTeffs = []
isochroneLogggs = []
teff2s = []
teff2errs = []
for isofile in isofiles:
    isochroneMass, isochroneLogTeff, isochroneLogg, kp = np.loadtxt(isofile, unpack=True, usecols = (1, 2, 3, 13))

#### Read in mass, logteff, magnitude from isochrone file (only for points with logg >= 4.1) ####

    m = np.where(isochroneLogg >= 4.1)    # this restricts the isochrone points used to 
                                          # those with logg >= 4.1 in effect, this cuts the 
                                          # isochrone down to only include the main sequence.
    isochroneMass = isochroneMass[m]
    isochroneLogTeff = isochroneLogTeff[m]
    isochroneLogg = isochroneLogg[m]
    kp = kp[m]                                         # magnitude in kepler bandpass

#### Calculate two new columns for the isochrone: radius squared and surface brightness ###

    isochroneRadiusSquared = isochroneMass/(10**isochroneLogg)
    isochroneSb = 10**(-kp/2.5)/isochroneRadiusSquared # isochrone surface brightness column

#### We want a more fine-grained set of logteff and surface brightness values, so we interpolate ###

    sb1_grid = interp1d(isochroneLogTeff, isochroneSb) # this returns a scipy.interp1d object
    
    #print(np.log10(aspcapTeff))
    #print(sb1_grid.x)
        
    sb1 = sb1_grid(np.log10(aspcapTeff)) # This estimates sb1, the surface brightness of 
                                         # the primary, which is assumed to correspond to
                                         # the ASPCAP Teff, because we've interpolated, 
                                         # we can now throw any logTeff we want at it and 
                                         # get a surface brightness back!

#### Now we want to find the location at which the flux ratio matches the isochroneSb / sb1. ####

    fit = np.argmin(np.abs(isochroneSb/sb1 - fluxRatio))
    
#### Use the 'fit' index to select the Teff from the isochrone and calculate a temperature ratio ####

    tratio = 10**isochroneLogTeff[fit]/aspcapTeff

#### Now that we have the ratio, we can calculate the temperature of the secondary ####
    teff2 = aspcapTeff * tratio
    
    
    
#### Attempting some imprecise error propagation #### 
    fluxRatiomax = np.array(fluxRatio) + np.array(fluxRatioErr)
    fluxRatiomin = np.array(fluxRatio) - np.array(fluxRatioErr)
    fitmax = np.argmax(np.abs(isochroneSb/sb1 - fluxRatiomax))
    fitmin = np.argmin(np.abs(isochroneSb/sb1 - fluxRatiomin))
    isochroneLogTeffs.append(isochroneLogTeff)
    isochroneLogggs.append(isochroneLogg)
    tratiomax = 10**isochroneLogTeff[fitmax]/aspcapTeff
    tratiomin = 10**isochroneLogTeff[fitmin]/aspcapTeff
    
    teff2s.append(teff2)
    teff2max = aspcapTeff * tratiomax
    teff2min = aspcapTeff * tratiomin
    isochroneLogTeffs.append(isochroneLogTeff)
    isochroneLogggs.append(isochroneLogg)
    
    teff1err = 162
    Teff1err = 0.434 * teff1err/aspcapTeff
    teff2err = (teff2max - teff2min)/2
    Teff2err = np.log10(0.434  * teff2err/teff2)
    
    
#### Sanity Check print statements ####

#    print('Teff1_err = ', teff1err, 'Teff2_err = ', (teff2err)) 
#    print(fit, fitmax, fitmin)
#    print(tratio, tratiomax, tratiomin)

    print('Using', isofile, 'and {0} +/- {1} K for star 1,'.format(aspcapTeff, Teff1err))
    print('T2/T1 is {0:.3f} which means teff 2 is {1:.0f} K'.format(tratio, teff2))

#### Plotting stuff NOTE this is still inside the loop over all the stars, but outside the isochrone loop. ####

#### First, plot all the isochrones: ####

for logteff, logg, label in zip(isochroneLogTeffs, isochroneLogggs, labels):
    plt.plot(logteff, logg, ls=':', label=label)

#### Sanity Check print statements ####
   
    #plt.axvline(np.log10(aspcapTeff), color='C2', label='Primary')
    #plt.axvline(np.log10(teff2s[0]), color='C0', ls=':', label='Secondary')
    #plt.axvline(np.log10(teff2s[1]), color='C1', ls=':', label='Secondary')

##########################################################################################

#### Now,  plot points and error bars for each star in the binary ####

plt.errorbar(np.log10(aspcapTeff), logg1, yerr=logg1_err, xerr=Teff1err, color='C2', ls='None', marker='o', label='Primary')
plt.errorbar(np.log10(teff2s[0]), logg2, yerr=logg2_err, color='C9', ls='None', marker='o', label='Secondary')
plt.gca().invert_yaxis()                   # Inverts Y axis (increasing downward)
plt.gca().invert_xaxis()                   # Inverts X axis (increasing to the left)
plt.ylabel('$\log g$')
plt.xlabel('$\log T_{\mathrm{eff}}$')
plt.title(starId)
plt.legend()
plt.savefig('6778289HR.png')
plt.savefig('6778289HR.jpeg')
plt.show()