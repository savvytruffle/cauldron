'''This program plots an HR diagram with isochrones for a single binary star system, in 
   this case KIC 5285607. It prints out a calculated temperature ratio and teff for the 
   secondary star. It also plots the isochrone read in with logg vs. logTeff axes.'''

import numpy as np
from scipy.interpolate import interp1d
from sys import argv
import matplotlib.pyplot as plt

####################################### 5285607 ###########################################

starId = 5285607                                        # Kepler Input Catalog 
aspcapTeff = 6495                                       # Teff (primary) of starId from ASPCAP
teff1err = 156                                          # Teff (primary) error from ASPCAP
logTeff1err = [2.193125]                                # log_10 of Teff error
logTeff2err = [1.78049]                                 # log(err) on Teff2 (calculated below)
BFfluxRatio = [0.6202419237]                            # flux ratios from BF areas
BFfluxRatioErr = [0.02710]                              # error on the above
KEBLATfluxratio = [0.258]                               # KEBLAT flux ratio *NOT A SURFACE BRIGHTNESS*
KEBLATfluxratio_err = [0.046]                           # error on the above
logg1 = [4.080171]                                      # log(g) of star 1 from BF flux ratio
logg1_err = [0.041974]                                  # error on the above
logg2 = [4.1503056]                                     # log(g) of star 2 from BF flux ratio
logg2_err = [0.040514]                                  # error on the above 

print(' ')
print('Running analysis for star', starId)

######################################## 5285607 #########################################

#### Put filenames and corresponding labels for the isochrones you want to plot here ####

isofiles = ['isochrones/fehp00afep4_age1.txt', 'isochrones/fehm00afep8_age1.txt', 
    'isochrones/fehp05afep2_age1.txt']
labels = ['1 Gyr $Z=0.01$', '1 Gyr $Z=0.09$', '1 Gyr $Z=0.57$']

isochroneLogTeffs = []
isochroneLogggs = []
teff2s = []
teff2errs = []
for isofile in isofiles:
    isochroneMass, isochroneLogTeff, isochroneLogg, kp = np.loadtxt(isofile, unpack=True, usecols = (1, 2, 3, 13))

#### Read in mass, logteff, magnitude from isochrone file (only for points with logg >= 4.1) ####

    m = np.where(isochroneLogg >= 4.1)    ### this restricts the isochrone points used to 
                                          ### those with logg >= 4.1 in effect, this cuts the 
                                          ### isochrone down to only include the main sequence.
    isochroneMass = isochroneMass[m]
    isochroneLogTeff = isochroneLogTeff[m]
    isochroneLogg = isochroneLogg[m]
    kp = kp[m]                                         # magnitude in kepler bandpass

#### Calculate two new columns for the isochrone: radius squared and surface brightness ###

    isochroneRadiusSquared = isochroneMass/(10**isochroneLogg)
    isochroneSb = 10**(-kp/2.5)/isochroneRadiusSquared # isochrone surface brightness column

#### We want a more fine-grained set of logteff and surface brightness values, so we interpolate ###

    sb1_grid = interp1d(isochroneLogTeff, isochroneSb) ### this returns a scipy.interp1d object
    
    #print(np.log10(aspcapTeff))
    #print(sb1_grid.x)
        
    sb1 = sb1_grid(np.log10(aspcapTeff)) ### This estimates sb1, the surface brightness of 
                                         ### the primary, which is assumed to correspond to
                                         ### the ASPCAP Teff, because we've interpolated, 
                                         ### we can now throw any logTeff we want at it and 
                                         ### get a surface brightness back!

#### Now we want to find the location at which the flux ratio matches the isochroneSb / sb1. ####

    fit = np.argmin(np.abs(isochroneSb/sb1 - KEBLATfluxratio)) #### Using the KEBLAT SB Ratio
    
#### Use the 'fit' index to select the Teff from the isochrone and calculate a temperature ratio ####

    tratio = 10**isochroneLogTeff[fit]/aspcapTeff  ####ONLY from KEBLAT and isochrones 

#### Now that we have the ratio, we can calculate the temperature of the secondary ####
    teff2 = aspcapTeff * tratio
    
#### Attempting some imprecise error propagation #### 

    

    KEBLATfluxratiomax = np.array(KEBLATfluxratio) + np.array(KEBLATfluxratio_err)
    KEBLATfluxratiomin = np.array(KEBLATfluxratio) - np.array(KEBLATfluxratio_err)
    fitmax = np.argmax(np.abs(isochroneSb/sb1 - KEBLATfluxratiomax))
    fitmin = np.argmin(np.abs(isochroneSb/sb1 - KEBLATfluxratiomin))
    isochroneLogTeffs.append(isochroneLogTeff)
    isochroneLogggs.append(isochroneLogg)
    tratiomax = 10**isochroneLogTeff[fitmax]/aspcapTeff
    tratiomin = 10**isochroneLogTeff[fitmin]/aspcapTeff

    
    teff2s.append(teff2)
    teff2max = aspcapTeff * tratiomax
    teff2min = aspcapTeff * tratiomin
    isochroneLogTeffs.append(isochroneLogTeff)
    isochroneLogggs.append(isochroneLogg)
    
    Teff1err = 0.434 * teff1err/aspcapTeff    
    teff2err = (teff2max - teff2min)/2
    Teff2err = np.log10(0.434  * teff2err/teff2)
    
    
#### Sanity Check print statements ####

    
#    print(fit, fitmax, fitmin)
#    print(tratio, tratiomax, tratiomin)

    print('Using', isofile, 'and {0} +/- {1} K for star 1,'.format(aspcapTeff, Teff1err))
    print('T2/T1 is {0:.3f} which means teff 2 is {1:.0f} K'.format(tratio, teff2))
    print('Teff1_err = ', teff1err, 'Teff2_err = ', (teff2err)) 
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
#plt.errorbar(np.log10(teff2s[0]), logg2, yerr=logg2_err, xerr=Teff2err, color='C3', ls='None', marker='o', label='Secondary')
plt.errorbar(np.log10(teff2s[0]), logg2, yerr=logg2_err, color='C3', ls='None', marker='o', label='Secondary')
plt.gca().invert_yaxis()                   # Inverts Y axis (increasing downward)
plt.gca().invert_xaxis()                   # Inverts X axis (increasing to the left)
plt.ylabel('$\log g$')
plt.xlabel('$\log T_{\mathrm{eff}}$')
plt.title(starId)
plt.legend()
plt.savefig('figures/HRDiagrams/5286507HRK.png')
plt.savefig('figures/HRDiagrams/5286507HRK.jpeg')
plt.show()