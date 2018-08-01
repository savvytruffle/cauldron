'''This program plots an HR diagram with isochrones for a binary star system.

It reads in at least one isochrone file as well as
teff for one star and a surface brightness ratio (both on the command line).

It prints out a calculated temperature ratio and teff for the secondary star.
It also plots the isochrone read in with logg vs. logTeff axes.

Implementation questions:
- what's the best way to propagate the uncertainties?
- what's our plan for choosing an appropriate isochrone?
- in order to plot logg for our stars, we need a mass (keblat is good) and a radius(???)
  
- how confident are we that all our systems really are main sequence stars?
- Goal: be able to plug in KEBLAT results from Matson paper into this (+ other analysis?) and see if 
  we get similar temperatures, radii, etc.

note: we are using 'frat', the keblat flux ratio, as a proxy for the surface brightness ratio

'''

import numpy as np
from scipy.interpolate import interp1d
from sys import argv
import matplotlib.pyplot as plt

starIds = [5285607, 6864859, 6778289, 6449358, 4285087, 6131659, 6781535]                                        # Kepler Input Catalog #
aspcapTeffs = [6495, 6417, 6572, 6237, 5664, 4845, 5749]                                                         # Teff of starId from ASPCAP
Tefferrs = [156, 159, 162, 179, 146, 98, 125]                                                                    # error on the above
fluxRatios = [0.6202419237, 0.7732692053, 0.4621566737, 0.3920984216, 0.9966328203, 0.6478124125, 0.8909515437]  # flux ratios from BF areas
fluxRatioErrs = [0.02710, 0.02062, 0.02931, 0.01441, 0.02241, 0.02824, 0.07272]                                  # error on the above
logg1s = [4.080171, 4.120715, 4.291286, np.nan, 4.487564, 4.813282, 4.355104]                                    # log(g) of star 1 from BF flux ratio
logg1_errs = [0.041974, 0.023629, 0.049574, np.nan, 0.021908, 0.05458, 0.064560]                                 # error on the above
logg2s = [4.1503056, 4.277946, 4.254075, np.nan, 4.457921, 4.380783, 4.508261]                                   # log(g) of star 2 from BF flux ratio
logg2_errs = [0.040514, 0.025099, 0.055558, np.nan, 0.020566, 0.038261, 0.072749]                                # error on the above 

for starId, aspcapTeff, Tefferr, fluxRatio, fluxRatioErr, logg1, logg1_err, logg2, logg2_err in zip(
    starIds, aspcapTeffs, Tefferrs, fluxRatios, fluxRatioErrs, logg1s, logg1_errs, logg2s, logg2_errs):

    print(' ')
    print('Running analysis for star', starId)

    # Put filenames and corresponding labels for the isochrones you want to plot here
    isofiles = ['afep0_age1.txt', 'fehm05afep0_age1.txt']
    labels = ['1 Gyr $Z=0$', '1 Gyr $Z=-0.5$']

    ###5738698 Matson et al test star
    #isofiles = ['fehm048afem2_age2.txt', 'fehm048afem2_age2p5.txt', 'fehm048afem2_age3.txt']
    #labels = ['2 Gyr Z=-0.48', '2.5 Gyr Z=-0.48', '3 Gyr Z=-0.48']

    isochroneLogTeffs = []
    isochroneLogggs = []
    teff2s = []
    teff2errs = []
    for isofile in isofiles:
        isochroneMass, isochroneLogTeff, isochroneLogg, kp = np.loadtxt(isofile, unpack=True, usecols = (1, 2, 3, 13))

        # Read in mass, logteff, magnitude from isochrone file (only for points with logg >= 4.1)
        m = np.where(isochroneLogg >= 4.1)  # this restricts the isochrone points used to those with logg >= 4.1
                                            # in effect, this cuts the isochrone down to only include the main sequence.
        isochroneMass = isochroneMass[m]
        isochroneLogTeff = isochroneLogTeff[m]
        isochroneLogg = isochroneLogg[m]
        kp = kp[m]                          # magnitude in kepler bandpass

        # Calculate two new columns for the isochrone: radius squared and surface brightness
        isochroneRadiusSquared = isochroneMass/(10**isochroneLogg)
        isochroneSb = 10**(-kp/2.5)/isochroneRadiusSquared  # isochrone surface brightness column

        # We want a more fine-grained set of logteff and surface brightness values, so we interpolate
        sb1_grid = interp1d(isochroneLogTeff, isochroneSb)  # this returns a scipy.interp1d object
        #print(np.log10(aspcapTeff))
        #print(sb1_grid.x)
        
        sb1 = sb1_grid(np.log10(aspcapTeff))
        # This estimates sb1, the surface brightness of the primary, which is assumed to correspond to the ASPCAP Teff
        # Because we've interpolated, we can now throw any logTeff we want at it and get a surface brightness back!

        # Now we want to find the location at which the flux ratio matches the isochroneSb / sb1.
        fit = np.argmin(np.abs(isochroneSb/sb1 - fluxRatio))
        
        # Use the 'fit' index to select the Teff from the isochrone and calculate a temperature ratio
        tratio = 10**isochroneLogTeff[fit]/aspcapTeff
    
        # Now that we have the ratio, we can calculate the temperature of the secondary
        teff2 = aspcapTeff * tratio
        
        # Attempting some imprecise error propagation... 
        fluxRatiomax = fluxRatio + fluxRatioErr
        fluxRatiomin = fluxRatio - fluxRatioErr
        fitmax = np.argmin(np.abs(isochroneSb/sb1 - fluxRatiomax))
        fitmin = np.argmin(np.abs(isochroneSb/sb1 - fluxRatiomin))
        tratiomax = 10**isochroneLogTeff[fitmax]/aspcapTeff
        tratiomin = 10**isochroneLogTeff[fitmin]/aspcapTeff
        teff2max = aspcapTeff * tratiomax
        teff2min = aspcapTeff * tratiomin
        teff2err = (teff2max - teff2min)/2
        #print(fluxRatio, fluxRatiomax, fluxRatiomin)
        #print(fit, fitmax, fitmin)
        #print(tratio, tratiomax, tratiomin)

        teff2s.append(teff2)
        teff2errs.append(teff2err)  # values are questionable at best
        isochroneLogTeffs.append(isochroneLogTeff)
        isochroneLogggs.append(isochroneLogg)
    
        # Now, we will find the error in the ratio and the teff of the secondary
        #this = (10**isochroneLogTeff[fit])
        #teff2StD = np.std(this) / (Tefferr)
        #print('teff2StD = ', teff2StD)

        print('Using', isofile, 'and {0} +/- {1} K for star 1,'.format(aspcapTeff, Tefferr))
        print('T2/T1 is {0:.3f} which means teff 2 is {1:.0f} K'.format(tratio, teff2))
        # TODO: calculate and print out appropriate uncertainties for tratio and teff2
        # (this means you'll need to input some uncertainties!!)

    ### Plotting stuff ###
    # NOTE this is still inside the loop over all the stars, but outside the isochrone loop.
    # First plot all the isochrones
    for logteff, logg, label in zip(isochroneLogTeffs, isochroneLogggs, labels):
        plt.plot(logteff, logg, ls=':', label=label)
        #plt.axvline(np.log10(aspcapTeff), color='C2', label='Primary')  # vertical lines
        #plt.axvline(np.log10(teff2s[0]), color='C0', ls=':', label='Secondary')
        #plt.axvline(np.log10(teff2s[1]), color='C1', ls=':', label='Secondary')
    ########################################################################################
    # Now plot points for each star in the binary

    plt.errorbar(np.log10(aspcapTeff), logg1, yerr=logg1_err, color='C2', ls='None', marker='o', label='Primary')
    plt.errorbar(np.log10(teff2s[0]), logg2, yerr=logg2_err, color='C3', ls='None', marker='o', label='Secondary')
    # may want to add xerr for temperature error bars too? Star 1 is easy (ASPCAP) but Star 2 is hard...
    # please note!!! just naively plotting teff2s[0] assumes the very first (zeroth) 
    # isochrone is the right one for calculating teff2!!
    
    plt.ylabel('$\log g$')
    plt.xlabel('$\log T_{\mathrm{eff}}$')
    plt.title(starId)
    plt.legend()
    plt.show()
