'''This program plots an HR diagram with isochrones for a binary star system.

It reads in at least one isochrone file as well as
teff for one star and a surface brightness ratio (both on the command line).

It prints out a calculated temperature ratio and teff for the secondary star.
It also plots the isochrone read in with logg vs. logTeff axes.

Implementation questions:
- what's the best way to propagate the uncertainties?
- what's our plan for choosing an appropriate isochrone?
- in order to plot logg for our stars, we need a mass (keblat is good) and a radius(???)

Science questions:
- look at Matson paper Fig 8 and what they did to choose isochrones and draw conclusions about
  evolutionary histories
  
  ###
  They formed a grid of isochrones over a range of logZ/Z_sun (log of the metallicity in 
  terms of the Sun) values and ages and determined which fit best with a chi squared good-
  ness of fit statistic. 
  ###
  
- how confident are we that all our systems really are main sequence stars?
- Goal: be able to plug in KEBLAT results from Matson paper into this (+ other analysis?) and see if 
  we get similar temperatures, radii, etc.
  
TO USE THIS PROGRAM:

python HR.py starId aspcapTeff Tefferr keblatSbRatio SbRatioError

e.g., python HR.py 123456 5000 200 0.5 0.05

note: we are using 'frat', the keblat flux ratio, as a proxy for the surface brightness ratio

'''

import numpy as np
from scipy.interpolate import interp1d
from sys import argv
import matplotlib.pyplot as plt

# Read starId, aspcapTeff, Tefferr, and keblatSbRatio from the command line
starId = argv[1]
aspcapTeff = float(argv[2])
Tefferr = float(argv[3])
keblatSbRatio = float(argv[4])
SbErr = float(argv[4])

print('Running analysis for star', starId)
print('You entered teff 1 as', aspcapTeff)

# Put filenames and corresponding labels for the isochrones you want to plot here
#isofiles = ['afep0.txt', 'fehm25afem2_age1.txt', 'fehp02afep0_age1.txt']
#labels = ['1 Gyr Z=0', '1 Gyr Z=-0.2', '1 Gyr Z=+0.2']

###5738698 Matson et al test star
isofiles = ['fehm048afem2_age2.txt', 'fehm048afem2_age2p5.txt', 'fehm048afem2_age3.txt']
labels = ['2 Gyr Z=-0.48', '2.5 Gyr Z=-0.48', '3 Gyr Z=-0.48']



teff2StDs = []
isochroneLogTeffs = []
isochroneLogggs = []
teff2s = []
teff2err = []
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
    sb1 = interp1d(isochroneLogTeff, isochroneSb)  # this returns a scipy.interp1d object
    sb1 = sb1(np.log10(aspcapTeff))
    # This estimates sb1, the surface brightness of the primary, which is assumed to correspond to the ASPCAP Teff
    # Because we've interpolated, we can now throw any logTeff we want at it and get a surface brightness back!

    # Now we want to find the location at which Keblat's SB ratio matches the isochroneSb / sb1.
    fit = np.argmin(np.abs(isochroneSb/sb1 - keblatSbRatio))
    
    # We also need a fit for the maxima and minima of the Keblat's SB ratio
    keblatSbRatiomax = keblatSbRatio + SbErr
    keblatSbRatiomin = keblatSbRatio + SbErr
    fitmax = np.argmin(np.abs(isochroneSb/sb1 - keblatSbRatiomax))
    fitmin = np.argmin(np.abs(isochroneSb/sb1 - keblatSbRatiomin))

    # Use the 'fit' index to select the Teff from the isochrone and calculate a temperature ratio
    
    tratio = 10**isochroneLogTeff[fit]/aspcapTeff
    tratiomax = 10**isochroneLogTeff[fitmax]/aspcapTeff
    tratiomin = 10**isochroneLogTeff[fitmin]/aspcapTeff
    
    # Now that we have the ratio, we can calculate the temperature of the secondary
    teff2 = aspcapTeff * tratio
    teff2max = aspcapTeff * tratiomax
    teff2min = aspcapTeff * tratiomin
    
    teff2err = teff2max - teff2
    
    teff2s.append(teff2)
    isochroneLogTeffs.append(isochroneLogTeff)
    isochroneLogggs.append(isochroneLogg)
    
    # Now, we will find the error in the ratio and the teff of the secondary
    #this = (10**isochroneLogTeff[fit])
    #teff2StD = np.std(this) / (Tefferr)
    #print('teff2StD = ', teff2StD)
    
    print('Using the isochrone file', isofile, 'for star 1,')
    print('We calculate the temperature ratio T2/T1 is', tratio)
    print('Which means teff 2 is', teff2, '+/-', teff2err)
    # TODO: calculate and print out appropriate uncertainties for tratio and teff2
    # (this means you'll need to input some uncertainties!!)

### Plotting stuff ###

# First plot all the isochrones
for logteff, logg, label in zip(isochroneLogTeffs, isochroneLogggs, labels):
    plt.plot(logteff, logg, label=label)
    #plt.axvline(np.log10(aspcapTeff), color='C1', label='Primary')
    #plt.axvline(np.log10(teff2), color='C2', label='Secondary')

# Now plot points for each star in the binary
# TODO: actually estimate logg because these hard-wired values are based on keblat radii and
#       don't have any error bars!

# 5285607
#plt.plot(np.log10(aspcapTeff), 3.925, color='C1', ls='None', marker='o', label='Primary')
#plt.plot(np.log10(teff2s[0]), 4.375, color='C2', ls='None', marker='o', label='Secondary')
#plt.errorbar(np.log10(teff2s[0]), 4.375, color='C2', ls='None', marker='o', label='Secondary')

# 6864859
#plt.plot(np.log10(aspcapTeff), 4.250, color='C1', ls='None', marker='o', label='Primary')
#plt.plot(np.log10(teff2s[0]), 4.149, color='C2', ls='None', marker='o', label='Secondary')

# 6778289
#plt.plot(np.log10(aspcapTeff), 4.132, color='C1', ls='None', marker='o', label='Primary')
#plt.plot(np.log10(teff2s[0]), 4.478, color='C2', ls='None', marker='o', label='Secondary')

# 4285087
#plt.plot(np.log10(aspcapTeff), 4.465, color='C1', ls='None', marker='o', label='Primary')
#plt.plot(np.log10(teff2s[0]), 4.479, color='C2', ls='None', marker='o', label='Secondary')

# 6131659
#plt.plot(np.log10(aspcapTeff), 4.496, color='C1', ls='None', marker='o', label='Primary')
#plt.plot(np.log10(teff2s[0]), 4.704, color='C2', ls='None', marker='o', label='Secondary')
#plt.plot(np.log10(teff2err[0]), 4.704, color='C3', ls='None', marker='o')

# 6781535
plt.plot(np.log10(aspcapTeff), 4.465, color='C1', ls='None', marker='o', label='Primary')
plt.plot(np.log10(teff2s[0]), 4.479, color='C2', ls='None', marker='o', label='Secondary')
plt.errorbar(np.log10(teff2s[0]), 4.479, yerr = None, xerr = None, ecolor='C2')

#plt.grid()
plt.ylabel('$\log g$')
plt.xlabel('$\log T_{\mathrm{eff}}$')
plt.title(starId)
#plt.axis('tight')
plt.legend()
plt.show()
#plt.savefig('5738698fehm0.48HR.png')
