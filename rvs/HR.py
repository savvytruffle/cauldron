'''This program does HR diagram stuff for a binary star system.

It reads in an isochrone file (hardcoded for now) as well as
teff for one star and a surface brightness ratio (both on the command line).

It prints out a calculated temperature ratio and teff for the secondary star.
It also plots the isochrone read in with logg vs. logTeff axes.

Outstanding questions:
- what to do once we have temperature ratio and teff2
- what motivates us to choose one isochrone over another (why 1Gyr? why z=0?)
- wasn't a goal to better constrain radii? How does that get us this?
    --> maybe because we can now read logg from the y-axis and plus in our known mass and get R?
        but is that really better?
- look at Matson paper Fig 8 and what they did to choose isochrones and draw conclusions about
  evolutionary histories
- how confident are we that all our systems really are main sequence stars?
- Goal: be able to plug in KEBLAT results from Matson paper into this (+ other analysis?) and see if 
  we get similar temperatures, radii, etc. How does Keivan envision this happening?

'''

import numpy as np
from scipy.interpolate import interp1d
from sys import argv
import matplotlib.pyplot as plt

# Read in the Dartmouth isochrone file
# For now, default to age of 1 Gyr and metallicity of zero
# TODO: generalize to more ages and metallicities by reading in different files
isofile = 'afep0.txt'
isochroneMass, isochroneLogTeff, isochroneLogg, kp = np.loadtxt(isofile, unpack=True, usecols = (1, 2, 3, 13))
starId = argv[1]
aspcapTeff = float(argv[2])
keblatSbRatio = float(argv[3])

# Read in mass, logteff, magnitude from isochrone file (only for points with logg >= 4.1)
m = np.where(isochroneLogg >= 4.1)  # this restricts the isochrone points used to those with logg >= 4.1
                           # in effect, this cuts the isochrone down to only include the main sequence.
isochroneMass = isochroneMass[m]
isochroneLogTeff = isochroneLogTeff[m]
isochroneLogg = isochroneLogg[m]
kp = kp[m]  # magnitude in kepler bandpass

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

# Use the 'fit' index to select the Teff from the isochrone and calculate a temperature ratio
tratio = 10**isochroneLogTeff[fit]/aspcapTeff

# Now that we have the ratio, we can calculate the temperature of the secondary
teff2 = aspcapTeff * tratio

print('Running analysis for star', starId)
print('You entered teff 1 as', aspcapTeff)
print('Using the isochrone file', isofile, 'for star 1,')
print('We calculate the temperature ratio T2/T1 is', tratio)
print('Which means teff 2 is', teff2)
print('Here\'s a plot of the isochrone')

# Plotting stuff
#log(g) vs. log(teff) complete with Dartmouth evolutionary track
logt1 = np.log10(aspcapTeff)
logt2 = np.log10(teff2)
plt.plot(isochroneLogTeff, isochroneLogg, label='Isochrone used for Primary')
plt.axvline(np.log10(aspcapTeff), color='C1', label='Primary')
plt.axvline(np.log10(teff2), color='C2', label='Secondary')
#plt.grid()
plt.ylabel('$\log g$')
plt.xlabel('$\log T_{\mathrm{eff}}$')
plt.title(starId)
#plt.axis('tight')
plt.legend()
plt.show()
