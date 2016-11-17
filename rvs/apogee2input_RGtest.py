from __future__ import print_function
import numpy as np
import apogee.tools.read as apread
from apogee.tools import apStarWavegrid
import apogee.spec.plot as splot
from apogee.spec import continuum
import matplotlib.pyplot as plt
import apogee as apg
from apogee.modelspec import ferre
'''
This program uses jobovy/apogee to make text files for APOGEE visit spectra and a model spectrum.
The model spectrum should have stellar parameters similar to the target star for use with BF_python.
All the final spectra are continuum normalized.
'''

# define useful variables upfront here so they're easy to change in one place.
KIC = 'RGtest'
plate = 4464
ID = '2M19315429+4232516'
Nvisits = 27
Teff = 4516
logg = 2.50
FeH = -0.34
lib = 'GK'
modelfileout = 'data/' + str(KIC) + '/model_' + str(Teff) + '_' + str(logg) + '_' + str(FeH) + '.txt'

visits = [N for N in range(2, 2+Nvisits)]
for visit in visits:

    # Define outfile for visit spectrum
    specfileout = 'data/' + str(KIC) + '/obsspecnorm' + str(visit) + '.txt'

    # Read in spectra:::[0],[1] is the combined spectra:::[2]-[n] are each visit###
    mystar = apread.apStar(plate, ID, ext=1, header=False)[visit]#KIC6449358

    # Read in error###
    mystarerr = apread.apStar(plate, ID, ext=2, header=False)[visit]

    # Reshape the arrays following the example in jobovy/apogee README
    mystar = np.reshape(mystar, (1, len(mystar)))
    mystarerr = np.reshape(mystarerr, (1, len(mystarerr)))

    # Define Flux and Wavelength
    fluxdata = mystar[0]
    wavedata = apStarWavegrid()

    # Option to print Flux and Wavelength
    #print(fluxdata)
    #for item in fluxdata:
    #    print(item)
    #print(wavedata)

    # Fit the continuum
    fitcontinuum = continuum.fit(mystar, mystarerr, type='aspcap')
    fitcontinuumentry = fitcontinuum[0] 
    #            if aspcap does not work well, it is customizable yo

    # Divide by the continuum (without dividing by zero)
    fluxnorm = [fluxdataentry/fitcontinuumentry if fitcontinuumentry != 0 else np.nan for 
        fluxdataentry, fitcontinuumentry in zip(fluxdata, fitcontinuum[0])]

    # fluxnorm is great, but it's still spikey - fix that
    threshold = 1.03 # value above which to cut off spikes
    fluxnorm = [flux if flux < threshold else 1.0 for flux in fluxnorm]

    # Print visit spectra to txt files
    realstar = open(specfileout, 'w') 
    for wave, flux in zip(wavedata, fluxnorm): 
        if flux > 0 and flux != np.nan: # only print positive fluxes that aren't nan
            print(wave, flux, file=realstar)
    realstar.close()

# Create an appropriate model spectrum by interpolating with FERRE 
modelspec = ferre.interpolate(Teff,  logg, FeH,  0.,      0.,  0., lib=lib)
#                            (Teff., logg, Fe/H, alphafe, nfe, cfe, 'F' or 'GK')

if -1000 in modelspec:
    print('There was a problem making the model. Ensure parameters are fall inside grid boundaries.')
    
else:
    # Print the model spectrum to a txt file
    templatedata = open(modelfileout, 'w')
    for wave, flux in zip(wavedata, modelspec):
        if flux > 0 and flux != np.nan: # only print positive fluxes that aren't nan
            print(wave, flux, file=templatedata)
    templatedata.close()

# Plot the model on top of the last visit spectrum
plt.plot(wavedata, fluxnorm)
plt.plot(wavedata, modelspec, color='r')
plt.show()

