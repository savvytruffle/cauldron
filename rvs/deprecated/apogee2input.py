from __future__ import print_function
import numpy as np
import apogee.tools.read as apread
from apogee.tools import apStarWavegrid
import apogee.spec.plot as splot
from apogee.spec import continuum
import matplotlib.pyplot as plt
import apogee as apg
from apogee.modelspec import ferre
import os
from astropy.io import fits
import csv
'''
This program uses jobovy/apogee to make text files for APOGEE visit spectra and a model spectrum.
The model spectrum should have stellar parameters similar to the target star for use with BF_python.
All the final spectra are continuum normalized.
'''
### define useful variables upfront here so they're easy to change in one place. ###
KIC = 3848919
#ApogeeID = '2M19432016+3957081'

locIDs, mjds, fiberIDs = np.loadtxt('data/' + str(KIC) +'/' + str(KIC) + 'Visitlist.txt', 
    usecols=(1, 2, 3), unpack=True, delimiter=',')

infilelist = []; HJDlist = []; BCVlist = []

for locID, mjd, fiberID in zip(locIDs, mjds, fiberIDs):
    # the three arguments are location ID, MJD, and fiber ID, and they need to be integers
    locID = int(locID)
    mjd = int(mjd)
    fiberID = int(fiberID)
#    print('locID, mjd, fiberID: ', locID, mjd, fiberID) # test to make sure the values are correct
    
### Trying to fix the issue with spectra not being found on the server ###
    locID = int(locID)
    fitsfilepath = 'data/'+str(KIC)+'/apVisit-r5-'+ str(locID) + '-' + str(mjd) + '-' + str(fiberID) + '.fits'
    header = fits.open(fitsfilepath)[0].header
    #HJDlist.append(float('24'+str(header['HJD'])))
    #BCVlist.append(header['BC'])

### Make a unique outfile for each visit spectrum ###
    specfileout = 'data/'+str(KIC)+'/apVisitnorm-'+str(locID)+'-'+str(mjd)+'-'+str(fiberID)+'.txt'
    infilelist.append(specfileout)

    fluxes = apread.apVisit(locID, mjd, fiberID, ext=1, header=False)
    fluxerrs = apread.apVisit(locID, mjd, fiberID, ext=2, header=False)
    waves = apread.apVisit(locID, mjd, fiberID, ext=4, header=False)
    #header = apread.apVisit(locID, mjd, fiberID, ext=1, header=True) # only gives 1st extension header data
    
    # Manually access the file you just downloaded and read the header from it (ugh, sorry)
    #SDSSland = os.environ.get('SDSS_LOCAL_SAS_MIRROR')
    #fitsfilepath = (SDSSland + '/dr12/apogee/spectro/redux/r5/apo25m/' + str(locID) + '/' + str(mjd) + 
    #               '/apVisit-r5-'+ str(locID) + '-' + str(mjd) + '-' + str(fiberID) + '.fits'
    
### Trying to fix the issue with spectra not being found on the server ###
    fitsfilepath = 'data/'+str(KIC)+'/apVisit-r5-'+ str(locID) + '-' + str(mjd) + '-' + str(fiberID) + '.fits'
    header = fits.open(fitsfilepath)[0].header
    HJDlist.append(float('24'+str(header['HJD'])))
    BCVlist.append(header['BC'])

### This part normalizes the spectrum ###
    contspec = continuum.fitApvisit(fluxes, fluxerrs, waves) #define continuum
    specnorm = fluxes/contspec #normalization is the spectra divided by the continuum

    # Plot visit spectrum with continuum fit to make sure it looks OK
    plt.plot(waves, fluxes)
    plt.plot(waves, contspec, lw=2, color='r')
    plt.show()
    plt.plot(waves, specnorm)
    plt.show()

### Save visit spectrum to text file ###
    realstar = open(specfileout, 'w')
    for wave, flux in zip(waves, specnorm): 
        print(wave, flux, file=realstar)
    realstar.close()
    print('Visit spectrum saved to {0}'.format(specfileout))

for file, HJD, BCV in zip(infilelist, HJDlist, BCVlist):
    print(file, HJD, BCV)
