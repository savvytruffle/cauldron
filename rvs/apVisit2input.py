from __future__ import print_function
import numpy as np
import apogee.tools.read as apread
from apogee.tools import apStarWavegrid
import apogee.spec.plot as splot
from apogee.spec import continuum
import matplotlib.pyplot as plt
import apogee as apg
from apogee.modelspec import ferre
import csv
'''
This program uses jobovy/apogee to make text files for APOGEE visit spectra and a model spectrum.
The model spectrum should have stellar parameters similar to the target star for use with BF_python.
All the final spectra are continuum normalized.
'''

## define useful variables upfront here so they're easy to change in one place.

KIC = 6781535
plate = 4464
ID = '2M19321788+4216489'
#visit = 3
#modelfileout = 'data/'+str(KIC)+'/modeltest.txt'
#specfileout = 'data/'+str(KIC)+'/obsspecnorm'+str(visit)+'.txt' #For da loop
specfileout = 'data/'+str(KIC)+'/obsspecnormTEST'+'.txt'

# reference for other stars
#(4263, '2M19432016+3957081', ext=1, header=False)[visit] #KIC4851217 (6 Visits):::
#(4263, '2M19390532+4027346', ext=1, header=False)[visit] #KIC5285607 (6 Visits):::
#(4263, '2M19355993+3813561', ext=1, header=False)[visit] #KIC3127817 (6 Visits):::
#(4464, '2M19353513+4149543', ext=1, header=False)[visit] #KIC6449358 (25 Visits):::
#(4263, '2M19373173+4027078', ext=1, header=False)[visit] #KIC5284133 (6 Visits):::
#(4464, '2M19282456+4215080', ext=1, header=False)[visit] #KIC6778289 (25 Visits):::
#(4464, '2M19321788+4216489', ext=1, header=False)[Visit] #KIC6781535 (25 Visits):::

# the three arguments are location ID, MJD, and fiber ID, defining them here is neater!
spec = apread.apVisit(7439, 56763, 207, ext=1, header=False)
specerr = apread.apVisit(7439, 56763, 207, ext=2, header=False)
wave = apread.apVisit(7439, 56763, 207, ext=4, header=False)
header = apread.apVisit(7439, 56763, 207, ext=1, header=True)[1]

weird_format_spec = apread.apVisit(7439, 56763, 207, ext=1, header=True)[0]
weird_format_wave = apread.apVisit(7439, 56763, 207, ext=4, header=True)[0]

#Read the (KICnumber)Visitlist.txt and normalize each spectra in it (hopefully)
#2MassID,PlateID,MJD,Fiber,RA,Dec,ReductionVersion,SN,RV 
Visitlist = csv.reader(open{'data/'+str(KIC)+'/4851217Visitlist'+'.txt','r',delimeter=','))
header = reader.next()
column = col(*reader)
for row in Visitlist
	locID = column[2]
	MJD = column[3]
	fiberID = column[4]

cont = continuum.fitApvisit(spec, specerr, wave) #define continuum
specnorm = spec/cont #normalization is the spectra divided by the continuum

plt.plot(wave, spec)
plt.plot(wave, cont, lw=2, color='r')
plt.show()
plt.plot(wave, specnorm)
plt.show()

realstar = open(specfileout, 'w') 
for wave, flux in zip(wave, specnorm): 
#	if flux > 0 and flux != np.nan: # only print positive fluxes that aren't nan
	print(wave, specnorm, file=realstar)
realstar.close()

'''
visits = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]

### BEGIN LOOP OVER VISITS
for visit in visits:

    # Define outfile for visit spectrum
    specfileout = 'data/'+str(KIC)+'/obsspecnorm'+str(visit)+'.txt'

    # Read in spectra:::[0],[1] is the combined spectra:::[2]-[n] are each visit###
    mystar = apread.apStar(plate, ID, ext=1, header=False)[visit]#KIC6449358

    # Read in error###
    mystarerr = apread.apStar(plate, ID, ext=2, header=False)[visit]

    # Reshape the arrays following the example in jobovy/apogee README
    mystar = np.reshape(mystar, (1, len(mystar)))
    mystarerr = np.reshape(mystarerr, (1, len(mystarerr)))

    ###Define Flux and Wavelength
    fluxdata = mystar[0]
    wavedata = apStarWavegrid()

    ###Option to print Flux and Wavelength
    #print(fluxdata)
    #for item in fluxdata:
    #    print(item)
    #print(wavedata)

    ###Fit the continuum
    fitcontinuum = continuum.fit(mystar, mystarerr, type='aspcap')
    fitcontinuumentry = fitcontinuum[0] 
    #            if aspcap does not work well, it is customizable yo

    ###Divide by the continuum (without dividing by zero)
    fluxnorm = [fluxdataentry/fitcontinuumentry if fitcontinuumentry != 0 else np.nan for 
        fluxdataentry, fitcontinuumentry in zip(fluxdata, fitcontinuum[0])]

    ### fluxnorm is great, but it's still spikey.
    threshold = 1.03 # value above which to cut off spikes
    fluxnorm = [flux if flux < threshold else 1.0 for flux in fluxnorm]

    # Meredith had an idea but it didn't turn out to be useful, I don't think
    #range = np.nanmax(fluxnorm) - np.nanmin(fluxnorm)
    #fluxnorm = [(flux-np.nanmin(fluxnorm))/range for flux in fluxnorm]

    ###Print visit spectra to txt files
    realstar = open(specfileout, 'w') 
    for wave, flux in zip(wavedata, fluxnorm): 
        if flux > 0 and flux != np.nan: # only print positive fluxes that aren't nan
            print(wave, flux, file=realstar)
    realstar.close()

## END LOOP OVER VISIT


###Create an appropriate model spectrum by interpolating with FERRE 
## don't forget to set Teff, logg, Fe/H and specify cooler (GK) or hotter (F) library

#modelspec = ferre.interpolate(4812., 4.5, 0.1, 0., 0., 0.)#KIC4851217
#modelspec = ferre.interpolate(6000., 5., 0.1, 0., 0., 0., lib='F')
#modelspec = ferre.interpolate(7000., 4.9, -1.0, 0., 0., 0., lib='F')#KIC5285607
modelspec = ferre.interpolate(6510., 4.6, 0.2, 0., 0., 0., lib='F')#KIC6449358
#                             (Teff., logg, metals, alphafe, nfe, cfe)

###When in doubt, print it out
#for item in modelspec: 
#    print(item)

###Plot the model on top of the last visit spectrum
plt.plot(wavedata, fluxnorm)
plt.plot(wavedata, modelspec, color='r')
plt.show()

###Print the model spectrum to a txt file
templatedata = open(modelfileout, 'w')
for wave, flux in zip(wavedata, modelspec):
    if flux > 0 and flux != np.nan: # only print positive fluxes that aren't nan
        print(wave, flux, file=templatedata)
templatedata.close()


# Old stuff is below

###Overlay model spectrum 
#print('4263','2M19390532+4027346')
#print(mystar['LOCATION_ID'],mystar['APOGEE_ID'])

###Is the data where we think it is?
#xdummy = np.arange(len(mystar[0]))
#plt.plot(xdummy, mystar[0])
#plt.show()

#print(mystar)
#splot.waveregions(mystar[0], labelLines=False, apStar=True)

#spec is the apogee spectra
#spec, hdr= apread.apStar(4102,'2M21353892+4229507',ext=1)
#                        (Plate,2massID,ext=1)
#mspec is the model spectra interpolated by ferre 
#mspec= ferre.interpolate(4750.,2.5,-0.1,0.1,0.,0.)
#                        (Teff.,logg,metals,alphafe,nfe,cfe.)
#apogee.spec.plot.waveregions(spec)
#apogee.spec.plot.waveregions(mspec)
'''