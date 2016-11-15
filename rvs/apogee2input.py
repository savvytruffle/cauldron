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

## TODO: define useful variables upfront here so they're easy to change in one place.
##       To include filenames, PlateID, APOGEEID, desired Teff/logg/FeH of model, etc.

# to get you started...
KIC = 5285607
visit = 7
modelfileout = 'data/'+str(KIC)+'/modeltest.txt'
specfileout = 'data/'+str(KIC)+'/obsspecnorm'+str(visit)+'.txt'

## TODO: BEGIN LOOP OVER VISIT
# visits = [2, 3, 4, 5, 6, 7]
# for visit in visits:
# ... etc.

###Read in spectra:::[0],[1] is the combined spectra:::[2]-[n] are each visit###

#mystar = apread.apStar(4263, '2M19390532+4027346', ext=1, header=False)[2]
#mystar = apread.apStar(4263, '2M19432016+3957081', ext=1, header=False)[7]#KIC4851217
mystar = apread.apStar(4263, '2M19390532+4027346', ext=1, header=False)[visit]	#KIC5285607
#mystar = apread.apStar(4464, '2M19353513+4149543', ext=1, header=False)[2]	#KIC6449358
#mystar = apread.apStar(4263, '2M19355993+3813561', ext=1, header=False)[2]	#KIC3127817Overlap
#mystar = apread.apStar(4263, '2M19373173+4027078', ext=1, header=False)[2]	#KIC5284133

###Read in error###

#mystarerr = apread.apStar(4263, '2M19390532+4027346', ext=2, header=False)[2]
#mystarerr = apread.apStar(4263, '2M19432016+3957081', ext=2, header=False)[7]#KIC4851217
mystarerr = apread.apStar(4263, '2M19390532+4027346', ext=2, header=False)[visit]	#KIC5285607
#mystarerr = apread.apStar(4464, '2M19353513+4149543', ext=2, header=False)[2]	#KIC6449358
#mystarerr = apread.apStar(4263, '2M19355993+3813561', ext=2, header=False)[2]	#KIC3127817Overlap
#mystarerr = apread.apStar(4263, '2M19373173+4027078', ext=2, header=False)[2]	#KIC5284133

###Reshape the arrays following the example in jobovy/apogee README
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
#			if aspcap does not work well, it is customizable yo

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

## TODO: END LOOP OVER VISIT


###Read in Model Spec 
####::::These numbers are not true values, don't believe in the lies::::####

#modelspec = ferre.interpolate(4812., 4.5, 0.1, 0., 0., 0.)#KIC4851217
#modelspec = ferre.interpolate(6000., 5., 0.1, 0., 0., 0., lib='F')
modelspec = ferre.interpolate(6438., 4.6, 0.1, 0., 0., 0., lib='F')#KIC5285607
#modelspec = ferre.interpolate(6510., 4.8, 0.1, 0., 0., 0.)#KIC3127817
#							 (Teff.,logg,metals,alphafe,nfe,cfe.)

###When in doubt, print it out
#for item in modelspec: 
#	print(item)

plt.plot(wavedata, fluxnorm)
plt.plot(wavedata, modelspec, color='r')
plt.show()
#data = apread.rcsample()
#indx = data['SNR'] > 200.
#data = data[indx]


###Print data to txt files###
realstar = open('data/spectest6.txt', 'w') 
templatedata = open('data/modeltest.txt', 'w')
for wave, flux in zip(wavedata, fluxnorm): 
	if flux > 0 and flux != np.nan: # only print positive fluxes that aren't nan
		print(wave, flux, file=realstar)
for wave, flux in zip(wavedata, modelspec):
    if flux > 0 and flux != np.nan: # only print positive fluxes that aren't nan
        print(wave, flux, file=templatedata)
realstar.close()
templatedata.close()

###Overlay model spectrum 
#print('4263','2M19390532+4027346')
#print(mystar['LOCATION_ID'],mystar['APOGEE_ID'])

###Is the data where we think it is?
#xdummy = np.arange(len(mystar[0]))
#plt.plot(xdummy, mystar[0])
#plt.show()


#print(mystar)
#splot.waveregions(mystar[0], labelLines=False, apStar=True)

#Non Continuum Normalized
#splot.waveregions(data[3512]['LOCATION_ID'],data[3512]['APOGEE_ID'],ext=1,
#                  apStar=True,labelID=data[3512]['APOGEE_ID'],
#                  labelTeff=data[3512]['TEFF'],
#                  labellogg=data[3512]['LOGG'],
#                  labelmetals=data[3512]['METALS'],
#                  labelafe=data[3512]['ALPHAFE'])

#Plot a whole detector 
#splot.detector(data[3512]['LOCATION_ID'],data[3512]['APOGEE_ID'],
#               'blue',ext=1,labelLines=False,
#               labelID=data[3512]['APOGEE_ID'],
#               labelTeff=data[3512]['TEFF'],
#               labellogg=data[3512]['LOGG'],
#               labelmetals=data[3512]['METALS'],
#               labelafe=data[3512]['ALPHAFE'])



#spec is the apogee spectra
#spec, hdr= apread.apStar(4102,'2M21353892+4229507',ext=1)
#						(Plate,2massID,ext=1)
#mspec is the model spectra interpolated by ferre 
#mspec= ferre.interpolate(4750.,2.5,-0.1,0.1,0.,0.)
#						(Teff.,logg,metals,alphafe,nfe,cfe.)
#apogee.spec.plot.waveregions(spec)
#apogee.spec.plot.waveregions(mspec)
 
#print(mspec[0])
#print(mspec[1])
#plt.show()