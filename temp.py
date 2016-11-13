from __future__ import print_function
import numpy as np
import apogee.tools.read as apread
from apogee.tools import apStarWavegrid
import apogee.spec.plot as splot
from apogee.spec import continuum
import matplotlib.pyplot as plt
import apogee as apg
from apogee.modelspec import ferre



###Read in spectra:::[0],[1] is the combined spectra:::[2][n] are each visit
mystar = apread.apStar(4263, '2M19390532+4027346', ext=1, header=False)[3]

###Read in error
mystarerr = apread.apStar(4263, '2M19390532+4027346', ext=2, header=False)[3]

###Reshape the arrays!
mystar = np.reshape(mystar, (1, len(mystar)))
mystarerr = np.reshape(mystarerr, (1, len(mystarerr)))

###Define Flux and Wavelength
fluxdata = mystar[0]
wavedata = apStarWavegrid()

###Show Flux and Wavelength
#print(fluxdata)
#print(wavedata)

###Normalizing to continuum
fitcontinuum = continuum.fit(mystar, mystarerr, type='aspcap')
###is aspcap does not work well, it is customizable yo


###Divide by continuum (done in plot because dividing by zero isnan)
#fluxnorm = fluxdata/fitcontinuum[0]

###Read in Model Spec A4851217
####::::These numbers are not true values, don't believe in the lies::::####
modelspec = ferre.interpolate(4812., 2.5, 0.1, 0., 0., 0.)
#							 (Teff.,logg,metals,alphafe,nfe,cfe.)

###When in doubt, print it out
#for item in modelspec: 
#	print(item)

plt.plot(wavedata, fluxdata/fitcontinuum[0])
plt.plot(wavedata, modelspec, color='r')
plt.show()
#data = apread.rcsample()
#indx = data['SNR'] > 200.
#data = data[indx]


###Print data to txt files 
apstar = file.open('obsspec1.txt', 'w')
templatedata = file.open('modelspec.txt', 'w')
for wave, flux in zip(wavedata, fluxnorm): # need to make fluxnorm without any zero division problems!
    print(wave, flux, file=apstar)
for wave, flux in zip(wavedata, modelspec):
    print(wave, flux, file=templatedata)
apstar.close()
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