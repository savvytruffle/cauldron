#from __future__ import print_function # good idea for python < 3
#from __future__ import division # good idea for python < 3
import copy
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from apogee.spec.continuum import _fit_aspcap
'''
reads in a PHOENIX model spectrum (FITS format)
makes a continuum flattened (normalized) spectrum over the APOGEE wavelength range
saves the result to a text file

'specfile' is the FITS HiRes PHOENIX spectrum file which lives in directory 'dir'
'outtxt' is the new text file that will be written

you can get PHOENIX spectra from here: http://phoenix.astro.physik.uni-goettingen.de
(click "version 2.0 of the spectral library" and login as a guest when prompted)
'''

def fitPhoenix(spec, specerr, wave, deg=4, niter=10, usigma=3., lsigma=0.1, cont_pixels=None):
    """
    Continuum fitting routine for apVisit spectra (one spectrum at a time; aspcap method only)
    INPUT:
       spec - single spectrum to fit
       specerr - error on the spectrum; assume no covariances
       wave - wavelength grid corresponding to spec and specerr; must have length 12288
       keywords:
          deg = (4) degree of the polynomial
          niter = (10) number of sigma-clipping iterations to perform
          usigma, lsigma = (3., 0.1) upper and lower sigmas for sigma clipping
    OUTPUT:
       continuum
    Added by Meredith Rawls, 2016-11
       TODO: -Generalize to work for wavelength grid of any length
             -Allow multiple apVisits to be continuum-normalized at once, like the regular 'fit'
    """
    # Parse the input
    tspec = copy.copy(spec)
    tspecerr = specerr
    cont = np.empty_like(tspec)
    cont = _fit_aspcap(wave, tspec, tspecerr, deg, niter, usigma, lsigma)
    return cont


### MAIN PROGRAM BEGINS HERE ###

#######################
## edit values below ##
#######################
# directory where PHOENIX wavelength file lives:
wavedir = 'data/PHOENIX/'
# directory where PHOENIX spectra live:
specdir = 'data/PHOENIX/Z-0.0/'
# spectrum filename:
specfile = 'lte05500-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits' #example
#specfile = 'lte06900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits' #4851217
#specfile = 'lte06400-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits' #5285607
#specfile = 'lte06500-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits' #6449358
#specfile = 'lte05700-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits' #6781535
#specfile = 'lte04200-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits' #5284133
#specfile = 'lte06300-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits' #6778289
#specfile = 'lte06800-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits' #6864859

#######################
## edit values above ##
#######################

# the following are defined based on the names above and the fact that we're running this for apogee
wavefile = 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'
wavepath = wavedir + wavefile
specpath = specdir + specfile
outtxt = specfile[:-5] + '_apogee_norm.txt'
wavestart = 15100
waveend = 17000

# read in PHOENIX FITS spectrum and corresponding wavelength FITS file
hdu = fits.open(specpath)
spec = hdu[0].data
hduwave = fits.open(wavepath)
wave = hduwave[0].data

# put the wavelength array in goddamn air units
# reference: Huswer et al. 2013
wave = wave / (1.0 + 0.05792105/(238.0185-(1e4/wave)**2) + 0.00167917/(57.362-(1e4/wave)**2))

# truncate spectrum to the range you want to fit
idxstart = np.where(wave > wavestart)[0][0]
idxend = np.where(wave > waveend)[0][0]
wavechunk = wave[idxstart:idxend]
specchunk = spec[idxstart:idxend]
errdummy = np.zeros(len(specchunk))
errdummy.fill(1.0)

# measure the continuum level and divide by it to flatten the spectrum
print('Flattening spectrum...')
contspec = fitPhoenix(specchunk, errdummy, wavechunk, deg=3, niter=20, usigma=3., lsigma=0.1, cont_pixels=None)
specnorm = specchunk/contspec

# plot the normalization and the final spectrum
plt.plot(wavechunk, specchunk)
plt.plot(wavechunk, contspec, lw=2, color='r')
plt.show()
plt.plot(wavechunk, specnorm)
plt.show()

# write the result to a text file
f = open(outtxt, 'w')
for wentry, sentry in zip(wavechunk, specnorm):
    print(wentry, sentry, file=f)
f.close()
print('New spectrum written to {0}'.format(outtxt))