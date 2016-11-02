import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# Point to the files
fitsfile = '../../../joni_EBs/apStar-r5-2M18413657-0334237-rvstd.fits'
outfile = 'template_test.txt'

# Get data out of the FITs file
hdu = fits.open(fitsfile)
head = hdu[0].header
fluxes = hdu[1].data[0] # special case for apStar! this gives us the combined spectrum (zeroth array)

# Normalize the fluxes in the easiest way possible
fluxes = fluxes / np.median(fluxes)

# Make the wavelength scale from header info
# (it is log10 by default!)
headerdwave = head['CDELT1']
headerwavestart = head['CRVAL1']
headerwavestop = headerwavestart + headerdwave*len(fluxes)
logwaves = np.arange(headerwavestart, headerwavestop, headerdwave)
waves = np.power(10, logwaves)

# Make sure you have what you think you have
#print(waves)
#print(fluxes)
#print(len(waves), len(fluxes))

# Write wavelength, flux out to a file
fout = open(outfile, 'w')
for wave, flux in zip(waves, fluxes):
    print(wave, flux, file=fout)

# Make a lovely plot
plt.plot(waves, fluxes)
plt.xlabel('Wavelength (A)')
plt.ylabel('Flux')
plt.show()