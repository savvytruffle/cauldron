from __future__ import print_function
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from PyAstronomy import pyasl

# read in wavelengths and fluxes for a list of APOGEE spectra
# this currently works for apVisit spectra only!!! (not apStar)

wavelist = []
speclist = []
#filelist = '/Users/revhalzoo/SDSS/A4851217/A4851217infiles.txt'
#filelist = '/Users/revhalzoo/SDSS/B5285607/B5285607infiles.txt'
#filelist = '/Users/revhalzoo/SDSS/C6449358/C6449358infiles.txt'
#filelist = '/Users/revhalzoo/SDSS/DCA/DChoAinfiles.txt'
#filelist = '/Users/revhalzoo/SDSS/A4851217/A4851217phx.txt'
filelist = 'data/4851217/A4851217infiles.txt'
f1 = open(filelist)
infilelist = [] # for use later to make outfilelist

# option to just turn the fits file into a text file without any despiking
doDespike = True
isFits = False

if isFits == True: #it's a FITS file
    for line in f1:
        infile = line.rstrip()
        infilelist.append(infile)
        print (infile)
        hdu = fits.open(infile)
        # APOGEE! the data is in a funny place and backwards
        wave = hdu[4].data
        wave = wave.flatten()
        wave = wave[::-1]
        spec = hdu[1].data
        spec = spec.flatten()
        spec = spec[::-1]
        spec = spec / np.median(spec) # put the continuum roughly near 1
        wavelist.append(wave)
        speclist.append(spec)
    f1.close()

else: # it's a text file
    for line in f1:
        infile = line.rstrip()
        infilelist.append(infile)
        wave, spec = np.loadtxt(infile, usecols=(0,1), unpack=True)
        wavelist.append(wave)
        speclist.append(spec)
    f1.close()

# do something to 'spec' to cut out spikes, call it newspec
newwavelist = []
newspeclist = []
for wave, spec in zip(wavelist, speclist):
    if doDespike == True:
        r = pyasl.pointDistGESD(spec, maxOLs=1000, alpha=5000)
        # r[0] is number of outliers found, r[i] is indices of outliers
        # maxOLs is max number of outliers that may be identified; increase alpha to find more
        #print(r[0], 'outliers found')
        newwave, newspec = np.delete(wave, r[1]), np.delete(spec, r[1])
    else:
        newspec = spec
        newwave = wave
    # option to plot the result
    plt.plot(wave, spec)
    plt.plot(newwave, newspec, color='r')
    plt.show()
    newwavelist.append(newwave)
    newspeclist.append(newspec)
    
# write out a set of two-column text files,
# each containing one element of wavelist and one element of newspeclist
for file, wave, newspec in zip(infilelist, wavelist, newspeclist):
    # create outfile based on infile name
    outfile = file[0:-4]+'_despiked.txt'
    print(outfile)
    f = open(outfile, 'w')
    for wentry, sentry in zip(wave, newspec):
        print(wentry, sentry, file=f)
