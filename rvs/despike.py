from __future__ import print_function
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from PyAstronomy import pyasl
import argparse
import os
'''
Removes spikes (tellurics) from continuum-normalized apVisit APOGEE spectra.

Typically you will run this after apVisit2input.py, which finds the apVisit
spectra downloaded via the python apogee module for one target and 
continuum-normalizes them. If you somehow have an apVisit-style FITS file
that is continuum-normalized already (or close), this can despike it too.

First, read_infiles to reads in wavelengths and fluxes for a list of spectra
Then, despike_spectra despike the spectra with one of two techniques

A plot will pop up for each spectrum so you can inspect the before/after despiking
(blue/red, respectivly) and the threshold levels selected (green dotted lines)

The new, despiked spectra are written to two-column (wavelength, flux) files
with similar names as the original, but they now end in '_despiked.txt'.

USAGE:
python despike.py -d datapath -i infilelist

datapath: /path/to/infilelist/and/files/listed/in/infilelist

infilelist: single-column text file with list of *continuum-normalized* single-visit
            APOGEE files you want to despike
            the files in this list can either be two-column text files (wave, flux)
            or FITS files in the apVisit format
'''


def main():
    '''
    Read a list of files, do despiking, and write the result to outfiles
    '''
    infilelist, wavelist, speclist = read_infiles(datapath, filelist)
    newwavelist, newspeclist = despike_spectra(wavelist, speclist)

    # write out a set of two-column text files,
    # each containing one element of newwavelist and one element of newspeclist
    for file, newwave, newspec in zip(infilelist, newwavelist, newspeclist):
        # create outfile based on infile name
        outfile = os.path.splitext(file)[0] + '_despiked.txt'
        with open (outfile, 'w') as f:
            for wentry, sentry in zip(newwave, newspec):
                print(wentry, sentry, file=f)
    return


def read_infiles(datapath, filelist, isFits=False):
    '''
    Read in filelist, a text file containing a list of continuum-normalized spectra
    
    The list file and the spectrum files must both live in datapath
    '''
    print(datapath, filelist)
    wavelist = []; speclist = []
    with open(os.path.join(datapath, filelist)) as f1:
        infilelist = [] # for use later to make outfilelist
        if isFits: #it's a FITS file
            for line in f1:
                infile = line.rstrip()
                infile = os.path.join(datapath, infile)
                infilelist.append(infile)
                with fits.open(infile) as hdu:
                    # APOGEE! the data is in a funny place and backwards
                    wave = hdu[4].data.flatten()
                    wave = wave[::-1]
                    spec = hdu[1].data
                spec = spec.flatten()
                spec = spec[::-1]
                spec = spec / np.median(spec) # put the continuum roughly near 1
                wavelist.append(wave)
                speclist.append(spec)
        else: # it's a text file
            for line in f1:
                infile = line.rstrip()
                infile = os.path.join(datapath, infile)
                infilelist.append(infile)
                wave, spec = np.loadtxt(infile, usecols=(0,1), unpack=True)
                wavelist.append(wave)
                speclist.append(spec)
    return infilelist, wavelist, speclist


def simpledespike(wave, spec, delwindow=6, stdfactorup=0.7, stdfactordown=3, plot=True):
    '''
    Implement a simple despiking routine based on the stdev of 1D fluxes
    
    This function is called by despike_spectra, which is where you can adjust the
    values of these arguments as desired:
    
    stdfactorup: Outliers (spikes) are identified as exceeding 
                 stdfactorup*sigma above the continuum
    stdfactordown: Additional outliers (downward spikes) must exceed 
                   stdfactordown*sigma below the continuum
    delwindow: Around each outlier, adjacent points in a window 
               of +/- delwindow are also flagged as outliers
    
    All the outlier points are deleted from input wave, spec to yield newwave, newspec
    '''
    pointstodelete = []
    outliers = (np.where((spec > 1.0 + stdfactorup*np.std(spec)) | 
                         (spec < 1.0 - stdfactordown*np.std(spec))))[0]
    for point in outliers: # add +/- delwindow points around each outlier
        pointstodelete.extend(range(point-delwindow, point+delwindow+1))
    pointstodelete = [point for point in pointstodelete if point >= 0]
    newwave, newspec = np.delete(wave, pointstodelete), np.delete(spec, pointstodelete)
    if plot:
        plt.plot(wave, spec)
        plt.plot(newwave, newspec, color='r')
        plt.xlabel('Wavelength ({\AA})')
        plt.ylabel('Normalized flux')
        plt.axhline(y = (1 + stdfactorup*np.std(spec)), ls=':', color='g')
        plt.axhline(y = (1 - stdfactordown*np.std(spec)), ls=':', color='g')
        plt.show()
    return newwave, newspec


def generalizedESDdespike(wave, spec, maxOLs=1000, alpha=5000):
    '''
    Implement a complicated despiking routine from PyAstronomy
    
    (this function is not tested)
    '''
    r = pyasl.pointDistGESD(spec, maxOLs, alpha)
    # r[0] is number of outliers found, r[i] is indices of outliers
    # maxOLs is max outliers that may be identified; increase alpha to find more
    newwave, newspec = np.delete(wave, r[1]), np.delete(spec, r[1])
    return newwave, newspec


def despike_spectra(wavelist, speclist, type='simple', plot=True):
    '''
    Do one of two things to remove spikes
    
    type='simple' is recommended (see simpledespike)
    type=<anything not 'simple'> will use generalizedESDdespike instead
    
    INPUT
    wavelist: input list of wavelength arrays
    speclist: input list of corresponding flux arrays (for 1D spectra)
    
    OUTPUT
    newwavelist: output file of wavelength arrays without any spikes
    newspeclist: output file of corresponding flux arrays without any spikes
    '''
    newwavelist = []; newspeclist = []
    for wave, spec in zip(wavelist, speclist):
        if type == 'simple':
            delwindow = 6
            stdfactorup = 0.7
            stdfactordown = 3
            newwave, newspec = simpledespike(wave, spec, 
                                             delwindow=delwindow, 
                                             stdfactorup=stdfactorup, 
                                             stdfactordown=stdfactordown,
                                             plot=plot)
        else:
            newwave, newspec = generalizedESDdespike(wave, spec, maxOLs=1000, alpha=5000)
        newwavelist.append(newwave)
        newspeclist.append(newspec)
    return newwavelist, newspeclist


if __name__ == '__main__':
    # parse command line arguments with argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='filelist', required=True,
                        help='text file containing a list of spectra')
    parser.add_argument('-d', dest='datapath', required=True,
                        help='path to where filelist and files in filelist live')
    args = parser.parse_args()
    filelist = args.filelist
    datapath = args.datapath
    if not os.path.exists(os.path.join(datapath, filelist)):
            raise argparse.ArgumentTypeError("{0} does not exist".format(filelist))
    main()
