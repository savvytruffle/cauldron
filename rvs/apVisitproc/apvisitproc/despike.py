from __future__ import print_function

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from PyAstronomy import pyasl

'''
Remove spikes (tellurics) from continuum-normalized apVisit APOGEE spectra.

Typically you will run this after apVisit2input.py, which finds the apVisit
spectra downloaded via the python apogee module for one target and
continuum-normalizes them. If you somehow have an apVisit-style FITS file
that is continuum-normalized already (or close), this can despike it too.

First, read_infiles reads in wavelengths and fluxes for a list of spectra;
Then, despike_spectra despikes the spectra with one of two techniques.

Usage
-----
python despike.py -d datapath -i filelist

datapath: /path/to/filelist/and/files/listed/in/filelist

filelist: Single-column text file with list of *continuum-normalized*
          single-visit APOGEE files you want to despike.
          The files in this list can either be two-column text files
          (wave, flux) or FITS files in the apVisit format.

Result
------
The new despiked spectra are written to two-column (wavelength, flux) files
with similar names as the original, but they now end in '_despiked.txt'.
'''


def main():
    '''
    Parse arguments, despike a set of spectra, and write the results to files.
    '''
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

    infilelist, wavelist, speclist = read_infiles(datapath, filelist)
    newwavelist, newspeclist = despike_spectra(wavelist, speclist)

    # write out a set of two-column text files,
    # each containing one element of newwavelist and one element of newspeclist
    for file, newwave, newspec in zip(infilelist, newwavelist, newspeclist):
        # create outfile based on infile name
        outfile = os.path.splitext(file)[0] + '_despiked.txt'
        with open(outfile, 'w') as f:
            for wentry, sentry in zip(newwave, newspec):
                print(wentry, sentry, file=f)
    return


def read_infiles(datapath, filelist, isFits=False):
    '''
    Load a text file containing a list of continuum-normalized spectra

    Parameters
    ----------
    datapath: `str`
        Path to the directory containing both filelist and the files therein
    filelist: `str`
        Name of a text file containing a list of continuum-normalized apVisit
        spectra you want to despike. The files in this list can either be
        two-column text files (wave, flux) or FITS files in the apVisit format.
    isFits: `bool`
        True if the files listed in filelist are FITS files, else False.

    Returns
    -------
    infilelist: `list`
        The full path to each spectrum file in filelist
    wavelist: `list`
        A list of lists containing wavelength values for each spectrum
    speclist: `list`
        A list of lists containing flux values for the same spectra
    '''
    print(datapath, filelist)
    wavelist = []
    speclist = []
    with open(os.path.join(datapath, filelist)) as f1:
        infilelist = []  # for use later to make outfiles
        if isFits:  # it's a FITS file
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
                spec = spec / np.median(spec)  # put the continuum roughly near 1
                wavelist.append(wave)
                speclist.append(spec)
        else:  # it's a text file
            for line in f1:
                infile = line.rstrip()
                infile = os.path.join(datapath, infile)
                infilelist.append(infile)
                wave, spec = np.loadtxt(infile, usecols=(0, 1), unpack=True)
                wavelist.append(wave)
                speclist.append(spec)
    return infilelist, wavelist, speclist


def simpledespike(wave, spec, delwindow=6, stdfactorup=0.7, stdfactordown=3, plot=True):
    '''
    Implement a simple despiking routine based on the stdev of 1D fluxes

    All outlier points are deleted from wave, spec to yield newwave, newspec.

    Parameters
    ----------
    wave: `list`
        A 1D list of wavelength values for one spectrum
    spec: `list`
        A 1D list of flux values for the same spectrum
    delwindow: `int`
        Around each outlier (upward spike), adjacent points in a window of +/-
        delwindow are also flagged as outliers
    stdfactorup: `float`
        Outliers (upward spikes) are identified as exceeding stdfactorup*sigma
        above the continuum
    stdfactordown: `float`
        Additional outliers (downward spikes) are identified as exceeding
        stdfactordown*sigma below the continuum
    plot: `bool`
        True = show an interactive plot of each spectrum as it is despiked

    Returns
    -------
    newwave: `list`
        A 1D list of wavelength values for one despiked spectrum
    newspec: `list`
        A 1D list of flux values for the same despiked spectrum
    '''
    pointstodelete = []
    outliers = (np.where((spec > 1.0 + stdfactorup*np.std(spec)) |
                         (spec < 1.0 - stdfactordown*np.std(spec))))[0]
    for point in outliers:  # add +/- delwindow points around each outlier
        pointstodelete.extend(range(point-delwindow, point+delwindow+1))
    pointstodelete = [point for point in pointstodelete if point >= 0]
    newwave, newspec = np.delete(wave, pointstodelete), np.delete(spec, pointstodelete)
    if plot:
        plt.plot(wave, spec)
        plt.plot(newwave, newspec, color='r')
        plt.xlabel('Wavelength ({\AA})')
        plt.ylabel('Normalized flux')
        plt.axhline(y=(1 + stdfactorup*np.std(spec)), ls=':', color='g')
        plt.axhline(y=(1 - stdfactordown*np.std(spec)), ls=':', color='g')
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
    Use one of two techniques to remove spikes from a series of spectra.

    Parameters
    ----------
    wavelist: `list`
        Input list of wavelength arrays
    speclist: `list`
        Input list of corresponding flux arrays (for 1D spectra)
    type: `str`
        type='simple' is recommended (see simpledespike)
        type=<anything not 'simple'> will use generalizedESDdespike instead
    plot: `bool`
        True = show an interactive plot of each spectrum as it is despiked

    Returns
    -------
    newwavelist: `list`
        A list of lists containing wavelength values for each despiked spectrum
    newspeclist: `list`
        A list of lists containing flux values the same despiked spectra
    '''
    newwavelist = []
    newspeclist = []
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
    main()
