from __future__ import print_function

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
import apogee.tools.read as apread
from apogee.spec import continuum

'''
Create continuum-normalized text file spectra from a set of APOGEE visit spectra.

This program assumes you are planning to do a broadening function analysis on
the spectra to measure radial velocities, and thus refers to a single "model"
spectrum and a time series of "target" spectra.

Written by Meredith Rawls

Usage
------
$ python apVisit2input.py -d datadir -k KIC
This assumes the apVisit metadata file is datadir/KIC/KICVisitlist.txt, e.g.,
$ python apVisit2input.py -d data -k 1234567

Result
------
A list of the new text file spectra created is printed to the terminal.
Each spectrum's date of observation (HJD) and barycentric velocity (BCV) is
also printed out for easy reference.
'''


def main():
    '''
    Parse arguments, create new spectrum text files, and print useful info.

    More specifically...
    - Read datadir and KIC from command line arguments
    - Load data for all the visits for the target specified
    - Loop over all the visits, normalize the spectra, and save to file
    - Print the names of the new files with their HJDs and BCVs
    '''
    # parse command line arguments with argparse
    parser = argparse.ArgumentParser(description='Run with python apVisit2input.py -d datadir -k KIC')
    parser.add_argument('-d', dest='datadir', required=True,
                        help='directory containing KIC subdirectories')
    parser.add_argument('-k', dest='KIC', required=True,
                        help='KIC of target')
    args = parser.parse_args()
    datadir = args.datadir
    KIC = args.KIC

    locIDs, mjds, fiberIDs = load_allvisitinfo(datadir, KIC)
    print('New spectrum file (despiked), HJD, and BCV:')

    # loop over all visit spectra
    for locID, mjd, fiberID in zip(locIDs, mjds, fiberIDs):
        fitsfilepath, specfileout, wave, flux, fluxerr = load_apVisit(datadir, KIC, locID, mjd, fiberID)
        specnorm, specnormerr = normalize_spec(wave, flux, fluxerr, plot=True)
        with open(specfileout, 'w') as f:
            for wavepoint, specpoint, specnormpoint in zip(wave, specnorm, specnormerr):
                print(wavepoint, specpoint, specnormpoint, file=f)
        HJD, BCV = make_BFinfile(fitsfilepath)
        print(specfileout, HJD, BCV)
    return


def load_allvisitinfo(datadir, KIC):
    '''
    Retrieve necessary metadata from an apVisit metadata file for a target

    NOTE: the user must have an apVisit metadata file already (see README)
    An example file of the correct format is in sample_Visitlist.txt

    Parameters
    ----------
    datadir: `str`
        The directory which must contain 'KIC/KICVisitlist.txt'
    KIC: `str`
        The ID or name of the target. Suggested to be the 7-8 digit Kepler
        identifier, but may be any string that is part of datadir

    Returns
    -------
    locIDs: `list`
        The location IDs for each spectrum, typically 4 digits
    mjds: `list`
        The date of observation in MJD of each spectrum, typically 5 digits
    fiberIDs: `list`
        The ID of the fiber used for each spectrum, typically 3 digits
    '''
    visitlist = os.path.join(datadir, KIC, KIC + 'Visitlist.txt')
    if not os.path.exists(visitlist):
        raise argparse.ArgumentTypeError("{0} does not exist".format(visitlist))
    locIDs, mjds, fiberIDs = np.loadtxt(visitlist, usecols=(1, 2, 3),
                                        unpack=True, delimiter=',')
    return locIDs, mjds, fiberIDs


def load_apVisit(datadir, KIC, locID, mjd, fiberID):
    '''
    Returns original and new filenames plus raw data for one apVisit spectrum

    NOTE: currently only works for data releases DR12 and DR13

    Parameters
    ----------
    datadir: `str`
        The directory which must contain 'KIC/KICVisitlist.txt'
    KIC: `str`
        The ID or name of the target. Suggested to be the 7-8 digit Kepler
        identifier, but may be any string that is part of datadir
    locID: `float`
        The location IDs for one spectrum, typically 4 digits
    mjd: `float`
        The date of observation in MJD of one spectrum, typically 5 digits
    fiberID: `float`
        The ID of the fiber used of one spectrum, typically 3 digits

    Returns
    -------
    fitsfilepath: `str`
        The path to the original apVisit file, as determined by the apogee module
    specfileout: `str`
        The path to the new spectrum text file that will be created
    wave: `list`
        A 1D list of wavelength values for one spectrum
    flux: `list`
        A 1D list of flux values for the same spectrum
    fluxerr: `list`
        A 1D list of flux error values for the same spectrum
    '''
    locID = str(int('{:04d}'.format(int(locID))))
    mjd = str('{:05d}'.format(int(mjd)))
    fiberID = str('{:03d}'.format(int(fiberID)))
    specfileout = os.path.join(datadir, KIC, 'apVisitnorm-'+locID+'-'+mjd+'-'+fiberID+'.txt')
    wave = apread.apVisit(int(locID), mjd, fiberID, ext=4, header=False)
    flux = apread.apVisit(int(locID), mjd, fiberID, ext=1, header=False)
    fluxerr = apread.apVisit(int(locID), mjd, fiberID, ext=2, header=False)
    SDSS_PATH = os.environ.get('SDSS_LOCAL_SAS_MIRROR')
    SDSS_VERSION = os.environ.get('RESULTS_VERS')
    if SDSS_PATH is None:
        raise RuntimeError('You haven\'t defined the environment variable SDSS_LOCAL_SAS_MIRROR')
    if SDSS_VERSION == 'v603':
        drnum = 'dr12'
        rnum = 'r5'
    elif SDSS_VERSION == 'l30e.2':
        drnum = 'dr13'
        rnum = 'r6'
    else:
        raise RuntimeError('You don\'t appear to be using DR12 or DR13, cannot proceed')
    fitsfile = 'apVisit-' + rnum + '-' + locID + '-' + mjd + '-' + fiberID + '.fits'
    fitsfilepath = os.path.join(SDSS_PATH, drnum, 'apogee', 'spectro', 'redux', rnum,
                                'apo25m', locID, mjd, fitsfile)
    return fitsfilepath, specfileout, wave, flux, fluxerr


def normalize_spec(wave, flux, fluxerr, plot=True):
    '''
    Continuum normalize a single spectrum using the apogee module normalizer

    Parameters
    ----------
    wave: `list`
        A 1D list of wavelength values for one spectrum
    flux: `list`
        A 1D list of flux values for the same spectrum
    fluxerr: `list`
        A 1D list of flux error values for the same spectrum
    plot: `bool`
        Choose whether to show an interactive plot of the continuum
        normalization process for visual inspection

    Returns
    -------
    specnorm: `list`
        A 1D list of new normalized flux values for a spectrum (goes with wave)
    specnormerr: `list`
        A 1D list of normalized flux errors for the same spectrum
    '''
    contspec = continuum.fitApvisit(flux, fluxerr, wave)
    specnorm = flux/contspec
    specnormerr = fluxerr/contspec
    if plot:
        plt.plot(wave, flux)
        plt.plot(wave, contspec, lw=2, color='r')
        plt.show()
        plt.plot(wave, specnorm)
        plt.show()
    return specnorm, specnormerr


def make_BFinfile(fitsfilepath):
    '''
    Return info about one spectrum that is useful for the broadening function

    Parameters
    ----------
    fitsfilepath: `str`
        The path to an apVisit file, as determined by the apogee module

    Returns
    -------
    HJD: `float`
        Time of observation from the FITS header, in heliocentric julian days
    BCV: `float`
        Barycentric velocity from the FITS header, in km/s
    '''
    header = fits.open(fitsfilepath)[0].header
    HJD = float('24'+str(header['HJD']))
    BCV = header['BC']
    return HJD, BCV


if __name__ == '__main__':
    main()
