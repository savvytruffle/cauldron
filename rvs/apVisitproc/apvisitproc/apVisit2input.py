from __future__ import print_function
import numpy as np
import apogee.tools.read as apread
from apogee.spec import continuum
import matplotlib.pyplot as plt
import os
from astropy.io import fits
import argparse
'''
This program uses jobovy/apogee to make text files for APOGEE visit spectra and a model spectrum.
The model spectrum should have stellar parameters similar to the target star for use with BF_python.
All the final spectra are continuum normalized.

USAGE:
$ python apVisit2input.py -d datadir -k KIC
This assumes that the apVisit metadata file is datadir/KIC/KICVisitlist.txt
e.g.,
$ python apVisit2input.py -d data -k 1234567
'''


def main():
    '''
    Read datadir and KIC from command line
    Load data for all the visits for the target specified
    Loop over all the visits, normalize the spectra, and save to file
    Print the names of the new files with their HJDs and BCVs
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
    runs on a list of visit spectra
    returns lists of fits files on disk, locIDs, mjds, and fiberIDs for each visit
    '''
    visitlist = os.path.join(datadir, KIC, KIC + 'Visitlist.txt')
    if not os.path.exists(visitlist):
        raise argparse.ArgumentTypeError("{0} does not exist".format(visitlist))
    locIDs, mjds, fiberIDs = np.loadtxt(visitlist, usecols=(1, 2, 3),
                                        unpack=True, delimiter=',')
    return locIDs, mjds, fiberIDs


def load_apVisit(datadir, KIC, locID, mjd, fiberID):
    '''
    runs on a single visit spectrum
    returns the new outfile name and the raw data from the visit spectrum
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
    runs on a single visit spectrum
    returns a continuum-normalized spectrum
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
    runs on a single visit spectrum
    returns metadata for use in making a BF_python infile
    '''
    header = fits.open(fitsfilepath)[0].header
    HJD = float('24'+str(header['HJD']))
    BCV = header['BC']
    return HJD, BCV


if __name__ == '__main__':
    main()
