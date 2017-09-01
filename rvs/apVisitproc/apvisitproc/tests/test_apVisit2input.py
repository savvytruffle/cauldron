import numpy as np
from apvisitproc import apVisit2input as ap
import os

DATAPATH = os.path.dirname(__file__)
KIC = 'testKIC'
FITSFILEPATH = os.path.join(DATAPATH, 'apVisit-r5-7125-56557-285.fits')
# One set of locID, mjd, and fiberID that exist in the test Visitlist
locID = 5215
mjd = 55840
fiberID = 277


def test_load_allvisitinfo():
    '''
    Test loading locID, mjd, and fiberID values from a test Visitlist file
    '''
    locIDs, mjds, fiberIDs = ap.load_allvisitinfo(DATAPATH, KIC)
    assert len(locIDs) > 0
    assert len(locIDs) == len(mjds)
    assert len(mjds) == len(fiberIDs)
    assert float(locID) in locIDs
    assert float(mjd) in mjds
    assert float(fiberID) in fiberIDs


def test_load_apVisit():
    '''
    Test reading data from an apVisit file on disk (uses apogee module tools)
    '''
    result = ap.load_apVisit(DATAPATH, KIC, locID, mjd, fiberID)
    # result contains: fitsfilepath, specfileout, wave, flux, fluxerr
    assert os.path.isfile(result[0])
    wave = result[2]
    flux = result[3]
    fluxerr = result[4]
    assert len(wave) > 0
    assert len(wave) == len(flux)
    assert len(flux) == len(fluxerr)


def test_normalize_spec():
    '''
    Test spectrum normalization for one spectrum
    '''
    result = ap.load_apVisit(DATAPATH, KIC, locID, mjd, fiberID)
    wave = result[2]
    flux = result[3]
    fluxerr = result[4]
    specnorm, specnormerr = ap.normalize_spec(wave, flux, fluxerr, plot=False)
    assert len(specnorm) == len(specnormerr)
    assert len(specnorm) == len(wave)
    # The median of the normalized fluxes should be closer to 1 than before
    assert np.abs(1 - np.median(specnorm)) < np.abs(1 - np.median(flux))
    # The median of the normalized fluxes should be more similar at the beginning
    # and end of the wavelength range than the original fluxes
    assert ((np.median(specnorm[0:50]) - np.median(specnorm[-50:-1])) <
            (np.median(flux[0:50]) - np.median(flux[-50:-1])))


def test_make_BFinfile():
    '''
    Test that HJD and BCV header values can be accessed and contain
    the expected values returned by fitsheader (plus 2400000 for HJD)
    '''
    HJD, BCV = ap.make_BFinfile(FITSFILEPATH)
    assert HJD == 2456557.7326138
    assert BCV == -10.9548025967
