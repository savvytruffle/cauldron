import numpy as np
from apvisitproc import despike
import pytest
import os

DATAPATH = os.path.dirname(__file__)
FILELIST1 = os.path.join(DATAPATH, 'list_of_txt_spectra.txt')
FILELIST2 = os.path.join(DATAPATH, 'list_of_fits_spectra.txt')


@pytest.fixture
def wave_spec_generate():
    '''
    Read in three small chunks of spectra for testing purposes

    'wavelist_speclist_generate' can be used as input to any other test function
    that needs access to the variables it returns!

    wave1, spec1 are a single chunk of 1d spectrum
    wavelist, speclist are lists of three chunks of 1d spectrum
    '''
    wave1, spec1 = np.loadtxt(os.path.join(DATAPATH, 'spec1test.txt'), unpack=True)
    wave2, spec2 = np.loadtxt(os.path.join(DATAPATH, 'spec2test.txt'), unpack=True)
    wave3, spec3 = np.loadtxt(os.path.join(DATAPATH, 'spec3test.txt'), unpack=True)
    wavelist = [wave1, wave2, wave3]
    speclist = [spec1, spec2, spec3]
    return wave1, spec1, wavelist, speclist


@pytest.mark.parametrize('filelist, cond', [
    (FILELIST1, False),
    (FILELIST2, True),
    ])
def test_read_infiles(filelist, cond):
    '''
    Test reading in both text and fits files

    Each resulting wavelength array should be sorted in ascending order
    '''
    infilelist, wavelist, speclist = despike.read_infiles(DATAPATH, filelist, cond)
    assert len(infilelist) > 0
    assert len(infilelist) == len(wavelist)
    assert len(wavelist) == len(speclist)
    for wave in wavelist:
        assert all(value >= 0 for value in wave)
        assert list(np.sort(wave)) == list(wave)
        assert all(np.equal(np.sort(wave), wave))


def test_simpledespike(wave_spec_generate):
    '''
    spike condition is met at pixels 15, 16, 17 and 18
    so indices 9 through 24, inclusive, should be removed
    '''
    wave, spec = wave_spec_generate[0], wave_spec_generate[1]
    newwave, newspec = despike.simpledespike(wave, spec, delwindow=6,
                                             stdfactorup=0.7, stdfactordown=3,
                                             plot=False)
    assert len(newwave) == len(newspec)
    assert len(newwave) <= len(wave)
    assert len(newspec) <= len(spec)
    assert all(np.equal(np.hstack((wave[0:9], wave[25:])), newwave))
    assert all(np.equal(np.hstack((spec[0:9], spec[25:])), newspec))


def test_despike_spectra(wave_spec_generate):
    '''
    Test that new spectra are shorter than the original because the outliers are gone
    '''
    wavelist, speclist = wave_spec_generate[2], wave_spec_generate[3]
    newwavelist, newspeclist = despike.despike_spectra(wavelist, speclist, type='simple', plot=False)
    assert len(newwavelist) == len(wavelist)
    assert len(newspeclist) == len(speclist)
