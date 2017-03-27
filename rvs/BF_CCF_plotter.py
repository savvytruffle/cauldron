import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
'''
This program makes awesome plots of APOGEE CCFs compared to our BFs.
The idea is for you to run it once for each target.
'''

# Where the apStar file lives
dir = 'data/6864859/'
# The name of the apStar file for a single target
ccffile = 'apStar-r5-2M19292405+4223363.fits'
ccfinfile = dir + ccffile

# Read in relevant CCF info from apStar file
hdu = fits.open(ccfinfile)
# The CCF information lives in the 9th HDU, because of course it does
hdu9 = hdu[9].data
CCFvalues = hdu9['CCF'][0]
CCFerrors = hdu9['CCFERR'][0]
CCFxaxis = hdu9['CCFLAG'][0]  # pixels, needs to be turned into RVs
CCF_w0 = hdu9['CCFW0'][0]  # starting wavelength in Angstroms ??
CCF_delta = hdu9['CCFDW'][0]  # Delta (log_10 lambda)
# log lambda spacing per lag step (6.d-6 corresponds to 4.145 km/s)

# Get the systemic RV for each visit according to APOGEE
## THIS IS NON-TRIVIAL. APOGEE really likes to remove all RV information.
## There are VREL1 - VRELNVISIT in HDU0, but they're just BC + VREL = VHELIO.
## This is not helpful.
## TODO: figure out actual systemic RVs according to APOGEE
## (or maybe just move our BFs so they're symmetric around 0 and move on with our lives?)
#systemics = hdu9['VREL'][0]
systemics = hdu9['SYNTHVREL'][0]
#bcvs = hdu9['BC'][0]
#print(systemics) # these differ from what we think the systemic RV is
# this would work if the values in systemics were accurate, but they're not
CCF_rvaxis = []
for rv in systemics:
    CCF_rvaxis.append(CCFxaxis * 4.145 + rv)
# for now, let's just use:
CCF_rvaxis = CCFxaxis * 4.145

# Set up the main figure
fig = plt.figure(1, figsize=(12,8))
windowrows = 7
windowcols = 4
fig.text(0.5, 0.04, 'Arbitrary RV (km s$^{-1}$)', ha='center', va='center', size='small')
fig.text(0.07, 0.6, 'CCF or BF amplitude', ha='center', va='center', size='small', rotation='vertical')

# Loop over and plot CCF data
for idx, CCFdata in enumerate(CCFvalues):
    if idx >= 2: # only plot visits, not 0th and 1st columns (combined spectra)
        ax = fig.add_subplot(windowrows, windowcols, idx-1)
        plt.axis([-150, 150, -0.1, 0.5])
        ax.set_xticklabels(()) # add some of these back in later for certain idx
        ax.set_yticklabels(()) # add some of these back in later for certain idx
        plt.plot(CCF_rvaxis, CCFvalues[idx] - 0.2) 
        plt.subplots_adjust(wspace=0, hspace=0)   

# Read in relevant BF info for the same target
bffile = '6864859BFOut.txt' # WARNING: we omitted some visits from this BF run!
bfinfile = dir + bffile

# Read in relevant BF info from the BF infile
bfdata = np.loadtxt(bfinfile, comments='#', usecols=(0,1,2), unpack=True)
BFrvaxis = bfdata[0]
BFvalues = bfdata[1]
BFerrors = bfdata[2]

# Save the indices that separate each BF's visit in the input file
visitidx = 0
BFindices = []
for idx, (rv, value, error) in enumerate(zip(BFrvaxis, BFvalues, BFerrors)):
    if np.abs(BFrvaxis[idx-1] - rv) > 100:
        visitidx = visitidx + 1
        BFindices.append(idx)
print(BFindices)

# Loop over and plot BF data
## You will get several WARNING statements printed if there are fewer BFs than CCFs,
## (and you will also notice that not every plot window has a BF on top of the CCF)
xoffset = 100 # arbitrary systemic RV offset to make CCF and BF line up
yamp = 8 # arbitrary factor to stretch BF values to make CCF and BF line up
for idx in range(0, len(CCFvalues)-2):
    ax = fig.add_subplot(windowrows, windowcols, idx+1)
    plt.axis([-150, 150, -0.1, 0.5])
    ax.set_xticklabels(())
    ax.set_yticklabels(())
    try:
        plt.plot(BFrvaxis[BFindices[idx]:BFindices[idx+1]] - xoffset, 
                 yamp*BFvalues[BFindices[idx]:BFindices[idx+1]])
    except:
        print('WARNING: BF index {0} is out of range, skipping'.format(idx))
        continue

plt.show()

