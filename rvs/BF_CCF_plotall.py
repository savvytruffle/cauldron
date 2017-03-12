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

# Get the systemic RV for each visit
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

# Read in relevant BF info for the same target
bffile = '6864859BFOut.txt'
bfinfile = dir + bffile

# Read in relevant BF info from the BF infile
bfdata = np.loadtxt(bfinfile, comments='#', usecols=(0,1,2), unpack=True)
BFrvaxis = bfdata[0]
BFvalues = bfdata[1]
BFerrors = bfdata[2]

# Starter loop, NOT FINISHED YET. Need to save the index of where each visit starts/ends.
visitidx = 0
for idx, (rv, value, error) in enumerate(zip(BFrvaxis, BFvalues, BFerrors)):
    if np.abs(BFrvaxis[idx-1] - rv) > 100:
   # if np.abs(BFrvaxis[idx-1] - rv) > 400:
        visitidx = visitidx + 1
        print(visitidx, idx)
        # save idx at this point; e.g., it should be 400 on the first time through
        #408 is where "###" happens, but 414 is where the next pocket of data starts

# Meanwhile, for just one visit, we can hardwire [0:400] as the start:end.
# And we can also adjust the RV zeropoints and BF amplitudes arbitrarily.
plt.plot(BFrvaxis[0:400] - 107, 10*BFvalues[0:400])
plt.plot(CCF_rvaxis, CCFvalues[2] - 0.2)
plt.axis([-150, 150, -0.1, 0.5])
plt.xlabel('Arbitrary radial velocity')
plt.ylabel('CCF or BF amplitude')
plt.show()

