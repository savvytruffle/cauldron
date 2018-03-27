import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from astropy.io import fits
from astropy.time import Time
'''
This program makes awesome plots of APOGEE CCFs compared to our BFs.
The idea is for you to run it once for each target.
'''

# Where the apStar file lives
#dir = 'data/5285607/'
#dir = 'data/6864859/' # promising
#dir = 'data/6778289/'
#dir = 'data/6449358/'
#dir = 'data/4285087/'
#dir = 'data/6781535/'
dir = 'data/6131659/' # promising

# The name of the apStar file for a single target
#ccffile = 'apStar-r5-2M19390532+4027346.fits' #5285607
#ccffile = 'apStar-r6-2M19292405+4223363.fits' #6864859
#ccffile = 'apStar-r5-2M19282456+4215080.fits' #6778289
#ccffile = 'apStar-r5-2M19353513+4149543.fits' #6449358
#ccffile = 'apStar-r5-2M19463571+3919069.fits' #4285087
#ccffile = 'apStar-r5-2M19321788+4216489.fits' #6781535
ccffile = 'apStar-r6-2M19370697+4126128.fits' #6131659

# The name of the file with relevant BF info for the same target
#bffile = '5285607BFOutAPstar.txt' # in the same order as APOGEE even though it's not chronologic
#bffile = '6864859BFOutALL.txt'
#bffile = '6778289BFOut.txt'
#bffile = '6449358BFOut.txt'
#bffile = '4285087BFOut.txt'
#bffile = '6781535BFOut.txt'
bffile = '6131659BFOutALL.txt'

# The name of the bjdinfile used with BF_python
# (heliocentric/barycentric velocities in col=(2,); note the top row is for the template)
#BCVfile = '5285607bjdinfile.txt'
#BCVfile = '6864859bjdinfileALL.txt'
BCVfile = '6131659bjdinfileALL.txt'

ccfinfile = dir + ccffile
bfinfile = dir + bffile
bjdinfile = dir + BCVfile


# Read in relevant CCF info from apStar file
hdu = fits.open(ccfinfile)
# The CCF information lives in the 9th HDU, because of course it does
hdu9 = hdu[9].data
CCFvalues = hdu9['CCF'][0][2:]
CCFerrors = hdu9['CCFERR'][0][2:]
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

# Get the timestamp for each CCF visit from the 0th HDU
hdu0 = hdu[0].header
ccftimes = []
for idx in range(1, len(CCFvalues)+1):
    headercard = 'HJD' + str(idx)
    ccftimes.append(hdu0[headercard] + 2400000.)
#print(ccftimes)

ccftimesAstropy = []
for ccftime in ccftimes:
    ccftimesAstropy.append(Time(ccftime, scale='utc', format='jd'))

# Set up the main figure
fig = plt.figure(1, figsize=(15,10))
windowrows = 4
windowcols = 7
fig.text(0.5, 0.04, 'Radial Velocity (km s$^{-1}$)', ha='center', va='center', size='small')
fig.text(0.07, 0.5, 'Arbitrary CCF or BF amplitude', ha='center', va='center', size='small', rotation='vertical')
axlims = [-140, 140, -0.06, 0.56]
ccfyoffset = -0.12
#bfxoffset = -100 # arbitrary systemic RV offset to make CCF and BF line up
#bfxoffset = -61.7 #5285607 visit1
#bfxoffset = -97.4 #6864859 visit1
bfxoffset = -100 #613 estimate
#manualxoffsets = [0, 0, 0, 8, -9, 0, 0,  # 686
#                  0, 0, 0, 0,  0, 0, 0,
#                  0, 0, 0, 0,  0, 0, 0,
#                  0, 0, 0, 0,  0, 0, 0]
manualxoffsets = [14, 0, 12,  13, 12,  0, -2,  # 613
                  0,  0, -7,   0,  0,  0, 15,
                  17, 17, 17, 17, 17, -4, -4,
                  0,  0,   5, 10, 11, 17]
#bfxoffset = 0
#bfyamp = 7 # arbitrary factor to stretch BF values to make CCF and BF line up 686
bfyamp = 9 #613

# Loop over and plot CCF data
for idx, CCFdata in enumerate(CCFvalues):
    ax = fig.add_subplot(windowrows, windowcols, idx+1)
    plt.axis(axlims)
    plt.plot(CCF_rvaxis, CCFvalues[idx] + ccfyoffset)
    plt.text(23, 0.5, ccftimesAstropy[idx].iso[0:10], size=10)
#    plt.text(25, 0.5, '{0:.3f}'.format(ccftimes[idx] - 2450000.), size=10, color='C0')
#    plt.text(25, 0.3, idx, size=10, color='C0')
    plt.subplots_adjust(wspace=0, hspace=0)   

# Read in relevant BF info from the BF infile
bfdata = np.loadtxt(bfinfile, comments='#', usecols=(0,1,2), unpack=True)
BFrvaxis = bfdata[0]
BFvalues = bfdata[1]
BFerrors = bfdata[2]

# Get the timestamp for each BF
with open(bfinfile) as bfinfo:
    bftimes = [float(line[13:]) for line in bfinfo if 'timestamp' in line]
    
# Save the indices that separate each BF's visit in the input file
visitidx = 0
BFindices = []
for idx, (rv, value, error) in enumerate(zip(BFrvaxis, BFvalues, BFerrors)):
    if np.abs(BFrvaxis[idx-1] - rv) > 100:
        visitidx = visitidx + 1
        BFindices.append(idx)

# Read in barycentric velocity correction info from the BJD infile
BCVdata = np.loadtxt(bjdinfile, comments='#', usecols=(2,), unpack=True)
BCVdata = BCVdata[1:]
#print(BCVdata)

# Loop over and plot BF data
for idx in range(0, len(bftimes)):
    ax = fig.add_subplot(windowrows, windowcols, idx+1)
    plt.axis(axlims)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    if (idx!= 18 and idx!=19 and idx!=20 and idx!=21 and idx!=22 and idx!=23 and idx!=24 and idx!=25 and idx!=26):
        ax.set_xticklabels(())
    if (idx!=0 and idx!=7 and idx!=14 and idx!=21):
        ax.set_yticklabels(())
#    plt.text(25, 0.4, '{0:.3f}'.format(bftimes[idx] - 2450000.), size=10, color='C1')
#    plt.text(25, 0.2, idx, size=10, color='C1')
    try:
        plt.plot(BFrvaxis[BFindices[idx]:BFindices[idx+1]] + bfxoffset + BCVdata[idx] + manualxoffsets[idx], 
                 bfyamp*BFvalues[BFindices[idx]:BFindices[idx+1]])
    except: # handle the final case where there is no idx+1
        plt.plot(BFrvaxis[BFindices[idx]::] + bfxoffset + BCVdata[idx] + manualxoffsets[idx], 
                 bfyamp*BFvalues[BFindices[idx]::])

plt.show()
#fig.savefig('6136159BFCCF.png')
