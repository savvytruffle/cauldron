import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
'''
This program makes awesome plots of APOGEE CCFs compared to our BFs.
The idea is for you to run it once for each target.
'''

# Where the apStar file lives
#dir = 'data/6864859/'
dir = 'data/5285607/'
#dir = 'data/5284133/'
#dir = 'data/4285087/'
#dir = 'data/6449358/'
#dir = 'data/6781535/'
#dir = 'data/6778289/'
#dir = 'data/4851217/'
#dir = 'data/6610219/'
# The name of the apStar file for a single target
#ccffile = 'apStar-r5-2M19245882+4714573.fits' #6864859
ccffile = 'apStar-r5-2M19390532+4027346.fits' #5285607
#ccffile = 'apStar-r5-2M19373173+4027078.fits' #5284133
#ccffile = 'apStar-r5-2M19321788+4216489.fits' #6781535
#ccffile = 'apStar-r5-2M19353513+4149543.fits' #6449358
#ccffile = 'apStar-r5-2M19463571+3919069.fits' #4285087
#ccffile = 'apStar-r5-2M19282456+4215080.fits' #6778289
#ccffile = 'apStar-r5-2M19320615+4200049.fits' #6610219
#ccffile = 'apStar-r5-2M19432016+3957081.fits' #4851217

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

# Get the timestamp for each CCF visit from the 0th HDU
hdu0 = hdu[0].header
ccftimes = []
for idx in range(1, len(CCFvalues)-1):
    headercard = 'HJD' + str(idx)
    ccftimes.append(hdu0[headercard] + 2400000.)

# Set up the main figure
fig = plt.figure(1, figsize=(12,8))
windowrows = 2
#windowrows = 7
windowcols = 4
fig.text(0.5, 0.04, 'Arbitrary RV (km s$^{-1}$)', ha='center', va='center', size='small')
fig.text(0.07, 0.6, 'CCF or BF amplitude', ha='center', va='center', size='small', rotation='vertical')

# Loop over and plot CCF data
for idx, CCFdata in enumerate(CCFvalues):
    ax = fig.add_subplot(windowrows, windowcols, idx+1)
    plt.axis([40, 150, -0.1, 0.5])
    ax.set_xticklabels(()) # add some of these back in later for certain idx
    ax.set_yticklabels(()) # add some of these back in later for certain idx
    plt.plot(CCF_rvaxis, CCFvalues[idx] - 0.2)
    if idx < 2:
        plt.text(20, 0.4, 'combined', size='10', color='C0')
    if idx >= 2:
        plt.text(20, 0.4, '{0:.5f}'.format(ccftimes[idx-2]), size='10', color='C0')
    plt.subplots_adjust(wspace=0, hspace=0)   

# Read in relevant BF info for the same target
#bffile = '6864859BFOut.txt' # 6864859 WARNING: we omitted some visits from this BF run!
#bffile = '6864859BFOutALL.txt' #6864859 ALL Visits!
bffile = '5285607BFOutAPstar.txt' 
#bffile = '5284133BFOut.txt' 
#bffile = '6449358BFOut.txt'
#bffile = '6781535BFOut.txt'
#bffile = '4851217BFOut.txt'
#bffile = '6778289BFOut.txt'
#bffile = '6610219BFOut.txt'
#bffile = '4285087BFOut.txt'
bfinfile = dir + bffile

#Read in heliocentric/barycentric velocities for the same target
HCV = '5285607bjdinfileApstar.txt'
#HCV = '6864859bjdinfile.txt

bjdinfile = dir + HCV

# Read in relevant BF info from the BF infile
bfdata = np.loadtxt(bfinfile, comments='#', usecols=(0,1,2), unpack=True)
BFrvaxis = bfdata[0]
BFvalues = bfdata[1]
BFerrors = bfdata[2]

# Read in relevant BF info from the BF infile and add barycentric velocity correction
HCVdata = np.loadtxt(bjdinfile, comments='#', usecols=(2,), unpack=True)
HCVx = HCVdata[3]

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

# Loop over and plot BF data
#xoffset = 100 # arbitrary systemic RV offset to make CCF and BF line up
#xoffset = 61.7 #5285607
#xoffset = 97.4 #6864859
yamp = 8 # arbitrary factor to stretch BF values to make CCF and BF line up
for idx in range(0, len(bftimes)):
    ax = fig.add_subplot(windowrows, windowcols, idx+3)
    plt.axis([-150, 150, -0.1, 0.5])
    ax.set_xticklabels(())
    ax.set_yticklabels(())
    plt.text(20, 0.28, bftimes[idx], size='10', color='C1')
    try:
        #plt.plot(BFrvaxis[BFindices[idx]:BFindices[idx+1]] - xoffset, 
        #         yamp*BFvalues[BFindices[idx]:BFindices[idx+1]])
        plt.plot(BFrvaxis[BFindices[idx]:BFindices[idx+1]] - HCVx[idx+1], 
                 yamp*BFvalues[BFindices[idx]:BFindices[idx+1]])
    except: # handle the final case where there is no idx+1
        plt.plot(BFrvaxis[BFindices[idx]::] - HCVx[idx+1], 
                 yamp*BFvalues[BFindices[idx]::])

plt.show()
#fig.savefig('5286507BFCCF.png')
