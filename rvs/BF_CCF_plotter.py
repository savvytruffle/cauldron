import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from astropy.io import fits
from astropy.time import Time
'''
This program makes awesome plots of APOGEE CCFs compared to our BFs.
The idea is for you to run it once for each target.
'''

#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#plt.rc('text', usetex=True)

#KIC = '5285607'; windowrows = 1; windowcols = 6
#KIC = '6864859'; windowrows = 4; windowcols = 7
#KIC = '6778289'; windowrows = 4; windowcols = 7
#KIC = '6449358'; windowrows = 4; windowcols = 7
#KIC = '4285087'; windowrows = 1; windowcols = 6
#KIC = '6781535'; windowrows = 4; windowcols = 7
KIC = '6131659'; windowrows = 4; windowcols = 7

dir = os.path.join('data', KIC)

# The name of the apStar file for a single target
#ccffile = 'apStar-r8-2M19390532+4027346.fits' #5285607
#ccffile = 'apStar-r6-2M19292405+4223363.fits' #6864859  #r8 is missing one visit
#ccffile = 'apStar-r8-2M19282456+4215080.fits' #6778289
#ccffile = 'apStar-r8-2M19353513+4149543.fits' #6449358
#ccffile = 'apStar-r8-2M19463571+3919069.fits' #4285087
#ccffile = 'apStar-r8-2M19321788+4216489.fits' #6781535
ccffile = 'apStar-r8-2M19370697+4126128.fits' #6131659

# The name of the file with relevant BF info for the same target
#bffile = '5285607BFOutAPstar.txt' # in the same order as APOGEE even though it's not chronologic
#bffile = '6864859BFOutALL.txt'
#bffile = '6778289BFOutJCAllVisits.txt'
#bffile = '6449358BFOutALL.txt'
#bffile = '4285087BFOut.txt'
#bffile = '6781535BFOutALL.txt'
bffile = '6131659BFOutALL.txt'

# The name of the bjdinfile used with BF_python
# (heliocentric/barycentric velocities in col=(2,); note the top row is for the template)
#BCVfile = '5285607bjdinfile.txt'
#BCVfile = '6864859bjdinfileALL.txt'
#BCVfile = '6778289bjdinfiles.txt'
#BCVfile = '6449358bjdinfile.txt'
#BCVfile = '4285087bjdinfile.txt'  ## apogee visits aren't chronological WHY
#BCVfile = '6781535bjdinfileALL.txt'
BCVfile = '6131659bjdinfile.txt'

#ccfinfile = os.path.join(dir, ccffile)
ccfinfile = os.path.join('data', ccffile)
bfinfile = os.path.join(dir, bffile)
bjdinfile = os.path.join(dir, BCVfile)

# Read in relevant CCF info from apStar file
hdu = fits.open(ccfinfile)
# The CCF information lives in the 9th HDU, because of course it does
hdu9 = hdu[9].data
CCFvalues = hdu9['CCF'][0][2:]
CCFerrors = hdu9['CCFERR'][0][2:]
CCFxaxis = hdu9['CCFLAG'][0]  # pixels, needs to be turned into RVs
CCF_delta = hdu9['CCFDW'][0]  # Delta (log_10 lambda)
#print(hdu9['VHELIO'][0] - hdu9['VREL'][0])
# (the data model website says this log lambda spacing per lag step of
# 6.d-6 corresponds to 4.145 km/s)

# Get the systemic RV of the star according to APOGEE
apSystemic = 82  ## TODO: FIND ACTUAL SYSTEMIC RV PARAMETER !!!
# Get the barycentric velocities of each visit according to APOGEE
# This is called VHELIO for super obvious reasons, Jen Sobeck private communication, for reals
apBCVs = hdu9['VHELIO'][0]
CCF_rvaxis = [CCFxaxis * 4.145 + bcv for bcv in apBCVs]

# Get the timestamp for each CCF visit from the 0th HDU
hdu0 = hdu[0].header
ccftimes = []
for idx in range(1, len(CCFvalues)+1):
    headercard = 'HJD' + str(idx)
    ccftimes.append(hdu0[headercard] + 2400000.)

ccftimesAstropy = []
for ccftime in ccftimes:
    ccftimesAstropy.append(Time(ccftime, scale='utc', format='jd'))

# Set up the main figure
fig = plt.figure(1, figsize=(15,10))
ax = fig.add_subplot(111)
ax.tick_params(top='off', bottom='off', left='off', right='off')
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.set_xlabel('Radial Velocity (km s$^{-1}$)', labelpad=20, size='x-large')
ax.set_ylabel('Arbitrary CCF or BF amplitude', labelpad=20, size='x-large')
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
plt.title('KIC ' + KIC, size='x-large')

#windowrows = 4
#windowcols = 7
#fig.text(0.5, 0.04, 'Radial Velocity (km s$^{-1}$)', ha='center', va='center', size='x-large')
#fig.text(0.07, 0.5, 'Arbitrary CCF or BF amplitude', ha='center', va='center', size='x-large', rotation='vertical')
axlims = [-140, 140, -0.06, 0.56]
ccfyoffset = -0.12
bfyamp = 9 # arbitrary factor to stretch BF values for clarity
# used 7 for 686 and 9 for 613

# Loop over and plot CCF data
for idx, CCFdata in enumerate(CCFvalues):
    ax0 = fig.add_subplot(windowrows, windowcols, idx+1)
    plt.axis(axlims)
    plt.plot(CCF_rvaxis[idx], CCFvalues[idx] + ccfyoffset)
    plt.text(23, 0.5, ccftimesAstropy[idx].iso[0:10], size=10)
    ## TURN ON FOR TIMESTAMP TROUBLESHOOTING
    #plt.text(25, 0.5, '{0:.3f}'.format(ccftimes[idx] - 2450000.), size=10, color='C0')
    #plt.text(25, 0.3, idx, size=10, color='C0')
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

# Loop over and plot BF data
for idx in range(0, len(bftimes)):
    ax1 = fig.add_subplot(windowrows, windowcols, idx+1)
    plt.axis(axlims)
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    if (idx < 20):  # set to 18 for 686 and 20 for 613
        ax1.set_xticklabels(())
    if (idx!=0 and idx!=7 and idx!=14 and idx!=21):
        ax1.set_yticklabels(())
    ## TURN ON FOR TIMESTAMP TROUBLESHOOTING
    #plt.text(25, 0.4, '{0:.3f}'.format(bftimes[idx] - 2450000.), size=10, color='C1')
    #plt.text(25, 0.2, idx, size=10, color='C1')
    try:
        plt.plot(BFrvaxis[BFindices[idx]:BFindices[idx+1]] - apSystemic + BCVdata[idx], 
                 bfyamp*BFvalues[BFindices[idx]:BFindices[idx+1]])
    except: # handle the final case where there is no idx+1
        try:
            plt.plot(BFrvaxis[BFindices[idx]::] - apSystemic + BCVdata[idx], 
                    bfyamp*BFvalues[BFindices[idx]::])
        except:
            print('You\'re missing a BF where an APOGEE CCF exists')
            continue

plt.show()
