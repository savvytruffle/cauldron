import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
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

windowcols = 3 #The plot with all the subplots should have three columns
windowrows = 6 #The number of rows the subplots take up, manually set for now.
# CUSTOMIZED BF WIDTH AND PLOT LIMITS
widlimits = [0,15, 0,15]; rvneg = -100; rvpos = 100; ymin = -0.15; ymax = 1.19 # good starting default
xmin = rvneg
xmax = rvpos
fig = plt.figure(1, figsize=(15,10))

# Starter loop, NOT FINISHED YET. Need to save the index of where each visit starts/ends.
visitidx = 0
for idx, (rv, value, error) in enumerate(zip(BFrvaxis, BFvalues, BFerrors)):
    ax = fig.add_subplot(windowrows, windowcols, 1) # removed "i"
    ax.yaxis.set_major_locator(MultipleLocator(0.4))
    #if windowcols == 4 and (i!=1 and i!=5 and i!=9 and i!=13 and i!=17 and i!=21 and i!=25):
    #    ax.set_yticklabels(())
    #if windowcols == 3 #and (i!=1 and i!=4 and i!=7 and i!=10 and i!=13 and i!=16 and i!=19 and i!=22 and i!=25):
    #    ax.set_yticklabels(())
    #if i!=20 and i!=21 and i!=22 and i!=23 and i!=24 and i!=25:
    #if i < visit-windowrows:
    #if i!=13 and i!=14 and i!=15 and i!=16:
    #    ax.set_xticklabels(())
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.axis([xmin, xmax, ymin, ymax])
    plt.tick_params(axis='both', which='major')
    if np.abs(BFrvaxis[idx-1] - rv) > 400:
        visitidx = visitidx + 1
        print(visitidx, idx)
        # save idx at this point; e.g., it should be 400 on the first time through

# Meanwhile, for just one visit, we can hardwire [0:400] as the start:end.
# And we can also adjust the RV zeropoints and BF amplitudes arbitrarily.
plt.plot(BFrvaxis[0:400] - 107, 10*BFvalues[0:400])
plt.plot(BFrvaxis[401:801] - 107, 10*BFvalues[401:801])
plt.plot(BFrvaxis[801:1201] - 107, 10*BFvalues[801:1201])
plt.plot(CCF_rvaxis, CCFvalues[2] - 0.2)
plt.axis([-150, 150, -0.1, 0.5])
plt.xlabel('Arbitrary radial velocity')
plt.ylabel('CCF or BF amplitude')
plt.show()

