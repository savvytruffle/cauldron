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

######                    This is where the plotting happens                     ######

###Reading this in so that the plot can have the appropriate number of subplots###
infiles = 'data/4285087/4285087infiles.txt'
specdata = bff.read_specfiles(infiles)
nspec = specdata[0]; filenamelist = specdata[1]
windowcols = 3 # 4                             # how many columns the plot should have
#windowrows = 6                                # manually set number of plot rows here, or automatically below
windowrows = int([np.rint((nspec-1)/windowcols) if (np.float(nspec-1)/windowcols)%windowcols == 0 else np.rint((nspec-1)/windowcols)+1][0])
xmin = rvneg
xmax = rvpos
fig = plt.figure(1, figsize=(15,10))
fig.text(0.5, 0.04, 'Arbitrary Radial Velocity (km s$^{-1}$)', ha='center', va='center', size='large')
#########0.5, 0.04
fig.text(0.07, 0.6, 'Broadening Function/Cross Correlation Function Amplitude', ha='center', va='center', size='large', rotation='vertical')
#########0.07, 0.5

# Starter loop, NOT FINISHED YET. Need to save the index of where each visit starts/ends.
visitidx = 0
for idx, (rv, value, error) in enumerate(zip(BFrvaxis, BFvalues, BFerrors)):
    if np.abs(BFrvaxis[idx-1] - rv) > 100:
   # if np.abs(BFrvaxis[idx-1] - rv) > 400:
        visitidx = visitidx + 1
        print(visitidx, idx)
        # save idx at this point; e.g., it should be 400 on the first time through
        # 408 is where "###" happens, but 414 is where the next pocket of data starts
        
    plt.plot(BFrvaxis[0:400] - 107, 10*BFvalues[0:400], color=colors[14], lw=1.5,ls='-', label='BF')

for i in range (1,nspec):
    ax = fig.add_subplot(windowrows, windowcols, i) # out of range if windowcols x windowrows < nspec
    ax.yaxis.set_major_locator(MultipleLocator(0.4))
    if windowcols == 4 and (i!=1 and i!=5 and i!=9 and i!=13 and i!=17 and i!=21 and i!=25):
        ax.set_yticklabels(())
    if windowcols == 3 and (i!=1 and i!=4 and i!=7 and i!=10 and i!=13 and i!=16 and i!=19 and i!=22 and i!=25):
        ax.set_yticklabels(())
    #if i!=20 and i!=21 and i!=22 and i!=23 and i!=24 and i!=25:
    if i < nspec-windowrows:
    #if i!=13 and i!=14 and i!=15 and i!=16:
        ax.set_xticklabels(())
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.axis([xmin, xmax, ymin, ymax])
    plt.tick_params(axis='both', which='major')
    #plt.text(xmax - 0.16*(np.abs(xmax-xmin)), 0.75*ymax, '%.3f $\phi$' % (phase[i]), size='small')
    #plt.text(xmax - 0.26*(np.abs(xmax-xmin)), 0.55*ymax, '%s' % (datetimelist[i].iso[0:10]), size='small')
    #plt.plot(bf_ind, bfsmoothlist[i], color=colors[14], lw=1.5, ls='-', label='Smoothed BF')
    
###########        This is where we tell it to plot the BF and CCF, but can they both be in the same for loop?      #####
    
    plt.plot(CCF_rvaxis, CCFvalues[2] - 0.2, color=colors[0], lw=3, ls='-', label='CCF')

    # OPTION TO PLOT VERTICAL LINE AT ZERO
    #plt.axvline(x=0, color=colors[15])    
    # print legend
    if i==nspec-1: ax.legend(bbox_to_anchor=(2.6,0.7), loc=1, borderaxespad=0., 
                        frameon=False, handlelength=3, prop={'size':20})
plt.show()
#######

