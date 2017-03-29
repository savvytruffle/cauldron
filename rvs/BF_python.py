from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astropy.io import fits
from astropy.time import Time
from PyAstronomy import pyasl
from scipy import ndimage
import pandas as pd
import gaussfitter as gf
import BF_functions as bff

'''
Program to extract radial velocities from a double-lined binary star spectrum.
Uses the Broadening Function technique.

Meredith Rawls
2014-2015

Based loosely on Rucinski's BFall_IDL.pro, and uses the PyAstronomy tools.
http://www.astro.utoronto.ca/~rucinski/BFdescription.html
http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/svd.html

In practice, you will run this twice: once to do the initial BF, and then again
to properly fit the peaks of each BF with a Gaussian.

INPUT
infiles:    single-column file with one FITS or TXT filename (w/ full path) per line
            1st entry must be for the template star (e.g., arcturus or phoenix model)
            (the same template is used to find RVs for both stars)
            NO comments are allowed in this file
            FUN FACT: unless APOGEE, these should be continuum-normalized to 1 !!!
bjdinfile:     columns 0,1,2 must be filename, BJD, BCV (e.g., from IRAF bcvcorr)
            top row must be for the template star (e.g., arcturus)
            (the 0th column is never used, but typically looks like infiles_BF.txt)
            one line per observation
            comments are allowed in this file using #
gausspars:    your best initial guesses for fitting gaussians to the BF peaks
            the parameters are [amp1, offset1, width1, amp2, offset2, width2]
            the top line is ignored (template), but must have six values
            one line per observation
            comments are allowed in this file using #

OUTPUT
outfile:    a file that will be created with 8 columns: BJD midpoint, orbital phase,
            Kepler BJD, RV1, RV1 error, RV2, RV2 error
bfoutfile:  a file that contains all the BF function data (raw RV, BF, gaussian model)

IMMEDIATELY BELOW, IN THE CODE
You need to specify whether you have APOGEE (near-IR) or "regular" (e.g., ARCES)
spectra with the 'isAPOGEE' flag. You also need to set the binary's PERIOD and BJD0,
both in days, and the constant RV and BCV of whatever template you are using.
'''

##########
# YOU NEED TO HAVE THESE INPUT FILES !!!
# THE OUTPUT FILE WILL BE CREATED FOR YOU

# EXAMPLE INFILES AND OUTFILES
#infiles =   'infiles.txt'; bjdinfile = 'bjdinfile.txt'
#gausspars = 'gausspars.txt'
#outfile =   'rvoutfile.txt'; bfoutfile = 'bfoutfile.txt'

#4851217
#infiles =   'data/4851217/4851217infiles.txt'; bjdinfile = 'data/4851217/4851217bjdinfile.txt'
#gausspars = 'data/4851217/4851217gausspars.txt'
#outfile =   'data/4851217/4851217Outfile.txt'

#5285607
#infiles =   'data/5285607/5285607infiles.txt'; bjdinfile = 'data/5285607/5285607bjdinfile.txt'
#gausspars = 'data/5285607/5285607gausspars.txt'
#outfile =   'data/5285607/5285607Outfile-Meredith.txt'; bfoutfile = 'data/5285607/5285607BFdata.txt'

#4075064
#infiles = 'data/4075064/4075064infiles.txt'; bjdinfile = 'data/4075064/4075064bjdinfile.txt'
#gausspars = 'data/4075064/4075064gausspars.txt'
#outfile = 'data/4075064/4075064outfile.txt'; bfoutfile = 'data/4075064/4075064BFdata.txt'

#3848919
#infiles = 'data/3848919/3848919infiles.txt'; bjdinfile = 'data/3848919/3848919bjdinfile.txt'
#gausspars = 'data/3848919/3848919gausspars.txt'
#outfile = 'data/3848919/3848919outfile.txt'; bfoutfile = 'data/3848919/3848919BFdata.txt'

#6610219
#infiles = 'data/6610219/6610219infiles.txt'; bjdinfile = 'data/6610219/6610219bjdinfile.txt'
#gausspars = 'data/6610219/6610219gausspars1.txt'
#outfile = 'data/6610219/6610219outfile.txt'; bfoutfile = 'data/6610219/6610219BFOut.txt'

#4285087
#infiles = 'data/4285087/4285087infiles.txt'; bjdinfile = 'data/4285087/4285087bjdinfile.txt'
#gausspars = 'data/4285087/4285087gausspars.txt'
#outfile = 'data/4285087/4285087outfile.txt'; bfoutfile = 'data/4285087/4285087BFdata.txt'

#6449358
infiles =   'data/6449358/6449358infiles1.txt'; bjdinfile = 'data/6449358/6449358bjdinfile1.txt'
gausspars = 'data/6449358/6449358gausspars.txt'
outfile =   'data/6449358/6449358Outfile1.txt'; bfoutfile = 'data/6449358/6449358BFOut.txt'

#5284133
#infiles =   'data/5284133/5284133infiles.txt'; bjdinfile = 'data/5284133/5284133bjdinfile.txt'
#gausspars = 'data/5284133/5284133gausspars.txt'
#outfile =   'data/5284133/5284133Outfile.txt'; bfoutfile = 'data/5284133/5284133BFOut.txt'

#6778289
#infiles =   'data/6778289/6778289infiles.txt'; bjdinfile = 'data/6778289/6778289bjdinfile.txt'
#gausspars = 'data/6778289/6778289gausspars.txt'
#outfile =   'data/6778289/6778289Outfile.txt'; bfoutfile = 'data/6778289/6778289BFOut.txt'

#6781535
#infiles =   'data/6781535/6781535infiles.txt'; bjdinfile = 'data/6781535/6781535bjdinfile.txt'
#gausspars = 'data/6781535/6781535gausspars.txt'
#outfile =   'data/6781535/6781535Outfile.txt'; bfoutfile = 'data/6781535/6781535BFOut.txt'

#6864859
#infiles =   'data/6864859/6864859infiles1.txt'; bjdinfile = 'data/6864859/6864859bjdinfile1.txt'
#gausspars = 'data/6864859/6864859gausspars1.txt'
#outfile =   'data/6864859/6864859Outfile-Meredith.txt'; bfoutfile = 'data/6864859/6864859BFOut.txt'

# ORBITAL PERIOD AND ZEROPOINT !!!
#period = 2.47028; BJD0 = 2455813.69734 # 4851217
#period = 3.8994011; BJD0 = 2454959.576010 # 5285607
period = 5.7767904; BJD0 = 2456760.90580 # 6449358
#period = 8.7845759; BJD0 = 245800.46231 #5284133
#period = 30.13015; BJD0 = 2456557.73097 #6778289
#period = 9.1220856; BJD0 = 2456557.733 #6781535
#period = 40.8778427; BJD0 = 2454955.556300 #6864859
#period = 61.4228063; BJD0 = 2455813.69734 #4075064
#period = 1.0472603; BJD0 = 2455811.61005 #3848919
#period = 11.3009948; BJD0 = 2456557.73097 #6610219
#period = 4.4860312; BJD0 = 2455813.69734

# STUFF YOU NEED TO DEFINE CORRECTLY !!!
isAPOGEE = True        # toggle to use near-IR stuff, or not
SpecPlot = True         # toggle to plot spectra before BFs, or not
bjdoffset = 2454833.    # difference between real BJDs and 'bjdfunny' (truncated BJDs)
amplimits = [0,1.2, 0,1.2] # limits for gaussian normalized amplitude [min1,max1,min2,max2]
threshold = 10             # margin for gaussian position (raw RV in km/s)
#widlimits = [0,25, 0,22]   # limits for gaussian width (km/s) [min1,max1,min2,max2]
# ^^^ widlimits IS NOW SPECIFIED ON A PER-STAR BASIS BELOW

# RADIAL VELOCITY AND BCV INFO FOR TEMPLATE (km/s; set both to 0 if using a model !!!)
rvstd = 0; bcvstd = 0 # model template

# PARAMETERS FOR THE BROADENING FUNCTION (IMPORTANT PAY ATTENTION !!!)
smoothstd = 1.5      # stdev of Gaussian to smooth BFs by (~slit width in pixels)
#w00 = 5400          # starting wavelength for new grid
#n = 38750           # number of wavelength points for new grid
#stepV = 1.7         # roughly 3e5 / (max_wavelength / wavelength_step) km/s, rounded down
m = 401              # length of the BF (must be longer if RVs are far from 0)
## good values for APOGEE:
#w00 = 15170; n = 32000; stepV = 1.0 # all of APOGEE, (too) high res
w00 = 15170; n = 22000; stepV = 1.5 # all of APOGEE, still pretty high res
#w00 = 15170; n = 2000; stepV = 4.0 # a little piece of APOGEE (lower res, apStar)

# CUSTOMIZED BF WIDTH AND PLOT LIMITS
widlimits = [0,15, 0,15]; rvneg = -100; rvpos = 100; ymin = -0.15; ymax = 1.19 # good starting default
#widlimits = [0,15, 0,15]; rvneg = -70; rvpos = 270; ymin = -0.15; ymax = 1.19 # 5285607
#widlimits = [0,5, 0,5]; rvneg = 0; rvpos = 200; ymin = -0.15; ymax = 1.1 #6449358
#widlimits = [0,5, 0,5]; rvneg = 0; rvpos = 200; ymin = -0.15; ymax = 1.1 #6778289
#widlimits = [0,9, 0,9]; rvneg = 30; rvpos = 170; ymin = -0.15; ymax = 1.19 # 6864859
#widlimits = [0,9, 0,9]; rvneg = -150; rvpos = 50; ymin = -0.15; ymax = 1.19 # 6610259a
#widlimits = [0,15, 0,15]; rvneg = -50; rvpos = 10; ymin = -0.15; ymax = 1.19 # 6610219b


colors = bff.user_rc()

print('Welcome to the Broadening Function party!')
print('')
print('MAKE SURE THIS IS WHAT YOU WANT:')
print('You set Porb = {0} days, BJD0 = {1} days'.format(period, BJD0))

# CREATE NEW SPECTRUM IN LOG SPACE
# This uses w00, n, and stepV, defined above. The new wavelength grid is w1.
# The BF will be evenly spaced in velocity with length m.
# The velocity steps are r (km/s/pix).
w1, m, r = bff.logify_spec(isAPOGEE, w00, n, stepV, m)

# READ IN ALL THE THINGS
specdata = bff.read_specfiles(infiles, bjdinfile, isAPOGEE)
nspec = specdata[0]; filenamelist = specdata[1]
datetimelist = specdata[2]; wavelist = specdata[3]; speclist = specdata[4]

# INTERPOLATE THE TEMPLATE AND OBJECT SPECTRA ONTO THE NEW LOG-WAVELENGTH GRID
# OPTION TO PLOT THIS
newspeclist = []
yoffset = 0
if SpecPlot == True:
    plt.axis([w1[0], w1[-1], 0, nspec+3])
    plt.xlabel(r'Wavelength ({\AA})')
for i in range (0, nspec):
    newspec = np.interp(w1, wavelist[i], speclist[i])
    newspeclist.append(newspec)
    if SpecPlot == True:
        if i == 0: # plot template in red
            plt.plot(w1, newspec+yoffset, label=datetimelist[i].iso[0:10], color=colors[6], marker='.')
        else: # plot the rest in blue
            plt.plot(w1, newspec+yoffset, label=datetimelist[i].iso[0:10], color=colors[0], marker='.')
    yoffset = yoffset + 1
if SpecPlot == True:
    ##plt.legend()
    plt.show()

# BROADENING FUNCTION TIME
svd = pyasl.SVD()
# Single Value Decomposition
svd.decompose(newspeclist[0], m)
singularvals = svd.getSingularValues()
bflist = []
bfsmoothlist = []
for i in range (0, nspec):
    # Obtain the broadening function
    bf = svd.getBroadeningFunction(newspeclist[i]) # this is a full matrix
    bfarray = svd.getBroadeningFunction(newspeclist[i], asarray=True)
    # Smooth the array-like broadening function
    # 1ST LINE - python 2.7 with old version of pandas; 2ND LINE - python 3.5 with new version of pandas
    #bfsmooth = pd.rolling_window(bfarray, window=5, win_type='gaussian', std=smoothstd, center=True)
    bfsmooth = pd.Series(bfarray).rolling(window=5, win_type='gaussian', center=True).mean(std=smoothstd)
    # The rolling window makes nans at the start because it's a punk.
    for j in range(0,len(bfsmooth)):
        if np.isnan(bfsmooth[j]) == True:
            bfsmooth[j] = 0
        else:
            bfsmooth[j] = bfsmooth[j]
    bflist.append(bf)
    bfsmoothlist.append(bfsmooth)
    
bfnormlist = []
for a in bfsmoothlist:
    bfnormlist.append((a-np.min(a))/(np.max(a)-np.min(a)))

# Obtain the indices in RV space that correspond to the BF
bf_ind = svd.getRVAxis(r, 1) + rvstd - bcvstd

# OPTION TO PLOT THE SINGULAR VALUES TO SEE WHERE THEY AREN'T A MESS
# this probably isn't important, because instead of choosing which values to throw out,
# we use "Route #2" as described by Rucinski and just use the final row of the BF array
# and smooth it with a Gaussian to get rid of noise problems.
# for more info, seriously, read http://www.astro.utoronto.ca/~rucinski/SVDcookbook.html
##plt.figure(2)
#plt.semilogy(singularvals, 'b-')
#plt.xlabel('BF Index')
#plt.ylabel('Singular Values')
#plt.show()

# OPTION TO PLOT THE SMOOTHED BFs
plt.axis([rvneg, rvpos, -0.2, float(nspec)+1])
plt.xlabel('Radial Velocity (km s$^{-1}$)')
plt.ylabel('Broadening Function (arbitrary amplitude)')
yoffset = 0.0
for i in range(1, nspec):
    plt.plot(bf_ind, bfnormlist[i]+yoffset, color=colors[0], marker='.')
    plt.axhline(y=yoffset, color=colors[15], ls=':')
    yoffset = yoffset + 1.0
plt.show()

# FIT THE SMOOTHED BF PEAKS WITH TWO GAUSSIANS
# you have to have pretty decent guesses in the gausspars file for this to work.
bffitlist = bff.gaussparty(gausspars, nspec, filenamelist, bfnormlist, bf_ind, amplimits, threshold, widlimits)
rvraw1 = []; rvraw2 = []; rvraw1_err = []; rvraw2_err = []
rvraw1.append(0), rvraw2.append(0), rvraw1_err.append(0), rvraw2_err.append(0)
for i in range(1, len(bffitlist)):
    rvraw1.append(bffitlist[i][0][1]) # [0,1,2] is amp,rv,width for star 1; [4,5,6] is same for star2
    rvraw2.append(bffitlist[i][0][4])
    rvraw1_err.append(bffitlist[i][2][1])
    rvraw2_err.append(bffitlist[i][2][4])

# CALCULATE ORBITAL PHASES AND FINAL RV CURVE
rvdata = bff.rvphasecalc(bjdinfile, bjdoffset, nspec, period, BJD0, rvraw1, rvraw1_err, rvraw2, rvraw2_err, rvstd, bcvstd)
phase = rvdata[0]; bjdfunny = rvdata[1]
rv1 = rvdata[2]; rv2 = rvdata[3]
rv1_err = rvdata[4]; rv2_err = rvdata[5]
g2 = open(outfile, 'w')
print('# RVs calculated with BF_python.py', file=g2)
print('#', file=g2)
print('# Porb = {0} days, BJD0 = {1} days'.format(period, BJD0), file=g2)
print('# Wavelength axis = [{0} - {1}] Angstroms'.format(w1[0], w1[-1]), file=g2)
print('#', file=g2)
print('# Template spectrum (line 0 of infiles):  {0}'.format(filenamelist[0]), file=g2)
print('# RV of template, BCV of template (km/s): {0}, {1}'.format(rvstd, bcvstd), file=g2)
print('#', file=g2)
print('# List of all input spectra (infiles): {0}'.format(infiles), file=g2)
print('# Target BJD and BCV info (bjdinfile): {0}'.format(bjdinfile), file=g2)
print('# Gaussian fit guesses (gausspars):    {0}'.format(gausspars), file=g2)
print('#', file=g2)
print('# BF parameters: w00 = {0}; n = {1}; stepV = {2}'.format(w00, n, stepV), file=g2)
print('# BF parameters: smoothstd = {0}; m = {1}'.format(smoothstd, m), file=g2)
print('# gaussfit: amplimits = {0}; threshold = {1}, widlimits = {2}'.format(amplimits, threshold, widlimits), file=g2)
print('#', file=g2)
print('# time, phase, adjusted_time, RV1 [km/s], error1 [km/s], RV2 [km/s], error2 [km/s]', file=g2)
print('#', file=g2)
for i in range(1, nspec):
    print ('%.9f %.9f %.9f %.5f %.5f %.5f %.5f' % (bjdfunny[i] + bjdoffset, phase[i], bjdfunny[i], 
            rv1[i], rv1_err[i], rv2[i], rv2_err[i]), file=g2)
g2.close()
print('BJD, phase, and RVs written to %s.' % outfile)
print('Use rvplotmaker.py to plot the RV curve.')

try:
    bfout = open(bfoutfile, 'w')
    for idx in range(1, nspec):
        print('###', file=bfout)
        print('# timestamp: {0}'.format(datetimelist[idx]), file=bfout)
        print('# Gaussian 1 [amp, RV +/- err, wid]: [{0:.2f}, {1:.2f} +/- {2:.2f}, {3:.2f}]'.format(bffitlist[i][0][0], rvraw1[i], rvraw1_err[i], bffitlist[i][0][2]), file=bfout)
        print('# Gaussian 2 [amp, RV +/- err, wid]: [{0:.2f}, {1:.2f} +/- {2:.2f}, {3:.2f}]'.format(bffitlist[i][0][3], rvraw2[i], rvraw2_err[i], bffitlist[i][0][5]), file=bfout)
        print('# Uncorrected_RV, BF_amp, Gaussian_fit', file=bfout)
        print('###', file=bfout)
        for vel, amp, modamp in zip(bf_ind, bfsmoothlist[idx], bffitlist[idx][1]):
            print(vel, amp, modamp, file=bfout)
    bfout.close()
except:
    print('No BF outfile specified, not saving BF data to file')
    
# handy little gaussian function maker
def gaussian(x, amp, mu, sig): # i.e., (xarray, amp, rv, width)
    return amp * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

# PLOT THE FINAL SMOOTHED BFS + GAUSSIAN FITS IN INDIVIDUAL PANELS
# manually adjust this multi-panel plot based on how many spectra you have
#plt.figure(4)
windowcols = 3 # 4                             # how many columns the plot should have
#windowrows = 6                                # manually set number of plot rows here, or automatically below
windowrows = int([np.rint((nspec-1)/windowcols) if (np.float(nspec-1)/windowcols)%windowcols == 0 else np.rint((nspec-1)/windowcols)+1][0])
xmin = rvneg
xmax = rvpos
fig = plt.figure(1, figsize=(15,10))
fig.text(0.5, 0.04, 'Uncorrected Radial Velocity (km s$^{-1}$)', ha='center', va='center', size='large')
#########0.5, 0.04
fig.text(0.07, 0.6, 'Broadening Function', ha='center', va='center', size='large', rotation='vertical')
#########0.07, 0.5
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
    plt.text(xmax - 0.16*(np.abs(xmax-xmin)), 0.75*ymax, '%.3f $\phi$' % (phase[i]), size='small')
    plt.text(xmax - 0.26*(np.abs(xmax-xmin)), 0.55*ymax, '%s' % (datetimelist[i].iso[0:10]), size='small')
    #plt.plot(bf_ind, bfsmoothlist[i], color=colors[14], lw=1.5, ls='-', label='Smoothed BF')
    plt.plot(bf_ind, bfnormlist[i], color=colors[14], lw=1.5, ls='-', label='Normalized Smoothed BF')
    plt.plot(bf_ind, bffitlist[i][1], color=colors[0], lw=3, ls='-', label='Two Gaussian fit')
    gauss1 = gaussian(bf_ind, bffitlist[i][0][0], bffitlist[i][0][1], bffitlist[i][0][2])
    gauss2 = gaussian(bf_ind, bffitlist[i][0][3], bffitlist[i][0][4], bffitlist[i][0][5])
    plt.plot(bf_ind, gauss1, color=colors[6], lw=3, ls='--')#, label='Gaussian fit 1')
    plt.plot(bf_ind, gauss2, color=colors[2], lw=3, ls='--')#, label='Gaussian fit 2')
    # OPTION TO PLOT VERTICAL LINE AT ZERO
    #plt.axvline(x=0, color=colors[15])    
    # print legend
    if i==nspec-1: ax.legend(bbox_to_anchor=(2.6,0.7), loc=1, borderaxespad=0., 
                        frameon=False, handlelength=3, prop={'size':20})
plt.show()