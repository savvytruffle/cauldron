from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import RegularPolyCollection
from BF_functions import user_rc
'''
Radial velocity plotter!
Makes a plot with two panels: top is RV vs. time, bottom is RV vs. orbital phase

You need to have these columns in your input file: TIME, PHASE, RV1, RV1_ERR, RV2, RV2_ERR

Update September 2015:
Has two flag options
1. apply some shift to the RVs before plotting them (doShift)
2. read in another set of calculated RVs from cols 8,9,10,11 and plot RV1 = RV_col8-RV_col3
   and RV2 = RV_col10-RV_col5, plus a line at RV = 0 (compareRVs)
    
Update June 2016:
Simplified some options; no longer manually sets point shape as a function of "source" string.
If you want that functionality, use an older version of this code... it was messy.
**NOTE that any RV value with an error bar = 0 is not plotted!**
'''

dateoffset = 2454833. # this value will be subtracted from bjds in pane vs. time

#sysname = '5285607'; filename = '5285607Outfile_take2.txt'
#timestart = 980; timeend = 1020
#phasemin = 0.5; phasemax = 1.5
#RVmin = -45; RVmax = 180

sysname = '6449358'; filename = 'data/6449358/6449358Outfile.txt'
timestart = 1700; timeend = 2000
phasemin = 0.5; phasemax = 1.5
RVmin = 0; RVmax = 140

colors = user_rc()

#sysname = '5285607'; filename = 'data/5285607/5285607OutfileJC.txt'
#timestart = 975; timeend = 1040
#phasemin = 0.5; phasemax = 1.5
#RVmin = -50; RVmax = 200 # 5285607

#sysname = '6778289'; filename = 'data/6778289/6778289OutfileJC.txt'
#timestart = 1720; timeend = 1990
#phasemin = 0.5; phasemax = 1.12
#RVmin = -50; RVmax = 150 # 6778289

#sysname = '6864859'; filename = 'data/6864859/6864859OutfileJC.txt'
#timestart = 1720; timeend = 1990
#phasemin = 0.53; phasemax = 1.53
#RVmin = 40; RVmax = 140 # 6864859

# Other useful definitions
red = '#e34a33' # red, star 1
yel = '#fdbb84' # yellow, star 2

# usecols=(0,1,3,4,5,6) # this is the default, with RVs in 3,4,5,6 not 8,9,10,11
bjd, phase, rv1, rverr1, rv2, rverr2 = np.loadtxt(filename, comments='#', 
    dtype={'names': ('bjd', 'phase', 'rv1', 'rverr1', 'rv2', 'rverr2'),
    'formats': (np.float64, np.float64, np.float64, np.float64, np.float64, np.float64)},
    usecols=(0,1, 3,4,5,6), unpack=True)

# Skip any RV values that have 0 for error bars
for idx, err in enumerate(rverr1):
    if err == 0:
        rv1[idx] = None
        rverr1[idx] = None
for idx, err in enumerate(rverr2):
    if err == 0:
        rv2[idx] = None
        rverr2[idx] = None
rv1mask = np.isfinite(rv1)
rv2mask = np.isfinite(rv2)

# Double the arrays so we can plot any phase from 0 to phase 2... assuming phase is in range (0,1)
rv1_double = np.concatenate((rv1,rv1), axis=0)
rv2_double = np.concatenate((rv2,rv2), axis=0)
phase_double = np.concatenate((np.array(phase),np.array(phase)+1.0), axis=0)
rverr1_double = np.concatenate((rverr1,rverr1), axis=0)
rverr2_double = np.concatenate((rverr2,rverr2), axis=0)

# Set up the figure
fig = plt.figure(1, figsize=(13,9))

# Unfolded RV vs time (BJD-2454833)
ax2 = plt.subplot(2,1,1)
plt.axis([timestart, timeend, RVmin, RVmax])
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.xaxis.set_ticks_position('bottom')
ax2.yaxis.set_ticks_position('left')
plt.tick_params(axis='both', which='major')
# dotted lines to guide the eye
plt.plot(bjd[rv1mask]-dateoffset, rv1[rv1mask], color=colors[15], mfc=None, mec=None, lw=1.5, ls=':')
plt.plot(bjd[rv2mask]-dateoffset, rv2[rv2mask], color=colors[15], mfc=None, mec=None, lw=1.5, ls=':')
for idx, date in enumerate(bjd):
    plt.errorbar(date-dateoffset, rv1[idx], yerr=rverr1[idx], fmt='ko', color=colors[15], mfc=colors[6], mec=colors[14], ms=10, lw=1.5)
    plt.errorbar(date-dateoffset, rv2[idx], yerr=rverr2[idx], fmt='ko', color=colors[15], mfc=colors[2], mec=colors[14], ms=10, lw=1.5)
plt.xlabel("Time (BJD -- {0:.0f})".format(dateoffset))

# Folded RV vs phase
ax1 = plt.subplot(2,1,2)
plt.axis([phasemin, phasemax, RVmin, RVmax])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
plt.tick_params(axis='both', which='major')
for idx, ph in enumerate(phase_double):
    plt.errorbar(phase_double[idx], rv1_double[idx], yerr=rverr1_double[idx], marker='o', color=colors[6], mec=colors[14], ecolor=colors[6], ms=10, ls='None', lw=1.5)
    plt.errorbar(phase_double[idx], rv2_double[idx], yerr=rverr2_double[idx], marker='o', color=colors[2], mec=colors[14], ecolor=colors[2], ms=10, ls='None', lw=1.5)
plt.xlabel("Orbital Phase")

# Draw vertical lines at phase = 0.5
#plt.axvline(x=0.5, ymin=-59, ymax=45, color='k', ls=':')
#plt.axvline(x=1.5, ymin=-59, ymax=45, color='k', ls=':')

# Option for a legend and labels (note: for a legend you will need to add a label to the plt.errorbar commands)
#plt.legend(ncol=2, loc=1, numpoints=1, frameon=False, bbox_to_anchor=(1,2.35), columnspacing=0.7)
fig.text(0.07, 0.5, 'Radial Velocity (km s$^{-1}$)', ha='center', va='center', rotation='vertical', size='large')
fig.text(0.14, 0.13, 'Folded')
fig.text(0.14, 0.55, 'Unfolded')
fig.text(0.14, 0.9, sysname, size='large')

plt.show()
#fig.savefig('6449358rv.png')