import os
import numpy as np
import matplotlib.pyplot as plt

"""Plotting script for BF panels, with optional multiple RV2 cases.

Needed Stuff
------------
Files named data/KICBFOut_ID.txt
where KIC is the star name and ID is from 1-N visits
the columns should be uncorrected RV (km/s), BF amplitude, Gaussian fit

File named data/KIC_bcv.txt
where KIC is the star name
the columns should be timestamp, BCV (km/s)

Optional Stuff
--------------
File named data/KIC_fowardmodel.rvs
where KIC is the star name
the columns should be three sets of RV2s (m/s)


Originally by Diana, edited/adapted by Meredith
"""

maxVisit = 100  # there can't be any more visits than this
KIC = '6449358'
bfoutpath = os.path.join('data', KIC+'BFOut_visit')
bcvpath = os.path.join('data', KIC+'_bcv.txt')
multipleRV2data = os.path.join('data', KIC+'_forwardmodel.rvs')

multipleRV2s = True
try:
    rv2_case1, rv2_case2, rv2_case3 = np.loadtxt(multipleRV2data, usecols=(0,1,2), unpack=True)
except FileNotFoundError:
    multipleRV2s = False

bfoutdata = []
timestamps = []
gaussfit1 = []
gaussfit2 = []
gaussfit3 = []
for visit in range(maxVisit):
    try:
        filename = bfoutpath + '{}.txt'.format(visit+1)
        bfoutdata.append(np.loadtxt(filename))
    except FileNotFoundError:
        pass
    else:
        with open(filename) as f:
            for idx, line in enumerate(f):
                if idx == 0:  # 1st line of file
                    timestamp = np.float(line[15:].strip())
                    timestamps.append(timestamp)
                elif idx == 1:
                    gauss1 = line.split()
                elif idx == 2:
                    gauss2 = line.split()
                elif idx == 3:
                    if line[2:12] == 'Gaussian 3':  # only if info for a 3rd Gaussian exists
                        gauss3 = line.split()
            gaussfit1.append([float(value) for value in gauss1[9:13]])
            gaussfit2.append([float(value) for value in gauss2[9:13]])
            if gauss3:
                gaussfit3.append([float(value) for value in gauss3[9:13]])

#print(gaussfit1)  # four values: amplitude, RV, RVerr, width
#print(gaussfit2)

bfoutdata = np.array(bfoutdata)
timestamps_bcv, bcvs = np.loadtxt(bcvpath, usecols=(0,1), unpack=True)

#print(timestamps, len(timestamps))
#print(timestamps_bcv, len(timestamps_bcv))

plt.figure()

# TODO: 
#       make subplots fancier, add labels, adjust font sizes, etc.
#       be clever about cross-checking the bfout timestamps and the bcv timestamps
#       plot locations of the fit peak(s) ... use gaussfit1 and gaussfit2 !

for visit in range(len(timestamps)):                                                                                                                                                                 
    plt.subplot(4, 4, visit+1)
    plt.plot(bfoutdata[visit, :, 0], bfoutdata[visit, :, 1])
    #plt.text(-244, -0.016, r'$\phi$='+str((timestamps[visit]-keblat.tpe)%keblat.period/keblat.period)[:5])
    plt.xlim((-250, 250))
    if multipleRV2s:
        plt.axvline(rv2_case1[visit]*1e-3-bcvs[visit], color='green', ls='--')
                    #label=str(np.round(rvpars_case1[0], 1))+', '+str(np.round(rvpars_case1[1], 1)))
        plt.axvline(rv2_case2[visit]*1e-3-bcvs[visit], color='red', ls='--')
                    #label=str(np.round(rvpars_case2[0], 1)) +', '+str(np.round(rvpars_case2[1], 1)))
        plt.axvline(rv2_case3[visit]*1e-3-bcvs[visit], color='orange', ls='--')
                    #label=str(np.round(rvpars_case3[0], 1)) +', '+str(np.round(rvpars_case3[1], 1)))
    plt.legend(loc='upper left')

plt.show()