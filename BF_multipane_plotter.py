import os
import numpy as np
import matplotlib.pyplot as plt

"""Plotting script for BF panels with multiple RV2 cases.

Needed Stuff
------------
Files named data/KICBFOut_ID.txt
where KIC is the star name and ID is from 1-N visits
the columns should be uncorrected RV (km/s), BF amplitude (arbitrary), Gaussian fit

File named data/KIC_bcv.txt
where KIC is the star name
the columns should be timestamp, BCV (km/s)


Originally by Diana, edited/adapted by Meredith
"""

maxVisit = 100  # there can't be any more visits than this
KIC = '6449358'
bfoutpath = os.path.join('data', KIC+'BFOut_')
bcvpath = os.path.join('data', KIC+'_bcv.txt')

rv2_case1 = None
rv2_case2 = None
rv2_case3 = None

bfoutdata = []
timestamps = []
gauss1 = []
gauss2 = []
gauss3 = []
for visit in range(maxVisit):
    try:
        filename = bfoutpath + '{}.txt'.format(visit+1)
        bfoutdata.append(np.loadtxt(filename))
    except FileNotFoundError:
        pass
    else:
        with open(filename) as f:
            for idx, line in enumerate(f):
                if idx == 1:  # 2nd line of file
                    timestamp = np.float(line[13:])
                    timestamps.append(timestamp)
                elif idx == 2:  # 3rd line of file
                    gauss1.append(line)
                elif idx == 3:  # 4th line of file
                    gauss2.append(line)
                elif idx == 4:  # 5th line of file
                    if line[2:12] == 'Gaussian 3':  # only if info for a 3rd Gaussian exists
                        gauss3.append(line)

# Not currently using gauss1, gauss2, or gauss3 for anything because the data in the file is wrong
# and it's a huge pain to do string parsing to get rid of the []s WHOSE IDEA WAS THAT MEREDITH

bfoutdata = np.array(bfoutdata)
timestamps_bcv, bcvs = np.loadtxt(bcvpath, usecols=(0,1), unpack=True)

print(timestamps, len(timestamps))
print(timestamps_bcv, len(timestamps_bcv))

plt.figure()

# TODO: read in and plot rv2_case1, case2, and case3
#       make subplots fancier
#       be clever about cross-checking the bfout timestamps and the bcv timestamps
#       fix bug that writes BFOutfile with the same damn gaussian info at the top of each chunk
#       plot locations of the fit peak(s) ... will need to fix ^^ first

for visit in range(len(timestamps)):                                                                                                                                                                 
    plt.subplot(4, 4, visit+1)
    plt.plot(bfoutdata[visit, :, 0], bfoutdata[visit, :, 1])
    #plt.text(-244, -0.016, r'$\phi$='+str((timestamps[visit]-keblat.tpe)%keblat.period/keblat.period)[:5])
    plt.xlim((-250, 250))
    #plt.axvline(rv2_case1[visit]*1e-3-bcvs[visit], color='green', ls='--', 
    #            label=str(np.round(rvpars_case1[0], 1))+', '+str(np.round(rvpars_case1[1], 1)))
    #plt.axvline(rv2_case2[visit]*1e-3-bcvs[visit], color='red', ls='--', 
    #            label=str(np.round(rvpars_case2[0], 1)) +', '+str(np.round(rvpars_case2[1], 1)))
    #plt.axvline(rv2_case3[visit]*1e-3-bcvs[visit], color='orange', ls='--', 
    #            label=str(np.round(rvpars_case3[0], 1)) +', '+str(np.round(rvpars_case3[1], 1)))
    plt.legend(loc='upper left')

plt.show()