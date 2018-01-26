from __future__ import print_function

import numpy as np
from scipy.integrate import simps
from numpy import trapz


### Read in all the things! ###

areaout = 'data/5285607/5285607BFArea.txt'

#PAmp, Perr, PWidth, Samp, Serr, SWidth = np.loadtxt('data/5285607/5285607Gin.txt',
#	usecols=(0,1,2,3,4,5),unpack=True)

PAmp, Perr, PWidth, Samp, Serr, SWidth = np.loadtxt('data/5285607/5285607Gin.txt',
	usecols=(0,1,2,3,4,5),unpack=True)

PArea = (PAmp*PWidth)/(2.35*0.3984)
PAAve = np.mean(PArea)
SArea = (Samp*SWidth)/(2.35*0.3984)
SAAve = np.mean(SArea)

AreaRat = (PArea/SArea)
AvAreaRat = (PAAve/SAAve)
test = (SAAve/PAAve)
print(AvAreaRat)

dataout = np.vstack((PArea,SArea))
np.savetxt('5285607BFAreaout.txt', dataout.T, fmt = '%.5f')

#BFAreaout = open(areaout, 'w')
#print(PArea, file=BFAreaout)
    #print('Secondary BF Area: {0} +/- {1} width {2} xmax {3}'.format(bffitlist[idx][0][3], bffitlist[idx][2][3], bffitlist[idx][0][5], bffitlist[idx][0][4]), file=gout)

#BFAreaout.close()
    
print('DONE! test =', test, 'Primary/Secondary ='AvAreaRat)