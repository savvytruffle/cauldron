from __future__ import print_function

import numpy as np
from scipy.integrate import simps
from numpy import trapz
from sys import argv

#Calculate the Area underneath the BF curves 

starId = argv[1] 
KTRatio = float(argv[2])
KRsum = float(argv[3])

### Read in all the things! ###

areaout = 'data/6778289/6778289BFArea.txt'

#PAmp, Perr, PWidth, Samp, Serr, SWidth = np.loadtxt('data/5285607/5285607Gin.txt',
#	usecols=(0,1,2,3,4,5),unpack=True)

PAmp, Perr, PWidth, Samp, Serr, SWidth = np.loadtxt('data/6778289/6778289Gin.txt',
	usecols=(0,1,2,3,4,5),unpack=True)

PArea = (PAmp*PWidth)/(2.35*0.3984)
PAAve = np.mean(PArea)
SArea = (Samp*SWidth)/(2.35*0.3984)
SAAve = np.mean(SArea)

PAmed = np.median(PArea)
SAmed = np.median(PArea)

AreaRat = (PArea/SArea)
AvAreaRat = (PAAve/SAAve)
MedAreaRat= (SAmed/PAmed)
SoP = (SAAve/PAAve)
StD = np.std(SArea)/np.std(PArea)
#print(AvAreaRat)

dataout = np.vstack((PArea,SArea))
np.savetxt('6778289BFAreaout.txt', dataout.T, fmt = '%.5f')

#BFAreaout = open(areaout, 'w')
#print(PArea, file=BFAreaout)
    #print('Secondary BF Area: {0} +/- {1} width {2} xmax {3}'.format(bffitlist[idx][0][3], bffitlist[idx][2][3], bffitlist[idx][0][5], bffitlist[idx][0][4]), file=gout)

#BFAreaout.close()
    
#print(SAmed, PAmed)
#print('SArea = ', SArea)
#print('PArea = ', PArea)

print('S/P =', SoP, 'Std =', StD, 'Median =' ,MedAreaRat, 'Mean =', AvAreaRat)

RRat = 1/( np.sqrt(SoP/KTRatio ** 4) )

print('R_1/R_2 = ', RRat)

R2 = KRsum / (1 + SoP)

print('R2 = ', R2)

###Now, calculate log(g) to put on our HR diagrams from KEBLAT masses and individual radii just found###


