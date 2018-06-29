from __future__ import print_function

import numpy as np
from scipy.integrate import simps
from numpy import trapz
from sys import argv

# python BFArea.py 
#Calculate the Area underneath the BF curves 

#starId = argv[1] 
#KTRatio = float(argv[2])
#KRsum = float(argv[3])

### Read in all the things! ###

#5285607
areaout = 'data/5285607/5285607BFArea.txt'
starId = 5285607
KTRatio = 0.4922
KRsum = 3.489
LCM1 = 1.554
LCM1err = 0.023
LCM2 = 1.333
LCM2err = 0.020

'''
#6864859
areaout = 'data/6864859/6864859BFArea.txt'
starId = 6864859
KTRatio = 0.4965
KRsum = 3.104
LCM1 = 1.354
LCM1err = 0.029
LCM2 = 1.411
LCM2err = 0.028

#6778289
areaout = 'data/6778289/6778289BFArea.txt'
starId = 6778289
KTRatio = 0.4865
KRsum =2.745
LCM1 = 1.510
LCM1err = 0.022
LCM2 = 1.091
LCM2err = 0.018

#4285607
areaout = 'data/4285607/4285607BFArea.txt'
starId = 4285607
KTRatio = 0.969
KRsum =2.033
LCM1 = 1.137
LCM1err = 0.013
LCM2 = 1.103
LCM2err = 0.014

#6131659
areaout = 'data/6131659/6131659BFArea.txt'
starId = 6131659
KTRatio = 0.6599
KRsum =1.5251
LCM1 = 0.9422
LCM1err = 0.0093
LCM2 = 0.7028
LCM2err = 0.0078

#6781535
areaout = 'data/6781535/6781535BFArea.txt'
starId = 6781535
KTRatio = 0.5562
KRsum = 2.0408
LCM1 = 1.0057
LCM1err = 0.0327
LCM2 = 1.0346
LCM2err = 0.0330
'''

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

R1 = KRsum - R2

print('R_1 = ', R1, 'R2 = ', R2)

###Now, calculate log(g) to put on our HR diagrams from KEBLAT masses and individual radii just found###

g1 = ( LCM1 / R1**2 ) * 27400
logg1 = np.log10(g1)
logg1err = np.log10( (LCM1err/LCM1)**2 + (2*(R1err/R1))**2)

g2 = ( LCM2 / R2**2 ) * 27400
logg2 = np.log10(g2)
logg2err = np.log10(np.std(g2))

print('logg1 = ', logg1, '+/-', log1err, 'logg2 = ', logg2, '+/-', logg2err)
