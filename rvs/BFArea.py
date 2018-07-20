from __future__ import print_function

import numpy as np
from scipy.integrate import simps
from numpy import trapz
from sys import argv
### BFArea.py
### calculate the Area underneath the BF curves of the primary and secondary members of a 
#   spectroscopic eclipsing binary. The ratio of these areas is directly proportional to 
#   the flux ratio of the binary (Bayless & Orosz 2006; Stassun et al. 2007).   


# python BFArea.py 

#starId = argv[1] 
#KTRatio = float(argv[2])
#KRsum = float(argv[3])

### Read in all the things! ###

#5285607
#areaout = 'data/5285607/5285607BFArea.txt'
#starId = 5285607
#KTRatio = 0.4922
#KTRaterr = 0.1086
#KRsum = 3.489
#KRsumerr = 0.051
#LCM1 = 1.554
#LCM1err = 0.023
#LCM2 = 1.333
#LCM2err = 0.020


#6864859
#areaout = 'data/6864859/6864859BFArea.txt'
#starId = 6864859
#KTRatio = 0.4965
#KTRaterr = 0.1143
#KRsum = 1.042
#KRsumerr = 0.016
#LCM1 = 1.354
#LCM1err = 0.029
#LCM2 = 1.411
#LCM2err = 0.028


#6778289
#areaout = 'data/6778289/6778289BFArea.txt'
#starId = 6778289
#KTRatio = 0.4865
#KTRaterr = 0.2673
#KRsum =0.7226
#KRsumerr = 0.0086
#LCM1 = 1.510
#LCM1err = 0.022
#LCM2 = 1.091
#LCM2err = 0.018


#4285607
#areaout = 'data/4285607/4285607BFArea.txt'
#starId = 4285607
#KTRatio = 0.969
#KTRaterr = 0.1405
#KRsum =0.9696
#KRsumerr = 0.0080
#LCM1 = 1.137
#LCM1err = 0.013
#LCM2 = 1.103
#LCM2err = 0.014


#6131659
#areaout = 'data/6131659/6131659BFArea.txt'
#starId = 6131659
#KTRatio = 0.6599
#KTRaterr = 0.2192
#KRsum =0.7458
#KRsumerr = 0.0052
#LCM1 = 0.9422
#LCM1err = 0.0093
#LCM2 = 0.7028
#LCM2err = 0.0078


#6781535
areaout = 'data/6781535/6781535BFArea.txt'
starId = 6781535
KTRatio = 0.5562
KTRaterr = 0.3508
KRsum = 1.0283
KRsumerr = 0.0213
LCM1 = 1.0057
LCM1err = 0.0327
LCM2 = 1.0346
LCM2err = 0.0330


#PAmp, Perr, PWidth, Samp, Serr, SWidth = np.loadtxt('data/5285607/5285607Gin.txt',
#	usecols=(0,1,2,3,4,5),unpack=True)

### Read in data from BF_python.py 
  # PAmp is the amplitude of the primary BF peak for each APOGEE visit spectra, similarly
  # Samp is the amplitude of the secondary BF peak for each APOGEE visit spectra. The PWidth
  # and SWidth are the width of the peaks, and the Perr and Serr are the associated errors. 
  
PAmp, Perr, PWidth, Samp, Serr, SWidth = np.loadtxt('data/6778289/6778289Gin.txt',
	usecols=(0,1,2,3,4,5),unpack=True)

### Calculate the area under the BF primary and secondary peaks for each visit, then take '
  # the average of all visits respective primary and secondary peaks, to find the average. 
  
PArea = (PAmp*PWidth)/(2.35*0.3984)
PAAve = np.mean(PArea)
SArea = (Samp*SWidth)/(2.35*0.3984)
SAAve = np.mean(SArea)

PAmed = np.median(PArea)
SAmed = np.median(PArea)

### So this ratio of areas for the BF primary and secondary peaks is proportional to the 
  # flux ratio of the binary.
  
AreaRat = (PArea/SArea)
AvAreaRat = (PAAve/SAAve)
MedAreaRat= (SAmed/PAmed)
SoP = (SAAve/PAAve)

### This is where the standard deviation of the SoP is calculated using the numpy package 
  # np.std which takes the standard deviation of the array elements, in this case SArea and
  # PArea. These standard deviations are then divided, to find the error in the Area ratio.
  
StD = np.std(SArea)/np.std(PArea)
#print(AvAreaRat)

### This is where the Primary and Secondary Areas were once saved.
#dataout = np.vstack((PArea,SArea))
#np.savetxt('6778289BFAreaout.txt', dataout.T, fmt = '%.5f')

#BFAreaout = open(areaout, 'w')
#print(PArea, file=BFAreaout)
    #print('Secondary BF Area: {0} +/- {1} width {2} xmax {3}'.format(bffitlist[idx][0][3], bffitlist[idx][2][3], bffitlist[idx][0][5], bffitlist[idx][0][4]), file=gout)

#BFAreaout.close()
    
#print(SAmed, PAmed)
#print('SArea = ', SArea)
#print('PArea = ', PArea)

print('S/P =', SoP, 'Std =', StD, 'Median =' ,MedAreaRat, 'Mean =', AvAreaRat)

RRat = 1/( np.sqrt(SoP/KTRatio ** 4) )

### This is where I attempt to calculate the error in the radius ratio, reasoning as follows:
  # For powers and square roots, you multiply the relative standard error by the power, in 
  # our ratio calculation this meant multiplying by powers of (-1), (1/2), and (4). The (-1) 
  # power comes from the 1 over the function, which is equivalent to the function to the (-1)
  # power. The (1/2) power comes from the square root. Finally the power of (4) comes from the 
  # fourth power in the denominator of the Rrat calculation, which I am not sure I have treated 
  # correctly in this case. 
  
RRatstd = -2 * np.absolute( SoP / KTRatio ) * np.sqrt ( (StD / SoP ) ** 2 + ( KTRaterr / KTRatio ) ** 2 - 2 *  ( ( StD * KTRaterr ) / ( SoP * KTRatio ) ))

print('R_1/R_2 = ', RRat, '+/- = ', RRatstd)

R2 = KRsum / (1 + SoP)
R2err = np.abs (R2) * ( np.sqrt ( ( KRsumerr / KRsum ) ** 2 + ( StD / SoP ) ** 2 ) )

R1 = KRsum - R2
R1err = np.sqrt( KRsumerr ** 2 + R2err ** 2 + 2 * KRsumerr * R2err )
print('R_1 = ', R1, '+/-', R1err, 'R2 = ', R2, '+/-', R2err)

###Now, calculate log(g) to put on our HR diagrams from KEBLAT masses and individual radii just found###

g1 = ( LCM1 / R1**2 ) * 27400
logg1 = np.log10(g1)
logg1err = np.log10( (LCM1err/LCM1)**2 + (2*(R1err/R1))**2)

g2 = ( LCM2 / R2**2 ) * 27400
logg2 = np.log10(g2)
logg2err = np.log10(np.std(g2))

print('logg1 = ', logg1, '+/-', logg1err, 'logg2 = ', logg2, '+/-', logg2err)
