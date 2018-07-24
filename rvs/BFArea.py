from __future__ import print_function

import numpy as np
from scipy.integrate import simps
from numpy import trapz
from sys import argv
### BFArea.py
### Calculate the Area underneath the BF curves of the primary and secondary members of a 
#   spectroscopic eclipsing binary. The ratio of these areas is directly proportional to 
#   the flux ratio of the binary (Bayless & Orosz 2006; Stassun et al. 2007).
### Backs out the keblat Tratio from the keblat flux ratio and keblat radius ratio

starIds = [5285607, 6864859, 6778289, 6449358, 4285087, 6131659, 6781535] #KEPLER Input Catalog
kRsums =     [3.489, 3.104, 2.745, np.nan, 2.033, 1.5251, 2.0408]         #KEBLAT Radius Sums
kRsum_errs = [0.051, 0.016, 0.0086, np.nan, 0.0080, 0.0052, 0.0213]       #KEBLAT Radius Sum errors
kM1s =     [1.554, 1.354, 1.510, np.nan, 1.137, 0.9422, 1.0057]           #KEBLAT Mass_1
kM1_errs = [0.023, 0.029, 0.022, np.nan, 0.013, 0.0093, 0.0327]           #KEBLAT Mass_1 errors
kM2s =     [1.333, 1.411, 1.091, np.nan, 1.103, 0.7028, 1.0346]           #KEBLAT Mass_2
kM2_errs = [0.020, 0.028, 0.018, np.nan, 0.014, 0.0078, 0.0330]           #KEBLAT Mass_2 errors
kfluxRatios = [0.258, 1.407, 0.19138, 0, 0.901, 0.1483, 0.9201]           #KEBLAT Flux ratios 
kfluxRatioErrs = [0.046, 0.101, 2.6e-5, 0.080, 0.0017, 0.0524]            #KEBLAT Flux ratios errors 
kradRatios = [0.551, 1.149, 0.57093, np.nan, 0.969, 0.6799, 0.8641]       #KEBLAT Radius Ratios 
kradRatiosErrs = [0.048, 0.020, 0.013, np.nan, 0.0080, 0.0057, 0.0275]    #KEBLAT Radius Ratio errors
  
for starId, kRsum, kRsum_err, kM1, kM1_err, kM2, kM2_err, kfluxRatio, kfluxRatioErr, kradRatio, kradRatiosErr in zip(
    starIds, kRsums, kRsum_errs, kM1s, kM1_errs, kM2s, kM2_errs, kfluxRatios, kfluxRatioErrs, kradRatios, kradRatiosErrs):
    
    print(' ')
    print('###')
    print(starId)
    print('###')

### Read in data from BF_python.py 
  # PAmp is the amplitude of the primary BF peak for each APOGEE visit spectra, similarly
  # Samp is the amplitude of the secondary BF peak for each APOGEE visit spectra. The PWidth
  # and SWidth are the width of the peaks, and the Perr and Serr are the associated errors. 

    areaout = 'data/' + str(starId) + '/' + str(starId) + 'BFArea.txt'
    gaussfile = 'data/' + str(starId) + '/' + str(starId) + 'Gin.txt'

    PAmp, Perr, PWidth, Samp, Serr, SWidth = np.loadtxt(gaussfile, usecols=(0,1,2,3,4,5), unpack=True)

    ### Calculate the area under the BF primary and secondary peaks for each visit,
    ### then do all the statistics.
  
    PAreas = (PAmp*PWidth)/(2.35*0.3984)  # list of one value per visit
    SAreas = (Samp*SWidth)/(2.35*0.3984)  # list of one value per visit

    SoP = np.mean(SAreas)/np.mean(PAreas)  # secondary over primary ratio of mean BF areas
    
    Psum = 0
    Ssum = 0
    for PArea, SArea in zip(PAreas, SAreas):
        Pdiff = (PArea - np.mean(PAreas))**2
        Psum += Pdiff
        Sdiff = (SArea - np.mean(SAreas))**2
        Ssum += Sdiff
    
    P_std = np.sqrt(Psum / (len(PAreas) - 1))
    S_std = np.sqrt(Ssum / (len(SAreas) - 1))

    P_stderror = P_std / np.sqrt(len(PAreas))
    S_stderror = S_std / np.sqrt(len(SAreas))
    
    SoP_stderror = SoP * np.sqrt((P_stderror/np.mean(PAreas))**2 + (S_stderror/np.mean(SAreas))**2)

    print('BF Area 1 mean', np.mean(PAreas), '+/-', P_std)
    print('BF Area 2 mean', np.mean(SAreas), '+/-', S_std)
    print('Flux Ratio (2/1) =', SoP, '+/-', SoP_stderror)  # this is definitely correct, hooray!

### We back out the keblat Tratio from the keblat flux ratio and keblat radius ratio.
    kTemprat = (kfluxRatio / kradRatio**2)**(1/4)
    kTempraterr = np.sqrt(((kTemprat / 4 * kfluxRatio) * kradRatiosErr)**2 + (-kTemprat / 2 * kradRatiosErr) **2)
### The following lines were used to compare the temperature ratios calculated from the BF flux 
  # and KEBLAT radii, to the temperature ratios calculated from the KEBLAT flux and radius ratios.
#    BFTemprat = (SoP / ((R1/R2)**2))**(1/4)
#    Tempratdif = ((kTemprat - BFTemprat) / ((kTemprat + BFTemprat)/2)) * 100
#    print('Temperature Ratio from Keblat = ', kTemprat)
#    print('Temperature Ratio from BF =', BFTemprat)
#    print('So, the difference between them is', Tempratdif, '%')

### Next, we will calculate the Radii of the stars using the flux ratio from the BF (SoP) 
  # and the KEBLAT temperature ratio (kTemprat)

    BFKRadrat = 1 / (np.sqrt(SoP / kTemprat**4))
    BFKRadraterr = np.sqrt(((-1 / (2 * kTemprat**4 * (SoP / kTemprat**4))**(3/2)) * SoP_stderror)**2 + (2 / (kTemprat * (SoP / kTemprat**4)) * kTempraterr)**2)
    R2 = kRsum / (1 + BFKRadrat)
    R2err = R2 * np.sqrt((kRsum_err / kRsum)**2 + (SoP_stderror / SoP)**2)
    R1 = kRsum - R2
    R1err = np.sqrt(kRsum_err**2 + R2err**2)
    
    print('BFKRadrat = {0:.2f} +/- {1:.2f}'.format(BFKRadrat, BFKRadraterr))
#    print('WARNING these radius and logg values are WRONG (but the error propagation is right, dangit)')
    print('R1 = {0:.3f} +/- {1:.3f}, R2 = {2:.3f} +/- {3:.3f}'.format(R1, R1err, R2, R2err))
    

### Now, calculate log(g) to put on our HR diagrams from KEBLAT masses and individual radii just found###

    g1 = (kM1 / R1**2) * 27400
    logg1 = np.log10(g1)
    g1err = np.sqrt((27400/(R1**2) * kM1_err)**2 + ((-54800*kM1/(R1**3)) * R1err)**2)
    logg1err = 0.434*(g1err/g1)

    g2 = (kM2 / R2**2) * 27400
    logg2 = np.log10(g2)
    g2err = np.sqrt((27400/(R2**2) * kM2_err)**2 + ((-54800*kM2/(R2**3)) * R2err)**2)
    logg2err = 0.434*(g2err/g2)

    print('logg1 = {0:.3f} +/- {1:.3f}, logg2 = {2:.3f} +/- {3:.3f}'.format(logg1, logg1err, logg2, logg2err))
