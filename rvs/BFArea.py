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

starIds = [5285607, 6864859, 6778289, 6449358, 4285087, 6131659, 6781535]
kRsums =     [3.489, 3.104, 2.745, np.nan, 2.033, 1.5251, 2.0408]
kRsum_errs = [0.051, 0.016, 0.0086, np.nan, 0.0080, 0.0052, 0.0213]
kM1s =     [1.554, 1.354, 1.510, np.nan, 1.137, 0.9422, 1.0057]
kM1_errs = [0.023, 0.029, 0.022, np.nan, 0.013, 0.0093, 0.0327]
kM2s =     [1.333, 1.411, 1.091, np.nan, 1.103, 0.7028, 1.0346]
kM2_errs = [0.020, 0.028, 0.018, np.nan, 0.014, 0.0078, 0.0330]

### Read in data from BF_python.py 
  # PAmp is the amplitude of the primary BF peak for each APOGEE visit spectra, similarly
  # Samp is the amplitude of the secondary BF peak for each APOGEE visit spectra. The PWidth
  # and SWidth are the width of the peaks, and the Perr and Serr are the associated errors. 
  
for starId, kRsum, kRsum_err, kM1, kM1_err, kM2, kM2_err in zip(
    starIds, kRsums, kRsum_errs, kM1s, kM1_errs, kM2s, kM2_errs):
    
    print(' ')
    print('###')
    print(starId)
    print('###')

    areaout = 'data/' + str(starId) + '/' + str(starId) + 'BFArea.txt'
    gaussfile = 'data/' + str(starId) + '/' + str(starId) + 'Gin.txt'

    PAmp, Perr, PWidth, Samp, Serr, SWidth = np.loadtxt(gaussfile, usecols=(0,1,2,3,4,5), unpack=True)

    ### Calculate the area under the BF primary and secondary peaks for each visit,
    ### then do all the statistics.
  
    PAreas = (PAmp*PWidth)/(2.35*0.3984)  # list of one value per visit
    SAreas = (Samp*SWidth)/(2.35*0.3984)

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

    # Forging ahead to calculate radii and logg... 

    R2 = kRsum / (1 + SoP)  # this is NOT correct... !
    #RRatio = 1 / (np.sqrt(SoP / TRatio**4))
    #R2 = kRsum / (1 + RRatio)  # this IS correct, but we need Tratio from HR.py to have RRatio!
    R2err = R2 * np.sqrt((kRsum_err / kRsum)**2 + (SoP_stderror / SoP)**2)

    R1 = kRsum - R2
    R1err = np.sqrt(kRsum_err**2 + R2err**2)
    print('WARNING these radius and logg values are WRONG (but the error propagation is right, dangit)')
    print('R1 = {0:.3f} +/- {1:.3f}, R2 = {2:.3f} +/- {3:.3f}'.format(R1, R1err, R2, R2err))

    ###Now, calculate log(g) to put on our HR diagrams from KEBLAT masses and individual radii just found###

    g1 = (kM1 / R1**2) * 27400
    logg1 = np.log10(g1)
    g1err = np.sqrt((27400/(R1**2) * kM1_err)**2 + ((-54800*kM1/(R1**3)) * R1err)**2)
    logg1err = 0.434*(g1err/g1)

    g2 = (kM2 / R2**2) * 27400
    logg2 = np.log10(g2)
    g2err = np.sqrt((27400/(R2**2) * kM2_err)**2 + ((-54800*kM2/(R2**3)) * R2err)**2)
    logg2err = 0.434*(g2err/g2)

    print('logg1 = {0:.3f} +/- {1:.3f}, logg2 = {2:.3f} +/- {3:.3f}'.format(logg1, logg1err, logg2, logg2err))
