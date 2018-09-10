from __future__ import print_function

import numpy as np
from scipy.integrate import simps
from numpy import trapz
from sys import argv

##################################### BFArea.py ##########################################
'''Calculate the Area underneath the BF curves of the primary and secondary members of a''' 
'''spectroscopic eclipsing binary. The ratio of these areas is directly proportional to '''
'''the flux ratio of the binary (Bayless & Orosz 2006; Stassun et al. 2007). This version'''
'''of BFArea.py also backs out the KEBLAT Tratio from the keblat flux ratio and KEBLAT ''' 
'''radius ratio and compares this Tratio to the Tratio calculated using the BF flux    '''
'''and the KEBLAT radius ratio.                                                        '''
##########################################################################################


starIds = [5285607, 6864859, 6778289, 6449358, 4285087, 6131659, 6781535]    #KEPLER Input Catalog
ASPCAPTeffs = [6495, 6417, 6572, 6237, 5664, 4845, 5749]
ASPCAPTeff_errs = [156, 159, 162, 179, 146, 98, 125]
kRsums =     [3.489, 3.104, 2.745, np.nan, 2.033, 1.5251, 2.0408]            #KEBLAT Radius Sums [R_sun]
kRsum_errs = [0.051, 0.016, 0.0086, np.nan, 0.0080, 0.0052, 0.0213]          #KEBLAT Radius Sum errors
kM1s =     [1.554, 1.354, 1.510, np.nan, 1.137, 0.9422, 1.0057]              #KEBLAT Mass_1
kM1_errs = [0.023, 0.029, 0.022, np.nan, 0.013, 0.0093, 0.0327]              #KEBLAT Mass_1 errors
kM2s =     [1.333, 1.411, 1.091, np.nan, 1.103, 0.7028, 1.0346]              #KEBLAT Mass_2
kM2_errs = [0.020, 0.028, 0.018, np.nan, 0.014, 0.0078, 0.0330]              #KEBLAT Mass_2 errors
kfluxRatios = [0.258, 1.407, 0.19138, np.nan, 0.901, 0.1483, 0.9201]         #KEBLAT Flux ratios 
kfluxRatioErrs = [0.046, 0.101, 2.6e-5, np.nan, 0.080, 0.0017, 0.0524]       #KEBLAT Flux ratios errors 
kradRatios = [0.551, 1.149, 0.57093, np.nan, 0.969, 0.6799, 0.8641]          #KEBLAT Radius Ratios 
kradRatiosErrs = [0.048, 0.020, 0.013, np.nan, 0.0080, 0.0057, 0.0275]       #KEBLAT Radius Ratio errors
GAIAparallaxs = [1.2504, 1.4897, 0.9093, 1.1974, 1.619, np.nan, np.nan]
GAIAparallax_errs = [0.0216, 0.0241, 0.0222, 0.0264, np.nan, np.nan]
  
for starId, ASPCAPTeff, ASPCAPTeff_err, kRsum, kRsum_err, kM1, kM1_err, kM2, kM2_err, kfluxRatio, kfluxRatioErr, kradRatio, kradRatiosErr, GAIAparallax, GAIAparallax_err in zip(
    starIds, ASPCAPTeffs, ASPCAPTeff_errs, kRsums, kRsum_errs, kM1s, kM1_errs, kM2s, kM2_errs, kfluxRatios, kfluxRatioErrs, kradRatios, kradRatiosErrs, GAIAparallaxs, GAIAparallax_errs):
    
    print(' ')
    print('##############################################################')
    print('Running analysis for star', starId)

#### Read in data from BF_python.py 
   # PAmp is the amplitude of the primary BF peak for each APOGEE visit spectra, similarly
   # Samp is the amplitude of the secondary BF peak for each APOGEE visit spectra. The PWidth
   # and SWidth are the width of the peaks, and the Perr and Serr are the associated errors. 

    areaout = 'data/' + str(starId) + '/' + str(starId) + 'BFArea.txt'
    gaussfile = 'data/' + str(starId) + '/' + str(starId) + 'Gin.txt'

    PAmp, Perr, PWidth, Samp, Serr, SWidth = np.loadtxt(gaussfile, usecols=(0,1,2,3,4,5), unpack=True)

#### Calculate the area under the BF primary and secondary peaks for each visit ####
   # The areas under the primary or secondary are (respectively) roughly equal for each visit,
   # the mean is taken for the primary and secondary (respectively)  #### 
   
    PAreas = (PAmp*PWidth)/(2.35*0.3984)                                    #list of one value per visit (primary)
    SAreas = (Samp*SWidth)/(2.35*0.3984)                                    #list of one value per visit (secondary)

    SoP = np.mean(SAreas)/np.mean(PAreas)  # secondary over primary ratio of mean BF areas
    
#### Calculate the std error in both the primary and secondary area averages, then propagate 
   # through to the SoP ratio.
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
    print('BF Flux Ratio (2/1) =', SoP, '+/-', SoP_stderror) 

#### We back out the keblat Tratio from the KEBLAT flux ratio and KEBLAT radius ratio.
    
#    KEBLATTemprat = (kfluxRatio / kradRatio**2)**(1/4)
#    kTempraterr = np.sqrt(((KEBLATTemprat / 4 * kfluxRatio) * kradRatiosErr)**2 + (-KEBLATTemprat / 2 * kradRatiosErr) **2)

#### Calculate R1 and R2 using the KEBLAT radius sum and the KEBLAT radius ratio

    R2 = kRsum / (1 + kradRatio)
    R2err = R2 * np.sqrt((kRsum_err / kRsum)**2)
    R1 = kRsum - R2
    R1err = np.sqrt(kRsum_err**2 + R2err**2)
    R2oR1 = R2/R1
    R1oR2 = R1/R2
    R2oR1err = R2oR1 * np.sqrt((R2err / R2)**2 + (R1err / R1)**2)
    R1oR2err = R1oR2 * np.sqrt((R2err / R2)**2 + (R1err / R1)**2)

#### Now, we find the sum of the KEBLAT fluxes using the ASPCAP temperature, the sum of the
   # squared radii and the distance to our targets from GAIA
    solartocm = 6.599e10
    pctocm = 3.086e18
    sigma = 5.6704e-5 
    GAIAdistance = (1 / GAIAparallax) * pctocm
#    GAIAdistance_err = np.sqrt((                                              # [erg*cm^-2*K-4]
    KEBLATFluxsum = ((sigma * ASPCAPTeff **4 * ((R1*solartocm)**2 + (R2*solartocm)**2)) / GAIAdistance **2)
#    KEBLATFluxsum = ((sigma * ASPCAPTeff **4 * ((R1)**2 + (R2)**2)) / (GAIAdistance*pctom) **2)
    Flux1KEBASP = (sigma * ASPCAPTeff **4 * 4 * np.pi * (R1*solartocm)**2) / GAIAdistance
    Flux2KEBASP = (sigma * ASPCAPTeff **4 * 4 * np.pi * (R2*solartocm)**2) / GAIAdistance
    F2overF1 = Flux1KEBASP/Flux2KEBASP
    
#    print('BFKRadrat = {0:.2f} +/- {1:.2f}'.format(BFKRadrat, BFKRadraterr))
    print('R1 = {0:.3f} +/- {1:.3f}, R2 = {2:.3f} +/- {3:.3f}'.format(R1, R1err, R2, R2err))
    print('R2/R1 = {0:.3f}, +/- {1:.3f}, R1/R2 = {2:.3f}, +/- {3:.3f}'.format(R2oR1, R2oR1err, R1oR2, R1oR2err))
    print('Flux1KEBASP = {0:.3f} Flux1KEBASP = {1:.3f} F2/F1 = {2:.3f}'.format(Flux1KEBASP, Flux2KEBASP, F2overF1))
    
#### Now, we will use the fluxes we just found, and assuming the ASPCAP temperature is the flux weighted sum 
   # of the system, and that for each system the primary contributes the majority of the light, we can 
   # calculate the temperatures of the components. 
    R1cm = R1 * solartocm
    R2cm = R2 * solartocm
    T1KEBASP = ((Flux1KEBASP * GAIAdistance**2) / (sigma * R2cm**2))**(1/4)
    T2KEBASP = ((Flux2KEBASP * GAIAdistance**2) / (sigma * R2cm**2))**(1/4)
    
#    print(KEBLATFluxsum, GAIAdistance, R1cm, R1cm)
    print('T1KEBASP = ', T1KEBASP,'T2KEBASP = ', T2KEBASP)

   
#### The following lines were used to compare the temperature ratios calculated from the BF flux 
   # and KEBLAT radii, to the temperature ratios calculated from the KEBLAT flux and radius ratios.
#    BFTemprat = (SoP / ((R1/R2)**2))**(1/4)
#    Tempratdif = ((kTemprat - BFTemprat) / ((kTemprat + BFTemprat)/2)) * 100
#    print('Temperature Ratio from Keblat = ', KEBLATTemprat)
#    print('Temperature Ratio from BF =', BFTemprat)
#    print('So, the difference between them is', Tempratdif, '%')
    
    ### Now, calculate log(g) to put on our HR diagrams from KEBLAT masses and individual radii just found###
    g1 = (kM1 / R1**2) * 27400
    logg1 = np.log10(g1)
    g1err = np.sqrt((27400/(R1**2) * kM1_err)**2 + ((-54800*kM1/(R1**3)) * R1err)**2)
    logg1err = 0.434*(g1err/g1)

    g2 = (kM2 / R2**2) * 27400
    logg2 = np.log10(g2)
    g2err = np.sqrt((27400/(R2**2) * kM2_err)**2 + ((-54800*kM2/(R2**3)) * R2err)**2)
    logg2err = 0.434*(g2err/g2)

    print('logg1 = {0:.6f} +/- {1:.6f}, logg2 = {2:.6f} +/- {3:.6f}'.format(logg1, logg1err, logg2, logg2err))
    
#### Next, we will confirm that the ASPCAP temperature is the flux weighted sum of the temperature components:
