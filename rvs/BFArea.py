from __future__ import print_function
import numpy as np
from scipy.integrate import simps
from numpy import trapz
from sys import argv
from astropy import units as u
from astropy import constants as const

##########################################################################################
##################################### BFArea.py ##########################################
'''Calculate the Area underneath the BF curves of the primary and secondary members of a''' 
'''spectroscopic eclipsing binary. The ratio of these areas is directly proportional to '''
'''the flux ratio of the binary (Bayless & Orosz 2006; Stassun et al. 2007). '''
##########################################################################################
##########################################################################################

starIds = [5285607, 6864859, 6778289, 6449358, 4285087, 6131659, 6781535]                 #KEPLER Input Catalog
ASPCAPTeffs = [6495, 6417, 6572, 6237, 5664, 4845, 5749] * u.K                            #ASPCAP Effective Temperature
ASPCAPTeff_errs = [156, 159, 162, 179, 146, 98, 125]                                      #Error on the above 
kRsums =     [3.489, 3.104, 2.745, np.nan, 2.033, 1.5251, 2.0408] * u.Rsun                #KEBLAT Radius Sums [R_sun]
kRsum_errs = [0.051, 0.016, 0.0086, np.nan, 0.0080, 0.0052, 0.0213]                       #KEBLAT Radius Sum errors
kM1s =     [1.554, 1.354, 1.510, np.nan, 1.137, 0.9422, 1.0057] * u.Msun                  #KEBLAT Mass_1
kM1_errs = [0.023, 0.029, 0.022, np.nan, 0.013, 0.0093, 0.0327]                           #KEBLAT Mass_1 errors
kM2s =     [1.333, 1.411, 1.091, np.nan, 1.103, 0.7028, 1.0346] * u.Msun                  #KEBLAT Mass_2
kM2_errs = [0.020, 0.028, 0.018, np.nan, 0.014, 0.0078, 0.0330]                           #KEBLAT Mass_2 errors
kfluxRatios = [0.258, 1.407, 0.19138, np.nan, 0.901, 0.1483, 0.9201]                      #KEBLAT Flux ratios 
kfluxRatioErrs = [0.046, 0.101, 2.6e-5, np.nan, 0.080, 0.0017, 0.0524]                    #KEBLAT Flux ratios errors 
kradRatios = [0.551, 1.149, 0.57093, np.nan, 0.969, 0.6799, 0.8641]                       #KEBLAT Radius Ratios 
kradRatiosErrs = [0.048, 0.020, 0.013, np.nan, 0.0080, 0.0057, 0.0275]                    #KEBLAT Radius Ratio errors
GAIAdistances = [781.979, 658.706, 1066.444, 815.946, 607.188, 3527.290, np.nan] * u.pc   #GAIA distances 
GAIAdistance_errs = [0.0216, 0.0241, 0.0222, 0.0264, np.nan, np.nan]                      #Error on the above

for starId, ASPCAPTeff, ASPCAPTeff_err, kRsum, kRsum_err, kM1, kM1_err, kM2, kM2_err, kfluxRatio, kfluxRatioErr, kradRatio, kradRatiosErr, GAIAdistance, GAIAdistance_err in zip(
    starIds, ASPCAPTeffs, ASPCAPTeff_errs, kRsums, kRsum_errs, kM1s, kM1_errs, kM2s, kM2_errs, kfluxRatios, kfluxRatioErrs, kradRatios, kradRatiosErrs, GAIAdistances, GAIAdistance_errs):

    print(' ')
    print('##############################################################')
    print('Running analysis for star', starId)

############################# Read in data from BF_python.py #############################
###PAmp is the amplitude of the primary BF peak for each APOGEE visit spectra similarly###
###Samp is the amplitude of the secondary BF peak for each APOGEE visit spectra PWidth ###
###and SWidth are the width of the peaks and the Perr and Serr are the associated error###
##########################################################################################

    areaout = 'data/' + str(starId) + '/' + str(starId) + 'BFArea.txt'
    gaussfile = 'data/' + str(starId) + '/' + str(starId) + 'Gin.txt'
    PAmp, Perr, PWidth, Samp, Serr, SWidth = np.loadtxt(gaussfile, usecols=(0,1,2,3,4,5), unpack=True)

############################ Calculate the area under the BF #############################
###The areas under the primary or secondary are (respectively) roughly equal for each  ###
### visit. The mean is taken for the primary and secondary (respectively)Then calculate###
### the std error in both the primary and secondary area averages, then propagate thru ###
### to the SoP ratio.                                                                  ###
##########################################################################################

    PAreas = (PAmp*PWidth)/(2.35*0.3984)           #one value per visit (primary)
    SAreas = (Samp*SWidth)/(2.35*0.3984)           #one value per visit (secondary)
    SoP = np.mean(SAreas)/np.mean(PAreas)          #secondary over primary ratio of mean BF areas
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

################################## Calculate R1 and R2 ###################################
### Using the KEBLAT radius sum and the KEBLAT radius ratio calculate R1 and R2 and its###
### error propagation                                                                  ###
##########################################################################################

    R2 = kRsum / (1 + kradRatio)
    R2err = R2 * np.sqrt((kRsum_err / kRsum)**2)
    R1 = kRsum - R2
    R1err = np.sqrt(kRsum_err**2 + R2err**2)
    R2oR1 = R2/R1
    R1oR2 = R1/R2
    R2oR1err = R2oR1 * np.sqrt((R2err / R2)**2 + (R1err / R1)**2)
    R1oR2err = R1oR2 * np.sqrt((R2err / R2)**2 + (R1err / R1)**2)

######################### Calculate the sum of the KEBLAT fluxes #########################
### Using the ASPCAP temperature, the sum of the squared radii and the distance to our ###
### targets from GAIA(Bailer-Jones et al 2018)and follow through with error propagation###
##########################################################################################

    R1.to(u.m)
    R2.to(u.m)
    GAIAdistance.to(u.m)
    KEBLATFluxsum = (((const.sigma_sb * ASPCAPTeff **4 * ((R1)**2 + (R2)**2)) / (GAIAdistance) **2))
    KEBLATFlux2 = KEBLATFluxsum / (1 + kfluxRatio)
    KEBLATFlux1 = KEBLATFluxsum - KEBLATFlux2
    KEBLATFluxsum = ((const.sigma_sb * ASPCAPTeff **4 * ((R1)**2 + (R2)**2)) / (GAIAdistance) **2)
    F2overF1 = KEBLATFlux2/KEBLATFlux1
    
    #print('BFKRadrat = {0:.2f} +/- {1:.2f}'.format(BFKRadrat, BFKRadraterr))
    print('R1 = {0:.3f} +/- {1:.3f}, R2 = {2:.3f} +/- {3:.3f}'.format(R1, R1err, R2, R2err))
    print('R2/R1 = {0:.3f}, +/- {1:.3f}, R1/R2 = {2:.3f}, +/- {3:.3f}'.format(R2oR1, R2oR1err, R1oR2, R1oR2err))
    print('KEBLATFlux1 = {0:.3f} KEBLATFlux2 = {1:.3f} F2/F1 = {2:.3f}'.format(KEBLATFlux1, KEBLATFlux2, F2overF1))
    #print('KEBLATFlux1 =', KEBLATFlux1, 'KEBLATFlux2 = ', KEBLATFlux2)

################################### Calculate T1 and T2 ##################################
###Now we will use the fluxes we just found, and assuming the ASPCAP temperature is the###
###flux weighted sum of the system and that for each system the primary contributes the###
###majority of the light, we can calculate the temperatures of the components.         ###
##########################################################################################

    T1KEBASP = (KEBLATFlux1 * ((GAIAdistance)**2) / (const.sigma_sb * (R1)**2))**(1/4)
    T2KEBASP = (KEBLATFlux2 * ((GAIAdistance)**2) / (const.sigma_sb * (R2)**2))**(1/4)
    T2oT1KEBASP = (T2KEBASP/T1KEBASP)
    print('T1KEBASP = {0:.3f} T2KEBASP = {1:.3f} T2oT1KEBASP = {2:.3f}'.format(T1KEBASP, T2KEBASP, T2oT1KEBASP))
    #print('T1KEBASP = ', T1KEBASP, 'T2KEBASP = ', T2KEBASP, 'T2/T1 =', T2oT1KEBASP)

#################################### Calculate log(g) #################################### 
### Calculate log(g) to put on our HR diagrams from KEBLAT masses and individual radii ###
### found above just found and propagate relevant errors                               ###
##########################################################################################

    kM1.to(u.g)
    kM1u = kM1 / (u.g)
    kM2.to(u.g)
    kM2u = kM2 / (u.g)
    
    R1u = R1 / (u.cm)
    R2u = R2 / (u.cm)
    
    
    g1 = (kM1u / (R1u**2))
    logg1 = np.log10(g1)
    g1err = np.sqrt(R1u**2 * kM1_err**2 + kM1u/R1u**3 * R1err**2)
    logg1err = 0.434*(g1err/g1)

    g2 = (kM2u / (R2u**2))
    logg2 = np.log10(g2)
    g2err = np.sqrt(R2u**2 * kM2_err**2 + kM2u/R2u**3 * R2err**2)
    logg2err = 0.434*(g2err/g2)

    #print('logg1 = {0:.6f} +/- {1:.6f}, logg2 = {2:.6f} +/- {3:.6f}'.format(logg1, logg1err, logg2, logg2err))
    
#### Next, we will confirm that the ASPCAP temperature is the flux weighted sum off the temperature components:
    FWSTemp = (T1KEBASP * KEBLATFlux1 + T2KEBASP * KEBLATFlux2) / (KEBLATFlux1 + KEBLATFlux2)
    FWSTempASPCAPdiff = ((np.abs(ASPCAPTeff - FWSTemp)) / ((ASPCAPTeff + FWSTemp) / 2) * 100)
    print('ASPCAP temperature  = ', ASPCAPTeff)
    print('Temperature from ASPCAP and KEBLAT calculation =', T1KEBASP)
    print('So, the difference between them is', FWSTempASPCAPdiff, '%')