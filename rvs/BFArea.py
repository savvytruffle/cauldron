from __future__ import print_function
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
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

    PAreas = (PAmp*PWidth)/(2.35*0.3984)           ### Primary (one value per visit) 
    SAreas = (Samp*SWidth)/(2.35*0.3984)           ### Secondary (one value per visit)
    SoP = np.mean(SAreas)/np.mean(PAreas)          ### secondary over primary ratio of mean BF areas
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

    print('BF Area 1 mean = {0:.3f} +/- {1:.3f}'.format(np.mean(PAreas), P_std))
    print('BF Area 2 mean = {0:.3f} +/- {1:.3f}'.format(np.mean(SAreas), S_std))
    print('BF Flux Ratio (2/1) = {0:.3f} +/- {1:.3f}'.format(SoP, SoP_stderror))

################################## Calculate R1 and R2 ###################################
### Using the KEBLAT radius sum and the KEBLAT radius ratio calculate R1 and R2 and its###
### error propagation                                                                  ###
##########################################################################################

    R1 = kRsum / (1 + kradRatio)
    #R1_err = (R1.value) * np.sqrt((kRsum_err / kRsum.value)**2)
    R1_err = np.sqrt((kRsum_err / kRsum.value)**2)
    R2 = kRsum - R1
    #R2_err = (R2.value) * np.sqrt(((kRsum_err/kRsum.value)**2) + (R1_err/R1.value)**2)
    R2_err = np.sqrt(((kRsum_err/kRsum.value)**2) + (R1_err/R1.value)**2)

    R2oR1 = R2/R1
    #R2oR1_err = R2oR1 * np.sqrt((R2_err/R2.value)**2 + (R1_err/R1.value)**2)

    R2oR1_err = np.sqrt((R2_err/R2.value)**2 + (R1_err/R1.value)**2)

######################### Calculate the sum of the KEBLAT fluxes #########################
### Using the ASPCAP temperature, the sum of the squared radii and the distance to our ###
### targets from GAIA(Bailer-Jones et al 2018)and follow through with error propagation###
##########################################################################################


    KEBLATFluxsum = ((const.sigma_sb * ASPCAPTeff**4 * ((R1.to(u.m))**2 +
                           (R2.to(u.m))**2)) / (GAIAdistance.to(u.m)) **2)
    KEBLATFluxsum_err = np.sqrt( ((ASPCAPTeff_err/ASPCAPTeff.value**2)*R1_err**2)+
                           ((R2_err/R2.value)**2)+(GAIAdistance_err/(GAIAdistance.value)**2))
    KEBLATFlux1 = KEBLATFluxsum / (1 + kfluxRatio)
    #KEBLATFlux1_err = (KEBLATFlux1.value) * np.sqrt((KEBLATFluxsum_err/(KEBLATFluxsum.value))**2 +
    #                    (kfluxRatioErr/kfluxRatio)**2)
    
    KEBLATFlux1_err = np.sqrt((KEBLATFluxsum_err/(KEBLATFluxsum.value))**2 +
                           (kfluxRatioErr/kfluxRatio)**2)
    KEBLATFlux2 = KEBLATFluxsum - KEBLATFlux1
    #KEBLATFlux2_err = (KEBLATFlux2.value) * np.sqrt((KEBLATFluxsum_err/(KEBLATFluxsum.value))**2 +
    #                       (KEBLATFlux1_err/KEBLATFlux1.value)**2)
    KEBLATFlux2_err = np.sqrt((KEBLATFluxsum_err/(KEBLATFluxsum.value))**2 +
                           (KEBLATFlux1_err/KEBLATFlux1.value)**2)
                           
    F2overF1 = KEBLATFlux2/KEBLATFlux1
    
    print('R1 = {0:.3f} +/- {1:.3f}, R2 = {2:.3f} +/- {3:.3f} R2/R1 = {4:.3f}, +/- {5:.3f}'
                            .format(R1, R1_err, R2, R2_err, R2oR1, R2oR1_err))
    print('KEBLATFlux1 = {0:.3e} +/-  {1:.3f} KEBLATFlux2 = {2:.3e} +/- = {3:.3f} F2/F1 = {4:.3f}'
                            .format(KEBLATFlux1, KEBLATFlux1_err, KEBLATFlux2, KEBLATFlux2_err, F2overF1))

################################### Calculate T1 and T2 ##################################
###Now we will use the fluxes we just found, and assuming the ASPCAP temperature is the###
###flux weighted sum of the system and that for each system the primary contributes the###
###majority of the light, we can calculate the temperatures of the components.         ###
##########################################################################################

    T1KEBASP = ((KEBLATFlux1 * (GAIAdistance.to(u.m))**2) / (const.sigma_sb * 
                            (R1.to(u.m))**2))**(1/4)
    T1KEBASP_err = np.sqrt( ((T1KEBASP.value)/4*(KEBLATFlux1.value)*KEBLATFlux1_err)**2 +
                            ((T1KEBASP.value)/2*(GAIAdistance.value)*GAIAdistance_err)**2 +
                            ((-T1KEBASP.value)/2*(R1.value)*R1_err)**2 
                            +(const.sigma_sb.value * R1_err) )

    T2KEBASP = ((KEBLATFlux2 * (GAIAdistance.to(u.m))**2) / (const.sigma_sb * 
                            (R2.to(u.m))**2))**(1/4)
    T2KEBASP_err = np.sqrt( ((T2KEBASP.value)/4*(KEBLATFlux2.value)*KEBLATFlux2_err)**2 +
                            ((T2KEBASP.value)/2*(GAIAdistance.value)*GAIAdistance_err)**2 +
                            ((-T2KEBASP.value)/(2*(R2.value))*R2_err)**2 
                            + (const.sigma_sb.value * R2_err) )

    T2oT1KEBASP = (T2KEBASP/T1KEBASP)

    print('T1KEBASP = {0:.0f} +/- {1:.0f} T2KEBASP = {2:.0f} +/- {3:.0f} T2oT1KEBASP = {4:.3f}'
                            .format(T1KEBASP, T1KEBASP_err, T2KEBASP, T2KEBASP_err, T2oT1KEBASP))

#################################### Calculate log(g) #################################### 
### Calculate log(g) to put on our HR diagrams from KEBLAT masses and individual radii ###
### found above and propagate relevant errors                                          ###
##########################################################################################

    
    logg1 = u.Dex((const.G*kM1 / R1**2).cgs)
    logg1_err = 0.434*(u.Dex((const.G*kM1_err*u.Msun/(R1_err*u.Rsun)**2).cgs)).value/logg1.value

    logg2 = u.Dex((const.G*kM2 / R2**2).cgs)
    logg2_err = 0.434*(u.Dex((const.G*kM2_err*u.Msun / (R2_err*u.Rsun)**2).cgs)).value / logg2.value

    print('logg1 = {0:.4f} +/- {1:.4f}, logg2 = {2:.4f} +/- {3:.4f}'
                            .format(logg1, logg1_err, logg2, logg2_err))
    
#### Next, we will confirm that the ASPCAP temperature is the flux weighted sum off the temperature components:
    FWSTemp = (T1KEBASP * KEBLATFlux1 + T2KEBASP * KEBLATFlux2) / (KEBLATFlux1 + KEBLATFlux2)
    FWSTemp_err = np.sqrt(  ((KEBLATFlux1.value/KEBLATFluxsum.value)*T1KEBASP_err)**2 +
                            (((KEBLATFlux1.value*T1KEBASP.value-T2KEBASP.value) /
                            KEBLATFluxsum.value**2)*KEBLATFlux1_err)**2 +
                            ((KEBLATFlux1.value/KEBLATFluxsum.value)*T2KEBASP_err)**2 +
                            ((KEBLATFlux1.value*T2KEBASP.value-T1KEBASP.value /
                            KEBLATFluxsum.value**2)*KEBLATFlux1_err)**2 ) 
    diff = ASPCAPTeff - FWSTemp
    sum = ASPCAPTeff + FWSTemp
    FWSTempASPCAPdiff = np.abs(diff) / (sum/2) * 100

    #print('ASPCAP temperature  = {0:.3f} +/- {1:.3}'.format(ASPCAPTeff, ASPCAPTeff_err)
    #print('Temperature from our calculation using ASPCAP and light curve analysis = {0:.3f} +/- {1:.3f}'.format(T1KEBASP, T1KEBASP_err)
    #print('The flux weighted sum of the temperatures we calculated is {0:.3f} +/- {1:.3f}'.format(FWSTemp, FWSTemp_err)
    #print('So, the difference between them is', FWSTempASPCAPdiff, '%')

################################   HR Diagram Party!   ################################### 
### Now, we plot HR diagrams for our targets with evolutionary tracks from theDartmouth###
### Stellar Evolution database (Dotter et al. 2008). A fit finds where the magnitude in###
### the Kepler bandpass (from KEBLAT flux) is closest to the magnitude in the Kepler   ###
### bandpass reported in the Stellar Evolution Database. 
##########################################################################################

    makePlots = False                        ###       "True" for HR diagrams!         ###
    
    amagkep1=-2.5*np.log10(kfluxRatio**(-1)) ### We find the apparent magnitude in the ###
    amagkep2=2.5* np.log10(kfluxRatio**(-1)) ### Kepler band through the apparent mag  ###
                                             ### flux relation with KEBLAT flux ratios.###

    if starId == 5285607: 
        isofiles = ['isochrones/fehp00afep4_age1.txt', 'isochrones/fehm00afep8_age1.txt', 
        'isochrones/fehp05afep2_age1.txt', 'isochrones/fehp02afep0_age2p5.txt']
        labels = ['1 Gyr $Z=0.01$', '1 Gyr $Z=0.09$', '1 Gyr $Z=0.57$', '2.5 Gyr $Z=0.21$'] ### labels verified

    elif starId == 6864859:
        isofiles = ['isochrones/fehm00afep8_age1.txt', 
        'isochrones/fehp02afep0_age2p5.txt', 'isochrones/fehp05afep0_age1.txt']
        labels = ['1 Gyr $Z=0.09$', '2.5 Gyr $Z=0.21$', '1 Gyr $Z=0.56$']

    elif starId == 6778289:
        isofiles = ['isochrones/fehm00afem2_age1.txt', 
        'isochrones/fehm07afep0_age1.txt', 'isochrones/fehm07afep0_age2.txt']
        labels = ['1 Gyr $Z=0.6$', '1 Gyr $Z=0.7$', '2 Gyr $Z=0.7$']

    elif starId == 4285087:
        isofiles = ['isochrones/fehm048afep0p8_age1p5.txt', 'isochrones/fehm00afep8_age1.txt', 
        'isochrones/fehm00afep8_age1p25.txt', 'isochrones/fehm00afep8_age2.txt']
        labels = ['1.5 Gyr $Z=-0.48$', '1 Gyr $Z=0.9$',
        '1.25 Gyr $Z=0.9$', '2 Gyr $Z=0.09$']

    elif starId == 6131659:
       isofiles = ['isochrones/fehp02afep0_age15.txt', 'isochrones/fehp05afep0_age1.txt',
       'isochrones/fehp02afep0_age2p5.txt','isochrones/fehm07afep0_age2.txt']
       labels = ['1 Gyr $Z=0.56$','2.5 Gyr $Z=0.21$','2 Gyr $Z=0.7$']

    elif starId == 6781535:
       isofiles = ['isochrones/fehm05afep0_age1.txt', 'isochrones/fehm05afep0_age1p75.txt',
       'isochrones/fehm048afep0p8_age1.txt', 'isochrones/fehm00afem2_age1.txt', 
       'isochrones/fehm07afep0_age1.txt','isochrones/fehm07afep0_age2.txt', 
       'isochrones/fehm00afep8_age3.txt']
       labels = ['1 Gyr $Z=-0.56$', '1.75 Gyr $Z=-0.56$', '1 Gyr $Z=-0.48$', 
       '1 Gyr $Z=0.6','1 Gyr $Z=0.07', '2 Gyr $Z=0.07', '3 Gyr $Z=0.09']

    else: print('No Isochrone file specified')

    if makePlots == True:
        isochroneLogTeffs = []
        isochroneLogggs = []
        for isofile in isofiles:
            isochroneMass, isochroneLogTeff, isochroneLogg, kp = np.loadtxt(isofile, 
            unpack=True, usecols = (1, 2, 3, 13))
            m = np.where(isochroneLogg >= 4.1) ### this restricts the isochrone points ###
                                               ### used to those with logg >= 4.1 in   ###
                                               ### effect, this cuts the isochrone down###
                                               ### to only include the main sequence.. ###

            #m = np.where(isochroneLogg>= 3.85)###           TEST for 5285607          ###
                                               ### In its HR diagram, KIC 5285607 primary#
                                               ### had a log(g) < ~3.7. This test is to###
                                               ### see if the isochrones selected for  ###
                                               ### this target will follow the primary ###
                                               ### to its log(g).     (Positive)       ###

            isochroneMass = isochroneMass[m]
            isochroneLogTeff = isochroneLogTeff[m]
            isochroneLogg = isochroneLogg[m]
            kp = kp[m]                         ### apparent magnitude in kepler bandpass###
                                               ### Stellar evolution database.         ###

            isochroneLogTeffs.append(isochroneLogTeff)
            isochroneLogggs.append(isochroneLogg)
        #sns.set()
        sns.color_palette("colorblind", n_colors=7)        
        for logteff, logg, label in zip(isochroneLogTeffs, isochroneLogggs, labels):
            with sns.color_palette("colorblind", 5):
                plt.plot(logteff, logg, ls=':', label=label)

            #print('logteff =', logteff, 'logg =', logg, 'isochroneLogTeff =', isochroneLogTeff)

        with sns.color_palette("colorblind", 2):     ### Color-blindness color palette ###
            plt.plot(np.log10(T1KEBASP.value), logg1, color='C3', ls='None', marker='o', label='Primary')
            #plt.plot(np.log10(T2KEBASP.value), logg2, color='C2', ls='None', marker='o', label='Secondary')
            plt.plot(np.log10(T2KEBASP.value), logg2, color='C2', ls='None', marker='o', markersize=5,
            markeredgewidth=2, markeredgecolor='g', markerfacecolor='None', label='Secondary')


        #plt.errorbar(np.log10(T1KEBASP.value), logg1, yerr=logg1_err, xerr=T1KEBASP_err, color='C2', ls='None', marker='o', label='Primary')
        #plt.errorbar(np.log10(teff2s[0]), logg2, yerr=logg2_err, xerr=Teff2err, color='C3', ls='None', marker='o', label='Secondary')
        #plt.errorbar(np.log10(teff2s[0]), logg2, yerr=logg2_err, color='C3', ls='None', marker='o', label='Secondary')
        plt.gca().invert_yaxis()               ### Inverts Y axis (increasing downward)###
        plt.gca().invert_xaxis()               ### Inverts X axis (increasing to left) ###
        plt.ylabel('$\log g$')
        plt.xlabel('$\log T_{\mathrm{eff}}$')
        plt.title(starId)
        plt.legend()
        plt.show()

print('##############################################################')
print('Analysis Complete, see above')

###############################   Compare flux ratios   ################################## 
### Compare the flux ratios from BF and KEBLAT after applying Bolometric corrections   ###
### between the H and Kepler bandpasses (color).                                       ###
##########################################################################################