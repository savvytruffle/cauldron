import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from numpy import trapz
from sys import argv
from astropy import units as u
from astropy import constants as const
from cycler import cycler

##########################################################################################
##################################### BFArea.py ##########################################
'''Calculate the Area underneath the BF curves of the primary and secondary members of a''' 
'''spectroscopic eclipsing binary. The ratio of these areas is directly proportional to '''
'''the flux ratio of the binary (Bayless & Orosz 2006; Stassun et al. 2007). '''
##########################################################################################
##########################################################################################

starIds =         [5285607, 4285087, 6131659, 6449358, 6778289, 6781535, 6864859] #KEPLER Input Catalog
ASPCAPTeffs =     [6495,    5664,    4845,    6237,    6572,    5749,    6417]   * u.K    #ASPCAP Effective Temperature
ASPCAPTeff_errs = [156,     146,     98,      179,     162,     125,     159]    * u.K    #Error on the above 
kRsums =          [3.679,   2.060,   1.525,   2.80,    2.746,   2.641,   3.110]  * u.Rsun #KEBLAT Radius Sums [R_sun]
kRsum_errs =      [0.033,   0.008,   0.005,   0,       0.0003,  0.031,   0.020]  * u.Rsun #KEBLAT Radius Sum errors
R1s =             [2.001,   1.046,   0.873,   2.108,   1.748,   1.381,   1.655]  * u.Rsun #Calculated from the analysis below
R1_errs =         [0.018,   0.004,   0.003,   0.000,   0.000,   0.016,   0.011]  * u.Rsun #Calculated from the analysis below
R2s =             [1.678,   1.014,   0.652,   0.692,   0.998,   1.260,   1.455]  * u.Rsun #Calculated from the analysis below
R2_errs =         [0.021,   0.006,   0.003,   0,       0.000,   0.021,   0.013]  * u.Rsun #Calculated from the analysis below
kM1s =            [1.557,   1.135,   0.942,   1.93,    1.512,   1.006,   1.411]  * u.Msun #KEBLAT Mass_1
kM1_errs =        [0.038,   0.014,   0.010,   0,       0.022,   0.036,   0.028]  * u.Msun #KEBLAT Mass_1 errors
kM2s =            [1.346,   1.101,   0.703,   0.703,   1.092,   1.037,   1.354]  * u.Msun #KEBLAT Mass_2
kM2_errs =        [0.033,   0.014,   0.008,   0,       0.019,   0.036,   0.028]  * u.Msun #KEBLAT Mass_2 errors
kfluxRatios =     [0.6579,  0.9455,  0.1480,  0.107,   0.1916,  0.7610,  0.7256]          #KEBLAT Flux ratios 
kfluxRatioErrs =  [0.1025,  0.0396,  0.0010,  0,       0.00002, 0.2569,  0.0129]          #KEBLAT Flux ratios errors
kradRatios =      [0.839,   0.969,   0.746,   0.328,   0.5708,  0.913,   0.879]           #KEBLAT Radius Ratios 
kradRatiosErrs =  [0.067,   0.020,   0.003,   0,       0.0003,  0.128,   0.008]           #KEBLAT Radius Ratio errors
GAIAdistances =   [799.744, 617.665, np.nan,  835.143, 1099.75, np.nan,  671.276]* u.pc   #Gaia distances 
GAIAdistance_errs = [13.82, 11.903,  np.nan,  18.413,  26.8496, np.nan,  10.8597]* u.pc   #Error on distances

# Placeholder arrays for values we will calculate in the superloop
T1s = []
T1_errs = []
T2s = []
T2_errs = []
logg1s = []
logg1_errs = []
logg2s = []
logg2_errs = []

def isochroneParse(isoDir, isoFile, loggMin=3.9):
    '''
    Given an isochrone file and minimum logg, return columns (lists) for the
    isochrones mass (Msun), radius (Rsun), temp (K), and logg (cgs).
    
    Input
    -----
    isoDir : `str`
        Directory where the isochrone files live
    isoFile : `str`
        Text file for a single age Dartmouth isochrone
    loggMin : `float`
        Only compute stellar parameters for logg above this cutoff
        Default 3.9 is a good proxy for the main sequence
        
    Returns
    -------
    masses : `list`
        Stellar masses, units of grams
    radii : `list`
        Stellar radii, units of cm
    temps : `list`
        Stellar effective temperatures, units of K
    loggs : `list`
        Stellar log(surface gravities), cgs units (dex)
    '''
    isoPath = os.path.join(isoDir, isoFile)
    isoMass, isoLogTeff, isoLogg = np.loadtxt(isoPath, unpack=True, usecols = (1,2,3))
    m = np.where(isoLogg >= loggMin)
    masses = isoMass[m]*u.Msun
    loggs = u.Dex(isoLogg[m]).cgs
    temps = np.power(10, isoLogTeff[m])*u.K
    radii = np.sqrt(const.G*masses/(np.power(10, loggs.value)*(u.cm/(u.s*u.s))))
    return masses.cgs, radii.cgs, temps, loggs


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

    print('BF Area 1 mean = {0:.3f} +/- {1:.3f}'.format(np.mean(PAreas), P_stderror))
    print('BF Area 2 mean = {0:.3f} +/- {1:.3f}'.format(np.mean(SAreas), S_stderror))
    print('BF Flux Ratio (2/1) = {0:.3f} +/- {1:.3f}'.format(SoP, SoP_stderror))

################################## Calculate R1 and R2 ###################################
### Using the KEBLAT radius sum and the KEBLAT radius ratio calculate R1 and R2 and its###
### error propagation                                                                  ###
##########################################################################################

    R1 = kRsum / (1 + kradRatio)
    R1_err = R1 * np.sqrt((kRsum_err.value / (kRsum.value))**2)
    R2 = kRsum - R1
    R2_err = R2 * np.sqrt((kRsum_err.value/kRsum.value)**2 + (R1_err.value/R1.value)**2)
    R2oR1 = R2/R1
    R2oR1_err = R2oR1 * np.sqrt((R2_err.value/R2.value)**2 + ((R1_err.value/R1.value)**2))

######################### Calculate the sum of the KEBLAT fluxes #########################
### Using the ASPCAP temperature, the sum of the squared radii and the distance to our ###
### targets from GAIA(Bailer-Jones et al 2018)and follow through with error propagation###
##########################################################################################

    KEBLATFluxsum = (const.sigma_sb * ASPCAPTeff**4 * (R1**2 + R2**2)) / (GAIAdistance**2)
    
    # partial derivatives for error propagation
    dFsumdDist = -2 * (R1.to(u.cm).value**2 + R1.to(u.cm).value**2) * const.sigma_sb.cgs.value * ASPCAPTeff.value**4 / GAIAdistance.to(u.cm).value**3
    dFsumdR1 = 2 * R1.to(u.cm).value * const.sigma_sb.cgs.value * ASPCAPTeff.value**4 / GAIAdistance.to(u.cm).value**2
    dFsumdR2 = 2 * R2.to(u.cm).value * const.sigma_sb.cgs.value * ASPCAPTeff.value**4 / GAIAdistance.to(u.cm).value**2
    dFsumdTeff = 4 * (R1.to(u.cm).value**2 + R1.to(u.cm).value**2) * const.sigma_sb.cgs.value * ASPCAPTeff.value**3 / GAIAdistance.to(u.cm).value**2
    
    KEBLATFluxsum_err = np.sqrt( (GAIAdistance_err.to(u.cm).value * dFsumdDist)**2 + 
                                 (R1_err.to(u.cm).value * dFsumdR1)**2 +
                                 (R2_err.to(u.cm).value * dFsumdR2)**2 +
                                 (ASPCAPTeff_err.value * dFsumdTeff)**2 ) * u.erg/u.cm**2/u.s
    KEBLATFlux1 = KEBLATFluxsum / (1 + kfluxRatio)
    KEBLATFlux1_err = KEBLATFlux1.cgs * np.sqrt((KEBLATFluxsum_err/KEBLATFluxsum)**2 + 
                                                (kfluxRatioErr/kfluxRatio)**2)
    KEBLATFlux2 = KEBLATFluxsum - KEBLATFlux1
    KEBLATFlux2_err = KEBLATFlux2.cgs * np.sqrt((KEBLATFluxsum_err/KEBLATFluxsum)**2 +
                                                (KEBLATFlux1_err/KEBLATFlux1)**2)
    F2overF1 = KEBLATFlux2/KEBLATFlux1
    F2overF1_err = F2overF1 * np.sqrt((KEBLATFlux1_err.cgs.value/KEBLATFlux1.cgs.value)**2 + 
                                      (KEBLATFlux2_err.cgs.value/KEBLATFlux2.cgs.value)**2)

    #print(KEBLATFlux1.unit)
    print('R1 = {0:.3f} +/- {1:.3f}'.format(R1, R1_err))
    print('R2 = {0:.3f} +/- {1:.3f}'.format(R2, R2_err))
    print('R2/R1 = {0:.3f} +/- {1:.3f}'.format(R2oR1, R2oR1_err))
    print('KEBLATFlux1 = {0:.3e} +/-  {1:.3e}'.format(KEBLATFlux1.cgs, KEBLATFlux1_err.cgs))
    print('KEBLATFlux2 = {0:.3e} +/- {1:.3e}'.format(KEBLATFlux2.cgs, KEBLATFlux2_err.cgs))
    print('F2/F1 = {0:.5f} +/- {1:.5f}'.format(F2overF1, F2overF1_err))

################################### Calculate T1 and T2 ##################################
###Now we will use the fluxes we just found, and assuming the ASPCAP temperature is the###
###flux weighted sum of the system and that for each system the primary contributes the###
###majority of the light, we can calculate the temperatures of the components.         ###
##########################################################################################

    T1KEBASP = (KEBLATFlux1 * GAIAdistance**2 / (const.sigma_sb * R1**2))**(1/4)

    # partial derivatives for error propagation
    dT1dDist = T1KEBASP.value / (2*GAIAdistance.to(u.cm).value)
    dT1dFlux = T1KEBASP.value / (4*KEBLATFlux1.cgs.value)
    dT1dRad = T1KEBASP.value / (-2*R1.to(u.cm).value)

    T1KEBASP_err = np.sqrt( (GAIAdistance_err.to(u.cm).value * dT1dDist)**2 + 
                            (KEBLATFlux1_err.cgs.value * dT1dFlux)**2 +
                            (R1_err.to(u.cm).value * dT1dRad)**2 )

    T2KEBASP = (KEBLATFlux2 * GAIAdistance**2 / (const.sigma_sb * R2**2))**(1/4)

    # partial derivatives for error propagation
    dT2dDist = T2KEBASP.value / (2*GAIAdistance.to(u.cm).value)
    dT2dFlux = T2KEBASP.value / (4*KEBLATFlux2.cgs.value)
    dT2dRad = T2KEBASP.value / (-2*R2.to(u.cm).value)

    T2KEBASP_err = np.sqrt( (GAIAdistance_err.to(u.cm).value * dT2dDist)**2 + 
                            (KEBLATFlux2_err.cgs.value * dT2dFlux)**2 +
                            (R2_err.to(u.cm).value * dT2dRad)**2 )

    T2oT1KEBASP = (T2KEBASP/T1KEBASP)

    print('T1KEBASP = {0:.0f} +/- {1:.0f}'.format(T1KEBASP, T1KEBASP_err))
    print('T2KEBASP = {0:.0f} +/- {1:.0f}'.format(T2KEBASP, T2KEBASP_err))
    print('T2oT1KEBASP = {0:.3f}'.format(T2oT1KEBASP))
    
    # Actually save T1 and T2 so we can loop through them later
    T1s.append(T1KEBASP)
    T1_errs.append(T1KEBASP_err)
    T2s.append(T2KEBASP)
    T2_errs.append(T2KEBASP_err)

#################################### Calculate log(g) #################################### 
### Calculate log(g) to put on our HR diagrams from KEBLAT masses and individual radii ###
### found above and propagate relevant errors                                          ###
##########################################################################################

    logg1 = u.Dex((const.G * kM1 / R1**2).cgs)
    g1_err = np.sqrt( (kM1_err.to(u.g).value * const.G.cgs.value/(R1.to(u.cm).value)**2)**2 + 
                      (R1_err.to(u.cm).value * (-2*const.G.cgs.value*kM1.to(u.g).value/(R1.to(u.cm).value)**3))**2 )
    logg1_err = 0.434 * g1_err / 10**(logg1.value)

    logg2 = u.Dex((const.G * kM2 / R2**2).cgs)
    g2_err = np.sqrt( (kM2_err.to(u.g).value * const.G.cgs.value/(R2.to(u.cm).value)**2)**2 + 
                      (R2_err.to(u.cm).value * (-2*const.G.cgs.value*kM2.to(u.g).value/(R2.to(u.cm).value)**3))**2 )
    logg2_err = 0.434 * g2_err / 10**(logg2.value)

    print('logg1 = {0:.4f} +/- {1:.4f}'.format(logg1, logg1_err))
    print('logg2 = {0:.4f} +/- {1:.4f}'.format(logg2, logg2_err))
    
    # Actually save logg1 and logg2 so we can loop through them later
    logg1s.append(logg1)
    logg1_errs.append(logg1_err)
    logg2s.append(logg2)
    logg2_errs.append(logg2_err)
    
#### Next, we will confirm that the ASPCAP temperature is the flux weighted sum off the temperature components:
# MR says: welllll this is a unit disaster ...
#    FWSTemp = (T1KEBASP * KEBLATFlux1 + T2KEBASP * KEBLATFlux2) / (KEBLATFlux1 + KEBLATFlux2)
#    FWSTemp_err = np.sqrt(((KEBLATFlux1/(u.W/(u.m**2))/(KEBLATFluxsum/(u.W/(u.m**2))))*T1KEBASP_err)**2+((KEBLATFlux1/(u.W/(u.m**2))*
#    ((T1KEBASP/u.K)-(T2KEBASP/u.K)/(KEBLATFluxsum/(u.W/(u.m**2)))**2))*KEBLATFlux1_err)**2+(((KEBLATFlux1/(u.W/(u.m**2)))/(KEBLATFluxsum/(u.W/(u.m**2))))*T2KEBASP_err)**2
#    +((KEBLATFlux1/(u.W/(u.m**2))*((T2KEBASP/u.K)-(T1KEBASP/u.K)/(KEBLATFluxsum/(u.W/(u.m**2)))**2))*KEBLATFlux1_err)**2)
#    diff = ASPCAPTeff - FWSTemp
#    sum = ASPCAPTeff + FWSTemp
#    FWSTempASPCAPdiff = np.abs(diff) / (sum/2) * 100
    #print('ASPCAP temperature  = {0:.3f} +/- {1:.3}'.format(ASPCAPTeff, ASPCAPTeff_err)
    #print('Temperature from our calculation using ASPCAP and light curve analysis = {0:.3f} +/- {1:.3f}'.format(T1KEBASP, T1KEBASP_err)
    #print('The flux weighted sum of the temperatures we calculated is {0:.3f} +/- {1:.3f}'.format(FWSTemp, FWSTemp_err)
    #print('So, the difference between them is', FWSTempASPCAPdiff, '%')
    print('Analysis Complete for star', starId)
    print('##############################################################')

################################   HR Diagram Party!   ################################### 
### Now, we plot HR diagrams for our targets with evolutionary tracks from theDartmouth###
### Stellar Evolution database (Dotter et al. 2008). A fit finds where the magnitude in###
### the Kepler bandpass (from KEBLAT flux) is closest to the magnitude in the Kepler   ###
### bandpass reported in the Stellar Evolution Database. Note, you must do EITHER the  ###
### HR Diagrams (line 276:makeHRPlots = True) OR the Mass VS. Radius plot (line 327:   ###
### makeMVRPlots = True, but you cannot do both with one execution of the program.     ###
##########################################################################################
    isochroneLogTeffs = []
    isochroneLogggs = []
    
    makeHRPlot = True                      ###       "True" for HR diagrams!         ###
    amagkep1=-2.5*np.log10(kfluxRatio**(-1)) ### We find the apparent magnitude in the ###
    amagkep2=2.5* np.log10(kfluxRatio**(-1)) ### Kepler band through the apparent mag  ###
                                             ### flux relation with KEBLAT flux ratios.###
                                             

################################   Log(Teff) VS Log(g)   #################################
##################### ALL the targets on the same Mass-Radius Plot #######################
########################################################################################## 

plt.figure(figsize=(9,6))
plt.rc('text', usetex=True)
isofiles = ['fehp00afep0_age1.txt', 'fehm0p5afep0_age1.txt',
            'fehp00afep0_age3.txt', 'fehm1afep0_age1.txt']
colors = ['#e41a1c','#377eb8','#4daf4a','#000000','#984ea3','#ff7f00','#a65628']
linestyles = ['-', ':', '-.', '--']
#isolabels = ['1.00 Gyr, [Fe/H] $=0$', '1.50 Gyr, [Fe/H] $=0$', '1.00 Gyr, [Fe/H] $=-0.5$',
#             '0.75 Gyr, [Fe/H] $=-1$'] #, '1.00 Gyr, [Fe/H] $=-1$']
isolabels = [i[0:-9] for i in isofiles]  # for now

for idx, (starId,  T1,  T1_err,  T2,  T2_err,  logg1,  logg1_err,  logg2,  logg2_err) in enumerate(zip(
          starIds, T1s, T1_errs, T2s, T2_errs, logg1s, logg1_errs, logg2s, logg2_errs)):
    
    plt.errorbar(np.log10(T1.value), logg1.value, yerr=logg1_err, 
            xerr=(0.434*(T1_err/T1.value)), marker='o', markersize=8, 
            markeredgewidth=1, ls='None', label=starId, c=colors[idx])

    plt.errorbar(np.log10(T2.value), logg2.value, yerr=logg2_err, 
            xerr=(0.434*(T2_err/T2.value)), ls='None', marker='o', 
            markersize=8, markeredgewidth=1, markerfacecolor='None',
            label='_nolegend_', c=colors[idx])

for idx, isofile in enumerate(isofiles):     
    isoMasses, isoRadii, isoTemps, isoLoggs = isochroneParse('isochrones_plot', isofile)
    plt.plot(np.log10(isoTemps.value), isoLoggs, ls=linestyles[idx],
             lw=2, color='0.75', label=isolabels[idx])

#plt.gca().invert_yaxis()                 ### Y axis increasing downward ###
#plt.gca().invert_xaxis()                 ### X axis increasing to left ###


plt.xlim([3.9, 3.65])                     ### Zooms in plots 
plt.ylim([4.8, 3.9])                     ### Zooms in plots
#    plt.xlim([3.85, 3.78])                   ### Zooms in plots
#    plt.ylim([4.32, 4.0])                    ### Zooms in plots

plt.xlabel('$\log T_{\mathrm{eff}}$', size=16)
plt.ylabel('$\log g$',size=16)
#plt.title('$\log g$ VS $\log T_{\mathrm{eff}}$')
plt.legend(frameon=False, fontsize=12)
plt.show()
    

#     if starId == 5285607: 
#         isofiles = ['isochrones_plot/fehp00afep0_age1.txt', 'isochrones_plot/fehp00afep0_age2.txt', 
#         'isochrones_plot/fehp00afep0_age3.txt', 'isochrones_plot/fehp00afep0_age4.txt', 
#         'isochrones_plot/fehp00afep0_age5.txt']
#         labels = ['1 Gyr', '2 Gyr', '3 Gyr', '4 Gyr', '5 Gyr']
#  
#     elif starId == 6864859:
#         isofiles = ['isochrones_plot/fehp00afep0_age3.txt', 
#         'isochrones_plot/fehp00afep0_age4.txt', 'isochrones_plot/fehp00afep0_age4p5.txt', 
#         'isochrones_plot/fehp00afep0_age5.txt']
#         labels = ['3 Gyr', '4 Gyr', '4.5 Gyr', '5 Gyr']
# 
#     elif starId == 6778289:
#         isofiles = ['isochrones_plot/fehp00afep0_age1p5.txt',
#         'isochrones_plot/fehp00afep0_age2.txt', 'isochrones_plot/fehp00afep0_age3.txt', 
#         'isochrones_plot/fehp00afep0_age4.txt', 'isochrones_plot/fehp00afep0_age5.txt', 
#         'isochrones_plot/fehp00afep0_age6.txt','isochrones_plot/fehp00afep0_age7.txt']
#         labels = ['1.5 Gyr', '2 Gyr', '3 Gyr', '4 Gyr', '5 Gyr', '6 Gyr', '7 Gyr']
# 
#     elif starId == 6449358:
#         isofiles = ['isochrones_plot/fehp00afep0_age1.txt', 'isochrones_plot/fehp00afep0_age2.txt', 
#         'isochrones_plot/fehp00afep0_age3.txt', 'isochrones_plot/fehp00afep0_age5.txt', 
#         'isochrones_plot/fehp00afep0_age8.txt', 'isochrones_plot/fehm05afep0_age1.txt',
#         'isochrones_plot/fehm1afep0_age1.txt']
#         labels = ['1 Gyr', '2 Gyr', '3 Gyr', '5 Gyr', '8 Gyr', '-0.5, 1 Gyr', '-1.0, 1 Gyr']
# 
#     elif starId == 4285087:
#         isofiles = ['isochrones_plot/fehp00afep0_age5.txt', 
#         'isochrones_plot/fehp00afep0_age6.txt', 'isochrones_plot/fehp00afep0_age7.txt',
#         'isochrones_plot/fehp00afep0_age8.txt',
#         'isochrones_plot/fehp00afep0_age10.txt', 'isochrones_plot/fehp00afep0_age12.txt']
#         labels = ['5 Gyr', '6 Gyr', '7 Gyr', '8 Gyr', '10 Gyr', '12 Gyr']
# 
#     elif starId == 6131659:                   ###   No Gaia distances for this target   ###
#         isofiles = ['isochrones_plot/fehp00afep0_age1.txt', 'isochrones_plot/fehp00afep0_age2.txt', 
#         'isochrones_plot/fehp00afep0_age3.txt']
#         labels = ['1 Gyr', '2 Gyr', '3 Gyr']
#         print('No Gaia distances for this target')
# 
# 
#     elif starId == 6781535:                   ### No Gaia distances for this target ###
#         isofiles = ['isochrones_plot/fehp00afep0_age1.txt']
#         labels = ['Fe/H = 0.00']
#         print('No Gaia distances for this target')
# 
# 
#     else: print('No Isochrone file specified')
# 
#         for isofile in isofiles:
#             isochroneMass, isochroneLogTeff, isochroneLogg, kp = np.loadtxt(isofile, 
#             unpack=True, usecols = (1, 2, 3, 13))
#             #m = np.where(isochroneLogg >= 4.1)### this restricts the isochrone points ###
#                                                ### used to those with logg >= 4.1 in   ###
#                                                ### effect, this cuts the isochrone down###
#                                                ### to only include the main sequence.. ###
# 
#             m = np.where(isochroneLogg>= 3.9)###           TEST for 5285607            ###
#                                                ### In its HR diagram, KIC 5285607 primary#
#                                                ### had a log(g) < ~3.7. This test is to###
#                                                ### see if the isochrones selected for  ###
#                                                ### this target will follow the primary ###
#                                                ### to its log(g).     (Positive)       ###
# 
#             isochroneMass = isochroneMass[m]
#             isochroneLogTeff = isochroneLogTeff[m]
#             isochroneLogg = isochroneLogg[m]
#             kp = kp[m]                         ### apparent magnitude in kepler bandpass###
#                                                ### Stellar evolution database.         ###
# 
#             isochroneLogTeffs.append(isochroneLogTeff)
#             isochroneLogggs.append(isochroneLogg)
# 
#         #for logteff, logg, label in zip(isochroneLogTeffs, isochroneLogggs, labels):
#         #            plt.plot(logteff, logg, ls=':', label=label)
#         
#         plt.errorbar(np.log10(T1KEBASP.value), logg1.value, yerr=logg1_err, 
#                          xerr=(0.434*(T1KEBASP_err/T1KEBASP.value)), color='C3', ls='None', marker='o', 
#                          label='Primary')
#         
#         plt.errorbar(np.log10(T2KEBASP.value), logg2.value, yerr=logg2_err, 
#                         xerr=(0.434*(T2KEBASP_err/T2KEBASP.value)), color='C2', ls='None', marker='o', 
#                         markersize=6, markeredgewidth=1, markeredgecolor='C2', markerfacecolor='None', 
#                         label='Secondary')
#         
#         plt.gca().invert_yaxis()                 ### Y axis increasing downward ###
#         plt.gca().invert_xaxis()                 ### X axis increasing to left ###
#         plt.xlabel('$\log T_{\mathrm{eff}}$')
# 
# #        plt.xlim([3.78, 3.72])                  ### Zooms in plots 
# #        plt.ylim([4.6, 4.4])                    ### Zooms in plots
# #        plt.xlim([3.85, 3.78])                   ### Zooms in plots
# #        plt.ylim([4.32, 4.0])                    ### Zooms in plots
#         plt.ylabel('$\log g$')
#         plt.title(starId)
#         plt.legend(frameon=False)
# 
# plt.show()

###############################   m(M_sun) VS r(R_sun)   ################################# 
##########################################################################################

plt.figure(figsize=(6,5))
plt.rc('text', usetex=True)
#isolabels = ['1.00 Gyr, [Fe/H] $=0$', '1.50 Gyr, [Fe/H] $=0$', '1.00 Gyr, [Fe/H] $=-0.5$',
#             '0.75 Gyr, [Fe/H] $=-1$'] #, '1.00 Gyr, [Fe/H] $=-1$']

for idx, (starId,  M1,   R1,  M1_err,   R1_err,  M2,   R2,  M2_err,   R2_err) in enumerate(zip(
          starIds, kM1s, R1s, kM1_errs, R1_errs, kM2s, R2s, kM2_errs, R2_errs)):
    plt.errorbar(M1.value, R1.value, yerr=R1_err.value, xerr=M1_err.value,
                 marker='o', markersize=8, markeredgewidth=1, ls='None',
                 c=colors[idx], label=starId)

    plt.errorbar(M2.value, R2.value, yerr=R2_err.value, xerr=M2_err.value, 
                 ls='None', marker='o', markersize=8, markeredgewidth=1, 
                 markerfacecolor='None', c=colors[idx], label='_nolegend_')
            
#MVRisofiles = ['fehp00afep0_age1.txt', 'fehp00afep0_age1p5.txt',
#               'fehm0p5afep0_age1.txt', 'fehm1afep0_age0p75.txt',
#               'fehm1afep0_age1.txt']
MVRisofiles = isofiles  # intentionally use the same ones as in the logTeff vs logg plot
isolabels = [i[0:-9] for i in isofiles]  # for now

for idx, isofile in enumerate(MVRisofiles):
    isoMasses, isoRadii, isoTemps, isoLoggs = isochroneParse('isochrones_plot', isofile)
    plt.plot(isoMasses.to_value(u.Msun), isoRadii.to_value(u.Rsun), ls=linestyles[idx],
             lw=2, color='0.75', label=isolabels[idx])

#plt.xlim([0.6, 1.6])
#plt.ylim([0.6, 2.1])
plt.ylabel('Radius $R_{\odot}$', size=16)
plt.xlabel('Mass $M_{\odot}$', size=16)
plt.legend(frameon=False, fontsize=12)
#plt.title('Mass VS Radius')

plt.show()

###############################   Compare flux ratios   ################################## 
### Compare the flux ratios from BF and KEBLAT after applying Bolometric corrections   ###
### between the H and Kepler bandpasses (color).                                       ###
##########################################################################################        