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
6845
starIds =         [5285607, 6864859, 6778289, 6449358, 4285087, 6131659,  6781535] #KEPLER Input Catalog
#ASPCAPTeffs =     [6495,    6417,  6572,     6237,    5804,    4845,     5749]    * u.K    #ASPCAP Effective Temperature
ASPCAPTeff_errs = [156,     159,     162,     179,     1,       98,       125]     * u.K    #Error on the above 
ASPCAPTeffs =     [6845,    6497,    6822,    6737,    5689,    5195,     5849]    * u.K    #ASCAP Effective Temperature with El-Badry et al correction
kRsums =          [3.679,   3.110,   2.746,   2.8231,  2.060,   1.525,    2.641]   * u.Rsun #KEBLAT Radius Sums
kRsum_errs =      [0.033,   0.020,   0.013,   0.0010,  0.008,   0.005,    0.031]   * u.Rsun #KEBLAT Radius Sum errors
R1s =             [2.003,   1.655,   1.748,   2.1254,  1.033,   0.908,    1.382]   * u.Rsun 
R1_errs =         [0.062,   0.013,   0.009,   0.0007,  0.012,   0.003,    0.090]   * u.Rsun 
R2s =             [1.679,   1.455,   0.998,   0.6977,  1.026,   0.616,    1.262]   * u.Rsun 
R2_errs =         [0.087,   0.012,   0.005,   0.00005, 0.011,   0.003,    0.092]   * u.Rsun 
kM1s =            [1.557,   1.411,   1.512,   np.nan,  1.135,   0.942,    1.006]   * u.Msun #KEBLAT Mass_1
kM1_errs =        [0.038,   0.028,   0.022,   np.nan,  0.014,   0.010,    0.036]   * u.Msun #KEBLAT Mass_1 errors
kM2s =            [1.346,   1.354,   1.092,   np.nan,  1.101,   0.703,    1.037]   * u.Msun #KEBLAT Mass_2
kM2_errs =        [0.033,   0.028,   0.019,   np.nan,  0.014,   0.008,    0.036]   * u.Msun #KEBLAT Mass_2 errors
kfluxRatios =     [0.6579,  0.7256,  0.19155, 0.10752, 0.9455,  0.1480,   0.7610]           #KEBLAT Flux ratios 
kfluxRatioErrs =  [0.1025,  0.0139,  0.00002, 6.5e-6,  0.0396,  0.0010,   0.2569]           #KEBLAT Flux ratios errors
kradRatios =      [0.839,   0.879,   0.5708,  0.3283,  0.969,   0.679,    0.913]            #KEBLAT Radius Ratios 
kradRatiosErrs =  [0.067,   0.008,   0.0003,  0.00002, 0.020,   0.003,    0.128]            #KEBLAT Radius Ratio errors
GAIAdistances =   [799.744, 671.276, 1099.75, 835.143, 617.665, np.nan,   np.nan]  * u.pc   #Gaia distances 
GAIAdistance_errs = [13.82, 10.8597, 26.8496, 18.413,  11.903,  np.nan,   np.nan]  * u.pc   #Error on distances

# Placeholder arrays for values we will calculate in the superloop
T1s = []
T1_errs = []
T2s = []
T2_errs = []
logg1s = []
logg1_errs = []
logg2s = []
logg2_errs = []

# Naively try adjusting ASPCAP Teffs to see what happens
#ASPCAPTeffs = ASPCAPTeffs + 500*u.K

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
    
    # Naively try adjusting loggs to see what happens
    #logg1 = logg1.cgs + u.Dex(0.2)
    #logg2 = logg2.cgs + u.Dex(0.2)
    
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

### GLOBAL PROPERTIES FOR ALL THE PLOTS ###

plt.rc('text', usetex=True)

isofiles = ['fehp00afep0_age0p8.txt',
            'fehp00afep0_age1.txt', 'fehm05afep0_age1.txt', 'fehm1afep0_age1.txt',
            'fehp00afep0_age3.txt',
            'fehp00afep0_age5.txt', 'fehm05afep0_age5.txt']

isolabels = ['0.8 Gyr, [Fe/H] $=0$',
             '1 Gyr, [Fe/H] $=0$', '1 Gyr, [Fe/H] $=-0.5$', '1 Gyr, [Fe/H] $=-1$', 
             '3 Gyr, [Fe/H] $=-1$',
             '5 Gyr, [Fe/H] $=0$', '5 Gyr, [Fe/H] $=-0.5$']
 
isolines = ['-',
            '-', '--', ':',
            '-',
            '-', '--']

isocolors = ['#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525']
isocolors = list(reversed(isocolors))  # darker colors first

starcolors = ['#e41a1c', '#377eb8', '#4daf4a', '#a65628', '#984ea3', '#737373', '#ff7f00']

assert len(isofiles) == len(isolabels)
assert len(isofiles) == len(isolines)
assert len(isofiles) == len(isocolors)
assert len(starcolors) == len(starIds)

################################   Log(Teff) VS Log(g)   #################################
##################### ALL the targets on the same Mass-Radius Plot #######################
########################################################################################## 

plt.figure(figsize=(9,6))

for idx, (starId,  T1,  T1_err,  T2,  T2_err,  logg1,  logg1_err,  logg2,  logg2_err) in enumerate(zip(
          starIds, T1s, T1_errs, T2s, T2_errs, logg1s, logg1_errs, logg2s, logg2_errs)):
    
    if starId != 6131659 and starId != 6781535 and starId != 6449358:
        plt.errorbar(np.log10(T1.value), logg1.value, yerr=logg1_err, 
                xerr=(0.434*(T1_err/T1.value)), ls='None', marker='o', 
                markersize=8, markeredgewidth=1, label=starId,
                c=starcolors[idx], zorder=2)

        plt.errorbar(np.log10(T2.value), logg2.value, yerr=logg2_err, 
                xerr=(0.434*(T2_err/T2.value)), ls='None', marker='o', 
                markersize=8, markeredgewidth=1, markerfacecolor='None',
                label='_nolegend_', c=starcolors[idx], zorder=2)
        
        #plt.axvline(x=np.log10(ASPCAPTeffs[idx].value + 200), c=starcolors[idx])

for idx, isofile in enumerate(isofiles):     
    isoMasses, isoRadii, isoTemps, isoLoggs = isochroneParse('isochrones_plot', isofile)
    
    plt.plot(np.log10(isoTemps.value), isoLoggs, ls=isolines[idx],
             lw=2, color=isocolors[idx], label=isolabels[idx], zorder=1)

plt.legend(frameon=False, fontsize=12, loc=1)

plt.xlim([3.9, 3.65])
plt.ylim([4.8, 3.9])

plt.xlabel('$\log T_{\mathrm{eff}}$', size=16)
plt.ylabel('$\log g$',size=16)

plt.show()
plt.savefig('LTeffvLg.eps')

################################   Log(Teff) VS Log(g)   #################################
################################## ONE target per panel ##################################
########################################################################################## 

fig = plt.figure(figsize=(9,6))

for idx, (starId,  T1,  T1_err,  T2,  T2_err,  logg1,  logg1_err,  logg2,  logg2_err) in enumerate(zip(
          starIds, T1s, T1_errs, T2s, T2_errs, logg1s, logg1_errs, logg2s, logg2_errs)):

    # If you only want to plot a subset of isochrones for any star, customize that here
    # Otherwise, the default values assigned above will be plotted instead
    if starId == 5285607:  
        pltidx = 1
    elif starId == 6864859:
        pltidx = 2
    elif starId == 6778289:
        pltidx = 3
    elif starId == 4285087:
        pltidx = 4
    else:
        pltidx = 0
    
    if pltidx > 0:
        # Make a new subplot panel for each star
        ax = fig.add_subplot(2, 2, pltidx)
        plt.subplots_adjust(wspace=0, hspace=0)
        
        for isoidx, isofile in enumerate(isofiles):   
            isoMasses, isoRadii, isoTemps, isoLoggs = isochroneParse('isochrones_plot', isofile)
            
            # Plot all the desired isochrones for each star
            plt.plot(np.log10(isoTemps.value), isoLoggs, ls=isolines[isoidx], lw=2,
                     color=isocolors[isoidx], label=isolabels[isoidx], zorder=1)
    
        # Plot the (Teff, logg) for each star
        plt.errorbar(np.log10(T1.value), logg1.value, yerr=logg1_err, 
                         xerr=(0.434*(T1_err/T1.value)), c=starcolors[idx],
                         ls='None', marker='o', zorder=2)
        plt.errorbar(np.log10(T2.value), logg2.value, yerr=logg2_err, 
                         xerr=(0.434*(T2_err/T2.value)), ls='None', marker='o', 
                         markersize=8, markeredgewidth=1, markerfacecolor='None',
                         c=starcolors[idx], zorder=2)
        
        # Annotate and format plot panels
        ax.invert_yaxis()
        ax.invert_xaxis()

        plt.xlim([3.89, 3.73]) 
        plt.ylim([4.65, 3.95])
        
        if starId != 5285607 and starId != 6778289:
            plt.gca().set_yticklabels([])
        else:
            plt.ylabel('$\log g$', size=14)
        if starId != 4285087 and starId != 6778289:
            plt.gca().set_xticklabels([])            
        else:
            plt.xlabel('$\log T_{\mathrm{eff}}$', size=14)
        
        #plt.legend(bbox_to_anchor=(1.05, 0.8), loc=2, borderaxespad=0., frameon=False, fontsize=12)

        plt.text(3.88, 4.55, starId, size=14)

plt.show()
plt.savefig('Teffloggsubs.eps')

###############################   m(M_sun) VS r(R_sun)   #################################
##################### ALL the targets on the same Mass-Radius Plot #######################
##########################################################################################  

plt.figure(figsize=(9,6))

for idx, (starId,  M1,   R1,  M1_err,   R1_err,  M2,   R2,  M2_err,   R2_err) in enumerate(zip(
          starIds, kM1s, R1s, kM1_errs, R1_errs, kM2s, R2s, kM2_errs, R2_errs)):

    if starId != 6449358:
        plt.errorbar(M1.value, R1.value, yerr=R1_err.value, xerr=M1_err.value,
                     marker='o', markersize=8, markeredgewidth=2, ls='None',
                     c=starcolors[idx], label=starId)

        plt.errorbar(M2.value, R2.value, yerr=R2_err.value, xerr=M2_err.value, 
                     ls='None', marker='o', markersize=8, markeredgewidth=2, 
                     markerfacecolor='None', c=starcolors[idx], label='_nolegend_')

isofiles = ['fehp00afep0_age0p8.txt',
            'fehp00afep0_age1.txt', 'fehm05afep0_age1.txt']

isolabels = ['0.8 Gyr, [Fe/H] $=0$',
             '1 Gyr, [Fe/H] $=0$', '1 Gyr, [Fe/H] $=-0.5$', '1 Gyr, [Fe/H] $=-1$', 
             '3 Gyr, [Fe/H] $=-1$',
             '5 Gyr, [Fe/H] $=0$', '5 Gyr, [Fe/H] $=-0.5$']

for idx, isofile in enumerate(isofiles):
    
    isoMasses, isoRadii, isoTemps, isoLoggs = isochroneParse('isochrones_plot', isofile)
    
    plt.plot(isoMasses.to_value(u.Msun), isoRadii.to_value(u.Rsun), ls=isolines[idx],
             lw=2, c=isocolors[idx], label=isolabels[idx])

plt.xlim([0.2, 1.7])
plt.ylim([0.5, 2.1])
plt.ylabel('Radius ($R_{\odot}$)', size=16)
plt.xlabel('Mass ($M_{\odot}$)', size=16)
plt.legend(frameon=False, fontsize=12, loc=2)

plt.show()

###### A Mass vs Teff plot? Because, why not ######
newfig = plt.figure(figsize=(7,6))

for idx, (starId,  M1,   T1,  M1_err,   T1_err,  M2,   T2,  M2_err,   T2_err) in enumerate(zip(
          starIds, kM1s, T1s, kM1_errs, T1_errs, kM2s, T2s, kM2_errs, T2_errs)):

    if starId != 6131659 and starId != 6781535 and starId != 6449358:
        plt.errorbar(np.log10(T1.value), M1.value, xerr=(0.434*(T1_err/T1.value)), yerr=M1_err.value,
                     marker='o', markersize=8, markeredgewidth=2, ls='None',
                     c=starcolors[idx], label=starId, zorder=2)

        plt.errorbar(np.log10(T2.value), M2.value, xerr=(0.434*(T2_err/T2.value)), yerr=M2_err.value, 
                     ls='None', marker='o', markersize=8, markeredgewidth=2, 
                     markerfacecolor='None', c=starcolors[idx], label='_nolegend_', zorder=2)
    
for idx, isofile in enumerate(isofiles):
    
    isoMasses, isoRadii, isoTemps, isoLoggs = isochroneParse('isochrones_plot', isofile)
    
    plt.plot(np.log10(isoTemps.value), isoMasses.to_value(u.Msun), ls=isolines[idx],
             lw=2, c=isocolors[idx], label=isolabels[idx], zorder=1)

plt.xlabel('$\log T_{\mathrm{eff}}$', size=14)
plt.ylabel('Mass ($M_{\odot}$)', size=16)
plt.legend(frameon=False, fontsize=12)
plt.xlim([3.9, 3.6])
plt.ylim([0.5, 2.2])

plt.show()
plt.savefig('MvR.eps')
