from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.integrate import simps
from numpy import trapz
from sys import argv
from astropy import units as u
from astropy import constants as const
from sklearn.preprocessing import StandardScaler
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

##########################################################################################
##################################### BFArea.py ##########################################
'''Calculate the Area underneath the BF curves of the primary and secondary members of a''' 
'''spectroscopic eclipsing binary. The ratio of these areas is directly proportional to '''
'''the flux ratio of the binary (Bayless & Orosz 2006; Stassun et al. 2007). '''
##########################################################################################
##########################################################################################

starIds = [5285607, 6864859, 6778289, 6449358, 4285087, 6131659, 6781535]#KEPLER Input Catalog
ASPCAPTeffs = [6495, 6417, 6572, 6237, 5664, 4845, 5749] * u.K           #ASPCAP Effective Temperature
ASPCAPTeff_errs = [156, 159, 162, 179, 146, 98, 125] * u.K               #Error on the above 
kRsums =     [3.489, 3.104, 2.745, 2.80, 2.033, 1.5251, 2.0408] * u.Rsun #KEBLAT Radius Sums [R_sun]
R1s = [2.250, 1.444, 1.747, 2.11, 1.033, 0.908, 1.095] * u.Rsun          #Calculated from the analysis below
R1_errs = [0.033, 0.007,  0.005, 0, 0.004, 0.003, 0.011] * u.Rsun        #Calculated from the analysis below
R2s = [1.239, 1.660, 0.998, 0.692, 1.000, 0.617, 0.946] * u.Rsun         #Calculated from the analysis below
R2_errs = [0.026, 0.012, 0.004, 0, 0.006, 0.003, 0.014] * u.Rsun         #Calculated from the analysis below
kRsum_errs = [0.051, 0.016, 0.0086, 0, 0.0080, 0.0052, 0.0213] * u.Rsun  #KEBLAT Radius Sum errors
kM1s =     [1.554, 1.354, 1.510, 1.93, 1.137, 0.9422, 1.0057] * u.Msun   #KEBLAT Mass_1
kM1_errs = [0.023, 0.029, 0.022, 0, 0.013, 0.0093, 0.0327] * u.Msun      #KEBLAT Mass_1 errors
kM2s =     [1.333, 1.411, 1.091, 0.783, 1.103, 0.7028, 1.0346] * u.Msun  #KEBLAT Mass_2
kM2_errs = [0.020, 0.028, 0.018, 0, 0.014, 0.0078, 0.0330] * u.Msun      #KEBLAT Mass_2 errors
kfluxRatios = [0.258, 1.407, 0.19138, 0.107, 0.901, 0.1483, 0.9201]      #KEBLAT Flux ratios 
kfluxRatioErrs = [0.046, 0.101, 2.6e-5, 0, 0.080, 0.0017, 0.0524]        #KEBLAT Flux ratios errors 
kradRatios = [0.551, 1.149, 0.57093, 0.328, 0.969, 0.6799, 0.8641]       #KEBLAT Radius Ratios 
kradRatiosErrs = [0.048, 0.020, 0.013, 0, 0.0080, 0.0057, 0.0275]        #KEBLAT Radius Ratio errors
GAIAdistances =  [799.7441, 671.2761, 1099.7471, 835.1428, 617.665, np.nan, np.nan] * u.pc   #Gaia distances 
GAIAdistance_errs = [13.8152, 10.8597, 26.8496, 18.4130, 11.9031, np.nan, np.nan] * u.pc     #Error on distances

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
    R1_err = R1 * np.sqrt((kRsum_err.value / (kRsum.value))**2)
    R2 = kRsum - R1
    R2_err = R2 * np.sqrt((kRsum_err.value/kRsum.value)**2 + (R1_err.value/R1.value)**2)
    R2oR1 = R2/R1
    R2oR1_err = R2oR1 * np.sqrt((R2_err.value/R2.value)**2 + ((R1_err.value/R1.value)**2))

######################### Calculate the sum of the KEBLAT fluxes #########################
### Using the ASPCAP temperature, the sum of the squared radii and the distance to our ###
### targets from GAIA(Bailer-Jones et al 2018)and follow through with error propagation###
##########################################################################################

    # MR is not entirely sure which of these is right ...
    #KEBLATFluxsum = (const.sigma_sb * ASPCAPTeff**4 * kRsum**2) / (GAIAdistance**2)
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

    #print(KEBLATFlux1.unit)  # MR confirmed this is W/m^2 in SI and erg/cm^2/s in cgs, i.e., g/s^3
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

################################   HR Diagram Party!   ################################### 
### Now, we plot HR diagrams for our targets with evolutionary tracks from theDartmouth###
### Stellar Evolution database (Dotter et al. 2008). A fit finds where the magnitude in###
### the Kepler bandpass (from KEBLAT flux) is closest to the magnitude in the Kepler   ###
### bandpass reported in the Stellar Evolution Database. 
##########################################################################################
    isochroneLogTeffs = []
    isochroneLogggs = []
    fig, ax = plt.subplots()
    col=[]
    mrk=[]

    makeHRPlots = True                       ###       "True" for HR diagrams!         ###
    amagkep1=-2.5*np.log10(kfluxRatio**(-1)) ### We find the apparent magnitude in the ###
    amagkep2=2.5* np.log10(kfluxRatio**(-1)) ### Kepler band through the apparent mag  ###
                                             ### flux relation with KEBLAT flux ratios.###

for i in range(len(starIds)):
    if starId == 5285607: 
        isofiles = ['isochrones_plot/fehp00afep0_age1.txt', 'isochrones_plot/fehp00afep0_age2.txt', 
        'isochrones_plot/fehp00afep0_age3.txt', 'isochrones_plot/fehp00afep0_age4.txt', 
        'isochrones_plot/fehp00afep0_age5.txt']
        labels = ['1 Gyr', '2 Gyr', '3 Gyr', '4 Gyr', '5 Gyr']
        x1, x2, y1, y2 = 3.86, 3.77, 3.8, 4.4
        axins = zoomed_inset_axes(ax, 1, loc=9)
        MVRisofiles = ['isochrones_plot/fehp00afep0_age3.txt']
        MVRlabels = ['3 Gyr'] 
        col.append('firebrick')
        mrk.append('o') 

    elif starId == 6864859:
        isofiles = ['isochrones_plot/fehp00afep0_age3.txt', 
        'isochrones_plot/fehp00afep0_age4.txt', 'isochrones_plot/fehp00afep0_age4p5.txt', 
        'isochrones_plot/fehp00afep0_age5.txt']
        labels = ['3 Gyr', '4 Gyr', '4.5 Gyr', '5 Gyr']
        x1, x2, y1, y2 = 3.83, 3.77, 4.1, 4.3
        axins = zoomed_inset_axes(ax, 1, loc=9)
        MVRisofiles = ['isochrones_plot/fehp00afep0_age4.txt']
        MVRlabels = ['4 Gyr']
        col.append('chocolate')
        mrk.append('v') 


    elif starId == 6778289:
        isofiles = ['isochrones_plot/fehp00afep0_age1p5.txt',
        'isochrones_plot/fehp00afep0_age2.txt', 'isochrones_plot/fehp00afep0_age3.txt', 
        'isochrones_plot/fehp00afep0_age4.txt', 'isochrones_plot/fehp00afep0_age5.txt', 
        'isochrones_plot/fehp00afep0_age6.txt','isochrones_plot/fehp00afep0_age7.txt']
        labels = ['1.5 Gyr', '2 Gyr', '3 Gyr', '4 Gyr', '5 Gyr', '6 Gyr', '7 Gyr']
        x1, x2, y1, y2 = 3.85, 3.75, 4.0, 4.6
        axins = zoomed_inset_axes(ax, 1, loc=9)
        MVRisofiles = ['isochrones_plot/fehp00afep0_age2.txt']
        MVRlabels = ['2 Gyr']
        col.append('yellow')
        mrk.append('>')


    elif starId == 6449358:
        isofiles = ['isochrones_plot/fehp00afep0_age1.txt', 'isochrones_plot/fehp00afep0_age2.txt', 
        'isochrones_plot/fehp00afep0_age3.txt', 'isochrones_plot/fehp00afep0_age5.txt', 
        'isochrones_plot/fehp00afep0_age8.txt']
        labels = ['1 Gyr', '2 Gyr', '3 Gyr', '5 Gyr', '8 Gyr']
        x1, x2, y1, y2 = 3.83, 3.75, 4.0, 4.8
        axins = zoomed_inset_axes(ax, 1, loc=9)
        MVRisofiles = ['isochrones_plot/fehp00afep0_age5.txt']
        MVRlabels = ['5 Gyr']
        col.append('darkgreen')
        mrk.append('^')

    elif starId == 4285087:
        isofiles = isofiles = ['isochrones_plot/fehp00afep0_age5.txt', 
        'isochrones_plot/fehp00afep0_age6.txt', 'isochrones_plot/fehp00afep0_age7.txt',
        'isochrones_plot/fehp00afep0_age8.txt',
        'isochrones_plot/fehp00afep0_age10.txt', 'isochrones_plot/fehp00afep0_age12.txt']
        labels = ['5 Gyr', '6 Gyr', '7 Gyr', '8 Gyr', '10 Gyr', '12 Gyr']
        x1, x2, y1, y2 = 3.77, 3.74, 4.4, 4.6
        axins = zoomed_inset_axes(ax, 3, loc=9)
        MVRisofiles = ['isochrones_plot/fehp00afep0_age8.txt']
        MVRlabels = ['8 Gyr']
        col.append('teal') 
        mrk.append('<')

    elif starId == 6131659:                  ###   No Gaia distances for this target   ###
        isofiles = ['isochrones_plot/fehp00afep0_age1.txt', 'isochrones_plot/fehp00afep0_age2.txt', 
        'isochrones_plot/fehp00afep0_age3.txt']
        labels = ['1 Gyr', '2 Gyr', '3 Gyr']
#        x1, x2, y1, y2 = 5, 5, 5, 5
        MVRisofiles = ['isochrones_plot/fehp00afep0_age8.txt']
        col.append('dodgerblue')
        mrk.append('d')
        print('No Gaia distances for this target')


    elif starId == 6781535:                  ###   No Gaia distances for this target   ###
        isofiles = ['isochrones_plot/fehp00afep0_age1.txt']
        labels = ['Fe/H = 0.00']
#        x1, x2, y1, y2 = 5, 5, 5, 5
        MVRisofiles = ['isochrones_plot/fehp00afep0_age8.txt']
        col.append('indigo')
        mrk.append('d')
        print('No Gaia distances for this target')


    else: print('No Isochrone file specified')

    if makeHRPlots == True:

        for isofile in isofiles:
            isochroneMass, isochroneLogTeff, isochroneLogg, kp = np.loadtxt(isofile, 
            unpack=True, usecols = (1, 2, 3, 13))
            #m = np.where(isochroneLogg >= 4.1) ### this restricts the isochrone points ###
                                               ### used to those with logg >= 4.1 in   ###
                                               ### effect, this cuts the isochrone down###
                                               ### to only include the main sequence.. ###

            m = np.where(isochroneLogg>= 3.9)###           TEST for 5285607          ###
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
            isochroneRadius = (np.sqrt(((6.67e-8)*(isochroneMass*(1.989e33)))/(10**isochroneLogg)))/(6.955e10)

 ###############################   Log(Teff) VS Log(g)   ################################# 
        for logteff, logg, label in zip(isochroneLogTeffs, isochroneLogggs, labels):
            ax.plot(logteff, logg, ls=':')
            axins.plot(logteff, logg, ls=':')
#        axins.invert_yaxis()               ### Inverts Y axis (increasing downward)###
#        axins.invert_xaxis()               ### Inverts X axis (increasing to left) ###
        
            ax.errorbar(np.log10(T1KEBASP.value), logg1.value, yerr=logg1_err, 
                 xerr=(0.434*(T1KEBASP_err/T1KEBASP.value)), color='C3', ls='None', marker='o', 
                 label='Primary')

            ax.errorbar(np.log10(T2KEBASP.value), logg2.value, yerr=logg2_err, 
                xerr=(0.434*(T1KEBASP_err/T2KEBASP.value)), color='C2', ls='None', marker='o', 
                markersize=6, markeredgewidth=1, markeredgecolor='C2', markerfacecolor='None', 
                label='Secondary')

            axins.errorbar(np.log10(T1KEBASP.value), logg1.value, yerr=logg1_err, 
                 xerr=(0.434*(T1KEBASP_err/T1KEBASP.value)), color='C3', ls='None', marker='o')
    
            axins.errorbar(np.log10(T2KEBASP.value), logg2.value, yerr=logg2_err, 
                xerr=(0.434*(T1KEBASP_err/T2KEBASP.value)), color='C2', ls='None', marker='o', 
                markersize=6, markeredgewidth=1, markeredgecolor='C2', markerfacecolor='None')

            ax.invert_yaxis()               ### Inverts Y axis (increasing downward)###
            ax.invert_xaxis()               ### Inverts X axis (increasing to left) ###

            axins.set_xlim(x1, x2)
            axins.set_ylim(y1, y2)
        ax.set_xlabel('$\log T_{\mathrm{eff}}$')
        ax.set_ylabel('$\log g$')
        ax.set_title(starId)
        plt.yticks(visible=False)
        plt.xticks(visible=False)
        mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")
        plt.legend()
        plt.show()

###############################   m(M_sun) VS r(R_sun)   ################################# 

for i in range(len(starIds)):
    if starId == 5285607:  
        col.append('firebrick')
        mrk.append('o') 
        MVRisofiles = ['isochrones_plot/fehp00afep0_age3.txt']
        MVRlabels = ['3 Gyr'] 
    elif starId == 6864859:
        col.append('chocolate')
        mrk.append('v') 
        MVRisofiles = ['isochrones_plot/fehp00afep0_age4.txt']
        MVRlabels = ['4 Gyr']
    elif starId == 6778289:
        col.append('yellow')
        mrk.append('>')
        MVRisofiles = ['isochrones_plot/fehp00afep0_age2.txt']
        MVRlabels = ['2 Gyr']
    elif starId == 6449358:
        col.append('darkgreen')
        mrk.append('^')
        MVRisofiles = ['isochrones_plot/fehp00afep0_age5.txt']
        MVRlabels = ['5 Gyr']
    elif starId == 4285087:
        col.append('teal') 
        mrk.append('<')
        MVRisofiles = ['isochrones_plot/fehp00afep0_age8.txt']
        MVRlabels = ['8 Gyr']
    elif starId ==6131659:
        col.append('dodgerblue')
        mrk.append('d')
    elif starId == 6781535:
        col.append('indigo')
        mrk.append('d')

    plt.errorbar(kM1s.value[i], R1s.value[i], yerr=R1_errs.value[i], xerr=kM1_errs.value[i],
                    ls='None') 
        
    plt.errorbar(kM2s.value[i], R2s.value[i], yerr=R2_errs.value[i], xerr=kM2_errs.value[i], 
                    ls='None', markersize=6, markeredgewidth=1, 
                    markerfacecolor='None')
#
    for MVRisofile in MVRisofiles:
            isochroneMass, isochroneLogTeff, isochroneLogg, kp = np.loadtxt(MVRisofile, 
            unpack=True, usecols = (1, 2, 3, 13))
                #m = np.where(isochroneLogg >= 4.1) ### this restricts the isochrone points ###
                                               ### used to those with logg >= 4.1 in   ###
                                               ### effect, this cuts the isochrone down###
                                               ### to only include the main sequence.. ###

            m = np.where(isochroneLogg>= 3.9)###           TEST for 5285607          ###
                                               ### In its HR diagram, KIC 5285607 primary#
                                               ### had a log(g) < ~3.7. This test is to###
                                               ### see if the isochrones selected for  ###
                                               ### this target will follow the primary ###
                                               ### to its log(g).     (Positive)       ###

    isochroneMass = isochroneMass[m]
    isochroneLogTeff = isochroneLogTeff[m]
    isochroneLogg = isochroneLogg[m]
                                               ### Stellar evolution database.         ###

    isochroneLogTeffs.append(isochroneLogTeff)
    isochroneLogggs.append(isochroneLogg)
    isochroneRadius = (np.sqrt(((6.67e-8)*(isochroneMass*(1.989e33)))/(10**isochroneLogg)))/(6.955e10)

 

    plt.plot(isochroneMass[i], isochroneRadius[i], ls=':', color = col[i])

    plt.ylabel('Radius $R_{\odot}$')
    plt.xlabel('Mass $M_{\odot}}$')
    plt.title('Mass VS Radius')
    plt.legend()
plt.show()



print('##############################################################')
print('Analysis Complete, see above')

###############################   Compare flux ratios   ################################## 
### Compare the flux ratios from BF and KEBLAT after applying Bolometric corrections   ###
### between the H and Kepler bandpasses (color).                                       ###
##########################################################################################