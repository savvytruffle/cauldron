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
##########################################################################################
isochroneLogTeffs = []
isochroneLogggs = []
col=[]
mrk=[]

for starId in zip(starIds):

    for i in (range(len(starIds))):

        if starId == 5285607:  
            col.append('firebrick')
            mrk.append('o') 
            MVRisofiles = ['isochrones_plot/fehp00afep0_age3.txt']
            labels = ['3 Gyr']
        elif starId == 6864859:
            col.append('chocolate')
            mrk.append('v') 
            MVRisofiles = ['isochrones_plot/fehp00afep0_age4.txt']
            labels = ['4 Gyr']
        elif starId == 6778289:
            col.append('yellow')
            mrk.append('>')
            MVRisofiles = ['isochrones_plot/fehp00afep0_age2.txt']
            labels = ['2 Gyr']
        elif starId == 6449358:
            col.append('darkgreen')
            mrk.append('^')
            MVRisofiles = ['isochrones_plot/fehp00afep0_age5.txt']
            labels = ['5 Gyr']
        elif starId == 4285087:
            col.append('teal') 
            mrk.append('<')
            MVRisofiles = ['isochrones_plot/fehp00afep0_age8.txt']
            labels = ['8 Gyr']
        elif starId ==6131659:
            col.append('dodgerblue')
            mrk.append('o')
            
        elif starId == 6781535:
            col.append('indigo')
            mrk.append('d')

            for MVRisofile in MVRisofiles:        
                isochroneMass, isochroneLogg = np.loadtxt(MVRisofile, unpack=True, usecols = (1, 3))
                #m = np.where(isochroneLogg >= 4.1)### this restricts the isochrone points ###
                                               ### used to those with logg >= 4.1 in   ###
                                               ### effect, this cuts the isochrone down###
                                               ### to only include the main sequence.. ###

                m = np.where(isochroneLogg>= 3.9)  ###           TEST for 5285607          ###
                                               ### In its HR diagram, KIC 5285607 primary#
                                               ### had a log(g) < ~3.7. This test is to###
                                               ### see if the isochrones selected for  ###
                                               ### this target will follow the primary ###
                                               ### to its log(g).     (Positive)       ###

                isochroneMass = isochroneMass[m]
                isochroneLogg = isochroneLogg[m]
                isochroneLogggs.append(isochroneLogg)
                isochroneRadius = (np.sqrt(((6.67e-8)*(isochroneMass*(1.989e33)))/(10**isochroneLogg)))/(6.955e10)

                plt.plot(isochroneMass.value, isochroneRadius.value, ls=':', color = col)
    

        col=[]
        mrk=[]
        if starId == 5285607:  
            col.append('firebrick')
            mrk.append('o') 
            labels = ['3 Gyr']
        elif starId == 6864859:
            col.append('chocolate')
            mrk.append('v') 
            labels = ['4 Gyr']
        elif starId == 6778289:
            col.append('yellow')
            mrk.append('>')
            labels = ['2 Gyr']
        elif starId == 6449358:
            col.append('darkgreen')
            mrk.append('^')
            labels = ['5 Gyr']
        elif starId == 4285087:
            col.append('teal') 
            mrk.append('<')
            labels = ['8 Gyr']
        elif starId ==6131659:
            col.append('dodgerblue')
            mrk.append('o')
        elif starId == 6781535:
            col.append('indigo')
            mrk.append('d')
        
    plt.errorbar(kM1s.value, R1s.value, yerr=R1_errs.value, xerr=kM1_errs.value,
                ls='None', marker=mrk[i], color=col[i]) 

    plt.errorbar(kM2s.value, R2s.value, yerr=R2_errs.value, xerr=kM2_errs.value, 
                ls='None', marker=mrk[i], color=col[i], markersize=6, markeredgewidth=1, 
                markerfacecolor='None')

plt.ylabel('Radius $R_{\odot}$')
plt.xlabel('Mass $M_{\odot}}$')
plt.title('Mass VS Radius')
plt.legend()
plt.show()