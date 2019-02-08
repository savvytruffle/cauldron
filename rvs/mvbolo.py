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
##################################### mvbolo.py ##########################################
'''Calculate the bolometric luminosity for a subset of Kepler/APOGEE/GAIA SEB from the '''
'''temperature estimates completed with BF_area.py and mass, radius estimates from LC  '''
'''analysis. Then plot this against masses derived from LC analysis. '''
##########################################################################################
##########################################################################################

starIds =         [5285607, 4285087, 6131659, 6449358, 6778289, 6781535, 6864859] #KEPLER Input Catalog
BFLCTeff1s =       [6892,    5684,    np.nan,  6737,    7007,    np.nan,  6541]       #Combined CAULDRON analysis
BFLCTeff1_errs =   [341,     100,     np.nan,  358,     296,     np.nan,  202]        #Combined CAULDRON analysis
BFLCTeff2s =       [6776,    5694,    np.nan,  6733,    6135,    np.nan,  6439]       #Combined CAULDRON analysis
BFLCTeff2_errs =   [389,     116,     np.nan,  500,     358,     np.nan,  274]        #Combined CAULDRON analysis
R1s =             [2.003,   1.033,   0.908,   2.1254,  1.748,   1.382,   1.655]  #KEBLAT R1
R1_errs =         [0.062,   0.012,   0.003,   0.0007,  0.009,   0.090,   0.013]  #KEBLAT R1_err
R2s =             [1.679,   1.026,   0.616,   0.6977,  0.998,   1.262,   1.455]  #KEBLAT R2
R2_errs =         [0.087,   0.011,   0.003,   0.00005, 0.005,   0.092,   0.012]  #KEBLAT R2_err

bololum1 = []
bololum2 = []


#### For the isochrone tracks that we overlay ####
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

for starId, BFLCTeff1, BFLCTeff1_err, BFLCTeff2, BFLCTeff2_err, R1, R1_err, R2, R2_err in zip(
    starIds, BFLCTeff1s, BFLCTeff1_errs, BFLCTeff2s, BFLCTeff2_errs, R1s, R1_errs, R2s, R2_errs):

    print(' ')
    print('##############################################################')
    print('Running analysis for star', starId)
    
    bololum1 = [4*np.pi*(R1*const.R_sun)**2*const.sigma_sb*BFLCTeff1**4] / const.L_bol0
    bololum2 = [4*np.pi*(R2*const.R_sun)**2*const.sigma_sb*BFLCTeff2**4] / const.L_bol0
    
####Error propagation in bolometric luminosity####
    dbolodR1 = [(8*np.pi*(R1*const.R_sun)*BFLCTeff1**4)/(R1*u.R_sun)]**2
    dbolodR2 = [(8*np.pi*(R2*const.R_sun)*BFLCTeff2**4)/(R2*u.R_sun)]**2
    dbolodT1 = [(16*np.pi*(R1*const.R_sun)**2*BFLCTeff1**3)/(BFLCTeff1)]**2
    dbolodT2 = [(16*np.pi*(R2*const.R_sun)**2*BFLCTeff2**3)/(BFLCTeff2)]**2

    bololum1_err = sqrt(dbolodR1+dbolodT1)
    bololum2_err = sqrt(dbolodR2+dbolodT2)

    print('Bolometric Luminosity Star 1 =', bololum1, '+/-', bololum1_err)
    print('Bolometric Luminosity Star 2 =', bololum2, '+/-', bololum2_err)

##########################################################################################
##################################### mvbolo.py ##########################################

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

newfig = plt.figure(figsize=(7,6))

for idx, (starId,  M1,   T1,  M1_err,   T1_err,  M2,   T2,  M2_err,   T2_err) in enumerate(zip(
          starIds, kM1s, T1s, kM1_errs, T1_errs, kM2s, T2s, kM2_errs, T2_errs)):

    if starId != 6131659 and starId != 6781535 and starId != 6449358:
        plt.errorbar(bololum1.value, M1.value, xerr=(bololum1_err.value), yerr=M1_err.value,
                     marker='o', markersize=8, markeredgewidth=2, ls='None',
                     c=starcolors[idx], label=starId, zorder=2)

        plt.errorbar(bololum2.value, M2.value, xerr=(bololum2_err.value), yerr=M2_err.value, 
                     ls='None', marker='o', markersize=8, markeredgewidth=2, 
                     markerfacecolor='None', c=starcolors[idx], label='_nolegend_', zorder=2)
    
for idx, isofile in enumerate(isofiles):
    
    isoMasses, isoRadii, isoTemps, isoLoggs = isochroneParse('isochrones_plot', isofile)
    
    plt.plot(np.log10(isoTemps.value), isoMasses.to_value(u.Msun), ls=isolines[idx],
             lw=2, c=isocolors[idx], label=isolabels[idx], zorder=1)

plt.xlabel('Bolometric Luminosity ($L_{\odot}$)', size=14)
plt.ylabel('Mass ($M_{\odot}$)', size=16)
plt.legend(frameon=False, fontsize=12)
plt.xlim([3.9, 3.6])
plt.ylim([0.5, 2.2])

plt.show()
