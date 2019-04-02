import numpy as np
import astropy.units as u
import astropy.constants as const

'''
Quick luminosity calculation for 4 SB2 EBs, with error propagation

Luminosity = 4 pi R^2 sigma T^4
'''

starids = ['528_1', '528_2', '686_1', '686_2', '677_1', '677_2', '428_1', '428_2']

teffs = [6845., 6776., 6497., 6439., 6822., 6135., 5684., 5694.] * u.K
radii = [2.003, 1.679, 1.655, 1.455, 1.748, 0.998, 1.033, 1.026] * u.solRad

tefferrs = [328., 478., 159., 274., 162., 358., 146., 243.] * u.K
radiierrs = [0.062, 0.087, 0.013, 0.012, 0.009, 0.005, 0.012, 0.011] * u.solRad

lums = []
lumerrs = []

for starid, rad, teff, raderr, tefferr in zip(starids, radii, teffs, radiierrs, tefferrs):
    ddr = 8 * np.pi * rad.to(u.cm).value * const.sigma_sb.cgs.value * teff.value**4
    ddt = 16 * np.pi * (rad.to(u.cm).value)**2 * const.sigma_sb.cgs.value * (teff.value)**3
    
    lums.append(4*np.pi*rad*rad*const.sigma_sb*(teff**4))
    lumerrs.append(np.sqrt( (raderr.to(u.cm).value * ddr)**2 + (tefferr.value * ddt)**2 ))


for starid, lum, lumerr in zip(starids, lums, lumerrs):
    print('{0}: {1:.3e} +/- {2:.3e} erg/s'.format(starid, lum.cgs.value, lumerr))

for starid, lum, lumerr in zip(starids, lums, lumerrs):
    print('{0}: {1:.3f} +/- {2:.3f} [log L in Lsun units]'.format(starid, np.log10(lum/u.solLum), 0.434*lumerr/lum.cgs.value))
