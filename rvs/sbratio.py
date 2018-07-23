'''This program is a translation of Keivan Stassun's IDL code to compute Temp ratios 
with surface brightnesses for Spectroscopic Eclipsing Binaries. To run this program,
execute it with the following input arguments: 
python sbratio.py 0 6237 0.3920984216 0 '''
##################sb ratio :: teff 1 APOCASK :: SBRatio :: Tratio???
import numpy as np
from scipy.interpolate import interp1d
from sys import argv

###   Read in the Dartmouth isochrone file   ###
mass, logteff, logg, kp = np.loadtxt('afep0.txt', unpack=True, usecols = (1, 2, 3, 13))
sbratio_Kepler = float(argv[1])
teff1 = float(argv[2])
sbratio = float(argv[3])
tratio = float(argv[4])
if len(argv)==6:
	age1 = float(argv[5])
if len(argv)==7:
	met = float(argv[6]) 
	age1 = float(argv[5])

###   Define all the things   ###
age1 = 1000
met = 0                                               #Metallicity
feh = str(abs(met)).zfill(3)
#teff1 = 4000 #for each SEB
#feh = 
a = str(age1).zfill(5)
m = np.where(logg >= 4.1)
mass = mass[m]
logteff = logteff[m]
logg = logg[m]
kp = kp[m]                                            #magnitude
r2 = mass/(10**logg)
sb = 10**(-kp/2.5)/r2                                 #r^2 came from Gm/r^2 = g
#teff1 = 10**(logteff[0])
#sbratio = 1.4                                        #Which is on top P or S???

###   Math Party!   ###
sb1 = interp1d(logteff, sb)
sb1 = sb1(np.log10(teff1))

###   Fit Party!   ###
#fit = np.where(abs(sb/sb1-sbratio) == min(sb/sb1-sbratio))
fit = np.argmin(np.abs(sb/sb1-sbratio))
#fit = fit[0]

###   Temp Ratio Party!   ###
tratio = 10**logteff[fit]/teff1     #Secondary/Primary
#teff2 = aspcapTeff * tratio
print('tratio =', tratio)
print('teff1 =', teff1)
#print('teff2 =', teff2)
###   Plotting Party!   ###
