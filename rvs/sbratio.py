'''This program is a translation of Keivan Stassun's IDL code to compute Temp ratios 
with surface brightnesses for Spectroscopic Eclipsing Binaries. To run this program,
execute it with the following input arguments: 
python sbratio.py 0 6237 0.3920984216 0 '''

## More instructive: python sbratio.py teff1 sbratio
## then explain WTF are teff1 and sbratio, where should they come from?
## (note that sbratio_Kepler and tratio are never used in the program so I omitted them)
## Finally explain what the program does (prints some values out? saves them to file? makes a plot?)

import numpy as np
from scipy.interpolate import interp1d
from sys import argv

###   Read in the Dartmouth isochrone file   ###
isofile = 'afep0.txt'
mass, logteff, logg, kp = np.loadtxt(isofile, unpack=True, usecols = (1, 2, 3, 13))
##sbratio_Kepler = float(argv[1])  ## never used
teff1 = float(argv[1])
sbratio = float(argv[2])
##tratio = float(argv[4])  ## never used

## the result of the if statements below are also never used since you hardwire age and met below
## (do you want to read in age1 and met or have them hardwired? why?)
##if len(argv)==6:
##	age1 = float(argv[5])
##if len(argv)==7:
##	met = float(argv[6]) 
##	age1 = float(argv[5])

###   Define all the things   ###  ## none of these things are ever used, why are they here?
##age1 = 1000
##met = 0                                               #Metallicity
##feh = str(abs(met)).zfill(3)  ## never used
##a = str(age1).zfill(5)  ## never used

## Read things from the isochrone file
m = np.where(logg >= 4.1)  # what is this? what do you expect if you add a print(m) statement?
mass = mass[m]
logteff = logteff[m]
logg = logg[m]
kp = kp[m]                                            #magnitude
r2 = mass/(10**logg)
sb = 10**(-kp/2.5)/r2                                 #r^2 came from Gm/r^2 = g
#teff1 = 10**(logteff[0])  ## you have this line commented out! why? it changes the value of teff1...
#sbratio = 1.4                                        #Which is on top P or S???

## print some stuff out for funsies
print('index, mass, logteff, logg, kp')
print('------------------------------')
for idx, (midx, t, g, kep) in enumerate(zip(mass, logteff, logg, kp)):
    print(idx, midx, t, g, kep)
print(' ')

###   Math Party!   ###
sb1 = interp1d(logteff, sb)  ## what is this interpolating and why?
sb1 = sb1(np.log10(teff1))

###   Fit Party!   ###
#fit = np.where(abs(sb/sb1-sbratio) == min(sb/sb1-sbratio))
fit = np.argmin(np.abs(sb/sb1 - sbratio))  ## what are you fitting and why? CONVINCE ME :D
#fit = fit[0]

###   Temp Ratio Party!   ###
tratio = 10**logteff[fit]/teff1     #Secondary/Primary

teff2 = teff1 * tratio

print('You entered teff 1 as ', teff1)
print('You entered surface brightness ratio as ', sbratio)
print('Using the isochrone file ', isofile)
print('The calculated temperature ratio (T2/T1) is ', tratio)
print('This means teff 2 is ', teff2)

###   Plotting Party!   ###

## this isn't a very interesting party just yet... ;)
