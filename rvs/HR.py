'''This program does HR diagram stuff'''

import numpy as np
from scipy.interpolate import interp1d
from sys import argv
import matplotlib.pyplot as plt

###   Read in the Dartmouth isochrone file   ###
isofile = 'afep0.txt'
mass, logteff, logg, kp = np.loadtxt(isofile, unpack=True, usecols = (1, 2, 3, 13))
teff1 = float(argv[1])
sbratio = float(argv[2])

mone = float(argv[3])
mtwo = float(argv[4])

rone = float(argv[5])
rtwo = float(argv[6])

########################### Read things from the isochrone file ##############################
m = np.where(logg >= 4.1)  # a print(m) statement will print all of the masses where the 
                           # log(g) is greater than 4.1. This isolates members of the main
                           # sequence. 
mass = mass[m]
logteff = logteff[m]
logg = logg[m]
kp = kp[m]                                            #magnitude
r2 = mass/(10**logg)
sb = 10**(-kp/2.5)/r2                                 #r^2 came from Gm/r^2 = g

# print some stuff out for funsies #
print (r2)

###   Math Party!   ###
sb1 = interp1d(logteff, sb)# this interpolates the relationship between the surface bright-
                           # ness and the effective temperature (mass-luminosity relation)
                           # which only applies to main sequence stars, which is one of the 
                           # reasons we need line 21



sb1 = sb1(np.log10(teff1)) # this calculates the surface brightness of the primary, which is 
                           # assumed to be the surface brightness of the object, because the 
                           # secondary does not contribute much to the flux observed

###   Fit Party!   ###
fit = np.argmin(np.abs(sb/sb1 - sbratio))  # this fits a line to the points at which the 
                                           # the surface brightness calculated for our target
                           # and the surface brightness ratio that was read in are closest
                           # to one another. This is then plotted, as it is the evolutionary
                           # track associated with our particular target. 


#fit = fit[0]

###   Ratio Party!   ###
tratio = 10**logteff[fit]/teff1            # this calculates the temperature ratio using 10
                                           # to the power of (log(Teff[fit]) using the fit 
                           # discussed above. This ratio comes from our calculation of the 
                           # surface brightness, which comes from our calculation of radius 
                           # of the secondary. This ratio is therefore (Teff of the secondary/
                           # Teff of the primary.)

teff2 = teff1 * tratio     # this line then calculates the effective temperature of the 
                           # secondary using the primary Teff and the ratio calculated above. 

loggone = 4.391            # these are for my troubleshooting, attempting to plot log(g) but 
                           # our calculated points do not fit on the HR diagram, but when we 
                           # declare loggone as done here, with the log(g)from the Dartmouth 
                           # file the point is on the track, therefore there is a problem 
                           # with how radius is calculated or how surface gravity is calculated. 

loggtwo = 4.391
print('You entered teff 1 as ', teff1)
print('You entered radius as ', rone)
print('Using the isochrone file ', isofile)
#print('r is ', r)

###   Plotting Party!   ###
#logt1 = np.log10(teff1)
#logt2 = np.log10(teff2)
#plt.plot(logt1, rone, 'r.', label='',) 
#plt.plot(logt2, rtwo, 'b.', label='',)
#plt.plot(logteff, r, label='',)
#plt.grid()
#plt.ylabel('R/$R_{\odot}$')
#plt.xlabel('log Teff')
#plt.title('5285607')
#plt.axis('tight')
#plt.show()
#plt.savefig('5285607HR.pdf')

#log(g) vs. log(teff) complete with Dartmouth evolutionary track
logt1 = np.log10(teff1)
logt2 = np.log10(teff2)
plt.plot(logt1, loggone, 'r.', label='',) 
plt.plot(logt2, loggtwo, 'b.', label='',)
plt.plot(logteff, logg, label='',)
plt.grid()
plt.ylabel('log(g)$')
plt.xlabel('log Teff')
plt.title('5285607')
plt.axis('tight')
plt.show()
plt.savefig('5285607HR.pdf')


