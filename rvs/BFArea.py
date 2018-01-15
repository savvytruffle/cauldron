from __future__ import print_function

import numpy as np
from scipy.integrate import simps
from numpy import trapz


# The y values.  A numpy array is used here,
# but a python list could also be used.
y = np.array([5, 20, 4, 18, 19, 18, 7, 4])

### Read in all the things! ###

uncorrRV, BFamp, Gaussian = np.loadtxt('data/5285607/5285607BFOut_visit1.txt',
	usecols=(0,1,2),unpack=True)

### Input file alternates between primary and secondary peak	
	#with open('data/5285607/5285607BFOut.txt', 'r') as f:
    #for count, line in enumerate(f, start=1):
    #    if count % 2 == 0:
    #        print(line)
        

# Compute the area using the composite trapezoidal rule.
area = trapz(BFamp, dx=5)
print("Trapezoidal area =", area)

# Compute the area using the composite Simpson's rule.
area = simps(BFamp, dx=5)
print("Simpson area =", area)