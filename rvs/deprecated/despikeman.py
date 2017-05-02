import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


#Where the despiked files live

despikelist = 'data/5284133/5284133infilesdespike.txt'
manuallist = 'data/5284133/5284133infilesman.txt'
f1 = open(despikelist)
f2 = open(manuallist)
wavedespikelist = despikelist[0]; specdespikelist = despikelist[1]
yoffset = 0

for line in f1:
#for i in range (0, 3):
#wave, spec of despiked .txt
    despikelist = line.rstrip()
    #despikespec = np.interp(wavedespikelist[i], specdespikelist[i])
    #despikelist.append(despikespec)
    #infilelist.append(infile)
    yoffset = yoffset + 1
    wavedespike, specdespike = np.loadtxt(despikelist, usecols=(0,1), unpack=True)
    #print(wavedespike, specdespike)
f1.close()
    
#wave, spec of manual despiked .txt
for line in f2:
    manuallist = line.rstrip()
    waveman, specman = np.loadtxt(manuallist, usecols=(0,1), unpack=True)
f2.close()
    
#plt.plot(despikespec+yoffset, color ='r')
plt.plot(wavedespike, specdespike, color='r')
plt.plot(wavedespike+yoffset, specdespike+yoffset, color='r')
plt.plot(waveman, specman, color='b')
plt.plot(waveman+yoffset, specman+yoffset, color='b')
plt.show()


