import os
import sys
import numpy as np

# set name of star
try:
    KIC = sys.argv[1]
except:
    KIC = '5285607'

# read in the big data file where each row is one long string
file = os.path.join('data', KIC, KIC + 'BFOut.txt')
bulk_data = np.loadtxt(file, dtype=str, delimiter=',', comments='# ')

# save other useful info from the BFOut file
metadata = []
with open(file) as f:
    for idx, line in enumerate(f):
        if '# ' in line:
            metadata.append(line)

# printing it out so you can see what it looks like
#print(bulk_data)

# identify indices of ### that separate each visit
visit_stops = np.where(bulk_data == bulk_data[0])[0]

# print out the indices... b/c of the way the outfile is written, each visit will have 2 rows of ###
#print(visit_stops)

visit_stops = np.append(visit_stops[::2], len(bulk_data))

nvisit = 0
# loop through each visit chunk, where i,j are the indices corresponding to the start and end of visit
idx = 0
for i,j in zip(visit_stops[:-1], visit_stops[1:]):
	idx1 = idx
	idx2 = idx+1
	idx3 = idx+2
	idx4 = idx+3
	nvisit += 1
	idx += 4
	# you can change the outfile naming convention, etc...
	print('File created: ', KIC + 'BFOut_visit{}.txt'.format(nvisit))
	header = str(metadata[idx1] + metadata[idx2] + metadata[idx3] + metadata[idx4])
	np.savetxt(KIC + 'BFOut_visit' + str(nvisit) + '.txt', bulk_data[i:j], fmt='%s', header=header)
	
# now you can go to the folder that you've saved everything in, and open up each visit individually, and
# read in the files individually for what you need to do
