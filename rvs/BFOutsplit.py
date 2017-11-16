import numpy as np

# read in the big data file where each row is one long string

bulk_data = np.loadtxt('data/5285607/5285607BFOut.txt', dtype=str, delimiter=',', comments='# ')

# printing it out so you can see what it looks like
print(bulk_data)

# identify indices of ### that separate each visit
visit_stops = np.where(bulk_data == '###')[0]

# print out the indices... b/c of the way the outfile is written, each visit will have 2 rows of ###
print(visit_stops)

visit_stops = np.append(visit_stops[::2], len(bulk_data))

nvisit=0
# loop through each visit chunk, where i,j are the indices corresponding to the start and end of visit
for i,j in zip(visit_stops[:-1], visit_stops[1:]):
	nvisit+=1
	# you can change the outfile naming convention, etc...
	np.savetxt('5285607BFOut_visit'+str(nvisit)+'.txt', bulk_data[i:j], fmt='%s')
	
# now you can go to the folder that you've saved everything in, and open up each visit individually, and
# read in the files individually for what you need to do