# quick script to count the numbers of kmers in a sliding window across multiple genomes. Saves results for each genome in a separate file in the outFolder

import sys
import os
from kmerAnalysis import doKmerCountWindows, doKmerCount

if len(sys.argv) != 6:
	sys.exit("USAGE: kmerCountWindows.py fastaMap windowSize stepSize k outFolder")

fastaMap = sys.argv[1]
windowSize = int(sys.argv[2])
stepSize = int(sys.argv[3])
k = int(sys.argv[4])
outFolder = sys.argv[5]

# read fasta information
names=[]
fnames=[]
with open(fastaMap, 'r') as nf:
	line = nf.readline().strip().split(',')
	while line != ['']:
		names.append(line[0])
		fnames.append(line[1])
		line = nf.readline().strip().split(',')


for name, fname in zip(names, fnames):
	data = doKmerCountWindows(fname, k, windowSize, stepSize)
	with open(os.path.join(outFolder, name+'_'+str(windowSize)+'_'+str(stepSize)+'_k'+str(k) + '.csv'), 'w') as of:
		# write header of file - names of tetras
		of.write(','.join(data[0][0].keys())+'\n')
		for windowWrap in data:
			of.write(str(windowWrap[1][0]) + ':' + str(windowWrap[1][1]) + ',' + ','.join([str(x) for x in windowWrap[0].values()]) + '\n')
