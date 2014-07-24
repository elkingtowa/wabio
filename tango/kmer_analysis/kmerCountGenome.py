# quick script to count the numbers of kmers across multiple genomes. Saves results for all in outFile

import sys
import os
from kmerAnalysis import doKmerCount

if len(sys.argv) != 4:
	sys.exit("USAGE: kmerCountWindows.py fastaMap k outFile")

fastaMap = sys.argv[1]
k = int(sys.argv[2])
outFile = sys.argv[3]

# read fasta information
names=[]
fnames=[]
with open(fastaMap, 'r') as nf:
	line = nf.readline().strip().split(',')
	while line != ['']:
		names.append(line[0])
		fnames.append(line[1])
		line = nf.readline().strip().split(',')

with open(outFile, 'w') as of:
	first = True
	for name, fname in zip(names, fnames):
		data = doKmerCount(fname, k)
		if first:
			# write header of file - names of tetras
			of.write(','.join(data.keys())+'\n')
		# write kmer counts 
		of.write(name+','+','.join([str(x) for x in data.values()]) + '\n')
		first = False