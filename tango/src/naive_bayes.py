# Building a naive Bayes classifier for clusters!

import scipy.stats
import math
import os
#Files we need 
meanFile = 'C:/Users/Admin/Documents/GitHub/tango/data/with_reverse_complement/NB_clusterMean.tsv'
varFile = 'C:/Users/Admin/Documents/GitHub/tango/data/with_reverse_complement/NB_clusterVar.tsv'
freqFile = 'C:/Users/Admin/Documents/GitHub/tango/data/with_reverse_complement/NB_clusterNums.tsv'

#read each file 
clusterNames = []
tetraNames = []
clusterMeans = []
clusterVars = []
clusterFreqs = []

with open(meanFile,'r') as mf, open(varFile,'r') as vf, open(freqFile,'r') as ff:
	mline = mf.readline().strip().split('\t')
	vline = vf.readline().strip().split('\t')
	fline = ff.readline().strip().split('\t')
	tetraNames = mline
	mline = mf.readline().strip().split('\t')
	vline = vf.readline().strip().split('\t')
	fline = ff.readline().strip().split('\t')
	while mline != ['']:
		clusterNames.append(mline[0])
		clusterMeans.append([float(m) for m in mline[1:]])
		clusterVars.append([float(v) for v in vline[1:]])
		clusterFreqs.append(int(fline[1]))
		mline = mf.readline().strip().split('\t')
		vline = vf.readline().strip().split('\t')
		fline = ff.readline().strip().split('\t')


seq='TGCGGCTGCCAGATCGTGTACGGGTTTGGAAGTCGACGGAGGGAACAGCGCGGGCCTAGAAGGCCCCGTAATGCCCCCTGAGAGCCCCGTAGACGGACGAACGGTGCGGATCGATAGATGGCACCGGAGACAAGCGAAGACGGCCGCAGAGCCGTCGCCGGCTGACGCCCGCGTAGGAAGATACTCATGTGAAGTGCGTCACATTCTACGGGTGAAACGCGAAAGTGGAAGGTTCCTTACCTATGGAGGGGTAAGGGAGCGAGCTCCAGCGAGCGACCGCACCCCGACATAGGTTCTTGTCGGGGTAGTCGAACGGAGAGAGACTACCCCTTTTAGCGACCTCCGGTCGCCCAGGTAGGTACCGAACGATGAGTGAGGTACCAGACCGTACAGGCCGGGGTTTATCCCCCGGCCGATACAGCATGGTCGTTTTGGGTAGGTACGTTACGTAAGCATCACTCACCAACAGAACCAGTGGTTACGTAACCGGGTACGTTACGTACGACTAGATACGTAACAGAACCACTAAACCCGTGGCCGCCCGGAAGGCGGCCCGCAGCGGGTTACGTTTGTCAGGTATGTCACCTTGACACGTAACAGATCGAGATTTCTGCCGGACTCTGTGGGTCCGGTAGCGCCCCTGGCGGCTAGCTACCGCGCGAGGGGCCTCTAACCATCCACAGAATGGAGCTCCGCATGGAGAAGAAGTGTTCTTTCGAAGGCTGTCCCAAGCCAGTACGCACTCGAGGCTGGTGCAACGGTCACTACCTGCAGCAGTGGTCAGGTAAGCCGATGATGCCTCTCAAGAACGTCCACCCTGAAACCTGCACCCTCGACTTCTGCGACCGACCCTATCGGTCGAACGGGCTGTGCGAAGGCCACTACTACCAGGCGCGTCGTGGGCGGCCGCTCAAGCCGCTACGTATCACCAGGCATCACCAGCCAGTACCTCCTTGCACGGTGGTCGGCTGCGATATCGAGGCCGTGTCTCACGATCGAC'
os.chdir('C:/Users/Admin/Documents/GitHub/tango/src')
from tetranucleotideAnalysis import *
a=TUDFromString(seq,4,enumerateKmers(4))
at=a.values()
# for b in a.values():
# 	if b>0:
# 		at.append(math.log(b,2))
# 	else:
# 		at.append(-99)

test = clusterMeans[0]

probabilities = []
for cluster in range(len(clusterNames)):
	cprob = 0
	for tetra, value in zip(range(len(tetraNames)), at):
		if clusterMeans[cluster][tetra] > -10:
			if scipy.stats.norm(clusterMeans[cluster][tetra], math.sqrt(clusterVars[cluster][tetra])).pdf(value) >0.0:
				cprob = cprob + math.log(scipy.stats.norm(clusterMeans[cluster][tetra], math.sqrt(clusterVars[cluster][tetra])).pdf(value))
		#print "c: " + str(clusterNames[cluster]) + "t: " + tetraNames[tetra] + str(cprob)
	probabilities.append(cprob*(clusterFreqs[cluster]/(float(sum(clusterFreqs)))))

