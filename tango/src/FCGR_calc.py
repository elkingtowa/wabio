# FCGR plots show the count of each kmer in the genome as a pixel in the heatmap
# calculate for each phage here and then probably plot in R
from tetranucleotideAnalysis import *
import numpy as np 
import matplotlib.pyplot as plt

# get phage and fastas
nameFile = 'C:/Users/Admin/Documents/GitHub/tango/data/phagesDB/sequenced_phage_map_windows.txt'
with open(nameFile,'r') as nf:
	names=[]
	fileNames=[]
	line = nf.readline().strip().split('\t')
	while line != ['']:
		names.append(line[0])
		fileNames.append(line[1])
		line = nf.readline().strip().split('\t')

# do for each phage, save results to file
fOutFile = 'C:/Users/Admin/Documents/GitHub/tango/data/with_reverse_complement/FCGR_all_frequency.tsv'
pOutFile = 'C:/Users/Admin/Documents/GitHub/tango/data/with_reverse_complement/FCGR_all_probability.tsv'
with open(fOutFile,'w') as of:
	# write header
	kmerList = enumerateKmers(4)
	kmerString = ''
	for kmer in kmerList:
		kmerString += kmer + '\t'
	of.write(kmerString+'\n')

	#write data for each file
	for name, fileName in zip(names, fileNames):
		print 'working with ' +name
		fDict = dokmerCount(fileName, 4, RC=True)
		valString = name
		for kmer in kmerList:
			valString += '\t'+str(fDict[kmer]) 
		of.write(valString+'\n')

with open(pOutFile,'w') as of:
	# write header
	kmerString = ''
	for kmer in kmerList:
		kmerString += kmer + '\t'
	of.write(kmerString+'\n')

	#write data for each file
	for name, fileName in zip(names, fileNames):
		print 'working with ' +name
		pDict = dokmerCount(fileName, 4, probability=True, RC=True)
		valString = name
		for kmer in kmerList:
			valString += '\t'+str(pDict[kmer]) 
		of.write(valString+'\n')
