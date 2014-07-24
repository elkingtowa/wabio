# compute kmer counts for a list of files, save each result to an individual matrix file
from tetranucleotideAnalysis import *
import numpy as np 
import matplotlib.pyplot as plt
import sys

#make into a command line script
assert len(sys.argv) ==6, "Usage: python window_kmer_count.py nameFile outFile k windowSize stepSize \n \
 						  nameFile is a tab delimited database of phage names and fasta locations \n \
 						  outFile is the filename to save the resulting statistics \n \
 						  k is the size of the oligonucleotide to count \n \
 						  windowSize is the length of the window to count in. Set to 0 for entire genome \n \
 						  stepSize is the number of bases to increment the window. Ignored if windowSize is 0"
nameFile = sys.argv[1]
outFile = sys.argv[2]
k = int(sys.argv[3])
windowSize = int(sys.argv[4])
stepSize = int(sys.argv[5])

def window_kmer_count(nameFile, outFile, k, windowSize, stepSize):
	# read names and fasta locations 
	with open(nameFile,'r') as nf:
		names=[]
		fileNames=[]
		line = nf.readline().strip().split('\t')
		while line != ['']:
			names.append(line[0])
			fileNames.append(line[1])
			line = nf.readline().strip().split('\t')
	# count kmers and write each line to outFile 
	with open(outFile,'w') as of:
		# write header of kmer names
		kmerList = enumerateKmers(k)
		kmerString = ''
		for kmer in kmerList:
			kmerString += '\t' + kmer 
		of.write(kmerString+'\n')
		# do for each phage
		for name,fileName in zip(names, fileNames):
			print "Working with: " + name
			# compute kmers in the entire genome
			if windowSize ==0:
				data=doKmerCount(fileName, k, probability=False, RC=False)
				toWrite = name
				for kmer in kmerList:
					toWrite += '\t' + str(data[kmer])
				of.write(toWrite+'\n')

			# compute kmers in sliding window
			else:
				data=doKmerCountWindows(fileName, k, windowSize, stepSize, probability=False)
				for dataDict, window in data:
					toWrite = name+'_'+str(window[0])+':'+str(window[1])
					for kmer in kmerList:
						toWrite += '\t'+str(dataDict[kmer])
					of.write(toWrite+'\n')

def main():
	window_kmer_count(nameFile, outFile, k, windowSize, stepSize)

if __name__ == '__main__':
	main()

# get phage and fastas
#nameFile = 'C:/Users/Admin/Documents/GitHub/tango/data/TDI_individual_clusters/windows/sequenced_phage_map_M1.txt'
#nameFile = 'C:/Users/Admin/Documents/GitHub/tango/data/phagesDB/sequenced_phage_map_windows.txt'


# SET WINDOW AND STEP
#windowSize=5000
#stepSize=2500

#save all phage windows to one file
#allOutFile = 'C:/Users/Admin/Documents/GitHub/tango/data/window_kmer_count/all_freqs_2000_k4.tsv'
