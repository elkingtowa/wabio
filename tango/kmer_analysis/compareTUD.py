#parameters 
import argparse

parser = argparse.ArgumentParser(description="This script calculates Tetranucleotide Usage Deviation (TUD) for multiple phage genomes. Originally it could only do 4-mers, but has now been generalized to all nucleotide lengths. The defualt function of this script is to save a nexus distance file that can be used for tree buliding with splitstree and other programs. Be aware that the nexus format is picky about special characters in the names of dataa, like parenthases. If you want to specify a cluster as part of the name, use a dash (ie \"Dante-F1\") The only required inputs are the location of a configuration file and location to save the resulting nexus distance file. The configuration file is a comma separated file of phage info with 2 necessary fields and 2 optional fields: \n name,fastaPath,subsetStart,subsetEnd\nEach phage to be compared should be separated by a new line. Additional arguments provide more control over the configuration the underlying calculations. If you want to save a csv file of the resulting usage deviation data, specify one with the --s option. A distance matrix can be saved by specifying a file name with the --d option.")
parser.add_argument("fastaMap", action="store", metavar="fastaMap", help="The file name of the comma separated file defining phage information. The name and fastaPath fields are required, subsets can be defined with integers in the next two fields.")
parser.add_argument("nexusFile", action="store", metavar="nexusFile", help="The file to save the resulting nexus to.")
parser.add_argument("--subset", action="store", metavar="subset", default="False", help="set to True to plot only calculate for the regions defined in the input file. If True, each line must have integers in fields 3 and 4 that represent the genomic region of each phage to compare. This allows you to subset on regions containing genes, etc for each phage.")
parser.add_argument("--k", action="store", metavar="k", default="4", help="Can also use this to compute 2,3,5-mers, etc.")
parser.add_argument("--s", action="store", metavar="s", default='', help="specify a filename to save resulting usage deviation data to. Be careful because these files can get big for higher values of k!")
parser.add_argument("--d", action="store", metavar="d", default='', help="specify a filename to save resulting distance matrix to")
parser.add_argument("--r", action="store", metavar="r", default="True", help="By default sequences are extended by their reverse complement before counting kmers and devation. Set to false to change this behavior.")
args=parser.parse_args()

fastaMap=args.fastaMap
nexusFile=args.nexusFile
if args.subset == "True": subset=True
else: subset=False
k=int(args.k)
deviationFile = args.s
distanceFile = args.d
if args.r == "True": RC=True
else: RC=False

#### DONE WITH ARGUMENTS | ACTUAL CODE BELOW ####

from kmerAnalysis import doKmerCount, doZeroOrderExpected, nexusWriter
from scipy.spatial.distance import pdist
from collections import OrderedDict 
import os

def compareTUD(fastaMap, nexusFile, subset, k, deviationFile, distanceFile, RC):
	# read information in fastaMap
	names = []
	fnames = []
	subsets = []
	with open(fastaMap, 'r') as nf:
		line = nf.readline().strip().split(',')
		while line != ['']:
			if subset and (len(line) < 4):
				print "Specified to subset but didn't include start and end position"
				print "Subsetting will be disabled"
				subset = False
			names.append(line[0])
			fnames.append(line[1])
			# parse subset information
			if subset:
				subsets.append([line[2],line[3]])
			line = nf.readline().strip().split(',')

	# compute deviation for each defined phage
	devDict = dict()
	if subset:
		for fname,name,subset in zip(fnames, names, subsets):
			obs = doKmerCount(fname, k, subset, RC)
			exp = doZeroOrderExpected(fname, k, subset) 
			devDict[name] = [float(obs[kmer])/exp[kmer] for kmer in obs.keys()]
	else:
		for fname, name in zip(fnames, names):
			obs = doKmerCount(fname,k, RC=RC)
			exp = doZeroOrderExpected(fname, k, RC=RC) 
			devDict[name] = [float(obs[kmer])/exp[kmer] for kmer in obs.keys()]

	# commpute distances and convert list to a symmetric distance matrix 
	# have to make sure order of names is preserved
	toDist = []
	for name in names:
		toDist.append(devDict[name])
	distances = pdist(toDist, 'euclidean')
	distanceMat = [[0 for x in range(len(names))] for x in range(len(names))]
	count = 0
	for i in range(len(names)):
		for j in range(i+1, len(names)):
			distanceMat[i][j] = distances[count]
			distanceMat[j][i] = distances[count]
			count+=1

	#done with calculations, start to write data
	# we always write a nexus file. 
	if os.path.isfile(nexusFile):
		print 'Warning: Nexus file exists and will be overwritten'
	nexusWriter('d', nexusFile, dataDict=OrderedDict((name, distance) for (name, distance) in zip(names, distanceMat)))

	# save deviation information if desired
	if deviationFile != '':
		if os.path.isfile(deviationFile):
			print 'Warning: deviation file exists and will be overwritten'		
		with open(deviationFile, 'w') as df:
			# wrtie names of kmers on first line
			kmerLine = '' 
			for kmer in obs.keys():
				kmerLine += kmer + ','
			df.write(kmerLine[:-1] + '\n')

			# write line for each phage
			for (name, deviations) in devDict.items():
				line = name
				for dev in deviations:
					line += ',' + str(dev)
				df.write(line+'\n')

	# save distance matrix if desired
	if distanceFile != '':
		if os.path.isfile(distanceFile):
			print 'Warning: distance file exists and will be overwritten'
		with open(distanceFile, 'w') as df:
			# write names on first line
			nameLine = ''
			for name in names:
				nameLine += name +','
			df.write(nameLine[:-1]+'\n')

			# write each line of distance matrix
			for name,distances in zip(names, distanceMat):
				toWrite = name
				for distance in distances:
					toWrite += ',' + str(distance)
				df.write(toWrite+'\n')

	print 'Done!   :)'

if __name__ == '__main__':
	compareTUD(fastaMap, nexusFile, subset, k, deviationFile, distanceFile, RC)