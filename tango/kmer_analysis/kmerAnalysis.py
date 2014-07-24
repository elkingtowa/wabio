#import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.spatial.distance import pdist

## TOOLS FOR THE ANALYSIS OF TETRANUCLEOTIDE USAGE DEVAIATION IN GENOMES
## Designed primarily for phage genomes. May not be the fastest way, but it gets things done!

#parseFasta(fileName): parses a fasta file defined by fileName.
# separations between sequences must be delimited by a line starting with  '>'
# if only one sequence is definied in the file, a string of the sequence is returned.
# if multiple sequencces are defined, a list of strings is returned. 
def parseFasta(fileName):
	sequences = []
	with open(fileName, 'r') as ff:
		line = ff.readline().strip()
		currentSequence = ''
		# read through file
		while line != '':
			if line[0] != '>':
				#append to current seq
				currentSequence += line
				line = ff.readline().strip()
			else:
				#new sequence
				if currentSequence != '':
					sequences.append(currentSequence)
					line = ff.readline().strip()
					currentSequence	= ''
				else: 
					#case in the first sequence of the file
					line = ff.readline().strip()
		# append last seq
		sequences.append(currentSequence)
	if len(sequences) ==1:
		return sequences[0]
	else: 
		return sequences

# reverseComplement(pattern): returns the reverse complement of DNA string pattern. Equivalent to an inversion. 
def reverseComplement(pattern):
	chars = list(pattern)
	complement = ''
	lookup = dict({("A","T"),("T","A"),("C","G"),("G","C"),("N","N")})
	for base in chars:
		complement += lookup[base]
	rev_complement = complement[::-1]
	return rev_complement 

#enumerateKmers(k): construct list of all kmers over alphabet ATCG
def enumerateKmers(k):
	alphabet = ['A','T','C','G']
	kmerList = ['']
	while k > 0:
		working = kmerList[:]
		kmerList = []
		for kmer in working:
			for letter in alphabet:
				kmerList.append(kmer+letter)
		k -= 1 
	return kmerList

#TDI(): computes the k-mer difference index for a given genome. Takes in a fileName, size of window to compute TDI in, 
# step size to move along the genome by, k, a list of kmers of length k, and a boolean plot. 
# if plot is true, plot the resulting figure using matplotlib. 
# always returns a list of window positions and Zscores
# if subset is defined, only return information about the sequence between the two positions. 
def TDI(sequence, windowSize, stepSize, k):
	#get TUD dict to use in TDI comparison
	obs = kmerCount(sequence,k)
	exp = zeroOrderExpected(sequence, k) 
	tudDict = {kmer:float(obs[kmer])/exp[kmer] for kmer in obs.keys()}
	# slide along the genome in windows defined by windowSize with steps defined by stepSize
	start = 0
	end = windowSize
	windowIndex= []
	differencesWindows = []
	while end < len(sequence):
		if end > len(sequence):
			end = len(sequence)
		window = sequence[start:end]
		#record list of windows
		windowIndex.append([start,end])
		
		# calculate observed and expected in the window 
		# returns a dictionary 
		obs = kmerCount(window, k)
		exp = zeroOrderExpected(window, k)

		#compute difference sum for all kmers
		differenceSum = 0
		for kmer in obs.keys():
			differenceSum += abs((obs[kmer]/exp[kmer])-tudDict[kmer])
		differencesWindows.append(differenceSum)

		#slide along the window by stepSize
		start += stepSize
		end += stepSize

	# compute Z-score   Z=(x-mu)/sigma
	Zscores = []
	mu = np.mean(differencesWindows)
	sigma = np.std(differencesWindows)
	for w in range(len(differencesWindows)):
		Zscores.append((differencesWindows[w] - mu)/sigma)

	xaxis = [x for x,y in windowIndex]
	return [xaxis,Zscores]

#doTDI(fileName, windowSize, stepSize, k, subset=None, plot=False): a wrapper
# for the calculation of tetranucleotide deviation index in a sequence. 
# reads the single sequence defined in fileName, calculates TDI for k-mers in a sliding window
# defined by windowSize and stepSize. 
# can extract a subset of the sequence by passing a list of two integers to subset
# if plot=True, plots an example graph to the graphics device
# returns [xaxis positions, zscores]
def doTDI(fileName, windowSize, stepSize, k, subset=False):
	#parse fasta from fileName
	sequence = parseFasta(fileName)
	# pick out part defined in subset
	if subset != False:
		sequence = sequence[subset[0]:subset[1]]

	#do TDI calculation
	tdi = TDI(sequence,windowSize,stepSize,k)

	return [tdi[0],tdi[1]]


#zeroOrderExpected(sequence, k): computes the expected number of a kmer in a genome using a 
# zero order markov model. Takes as input a sequence and a value for k
# returns a dictionary with the kmers as keys and the expected number as values
def zeroOrderExpected(sequence, k):
	#construct a dictionary of all kmers
	kmerDict = dict((key, 0 ) for key in enumerateKmers(k))
	nucleotides = kmerCount(sequence, 1, probability=True)

	# calc expeccted number for each kmer in the dict
	for kmer in kmerDict.keys():
		kA = kmer.count('A')
		kT = kmer.count('T')
		kC = kmer.count('C')
		kG = kmer.count('G')
		eKmer = (math.pow(nucleotides['A'],kA) * \
				 math.pow(nucleotides['T'],kT) * \
				 math.pow(nucleotides['C'],kC) * \
				 math.pow(nucleotides['G'],kG))* \
				(len(sequence)-k+1)
		kmerDict[kmer] = eKmer
	return kmerDict

#doZeroOrderExpected(fileName, k, subset=None): a wrapper for zeroOrderExpected.
# computes the expected number of kmers given a zero order markov model. Reads from a 
# single sequence fasta file defined in fileName, computes expected for a given k, returns dict
# of kemr:expeced. a subset of the sequence can be cosen by specifying two integers in a list with the subset option
def doZeroOrderExpected(fileName, k, subset=None, RC=False):
	sequence = parseFasta(fileName)
	# append reverse complement if RC=True
	if RC:
		sequence += reverseComplement(sequence)
	if subset != None:
		sequence = sequence[subset[0]:subset[1]]
	return zeroOrderExpected(sequence, k)

#kmerCount(filename, k, probability=False): computes the number of each kmer in a sequence. 
# returns a dictionary of kmers as keys and counts as values
# if probability is true, returns the count/(n-k+1) as values
def kmerCount(sequence, k, probability=False):
	#construct a dictionary of all kmers
	kmerDict = dict((key, 0 ) for key in enumerateKmers(k))

	# slide a window of length k across the genome and add to dict
	for i in range(len(sequence)-k+1):
		if kmerDict.has_key(sequence[i:i+k]):
			kmerDict[sequence[i:i+k]] += 1
	# divide by number of kmers in the genome 
	if probability:
		for key in kmerDict.keys():
			kmerDict[key] = float(kmerDict[key])/(len(sequence)-k+1)
	return kmerDict

#doKmerCount(fileName, k, probability=False, RC=False, subset=None): calls kmerCount
# on a fasta file containing a single sequence definied in fileName
# returns a dictionary of kmers as keys and counts as values
# if probability is true, returns the count/(n-k+1) as values
# if RC is true, extends the sequcene by the reverse complement before counting
# a subset of the sequence can be cosen by specifying two integers in a list with the subset option
def doKmerCount(fileName, k, probability=False, RC=False, subset=None):
	sequence = parseFasta(fileName)
	# append reverse complement if RC=True
	if RC:
		sequence += reverseComplement(sequence)
	if subset != None:
		sequence = sequence[subset[0]:subset[1]]
	return kmerCount(sequence, k, probability)

#doKmerCountWindows: computes kemrCount at sliding windows across the genome. 
# returns (list of dictionaries one for each window, list of window start and end positions)
def doKmerCountWindows(fileName, k, windowSize, stepSize, probability=False):
	kmerList = enumerateKmers(k)
	#parse fasta from filename
	sequence = parseFasta(fileName)
	#compute number of windows
	windows = int(math.ceil(len(sequence)/float(windowSize)))
	
	start = 0 
	end = start+windowSize
	toReturn = []
	while end < len(sequence):
		#append data
		toReturn.append((kmerCount(sequence[start:end], k, probability),[start,end]))
		# get start and end positions
		start += stepSize
		end += stepSize
	#do one last calc for the last window 
	toReturn.append((kmerCount(sequence[start:end], k, probability),[start,end]))
	return toReturn

#GCcontent(sequence, windowSize, stepSize): computes GC content in a sliding window across the genome. 
# windowSize and stepSize set the size of the sliding window
# returns (list of gc percentages, list of window starts, average GC content)
def GCcontent(sequence, windowSize, stepSize):
	start = 0
	end = start+windowSize
	starts = []
	GCpercents = []
	while end < len(sequence):
		sub = sequence[start:end]
		GCpercents.append((sub.count('G') + sub.count('C'))/float(len(sub)))
		starts.append(start)

		start+=stepSize
		end += stepSize
	#do one last calc for the last window
	GCpercents.append((sub.count('G') + sub.count('C'))/float(len(sub)))
	starts.append(start)

	return(GCpercents, starts, (sequence.count('G') + sequence.count('C'))/float(len(sequence)))

#doGCcontent(fileName, windowSize, stepSize, plot=False, RC=False): computes GC content for a sequence in a fasta file.
# windowSize and stepSize set the size of the sliding window
# if plot=True, plots the resulting line graph
# if RC=True, passes the reverse complement of the sequence to the GC content function
def doGCcontent(fileName, windowSize, stepSize, plot=False, RC=False):
	sequence = parseFasta(fileName)
	if RC:
		sequence = reverseComplement(sequence)
	#do calc
	gc = GCcontent(sequence, windowSize, stepSize)
	if plot:
		import matplotlib.pyplot as plt
		plt.plot(gc[1],gc[0])
		plt.xlabel('genomic position')
		plt.ylabel('GC fraction')
		plt.show()
	return gc

#nexusWriter(mode, outFile, dataDict=None, inFile=None): a general purpose function for writing to nexus files.
# specify a function with the mode option: 'd' is stance mode. A distances block will be written in the 
# nexus file. 
# furure modes to add: 'c', write character data to the nexus file, 'b', write both
# Either takes in data in the format of a dictionary of names:list, or reads in a data 
# matrix from the file specified in inFile. 
# Needs testingto ensure comliance with nexus format
def nexusWriter(mode, outFile, dataDict=None, inFile=None):
	assert (dataDict!=None or inFile!=None), "You have to specify some kind of input data!"
	assert mode=='d', "Distance mode is the only valid option for now"
	
	if inFile != None:
		# read in input data from a matrix for each mode
		if mode =='d':
			with open(distanceFile, 'r') as tf:
				# first line is the tab separated names of the tetranucleotides
				tetraNames = tf.readline().strip().split(',')

				line = tf.readline().strip().split(',')
				dataDict = dict()
				while line != ['']:
					dataDict[line[0]] = line[1:]
					line = tf.readline().strip().split(',')
	# otherwise dataDict is already defined and we can move on to writing
	# write to nexus file 
	with open(outFile, 'w') as of:
		of.write('#NEXUS\n')
		#write Taxa block
		of.write('BEGIN TAXA;\n')
		of.write('DIMENSIONS NTAX=' + str(len(dataDict.keys())) + ';\n')
		taxlabels = 'TAXLABELS'
		for name in dataDict.keys():
			taxlabels += ' ' + name
		taxlabels += ';\n'
		of.write(taxlabels)
		of.write('END;\n')

		if mode == 'd':
			#write distance block
			of.write('BEGIN DISTANCES;\n')
			of.write('\tDIMENSIONS NTAX='+str(len(dataDict.keys()))+';\n')
			of.write('\tFORMAT\n') 
			of.write('\t\tTRIANGLE=BOTH\n')
			of.write('\t\tDIAGONAL\n')
			of.write('\t\tLABELS=LEFT\n')
			of.write('\t;\n')
			of.write('\tMATRIX\n')
			for name,distance in zip(dataDict.keys(), dataDict.values()):
				toWrite = '\t\t'+name
				for dist in distance:
					toWrite += ' '+ str(dist)
				of.write(toWrite+'\n')
			of.write('\t;\n')
			of.write('END;\n')


#tudToNexus(tudFile, outFile, parseClusters): converts a tab separated TUD file to the nexus format
# each tetranucleotide will be treated as a character. 
# if parseClusters is True, remove cluster designations that are in parenthases at the end of each label
# because nexus complains a lot!
# DEPRECIATED
def tudToNexus(tudFile, outFile, parseClusters):
	# read data  
	with open(tudFile, 'r') as tf:
		# first line is the tab separated names of the tetranucleotides
		tetraNames = tf.readline().strip().split('\t')

		line = tf.readline().strip().split('\t')
		names = []
		data = []
		while line != ['']:
			names.append(line[0])
			data.append(line[1:])
			line = tf.readline().strip().split('\t')
		# we can get the clusters from the names when they're defined within a  '()'
		if parseClusters:
			clusters = [n.split('(')[1].split(')')[0] for n in names]
			newNames = [n.split('(')[0] for n in names]
			names = newNames
	print 'read file'
	#start to write things to the nexus file
	with open(outFile, 'w') as of:
		of.write('#NEXUS\n')
		#write Taxa block
		of.write('BEGIN TAXA;\n')
		of.write('DIMENSIONS NTAX=' + str(len(names)) + ';\n')
		taxlabels = 'TAXLABELS'
		for name in names:
			taxlabels += ' ' + name
		taxlabels += ';\n'
		of.write(taxlabels)
		of.write('END;\n')

		#write character block
		of.write('BEGIN CHARACTERS;\n')
		of.write('\tDIMENSIONS [NTAX='+str(len(names))+'] NCAHR='+str(len(data[0]))+';\n')
		of.write('\t[FORMAT\n') 
		of.write('\t\t[DATATYPE=CONTINUOUS]\n')
		of.write('\t;]\n')
		of.write('\tMATRIX\n')
		for name,mat in zip(names, data):
			toWrite = '\t\t'+name
			for tud in mat:
				toWrite += ' '+ str(tud)
			of.write(toWrite+'\n')
		of.write('\t;\n')
		of.write('END;\n')

#distanceToNexus(distanceFile, outFile, parseClusters): converts a tab separated distance file to the nexus format
# pairwise distances between each element will be recorded in the distances block
# if parseClusters is True, remove cluster designations that are in parenthases at the end of each label
# because nexus complains a lot!
# DEPRECIATED
def distanceToNexus(distanceFile, outFile, parseClusters):
	# read data  
	with open(distanceFile, 'r') as tf:
		# first line is the tab separated names of the tetranucleotides
		tetraNames = tf.readline().strip().split('\t')

		line = tf.readline().strip().split('\t')
		names = []
		data = []
		while line != ['']:
			names.append(line[0])
			data.append(line[1:])
			line = tf.readline().strip().split('\t')
		# we can get the clusters from the names when they're defined within a  '()'
		if parseClusters:
			clusters = [n.split('(')[1].split(')')[0] for n in names]
			newNames = [n.split('(')[0] for n in names]
			names = newNames
	#start to write things to the nexus file
	with open(outFile, 'w') as of:
		of.write('#NEXUS\n')
		#write Taxa block
		of.write('BEGIN TAXA;\n')
		of.write('DIMENSIONS NTAX=' + str(len(names)) + ';\n')
		taxlabels = 'TAXLABELS'
		for name in names:
			taxlabels += ' ' + name
		taxlabels += ';\n'
		of.write(taxlabels)
		of.write('END;\n')

		#write character block
		of.write('BEGIN DISTANCES;\n')
		of.write('\tDIMENSIONS NTAX='+str(len(names))+';\n')
		of.write('\tFORMAT\n') 
		of.write('\t\tTRIANGLE=BOTH\n')
		of.write('\t\tDIAGONAL\n')
		of.write('\t\tLABELS=LEFT\n')
		of.write('\t;\n')
		of.write('\tMATRIX\n')
		for name,mat in zip(names, data):
			toWrite = '\t\t'+name
			for tud in mat:
				toWrite += ' '+ str(tud)
			of.write(toWrite+'\n')
		of.write('\t;\n')
		of.write('END;\n')

#converts a tab separated TUD file to the mega format
# each tetranucleotide will be treated as a character
# if parseClusters is True, eliminate the cluster designation within parenthases
# UNTESTED
def distance_to_mega(distanceFile, outFile, parseClusters):
	# read data  
	with open(distanceFile, 'r') as tf:
		# first line is the tab separated names of the tetranucleotides
		tetraNames = tf.readline().strip().split('\t')

		line = tf.readline().strip().split('\t')
		names = []
		data = []
		while line != ['']:
			names.append(line[0])
			data.append(line[1:])
			line = tf.readline().strip().split('\t')
		# we can get the clusters from the names when they're defined within a  '()'
		if parseClusters:
			clusters = [n.split('(')[1].split(')')[0] for n in names]
			newNames = [n.split('(')[0] for n in names]
			names = newNames
	print len(data)
	print len(data[0])
	#start to write things to the mega file
	with open(outFile, 'w') as of:
		of.write('#MEGA\n')
		of.write('!Title tetst;\n')
		of.write('!Description tetst;\n')
		#write Taxa block
		for name in names:
			of.write('#' + name + '\n')

		#write distances as lower half of matrix
		count=0
		for line in data:
			print count
			toWrite = ''
			for i in range(count):
				toWrite += str(line[i]) + ' '
			of.write(toWrite+'\n')
			count += 1
		


