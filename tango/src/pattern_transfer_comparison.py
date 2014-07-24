from tetranucleotideAnalysis import *

#Working to implement the statistics in Liu et al. 2011
# alignment-free sequence comparison under the pattern transfer model 

#nullProbability(kmer, seq1, seq2): computes the probabilty of seeing the word under a null model. 
# use zero order markov with the two sequences appended for now. 
def nullProbability(kmer, seq1, seq2):
	sequence = seq1+seq2
	# count frequency of each letter in sequence and kmer.
	eA = sequence.count('A')/float(len(sequence))
	eT = sequence.count('T')/float(len(sequence))
	eC = sequence.count('C')/float(len(sequence))
	eG = sequence.count('G')/float(len(sequence))
	kA = kmer.count('A')
	kT = kmer.count('T')
	kC = kmer.count('C')
	kG = kmer.count('G')
	pKmer = (math.pow(eA,kA)*math.pow(eT,kT)*math.pow(eC,kC)*math.pow(eG,kG))
	return pKmer

#First, calculate Dstar2
kmerList = enumerateKmers(4)
#Dstar2(seq1, seq2, k, kmerList): computes the Dstar2 statistic for relatedness
# between 2 sequences. Takes 2 sequences, integer k, and enumerated lis tof kmers
# Returns (float of relatedness of seq1 and seq2, XtwDict, YtwDict, PwDict) 
# use first element of return for statistics. Use other for calculating Cstar2 efficiently
def Dstar2(seq1, seq2, k, kmerList):
	statDict = {}
	XtwDict = {}
	YtwDict = {}
	PwDict = {}
	#calculate statistic for each kmer
	for kmer in kmerList:
		# get probability of kmer in null model 
		Pw = nullProbability(kmer, seq1, seq2)
		PwDict[kmer] = Pw
		# count occurances of kmer in each seq
		reg =  r'(?=('+kmer+'))'
		Xw = float(len(re.findall(reg, seq1)))
		Yw = float(len(re.findall(reg, seq2)))
		# compute centralized count
		Xtw = Xw -(len(seq1)-k+1)*Pw
		Ytw = Yw -(len(seq2)-k+1)*Pw
		XtwDict[kmer] = Xtw
		YtwDict[kmer] = Ytw
		# finally, write to statDict
		statDict[kmer]= (Xtw * Ytw)/((len(seq1)*len(seq2))*Pw)
	return (sum(statDict.values()), XtwDict, YtwDict, PwDict)

#Cstar2(seq1, seq2, k, kmerList): computes the normalized Cstar2 statistic 
# input and output same as Dstar2
def Cstar2(seq1, seq2, k, kmerList):
	D = Dstar2(seq1,seq2, k, kmerList)
	numerator = ((len(seq1)*len(seq2))*D[0])
	#sum X and Y terms in denominator
	Xsum = 0
	Ysum = 0
	for kmer in kmerList:
		Xsum += math.pow(D[1][kmer],2)/D[3][kmer]
		Ysum += math.pow(D[2][kmer],2)/D[3][kmer]
	denominator = math.sqrt(Xsum)*math.sqrt(Ysum)
	# divide and return!
	return numerator / denominator

#RstarSum(seq1, seq2, k, windowSize, stepSize): calcluates the RstarSum statistic
# from the paper. Windows of size windowSize are moved along by stepSize
def RstarSum(seq1, seq2, k, windowSize, stepSize):
	kmerList = enumerateKmers(k)
	windows1 = int(math.floor((len(seq1)-windowSize)/stepSize))
	windows2 = int(math.floor((len(seq2)-windowSize)/stepSize))
	# initialize array to hold values
	mArray = np.zeros(shape=(windows1,windows2))
	# compute Cstar2 for each pair of windows
	for i in range(windows1):
		for j in range(windows2):
			print "calculating for window i=" + str(i)+" j="+ str(j)
			mArray[i,j] = Cstar2(seq1[i*stepSize:i*stepSize+windowSize],seq2[j*stepSize:j*stepSize+windowSize],k,kmerList)

	#initalize Xstar and Ystar arrays
	Xarray = np.zeros(shape=windows1)
	Yarray = np.zeros(shape=windows2)
	#compute max and store in Xstar and Ystar arrays
	for i in range(windows1):
		Xarray[i] = max(mArray[i,:])
	for j in range(windows2):
		Yarray[j] = max(mArray[:,j])

	#RhatSum is sum of Xarray and Yarray
	RhatSum = sum(Xarray) + sum(Yarray)
	#normalize and return
	return RhatSum/math.floor((len(seq1)*len(seq2)/stepSize))




import os
os.chdir('C:/Users/Admin/Documents/GitHub/tango/src')
f1 = 'C:/Users/Admin/Desktop/Dropbox/College/Junior/Phage Hunters TA/fasta/phage/Abrogate.fasta'
f2 = 'C:/Users/Admin/Desktop/Dropbox/College/Junior/Phage Hunters TA/fasta/phage/Aeneas.fasta'
for seq_record in SeqIO.parse(f1, "fasta"):
	seq1 = seq_record.seq.tostring().upper()
for seq_record in SeqIO.parse(f2, "fasta"):
	seq2 = seq_record.seq.tostring().upper()
