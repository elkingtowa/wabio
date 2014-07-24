

import pickle

# enumerateKmers(k): construct list of all k-mers over alphabet ATCG
# returns a list of the kmers.
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

# generateDatabase(input, k, fastaFile=True): generates a BLAST database for the single sequence contained in input for the integer k
# finds all locations of a k-mer and stores the locations in a dictionary
# used for faster lookup in later parts of BLAST. 
# if fastaFile is true (default), input contains a fastafile to be read and interpreted. 
# if fastaFile is false, input contains a string that needs to be hashed. 
# if saveName != None, saves a pickle of the database to a file
def generateDatabase(inputSeq, k, fastaFile=True, saveName=None):
	if fastaFile:
		#read fasta into sequence
		with open(inputSeq,'r') as f:
			sequence=''
			line=f.readline().strip()
			while line != '':
				if line[0] != '>':
					sequence +=line
				line=f.readline().strip()
	else: sequence = inputSeq
	#initalize database. probably don't need to do this but might be helpful in avoiding key errors later!
	kmerList = enumerateKmers(k)
	database = dict((k, []) for k in kmerList)

	# slide a window of length k along the genome and mark every position in the database.
	for i in range(len(sequence)-k+1):
		database[sequence[i:i+k]].append(i)

	# Save database to file if desired
	if saveName != None:
		with open(saveName, 'wb') as of:
			pickle.dump(database, of)

	return database

#a = generateDatabase('C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\homework\hw8\data\BRCA1-dna.txt',4, True)

def main():
	import argparse
	parser = argparse.ArgumentParser(description="This script can be used to save a database of a subject sequence so it doesn't have to be computed each time. This is really only useful for very large databases.")
	parser.add_argument("subject", action="store", metavar="subject", help="The location of the fasta file containing the subject sequence")
	parser.add_argument("k", action="store", type=int, metavar="k", help="Length of words to compute locations for")
	parser.add_argument("outFile", action="store", metavar="outFile", help="Filename to save the pickle dictionary to")
	args=parser.parse_args()
	generateDatabase(args.subject, int(args.k), saveName=args.outFile)
if  __name__ =='__main__':main()