#parameters 
import argparse

parser = argparse.ArgumentParser(description="Script to test contamination between a file of sequences and a file of vectors")
parser.add_argument("seqs", action="store", metavar="seqs", help="The file name of sequences")
parser.add_argument("vectors", action="store", metavar="vectors", help="The file name of the vectors")
args=parser.parse_args()

#parameters and paths specified in section above
seqs=args.seqs
vectors=args.vectors

#local_alignment_contamination(seq1,seq2): looks for alignments of length 11 with at most one mismatch between seq1 and seq2. Returns True if one is found, False otherwise.
def local_alignment_contamination(seq1, seq2):
	n = len(seq1)
	m = len(seq2)
	#matrix for keeping track of number of alignments and number of mismatches
	match_matrix = [[[0,0] for i in range(m+1)] for j in range(n+1)]
	for i in range(1, n+1):
		for j in range(1, m+1):
			if seq1[i-1] == seq2[j-1]: #match. add to match counter
				match_matrix[i][j][0] = match_matrix[i-1][j-1][0] + 1
				match_matrix[i][j][1] = match_matrix[i-1][j-1][1]
			else: #mismatch. add to mismatch counter
				match_matrix[i][j][0] = match_matrix[i-1][j-1][0]
				match_matrix[i][j][1] = match_matrix[i-1][j-1][1] + 1
			#remove i-11th and j-11th value.
			if (match_matrix[i][j][0]+match_matrix[i][j][1])>11:
				if seq1[i-12] == seq2[j-12]:
					match_matrix[i][j][0] += -1
				else:
					match_matrix[i][j][1] += -1
			#Criteria for having contamination
			if (match_matrix[i][j][0] >= 10) and (match_matrix[i][j][1] <= 1):
				#useful printlines...
				#print "seq: " + seq1[i-11:i] + " vector: " + seq2[j-11:j] + " mismatches:	 " + str(match_matrix[i][j][1])
				#print "seq: " + str(i) + " vector: " + str(j) + " mismatches:	 " + str(match_matrix[i][j][1])
				return True			
	return False

#contamination(f1,f2): reads from a seq file and vector file. File format must be like the examples presented with the homework.
# Returns a matrix of booleans representing if a given sequence is contaminated with a given vector. sequences are on rows, vectors on columns.
def contamination(seqs,vectors):
	sfile=open(seqs,"r")
	vfile=open(vectors,"r")

	seqs1 = sfile.read().split("\n\n")
	seqs=[]
	for seq in seqs1:
		seqs.append(seq.replace('\n',''))
	vectors1 = vfile.read().split("\n\n")
	vectors=[]
	for vec in vectors1:
		vectors.append(vec.replace('\n',''))

	CONTAMINATION = [[False for i in range(len(seqs))] for j in range(len(vectors))]

	for i in range(len(seqs)):
		print "analyzing all vectors for sequence: " + str(i+1)
		for j in range(len(vectors)):

			CONTAMINATION[j][i] = local_alignment_contamination(seqs[i],vectors[j])

	for row in CONTAMINATION:
		print row

#Function Call
def main():
	contamination(seqs,vectors)

if  __name__ =='__main__':main()