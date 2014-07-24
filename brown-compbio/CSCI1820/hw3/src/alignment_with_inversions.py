#parameters 
import argparse

parser = argparse.ArgumentParser(description="Script compute the optimal local alignment with inversions of two sequences. Takes input from two fasta sequences.")
parser.add_argument("fasta1", action="store", metavar="fasta1", help="The file name of the first fasta sequence")
parser.add_argument("fasta2", action="store", metavar="fasta2", help="The file name of the second fasta sequence")
parser.add_argument("--match", action="store", metavar="match", default="1", help="score for a match")
parser.add_argument("--mismatch", action="store", metavar="mismatch", default="-1", help="score for a mismatch")
parser.add_argument("--gap", action="store", metavar="gap", default="-1", help="score for a gap")
parser.add_argument("--invert", action="store", metavar="invert", default="-1", help="score for a invert")
args=parser.parse_args()

#parameters and paths specified in section above
fasta1=args.fasta1
fasta2=args.fasta2
match=int(args.match)
mismatch=int(args.mismatch)
gap=int(args.gap)
invert=int(args.invert)



#global_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty, return_score): Helper function for local alignment with inversions. 
# returns the optimal global alignment of two sequences by using the given parameters.
# uses the needleman-wunsch algorithm for alignment. Depends on helper function single_score. 
# if return_score is true, simply return the highest score for the alignment. otherwise, return a list of the two aligned strings. 
def global_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty, return_score):
	n = len(seq1)
	m = len(seq2)
	score_matrix=[[0 for i in range(m+1)] for j in range(n+1)]
	trace_matrix = [[0 for i in range(m+1)] for j in range(n+1)]
	for i in range(1, n+1):
		for j in range(1, m+1):
			#arguments to maximization 
			args=[(score_matrix[i-1][j-1] + single_score(seq1[i-1], seq2[j-1], match_score, mismatch_score)),
				(score_matrix[i-1][j] + gap_penalty),
				(score_matrix[i][j-1] + gap_penalty)]

			#pick max
			score_matrix[i][j] = max(args)
			#argmax for traceback
			trace_matrix[i][j] = args.index(max(args))
	# For GLOBAL alignment we want to start at the LAST score in the matrix.
	# Traceback
	al1 = ''
	al2 = ''
	cur_i = n
	cur_j = m

	#this code is for traceback, but not important for finding the score
	while (cur_i!=0) or (cur_j!=0):
		# if at one end of the matrix
		if cur_i ==0:
			al1 += "-"
			al2 += seq2[cur_j-1]
			cur_j+=-1
		# if at one end of the matrix
		elif cur_j ==0:
			al2 += "-"
			al1 += seq1[cur_i-1]
			cur_i+=-1
		#if match or mismatch
		elif trace_matrix[cur_i][cur_j] == 0:
			al1 += seq1[cur_i-1]
			al2 += seq2[cur_j-1]
			cur_i += -1
			cur_j += -1
			
		#if gap in i 
		elif trace_matrix[cur_i][cur_j] == 1:
			al1 += seq1[cur_i-1]
			al2 += "-"
			cur_i += -1
		#if gap in j
		elif trace_matrix[cur_i][cur_j] == 2:
			al1 += "-"
			al2 += seq2[cur_j-1]
			cur_j += -1

	#only return score of global alignment
	if return_score:
		return score_matrix[n][m]
	else:
		#return the actual alignment 
		return[al1[::-1],al2[::-1]]

#single_score(char1, char2): defines the match and mismatch scores for a single base
def single_score(char1, char2, match_score, mismatch_score):
	if char1==char2:
		return match_score
	else:
		return mismatch_score

# reverse_complement(pattern): returns the reverse complement of DNA string pattern. Equivalent to an inversion. 
def reverse_complement(pattern):
	chars = list(pattern)
	complement = ''
	lookup = dict({("A","T"),("T","A"),("C","G"),("G","C")})
	for base in chars:
		complement += lookup[base]
	rev_complement = complement[::-1]
	return rev_complement 

#local_alignment_inversions(seq1,seq2,match_score, mismatch_score, gap_penalty,inversion_penalty): performs local alignment with inversions 
# uses the algorithm from the waterman paper. 
# prints the alignment with "*" highlighting the region to be reverse complemented. 
def local_alignment_inversions(seq1,seq2,match_score, mismatch_score, gap_penalty,inversion_penalty):
	n = len(seq1)
	m = len(seq2)
	#Initialize matricies for scoring, traceback, and Z
	score_matrix=[[0 for i in range(m+1)] for j in range(n+1)]
	trace_matrix = [[[-1,0] for i in range(m+1)] for j in range(n+1)]
	Z = [[[[inversion_penalty for j in range(m+1)] for i in range(n+1)] for h in range (m+1)] for g in range (n+1)]
	max_score = -999999999
	max_pos=[0,0]

	for i in range(1, n+1):
		for j in range(1, m+1):
			for g in range(i, 0,-1):
				for h in range(j, 0, -1):
					#really only need to do inversions when sequences are the same length
					if (i-g == j-h):
						#Calculate Zscore by doing global alignment of the sequence and the reverse complement
						Z[g-1][h-1][i-1][j-1] = global_alignment(seq1[g-1:i],reverse_complement(seq2[h-1:j]),match_score,mismatch_score,gap_penalty,True) + inversion_penalty

			#get Z values for i and j. keep length of inversion
			zvals = []
			lengths = []
			for g in range(1, i+1):
				for h in range(1, j+1):
					zvals.append(score_matrix[g-2][h-2] + Z[g-1][h-1][i-1][j-1])
					lengths.append(i-g+1)

			#store values for argmax implementation 
			args = [
				max(zvals),
				(score_matrix[i-1][j-1] + single_score(seq1[i-1],seq2[j-1],match_score,mismatch_score)),
				(score_matrix[i-1][j] + gap_penalty),
				(score_matrix[i][j-1] + gap_penalty),
				0]

			#pick max
			score_matrix[i][j] = max(args)
			#argmax for traceback
			if args.index(max(args)) ==0:
				#if inversion, store the length
				trace_matrix[i][j] = [0,lengths[zvals.index(max(zvals))]]
			else: trace_matrix[i][j] = [args.index(max(args)),0]

			#keep track of max score and max position so we don't have to find it later
			if max(args) >= max_score:
				max_score=max(args)
				max_pos=[i,j]

	#TRACEBACK
	#sequences
	al1 = ''
	al2 = ''
	#current positions
	current_score= max_score
	cur_i =max_pos[0]
	cur_j =max_pos[1]

	#testing prints.....
	for row in score_matrix:
		print row
	print " --- "
	for row in trace_matrix:
		print row
	print " ****** "

	while current_score > 0:
		# if we're in an inversion
		if trace_matrix[cur_i][cur_j][0] == 0:
			# get length of inversion. Append sequence of inversion with markers around what to reverse complement
			invlength = trace_matrix[cur_i][cur_j][1]
			s1 = seq1[cur_i-invlength:cur_j]
			al1 += " " + s1[::-1] + " "
			s2 = seq2[cur_j-invlength:cur_j]
			al2 += "*" + s2[::-1] + "*"
			cur_i += -invlength
			cur_j += -invlength
			current_score=score_matrix[cur_i][cur_j]
		#if match or mismatch
		if trace_matrix[cur_i][cur_j][0] == 1:
			al1 += seq1[cur_i-1]
			al2 += seq2[cur_j-1]
			cur_i += -1
			cur_j += -1
			current_score=score_matrix[cur_i][cur_j]
		#if gap in i 
		elif trace_matrix[cur_i][cur_j][0] == 2:
			al1 += seq1[cur_i-1]
			al2 += "-"
			cur_i += -1
			current_score=score_matrix[cur_i][cur_j]
		#if gap in j
		elif trace_matrix[cur_i][cur_j][0] == 3:
			al1 += "-"
			al2 += seq2[cur_j-1]
			cur_j += -1
			current_score=score_matrix[cur_i][cur_j]
		#if alignment is discontinued 
		elif trace_matrix[cur_i][cur_j][0] == 4:
			current_score=0


	print al1[::-1]
	print al2[::-1]

	#don't really need to return anything
	#return[al1[::-1],al2[::-1]]

#helper function to read fastas and call alignment with inversions
def do_alignment(fasta1, fasta2,match_score, mismatch_score, gap_penalty,inversion_penalty):
	f1=open(fasta1, 'r')
	f2=open(fasta2, 'r')
	seq1 = ''
	seq2 = '' 

	#Read fasta files into a sequence
	lines1 = f1.read().split('\n')
	for line in lines1:
		if (len(line) > 0):
			if (line[0] != ">"):
				seq1=seq1+line

	lines2 = f2.read().split('\n')
	for line in lines2:
		if (len(line) > 0):
			if (line[0] != ">"):
				seq2=seq2+line

	local_alignment_inversions(seq1,seq2,match_score, mismatch_score, gap_penalty,inversion_penalty)

#FUNCTION CALL
def main():
	do_alignment(fasta1,fasta2,match,mismatch,gap,invert)

if  __name__ =='__main__':main()