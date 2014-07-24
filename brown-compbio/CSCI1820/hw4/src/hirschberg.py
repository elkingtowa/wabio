# Code was never used :(

#hirschberg(X,Y): returns the optimal global alignment of sequences X and Y using hirschberg's algorithm. 
def hirschberg(X,Y, match_score, mismatch_score, gap_penalty):
	Z=""
	W=""
	if (len(X) == 0) or (len(Y) == 0):
		if len(X) == 0:
			for i in range(len(Y)):
				Z = Z + '-'
				W = W + Y[i]
		elif len(Y) == 0:
			for i in range(len(X)):
				Z = Z + X[i]
				W = W + '-'
	elif (len(X) == 1) or (len(Y) == 1):
		temp = global_alignment(X,Y,match_score, mismatch_score, gap_penalty, False)
		Z = temp[0]
		W = temp[1]
	else:
		xmid = len(X)/2
		print "xmid: " + str(xmid)
		scoreL = nwscore(X[:xmid],Y,match_score, mismatch_score, gap_penalty)
		scoreR = nwscore(X[xmid:][::-1],Y[::-1],match_score, mismatch_score, gap_penalty)
		print "scoreL: " + str(scoreL)
		print "scoreR: " + str(scoreR)

		scoreR.reverse()
		for s in range(len(scoreL)):
			scoreR[s] += scoreL[s]

		print "new scoreR: " + str(scoreR)
		ymid = 	scoreR.index(max(scoreR))
		print "ymid: " + str(ymid)
		H1 = hirschberg(X[:xmid],Y[:ymid],match_score, mismatch_score, gap_penalty)
		H2 = hirschberg(X[xmid:],Y[ymid:],match_score, mismatch_score, gap_penalty)
		Z = H1[0] + H2[0]
		W = H1[1] + H2[1]
	return [Z,W]

#nwscore(seq1, seq2, match_score, mismatch_score, gap_penalty): returns the optimal local alignment of two sequences by using the given parameters.
# uses the needleman-wunsch algorithm for alignment. Depends on helper function single_score. 
def nwscore(seq1, seq2, match_score, mismatch_score, gap_penalty):
	n = len(seq1)
	m = len(seq2)
	score_matrix=[[0 for i in range(m+1)] for j in range(n+1)]
	for i in range(n+1):
		score_matrix[i][0] = gap_penalty *i
	for j in range(m+1):
		score_matrix[0][j] = gap_penalty *j

	for i in range(1, n+1):
		for j in range(1, m+1):
			#arguments to smith-waterman maximization 
			args=[(score_matrix[i-1][j-1] + single_score(seq1[i-1], seq2[j-1], match_score, mismatch_score)),
				(score_matrix[i-1][j] + gap_penalty),
				(score_matrix[i][j-1] + gap_penalty)]

			#pick max
			score_matrix[i][j] = max(args)

	return score_matrix[len(seq1)]

#global_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty, return_score): Helper function for local alignment with inversions. 
# returns the optimal global alignment of two sequences by using the given parameters.
# uses the needleman-wunsch algorithm for alignment. Depends on helper function single_score. 
# if return_score is true, simply return the highest score for the alignment. otherwise, return a list of the two aligned strings. 
def global_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty, return_score):
	n = len(seq1)
	m = len(seq2)
	score_matrix=[[0 for i in range(m+1)] for j in range(n+1)]
	trace_matrix = [[0 for i in range(m+1)] for j in range(n+1)]
	for i in range(n+1):
		score_matrix[i][0] = gap_penalty *i
	for j in range(m+1):
		score_matrix[0][j] = gap_penalty *j
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

#print hirschberg("ATGC","ATGC",1,-1,-1)
#print hirschberg("AGTACGCA","TATGC",2,-1,-2)
#print hirschberg("TATGC","AGTACGCA",2,-1,-2)

#helper function to read fastas and call alignment with inversions
def do_alignment(fasta1, fasta2,match_score, mismatch_score, gap_penalty):
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
	print "read sequences"
	hirschberg(seq1,seq2,match_score, mismatch_score, gap_penalty)

# do_alignment(
# 	"C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\homework\hw4\seqs\APOL3.fasta",
# 	"C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\homework\hw4\seqs\APOL4.fasta",
# 	1,-1,-1)