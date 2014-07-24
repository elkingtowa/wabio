#local_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty, return_score): returns the optimal local alignment of two sequences by using the given parameters.
# uses the smith-waterman algorithm for alignment. Depends on helper function single_score. 
# if return_score is true, simply return the highest score for the alignment. otherwise, return a list of the two aligned strings. 
def local_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty, return_score):
	n = len(seq1)
	m = len(seq2)
	score_matrix=[[0 for i in range(m+1)] for j in range(n+1)]
	trace_matrix = [[0 for i in range(m+1)] for j in range(n+1)]
	max_score = -999999999
	max_pos=[0,0]
	for i in range(1, n+1):
		for j in range(1, m+1):
			#arguments to smith-waterman maximization 
			args=[(score_matrix[i-1][j-1] + single_score(seq1[i-1], seq2[j-1], match_score, mismatch_score)),
				(score_matrix[i-1][j] + gap_penalty),
				(score_matrix[i][j-1] + gap_penalty),
				0]

			#pick max
			score_matrix[i][j] = max(args)
			#argmax for traceback
			trace_matrix[i][j] = args.index(max(args))

			#keep track of max so we don't have to find it later
			if max(args) >= max_score:
				max_score=max(args)
				max_pos=[i,j]

	# For LOCAL alignment we want to start at the higest score in the matrix.
	# Traceback
	al1 = ''
	al2 = ''
	current_score= max_score
	cur_i =max_pos[0]
	cur_j =max_pos[1]
	while current_score > 0:
		#if match or mismatch
		if trace_matrix[cur_i][cur_j] == 0:
			al1 += seq1[cur_i-1]
			al2 += seq2[cur_j-1]
			cur_i += -1
			cur_j += -1
			current_score=score_matrix[cur_i][cur_j]
		#if gap in i 
		elif trace_matrix[cur_i][cur_j] == 1:
			al1 += seq1[cur_i-1]
			al2 += "-"
			cur_i += -1
			current_score=score_matrix[cur_i][cur_j]
		#if gap in j
		elif trace_matrix[cur_i][cur_j] == 2:
			al1 += "-"
			al2 += seq2[cur_j-1]
			cur_j += -1
			current_score=score_matrix[cur_i][cur_j]
		#if alignment is discontinued 
		elif trace_matrix[cur_i][cur_j] == 3:
			current_score=0
	if return_score:
		return max_score
	else:
		return[al1[::-1],al2[::-1]]



#single_score(char1, char2): defines the match and mismatch scores for a single base
def single_score(char1, char2, match_score, mismatch_score):
	if char1==char2:
		return match_score
	else:
		return mismatch_score

# Test cases
#print local_alignment("ATG","ATG",1,-1,-1,False)
#print local_alignment("ATGATG","ATGCATG",1,-1,-1,False)
#print local_alignment("MISMATCHINGATGATG","ATGCATG",1,-1,-1,False)
#print local_alignment("MPRSFLVRKPSDPRRKPNYSELQDACVEFTFQQPYDQAHLLAAIPPPEVLNPAASLPTLIWDSLLVPQVRPVAWATLPLRESPKAVELTSLSDEDSGKSSQPPSPPSPASSFSSTSASS","MPRSFLVRKPSDPRRKPNYSELQEVLNPAASLPTLIWDSLLVPQVRPVAWATLPLRESPKAVELTSLSDEDSGKSSQPPSPPSPAPSSFSSTSASS",2,-2,-1,False)
#print local_alignment("MALQGISVVELSGLAPGPFCAMVLADFGARVVRVDRPGSRYDVSRLGRGKRSLVLDLKQPRGAAVLRRLCKRSDVLLEPFRRGVMEKLQLGPEILQRENPRLIYARLSGFGQSGSFCRLAGHDINYLALSGVLSKIGRSGENPYAPLNLLADFAGGGLMCALGIIMALFDRTRTGKGQVIDANMVEGTAYLSSFLWKTQKSSLWEAPRGQNMLDGGAPFYTTYRTADGEFMAVGAIEPQFYELLIKGLGLKSDELPNQMSMDDWPEMKKKFADVFAKKTKAEWCQIFDGTDACVTPVLTFEEVVHHDHNKERGSFITSEEQDVSPRPAPLLLNTPAIPSFKRDPFIGEHTEEILEEFGFSREEIYQLNSDKIIESNKVKASL","MAGPLSGLRVVELAGIGPGPHAAMILGDLGADVVRIDRPSSVDGISRDAMLRNRRIVTADLKSDQGLELALKLIAKADVLIEGYRPGVTERLGLGPEECAKVNDRLIYARMTGWGQTGPRSQQAGHDINYISLNGILHAIGRGDERPVPPLNLVGDFGGGSMFLLVGILAALWERQSSGKGQVVDAAMVDGSSVLIQMMWAMRATGMWTDTRGANMLDGGAPYYDTYECADGRYVAVGAIEPQFYAAMLAGLGLDAAELPPQNDRARWPELRALLTEAFASHDRDHWGAVFANSDACVTPVLAFGEVHNEPHIIERNTFYEANGGWQPMPAPRFSRTASSQPRPPAATIDIEAVLTDWDG",1,-1,-1,False)


