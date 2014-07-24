import networkx as nx

#assemble_paths(G, score): Assembles subgraphs that have errors preventing them from being simplified
# and cases left over from tour bus. one case in particular, where a start node has two neighbors, one of which
# points to the other neighbor, happend a lot in tour bus and needed to be corrected for. I simply remove the node
# belonging to the longer path in this case.  
# paths with possible errors are merged if they have an overlap greater than score.
def assemble_paths(G, score):
	#find start of subgraph - that's the node with in degree of 0
	in_degrees = G.in_degree()
	out_degrees = G.out_degree()
	starts = [key for (key, value) in in_degrees.items() if value==0]	

	#start exploring form nodes that have in degree of 0
	for start in starts:
		#print "starting from " + start
		previous = start
		can_cont = True
		while G.neighbors(previous) != [] and can_cont:
			next_nodes = G.neighbors(previous)
			if len(next_nodes)==1:
				#if they overlap by more than score, merge them 
				#print "len previous "+ str(len(previous))
				#print "len next_nodes[0] " + str(len(next_nodes[0]))
				if (len(previous)>0) and (len(next_nodes[0])> 0):
					overlap_out = max_overlap(previous, next_nodes[0], float(G.node[previous]['num'])/len(previous),
				 	float(G.node[next_nodes[0]]['num'])/len(next_nodes[0]))
				else:
					print '  somehow we tried to overlap with an empty node'
					break
				if overlap_out[1] > score:
					# we merge the two. 
					G.add_node(overlap_out[0], num=(G.node[previous]['num']+ G.node[next_nodes[0]]['num']))
					for pred_of_prev in G.predecessors(previous):
						G.add_edge(pred_of_prev, overlap_out[0])
					for neigh_of_next in G.neighbors(next_nodes[0]):
						G.add_edge(overlap_out[0], neigh_of_next)
					G.remove_node(previous)
					G.remove_node(next_nodes[0])
					previous=overlap_out[0]
				else:
					previous=next_nodes[0]

			#do the case where a previous points to next but theres an alternative path that passes though one node
			# in this case, just take the shorter path if they can be merged
			elif len(next_nodes) ==2:
				#print "2 next nodes"
				if (G.neighbors(next_nodes[0]) != []):
					if G.neighbors(next_nodes[0])[0] == next_nodes[1]:
						#take path from previous to next_nodes[1]
						G.remove_node(next_nodes[0])
					else:
						can_cont=False
				elif (G.neighbors(next_nodes[1]) != []):
					if G.neighbors(next_nodes[1])[0] == next_nodes[0]:
						#take path from previous to next_nodes[0]
						G.remove_node(next_nodes[1])
					else:
						can_cont=False
				else: 
					can_cont =False
			else: 
				can_cont=False




#max_overlap(seq1, seq2): find the position where the overlap between the two sequences has the highest score.
# seq2 will always come after seq1.
# cov1 and cov2 correspond to the coverage of the two sequences and are used to decide the majority node. 
def max_overlap(seq1, seq2, cov1, cov2):
	min_length = min(len(seq1), len(seq2))
	overlap_scores=[0 for i in range(min_length)]

	#overlap by i bases
	for i in range(min_length):
		s1 = seq1[len(seq1)-i:]
		s2 = seq2[:i]
		overlap_scores[i] = count_match(s1, s2)

	#find the first position in the overlap score array that has the highest score
	pos = overlap_scores.index(max(overlap_scores))
	#use the sequence with higher coverage as majority for overlap
	if cov1 > cov2:
		to_output = seq1 + seq2[pos:]
	else: to_output = seq1[:len(seq1)-pos] +seq2
	return [to_output, max(overlap_scores)]


#count_match(seq1, seq2): helper function for max_overlap - counts the number of matches and mismatches for that overlap
def count_match(seq1, seq2):
	count = 0
	for pos in range(len(seq1)):
		if seq1[pos] == seq2[pos]:
			count += 1
		else:
			count -= 1 
	return count
	


	



