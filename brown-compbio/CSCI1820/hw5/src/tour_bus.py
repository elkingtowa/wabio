import networkx as nx
import matplotlib.pyplot as plt

from local_alignment import *
from de_bruijn_velvet import *
from simplify import *
from simulate_reads import *
from remove_tips import *

#tour_bus(G, score, alignment_params): conducts the tour bus algorithm on a graph to remove bubbles caused by sequencing errors
# G: a similified graph passed in from another step of the algorithm
# score: what is the minimum local alignment score that needs to be achieved for two sequences to be merged?
# alignment_params: what parameters should be passed to the alignment algorithm
#  needs to be a list of length 3 [match score, mismatch socre, gap penalty]
def tour_bus(G, score, alignment_params):
	#print 'starting tour '
	for start_node in G.nodes():
		if start_node in G.nodes():
			#print "starting tour bus from " + start_node
			#bfs exploration of the graph. Visit nodes in increasing distance from the origin (the first node in the graph for now)
			# Need to make sure we cover all possible nodes. 
			start = start_node
			previous = start
			distances = {}
			distances[start] = 0.0
			visited = [start]
			#have we found a bubble?
			have_merged=False
			while not have_merged:
				#unveil neighbors of previous 
				new_nodes = [a for a in G.neighbors(previous) if a != previous]
				for i in new_nodes:
					if i in visited:
						#print i + " visited before!"
						#found a node that we've been to before!
						# backtrack and find closest common ancestor. both predecessors of i will be in the visited set. will always have 2 elements.
						to_traceback = [a for a in G.predecessors(i) if a in visited]
						
						#yay for literally the worst way to catch an error
						if len(to_traceback) <2:
							break

						#traceback path 1. shouldn't encounter a node that has more than one predecessor in the visited set
						traceback_path_1 = []
						traceback_pred_1 = [to_traceback[0]]
						while traceback_pred_1 != []:
							traceback_path_1.append(traceback_pred_1[0])
							#print G.predecessors(traceback_pred_1[0])
							traceback_pred_1 = [a for a in G.predecessors(traceback_pred_1[0]) if (a in visited) and (a not in traceback_path_1)]

						#traceback path 2. 
						traceback_path_2 = []
						traceback_pred_2 = [to_traceback[1]]
						while traceback_pred_2 != []:
							traceback_path_2.append(traceback_pred_2[0])
							traceback_pred_2 = [a for a in G.predecessors(traceback_pred_2[0]) if (a in visited) and (a not in traceback_path_2)]
						#find first element that occurs in both
						#print "tpath1: " + str(traceback_path_1) + "   tpath2:" + str(traceback_path_2)
						
						overlap_first = [j for j in traceback_path_1 if j in traceback_path_2][0]

						#find first overlap between lists
						overlap_pos1 = traceback_path_1.index(overlap_first)
						overlap_pos2 = traceback_path_2.index(overlap_first)

						#get rid of path past the overlap
						traceback_path_1 = traceback_path_1[:overlap_pos1]
						traceback_path_2 = traceback_path_2[:overlap_pos2]
						forward_path_1 = traceback_path_1[::-1]
						#print "forward_path_1: " + str(forward_path_1)
						forward_path_2 = traceback_path_2[::-1]
						#print "forward_path_2: " + str(forward_path_2)

						#extract sequences from corresponding paths
						sequence_1 = '' 
						for j in forward_path_1:
							if sequence_1 != '':
								sequence_1 = overlap(sequence_1,j)
							else: 
								sequence_1 = j
						#print sequence_1
						sequence_2 = '' 
						for j in forward_path_2:
							if sequence_2 != '':
								sequence_2 = overlap(sequence_2, j)
							else: 
								sequence_2 = j
						#print sequence_2

						#sequence to continue with reached end node first.
						#that's the path that doesn't have 'previous' in it
						#print "Trying to find what sequence to merge. looking for " + previous
						if previous in forward_path_1:
							seq_to_merge = 2
						elif previous in forward_path_2:
							seq_to_merge =1 
						else: 
							print "should never get here!"
							return False

						#align the two sequences. use local alignment and a strict score
						#could later convert this to a relative percentage match, something 
						# like the 80% they use in velvet. 
						merge = False
						alignment_score = local_alignment(sequence_1,sequence_2, alignment_params[0], alignment_params[1],alignment_params[2],True)
						if alignment_score >  score:
							merge=True

						#Should the two sequences be merged?
						if merge:
							have_merged=True
							#do merging!
							if seq_to_merge == 1:
								#path 1 is being represented in the merged graph. start at the beginning, look at the corresponding node in the other path
								# align and decide to merge. then copy coverage information (addative) and edges	
								j = 0 
								while j < len(forward_path_1):
									k=j
									while k < len(forward_path_2):
										if (local_alignment(forward_path_1[j],forward_path_2[k],alignment_params[0], alignment_params[1],alignment_params[2],True) > score):
											#merge these and copy info
											G.node[forward_path_1[j]]['num'] += G.node[forward_path_2[k]]['num']
											#copy (incoming edges to forward_path_2[k]) to forward_path_1[j]
											for l in G.predecessors(forward_path_2[k]):
												G.add_edge(l,forward_path_1[j])
											#copy (outgoing edges from forward_path_2[k]) to forward_path_1[j]
											for l in G.neighbors(forward_path_2[k]):
												G.add_edge(forward_path_1[j],l)
											#delete forward_path_2[k]
											G.remove_node(forward_path_2[k])
											j=k
											k = len(forward_path_2)
										else: k+=1
									j+=1
							if seq_to_merge == 2:
								#path 2 is being represented in the merged graph. start at the beginning, look at the corresponding node in the other path
								# align and decide to merge. then copy coverage information (addative) and edges
								j = 0 
								while j < len(forward_path_2):
									k=j
									while k < len(forward_path_1):
										if (local_alignment(forward_path_2[j],forward_path_1[k],alignment_params[0], alignment_params[1],alignment_params[2],True) > score):
											#merge these and copy info
											G.node[forward_path_2[j]]['num'] += G.node[forward_path_1[k]]['num']
											#copy (incoming edges to forward_path_1[k]) to forward_path_2[j]
											for l in G.predecessors(forward_path_1[k]):
												G.add_edge(l,forward_path_2[j])
											#copy (outgoing edges from forward_path_1[k]) to forward_path_2[j]
											for l in G.neighbors(forward_path_1[k]):
												G.add_edge(forward_path_2[j],l)
											#delete forward_path_1[k]
											G.remove_node(forward_path_1[k])
											j=k
											k = len(forward_path_1)
										else: k+=1
									j+=1
							return True
					else: 
						#calculate distance to new nodes
						distances[i] = distances[previous] + (float(len(i)) / G.node[i]['num'])
						visited.append(i)
				
				if not have_merged:
					if len(distances)>1:
						del distances[previous]
						#find node with minimum distance
						previous = min(distances, key=distances.get)
					else: 
						can_merge = False
						have_merged=True
	return False

# wrapper function to run tour_bus multiple times
def run_tour_bus(G, score, alignment_params):
	run_tour = True
	while run_tour:
		run_tour = tour_bus(G, score, alignment_params)



#### TEST CASES ####
#G9 = construct_de_bruijn_velvet(make_kmers(
#	['ABCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEFGHIJKLM',
#	'ABCDEFGOIJKLM','ABCDEFGHRJKLM'],4),"True",'')
#simplify(G9, True)
#run_tour_bus(G9, 5,[2,-2,-2])
