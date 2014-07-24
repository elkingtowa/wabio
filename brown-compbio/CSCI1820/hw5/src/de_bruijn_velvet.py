import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

#construct_de_bruijn_velvet(kmers, draw, write): constructs the De Bruijn graph for the sequences in kmers. 
# kmers is a dictionary with keys being the sequences and values being the number of time the kmer is represented
# all keys in kmers must be of length k. All values must be integers >0
# If draw is "True", draw a diagram of the graph to the screen. Can be messy and take a long time if the graph is large!
# If outfile is a nonempty string, output an adjacency list format of the graph to the specified file.
def construct_de_bruijn_velvet(kmers, draw, outfile):
	#make list of k-1mers for quick edge construction
	k1mers = [x[:-1] for x in kmers.keys()]
	k1mers_array = np.array(k1mers)

	#find overlaps
	edge_list = []
	for kmer in kmers.keys():
		matches = np.where(k1mers_array==kmer[1:])
		for match in matches[0]:
			#print match
			edge_list.append((kmer, kmers.keys()[match]))

	#make graph
	G = nx.DiGraph()
	#add seq_kmers as nodes and overlaps as edges
	for kmer in kmers.items():
		G.add_node(kmer[0], num=kmer[1])
	G.add_edges_from(edge_list)

	# draw the graph if desired
	if draw == "True":
		nx.draw_spring(G)
		plt.show()

	#output adjacency list format of the graph if desired
	if outfile != "":
		nx.write_adjlist(G, outfile)

	return G


#make_kmers(reads, k): constructs a dictionary of kmers from an input of reads
def make_kmers(reads,k):
	kmer_dict = {}
	for read in reads:
		#throw out reads less then length k
		if len(read)>=k:
			for start in range(len(read)-k+1):
				kmer = read[start:start+k]
				#count number of occurances
				if kmer in kmer_dict:
					kmer_dict[kmer] += 1
				else:
					kmer_dict[kmer] = 1
	return kmer_dict

#helper function for finding new titles of nodes
def overlap(s1, s2):
	#how many bases of overlap do we have?
	min_len = min(len(s1), len(s2))
	for i in range(min_len,0,-1):
		if s1[len(s1)-i:] == s2[:i]:
			return s1 + s2[i:]
	#didnt find an overlap?
	return ''
