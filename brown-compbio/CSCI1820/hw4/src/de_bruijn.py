#parameters 
import argparse

parser = argparse.ArgumentParser(description="Script to produce a visualization and text output of the De Brujin graph of an input genome by splitting the sequence into fragments of length k.")
parser.add_argument("genomefile", action="store", metavar="genomefile", help="The file name of the fasta sequence containing the genome")
parser.add_argument("k", action="store", metavar="k", help="An integer value k for constructing the graph.")
parser.add_argument("outfile", action="store", metavar="outfile", help="The file name to write results to.")
parser.add_argument("--draw", action="store", metavar="draw", default="True", help="set to False to supress drawing the graph.")
args=parser.parse_args()

genomefile=args.genomefile
k=int(args.k)
outfile = args.outfile
draw = args.draw

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

#de_bruijn(fasta,k): returns a visualization of the De Brujin graph by splitting the sequence in fasta into mers of length k
def de_bruijn(genomefile, k, outfile, draw):
	#read fasta sequence
	f = open(genomefile,'r')
	seq = ''
	lines = f.read().split('\n')
	for line in lines:
		if (len(line) > 0):
			if (line[0] != ">"):
				seq=seq+line
	L = len(seq)
	f.close()

	#split into unique kmers. 
	seq_kmers = []
	for base in range(L-k+1):
		seq_kmers.append(seq[base:base+k])	

	#make k-1mers. don't care if ther are duplicates here. 
	k1mers = [x[:-1] for x in seq_kmers]
	k1mers_array = np.array(k1mers)

	#print seq_kmers
	#print k1mers_array

	#find overlaps
	edge_list = []
	for kmer in seq_kmers:
		matches = np.where(k1mers_array==kmer[1:])
		for match in matches[0]:
			#print match
			edge_list.append((kmer, seq_kmers[match]))

	#print edge_list
	#make graph
	G = nx.DiGraph()
	#add seq_kmers as nodes and overlaps as edges
	G.add_nodes_from(seq_kmers)
	G.add_edges_from(edge_list)

	# draw the graph if desired
	if draw == "True":
		nx.draw_spring(G)
		plt.show()

	#output adjacency list format of the graph
	nx.write_adjlist(G, outfile)

# function call
def main():
	de_bruijn(genomefile, k, outfile, draw)

if  __name__ =='__main__':main()