import networkx as nx
from de_bruijn_velvet import *

#resolve_repeats(G): attempts to find and resolve repetitive sequences in the graph G
# works by identifying cycles
# then finds input and output of those cycles and condenses the repeated sequence
# repeats the sequence according to a strict division of the average coverage in the repeat by the 
# average coverage of the incoming and outgoing nodes
def resolve_repeats(G): 
	#does the graph have any cycles?
	cycles = nx.simple_cycles(G)
	if len(list(cycles)) > 0:
		for cycle in cycles:
			#find if there is an input or output from the cycle
			inputs = []
			outputs = []
			for node in cycle:
				if node in G.nodes():
					node_pred = [a for a in G.predecessors(node) if a not in cycle]
					if len(node_pred) != 0:
						for n in node_pred:
							inputs.append(n)
					node_neigh = [a for a in G.neighbors(node) if a not in cycle]
					if len(node_neigh) != 0:
						for n in node_neigh:
							outputs.append(n)
				else: break
			#compute average coverage in cycle
			total_cov = 0
			for node in cycle[1:]:
				total_cov += G.node[node]['num']
			av_cov = float(total_cov)/float(len(cycle[1:]))

			#get average coverage of input an output nodes
			total_outside_cov = 0
			for node in inputs:
				total_outside_cov += G.node[node]['num']
			for node in outputs: 
				total_outside_cov += G.node[node]['num']
			av_outside_cov = float(total_outside_cov)/float(len(inputs)+len(outputs))

			#we should repeat the cycle as many times as the coverage inside the cycle divided by the coverage outside the cycle
			repeat_num = int(round(av_cov/av_outside_cov))
			#extract sequence from repeat
			# is it just the second item in the list because it's been compressed
			sequence = cycle[1]
			#repeat this repeat_num times
			final_sequence = ''
			for i in range(repeat_num):
				final_sequence += sequence
			if repeat_num > 0:
				final_coverage = av_cov / repeat_num
			else:
				final_coverage=1
			#add node, edges to inside and outside nodes
			G.add_node(final_sequence, num=int(round(final_coverage)))
			for output_seq in outputs:
				G.add_edge(final_sequence, output_seq)
			for input_seq in inputs:
				G.add_edge(input_seq, final_sequence)

			#delete nodes in cycle
			for node in cycle:
				if node in G.nodes():
					G.remove_node(node)
	return G



###### TESTING ######
# repeat longer than k length
#a = make_kmers(['ZABCDEFGABCDEFGABCDEFGY'], 4)
#G1 = construct_de_bruijn_velvet(a, "True", '')
# repeat shorter then k length
#b = make_kmers(['ZABCDEFGABCDEFGABCDEFGY'], 9)
#G2 = construct_de_bruijn_velvet(b, "False", '')
#simplify(G2, False)
#G21 = resolve_repeats(G)
