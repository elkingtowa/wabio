import networkx as nx
import matplotlib.pyplot as plt

from de_bruijn_velvet import overlap

#simplify(G): takes an input graph and merges nodes that are continuous. if draw is True, draw the output graph after simplification
# node coverage is the sum of coverage of merged nodes
# node sequence is the overlap between nodes
def simplify(G, draw):
	#can we simplify? Is there a node with out degree 1 who's neighbor has an in degree of 1?
	can_simplify = True
	have_simplified = False
	while can_simplify:
		have_simplified=False
		for (node, out) in G.out_degree().items():
			if not have_simplified:
				if (out ==1) and (G.in_degree(G.neighbors(node)).values() ==[1]):
					#check to make sure the pair isnt a loop, which we don't want to simplify	
					if G.neighbors(node)[0] not in G.predecessors(node):
						have_simplified = True

						head_num = G.node[node]['num']
						tail_seq = G.neighbors(node)[0]
						tail_num = G.node[tail_seq]['num']

						new_seq = overlap(node, tail_seq)
						G.add_node(new_seq, num=head_num+tail_num)

						for pred in G.predecessors(node):
							G.add_edge(pred,new_seq)
						for neighbor in G.neighbors(tail_seq):
							G.add_edge(new_seq, neighbor)

						G.remove_node(node)
						G.remove_node(tail_seq)

		if not have_simplified:
			print "Simplification complete"
			can_simplify=False
	if draw:
		nx.draw_spring(G)
		plt.show() 
		
	return G

