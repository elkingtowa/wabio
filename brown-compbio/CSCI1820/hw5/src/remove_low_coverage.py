#remove_low_coverage(G, cutoff,length): removes nodes in G with a coverage in the bottom percentrile of coverage or length less than the specified value.
# a cutoff value of 0.1 means nodes with coverage less than 10% of the max will be removed. 
# a length of 40 means nodes with length less than 40 will be removed. 
def remove_low_coverage(G, cutoff, length):
	degrees = G.degree()
	isolated_nodes = [p[0] for p in degrees.items() if degrees[p[0]]==0]
	#print isolated_nodes

	coverage = [a[1]['num'] for a in G.nodes(data=True)]
	max_coverage = max(coverage)

	for node in G.nodes():
		if G.node[node]['num'] < float(max_coverage) * float(cutoff):
			G.remove_node(node)
	for node in G.nodes():
		if len(node) < length:
			G.remove_node(node)
			