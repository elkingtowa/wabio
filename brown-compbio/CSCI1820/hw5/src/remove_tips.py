
#remove_tips(G, k): continues with error correction by removing tips in G
# tips are a chain that's disconnected at one end. 
# Removed if 1) length less than 2k
#			 2) another neighbor of predecessor has higher coverage -or- another predecessor of neighbor has higher coverage  
def remove_tips(G, k):
	remove_tips = []
	#I imagine this might need to happen more than once? Only stop if we've looked and haven't found any tips
	found_tips = True
	while found_tips:
		#check for tips with in degree 0
		for (tip, degree) in G.out_degree().items():
			#print "tip: " + tip + " degree: " + str(degree)
			if (degree == 0) and (len(tip) < 2*k):
				# must also have 1 predecessor
				if len(G.predecessors(tip)) == 1:
					#get coverage of tip
					tip_num = G.node[tip]['num']
					#get neighbors of predecessor 
					neighbors_of_pred = G.neighbors(G.predecessors(tip)[0])
					#get coverage of neighbors of predecessor
					neighbor_coverage =[]
					for neighbor in neighbors_of_pred:
						neighbor_coverage.append(G.node[neighbor]['num'])
					#if coverage less than another, add to the list to remove
					if max(neighbor_coverage)>tip_num:
						remove_tips.append(tip)
						#print tip + " is a tip!"
		#Do the same but for tips with in degree 0
		for (tip, degree) in G.in_degree().items():
			if (degree == 0) and (len(tip) < 2*k):
				if len(G.neighbors(tip)) == 1:
					#get coverage of tip
					tip_num = G.node[tip]['num']
					#get predecessors of neighbor
					predecessors_of_neighbor = G.predecessors(G.neighbors(tip)[0])
					#get coverage of predecessors of neighbor
					predecessor_coverage = []
					for pred in predecessors_of_neighbor:
						predecessor_coverage.append(G.node[pred]['num'])
					#if coverage less than another, add to the list to remove
					if max(predecessor_coverage)>tip_num:
						remove_tips.append(tip)
						#print tip + " is a tip!"

		#print remove_tips
		if remove_tips == []:
			found_tips = False

		G.remove_nodes_from(remove_tips)

		#nx.draw(G)
		#plt.show()
		remove_tips = []
	return G

