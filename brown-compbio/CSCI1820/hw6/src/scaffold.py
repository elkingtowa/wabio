import argparse
parser = argparse.ArgumentParser(description="Script to assemble contigs and mate pair information into larger scaffolds.")
parser.add_argument("contigfile", action="store", metavar="contigfile", help="The file name of the newline delimited contigs")
parser.add_argument("matefile", action="store", metavar="matefile", help="The file name of the newline delimited mate pairs.")
args=parser.parse_args()
contigfile=args.contigfile
matefile = args.matefile


#extract paths(component): helper function to extract the path through a supercontig
def extract_path(component):
	#find first node
	for node in component.in_degree().items():
		if node[1] ==0:
			start = node[0]

	path=[]
	current = start
	#assume one connection...
	while component.neighbors(current) != []:
		path.append(current)
		current = component.neighbors(current)[0]
	path.append(current)
	return path

#component_match(comp1, comp2, num_contigs): determines if comp1 and comp2 are reverse complement supercontigs of one another 
def component_match(comp1, comp2, num_contigs):
	matching = True
	for i in comp1:
		if i >= num_contigs:
			j = i - num_contigs
		else:
			j = i + num_contigs
		if j not in comp2:
			matching = False
	return matching

#reverse_complement(pattern): helper function for simulate_reads. returns the reverse complement of a DNA sequence
# assumes an alphabet of {a,t,c,g}
def reverse_complement(pattern):
	chars = list(pattern)
	complement = ''
	lookup = dict({("a","t"),("t","a"),("c","g"),("g","c")})
	for base in chars:
		complement += lookup[base]
	rev_complement = complement[::-1]
	return rev_complement 


# read_data(contig_file, mate_file): reads the data from contigs and mate pairs. returns a list of contigs(including reverse complements) and mate pair data 
def read_data(contig_file, mate_file):
	cf = open(contig_file, 'r')
	mf = open(mate_file, 'r')

	contigs = []
	mates = []

	cline = cf.readline()
	while cline != '':
		contigs.append(cline.strip())
		cline = cf.readline()
	#append reverse complement of all contigs 
	for i in range(len(contigs)):
		contigs.append(reverse_complement(contigs[i]))

	mline=mf.readline()
	while mline != '':
		mline = mline.strip().split('\t')
		mates.append([mline[0],mline[1],mline[2]])
		mline=mf.readline()

	mf.close()
	cf.close()

	return contigs, mates


from local_alignment import *
import math
#OLD CODE. using local alignment in python was too damn slow. Would produce the same results as what I have below but take 1000x as long
def align_mates2(contigs, mates):

	mate_matches = dict()
	for i in range(len(mates)):
		mate_matches[i] = [[],[]]

	for i in range(len(mates)):
		mate1 = mates[i][0]
		mate2 = mates[i][2]
		threshold_score = len(mate1) - (int(math.ceil(len(mate1)*0.02))*2)
		for j in range(len(contigs)):
			contig = contigs[j]
			alignment_score1 = local_alignment(contig, mate1, 1, -1, -99999999, True)
			alignment_score2 = local_alignment(contig, mate2, 1, -1, -99999999, True)

			match1 = alignment_score1[0] >= threshold_score
			match2 = alignment_score2[0] >= threshold_score

			if not(match1 and match2):
				if match1:
					mate_matches[i][0].append([j, alignment_score1[1][0]])
					print "mate " + str(i) + " left matched in contig " + str(j)
				elif match2:
					mate_matches[i][1].append([j, alignment_score2[1][0]])
					print "mate " + str(i) + " right matched in contig " + str(j)
	return mate_matches

#how many exact matches
import regex
import math
import networkx as nx
def create_scaffold(contigs, mates):
	exact_matches1 = [[] for a in range(len(mates))]
	exact_matches2 = [[] for a in range(len(mates))]
	d= 0
	for i in range(len(mates)):
		mate = mates[i]
		mate1=mate[0]
		mate2=mate[2]
		#print d
		d+=1
		for j in range(len(contigs)):
			contig = contigs[j]
			allowed_mismatches = int(math.ceil(len(mate1)*0.02))
			#match using a regex with a certain number of allowed mismatches
			regex1 = '(?:' + mate1 + '){s<='+str(allowed_mismatches)+'}'
			regex2 = '(?:' + mate2 + '){s<='+str(allowed_mismatches)+'}'

			a = regex.search(regex1, contig)
			b = regex.search(regex2, contig)
			if not (a != None and b !=None):
				if a != None:
					exact_matches1[i] = [j, a.start()]
				if b != None:
					exact_matches2[i] = [j, b.start()]

	# find set where each mate maps to a contig
	initial_set = []
	for j in range(len(exact_matches1)):
	    if (exact_matches1[j] != [] and exact_matches2[j] != []):
	    	initial_set.append([mates[j], exact_matches1[j], exact_matches2[j]])

	#compute a matrix with the number of interactions between two contigs
	contig_mat = [[0 for l in range(len(contigs))] for l in range(len(contigs))]
	for i in range(len(initial_set)):
		y=initial_set[i][1][0]
		x=initial_set[i][2][0]
		contig_mat[x][y] += 1

	num_contigs = len(contigs)/2
	G = nx.DiGraph()
	#create a graph to represent the paths through the contigs
	for i in range(len(contigs)):
		for j in range(len(contigs)):
			num = contig_mat[i][j]
			if num >=2:
				#check to make sure we're not going to create a small loop
				if (j,i) in G.edges():
					 num1=G.edge[j][i]['num']
					 if num > num1:
					 	# if we are going to make a loop, but it can be improved, remove the old edge
					 	G.add_edge(i,j,num=num)
					 	G.remove_edge(j,i)
				else:
					G.add_edge(i, j, num=num)

	#components of graph are the paths of the assembly
	component_graphs = nx.weakly_connected_component_subgraphs(G)
	components = nx.weakly_connected_components(G)
	#each will be repeated with normal and reverse orientation. 
	component_matches = []
	for i in range(len(components)):
		comp1 = components[i]
		for j in range(len(components)):
			comp2 = components[j]
			if component_match(comp1, comp2, num_contigs):
				component_matches.append([i,j])

	#pick unique components
	unique_components = []
	for i in component_matches:
		if i[::-1] not in unique_components:
			unique_components.append(i)

	#fid paths for unique components
	unique_component_paths= []
	for i in unique_components:
		unique_component_paths.append([extract_path(component_graphs[i[0]]), extract_path(component_graphs[i[1]])])

	#report just one of the unique component paths
	paths = [i[0] for i in unique_component_paths]
	#change numbers greater than num_contigs to negative representation
	# simply for nicer output
	for i in range(len(paths)):
		for j in range(len(paths[i])):
			if paths[i][j] >= num_contigs:
				paths[i][j] = (paths[i][j]-num_contigs) * -1

	#return the paths as our supercontigs!
	return paths

#function to wrap everything together
def do_scaffold(contigfile, matefile):
	data =read_data(contigfile, matefile)
	paths = create_scaffold(data[0],data[1])

	#print results
	print "SCAFFOLD RESULTS"
	for path in paths:
		line = ''
		for pos in path:
			line += str(pos) + "  ->  "
		print line[:-4]

#function call
def main():
	do_scaffold(contigfile,matefile)

if  __name__ =='__main__':main()