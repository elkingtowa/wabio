from de_bruijn_velvet import *
from simplify import *
from simulate_reads import *
from remove_tips import *
from remove_low_coverage import *
from tour_bus import *
from assemble_paths import *
from resolve_repeats import *

#assemble(reads,k,tour_bus_args,low_coverage_cutoff): Runs all the sub functions necessary to assemble a set of reads into a graph. Takes as input many parameters 
# reads: a set of reads generated or given as input
# k: the k at which to aseemble the de bruijn graph
# tour_bus_args: a list of the arguments to pass to the tour bus algorithm: [minimum merge score, [alignment match score, mismatch score, gap penalty]]
# low_coverage_cutoff: value to pass to remove_low_coverage
# outfile: name of the output file for writing contigs
# write_intermediate: should intermediate results be written to files? useful for debugging. 

def assemble(reads, k, tour_bus_args, low_coverage_cutoff, outfile, write_intermediate):
	print "construct kmers"
	kmers = make_kmers(reads, k)

	print "construct graph"
	G = construct_de_bruijn_velvet(kmers, "False", '')
	if write_intermediate:
		nx.write_gexf(G, 'initial_graph.gexf')

	print "simplify"
	simplify(G, False)
	if write_intermediate:
		nx.write_gexf(G, 'after_simplify_1.gexf')

	print "remove tips"
	remove_tips(G, k)

	print "simplify"
	simplify(G, False)
	if write_intermediate:
		nx.write_gexf(G, 'after_simplify_2.gexf')

	print "resolve repeats"
	resolve_repeats(G)

	print "tour bus"
	run_tour_bus(G, tour_bus_args[0],tour_bus_args[1])
	if write_intermediate:
		nx.write_gexf(G, 'after_tour_bus.gexf')

	print "simplify"
	simplify(G, False)
	if write_intermediate:
		nx.write_gexf(G, 'after_simplify_3.gexf')

	print 'assemble_paths'
	assemble_paths(G, 3)
	if write_intermediate:
		nx.write_gexf(G, 'after_assemble_paths.gexf')

	print "remove low coverage nodes"
	remove_low_coverage(G, low_coverage_cutoff, 2*k)
	if write_intermediate:
		nx.write_gexf(G, 'final_graph.gexf')

	contig_length = [len(a) for a in G.nodes()]
	mean_contig = float(sum(contig_length)) / len(G.nodes())

	### information about assembly 
	print " "
	print " ***  ASSEMBLY INFORMATION  *** "
	print "number of contigs: " + str(len(G.nodes()))
	print "contig lengths: " + str(contig_length)
	print "max contig length: " + str(max(contig_length))
	print "mean contig length: " + str(mean_contig)

	of = open(outfile, 'w')
	for node in G.nodes():
		of.write(node+ '\n')

	return G.nodes()

### TESTING ###
#1000bp genome
#reads1 = simulate_reads('C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\homework\hw5\\fasta\\1000_genome.fasta',12,50,0.01,0,'')
#G = assemble(reads1, 25, [46,[2,-2,-2]],0.1, 'contigs1.txt', False)

