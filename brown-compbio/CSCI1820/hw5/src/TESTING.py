# TESTING
# test each part of the algorithm separately
# small genome without errors - high coverage. test simplification. 
b = simulate_reads('C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\homework\hw5\\fasta\small_genome.fasta',100,10,0,0)
G1 = construct_de_bruijn_velvet(make_kmers(b,5))
G11 = simplify(G1.copy())

## test some tip removal - tip at end
G4 = construct_de_bruijn_velvet(make_kmers(['ABCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEFGHIJZLM'],4),"True",'')
G41 = simplify(G4.copy())
G42 = remove_tips(G41.copy(),4)
## test some tip removal - tip at beginning
G5 = construct_de_bruijn_velvet(make_kmers(['AZCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEFGHIJKLM'],4),"True",'')
G51 = simplify(G5.copy())
G52 = remove_tips(G51.copy(),4)

#testing on first 1000 base pairs from genome, normal coverage. No errors because I haven't coded that yet.
# seems to produce good contigs (isolated nodes at the end of simplification at this point)
d = simulate_reads('C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\homework\hw5\\fasta\\1000_genome.fasta',8,30,0,0)
G6 = construct_de_bruijn_velvet(make_kmers(d,15),"False",'')
G61 = simplify(G6.copy())
G62 = remove_tips(G61.copy(),15)
G63 = simplify(G62.copy())
contigs6 = G63.nodes()

#testing on small genome, normal coverage. No errors because I haven't coded that yet.
# seems to produce good contigs (isolated nodes at the end of simplification at this point)
e = simulate_reads('C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\homework\hw5\\fasta\\SigmaKappaOmega9_genome_reduced.fasta',8,50,0,0)
G7 = construct_de_bruijn_velvet(make_kmers(e,25),"False",'')
print 1
G71 = simplify(G7.copy())
print 2
G72 = remove_tips(G71.copy(),25)
print 3
G73 = simplify(G72.copy())
contigs6 = G73.nodes()
of = open('C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\homework\hw5\contig_outfile.txt','w')
for contig in contigs6:
	of.write(contig + '\n')
of.close()

# TEST TOUR BUS 
#single sequencing error example -produces a single bubble in the graph 
G8 = construct_de_bruijn_velvet(make_kmers(['ABCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEZGHIJKLM'],4),"True",'')
simplify(G8, True)
run_tour_bus(G8, 5,[2,-2,-2])

#multiple sequencing errors - this actually produces two bubbles which can't be resolved to a single node with my algorithm
G9 = construct_de_bruijn_velvet(make_kmers(
	['ABCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEZGHIJKLM','ABCDEFGHRJKLM'],4),"True",'')
simplify(G9, True)
run_tour_bus(G9, 5,[2,-2,-2])

G2 = construct_de_bruijn_velvet(make_kmers(['ABCDEFGHIJKLM','ABCDEFGHIJKLM','ABCDEZGHIJKLM'],4),"True",'')
#simplify(G2)
run_tour_bus(G2,2,[2,-2,-2])
simplify(G2)
nx.draw(G2)
plt.show()

#example with errors 
b = simulate_reads('C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\homework\hw5\\fasta\\100_genome.fasta',10,10,0.01,0)
G3 = construct_de_bruijn_velvet(make_kmers(b,8),"True", '')
simplify(G3, True)
remove_tips(G3,8)
nx.draw(G3)
plt.show()
simplify(G3, True)
run_tour_bus(G3, 10,[2,-2,-2])


#REQUIREMENTS FOR HOMEWORK HANDIN
# Due Thursday March 20 is (1) completed code; (2) a set freads, assembly, de Bruijn graph,
# and contig mapg for one run of the assembler on the first 10000bp of artificial bacterial genome
# with settings a = 12, r = 0, l = 50, e = :01; (3) a set fassembly and contig mapg for a run of the
# assembler on the reads we provide in the second week with settings a = 12, r = 1, l = 50, e = :01
# (we would ask for the de Bruijn graph except it will probably be too large; if you cannot create a
# contig map either that is okay);

reads1 = simulate_reads('C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\homework\hw5\\fasta\SigmaKappaOmega9_genome_reduced.fasta',12,50,0.01,0)
assemble(reads1, 25, [25, [2,-2,-2]], 0.1)