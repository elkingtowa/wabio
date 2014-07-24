import random
#random.seed(1)

#reverse_complement(pattern): helper function for simulate_reads. returns the reverse complement of a DNA sequence
# assumes an alphabet of {A,T,C,G}
def reverse_complement(pattern):
	chars = list(pattern)
	complement = ''
	lookup = dict({("A","T"),("T","A"),("C","G"),("G","C")})
	for base in chars:
		complement += lookup[base]
	rev_complement = complement[::-1]
	return rev_complement 

#simulate_reads(genome, coverage, length, error_rate): reads an input genome and simulates reads according to the parameters
# coverage: the desired coverage from sampling the genome
# length: the mean length of the reads. Actual lengths will be distributed as a normal distribution with a standard deviation of 3
# error_rate: number [0,1] describing the probability of a sequencing error at any given base pair in the read
# rc: boolean(0,1) representing whether to add the reverse complement of each sequence to the output. 
# outfile: if not an empty string, output a list of newline delimited reads to this output file. 
def simulate_reads(genome, coverage, length, error_rate, rc, outfile):
	#open genome file and read sequence into string
	f = open(genome, 'r')
	lines1 = f.read().split('\n')
	genome = ''
	for line in lines1:
		if (len(line) > 0):
			if (line[0] != ">"):
				genome=genome+line
	f.close()
	G = len(genome)

	#how many reads should we create?
	N = int(round(float(coverage)*float(G)/float(length)))

	#pick N starting positions over [0-G]
	starts = []
	for i in range(N):
		starts.append(random.randrange(G))
	#pick N read lengths from L~N(length, 3)
	lengths = []
	for i in range(N):
		lengths.append(int(round(random.normalvariate(length,3))))
	
	#pick reads from genome. Deal with edge effects by discarding reads that continue past the end of the genome.
	reads = []
	for i in range(N):
		if starts[i]+lengths[i] <= G:
			reads.append(genome[starts[i]:starts[i]+lengths[i]])

	#add errors to the reads as specified by error_rate
	for read in range(len(reads)):
		for pos in range(len(reads[read])):
			if random.random() < error_rate:
				letters = ['A','T','C','G']
				letters.remove(reads[read][pos])
				reads[read] = reads[read][0:pos] + random.choice(letters) + reads[read][pos+1:]

	#if rc=1, add the reverse complement of every read to the output. 
	# reverse_complement(pattern): returns the reverse complement of DNA string pattern. Equivalent to an inversion. 
	if rc ==1:
		rc_reads =[]
		for read in reads:
			rc_reads.append(reverse_complement(read))
		reads += rc_reads

	#write output if desired
	if outfile != '':
		of = open(outfile,'w')
		for read in reads:
			of.write(read+'\n')
		of.close()

	#return list of reads
	return reads
