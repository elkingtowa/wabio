#parameters 
import argparse

parser = argparse.ArgumentParser(description="Script to align RNAseq reads to a genome and report coordinates of suspected genes. Identifies genes as regions with greater than 't' alignments in a sliding window of 500bp.")
parser.add_argument("genomefile", action="store", metavar="genomefile", help="The file name of the fasta sequence containing the genome")
parser.add_argument("readfile", action="store", metavar="readfile", help="The file name of the fasta sequence containing the RNAseq reads, one read per line")
parser.add_argument("outfile", action="store", metavar="outfile", help="The file name to write results to.")
parser.add_argument("--t", action="store", metavar="significance", default="35", help="threshold for gene calling.")
args=parser.parse_args()

genomefile=args.genomefile
readfile=args.readfile
outfile=args.outfile
significance=int(args.t)

import re
import numpy
import matplotlib.pyplot as plt
from scipy.stats import binom

def RNAseq(genomefile,readfile, significance, outfile):
	#read genome and readfile
	#genome=open("C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\homework\hw4\seqs\SigmaKappaOmega9_genome.fasta",'r')
	genomef=open(genomefile,'r')
	lines1 = genomef.read().split('\n')
	genome = ''
	for line in lines1:
		if (len(line) > 0):
			if (line[0] != ">"):
				genome=genome+line
	L = len(genome)
	genomef.close()

	#reads=open("C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\homework\hw4\seqs\8mer_reads.fasta",'r')
	readf=open(readfile,'r')
	lines2 = readf.read().split('\n')
	reads = []
	for line in lines2:
		if (len(line) > 0):
			if (line[0] != ">"):
				#substitute for RNA
				reads.append(line.replace("U","T"))
	readf.close()
	print "Finished reading genomes"
	#are there duplicates in this read set?  YES
	#duplicates = set([x for x in reads if reads.count(x)>1])
	#num_duplicates = [reads.count(x) for x in duplicates]

	#find all alignments of each read. report indicies with regexp match
	read_alignments= {read:[] for read in reads}

	for read in reads:
		iterator = re.finditer(read, genome)
		read_alignments[read].append([m.start() for m in iterator])

	# mean number of alignments =5
	num_alignments = [len(read) for read in read_alignments.values()]
	average_alignments = numpy.mean(num_alignments)

	#distribution of alignments
	start_positions = []
	for val in read_alignments.values():
		for pos in val:
			for start in pos:
				start_positions.append(start)
	#alignments at each base pair - takes a long ass time
	print "calculating alignments --  this can take a minute..."
	bp_alignments = [start_positions.count(x) for x in range(L)]

	#number of alignments in 500bp sliding window
	sliding_alignments500 = []
	for l in range(L-500):
		sliding_alignments500.append(sum(bp_alignments[l:l+500]))
	#plt.plot(range(L-500),sliding_alignments500)

	#gene calling! 
	print "calling genes..."
	allowed_below_threshold = 75
	window_size = 500
	genes = []

	position = 0
	in_gene =  False
	gstart =0
	gend =0
	bases_below=0

	# look over sliding window data
	while position < len(sliding_alignments500):
		if in_gene:
			# If less than significance but not too many bases below it, keep counting
			if sliding_alignments500[position] <= significance:
				# If too many positions below significance, stop gene 
				if bases_below > allowed_below_threshold:
					gend = position
					in_gene=False
					genes.append([gstart+(window_size/2),gend+(window_size/2)-allowed_below_threshold])
					bases_below=0
				else:
					bases_below+=1
		elif sliding_alignments500[position] > significance:
				in_gene=True
				gstart=position
				bases_below=0
		position +=1 


	# get number of alignments and probability of alignments in each gene
	# P = probability of a given read aligning to a given position under the null model. 
	P = float(len(reads)) / 65536.0	

	for gene in genes:
		gene.append(sum(bp_alignments[gene[0]:gene[1]]))
		gene.append(1- binom(gene[1]-gene[0],P).cdf(gene[2]))

	#write tab-separated results to an output file. 
	#outfile="C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\homework\hw4\seqs\RNAseq.out"
	of = open(outfile, 'w')
	of.write('start\tend\talignments\tprobability\n')
	for gene in genes:
		of.write(str(gene[0]) + '\t' + str(gene[1]) + '\t' +str(gene[2]) + '\t' +str(gene[3])+'\n')
	of.close()

	fig=plt.figure()
	fig.suptitle("Alignments in 500bp sliding window. Gene positions highlighted.", fontsize=16)
	plt.plot(range(L-500),sliding_alignments500)
	for gene in genes:
		plt.plot([gene[0],gene[1]],[significance+5,significance+5],'k-',lw=5)
	plt.show()

#function call
def main():
	RNAseq(genomefile,readfile,significance,outfile

if  __name__ =='__main__':main()