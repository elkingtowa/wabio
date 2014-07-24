import argparse
import sys 
import string
from simulate_reads import *
from assemble import *
from contig_map import *

parser = argparse.ArgumentParser(description="My implementation of the VELVET assembler! Simulates reads from an input genome or takes as input a list of reads and attempts genome assembly. Pay careful attention to the parameters becasue they can change your output widely! First, specify genome or read mode with the first argument. ")
parser.add_argument("mode", action="store", metavar="mode", help="specify a mode. must be exactly genome or read")
parser.add_argument("inputfile", action="store", metavar="inputfile", help="The file name of the fasta sequence containing the genome, or text file containing one read per line (each read must be of the same length).")
parser.add_argument("k", action="store", metavar="k", type=int, help="An integer value k for constructing the graph.")
parser.add_argument("contigfile", action="store", metavar="contigfile", help="Name of text file to write contigs at the end of assembly.")
parser.add_argument("--coverage", action="store", metavar="coverage", type=int, help="An integer coverage for simulating reads. Necessary if in genome mode")
parser.add_argument("--length", action="store", metavar="length", type=int, help="An integer read length for simulating reads. Necessary in genome mode")
parser.add_argument("--error_rate", action="store", metavar="error_rate", type=float, help="An float error rate for simulating reads. Necessary in genome mode")
parser.add_argument("--rc", action="store", metavar="rc", type=int, help="An boolean integer for rverse complement in simulating reads. Necessary in genome mode")
parser.add_argument("--read_outfile", action="store", metavar="read_outfile", default='', help="Specify a file to output reads to here. By default, nothing is saved.")
parser.add_argument("--coverage_outfile", action="store", metavar="coverage_outfile", default='', help="Specify a file to output a coverage map. By default, nothing is saved.")


args=parser.parse_args()

#my_velvet(): organizes code and calls subpieces
def my_velvet():
	if args.mode =="genome":
		#simulate reads
		reads=simulate_reads(args.inputfile, args.coverage, args.length, args.error_rate, args.rc, args.read_outfile)

	elif args.mode == "read":
		#take in reads
		inf = open(args.inputfile, 'r')
		reads = []
		line = inf.readline()
		while line != '':
			reads.append(string.upper(line.strip()))
			line=inf.readline()
	else:
		sys.exit()


	#set default parameters for algorithms used. These can be changed here but I didn't feel it was necessary to make them part of the arguments
	# local alignment score required to merge sequences
	tour_bus_score = args.k - 2
	# parameters to pass to local alignment algorithm [match score, mismatch score, gap score]
	alignment_parameters = [1, -1, -1]
	#percentile cutoff for removing low coverage nodes at the last step
	coverage_cutoff = 0.1 


	# DO ASSEMBLY!
	contigs = assemble(reads, args.k, [tour_bus_score, alignment_parameters], coverage_cutoff, args.contigfile, False)

	if args.coverage_outfile != '':
		print "Constructing contig map"
		rdict = contig_dict(reads, contigs, args.error_rate)
		contig_map(rdict, contigs, args.coverage_outfile)

#function call
def main():
	my_velvet()

if  __name__ =='__main__':main()