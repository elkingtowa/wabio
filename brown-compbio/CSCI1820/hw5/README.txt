README		BEN SIRANOSIAN 		CSCI1820 HW5
--------------------------------------------
PROBLEM 1 - GENOME ASSEMBLY

All code can be found in the code folder.
This README documents the usage of code for this problem. For an explanation of the algorithms used and similar things, please see the report pdf. 

The only piece of code you need to interact with is my_velvet.py. This script takes several variables as input and completes the genome assembly process. 

my_velvet works in two modes: genome and read. 
In genome mode, reads are generated from a specified genome according to the parameters coverage, length, error rate, rc.
	reads are printed to a newline delimited output file. In this mode, you must specify the optional arguments --coverage, 
	--length, --error_rate, --rc
In read mode, reads are given as a newline delimited read file. 
In both modes, assembly proceeds at a given k. after assembly, contigs are printed to a newline delimited contig file.
You also have the option of creating a contig map image after the assembly is complete. The image will be saved to the specified file 
	in .png format. This can take a long time for a large graph!

USAGE: my_velvet.py [-h] [--coverage coverage] [--length length]
             [--error_rate error_rate] [--rc rc]
             [--read_outfile read_outfile]
             [--coverage_outfile coverage_outfile]
             mode inputfile k contigfile
HELP:
My implementation of the VELVET assembler! Simulates reads from an input
genome or takes as input a list of reads and attempts genome assembly. Pay
careful attention to the parameters becasue they can change your output
widely! First, specify genome or read mode with the first argument.

positional arguments:
  mode                  specify a mode. must be exactly genome or read
  inputfile             The file name of the fasta sequence containing the
                        genome, or text file containing one read per line
                        (each read must be of the same length).
  k                     An integer value k for constructing the graph.
  contigfile            Name of text file to write contigs at the end of
                        assembly.

optional arguments:
  -h, --help            show this help message and exit
  --coverage coverage   An integer coverage for simulating reads. Necessary if
                        in genome mode
  --length length       An integer read length for simulating reads. Necessary
                        in genome mode
  --error_rate error_rate
                        An float error rate for simulating reads. Necessary in
                        genome mode
  --rc rc               An boolean integer for rverse complement in simulating
                        reads. Necessary in genome mode
  --read_outfile read_outfile
                        Specify a file to output reads to here. By default,
                        nothing is saved.
  --coverage_outfile coverage_outfile
                        Specify a file to output a coverage map. By default,
                        nothing is saved.
