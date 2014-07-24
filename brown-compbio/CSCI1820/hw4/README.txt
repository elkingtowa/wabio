README    BEN SIRANOSIAN    CS1820 HW4
---------
PROBLEM 1 - RNAseq

code for this problem is in RNAseq.py
this code aligns RNA reads to a genome and calls areas that are potentially genes. Genes are identified by significant enrichment in a sliding region of 500bp.
the genome is read from genomefile, the reads from readfile. results are outputed in tab separated format to outfile.

USAGE
python RNAseq.py genomefile readfile outfile
OPTIONAL ARGUMENTS 	
--t 	sets the threshold for number of alignments in a sliding window of 400bp to be considered a gene. DEFAULT=35

Depends on the following python packages:
argparse
re
numpy
scipy 
matplotlib
---------
PROBLEM 2 - De Bruijn graphs

code for this problem is in de_brujin.py
this script reads an input fasta file and creates the De Bruijn representation of the sequence with the parameter k. By default the script tries to draw the graph to the screen, this should be supressed for large inputs because it wll take _forever_
The graph is outputed in the adjacency list format to the outfile specified

USAGE python de_bruijn.py genomefile k outfile
OPTIONAL ARGUMENTS
--draw	set to False to supress drawing the graph

Depends on the following python packages:
argparse
networkx
matplotlib
numpy