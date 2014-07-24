README  |  Ben Siranosian  |  CSCI1820 hw2

There are three python scripts included in this handin. I will explain arguments and usage of each. 


----------------------
repeated_substrings.py
This script is used to find all repeated substrings of at least length k in the input text. Default parameters: k=3

USAGE: python repeated_substrings.py TEXT 
OPTIONAL ARGUMENTS: --k: find substrings of at least this length. 

EXAMPLE: 
python repeated_substrings.py ATGATGA --k 3 
ATG
ATGA
TGA

---------------------------------
repeated_substrings_mismatches.py
This script will find all repeated substings of at least length k with at most d mismatches over an alphabet in the input text. Default parameters: k=3, d=1, alphabet= A,T,C,G

USAGE: python repeated_substrings_mismatches.py TEXT
OPTIONAL ARGUMENTS: --k: find substrings of at least this length
					--d: find substrings with at most this many mismatches
					--alphabet: use this alphabet for finding mismatches 

EXAMPLE:
python repeated_substrings_mismatches.py ATGATG
AGG
TTG
ATG
ACG
 .
 .

-----------
dot_plot.py
This script will produce a dot-plot representing nucleotide identity between two sequences determined in fasta files. The first sequence will be placed on the horizontal of the plot, the second on the vertical 

USAGE: python dot_plot.py horizontal_fasta vertical_fasta
EXAMPLE: python dot_plot.py hrv16.fasta hrv1A.fasta