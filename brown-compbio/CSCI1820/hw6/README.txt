README    BEN SIRANOSIAN    CS1820 MIDTERM
---------
PROBLEM 1 - Scaffolding

Code to solve the scaffolding probelm can be found in scaffold.py
USAGE
python scaffold.py contigfile matefile
OUTPUT 
Order and orientation of computed scaffolds, one per line. Contig numbers are zero indexed. A negative sign indicates the reverse complement of the contig is included in the scaffold. 

depends on the following python packages:
argparse
math
regex
networkx

NOTES
I originally designed this algorithm to use my local alignment script from a previous homework. However, computing all possible overlaps with local alignment took far too long. To speed things up, I use fuzzy regexp matching with the regex package. I know local alignment is the ideal solution, but using the regexp made it possible to actually test my code. I left the original local alignemnt code in the file just in case you want to see it! 