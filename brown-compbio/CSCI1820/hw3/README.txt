README    BEN SIRANOSIAN    CS1820 HW3
---------
PROBLEM 0

code for this problem is in alignment_with_inversions.py
this program reports the optimal local alignment with inversions. 
default parameters: match score = 1, mismatch score = -1, gap score = -1, inversion score = -1

USAGE
python alignment_with_inversions.py fasta1 fasta2
OPTIONAL ARGUMENTS
scores can be set with the following:
--match
--mismatch
--gap
--invert

---------
PROBLEM 1

Code for this problem is in contamination.py

USAGE
python contamination.py seqs vectors
where seqs is the file of sequences to test, in the format given on the website.
      vecttors is the file of vecotrs to test, in the format given on the website. 

This program returns to stdout a matrix of sequences and vectors that is true if contamination has occured. Sequences are rows, vectors are columns. 