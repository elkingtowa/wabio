#parameters 
import argparse

parser = argparse.ArgumentParser(description="Script to create a dot-plot from two FASTA sequences of DNA, RNA or protein")
parser.add_argument("fasta1", action="store", metavar="fasta1", help="The file name of the fasta sequence to go on the horizontal of the plot")
parser.add_argument("fasta2", action="store", metavar="fasta2", help="The file name of the fasta sequence to go on the vertical of the plot")
args=parser.parse_args()

#parameters and paths specified in section above
fasta1=args.fasta1
fasta2=args.fasta2

# MAIN CODE
from pylab import *

#dot_polt(fasta1, fasta2): creates a dot plot by comparing two sequences from fasta files.
# fasta1 and fasta2 must be fasta formatted files containing a single sequence. 
def dot_plot(fasta1, fasta2):
	f1=open(fasta1, 'r')
	f2=open(fasta2, 'r')
	hseq = ''
	vseq = '' 

	#Read fasta files into a sequence
	lines1 = f1.read().split('\n')
	for line in lines1:
		if (len(line) > 0):
			if (line[0] != ">"):
				hseq=hseq+line

	lines2 = f2.read().split('\n')
	for line in lines2:
		if (len(line) > 0):
			if (line[0] != ">"):
				vseq=vseq+line

	#Inititalize matrix to zero
	score_matrix = [[0 for foo in hseq] for foo2 in vseq]
	#Loop through the sequences and mark positions where bases match
	for i in range(len(vseq)):
		for j in range(len(hseq)):
			if vseq[i]==hseq[j]:
				score_matrix[i][j] =1

	#Plot the matrix with grayscale
	figure(1)
	imshow(score_matrix, interpolation ="nearest", cmap=cm.gray_r)
	show()
	
#Function call
def main():
	dot_plot(fasta1,fasta2)

if  __name__ =='__main__':main()