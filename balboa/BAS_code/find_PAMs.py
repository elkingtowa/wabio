#using some code from Eddie's dataParsingScripy.py

import csv
import os
from Bio import SeqIO
from Bio.Blast import NCBIWWW

balboa_dir = 'C:\Users\Admin\Documents\GitHub\\balboa\\'
fasta_dir = 'C:\Users\Admin\Documents\GitHub\\balboa\phages\\fasta\\'

test = 'Sparky.fasta'
fullname=fasta_dir+test

#Extracts a region from the fasta file in filename.
# If start is greater than end, the reverse complement of the sequence is returned. 
# surrounding number of bases are included in the output. 
# If only_surrounding, returns a list of the surrounding bases from the match
def extract_region(filename, start, end, surrounding, only_surrounding):
	for seq_record in SeqIO.parse(filename, "fasta"):
		sequence = seq_record.seq

	if not only_surrounding:
		if start < end:
			return sequence[start-surrounding-1:end+surrounding].tostring()
		else:
			return sequence[end-surrounding-1:start+surrounding].reverse_complement().tostring()
	else: 
		if start < end:
			return [sequence[start-surrounding-1:start-1].tostring(),sequence[end:end+surrounding].tostring()]
		else:
			return [sequence[start:start+surrounding].reverse_complement().tostring(),sequence[end-surrounding-1:end-1].reverse_complement().tostring()]


#open and parse match databas. return a list of lists of data
# first element of data will be the column headers 
# qseqid,sseqid,filename,cluster,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore

import csv
datafile = 'C:\Users\Admin\Documents\GitHub\\balboa\\blast\\blast_all_filtered_remove_duplicates.csv'
def read_match_db(filename):
	with open(filename, 'rb') as datafile:
		csvreader = csv.reader(datafile)
		data =[]
		for row in csvreader:
			data.append(row)
	return data

data = read_match_db(datafile)
match_sequences = []
matches_plus10 = []
plus10 = []
for line in data[1:]:
	print "filename: " + fasta_dir+line[2]
	print "start: " + str(line[10]) + " end: " + str(line[11])
	match_sequences.append(extract_region(fasta_dir+line[2], int(line[10]), int(line[11]), 0, False))
	matches_plus10.append(extract_region(fasta_dir+line[2], int(line[10]), int(line[11]), 10, False))
	plus10.append(extract_region(fasta_dir+line[2], int(line[10]), int(line[11]), 10, True))


