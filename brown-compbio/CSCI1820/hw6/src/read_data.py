import os
os.chdir('C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\midterm\code')
c = 'C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\midterm\data\contigs.txt'
m = 'C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\midterm\data\mate_pairs.txt'
m2 = 'C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\midterm\data\mate_pairs_100.txt'


#reverse_complement(pattern): helper function for simulate_reads. returns the reverse complement of a DNA sequence
# assumes an alphabet of {a,t,c,g}
def reverse_complement(pattern):
	chars = list(pattern)
	complement = ''
	lookup = dict({("a","t"),("t","a"),("c","g"),("g","c")})
	for base in chars:
		complement += lookup[base]
	rev_complement = complement[::-1]
	return rev_complement 


# read data from contigs and mate pairs 
def read_data(contig_file, mate_file):
	cf = open(contig_file, 'r')
	mf = open(mate_file, 'r')

	contigs = []
	mates = []

	cline = cf.readline()
	while cline != '':
		contigs.append(cline.strip())
		cline = cf.readline()
	#append reverse complement of all contigs 
	for i in range(len(contigs)):
		contigs.append(reverse_complement(contigs[i]))

	mline=mf.readline()
	while mline != '':
		mline = mline.strip().split('\t')
		mates.append([mline[0],mline[1],mline[2]])
		mline=mf.readline()

	mf.close()
	cf.close()

	return contigs, mates


#do gapless local alignment. Find where each end of the mate maps and in which direction?
from local_alignment import *
import math
def align_mates(contigs, mates):
	#record what mates match each contig in a dict
	contig_matches = {}

	for i in range(len(contigs)):
		contig = contigs[i]
		for j in range(len(mates)):
			mate1 = mates[j][0]
			mate2 = mates[j][2]
			alignment_score1 = local_alignment(contig, mate1, 1, -1, -99999999, True)
			alignment_score2 = local_alignment(contig, mate2, 1, -1, -99999999, True)
			#allow for 2% errors, 
			threshold_score = len(mate1) - int(math.ceil(len(mate1)*0.02))
			print "contig: " + str(i) + " mate1: " + str(j) + " al1: " + str(alignment_score1) + " thresh: " + str(threshold_score)
			print "al2: " + str(alignment_score2) + " thresh: " + str(threshold_score)
			match1 = alignment_score1 >= threshold_score
			match2 = alignment_score2 >= threshold_score

			#if both ends match the contig, discard the mate
			if not(match1 and match2):
				if match1:
					if i in contig_matches.keys():
						contig_matches[i].append([j,1])
					else:
						contig_matches[i]=[[j,1]]
				elif match2:
					if i in contig_matches.keys():
						contig_matches[i].append([j,2])
					else:
						contig_matches[i]=[[j,2]]
	return contig_matches


## Try again using regexps with a certain sumber of mismatches to speed things up. 

#generates a list of sequences with errors by substituting every letter in alphabet at each position
# does one error, run multiple times to get multiple errors. 
def generate_errors(seq_list, alphabet):
	return_list = set()
	for seq in seq_list:
		for pos in range(len(seq)):
			for letter in alphabet:
				return_list.add(seq[:pos]+letter+seq[pos+1:])
	return list(return_list)

import math
#do the opposite way... iterate mates first and record those
def align_mates2(contigs, mates):

	mate_matches = dict()
	for i in range(len(mates)):
		mate_matches[i] = [[],[]]

	for i in range(len(mates)):
		mate1 = mates[i][0]
		mate2 = mates[i][2]
		threshold_score = len(mate1) - (int(math.ceil(len(mate1)*0.02))*2)
		for j in range(len(contigs)):
			contig = contigs[j]
			alignment_score1 = local_alignment(contig, mate1, 1, -1, -99999999, True)
			alignment_score2 = local_alignment(contig, mate2, 1, -1, -99999999, True)

			match1 = alignment_score1[0] >= threshold_score
			match2 = alignment_score2[0] >= threshold_score

			if not(match1 and match2):
				if match1:
					mate_matches[i][0].append([j, alignment_score1[1][0]])
					print "mate " + str(i) + " left matched in contig " + str(j)
				elif match2:
					mate_matches[i][1].append([j, alignment_score2[1][0]])
					print "mate " + str(i) + " right matched in contig " + str(j)
	return mate_matches

#c = 'contigs.txt'
#m = 'mate_pairs.txt'
#data = read_data(c,m)
data2 = read_data(c,m2)
#matches = align_mates(data[0],data[1])
matches2 = align_mates2(data2[0],data2[1][0:10])

import csv
w = csv.writer(open("output.csv", "w"))
for key, val in dict.items():
	w.writerow([key, val])
w.close()


