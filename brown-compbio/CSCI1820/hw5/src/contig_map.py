from local_alignment import *
import string 
from PIL import Image, ImageDraw

#contig_map(readfile, contigfile):creates a dictionary of where reads map from the input reads and contigs. Uses the error_rate from generation for mapping
def contig_dict(reads, contigs, error_rate): 
	#where does each read map? do gapless local alignment. allow a number of mismatches based on error rate
	read_dict ={}
	for read in reads:
		for contig in contigs: 
			score = local_alignment(read, contig, 1, -1, -99, True)
			#print score
			alignment = local_alignment(read, contig, 1, -1, -99, False)
			#print alignment
			threshold_score = len(read) - 5 #(int(round(error_rate * len(contig)))+1)
			if score >= threshold_score:
				alignment_pos = string.find(contig, alignment[1])

				#print "STORING! " + str([contigs.index(contig), alignment_pos])
				if read in read_dict:
					#print "append case"
					read_dict[read].append([contigs.index(contig), alignment_pos])
				else:

					read_dict[read] = [[contigs.index(contig), alignment_pos]]
	return read_dict

#contig_map(read_dict, contigs): does the actual image processing for the contig map given the calcluated dictionary of reads and list of contigs
def contig_map(read_dict, contigs, outfile):
	total_contig_length = 0
	for contig in contigs:
		total_contig_length += len(contig)
	contig_ends=[]
	for i in range(len(contigs)):
		if i ==0:
			contig_ends.append(len(contigs[i])+10)
		else:
			contig_ends.append(contig_ends[i-1] + len(contigs[i]) +10)


	W = total_contig_length + (len(contigs)*10) + 10
	H = 500
	img = Image.new("RGB", (W, H), "white")
	draw = ImageDraw.Draw(img)

	for i in range(len(contigs)):
		#draw contigs in black at bottom
		if i == 0:
			start = 1
		else:
			start = contig_ends[i-1] +10 
		end = contig_ends[i]
		draw.line([start+5, 490, end+5, 490], fill="black",width=10)

	#calculate positions of reads
	read_xpos = []
	for i in read_dict.items():
		for j in i[1]:
			if j[0]==0:
				start = j[1] +5
			else:
				start = contig_ends[j[0]-1] + j[1] + 15
			end = start +len(i[0])
			read_xpos.append([start,end])

	# find what y pixel they should map to. see if the range is already occcupied
	# default to 40 available levels. ideally This would be dynamic
	levels = [[] for a in range(40)]
	for read in read_xpos:
		r = range(read[0],read[1])
		found_level = False
		current_level = 0
		while not found_level:
			found_overlap = False
			for pos in r:
				if pos in levels[current_level]:
					found_overlap = True
			levels[current_level] += r
			if not found_overlap:
				#we can plot on that level!
				found_level = True
			else:
				current_level += 1
		#draw the reads as lines
		draw.line([read[0], 480-(current_level*8),read[1],480-(current_level*8)],fill='red')
	
	#save the image
	img.save(outfile,"PNG")
