import csv
import os.path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Alphabet import IUPAC



#opens file, parses through every other line (the crispr name)
#checks if name matches list of CRISPR names appearing in phages


validCRISPRs = ['NC_019951_9_21',
'NC_019951_8_42',
'NC_019951_9_20',
'NC_019951_9_9',
'NC_000962_6_12|NC_008769_4_13|NC_009525_6_12|NC_012207_4_13|NC_012943_2_12|NC_015758_5_14|NC_016768_2_12|NC_016804_4_13|NC_017524_8_4|NC_017528_7_12|NC_018078_2_12|NC_018143_6_12|NC_020089_4_25|NC_020245_4_13|NC_020559_5_14|NC_021054_8_11',
'NC_000962_6_9|NC_002755_6_10|NC_002945_5_10|NC_008769_4_10|NC_009525_6_9|NC_009565_6_7|NC_012207_4_10|NC_012943_2_15|NC_015758_5_11|NC_016768_2_15|NC_016804_4_10|NC_017524_8_1|NC_017528_7_9|NC_018078_2_15|NC_018143_6_9|NC_020089_4_22|NC_020245_4_10|NC_021054_8_8',
'NC_017026_7_8',
'NC_002755_5_12|NC_002945_4_10|NC_008769_3_12|NC_012207_3_12|NC_015758_4_23|NC_016804_3_12|NC_017026_6_18|NC_020089_4_12|NC_020245_3_12|NC_020559_4_12',
'NC_015848_5_21',
'NC_019952_4_2',
'NC_019951_9_23',
'NC_019951_9_27',
'NC_015758_4_10',
'NC_015848_5_23',
'NC_002755_6_3|NC_002945_5_4|NC_008769_4_4|NC_012207_4_4|NC_015758_5_4|NC_016804_4_4|NC_020089_4_15|NC_020245_4_4|NC_020559_5_4',
'NC_015758_4_11',
'NC_021192_10_12',
'NC_020276_1_1',
'NC_019951_9_8',
'NC_019952_5_10',
'NC_019951_8_24',
'NC_019951_8_41',
'NC_019951_9_16',
'NC_019951_9_24',
'NC_019951_9_34',
'NC_019952_4_4',
'NC_019965_7_9',
'NC_017026_8_3',
'NC_019951_8_49',
'NC_019951_9_28',
'NC_019951_9_42',
'NC_019951_9_11',
'NC_015848_5_15|NC_019950_7_15',
'NC_015848_5_7|NC_019950_7_7',]

CRISPRdict = {'crisprName': 'sequence'}
lineNumber = 1
crisprName = 'nothing'
crisprSequence = 'nothing'

#parses through sequence file, adds them to dictionary, removes names
with open('blast/mycobacteria_spacer_sequences_inc_names.txt', 'r') as nameFile:
	for line in nameFile:
		if lineNumber % 2 == 0:
			crisprSequence = line.replace("\n", "")

			#pair completed, can add to dictionary
			CRISPRdict[crisprName] = crisprSequence
		else:
			crisprName = line.replace(">", "").replace("\n", "")
		lineNumber += 1

#goes through list of valid CRISPR sequences. 
#problem: this doesn't take the sequences from phage that matched the sequence, just directly from the CRISPR sequences. 
#attempting to do different reading frames may result in either losing information or having a non- multiple of 3 sequence
#a better solution would be to parse through the CSV file again, and look at the sequence directly from the phage. 
#then, the script could pull out sequences before and after the CRISPR. but those could be different across phages.
#so, a good course of action would be to parse through the CSV file, and keep track if a sequence pulled out is unique or not.
#then, blast all unique sequences at once. 
#this might not even need the dictionary? which would be kind of sad face, because the dictionary was cool. 
#for crisprSequence in validCRISPRs:
	#untranslatedSequence = CRISPRdict[crisprSequence] 
	#print "CRISPR sequence: " + untranslatedSequence
	#translates sequence:
	#uses several reference frames
	#codingSequenceRF1 = Seq(untranslatedSequence, IUPAC.ambiguous_dna)
	#print "Translated amino acid sequence in RF1: " + codingSequenceRF1.translate()

#array to hold all info about unique sequences
#many phages have the same CRISPR as a blast match, but with slightly different sequences.
#the program adds unique sequences with information about their origins. 
dataArray = []
dataArray.append(["CRISPR Name", "FASTA file name", "phage cluster", "gene sequence", "start codon", "stop codon", "nucleotides before/after start codon", "PAM", "protein in RF1", "protein in RF2", "protein in RF3", "protein in RF4", "protein in RF5", "protein in RF6"])

with open('blast/blast_all_filtered.csv', 'r') as datafile:
	csvreader = csv.reader(datafile)
	for row in csvreader:
		if row[0] != 'qseqid':
			sequencename = row[0]
			phagename = row[1]
			filename = row[2]
			cluster = row[3]
			startcodon = int(row[10])
			stopcodon = int(row[11])

			#pulls out sequence from fasta file
			fullname = 'phages/fasta/' + filename

			for seq_record in SeqIO.parse(fullname, "fasta"):
				print(seq_record.id)
				sequence = seq_record.seq

			#makes sure stop and start codon are in the correct order
			#pulls before and after the gene
			deltaN = 5
			PAM = ""

			if startcodon < stopcodon:
				gene = sequence[startcodon - deltaN : stopcodon + deltaN]
				PAM = sequence[startcodon - deltaN : startcodon]
			else:
				gene = sequence[stopcodon - deltaN :startcodon + deltaN]
				PAM = sequence[stopcodon - deltaN: stopcodon]

			

			#parse dataArray to make sure name is unique
			uniqueSequence = True
			repeated = False
			repeatNumber = 0
			for dataRow in dataArray:
				if str(gene) == str(dataRow[3]):
					uniqueSequence = False

			print sequencename + " in " + phagename + " with sequence: " + str(gene)
			if uniqueSequence:
				print "is unique sequence, processing proteins"

				#performs translation in 6 different reference frames:
				#starting from the beginning of the sequence (including the genes at the beginning, not part of the CRISPR match!):

				proteinRF1 = gene.translate()
				#second nucleotide
				proteinRF2 = Seq(str(gene)[1:], IUPAC.ambiguous_dna).translate()
				#third nucleotide
				proteinRF3 = Seq(str(gene)[2:], IUPAC.ambiguous_dna).translate()

				#backwards
				proteinRF4 = Seq(str(gene)[::-1], IUPAC.ambiguous_dna).translate()
				proteinRF5 = Seq(str(gene)[::-1][1:], IUPAC.ambiguous_dna).translate()
				proteinRF6 = Seq(str(gene)[::-1][2:], IUPAC.ambiguous_dna).translate()

				dataArray.append([sequencename, filename, cluster, gene, startcodon, stopcodon, deltaN, PAM, proteinRF1, proteinRF2, proteinRF3, proteinRF4, proteinRF5, proteinRF6])

#writes to data file:
#probably not the most efficient way to do it, I could write the row to the CSV file as it is created

fileName = 'proteinInfoNotBLASTED.csv'

print "writing to CSV file: " + fileName 
with open(fileName, 'wb') as f:
	writer = csv.writer(f)

	for dataRow in dataArray:
		writer.writerow(dataRow)

#protein blasts each protein:
#output file format: CRISPR sequence name, phage of origin, PAM, refrence frame


for dataRow in dataArray[1:]:
	for x in range (8, 14):
		proteinSequence = dataRow[x]
		#truncates CRISPR name. Can occasionally cause a file name too long error. 
		outputfilelocation = "PBLAST results/" + dataRow[0][:50] + " " + dataRow[1] + " " + str(dataRow[7]) + " "  + str(x - 7) + ".xml"


		isFile = os.path.isfile(outputfilelocation)

		if not isFile:
			print 'PBLASTing at NCBI database \nCRISPR Name: ' + dataRow[0] + ' from phage ' + dataRow[1] + ' protein in RF' + str(x-7)
			#blasting
			result_handle = NCBIWWW.qblast("blastp", "nr" , proteinSequence)

			save_file = open(outputfilelocation, 'w')
			save_file.write(result_handle.read())
			save_file.close()
			result_handle.close()
			print 'successfully wrote results to file: ' + outputfilelocation
		else:
			print 'file already exists'
		

