import csv
import os.path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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

#pulls out unique sequences, saves them to data array, then writes them to a FASTA file for BLASTing locally.
dataArray = list()

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
			for record in dataArray:
				if str(gene) == str(record.seq):
					uniqueSequence = False

			print sequencename + " in " + phagename + " with sequence: " + str(gene)
			#if unique, create a SeqRecord for writing to FASTA file
			if uniqueSequence:
				print "sequence is unique, creating sequence record"

				record = SeqRecord(Seq(str(gene), IUPAC.ambiguous_dna), id = sequencename, description = "CRISPR found in " + filename + " cluster: " + cluster + " with PAM: " + str(PAM) )
				dataArray.append(record)

filename = "unique_CRISPRs.fasta"

print "writing data to " + filename

#write dataArray to FASTA file

SeqIO.write(dataArray,filename, "fasta")