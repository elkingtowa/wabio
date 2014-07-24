import csv
import os.path
from Bio import SeqIO
from Bio.Blast import NCBIWWW

#open up CSV file to parse

with open('blast/blast_all_filtered.csv', 'rb') as datafile:
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
			print fullname

			for seq_record in SeqIO.parse(fullname, "fasta"):
				print(seq_record.id)
				sequence = seq_record.seq

			#makes sure stop and start codon are in the correct order
			if startcodon < stopcodon:
				gene = sequence[startcodon: stopcodon]
			else:
				gene = sequence[stopcodon:startcodon]

			if not gene:
				print "gene sequence empty"
			else:
				print 'gene sequence: ' + gene

				outputfilelocation = 'NCBI BLAST results/' + sequencename[:50] + " " + phagename + '.xml'

				isFile = os.path.isfile(outputfilelocation)

				if not isFile:
					print 'BLASTing at NCBI database'
					#blasting
					result_handle = NCBIWWW.qblast("blastn", "nt", gene)

					save_file = open(outputfilelocation, 'w')
					save_file.write(result_handle.read())
					save_file.close()
					result_handle.close()
					print 'successfully wrote results to file'
				else:
					print 'file already exists'