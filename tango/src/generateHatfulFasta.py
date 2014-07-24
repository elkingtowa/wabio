from Bio import SeqIO

#parses through hatful names, adds them to array

filepath = '/Users/edwardwilliams/Documents/tango/data/hatful60phageNames.txt'
allphage = '/Users/edwardwilliams/Documents/balboa/phages/allphage.fasta'
outputfilepath = '/Users/edwardwilliams/Documents/tango/data/hatful.fasta'

nameArray = []

with open(filepath, 'r') as hatfulFile:
	for line in hatfulFile: 
		nameArray.append(line.split("(")[0])

recordArray = []

#pulls out sequence records of hatful phages
#there are some duplicates in this, but I think BLAST will ignore it

for seq_record in SeqIO.parse(allphage, "fasta"):
	for string in seq_record.id.split(" "):
		if string in nameArray:
			recordArray.append(seq_record)
			break


#writes recordArray to FASTA
SeqIO.write(recordArray, 'hatful.fasta', "fasta")