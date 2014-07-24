#code that was originally in tetranucleotideAnalysis.py written by Eddie.
#probably depreciated due to changes at this point. 

#open a list of tab separated phage names (col1) and filenames (col2), one per line. 
#calculates TUD for each and writes results to outfile. 
def phageTUDCalc(phageFile, outfile, k,RC=False):
	#open and read in names and filenames
	with open(phageFile, 'r') as pf:
		phages = []
		line = pf.readline()
		while line != '':
			phages.append(line.strip().split('\t'))
			line=pf.readline()

	with open(outfile, 'w') as of:
		kmerList = enumerateKmers(k)
		#get one tud, write kmers at top
		oneTud = TUD(phages[0][1], k, kmerList,RC=RC)
		kmers = ''
		for kmer in oneTud.keys():
			kmers += '\t' + kmer
		of.write(kmers+'\n')
		#get TUD data for each phage
		for phage in phages:
			name = phage[0]
			print('working with ' + name)
			fname = phage[1]
			tud = TUD(fname, k, kmerList,RC=RC)
			#write each result to file
			toWrite = name
			for j in tud.values():
				toWrite += '\t' + str(j)
			of.write(toWrite+'\n')


#similar to phageTUDCalc, except uses BioPython to parse through a single FASTA file with
#a large number of mycobacteria 
#modified to produce CSV files to reduce R processing errors
def mycobacteriaTUDCalc(filename, outfile, k):

	with open(outfile, 'w') as of:
		kmerList = enumerateKmers(4)
		wroteKmersAtTop = False

		for seq_record in SeqIO.parse(filename, "fasta"):
			sequence = seq_record.seq.tostring().upper()
			name = seq_record.description.split("| ")[1].replace(" ", "_").replace(',', '-')

			tud = TUDFromString(sequence, k, kmerList)

			#if it's the first time, write kmers at the top of the file
			#else, skip over
			if not wroteKmersAtTop:
				kmers = ''
				for kmer in tud.keys():
					kmers += '\t' + kmer
				of.write(kmers+'\n')
				wroteKmersAtTop = True

			print('working with ' + name)
			tud = TUDFromString(sequence, k, kmerList)
			#write each result to file
			toWrite = name
			for j in tud.values():
				toWrite += '\t' + str(j)
			of.write(toWrite+'\n')
