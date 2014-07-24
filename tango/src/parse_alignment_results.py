# parse alignment files into a matrix of distances
# read phage names from file
phageNameF = 'C:\\Users\\Admin\\Documents\\GitHub\\tango\\data\\hatful60phageNames.txt'
names = []
namesCluster = []
with open(phageNameF, 'r') as nf:
	origline = nf.readline()
	line1 = origline.strip().split('(')
	line2 = origline.strip().split('\t')

	while line1 != ['']:
		names.append(line1[0])
		namesCluster.append(line2[0])
		origline = nf.readline()
		line1 = origline.strip().split('(')
		line2 = origline.strip().split('\t')

# gets similarity number from output of stretcher multiple alignment
def getSimilarity(stretcherFile):
	with open(stretcherFile,'r') as sf:
		line=sf.readline().strip()
		c = 0
		while c < 50: 
			#print line[0:13]
			if line[0:13] == '# Similarity:':
				return float(line.strip().split()[2].split('/')[0]) / float(line.strip().split()[2].split('/')[1])
			else: 
				line=sf.readline().strip()
			c+=1
alignmentDir = 'Z:\\projects\\phage_fasta\\stretcher_out\\'
similarityMatrix = [[1.0 for x in range(len(names))] for x in range(len(names))]

# grab numbers for each alignment 
for i in range(len(names)):
	for j in range(i+1, len(names)):
		stretcherFile = alignmentDir+str(names[i])+'_'+str(names[j])
		similarityMatrix[i][j] = getSimilarity(stretcherFile)
		similarityMatrix[j][i] = getSimilarity(stretcherFile)

#write results to an output matrix
outfile = 'C:\\Users\\Admin\\Documents\\GitHub\\tango\\data\\pairwise_alignment_similarity_60.tsv'
with open(outfile,'w') as of:
	header = ''
	for nameC in namesCluster:
		header += '\t'+nameC
	of.write(header+'\n')
	for row in range(len(similarityMatrix)):
		toWrite = namesCluster[row]
		for col in range(len(similarityMatrix)):
			toWrite += '\t'+ str(similarityMatrix[row][col])
		of.write(toWrite+'\n')

#write results to an output matrix as 1-similarity
outfile = 'C:\\Users\\Admin\\Documents\\GitHub\\tango\\data\\pairwise_alignment_distance_60.tsv'
with open(outfile,'w') as of:
	header = ''
	for nameC in namesCluster:
		header += '\t'+nameC
	of.write(header+'\n')
	for row in range(len(similarityMatrix)):
		toWrite = namesCluster[row]
		for col in range(len(similarityMatrix)):
			if row==col:
				toWrite += '\t'+ str(0.0)
			else:
				toWrite += '\t'+ str(1-similarityMatrix[row][col])
		of.write(toWrite+'\n')

