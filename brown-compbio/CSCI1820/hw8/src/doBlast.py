import argparse
parser = argparse.ArgumentParser(description="My implementation of BLAST. Computes the highest scoring maximum segment pairs by aligning a query sequence (defined in the query fasta file) with the subject database (defined in the subject fast file). Additionally, reports alignment statistics in terms of the k* parameter. The subject database can be specified in a pickle file so that we don't have to generate it each time BLAST is called.")
parser.add_argument("query", action="store", metavar="query", help="The location of the fasta file containing the query sequence")
parser.add_argument("subject", action="store", metavar="subject", help="The location of the fasta file containing the subject sequence")
parser.add_argument("--databaseFile", action="store", metavar="databaseFile", default=None, help="You can specift a pickle file containing the location of each word in the database here. Should be precomputed by my generateDatabase.py script.")
parser.add_argument("--k", action="store", metavar="k", type=int, default=4, help="Word length to break up the query string. Default is 4.")
parser.add_argument("--T", action="store", metavar="T", type=int, default=10, help="Minimum score to create search terms after breking query sequence. Default is 10.")
parser.add_argument("--S", action="store", metavar="S", type=int, default=None, help="Minimum score to report MSPs. You can specify a fixed value here, or use the default which is (0.6*5*|qurey|)")
parser.add_argument("--threshold", action="store", metavar="threshold", type=float, default=0.75, help="Keep extending for matches until score falls below this fraction of max score. default is 0.75")
parser.add_argument("--substitutionMatrix", action="store", metavar="substitutionMatrix", type=str, default='[5,-4]', help="Substitution matrix for scoring matches. Must be specified as a python formatted list, where the first number is the score for a match and the second number is the score for a mismatch. Default is [5,-4]")
args=parser.parse_args()

#doBlast - Wrapper for all the blast functions!
from blast import *
from generateDatabase import *
import ast

#parseFasta: return a sequence by reading lines in a fasta file
def parseFasta(inFile):
	with open(inFile,'r') as fastafile:
		sequence=''
		line=fastafile.readline().strip()
		while line != '':
			if line[0] != '>':
				sequence +=line
			line=fastafile.readline().strip()
	return sequence

#gets sequences from files
querySeq = parseFasta(args.query)
subjectSeq = parseFasta(args.subject)
#get database if defined
if args.databaseFile != None:
	with open(args.databaseFile, 'r') as inFile:
		database = pickle.load(inFile)
else: database = None
#adjust S if defined
if args.S != None:
	blastS=int(args.S)
else: 
	blastS='default'

#get blast results in MSP list
print 'computing alignments...'
mspList = blast(querySeq, subjectSeq, database, threshold=float(args.threshold), k=int(args.k), T=int(args.T), S=blastS, score=ast.literal_eval(args.substitutionMatrix))

# sometimes matches within another match are reported.
# eliminate matches that are completely contained within another match
uniqueList = []
for i in range(len(mspList)):
	#compare agains all others 
	keepMe = True
	for j in range(len(mspList)):
		if i!=j:
			#start1 greater than start2 and end1 less than end2 for both query and subject
			# means match is contained within another match
			if ((mspList[i][2][0]>=mspList[j][2][0]) and (mspList[i][2][1]<=mspList[j][2][1])) and ((mspList[i][3][0]>=mspList[j][3][0]) and (mspList[i][3][1]<=mspList[j][3][1])):
				keepMe = False 
	if keepMe:
		uniqueList.append(mspList[i])

#pretty print 
print 'Found ' + str(len(uniqueList)) +' unique alignments'
for u in uniqueList:
	matchString = ''
	for pos in range(len(u[0])):
		if u[0][pos] == u[1][pos]:
			matchString+=('|')
		else:
			matchString+=(' ')

	print 'Query\t'+str(u[2][0]) + '\t' + u[0] + '\t' + str(u[2][1])
	print '\t\t'+matchString
	print 'Sbjct\t'+str(u[3][0]) + '\t' + u[1] + '\t' + str(u[3][1])
	print 'SCORE: ' + str(u[4])
	#get statistics
	# compute lambda for our matrix 
	lambdaStar = approximateLambda(ast.literal_eval(args.substitutionMatrix))
	stats = doStats(lambdaStar, len(uniqueList), len(querySeq),len(subjectSeq),u[4])

	print 'E = ' + str(stats[1])
	print 'Ps,m,n('+str(len(uniqueList))+') = ' + str(stats[0])