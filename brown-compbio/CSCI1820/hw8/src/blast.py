import sys
import os
import math
from generateDatabase import *

#blast(query, subject, databse, k=4, T=10, S=0.60*5*len(query)): The main function for blast! 
# Algorithm is described in the homework 8 handout. A database computed by generateDatabase must be supplied. 
# The database isn't generated every time and speeds up runtime. Also nets some sweet extra credit points. 
# Similarity mattix can be changed with the score option. specify a 2 item list with match and mismatch scores. 
# Can also specify None for database to have that generated on the spot
def blast(query, subject, database, threshold=0.5, k=4, T=10, S='default', score=[5,-4]):
	#if no database, make one
	if database ==None:
		database = generateDatabase(subject, k, fastaFile=False)

	#define S in terms of query if default
	if S=='default' : S=0.60*5*len(query)

	# break query into words of length k.
	if len(query) < k:
		print 'Query shorter than k. Exiting.'
		sys.exit()
	
	kList = []
	for i in range(len(query)-k+1):
		 kList.append(query[i:i+k])

	# find all words that we can match to with a score > T
	# Basically, how many mismatches can the word have?
	mismatches = 0
	for i in range(1,k+1):
		if ((k-i)*score[0]) + (i*score[1]) > T:
			mismatches += 1
		else: break 
	# generate a list of words that will score above T from the number of mismatches:
	# keep these indexed by the kmer that generated them for now
	wDict = generateMismatches(kList, mismatches)
	
	# hash the query so we know where each original word of length K came from
	queryMap = generateDatabase(query, k, False)

	# find MSPs for each item in kList by looking for all the mismatches. 
	# keep these indexed by the k for now?
	mspDict = dict((kmer, []) for kmer in kList)
	
	#print " STARTING MSP FINDING"
	#print 'kList: ' + str(kList)

	# for each original kmer, compute msps for each of the w's if they have a position in the genome
	for kmer in kList:
		kmerPositions = queryMap[kmer]
		for word in wDict[kmer]:
			#print 'kmer: ' + kmer +' word: ' + word
			wordPositions = queryMap[word]
			# search in each place the word happens in query
			for qPos in kmerPositions:
				# search in each place the word happens in the subject
				for sPos in database[word]:
					#print "getMSP at" + str((qPos,sPos))
					mspDict[kmer].append(getMSP(query, subject, qPos, sPos, k, S, score, threshold))

	# combine individual MSPs into a list, eliminate duplicates
	mspList = []
	for msps in mspDict.values():
		for msp in msps:
			if msp not in mspList and msp!=None:
				mspList.append(msp)

	return mspList


#getMSP(query, subject, qPos, sPos, k, S, score): starts with a match between query and subject defined by qPos and sPos.
# extends in both directions to find the MSP. returns the coordinates of the MSP for query and subject 
# if it scores above S. Stop searching if we fall below some percentage of the max score. 
# substitution matrix defined in score. Threshold for how far below max score can fall defined in threshold
# have a minimum length to extend as well? 
def getMSP(query, subject, qPos, sPos, k, S, score, threshold):
	#initialize dynamic programming array to length of query. This is the maximum space we can search over. 
	search = [[0 for i in range(len(query))] for i in range(len(query))]
	# quick compute initial score
	#match
	if query[qPos] == subject[sPos]: initalScore = score[0]
	#mismatch
	else: initalScore = score[1]

	#print "initalScore: " + str(initalScore)
	# set first number in score matrix
	search[0][0] = initalScore 
	maxScore = initalScore
	currentScore = initalScore
	maxPos = [0,0]

	#we can go both ways at the start
	extendRight = True
	extendLeft = True

	# current position we're extending to. 
	l = 0
	# while we can extend left and we're not at the start of the sequence
	while extendLeft and qPos-l>=0 and sPos-l>=0:
		extendRight = True
		#bases we've extended right
		r = maxPos[1]
		# while we can extend right and we're not at the end of the sequence.
		while extendRight and len(query)>qPos+r and len(subject)>sPos+r:
			#print 'l,r: ' + str((l,r))
			# extend to the right
			#pick next base in query and subject 
			qBase = query[qPos+r]
			sBase = subject[sPos+r]

			#match
			if qBase == sBase: newScore = search[l][r-1] + score[0]
			#mismatch
			else: newScore = search[l][r-1] + score[1]

			#make sure we dont fall too far below the max score
			if newScore > maxScore*threshold:
				# make sure we didn't update this value below... a bit of a hack
				if search[l][r] ==0:
					search[l][r] = newScore
					#print 'search['+str(l)+']['+str(r)+'] = ' + str(newScore)
				# if we've found a new maximum position, update
				if newScore > maxScore:
					maxScore=newScore
					maxPos = [l,r]
					#print "maxPos updated: " + str(maxPos)
					#print "maxScore: " + str(maxScore)
				#keep extending to the right	
				r+=1
			#if we fall below, stop extending
			else: extendRight = False
		
		# extend to the left if we're not at the start
		if qPos-1>0 and sPos-l>0:
			l += 1
			
			#reset r to 0
			r = maxPos[1]
			#print 'l,r: ' + str((l,r))
			#pick next base in query and subject 
			qBase = query[qPos-l]
			#print ' l qPos: ' +str(qPos)
			#print ' l qBase: ' +qBase
			sBase = subject[sPos-l]
			#print ' l sPos: '+ str(sPos)
			#print ' l sBase: '+ sBase
			#match
			if qBase == sBase: newScore = search[l-1][r] + score[0]
			#mismatch
			else: newScore = search[l-1][r] + score[1]
			#print ' l newScore: ' + str(newScore)
			#check if above threshold
			if newScore > maxScore*threshold:
				search[l][r] = newScore
				#print 'search['+str(l)+']['+str(r)+'] = ' + str(newScore)
				#try extending to the right to maximize
				extendRight = True
				if newScore > maxScore:
					maxScore=newScore
					maxPos = [l,r]
					#print "maxPos updated: " + str(maxPos)
					#print "maxScore: " + str(maxScore)
			# otherwise, we're done
			else: extendLeft = False
		# otherwise, we're done
		else: extendLeft = False

	# find the MSP
	# coordinates in query and subject
	qStart = qPos - maxPos[0]
	qEnd = qPos + maxPos[1]
	sStart = sPos - maxPos[0]
	sEnd = sPos + maxPos[1]
	#sequences
	qSeq = query[qStart:qEnd+1]
	sSeq = subject[sStart:sEnd+1]

	#for line in search:
	#	print line

	#only return if score >= S
	if maxScore >= S:
		return [qSeq, sSeq, (qStart,qEnd), (sStart, sEnd), maxScore]
	#else:
		#print "no MSP found with score >= "+str(S)

#generateMismatches(kList, mismatches):
# generateMismatches returns a dict form kList to words that will match with >T
# from there we can map each word back to where it came from in the genome.
def generateMismatches(kList, mismatches):
	misDict = dict((k, [k]) for k in kList)
	#do each round of mismatch generation
	while mismatches >0:
		# for each key in the dict...
		for originalK in misDict.keys():
			# make a copy of the list
			temp = misDict[originalK][:]
			for k in temp:
				# substitute at each base
				for base in range(len(k)):
					# substitute with each letter
					for letter in ['A','T','C','G']:
						sub = k[0:base]+letter+k[base+1:]
						# only add if not already there
						if sub not in misDict[originalK]:
							misDict[originalK].append(sub)
							#print 'append ' + sub
		mismatches += -1
	return misDict

#approximateLambda(S=[5,-4]): approximates the calculation of lambda based on the scoring matrix
# assumes equal nucleotide frequencies. 
# simply descend on to something near the value from above
# not exact, but simple. Doesn't take too long.
def approximateLambda(S=[5,-4]):
	# start with a grossly overestimated guess!
	lambdaGuess = 5.0
	prevGuess = lambdaGuess +1
	keepGuessing = True
	guesses = 0

	while keepGuessing:
		# fill in list of values to sum over
		valueList = []
		for i in range(4):
			for j in range(4):
				if i ==j:
					valueList.append((0.25*0.25)*math.pow(math.e, lambdaGuess*S[0]))
				else:
					valueList.append((0.25*0.25)*math.pow(math.e, lambdaGuess*S[1]))
		# descend toward real value
		if sum(valueList) > 1:
			lambdaGuess = lambdaGuess*0.9999
		else:
			keepGuessing=False

		if prevGuess < lambdaGuess or guesses > 10000000:
			keepGuessing = False
		guesses+=1
	return lambdaGuess

#doStats(lambdaStar, c, m, n, S, Kstar=None): computes blast statistics based on parameters
# reports Ps,m,n(c) and E for the given parameters
# These values will be returned as strings because they are in terms of Kstar
# Especially P, thats super annoying because its  a big polynomial 
# if Kstar = None, compute statistics in terms of Kstar
# if Kstar != None, computes exact statistics (UNIMPLEMENTED)
def doStats(lambdaStar, c, m, n, S, Kstar=None):
	assert Kstar==None, "Haven't coded exact statistics yet"

	#compute Ps,m,n(c)
	term1 = m * n * math.pow(math.e, -1*lambdaStar*S)
	# make a horribly, horribly ugly string to output P
	P = '1-exp(-'+str(term1)+'K) *['
	for i in range(c):
		if i ==0:
			P+=str(1)
		else:
			P +=('(('+str(term1)+'K^'+str(i)+')/'+str(math.factorial(i))+')')
		if i!=c-1:
			P+='+'
	P+=']'

	#compute E - same as other calculations above
	E = str(term1)+'K'
	
	return (P, E)



def main():
	print "running tests"
	#### TESTING 


	#TESTING getMSP
	query1 =   'GCGATGCAATTGC'
	subject1 = 'AAAATGCAATTAA'
	qPos1 = 5
	sPos1 = 5
	k = 4
	S = 10
	score = [5,-4]
	threshold = 0.1
	t1= getMSP(query1,subject1,qPos1,sPos1,k,S,score,threshold)
	assert t1 == ['ATGCAATT', 'ATGCAATT', (3, 10), (3, 10),40]

	#test match at start
	query2 = 'AGCTAGCT'
	subject2 = 'AGCTAGCT'
	qPos2=0
	sPos2=0
	t2= getMSP(query2,subject2,qPos2,sPos2,k,S,score,threshold)
	t2

	#TESTING blast
	# simple case with one obvious match
	query1 =   'GCGATGCAATTGC'
	subject1 = 'AAAATGCAATTAA'
	b1=blast(query1,subject1,None)
	assert b1==[['ATGCAATT', 'ATGCAATT', (3, 10), (3, 10), 40]]
	# match with a mismatch in the middle. lower S score
	query2 =   'GCGATGCGATTGC'
	subject2 = 'AAAATGCAATTAA'
	b2=blast(query2,subject2,None,S=30)
	assert b2==[['ATGCGATT', 'ATGCAATT', (3, 10), (3, 10), 31]]

	#match at start of sequence is giving error
	query21 = 'AGCTAGCT'
	subject21 = 'AGCTAGCT'
	b21 = blast(query21, subject21,None,S=20)

	# two matches in different locations
	query3 = 'ATTACGATCTCA'
	subject3 = 'ATTACTCAGTAGCTAGCTAGCGATCTCA'
	b3 = blast(query3, subject3,None,S=20)


if  __name__ =='__main__':main()
