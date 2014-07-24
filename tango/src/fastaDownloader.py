#downloads all fasta files from phagesDB 
#process outline:
#parses the list of phage names downloaded from phagesDB
#for each name, attempts to download once
#if the request is successful (status_code != 404), go on to the next name
#else, try other stuff 

import requests
import os.path

nameList = open("sequenced_phage_clusters.txt")

#need to append phageName and '.fasta'
urlStart = 'http://phagesdb.org/media/fastas/'

#text file where uncooperative results will be stored
badNames = open("badNames.txt", "wb")

attempts = 0
successfullyDownloaded = 0
failedDownloads = 0

for unformattedName in nameList:
	attempts += 1

	#parses string until the first space 
	phageName = unformattedName.strip().split('\t')[0]
	fileName = phageName + '.fasta'

	#checks if file already exists
	isFile = os.path.isfile("results/" + fileName)

	if not isFile:
		url = urlStart + fileName
		fileRequest = requests.get(url)

		if (fileRequest.status_code != 200):
			string = "phage " + phageName + " did not work, HTTP error code: " + str(fileRequest.status_code)
			badNames.write(string + '\n')
			print string
			failedDownloads += 1
		else:
			#download file to directory /results
			successfullyDownloaded += 1
			print "downloading fasta file for phage " + phageName
			with open("results/" + fileName, 'wb') as outputFile:
				outputFile.write(fileRequest.content)
	else:
		print "phage " + phageName + " was downloaded previously"
		successfullyDownloaded +=1


#displays results
print "results: \n"
print str(attempts) + " downloads attempted"
print str(successfullyDownloaded) + " fasta files downloaded successfully"
print str(failedDownloads) + " downloads failed"