import argparse

parser = argparse.ArgumentParser(description="This script is used to download fasta sequences from phagesdb.org. It takes as input a list of phage names downloaded from phagesdb.org and gets the corresponding fasta files from http://phagesdb.org/media/fastas/name.fasta. Sometimes, the names provided in the phagesdb data download file don't accurately match the names on the website. In this case, a file called \"bad_names.txt\" will be placed in the output directory describing the phage names that refused to download. All resulting files will be stored in the specified directory. Additional options can be used to include the name of the phage cluster in the filename.")
parser.add_argument("phageData", action="store", metavar="phageData", help="The file name of tab delimited phage data. We recommend using the file downloaded from http://phagesdb.org/data/?set=seq&type=simple. Each phage must be defined on a new line. The first field must contain the name of the phage (and therefore the name of the fasta file on the website). The second field optionally contains the phage cluster that can be included in the saved file name with the --parseClusters option")
parser.add_argument("saveFolder", action="store", metavar="saveFolder", help="All resulting files will be placed here. Will be created if it does not exist.")
parser.add_argument("--parseClusters", action="store", metavar="parseClusters", default="", help="Set to a character to have the phage cluster included in the file name following this character. For example, using \"-\" will save the resulting fasta as \"224-E.fasta\"")
parser.add_argument("--header", action="store", metavar="header", default="1", help="Number of lines to skip while reading phageData. Change this only if you have a custom file")
args=parser.parse_args()

phageData=args.phageData
saveFolder=args.saveFolder
parseClusters=args.parseClusters
header=int(args.header)
## DONE PARSING ARGS ##


import requests
import os

#downloads all fasta files from phagesDB 
#process outline:
#parses the list of phage names downloaded from phagesDB
#for each name, attempts to download once
#if the request is successful (status_code != 404), go on to the next name
#else, try other stuff 
def fastaDownloader(phageData, saveFolder, parseClusters, header):
	# read data from phageData
	names = []
	clusters = []
	with open(phageData, 'r') as pd:
		# skip header lines and end with first to read
		for i in range(header+1):
			line = pd.readline().strip().split('\t')

		#loop through doc
		while line != ['']:
			names.append(line[0])
			if parseClusters:
				clusters.append(line[1])
			line = pd.readline().strip().split('\t')

	# define where we're going to get fastas from 
	baseURL = 'http://phagesdb.org/media/fastas/'
	endURL = '.fasta'

	# make output dir if it doesn't exist
	if not os.path.exists(saveFolder):
		os.makedirs(saveFolder)

	# cheat to make iterating easier
	if len(names) != len(clusters):
		clusters = ['' for n in names]

	# begin downloading phage
	attempts = 0
	successfullyDownloaded = 0
	failedDownloads = 0
	notAttempted = 0
	badNames = []
	for name,cluster in zip(names, clusters):
		attempts += 1

		# WE KNOW SOME NAMES ARE BAD ALREADY
		if name == "GUmbie":
			name = "Gumbie"
		elif name == "Numberten":
			name = "NumberTen"
		elif name == "Seabiscuit":
			name = "SeaBiscuit"
		elif name == "Temp1":
			name = "BabyRatchet"

		#check to see if file already exists. If so, don't download
		if not os.path.isfile(os.path.join(saveFolder,name+endURL)):
			# def url
			url = baseURL + name + endURL
			# look for file
			fileRequest = requests.get(url)

			# display  error code
			if (fileRequest.status_code != 200):
				print "phage " + name + " did not work, HTTP error code: " + str(fileRequest.status_code)
				badNames.append(name)
				failedDownloads += 1
			#successful download
			else: 
				print "downloading fasta file for phage " + name
				with open(os.path.join(saveFolder,name+parseClusters+cluster+'.fasta'), 'wb') as ff:
					ff.write(fileRequest.content)
				successfullyDownloaded += 1
		# file exists, don't download
		else: 
			print name+endURL + ' exists, not downloading'
			notAttempted += 1

	# write badnames to file
	with open(os.path.join(saveFolder, 'badNames.txt'),'wb') as bn:
		for name in badNames:
			bn.write(name+'\n')

	#displays results
	print "results: \n"
	print str(attempts) + " downloads attempted"
	print str(successfullyDownloaded) + " fasta files downloaded successfully"
	print str(failedDownloads) + " downloads failed"
	print str(notAttempted) + " downloads not attempted"

if __name__ == '__main__':
	fastaDownloader(phageData, saveFolder, parseClusters, header)
