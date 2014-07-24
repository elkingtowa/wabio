import csv
import scipy.stats

#parses through file, converts the TUD measurements into an array of arrays

tudarray = []
subclusterdict = {}
clusterdict= {}
filepath = "/Users/edwardwilliams/Documents/tango/data/without_reverse_complement/all_phages_TUD.tsv"


with open(filepath, 'r') as tudfile:
	for line in tudfile:
		tudarray.append(line.split("	")) #splits tab separated values into array


#double for loop goes through every possible combination of phage, finds spearman correlation, saves to dictionaries sorted by cluster

totalwithsingleton = 0
corrnumberwithsingleton = 0
totalnosingleton = 0
corrnumbernosingleton = 0

for i in range(1, len(tudarray)):
	for j in range(1, len(tudarray)):
		corr = scipy.stats.spearmanr(tudarray[i][1:], tudarray[j][1:])
		totalwithsingleton += corr[0]
		corrnumberwithsingleton += 1
		phageone = tudarray[i][0]
		phagetwo = tudarray[j][0]
		#finding cluster of phage
		c1 = phageone.split("(")[1]
		clusterone = c1[:-1]
		c2 = phagetwo.split("(")[1]
		clustertwo = c2[:-1]
		print clusterone, clustertwo, corr
		#check if cluster is the same

		if clusterone != "Single" and clusterone !="Singleton" and clustertwo != "Single" and clustertwo != "Singleton": #excluding singleton phage
			totalnosingleton += corr[0]
			corrnumbernosingleton += 1

			if clusterone[0] == clustertwo[0] and clusterone != "Single" and clusterone != "Singleton":
				if clusterone[0] in clusterdict.keys():
					oldarray = clusterdict[clusterone[0]]
					newarray = oldarray.append(corr)
					clusterdict[clusterone[0]]
				else:
					clusterdict[clusterone[0]] = [corr]
				#comparing subcluster
				if (clusterone == clustertwo):
					if clusterone in subclusterdict.keys():
						oldarray = subclusterdict[clusterone]
						newarray = oldarray.append(corr)
						subclusterdict[clusterone]
					else:
						subclusterdict[clusterone] = [corr]


meancorrsingleton = totalwithsingleton/corrnumberwithsingleton

print "Mean correlation including singleton phage: " + str(meancorrsingleton)

meancorrnosingleton = totalnosingleton/corrnumbernosingleton

print "Mean correlation excluding singleton phage: " + str(meancorrnosingleton)

#finding average correlation vector between subclusters
mcorrsubdict = {}
mcorrclusdict= {}

for subcluster in subclusterdict.keys():
	corrarray = subclusterdict[subcluster]
	total = 0

	for corrvector in corrarray:
		total += corrvector[0]

	mean = total/len(corrarray)
	mcorrsubdict[subcluster] = mean

print mcorrsubdict

subclustotal = 0

for subcluster in mcorrsubdict.keys():
	corr = mcorrsubdict[subcluster]
	subclustotal += corr

meancorrsubcluster = subclustotal/len(mcorrsubdict)

print "Mean correlation between subclusters: " + str(meancorrsubcluster)

for cluster in clusterdict.keys():
	corrarray = clusterdict[cluster]
	total = 0

	for corrvector in corrarray:
		total += corrvector[0]

	mean = total/len(corrarray)
	mcorrclusdict[cluster] = mean

print mcorrclusdict

clustotal = 0

for cluster in mcorrclusdict.keys():
	corr = mcorrclusdict[cluster]
	clustotal += corr

meancorrcluster = clustotal/len(mcorrclusdict)

print "Mean correlation between clusters: " + str(meancorrcluster)
