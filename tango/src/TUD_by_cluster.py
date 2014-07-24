import csv

#computes the mean TUD for each sucluster, and compiles it into a TSV file. 

#unused
clusterNames = ["(A1)", "(A2)", "(A3)", "(A4)", "(A5)", "(A6)", "(A7)", "(A8)", "(A9)", "(A10)", "(A11)",
"(B1)", "(B2)", "(B3)", "(B4)", "(B5)", "(C1)", "(C2)"]

#scans through the name of each phage until it gets to the cluster name in parenthesis. 
#if the name is not found in dict.keys() then add it to the dictionary with an array containing the entire line.
#other lines are appended to the array containing the lines. 

path = '/Users/edwardwilliams/Documents/tango/data/all_phages_TUD.tsv'

clusterDictionary = {}
firstline = ""

#sort by subclusters
with open(path, 'r') as phagefile:
	for line in phagefile:
		#remove last parenthesis for ease of use
		lineArray = line.split("	")
		name = lineArray[0][:-1]
		#print name
		cluster = name.split("(")
		if len(cluster) > 1: #if not the first line
			subcluster = cluster[1]
			if subcluster == "Singleton": #add it to the dictionary with its name
				clusterDictionary[lineArray[0]] = lineArray
			elif subcluster in clusterDictionary.keys():
				oldarray = clusterDictionary[subcluster]
				newarray = oldarray.append(lineArray)
				clusterDictionary[subcluster]
			else:
				clusterDictionary[subcluster] = [lineArray]
		else:
			#saves the first line to write it to the output file
			firstLine = lineArray

#find average of each subcluster
#each array has 257 elements, 256 of which are TUD values

newfile = '/Users/edwardwilliams/Documents/tango/data/mean_tud_by_subcluster.csv'
with open(newfile, 'w') as outfile:
	csvwriter = csv.writer(outfile)
	csvwriter.writerow(firstLine)

	for subcluster in clusterDictionary.keys():

		#checks if it's a singleton (if its first element is a string and not an array)

		data = clusterDictionary[subcluster]

		if type(data[0]) is str:
			csvwriter.writerow(data)

		else:
			meantud = [subcluster]

			for i in range(1, 257):
				#sum TUD of each phage in this cluster for this tetranucleotide
				tudsum = 0
				phagecounter = 0
				for phage in data:
					phagecounter +=1
					tudsum += float(phage[i])
				mean = tudsum/phagecounter
				meantud.append(mean)
			csvwriter.writerow(meantud)



