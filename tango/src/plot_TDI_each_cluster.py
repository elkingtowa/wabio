clusterNames = ["A1","A10","A11","A2","A3","A4","A5","A6","A7","A8","A9",
				"B1","B2","B3","B4","B5","C1","C2","D1","D2","E","F1","F2","F3",
				"G","H1","H2","I1","I2","J","K1","K2","K3","K4","K5","L1","L2","L3",
				"M1","M2","N","O","P","Q","R","S","Singleton","T"]

filemap = "C:/Users/Admin/Documents/GitHub/tango/data/phagesDB/sequenced_phage_map.txt"
outputLoc = "C:/Users/Admin/Documents/GitHub/tango/data/TDI_individual_clusters/"
with open(filemap, 'r') as fm:
	lines=fm.readlines()

clusters=[]
for i in range(len(lines)):
	lines[i] = lines[i].strip()
	clusters.append(lines[i].split('\t')[0].split('(')[1][:-1])

cmap = {key:[] for key in clusterNames}
for i in range(len(lines)):
	cmap[clusters[i]].append(lines[i])

for key, values in cmap.items():
	outfile= outputLoc+"sequenced_phage_map_"+key+".txt"
	with open(outfile, 'w') as of:
		for value in values:
			of.write(value+'\n')

with open(outputLoc+"oscar_commands.sh","w") as of:
	for cluster in clusterNames:
		of.write("python compareTDI.py sequenced_phage_map_"+cluster+".txt ../../figures/TDI_individual_clusters/TDI_"+cluster+" \"TDI first 10 phage of "+cluster+'\"\n')