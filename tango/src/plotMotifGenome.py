from tetranucleotideAnalysis import *
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import numpy as np

#countMotif(sequence, motif, windowSize, stepSize): counts the occurances of a given motif in windows of windowSize, moving along sequence in intervals of stepSize
# returns a list of tuple: (position of start of region, occurances of the given motif for region) 
def countMotif(sequence, motif, windowSize, stepSize):
	start = 0
	end = start+windowSize
	counts =[]
	while end < len(sequence):
		window=sequence[start:end]
		#find all with regexp 
		reg =r'(?=('+motif+'))'
		count = len(re.findall(reg, window))
		counts.append((start, count))
		start += stepSize
		end += stepSize
	return counts


# code to plot the frequency of a given motif across the genome of multiple phage
def plotMotifGenome(nameFile, motif, windowSize, stepSize, title, saveName, maxNum, saveData=None, RC=False):
	# read information in nameFile
	names = []
	fnames = []
	with open(nameFile, 'r') as nf:
		line = nf.readline().strip().split('\t')
		while line != ['']:
			names.append(line[0])
			fnames.append(line[1])
			line = nf.readline().strip().split('\t')
	# do calculations for each sequence
	data=[]
	for fname, name in zip(fnames, names):
			#parse fasta from filename
		for seq_record in SeqIO.parse(fname, "fasta"):
			sequence = seq_record.seq.tostring().upper()

		if not RC:
			print "Computing for " + name
			d = countMotif(sequence, motif, windowSize, stepSize)
			data.append(d)

		## NEW CODE ## 
		# add data from counting across RC if RC=True
		if RC:
			rc = reverseComplement(sequence)
			f = countMotif(sequence, motif, windowSize, stepSize)
			r = countMotif(rc, motif, windowSize, stepSize)
			d = [(a[0],a[1]+b[1]) for a,b in zip(f,r)]
			data.append(d)

	# get from data into lists across genome
	ax = plt.subplot(1,1,1)
	# Code from SO to make colors distinguishable
	NUM_COLORS = maxNum
	cm = plt.get_cmap('gist_rainbow')
	cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
	scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
	ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(NUM_COLORS)])

	i = 0
	maxXaxis = []
	for phage,name in zip(data, names):
		xaxis = [start for (start, count) in phage]
		# find largest set of bins to save data with
		if len(xaxis) > len(maxXaxis):
			maxXaxis = xaxis
		counts  = [count for (start, count) in phage]
		if i < maxNum:
			color = cm(1.*i/NUM_COLORS)
			plt.plot(xaxis, counts, lw=1, label=name)
			i+=1
	print maxXaxis
	handles, labels = ax.get_legend_handles_labels()
	plt.title(title)
	plt.xlabel('genomic position')
	plt.ylabel('count of motif: '+ motif)
	plt.legend(prop={'size':5})
	plt.savefig(saveName, dpi=300)
	plt.clf()
	print "Plot \' " + title+ "\' saved to " + saveName 
	
	if saveData != None:
		with open(saveData, 'w') as of:
			#save data to output file
			# header of motif and genomic positions
			header = motif
			for pos in maxXaxis:
				header += ','+str(pos)
			of.write(header+'\n')
			for phage,name in zip(data, names):
				counts  = [count for (start, count) in phage]
				toWrite = name
				for i in range(len(maxXaxis)):
					if i < len(counts):
						toWrite += ','+str(counts[i])
					else:
						toWrite += ','

				of.write(toWrite+ '\n')
	print "Done  :)"


#plot this cool GATC motif for some clusters
# B3 is extremely elevated. look in comparison to a few others 
import os
# os.chdir('C:/Users/Admin/Documents/GitHub/tango/src')
# plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_B3.txt' , 'GATC', 5000,2500,'GATC frequency in B3 genomes. 1000/500', '../figures/without_reverse_complement/motif_plots/GATC_motif_B3.png',10,saveData='../data/without_reverse_complement/GATC_motif_B3.tsv',RC=False)
# plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_B2.txt' , 'GATC', 5000,2500,'GATC frequency in B2 genomes. 1000/500', '../figures/without_reverse_complement/motif_plots/GATC_motif_B2.png',10,saveData='../data/without_reverse_complement/GATC_motif_B2.tsv',RC=False)
# plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_B1.txt' , 'GATC', 5000,2500,'GATC frequency in B1 genomes. 1000/500', '../figures/without_reverse_complement/motif_plots/GATC_motif_B1.png',10,saveData='../data/without_reverse_complement/GATC_motif_B1.tsv',RC=False)

os.chdir('C:/Users/Admin/Documents/GitHub/tango/src')
plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_B3.txt' , 'GGATCC', 5000,2500,'GGATCC frequency in B3 genomes. 1000/500', '../figures/without_reverse_complement/motif_plots/GGATCC_motif_B3.png',10,saveData='../data/without_reverse_complement/GGATCC_motif_B3.tsv',RC=False)
plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_B2.txt' , 'GGATCC', 5000,2500,'GGATCC frequency in B2 genomes. 1000/500', '../figures/without_reverse_complement/motif_plots/GGATCC_motif_B2.png',10,saveData='../data/without_reverse_complement/GGATCC_motif_B2.tsv',RC=False)
plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_B1.txt' , 'GGATCC', 5000,2500,'GGATCC frequency in B1 genomes. 1000/500', '../figures/without_reverse_complement/motif_plots/GGATCC_motif_B1.png',10,saveData='../data/without_reverse_complement/GGATCC_motif_B1.tsv',RC=False)

# #cluster G have lots of TCGA. Compare to A1, B1
# plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_G.txt' , 'TCGA', 1000,900,'TCGA frequency in G genomes. 1000/900', '../figures/with_reverse_complement/motif_plots/TCGA_motif_G.png',10,saveData='../data/with_reverse_complement/TCGA_motif_G.tsv',RC=True)
# plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_A1.txt' , 'TCGA', 1000,900,'TCGA frequency in A1 genomes. 1000/900', '../figures/with_reverse_complement/motif_plots/TCGA_motif_A1.png',10,saveData='../data/with_reverse_complement/TCGA_motif_A1.tsv',RC=True)
# plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_B1.txt' , 'TCGA', 1000,900,'TCGA frequency in B1 genomes. 1000/900', '../figures/with_reverse_complement/motif_plots/TCGA_motif_B1.png',10,saveData='../data/with_reverse_complement/TCGA_motif_B1.tsv',RC=True)

# #cluter L have lots of CCTA. Compare to some others 
# plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_L1.txt' , 'CCTA', 1000,200,'CCTA frequency in L1 genomes. 1000/200', '../figures/with_reverse_complement/motif_plots/CCTA_motif_L1.png',10,saveData='../data/with_reverse_complement/CCTA_motif_L1.tsv',RC=True)
# plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_L2.txt' , 'CCTA', 1000,200,'CCTA frequency in L2 genomes. 1000/200', '../figures/with_reverse_complement/motif_plots/CCTA_motif_L2.png',10,saveData='../data/with_reverse_complement/CCTA_motif_L2.tsv',RC=True)
# plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_L3.txt' , 'CCTA', 1000,200,'CCTA frequency in L3 genomes. 1000/200', '../figures/with_reverse_complement/motif_plots/CCTA_motif_L3.png',10,saveData='../data/with_reverse_complement/CCTA_motif_L3.tsv',RC=True)
# plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_A1.txt' , 'CCTA', 1000,200,'CCTA frequency in A1 genomes. 1000/200', '../figures/with_reverse_complement/motif_plots/CCTA_motif_A1.png',10,saveData='../data/with_reverse_complement/CCTA_motif_A1.tsv',RC=True)
# plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_B1.txt' , 'CCTA', 1000,200,'CCTA frequency in B1 genomes. 1000/200', '../figures/with_reverse_complement/motif_plots/CCTA_motif_B1.png',10,saveData='../data/with_reverse_complement/CCTA_motif_B1.tsv',RC=True)
# #cluter L have lots of CTTA. Compare to some others 
# plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_L1.txt' , 'CTTA', 1000,200,'CTTA frequency in L1 genomes. 1000/200', '../figures/with_reverse_complement/motif_plots/CTTA_motif_L1.png',10,saveData='../data/with_reverse_complement/CTTA_motif_L1.tsv',RC=True)
# plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_L2.txt' , 'CTTA', 1000,200,'CTTA frequency in L2 genomes. 1000/200', '../figures/with_reverse_complement/motif_plots/CTTA_motif_L2.png',10,saveData='../data/with_reverse_complement/CTTA_motif_L2.tsv',RC=True)
# plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_L3.txt' , 'CTTA', 1000,200,'CTTA frequency in L3 genomes. 1000/200', '../figures/with_reverse_complement/motif_plots/CTTA_motif_L3.png',10,saveData='../data/with_reverse_complement/CTTA_motif_L3.tsv',RC=True)
# plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_A1.txt' , 'CTTA', 1000,200,'CTTA frequency in A1 genomes. 1000/200', '../figures/with_reverse_complement/motif_plots/CTTA_motif_A1.png',10,saveData='../data/with_reverse_complement/CTTA_motif_A1.tsv',RC=True)
# plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_B1.txt' , 'CTTA', 1000,200,'CTTA frequency in B1 genomes. 1000/200', '../figures/with_reverse_complement/motif_plots/CTTA_motif_B1.png',10,saveData='../data/with_reverse_complement/CTTA_motif_B1.tsv',RC=True)
