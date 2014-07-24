#parameters 
import argparse

parser = argparse.ArgumentParser(description="This script computes the Tetranucleotide Difference Index (TDI) plot for multiple phage genomes. Originally it could only do 4-mers, but has now been generalized to all nucleotide lengths. The only required inputs are the file name of a configuration file and a file name to save the resulting data. The configuration file is a comma separated file of phage info with 2 necessary fields and 2 optional fields: \n name,fastaPath,subsetStart,subsetEnd\nEach phage to be compared should be separated by a new line. To save a simple plot of the resulting TDI Z-score across the genome, specify a file with the --plotFile option. You will need to have the matplotlib library available for import to use this option. Additional arguments provide more control over the configuration of the underlying calculations and the plot.")
parser.add_argument("fastaMap", action="store", metavar="fastaMap", help="The file name of the comma separated file defining phage information. The name and fastaPath fields are required, subsets can be defined with integers in the next two fields.")
parser.add_argument("dataFile", action="store", metavar="dataFile", help="specify a filename to save resulting data to. Format will be CSV with the first row representing the start of each window. Data for each phage will be plottd on a new line. Z-scores for each window are the values in the file.")
parser.add_argument("--subset", action="store", metavar="subset", default="False", help="set to True to plot only the regions defined in the input file. If True, each line must have integers in fields 3 and 4 that represent the genomic region of each phage to compare")
parser.add_argument("--windowSize", action="store", metavar="windowSize", default="5000", help="The size of the window to compute TDI within. Default of 5000bp is used in Pride et. al 2006")
parser.add_argument("--stepSize", action="store", metavar="stepSize", default="1000", help="How many bases to move the window along the genome at each iteration. Default of 1000bp is used in Pride et. al 2006")
parser.add_argument("--k", action="store", metavar="k", default="4", help="Can also use this to compute 2,3,5-mers, etc.")
parser.add_argument("--plotFile", action="store", metavar="plotFile", default='', help="The file to save the resulting image to.")
parser.add_argument("--title", action="store", metavar="title", default='', help="The title for the resulting plot.")
parser.add_argument("--maxNum", action="store", metavar="maxNum", default="10", help="maximum number of phage to plot on the same figure. first maxNum of input file will be chosen. default: 10")
parser.add_argument("--xScale", action="store", metavar="xScale", default="False", help="Set to True to scale the x axis of the plot to relative genome position.")
args=parser.parse_args()

fastaMap=args.fastaMap
plotFile=args.plotFile
title=args.title
maxNum=int(args.maxNum)
if args.xScale == "True": xScale=True
else: xScale=False
if args.subset == "True": subset=True
else: subset=False
windowSize=int(args.windowSize)
stepSize=int(args.stepSize)
k=int(args.k)
dataFile = args.dataFile

from kmerAnalysis import doTDI
import numpy as np
import csv

#compareTDI(): computes tetranucleotide difference index z-scores for all genomes defined in fastaMap. 
# takes as input a  comma separated file of phage info with 2 necessary fields and 2 optional fields:
# 	name ,fastaPath,subsetStart,subsetEnd
# if xSacle is True, x axis is scaled to relative position along the genome to allow comparrison of different length sequences. 
# if subset is true, each sequence is shortened to the positions defined in the final two fields of the input.
# other arguments (windowSize, stepSize, k) are passed to the TDI function. 
# if plotFile != '', save a plot of the resulting data. Matplotlib must be installed.
# plot titled with title. saves the resulting plot to plotFile
# only first maxNum lines are plotted. the first 10 phage defined in the input file are plotted by default. 
def compareTDI(fastaMap, dataFile, xScale, subset, windowSize, stepSize, k, plotFile, title, maxNum):
	# read information in fastaMap
	names = []
	fnames = []
	subsets = []
	with open(fastaMap, 'r') as nf:
		line = nf.readline().strip().split(',')
		while line != ['']:
			if subset and (len(line) < 4):
				print "Specified to subset but didn't include start and end position"
				print "Subsetting will be disabled"
				subset = False
			names.append(line[0])
			fnames.append(line[1])
			# parse subset information
			if subset:
				subsets.append([line[2],line[3]])
			line = nf.readline().strip().split(',')
	# compute for each defined phage
	data = []
	if subset:
		calcNum=0
		for fname, name, subset in zip(fnames, named, subsets):
			if calcNum<maxNum:
				print "Computing for " + name
				d = TDI(fname, windowSize, stepSize, k,  subset=subset)
				data.append(d)
			calcNum+=1
	else:
		calcNum=0
		for fname, name in zip(fnames, names):
			if calcNum<maxNum:
				print "Computing for " + name
				d = doTDI(fname, windowSize, stepSize, k)
				data.append(d)
			calcNum+=1

	# Save resulting data
	with open(dataFile, 'wb') as sf:
		writer = csv.writer(sf)
		#find longest window list to write as header 
		xs = [d[0] for d in data]
		header = max(xs, key=len)
		writer.writerow(['name']+header)
		for row, name  in zip(data, names):
			writer.writerow([name]+row[1])

	#if we're plotting...
	if plotFile != '':
		import matplotlib.pyplot as plt
		import matplotlib.cm as mplcm
		import matplotlib.colors as colors

		#if xScale, rescale axes and plot on genome position scale
		if xScale:
			for i in range(len(data)):
				scaledX = []
				for j in range(len(data[i][0])):
					scaledX.append(data[i][0][j] / float(data[i][0][len(data[i][0])-1]))
				data[i][0] = scaledX
		
		ax = plt.subplot(1,1,1)
		# make colors distinguishable
		if maxNum > len(names): numColors = len(names)
		else: numColors = maxNum
		cm = plt.get_cmap('gist_rainbow')
		cNorm  = colors.Normalize(vmin=0, vmax=numColors-1)
		scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
		ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(numColors)])


		i = 0
		for [xaxis, Zscores],name in zip(data, names):
			color = cm(1.*i/numColors)
			plt.plot(xaxis, Zscores, lw=1, label=name)
			i+=1

		handles, labels = ax.get_legend_handles_labels()
		plt.title(title)
		plt.xlabel('genomic position')
		plt.ylabel('TDI Z-score')
		plt.legend(prop={'size':5})
		plt.savefig(plotFile, dpi=300)
		plt.clf()
		print "Plot \' " + title+ "\' saved to " + plotFile 



	print "Done  :)"

# function call
if __name__ == '__main__':
	compareTDI(fastaMap, dataFile, xScale, subset,windowSize,stepSize,k,plotFile,title,maxNum)