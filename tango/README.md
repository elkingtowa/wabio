Tango: alignment-free analysis of mycobacteriophage genomes
=====

This repository contains all the code, data, figures and results for the alignment-free sequence analysis project of the Brown University 2014 Phage Hunters class. 

##Try it for yourself!
###Quick setup
1. Install [Anaconda](http://continuum.io/downloads) (which includes Python, scipy, numpy, and matplotlib, dependencies you'll want)
2. [Download](https://github.com/bsiranosian/tango/archive/analysis-only.zip) our scripts—or [download everything](https://github.com/bsiranosian/tango/archive/master.zip) instead.  
3. Jump to [configuration files](https://github.com/bsiranosian/tango/blob/master/README.md#configuration-files) for instructions on running the scripts.  

###Detailed installation instructions
We've published two main scripts to analyze kmer usage in mycobacteriphage genomes. With these you can reconstruct the main results from our poster and presentation: a nexus file of distances between phage that can be used to build a neighbor joining phylogenetic tree and plots of genomic self-similarity in a sliding window that can be used to look for horizontal gene transfer. 

To use these scripts, you need to have [python 2.7](https://www.python.org/downloads/) installed on your system, as well as the following packages:

- [scipy](http://www.scipy.org/install.html) --  _for efficient calculation of distance matrices_
- [numpy](http://www.numpy.org/) --  _for some tricky math (comes with scipy)_
- [matplotlib](http://matplotlib.org/downloads.html) --  _if you chose to produce plots with the compareTDI script_

The actual scripts are contained within the _kmer\_analysis_ folder. 

_compareTUD.py_ is used to calculate the kmer usage deviation across a number of phage. The default output is a nexus file containing the distances between phage that can be used in a program like splitstree. It can also output raw usage deviation values and a raw distance matrix.

_comparreTDI.py_ calculates genomic self-similarity in a sliding window across the genome for a number of phage. The default output is a file containing the z-scored deviation for each window in each phage. It can also produce a plot of the z-score across each genome for quick comparisons. 

_fastaDownloader.py_ is an accessory script used to download all the genomes from phagesdb.org. We've included the database of mycobacteriophage fasta files as of 2014-06-20 in this repo, but you can download them again or keep the folder up to date with this script.

###configuration files
The first argument to each script is a comma separated configuration file that contains at least the name and fasta file location for each phage that is going to be  compared. Each phage should be defined on a separate line. The format is as follows: 

    name,fasta location,(optional) subset start, (optional) subset end
For an example configuration file that has the names of the 60 phage used in Graham Hatfull's 2010 comparative phage genomics paper, look at _kmer\_analysis\\examples\\fasta\_map\_Hatfull60.txt_. Note that locations defined in this file are specific to my computer and will have to be changed. Fields 3 and 4 can be used to specify a subset of the seqence to be compared with start and end coordinates (zero indexed). 

###compareTUD.py
Your command prompt must be in the _kmer\_analysis_ folder to use this script, as I haven't made it into an installable package yet.
 
	minimal usage: python compareTUD.py fastaMap nexusFile

Where _fastaMap_ is the location of the configuration file and _nexusFile_ is the desired output location. 
 
	detailed usage: compareTUD.py [-h] [--subset subset] [--k k] [--s s] [--d d] [--r r]
                    fastaMap nexusFile

    This script calculates Tetranucleotide Usage Deviation (TUD) for multiple
    phage genomes. Originally it could only do 4-mers, but has now been
    generalized to all nucleotide lengths. The default function of this script is
    to save a nexus distance file that can be used for tree building with
    splitstree and other programs. Be aware that the nexus format is picky about
    special characters in the names of data, like parentheses. If you want to
    specify a cluster as part of the name, use a dash (ie "Dante-F1") The only
    required inputs are the location of a configuration file and location to save
    the resulting nexus distance file. The configuration file is a comma separated
    file of phage info with 2 necessary fields and 2 optional fields:
    name,fastaPath,subsetStart,subsetEnd Each phage to be compared should be
    separated by a new line. Additional arguments provide more control over the
    configuration the underlying calculations. If you want to save a csv file of
    the resulting usage deviation data, specify one with the --s option. A
    distance matrix can be saved by specifying a file name with the --d option.

    positional arguments:
      fastaMap         The file name of the comma separated file defining phage
                       information. The name and fastaPath fields are required,
                       subsets can be defined with integers in the next two
                       fields.
      nexusFile        The file to save the resulting nexus to.

    optional arguments:
      -h, --help       show this help message and exit
      --subset subset  set to True to plot only calculate for the regions defined
                       in the input file. If True, each line must have integers in
                       fields 3 and 4 that represent the genomic region of each
                       phage to compare. This allows you to subset on regions
                       containing genes, etc for each phage.
      --k k            Can also use this to compute 2,3,5-mers, etc.
      --s s            specify a filename to save resulting usage deviation data
                       to. Be careful because these files can get big for higher
                       values of k!
      --d d            specify a filename to save resulting distance matrix to
      --r r            By default sequences are extended by their reverse
                       complement before counting kmers and devation. Set to false
                       to change this behavior.

###compareTDI.py
Your command prompt must be in the _kmer\_analysis_ folder to use this script, as I haven't made it into an installable package yet.
 
	minimal usage: python compareTUD.py fastaMap saveFile

Where _fastaMap_ is the location of the configuration file and _nexusFile_ is the desired output location for comma separated deviation data.

    detailed usage: compareTDI.py [-h] [--subset subset] [--windowSize windowSize]
                    [--stepSize stepSize] [--k k] [--plotFile plotFile]
                    [--title title] [--maxNum maxNum] [--xScale xScale]
                    fastaMap saveFile

    This script computes the Tetranucleotide Difference Index (TDI) plot for
    multiple phage genomes. Originally it could only do 4-mers, but has now been
    generalized to all nucleotide lengths. The only required inputs are the file
    name of a configuration file and a file name to save the resulting data. The
    configuration file is a comma separated file of phage info with 2 necessary
    fields and 2 optional fields: name,fastaPath,subsetStart,subsetEnd Each phage
    to be compared should be separated by a new line. To save a simple plot of the
    resulting TDI Z-score across the genome, specify a file with the --plotFile
    option. You will need to have the matplotlib library available for import to
    use this option. Additional arguments provide more control over the
    configuration of the underlying calculations and the plot.

    positional arguments:
      fastaMap              The file name of the comma separated file defining
                            phage information. The name and fastaPath fields are
                            required, subsets can be defined with integers in the
                            next two fields.
      saveFile              specify a filename to save resulting data to. Format
                            will be CSV with the first row representing the start
                            of each window. Data for each phage will be plottd on
                            a new line. Z-scores for each window are the values in
                            the file.

    optional arguments:
      -h, --help            show this help message and exit
      --subset subset       set to True to plot only the regions defined in the
                            input file. If True, each line must have integers in
                            fields 3 and 4 that represent the genomic region of
                            each phage to compare
      --windowSize windowSize
                            The size of the window to compute TDI within. Default
                            of 5000bp is used in Pride et. al 2006
      --stepSize stepSize   How many bases to move the window along the genome at
                            each iteration. Default of 1000bp is used in Pride et.
                            al 2006
      --k k                 Can also use this to compute 2,3,5-mers, etc.
      --plotFile plotFile   The file to save the resulting image to.
      --title title         The title for the resulting plot.
      --maxNum maxNum       maximum number of phage to plot on the same figure.
                            first maxNum of input file will be chosen. default: 10
      --xScale xScale       Set to True to scale the x axis of the plot to
                            relative genome position. 

###fastaDownloader.py
This script takes in a tab delimited list of phage names and possibly clusters and gets fasta files from phagesdb.org for each one. We recommend using the "simple data download" from phagesdb to get this name information, available [here](http://phagesdb.org/data/?set=seq&type=simple). Beware that sometimes the names in this file don't correspond to the names of the fasta files on the database - sometimes there are small changes in capitalization and etc. We've taken this in to account for the current database, but the names of any failed downloads will be written to a _badNames.txt_ file in the output folder. 
	
	minimal usasge: python fastaDownloader.py phageData saveFolder

By default, the script reads in phage names from the file defined in phageData and saves resulting output files to saveFolder. To include the cluster names in the saved files, see the --parseClusters option in detailed usage.


	detailed usage: fastaDownloader.py [-h] [--parseClusters parseClusters]
    		                           [--header header]
             		                   phageData saveFolder

    This script is used to download fasta sequences from phagesdb.org. It takes as
    input a list of phage names downloaded from phagesdb.org and gets the
    corresponding fasta files from http://phagesdb.org/media/fastas/name.fasta.
    Sometimes, the names provided in the phagesdb data download file don't
    accurately match the names on the website. In this case, a file called
    "bad_names.txt" will be placed in the output directory describing the phage
    names that refused to download. All resulting files will be stored in the
    specified directory. Additional options can be used to include the name of the
    phage cluster in the filename.

    positional arguments:
      phageData             The file name of tab delimited phage data. We
                            recommend using the file downloaded from
                            http://phagesdb.org/data/?set=seq&type=simple. Each
                            phage must be defined on a new line. The first field
                            must contain the name of the phage (and therefore the
                            name of the fasta file on the website). The second
                            field optionally contains the phage cluster that can
                            be included in the saved file name with the
                            --parseClusters option
      saveFolder            All resulting files will be placed here. Will be
                            created if it does not exist.

    optional arguments:
      -h, --help            show this help message and exit
      --parseClusters parseClusters
                            Set to a character to have the phage cluster included
                            in the file name following this character. For
                            example, using "-" will save the resulting fasta as
                            "224-E.fasta"
      --header header       Number of lines to skip while reading phageData.
                            Change this only if you have a custom file


##Directories 
- *kmer_analysis* has the scripts that users should interact with
- *data* contains input/output data
- *src* contains code that is still in development and not ready for release. Lots of ideas and bugs here. Enter if you dare!
- *figures* contains pretty pictures, most of which are just automated outputs from scripts. 
- *presentation* contains files for the presentation. Chen has a writeup of his presentation procedure on his [website](http://yeesus.com/tangoSEA/).
- *writing* contains pretty abstract prose. And paper stuff eventually too.

##Changelog

_2014-06-20_
Getting ready for the alpha release! Spicing up this readme... organizing files, explaining how people can interact with the two main scripts

_2014-06-16_
The second coming.  

_2014-06-15—2014-04-18_
Jeez I have no idea.  Review the commit log.  Actually please don't. 

_2014-04-17_
Added inital code for TUD calculation and figure generation. Two output files are created from the data in the Hatfull et. al 2010 comparative phage genomics paper. This subset of phages will be used for inital tests of our algorithms. 
