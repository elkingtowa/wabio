README		BEN SIRANOSIAN		BLAST
-------------------------------------
Problem 1: BLAST

The user should interface with two python scripts for this project. doBlast.py will be used for all blast functions. generateDatabase.py can be used to save a database (in pickle format) of the subject sequence so that it doesn't have to be computed every time. 

USAGE for generateDatabase.py
python generateDatabase.py subject k outfile 
PARAMETERS
subject: the location of a fasta file containing the subject sequence
k: integer word size to break up the sequence and compute match locations
outfile: filename to store the resulting database

USAGE for doBlast.py
python doBlast.py query subject
REQUIRED PARAMETERS
query: the location of a fasta file containing the query sequence
subject: the location of a fasta file containing the subject sequence
OPTIONAL PARAMETERS
--databaseFile: you can specify a databse that has been saved with generateDatabase.py with this option. If unspecified, doBlast will create a database in memory each time it runs. 
--k: Word length to break up the query string. Default is 4. must match the k used to compute the database
--T: Minimum score to create search terms after breking query sequence. Default is 10
--S: Minimum score to report MSPs. You can specify a fixed value here, or use the default which is (0.6*5*|qurey|)
--threshold: Keep extending for matches until score falls below this fraction of max score. default is 0.75
 --substitutionMatrix: Substitution matrix for scoring matches. Must be specified as a python formatted list, where the first number is the score for a match and the second number is the score for a mismatch. Default is [5,-4]


RESULTS FOR BLASTING BRCA1 EXONS
MOUSE
Query	0		GATGCGGAGTTTGTGTGTGAGCGGACACTGAAATATTTTCTGGGCATTGCAGGAGGAAAGTGGATAGTTAGCTATTCATG	79
				||||| |||||||||||||| |||||||||||||||||||| || ||||| |||| ||| ||| ||| ||||||||  ||
Sbjct	12042	GATGCTGAGTTTGTGTGTGAACGGACACTGAAATATTTTCTCGGAATTGCGGGAGCAAAATGGGTAGATAGCTATTTCTG	12121
SCORE: 301
E = 1.00426372836e-19K
Ps,m,n(1) = 1-exp(-1.00426372836e-19K) *[1]

DOG
Query	3	GTCTGGAGTTGATCAAAGAGCCTGTTTCTACAAAGTGTGATCACATATTTTGCAAGT	59
			|||||||||||||||| || ||||| || ||||||||||| |||||| ||||||| |
Sbjct	405	GTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATAATTTGCAAAT	461
SCORE: 222
E = 2.85263705761e-13K
Ps,m,n(1) = 1-exp(-2.85263705761e-13K) *[1]

COW
Query	1		CAGATCCTGAGTTTGTGTGTGAACGGACACTGAAATATTTCCTGGGAATTGCAGGAGGAAAATGGGTAGTTAGCTATTTTT	81
				||||| |||||||||||||||||||||||||||||||||| || |||||||| |||| ||||||||||| ||||||||| |
Sbjct	12040	CAGATGCTGAGTTTGTGTGTGAACGGACACTGAAATATTTTCTCGGAATTGCGGGAGCAAAATGGGTAGATAGCTATTTCT	12120
SCORE: 342
E = 4.00150492466e-23K
Ps,m,n(1) = 1-exp(-4.00150492466e-23K) *[1]

RAT
Query	0	GTTTGGAACTGATCAAAGAACCGGTTTCCACAAAGTGCGACCACATATTTTGCAA	54
			|| ||||  ||||||| ||||| || ||||||||||| ||||||||| |||||||
Sbjct	405	GTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATAATTTGCAA	459
SCORE: 203
E = 9.96575062155e-12K
Ps,m,n(1) = 1-exp(-9.96575062155e-12K) *[1]


BLAST COMPARED TO SW LOCAL ALIGNMENT ON RAT EXON
BLAST:
time python doBlast.py ../data/BRCA1-exons4.fasta ../data/BRCA1-dna.fasta --databaseFile ../data/BRCA_database.pkl
real	0m5.237s
user	0m5.223s
sys	0m0.008s

LOCAL ALIGNMENT:
time python local_alignment.py ../data/BRCA1-exons4.fasta ../data/BRCA1-dna.fasta 
['GTTTGGAACTGATCAAAGAACCGGTTTCCACAAAGTGCGACCACATATTTTGCAA', 'GTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATAATTTGCAA']

real	0m0.765s
user	0m0.756s
sys	0m0.011s

WAIT, so local alignmet is much faster. That isn't good. I think there's something funny in the code that is making BLAST take way too long. 