import sys
#### CONFIGURATION - CHANGE THESE VALUES TO WHAT YOU DESIRE!  ####
## DEFINE TRANSITION MATRIX
tmat = [[0.180,0.274,0.426,0.120,0.00011,0.00011,0.00011,0.00011],
        [0.171,0.368,0.274,0.188,0.00011,0.00011,0.00011,0.00011],
        [0.161,0.339,0.375,0.125,0.00011,0.00011,0.00011,0.00011],
        [0.079,0.355,0.384,0.182,0.00011,0.00011,0.00011,0.00011],
        [0.000044,0.000044,0.000044,0.000044,0.300,0.205,0.285,0.210],
        [0.000044,0.000044,0.000044,0.000044,0.322,0.298,0.078,0.302],
        [0.000044,0.000044,0.000044,0.000044,0.248,0.246,0.298,0.208],
        [0.000044,0.000044,0.000044,0.000044,0.177,0.239,0.292,0.292]]
## STATES 
states = [0,1,2,3,4,5,6,7]
## DEFINE INITIAL STATE PROBABILITIES 
imat = [0.14, 0.86]
## DEFINE EMISSION PROBABILITIES
emat = [[1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1],
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1]]
## READ FILE DEFINED BY SYS.ARGV[1]
filename = sys.argv[1]
#'C:\Users\Admin\Desktop\Dropbox\College\Junior\CSCI1820\homework\hw7\chr1_1000kb_1050kb_cpg.fasta' 

import os
import math
import random
import matplotlib.pyplot as plt
def viterbi(sequence, states, transitionMatrix, emissionProbabilities, initialProbabilities):
    #map base to emissions
    baseMap = dict([('A',0),('C',1),('G',2),('T',3)])
    L = len(sequence)
    ## initialize dynamic programming array
    vArray = [[0 for x in range(len(states))] for x in range(L+1)]
    ## initialize predicted state array
    ptr = [[0 for x in range(len(states))] for x in range(L+1)]
    ## calculate the initial probabilities given the first base
    for i in range(8):
        vArray[0][i]=-99999
    firstbase = baseMap[sequence[0]]
    vArray[0][firstbase] = math.log(initialProbabilities[0])
    vArray[0][firstbase+4] = math.log(initialProbabilities[1])
    ## recursive calls
    for i in range(1, L+1):
        for state in states:
            temp =[]
            for k in states:
                temp.append((math.log(transitionMatrix[state][k])+vArray[i-1][k]))
            #avoid log(0) - not that you really need to use emission probabilities anyway
            if emissionProbabilities[state][baseMap[sequence[i-1]]] ==0:
                vArray[i][state] = -999999
            else:
                vArray[i][state] =   max(temp) + math.log(emissionProbabilities[state][baseMap[sequence[i-1]]])
            ptr[i][state] = temp.index(max(temp))
    #termination case for pointers 
    pistar = vArray[L].index(max(vArray[L]))
    #initialize state array
    stateArray = [0 for x in range(L)]
    stateArray[L-1] = pistar
    #populate state array by tracing back pointers
    for x in range(L-1):
        stateArray[L-2-x] = ptr[L-1-x][stateArray[L-1-x]]

    # convert 0-3 to island, 4-7 to ocean 
    sa2=[]
    for s in stateArray:
        if s > 3:
            sa2.append(0)
        else: sa2.append(1)

    # plot to a scatter - can be difficult with a long sequence. 
    plt.scatter(range(len(sequence)),sa2, marker=".")
    plt.title("CpG island prediction. Islands are at 1, oceans at 0")
    plt.xlabel("Genomic position")
    plt.ylabel("Annotation of base")
    plt.show()
    return stateArray, vArray,ptr

# READ FILE INTO SEQUENCE
with open(filename,'r') as fastafile:
    sequence=''
    line=fastafile.readline().strip()
    while line != '':
        if line[0] != '>':
            sequence +=line
        line=fastafile.readline().strip()

#function call 
def main():
    data = viterbi(sequence, states, tmat, emat, imat)

if  __name__ =='__main__':main()


#TEST CASES I USED 
# simple island test case
seq1='CGCGCGCG'
# simple ocean test case
seq2='ATATATATAT'
# simply switching between occean and island
seq3= 'ATATATATATATATATATATCGCGCGCGCGCGCGCGATTCGCGCGCGCGCGCGCGATATATATATATATATATATATATATATATCGCGCGCGCGCGCGCGCGCG'
# first 200 bases of sequece, should be all island
seq4='GGTGGAGCGCGCCGCCACGGACCACGGGCGGGCTGGCGGGCGAGCGGCGAGCGCGCGGCGATCCGAGCCCCTAGGGCGGATCCCGGCTCCAGGCCCGCGCGCGCCTCAGGCCGTTTCCCTATTTAAGGCCTCGCCGCCGCGGGGTGTGTGAACCCGGCTCCGCATTCTTTCCCACACTCGCCCCAGCCAATCGACGGCCGCGCTCCTCCC'
