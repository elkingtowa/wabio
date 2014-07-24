import csv
import sys
import getopt
import json

#Given a .tsv file outputted by R, this script processes it into a textfile
#readable by Phylip
#Syntax is RtoPhylip.py -i inputFile

#the number of digits to export in matrix values
DIGITS = 16
#Specify disallowed characters and what to replace them with here
DISALLOWEDCHARS = {'(':'{', ')':'}', ':':'-', ';':'&', '[':'{', ']':'}'}

def openTSV(filename):
    #import the tsv
    with open(filename, 'r') as data:
        dataReader = csv.reader(data, dialect='excel-tab')
        data.close
        #convert into a two dimensional list
        return list(dataReader)

#replaces all chars in text that match a dic key with the dic value
def replace_all(text, dic):
    for i, j in dic.items():
        text = text.replace(i, j)
    return text


def process(matrix):
    #do ridiculous text preprocessing

    for row in matrix[1:]:
        #clip and pad names to 10 chars
        name = replace_all(row[0],DISALLOWEDCHARS)        
        row[0] = name[0:10].ljust(10)
        for i in range(1,len(row)):
            row[i] = row[i][0:DIGITS].zfill(DIGITS)
    #output a processed file
    stringrow = []
    for row in matrix[1:]:
        stringrow.append('  '.join(row))
    processed = '\n'.join(stringrow)
    return str(len(matrix) - 1) + '\n' + processed

def main(argv):
   inputfile = 'data/import.txt'
   outputfile = 'data/output.txt'
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["help", "ifile=","ofile="])
   except getopt.GetoptError:
      print('RtoPhylip.py -i <inputfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('RtoPhylip.py -i inputPath -o outputPath')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
   processed = process(openTSV(inputfile))
   output = open(outputfile, "w")
   output.write(processed)
   output.close

if __name__ == "__main__":
   main(sys.argv[1:])