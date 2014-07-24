#parameters 
import argparse

parser = argparse.ArgumentParser(description="Script to find all repeated substrings of at least length k. Default parameters unless otherwise specified: k=3")
parser.add_argument("text", action="store", metavar="text", help="String of letters to search for repeated substrings")
parser.add_argument("--k", action="store", metavar="k", default="3", help="Find strings of at least this length")
args=parser.parse_args()
#parameters and paths specified in section above
text=args.text
k=int(args.k)

# repeated_substrings(text, k): output all repeated substrings in text of at least length k. 
def repeated_substrings(text, k):
	# generate substrings length >= k 
	# add to dictionary if unique. increase count if already present
	kmer = {}
	for i in range(k, len(text)):
		for j in range(len(text)-i+1):
			sub = text[j:j+i]
			if sub in kmer:
				kmer[sub] = kmer[sub] + 1
			else:
				kmer[sub] = 1
	# return all substrings occuring 2 or more times
	keys = [key for key,value in kmer.items() if value>=2]
	for key in keys:
		print key

# function call
def main():
	repeated_substrings(text,k)

if  __name__ =='__main__':main()