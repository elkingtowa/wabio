#parameters 
import argparse

parser = argparse.ArgumentParser(description="Script to find all repeated substrings of at least length k with at most d mismatches. Default parameters unless otherwise specified: k=3, d=1, alphabet=['A','T','C','G']")
parser.add_argument("text", action="store", metavar="text", help="String of letters to search for repeated substrings")
parser.add_argument("--k", action="store", metavar="k", default="3", help="Find strings of at least this length")
parser.add_argument("--d", action="store", metavar="d", default="1", help="Find strings with at most this many mismatches")
parser.add_argument("--alphabet", action="store", metavar="alphabet", default="ATCG", help="The alphabet to iterate over for mismatches")
args=parser.parse_args()
#parameters and paths specified in section above
text=args.text
k=int(args.k)
d=int(args.d) 
alphabet=args.alphabet

# repeated_substrings_mismatches(text, k, d, alphabet): output all repeated substrings of at least length k with at most d mismatches over alphabet

# code outline: find each kmer, generate possible mismatches, add occurances to dictionary
def repeated_substrings_mismatches(text, k, d, alphabet):
	#find each kmer 
	kmer = {}
	for i in range(k, len(text)):
		for j in range(len(text)-i+1):
			sub = text[j:j+i]
			#generate list of mismatches for a given kmer
			mis_list = generate_mismatches(sub, d, alphabet)
			for mismatch in mis_list:
				if mismatch in kmer:
					kmer[mismatch] = kmer[mismatch] + 1
				else:
					kmer[mismatch] = 1
	# return all substrings occuring 2 or more times
	keys = [key for key,value in kmer.items() if value>=2]
	for key in keys:
		print key

# generate_mismatches(kmer, d): returns a list of strings at most d mismatches from kmer, using a specified alphabet
# this is a helper function for repeated_substrings_mismatches and should not be called on its own. 
def generate_mismatches(kmer, d, alphabet):
	mis_list =[]
	if d == 1:
		for letter in range(len(kmer)):
			for sub in alphabet:
				misk = kmer[0:letter] + sub + kmer[letter+1:len(kmer)]
				if misk not in mis_list:
					mis_list.append(misk)
	else:
		print "well, I haven't generalized this to greater numbers of mismatches yet!"
	return mis_list

# function call
def main():
	repeated_substrings_mismatches(text, k, d, alphabet)

if  __name__ =='__main__':main()