#!/usr/bin/python

import sys
import argparse
import gzip
import re
 
# Create parser to collect inputs
def createParser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-file1', action="store", help="Sample 1 cpgs.bed.gz file")
	parser.add_argument('-file2', action="store", help="Sample 2 cpgs.bed.gz file")
	return parser

def createID(fields):
	newid = "_".join(fields[0:3])
	return newid	

def convertFormat(fields):
	total_reads = int(fields[4])
	methylated_reads = 0
	try:
		Creads = int(re.findall('([0-9]+),G', fields[5])[0])
		Cbeta = float(re.findall('C:([0-9\\.]+):', fields[5])[0])
		methylated_reads = total_reads - (Creads - round(Cbeta * Creads))
	except:
		methylated_reads = 0
	formatted_string = str(methylated_reads) + "," + str(total_reads)
	return formatted_string
	

# No idea how this works..
def gen_line(filename):
	if (".gz" in filename):
		with gzip.open(filename) as f:
			for line in f:
				yield str(line.strip().decode())
	else:
		with open(filename) as f:
			for line in f:
				yield str(line.strip())

# Read the provided files
def readFiles(filea, fileb):
	print("cgID\tsample_1\tsample_2")
	# Not sure how the next 3 lines work. Pulled from StackOverflow...
	filenames = [filea, fileb]
	gens = [gen_line(n) for n in filenames]
	for linea, lineb in zip(*gens):
		fieldsa = linea.split("\t")
		fieldsb = lineb.split("\t")
		cpgid = createID(fieldsa)
		stringa = convertFormat(fieldsa)
		stringb = convertFormat(fieldsb)	
		outstring = "\t".join([cpgid, stringa, stringb])
		print(outstring)


# Main driver
def main(args):
	parser = createParser()
	args = parser.parse_args()
	readFiles(args.file1, args.file2)


# Call main() and exit after finishing
if __name__ == '__main__':
	sys.exit(main(sys.argv))
