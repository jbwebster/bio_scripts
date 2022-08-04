#!/usr/bin/python

import sys
import argparse
import gzip
import re
 
# Create parser to collect inputs
def createParser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-file1', action="store", help="Sample 1 cpgs.bed.gz file")
	return parser

def createChrBase(fields):
	newid = fields[0] + "." + fields[1]
	return newid	

def getReadInfo(fields):
	totalReads = int(fields[4])
	cfreq = 0.0
	tfreq = 0.0
	try:
		Creads = int(re.findall('([0-9]+),G', fields[5])[0])
		cfreq = (1.0 * Creads / totalReads) * 100
		tfreq = 100 - cfreq
	except:
		cfreq = 0.0
		treq = 0.0
	return totalReads, cfreq, tfreq

def doStuff(fields):
	chrbase = createChrBase(fields)	
	chrm = fields[0]
	base = fields[1]
	strand = "F" # Probably?
	totalReads, cfreq, tfreq = getReadInfo(fields)
	outstring = "\t".join([chrbase,chrm,base,strand,str(totalReads),str(cfreq),str(tfreq)])
	return outstring

# Read the provided files
def readFile(filename):
	print("chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT")

	if(".gz" in filename):
		with gzip.open(filename, "r") as f:
			for line in f:
				fields = line.split("\t")
				print(doStuff(fields))
		f.close()
	else:
		with open(filename, "r") as f:
			for line in f:
				fields = line.split("\t")
				print(doStuff(fields))
		f.close()
				

# Main driver
def main(args):
	parser = createParser()
	args = parser.parse_args()
	readFile(args.file1)


# Call main() and exit after finishing
if __name__ == '__main__':
	sys.exit(main(sys.argv))
