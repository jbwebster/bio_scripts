#!/usr/bin/python

import sys
import argparse
import gzip
import re

 
# Create parser to collect inputs
def createParser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-gtf', action="store", help="gtf file")
	parser.add_argument('-output', action="store", help="Output file")
	parser.add_argument('-log', action="store", help="Log file")
	return parser


def gtfToBed(gtf, bed):
	append = False
	with open(gtf, "r") as f:
		for line in f:
			if line[0] != "#":
				fields = line.split("\t")
				chrm = fields[0]
				pos = fields[3]
				endpos = fields[4]
				lastfield = fields[8].rstrip()
				lastfield = lastfield.replace(" ", "_")
				s = chrm + "\t" + pos + "\t" + endpos + "\t" + chrm + "_" + pos + "_" + endpos + ":" + lastfield + "\n"
				if (fields[2] != "gene"):
					writeOutput(bed, s, append)
					append = True
		
					
	f.close()	

def writeOutput(filename, string, append):
	if(append):
		with open(filename, "a") as f:
			f.write(string)
		f.close()
	else:
		with open(filename, "w") as f:
			f.write(string)
		f.close()


# Main driver
def main(args):
	parser = createParser()
	args = parser.parse_args()
	
	gtfToBed(args.gtf, args.output)



# Call main() and exit after finishing
if __name__ == '__main__':
	sys.exit(main(sys.argv))
