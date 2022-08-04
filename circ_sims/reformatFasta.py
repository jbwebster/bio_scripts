#!/usr/bin/python

import sys
import argparse
import gzip
import re

 
# Create parser to collect inputs
def createParser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-fasta', action="store", help="fasta file")
	parser.add_argument('-output', action="store", help="Output file")
	parser.add_argument('-log', action="store", help="Log file")
	return parser

def parseFasta(infile):
	print("Parsing fasta")
	isHeader = False
	append = False
	currHeader = ""
	currSeq = ""
	with open(infile, "r") as f:
		print("Fasta opened")
		for line in f:
			if line[0] == ">" and isHeader == False:
				if (currHeader == ""):
					isHeader = True
					currHeader = line
				elif(currHeader != ""):
					writeToOutput(currHeader, append)
					append = True
					writeToOutput(currSeq.rstrip() + "\n", append)	
					isHeader = True
					currHeader = line
					currSeq = ""
				writeToLog("Parsing : " + currHeader, True)
			elif line[0] == ">" and isHeader == True:
				# Sequence for previous header was missing
				writeToLog("Missing sequence for: " + currHeader, True)
				currHeader = line
				currSeq = ""
			else:
				isHeader = False
				fullseq = line.rstrip()
				seqlength = len(fullseq)
				newNumLines = int(seqlength / 50) # Truncates float
				currLineCount = 1
				currSeq = ""
				while(currLineCount <= newNumLines):
					startValue = (currLineCount - 1) * 50  # line 1 -> 0, line 2 -> 50
					endValue = startValue + 50 # line 1 -> 50, line 2 -> 100
					currLineSeq = fullseq[startValue:endValue] + "\n"
					currSeq = currSeq + currLineSeq
					currLineCount += 1
				extra = seqlength % 50
				if(extra > 0):
					currSeq = currSeq + fullseq[(len(fullseq) - extra):]
		writeToOutput(currHeader, append)
		writeToOutput(currSeq.rstrip() + "\n", append)
	f.close()
	print("Done parsing")
				
				
				

def writeToOutput(string, append):
	if (append):
		with open(outfile, "a") as f:
			f.write(string)
		f.close()
	else:
		with open(outfile, "w") as f:
			f.write(string)
		f.close()

def writeToLog(string, append):
	if (append):
		with open(logfile, "a") as f:
			f.write(string)
		f.close()	
	else:
		with open(logfile, "w") as f:
			f.write(string)
		f.close()

# Main driver
def main(args):
	parser = createParser()
	args = parser.parse_args()
	global logfile
	logfile = args.log
	global outfile
	outfile = args.output
	cmd_str = "python3 reformatFasta.py -fasta " + args.fasta + " -output " + args.output + " -log " + args.log + "\n"
	writeToLog(cmd_str, False)
	
	parseFasta(args.fasta)

	


# Call main() and exit after finishing
if __name__ == '__main__':
	sys.exit(main(sys.argv))
