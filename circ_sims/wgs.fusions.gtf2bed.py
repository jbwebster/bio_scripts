#!/usr/bin/python

import sys
import argparse
import gzip
import re
import random

class Fusion(object):
	gene = ""
	exonBreak = 0
	orientation = 0
	strand = ""
	start = 0
	end = 0
	preBreak = 0
	postBreak = 0
	breakpoint = 0
	chrm = ""

	def __init__(self, gene, exonBreak, orientation):
		self.gene = gene
		self.exonBreak = exonBreak
		self.orientation = orientation
		self.strand = None
		self.start = None
		self.end = None
		self.preBreak = None
		self.postBreak = None
		self.breakpoint = None
		self.chrm = None

	def __str__(self):
		return "_".join([str(self.chrm), str(self.gene), str(self.orientation), str(self.strand), str(self.preBreak), str(self.postBreak), str(self.breakpoint)])
	
 
# Create parser to collect inputs
def createParser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-gtf', action="store", help="Original gtf file")
	parser.add_argument('-fusions', action="store", help="fusions")
	parser.add_argument('-output', action="store", help="Genomic output file")
	parser.add_argument('-log', action="store", help="Log file")
	return parser

def parseFusions(infile):
	fusions = {}
	counter = 0
	transcripts = []
	with open(infile, "r") as f:
		for line in f:
			counter += 1
			if(counter > 2):
				fields = line.split("\t")
				geneA = fields[3]
				fus = Fusion(geneA, int(fields[4]), 5)
				fusions[geneA] = fus
				geneB = fields[6]
				fusB = Fusion(geneB, int(fields[7]), 3)
				fusions[geneB] = fusB
				transcripts.append(fields[2])
				transcripts.append(fields[5])
	return fusions, transcripts

def getExonEnds(gtf, output, fusions, transcripts):
	geneStart = 0
	geneEnd = 0
	geneStrand = " "
	atBreak = False
	with open(gtf, "r") as f:
		for line in f:
			if line[0] != "#":
				fields = line.split("\t")
				gene = re.findall('gene_name "(.*)"; gene_source', fields[8])[0]
				if(fields[2]=="gene" and gene in fusions.keys()):
					fus = fusions[gene]
					geneStart = int(fields[3])
					geneEnd = int(fields[4])
					geneStrand = fields[6]
					fus.start = geneStart - 500
					fus.end = geneEnd + 500
					fus.strand = geneStrand
					fus.chrm = "chr" + fields[0]
				if(fields[2]=="exon" and gene in fusions.keys()):
					transcript = re.findall('transcript_id "(.*)"; transcript_version', fields[8])[0]
					if transcript in transcripts:
						currExon = int(re.findall('exon_number "(.*)"; gene_name', fields[8])[0])
						fus = fusions[gene]

						if(fus.orientation==5):
							if fus.strand == "+":
								if currExon==fus.exonBreak:
									fus.preBreak = int(fields[4])
								elif currExon==fus.exonBreak + 1:
									fus.postBreak = int(fields[3])
							if fus.strand == "-":
								if currExon==fus.exonBreak:
									fus.postBreak = int(fields[3])
								elif currExon==fus.exonBreak + 1:
									fus.preBreak = int(fields[4])
						if(fus.orientation==3):
							if fus.strand == "+":
								if currExon==fus.exonBreak:
									fus.postBreak = int(fields[3])
								elif currExon==fus.exonBreak - 1:
									fus.preBreak = int(fields[4])
							if fus.strand == "-":
								if currExon==fus.exonBreak:
									fus.preBreak = int(fields[4])
								if currExon==fus.exonBreak - 1:
									fus.postBreak = int(fields[3])
						fusions[gene] = fus
						writeLog("Curr exon: " + str(currExon) + " Searching for: " + str(fus.exonBreak))
						writeLog(str(fus))
	f.close()	
	return fusions

def createBreakpoints(fusions):
	for key in fusions.keys():
		fus = fusions[key]
		dist = fus.postBreak - fus.preBreak
		bp = 0
		if fus.orientation == 5 and fus.strand == "+":
			bp = random.randint(fus.preBreak + int(0.1 * dist), int(fus.preBreak + (0.4 * dist)))
		elif fus.orientation == 5 and fus.strand == "-":
			bp = random.randint(int(fus.preBreak + (0.6 * dist)), fus.postBreak - int(0.1 * dist))
		elif fus.orientation == 3 and fus.strand == "+":
			bp = random.randint(int(fus.preBreak + (0.6 * dist)), fus.postBreak - int(0.1 * dist))
		elif fus.orientation == 3 and fus.strand == "-":
			bp = random.randint(fus.preBreak + int(0.1 * dist), int(fus.preBreak + (0.4 * dist)))
		fus.breakpoint = bp
	return fusions
		
		

def writeOutput(filename, fusions):
	with open(filename, "w") as f:
		for key in fusions.keys():
			fus = fusions[key]
			s = ""
			start = str(fus.start)
			breakpoint = str(fus.breakpoint)
			end = str(fus.end)
			if fus.orientation == 5 and fus.strand == "+":
				s = fus.chrm + "\t" + start + "\t" + breakpoint + "\t" + str(fus) + "\n"
			elif fus.orientation == 5 and fus.strand == "-":
				s = fus.chrm + "\t" + breakpoint + "\t" + end + "\t" + str(fus) + "\n"
			elif fus.orientation == 3 and fus.strand == "+":
				s = fus.chrm + "\t" + breakpoint + "\t" + end + "\t" + str(fus) + "\n"
			elif fus.orientation == 3 and fus.strand == "-":
				s = fus.chrm + "\t" + start + "\t" + breakpoint + "\t" + str(fus) + "\n"
			f.write(s)
	f.close()

def printFusions(fusions):
	s = ""
	for key in fusions.keys():
		s += str(fusions[key]) + "\n"
	return s

def validateFusions(fusions):
	s = ""
	for key in fusions.keys():
		fus = fusions[key]
		v = "FALSE"
		if fus.preBreak < fus.postBreak:
			v = "TRUE"
		s += str(fus) + ":" + v + ": " + str(fus.postBreak - fus.preBreak) + "\n"
	return s

def writeLog(s):
	if(append):
		with open(logfile, "a") as f:
			f.write(s)
		f.close()
	else:
		with open(logfile, "w") as f:
			f.write(s)
		f.close()

# Main driver
def main(args):
	random.seed(42)
	parser = createParser()
	args = parser.parse_args()
	global logfile
	logfile = args.log
	global append
	append = False
	
	print("Parsing selected fusions")
	fusions, transcripts = parseFusions(args.fusions)	
	writeLog(",".join(transcripts))
	print("Finding exon ends")
	fusions = getExonEnds(args.gtf, args.output, fusions, transcripts)
	writeLog(printFusions(fusions))
	print("Creating breakpoints")
	writeLog(validateFusions(fusions))
	fusions = createBreakpoints(fusions)
	print("Writing output")
	writeOutput(args.output, fusions)



# Call main() and exit after finishing
if __name__ == '__main__':
	sys.exit(main(sys.argv))
