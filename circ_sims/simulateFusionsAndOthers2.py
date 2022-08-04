#!/usr/bin/python

import sys
import argparse
import gzip
import re
import random

class Fusion(object):
	fusion_id = ""
	gene5 = ""
	gene3 = ""
	lastExon5Gene = ""
	firstExon3Gene = ""

	def __init__(self, fid, gene5, gene3, lastExon5Gene, firstExon3Gene):
		self.fusion_id = fid
		self.gene5 = gene5
		self.gene3 = gene3
		self.lastExon5Gene = lastExon5Gene
		self.firstExon3Gene = firstExon3Gene

	def __str__(self):
		return "\t".join([self.fusion_id, self.gene5, self.gene3, self.lastExon5Gene, self.firstExon3Gene])
	
 
# Create parser to collect inputs
def createParser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-fasta', action="store", help="Input fasta file")
	parser.add_argument('-fusions', action="store", help="Cosmic tsv of fusions")
	parser.add_argument('-output', action="store", help="Output file")
	parser.add_argument('-log', action="store", help="logfile")
	return parser

def getFusionInfo(fields):
	fid = re.findall("[0-9]+", fields[1])[0]
	gene5 = fields[2]
	gene3 = fields[5]
	lastExon5Gene = fields[4]
	firstExon3Gene = fields[7]
	fusion = Fusion(fid, gene5, gene3, lastExon5Gene, firstExon3Gene)
	writeToLog("\t".join([fid, gene5, lastExon5Gene, gene3, firstExon3Gene]) + "\n")
	return fusion
	

def collectTopFusions(filename):
	fusions = {}
	counter = 0
	with open(filename, "r") as f:
		for line in f:
			counter += 1
			if (counter > 2):
				fields = line.split("\t")
				fus = getFusionInfo(fields)
				fusionid = re.findall("[0-9]+", fields[1])[0]
				fusions[fusionid] = fus
	f.close()
	return fusions	

def extractTranscriptIds(fusionDict):
	transcripts = []
	for key in fusionDict.keys():
		transcripts.append(fusionDict[key].gene5)
		transcripts.append(fusionDict[key].gene3)
	return transcripts

def keepEntry(transcriptId, fields, topFusions):
	exon_number = re.findall('exon_number "([0-9]*)";', fields[8])[0]
	for key in topFusions.keys():
		fusion = topFusions[key]
		if (fusion.gene5 == transcriptId):
			if(exon_number <= fusion.lastExon5Gene):
				return True
		elif (fusion.gene3 == transcriptId):
			if(exon_number >= fusion.firstExon3Gene):
				return True
	return False

def keepEntry2():
	return True
	
			
def readInFasta(fasta):
	fasta_seqs = {}
	header = ""
	fullseq = ""
	exon = ""
	transcript = ""
	key = ""	
	with open(fasta, "r") as f:
		for line in f:
			if line[0] == ">":
				if header != "":
					fasta_seqs[key] = fullseq
					header = ""
					transcript = ""
					exon = ""
					key = ""
					fullseq = ""
				header = line.strip();
				transcript = re.findall('transcript_id_"(ENST[0-9]*)', header)[0]
				exon = re.findall('exon_number_"([0-9]*)', header)[0]
				key = transcript + "_" + exon
			else:
				fullseq += line
	f.close()
	return fasta_seqs

def simulateReportedFusions(seqs, fusions):
	sims = {}
	for key in fusions:
		writeToLog("Creating transcript for: " + str(fusions[key]))	
		fusion = fusions[key]
		gene5 = fusion.gene5
		lastExon = int(fusion.lastExon5Gene)
		currExon = 1
		fullseq = ""
		while (currExon <= lastExon):
			seq_key = gene5 + "_" + str(currExon)
			writeToLog("\tCurrKey: " + seq_key)
			fullseq += seqs[seq_key]
			currExon += 1
		gene3 = fusion.gene3
		currExon = int(fusion.firstExon3Gene)
		exonExists = True
		while (exonExists):
			seq_key = gene3 + "_" + str(currExon)
			writeToLog("\tCurrKey: " + seq_key)
			if (seq_key in seqs.keys()):
				fullseq += seqs[seq_key]
			else:
				exonExists = False
			currExon += 1
		header = ">" + gene5 + " to exon " + str(lastExon) + " >> " + gene3 + " at exon " + fusion.firstExon3Gene + "\n"
		sims[header] = fullseq
	return sims
	

def simulateFcircJunctions(seqs, fusions):		
	random.seed(42)
	sims = {}
	hasfcircRNA = False
	count = 0
	for key in fusions:
		v = random.randrange(100)
		fusion = fusions[key]
		if (v > 10):
			hasfcircRNA = True
		else:
			hasfcircRNA = False
		if(hasfcircRNA):
			s = "Attempting to create fcircRNA for fusion " + fusion.gene5 + ">>" + fusion.gene3 + ":"
			header = ">fcircRNA from " + fusion.gene5 + " to " + fusion.gene3 + " " + fusion.gene3 + ":"
			backspliceStart = int(fusion.firstExon3Gene) + 3 + int(random.randrange(5))
			backspliceEnd = int(fusion.lastExon5Gene) - 3 - int(random.randrange(5))
			threepgeneseq = ""
			fivepgeneseq = ""
			hasValidResult = True
			if (backspliceEnd < 1):
				s += "Failed. Invalid backspliceEnd in 5' gene"	
			else:
				currExon = int(fusion.firstExon3Gene)
				isValid = True
				
				while(isValid):
					seq_key = fusion.gene3 + "_" + str(currExon)
					if(seq_key in seqs.keys() and currExon <= backspliceStart):
						threepgeneseq += seqs[seq_key]
						header += str(currExon) + ","
						currExon += 1
					else:
						isValid = False
					if(currExon >= backspliceStart):
						s += "Included final exon of 3' gene. "
				if(threepgeneseq == ""):
					hasValidResult = False	
				isValid = True
				currExon = backspliceEnd
				header += " " + fusion.gene5 + ":"
				while(isValid):
					seq_key = fusion.gene5 + "_" + str(currExon)
					if(seq_key in seqs.keys() and currExon <= int(fusion.lastExon5Gene)):
						fivepgeneseq += seqs[seq_key]
						header += str(currExon) + ","
						currExon += 1
					else:
						isValid = False
					if(currExon == 1):
						s += "Included first exon of 5' gene. "
				if(fivepgeneseq == ""):
					hasValidResult = False
				fullseq = threepgeneseq + fivepgeneseq + threepgeneseq
				header += "\n"
				if(hasValidResult):
					sims[header] = fullseq
					count += 1
					s += "Success"
				else:
					s += " Invalid exons. Failed"
			writeToLog(s)
	print("Created fcircRNA for " + str(count) + " out of " + str(len(fusions.keys())) + "fusions")
	return sims	
					
		


def writeOutput(filename, sims):
	with open(filename, "w") as f:
		for key in sims.keys():
			f.write(key)
			f.write(sims[key])
	f.close()

def writeToLog(s):
	if (append):
		with open(logfile, "a") as f:
			f.write(s + "\n")
		f.close()
	else:
		with open(logfile, "w") as f:
			f.write(s + "\n")
		f.close()
	


# Main driver
def main(args):
	parser = createParser()
	args = parser.parse_args()
	global logfile
	logfile = args.log
	global append
	append = False
	writeToLog("Creating published fusion transcripts")
	append = True
	print("Collecting top fusions")
	topFusions = collectTopFusions(args.fusions)
	print("Reading in fasta")
	seqs = readInFasta(args.fasta) # Returns dict
	print("Simulating fusions")
	sims = simulateReportedFusions(seqs, topFusions)
	writeOutput(args.output, sims)
	print("Simulating fcircRNA")
	fcircSims = simulateFcircJunctions(seqs, topFusions)
	writeOutput(args.output + ".fcirc.fa", fcircSims)
	


# Call main() and exit after finishing
if __name__ == '__main__':
	sys.exit(main(sys.argv))
