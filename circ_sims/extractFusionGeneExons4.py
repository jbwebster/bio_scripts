#!/usr/bin/python

import sys
import argparse
import gzip
import re

class Fusion(object):
	fusion_id = ""
	gene5 = ""
	gene3 = ""
	lastExon5Gene = ""
	firstExon3Gene = ""
	gene5name = ""
	gene3name = ""
	count = 0

	def __init__(self, fid, gene5, gene3, lastExon5Gene, firstExon3Gene, gene5name, gene3name):
		self.fusion_id = fid
		self.gene5 = gene5
		self.gene3 = gene3
		self.lastExon5Gene = lastExon5Gene
		self.firstExon3Gene = firstExon3Gene
		self.gene5name = gene5name
		self.gene3name = gene3name

	def __str__(self):
		return "\t".join([self.fusion_id, self.gene5, self.gene3, self.lastExon5Gene, self.firstExon3Gene, self.gene5name, self.gene3name, str(self.count)])
	
	def increment(self):
		self.count += 1
 
# Create parser to collect inputs
def createParser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-gtf', action="store", help="Original gtf file")
	parser.add_argument('-fusions', action="store", help="Cosmic tsv of fusions")
	parser.add_argument('-output', action="store", help="Output file")
	parser.add_argument('-log', action="store", help="Log file")
	return parser

def getFusionInfo(fields):
	fid = re.findall("[0-9]+", fields[1])[0]
	gene5 = fields[2]
	gene3 = fields[5]
	lastExon5Gene = fields[4]
	firstExon3Gene = fields[7]
	gene5name = fields[3]
	gene3name = fields[6]
	fusion = Fusion(fid, gene5, gene3, lastExon5Gene, firstExon3Gene, gene5name, gene3name)
	writeToLog("\t".join([fid, gene5name, gene5, lastExon5Gene, gene3name, gene3, firstExon3Gene]) + "\n", True)
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
	writeToLog("Keeping entries from " + str(len(transcripts)) + " transcripts\n", True)
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
	
			
def parseGTF(filename, topFusions):

	transcripts = extractTranscriptIds(topFusions)
	writeToLog("Transcripts (" + str(len(transcripts)) + "):\n", True)
	for tra in transcripts:
		writeToLog(tra + "\n", True)
	keptLines = []	
	prevGene = ""
	geneAdded = False
	with open(filename, "r") as f:
		for line in f:
			if line[0] == "#":
				keptLines.append(line)
			else:
				fields = line.split("\t")
				if(fields[2] == "exon"):
					transcriptId = re.findall("ENST[0-9]*", fields[8])[0]
					if transcriptId in transcripts:
						if (not geneAdded):
							keptLines.append(prevGene)
							geneAdded = True
						#if(keepEntry(transcriptId, fields, topFusions)):
						if(keepEntry2()):
							fields[0] = "chr" + fields[0]
							outline = "\t".join(fields)
							keptLines.append(outline)
				if(fields[2] == "gene"):
					fields[0] = "chr" + fields[0]
					prevGene = "\t".join(fields)
					geneAdded = False
	f.close()
	return keptLines

def writeOutput(filename, entries):
	with open(filename, "w") as f:
		for entry in entries:
			f.write(entry)
	f.close()

def writeToLog(string, append):
	if(append):
		with open(logfile, "a") as f:
			f.write(string)
		f.close()	
	if(append is False):
		with open(logfile, "w") as f:
			f.write(string)
		f.close()

# Main driver
def main(args):
	parser = createParser()
	args = parser.parse_args()
	global logfile;
	logfile = args.log
	cmd_str = "python3 extractFusionGeneExons.py -gtf " + args.gtf + " -fusions " + args.fusions + " -output " + args.output + " -log " + args.log + "\n"
	writeToLog(cmd_str, False)
	topFusions = collectTopFusions(args.fusions)
	writeToLog("Num of fusions kept: " + str(len(topFusions.keys())) + "\n", True)

	# Collect information from gtf
	fusion_entries = parseGTF(args.gtf, topFusions)

	writeOutput(args.output, fusion_entries)
	


# Call main() and exit after finishing
if __name__ == '__main__':
	sys.exit(main(sys.argv))
