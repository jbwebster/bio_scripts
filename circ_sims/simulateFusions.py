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
	try:
		fid = fields[10]
		gene5 = re.findall("^ENST[0-9]*", fields[11])[0]
		gene3 = re.findall("::(ENST[0-9]*)", fields[11])[0]
		lastExon5Gene = fields[16]
		firstExon3Gene = fields[25]
		fusion = Fusion(fid, gene5, gene3, lastExon5Gene, firstExon3Gene)
		return fusion
	except:
		return None
	

def collectTopFusions(filename):
	# Count fusions
	fusion_ids = {}
	fusion_info = {}
	with open(filename, "r") as f:
		for line in f:
			fields = line.split("\t")
			if (fields[0]!="SAMPLE_ID"):
				fusion_id = fields[10]
				if(fusion_id in fusion_ids.keys()):
					fusion_ids[fusion_id] = fusion_ids[fusion_id] + 1
				else:
					fusion_ids[fusion_id] =  1
					fus = getFusionInfo(fields)
					fusion_info[fusion_id] = fus
	f.close()

	maxToKeep = 100
	numberKept = 0	
	keptGenes = []
	keptKeys = []
	keptFusionTmpIds = []
	# In later version of python, dicts keep their order
	sorted_dict = {k: v for k, v in sorted(fusion_ids.items(), key=lambda item: item[1])}
	sorted_keys = list(sorted_dict.keys())
	sorted_keys.reverse()
	for key in sorted_keys:
		if (key in fusion_info.keys()):
			curr = fusion_info[key]
			if(curr is not None):
				tmpId = curr.gene5 + "_" + curr.gene3
				if (tmpId not in keptFusionTmpIds and numberKept < maxToKeep and curr.gene5 not in keptGenes and curr.gene3 not in keptGenes):
					keptFusionTmpIds.append(tmpId)
					keptKeys.append(key)
					numberKept += 1
					keptGenes.append(curr.gene5)
					keptGenes.append(curr.gene3)
	
	topFusions = {}
	for key in fusion_info.keys():
		if key in keptKeys and key not in topFusions.keys():
			topFusions[key] = fusion_info[key]
	return topFusions

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
	


# Call main() and exit after finishing
if __name__ == '__main__':
	sys.exit(main(sys.argv))
