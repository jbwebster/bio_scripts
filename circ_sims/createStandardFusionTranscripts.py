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
		self.name = None
		self.sequence = None

	def __str__(self):
		return "\t".join([self.fusion_id, self.gene5, self.gene3, self.lastExon5Gene, self.firstExon3Gene])
	
 
# Create parser to collect inputs
def createParser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-gtf', action="store", help="gtf file")
	parser.add_argument('-fusions', action="store", help="Cosmic tsv of fusions")
	parser.add_argument('-output', action="store", help="Output file")
	parser.add_argument('-log', action="store", help="Log file")
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
					if(fus is not None):
						fusion_info[fusion_id] = getFusionInfo(fields)
	f.close()
	writeToLog("Number of fusion ids: " + str(len(fusion_ids.keys())) + "\n")
	writeToLog("Number of fusion infos: " + str(len(fusion_info.keys())) + "\n")

	maxToKeep = 100
	numberKept = 0	
	keptKeys = []
	keptFusionTmpIds = []
	# In later version of python, dicts keep their order
	sorted_dict = {k: v for k, v in sorted(fusion_ids.items(), key=lambda item: item[1])}
	sorted_keys = list(sorted_dict.keys())
	sorted_keys.reverse()
	for key in sorted_keys:
		if (key in fusion_info.keys()):
			curr = fusion_info[key]
			tmpId = curr.gene5 + "_" + curr.gene3
			if (tmpId not in keptFusionTmpIds and numberKept < maxToKeep):
				keptFusionTmpIds.append(tmpId)
				keptKeys.append(key)
				numberKept += 1
	
	topFusions = {}
	for key in fusion_info.keys():
		if key in keptKeys and key not in topFusions.keys():
			topFusions[key] = fusion_info[key]
	writeToLog("Fusions kept:\n")
	for key in topFusions.keys():
		f = topFusions[key]
		os = str(f)
		writeToLog(os + "\n")
	return topFusions

def extractTranscriptIds(fusionDict):
	transcripts = []
	for key in fusionDict.keys():
		transcripts.append(fusionDict[key].gene5)
		transcripts.append(fusionDict[key].gene3)
	writeToLog("Keeping entries from " + str(len(transcripts)) + " transcripts\n")
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
	
			
def parseGTF(filename, topFusions):

	transcripts = extractTranscriptIds(topFusions)
	writeToLog("Transcripts (" + str(len(transcripts)) + "):\n")
	for tra in transcripts:
		writeToLog(tra + "\n")
	keptLines = []	
	with open(filename, "r") as f:
		for line in f:
			if line[0] == "#":
				keptLines.append(line)
			else:
				fields = line.split("\t")
				if(fields[2] == "exon"):
					transcriptId = re.findall("ENST[0-9]*", fields[8])[0]
					if transcriptId in transcripts:
						if(keepEntry(transcriptId, fields, topFusions)):
							fields[0] = "chr" + fields[0]
							outline = "\t".join(fields)
							keptLines.append(outline)
	f.close()
	return keptLines

def writeOutput(filename, entries):
	with open(filename, "w") as f:
		for entry in entries:
			f.write(entry)
	f.close()

def writeToLog(string):
	with open(logfile, "a") as f:
		f.write(string)
	f.close()	

# Main driver
def main(args):
	parser = createParser()
	args = parser.parse_args()
	global logfile;
	logfile = args.log
	cmd_str = "python3 extractFusionGeneExons.py -gtf " + args.gtf + " -fusions " + args.fusions + " -output " + args.output + " -log " + args.log + "\n"
	writeToLog(cmd_str)
	topFusions = collectTopFusions(args.fusions)
	writeToLog("Num of fusions kept: " + str(len(topFusions.keys())) + "\n")

	# Collect information from gtf
	fusion_entries = parseGTF(args.gtf, topFusions)

	writeOutput(args.output, fusion_entries)
	


# Call main() and exit after finishing
if __name__ == '__main__':
	sys.exit(main(sys.argv))
