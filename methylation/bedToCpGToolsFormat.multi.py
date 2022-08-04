#!/usr/bin/python

import sys
import argparse
import gzip
import re
import itertools
 
# Create parser to collect inputs
def createParser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-files', action="store", help="Comma separated list of cpgs.bed.gz files")
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
	formatted_string = str(int(methylated_reads)) + "," + str(int(total_reads))
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
def readFiles(infiles):
	counter = 0
	files = []
	for filename in infiles:
		if(".gz" in filename):
			files.append(gzip.open(filename))
		else:
			files.append(open(filename))
	#files = [open(f) for f in infiles]
	print("cgID\t" + "\t".join(infiles))
	for lines in itertools.izip(*files):
		#Confirm that lines from files are all for same position in genome
		currPos = ""
		useLine = True
		for line in lines:
			fields = line.split("\t")
			cpgid = createID(fields)	
			if currPos=="":
				currPos = cpgid
			elif currPos!=cpgid:
				useLine = False
		# If lines descibe same location
		if useLine:
			outstring = currPos
			for line in lines:			
				samplestring = convertFormat(line.split("\t"))
				outstring += "\t" + samplestring	
			print(outstring)
		counter += 1
		if(counter>5):
			sys.exit()

	
	#print("cgID\tsample_1\tsample_2")
	# Not sure how the next 3 lines work. Pulled from StackOverflow...
	#filenames = [filea, fileb]
	#gens = [gen_line(n) for n in filenames]
	#for linea, lineb in zip(*gens):
	#	fieldsa = linea.split("\t")
	#	fieldsb = lineb.split("\t")
	#	cpgid = createID(fieldsa)
	#	stringa = convertFormat(fieldsa)
	#	stringb = convertFormat(fieldsb)	
	#	outstring = "\t".join([cpgid, stringa, stringb])
	#	print(outstring)


# Main driver
def main(args):
	parser = createParser()
	args = parser.parse_args()
	infiles = args.files.split(",")	
	readFiles(infiles)


# Call main() and exit after finishing
if __name__ == '__main__':
	sys.exit(main(sys.argv))
