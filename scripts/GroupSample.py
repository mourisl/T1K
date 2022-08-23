#!/usr/bin/env python3

#import sys
import argparse


if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(description = "Group samples into people-level")
	parser.add_argument("-l", help="list of genotyping results", dest="filelist", required=True)
	parser.add_argument("-q", help="ignore allels with less or equal quality scores", dest="qual", required=False, default=29)
	parser.add_argument("-d", help="number of HLA digits", dest="digits",
			required=False, default=2)

	args = parser.parse_args()

	fpList = open(args.filelist)
	geneSet = set({"HLA-A", "HLA-B", "HLA-C"})
	sampleSignature = {}
	badSamples = set({})
	signatureToSamples = {}
	sampleGroupId = {}
	for f in fpList:
		f = f.rstrip()
		fp = open(f)
		for line in fp:
			cols = line.rstrip().split("\t")
			if (cols[0] not in geneSet):
				continue
			if (f not in sampleSignature):
				sampleSignature[f] = set({})
			if (int(cols[1]) >= 1):
				qual = int(cols[4]) 
				sampleSignature[f].add(":".join(cols[2].split(',')[0].split(':')[0:(args.digits)]))
				if (qual <= args.qual):
					badSamples.add(f)
			if (int(cols[1]) >= 2):
				qual = int(cols[7])
				sampleSignature[f].add(":".join(cols[5].split(',')[0].split(':')[0:(args.digits)]))
				if (qual <= args.qual):
					badSamples.add(f)
	
	for s in sampleSignature:
		if (s in badSamples):
			sampleGroupId[s] = -1
			continue
		signature = tuple(sorted(list(sampleSignature[s])))
		if (signature not in signatureToSamples):
			signatureToSamples[signature] = []
		signatureToSamples[signature].append(s)
	
	i = 0
	for signature, samples in signatureToSamples.items():
		for s in samples:
			sampleGroupId[s] = i
		i += 1

	for s, groupid in sampleGroupId.items():
		print(s, groupid)
