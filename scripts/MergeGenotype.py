#!/usr/bin/env python3

# Merge the genotype results from multiple results, e.g: smart-seq

import sys
import argparse
import re

geneAlleles = {}
geneAlleleList = {}
files = []
finalAlleles = {}

if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(description = "Combine the genotyping results from mutiple files.")
	parser.add_argument("-l", help="list of genotyping results", dest="filelist", required=True)
	parser.add_argument("-n", help="number of alleles per gene", dest="numAllelePerGene", required=False, default=2)
	parser.add_argument("-q", help="ignore allels with less or equal quality scores", dest="qual", required=False, default=0)
	parser.add_argument("--tq", help="ignore allels with less or equal total quality scores", dest="totalQual", required=False, default=30)

	args = parser.parse_args()
	
	# Collect the allele information
	fpList = open(args.filelist)
	for f in fpList:
		f = f.rstrip()
		files.append(f)
		fp = open(f)
		for line in fp:
			cols = line.rstrip().split("\t")
			gene = cols[0]
			if (gene not in geneAlleles):
				geneAlleles[gene] = {}
			for i in [2, 5]:
				if (i < len(cols) and float(cols[i + 2]) > args.qual):
					equalAlleles = cols[i].split(",")

					for allele in equalAlleles[0:1]: # only use the first one for voting
						if (allele not in geneAlleles[gene]):
							geneAlleles[gene][allele] = 0
						geneAlleles[gene][allele] += float(cols[i + 2])
		fp.close()
	fpList.close()

	# Select the representative alleles 
	for gene in geneAlleles:
		for allele in sorted(geneAlleles[gene].keys(), key=lambda x:geneAlleles[gene][x], reverse=True)[0:int(args.numAllelePerGene)]:
			if (geneAlleles[gene][allele] >= float(args.totalQual)):
				finalAlleles[allele] = geneAlleles[gene][allele]
	
	# output the count matrix
	output = ["sample"]
	output += sorted(finalAlleles.keys())
	output += ["inconsistency"]
	print("\t".join(output))
	for f in files:
		fp = open(f)
		output = []
		sampleAlleles = {}
		inconsistentAlleles = []
		for allele in finalAlleles:
			sampleAlleles[allele] = 0
		for line in fp:
			cols = line.rstrip().split("\t")
			gene = cols[0]
			
			for i in [2, 5]:
				if (i < len(cols) and float(cols[i + 2]) > args.qual):
					conflict = True
					equalAlleles = cols[i].split(",")
					for allele in equalAlleles:
						if (allele in finalAlleles):
							sampleAlleles[allele] += float(cols[i + 1])
							conflict = False
							break 
					if (conflict):
						inconsistentAlleles.append("_".join(equalAlleles + cols[i+1:i+3]))
		sampleName = ".".join(f.split("/")[-1].split(".")[0:-1])
		if (re.search("_genotype$", sampleName)):
			sampleName = sampleName[0:-9]
		output = [sampleName]
		output += [str(sampleAlleles[allele]) for allele in sorted(sampleAlleles.keys())]
		output += [",".join(inconsistentAlleles)]
		print("\t".join(output))
		fp.close()				
