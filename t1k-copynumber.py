#!/usr/bin/env python3

import argparse
import math

geneList = []

# The likelihood of normal distribution without the constant factor
# 1/sigma * e^{-1/2*(x-mu)^2/sigma^2}
def NormalLikelihoodFactor(x, params):
	mu = params[0]
	sigma = math.sqrt(params[1])
	return math.exp(-0.5 * math.pow((x - mu)/sigma, 2)) / sigma

def LogNormalLikelihoodFactor(x, params):
	mu = params[0]
	sigma = math.sqrt(params[1])
	return -0.5 * math.pow((x - mu)/sigma, 2) - math.log(sigma)

def AbundTransform(x):
	#return x
	return math.sqrt(x)

if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(description = "Infer the allele copy number. Output directly to stdout.")
	parser.add_argument("-g", help="T1K's genotyping result file", dest="gfile", required=True)
	parser.add_argument("--nomissing", help="A comma separated list of genes that should be on every chromosome for inference one-copy parameters (will ignore the quantile options below)", dest="nomissing_list",
			required=False, default="")
	parser.add_argument("--upper-quantile", help="The upper quantile of alleles used to inference one-copy parameters (the alleles are from the predicted heterogzyous genes)", dest="upper_quantile", required=False, default=0.3)
	parser.add_argument("--lower-quantile", help="The upper quantile of alleles used to inference one-copy parameters", dest="lower_quantile", required=False, default=0)
	parser.add_argument("--adjust-var", help="Adjust variance by the given factor", dest="adjust_var", required=False, default=1.0)
	parser.add_argument("-q", help="ignore allels with less or equal quality scores", dest="qual", required=False, default=0)

	args = parser.parse_args()

	geneRank = {}
	geneToAlleles = {}
	alleleInfo = {}
	
	nomissingGenes = {}
	if (args.nomissing_list != ""):
		nomissingGenes = {g:1 for g in args.nomissing_list.split(",")}

	# Read in the allele information
	fp = open(args.gfile)
	geneIdx = 0
	alleleIdx = 0
	for line in fp:
		cols = line.rstrip().split()
		geneRank[cols[0]] = geneIdx
		geneToAlleles[cols[0]] = []
		geneIdx += 1 
		geneCopy = int(cols[1])
		for i in range(geneCopy):
			k = 2 if (i == 0) else 5 
			allele = cols[k]
			abund = float(cols[k + 1])
			quality = int(cols[k + 2])
			if (quality <= args.qual):
				continue
			alleleInfo[allele] = {}
			alleleInfo[allele]["abund"] = abund
			alleleInfo[allele]["rank"] = alleleIdx
			geneToAlleles[cols[0]].append(allele)
			alleleIdx += 1
	fp.close()
	
	abundances = []
	usedAlleles = 0
	if len(nomissingGenes) > 0:
		for g in nomissingGenes:
			if (g not in geneToAlleles):
				continue
			if (len(geneToAlleles[g]) > 1):
				for a in geneToAlleles[g]:
					abundances.append(AbundTransform(alleleInfo[a]["abund"])) 
			elif (len(geneToAlleles[g]) == 1):
				abundances.append(AbundTransform(alleleInfo[ geneToAlleles[g][0] ]["abund"]) / 2 ) 
			usedAlleles += len(geneToAlleles[g])

	start = int((len(alleleInfo) - usedAlleles) * float(args.lower_quantile))
	end = int((len(alleleInfo) - usedAlleles) * float(args.upper_quantile)) 
	#abundances = sorted([AbundTransform(alleleInfo[a]["abund"]) for a in alleleInfo])[start:end]
	
	heterAlleles = {}
	for g in geneToAlleles:
		if (g in nomissingGenes or len(geneToAlleles[g]) <= 1):
			continue
		for a in geneToAlleles[g]:
			heterAlleles[a] = 1
	abundances.extend( sorted([AbundTransform(alleleInfo[a]["abund"]) for a in heterAlleles])[start:end] )
	
	inspectAlleleCnt = len(abundances)
	# Infer the parameters 
	mean = sum(abundances)/inspectAlleleCnt
	var = sum([a*a for a in abundances]) / inspectAlleleCnt - mean * mean
	var *= float(args.adjust_var)
	#mean *= float(args.adjust_var)

	# Calculate the copy number
	for allele in alleleInfo:
		likelihoods = []
		for copy in range(8):
			likelihoods.append([copy + 1, NormalLikelihoodFactor(AbundTransform(alleleInfo[allele]["abund"]), [mean * (copy + 1), var * (copy + 1)])])
		sortLls = sorted(likelihoods, key=lambda x:x[1], reverse=True)
		alleleInfo[allele]["copy"] = sortLls[0][0]
		alleleInfo[allele]["ratio"] = sortLls[0][1] - sortLls[1][1]
	
	# Output final result
	for gene in sorted(geneRank.keys(), key=lambda x:geneRank[x]):
		line = gene + "\t" + str(len(geneToAlleles[gene]))
		for i in range(2):
			if (i < len(geneToAlleles[gene])):
				allele = geneToAlleles[gene][i]
				line += "\t%s\t%d\t%.2f"%(allele, alleleInfo[allele]["copy"], alleleInfo[allele]["ratio"])
			else:
				line += "\t.\t-1\t0"
		print(line)
