#!/usr/bin/env python3

import sys
import argparse

def GetGenotype(geneAbund):
	majorAlleleAbund = {}
	majorGeneAbund = {}
	allelesInGene = {}
	for gene in geneAbund:
		majorGene = gene.split("*")[0]
		allele = gene[0:(len(majorGene) + 4)]
		if (allele not in majorAlleleAbund):
			majorAlleleAbund[allele] = [0, 0]
		majorAlleleAbund[allele][0] += geneAbund[gene][0]
		majorAlleleAbund[allele][1] += geneAbund[gene][1]
		if (majorGene not in majorGeneAbund):
			majorGeneAbund[majorGene] = [0, 0]
			allelesInGene[majorGene] = {}
		majorGeneAbund[majorGene][0] += geneAbund[gene][0]
		majorGeneAbund[majorGene][1] += geneAbund[gene][1]
		if (allele not in allelesInGene[majorGene]):
			allelesInGene[majorGene][allele] = 0
		allelesInGene[majorGene][allele] += 1

	genotype = {} 
	for gene in majorGeneAbund:
		l = []
		for allele in allelesInGene[gene]:
			if (majorAlleleAbund[allele][1] >= majorGeneAbund[gene][1] * 0.1):
				l.append([allele, majorAlleleAbund[allele][1]])
		l = sorted(l, key = lambda x: x[1], reverse=True)
		genotype[gene] = []
		for i in range(0, min(2, len(l))):
			genotype[gene].append(l[i][0])
	return genotype


if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(description="Genotype from IPDKIR abundance estimation")
	parser.add_argument("-a", help="abunance estimation file", dest="abundfile", required=True)
	args = parser.parse_args()
	
	geneAbund = {}

	fp = open(args.abundfile)
	for line in fp:
		if ("KIR" not in line): # header or other genes
			continue
		cols = line.rstrip().split()
		geneAbund[cols[0]] = [float(cols[3]), float(cols[4])]
	fp.close()

	genotype = GetGenotype(geneAbund)
	for gene in sorted(list(genotype.keys())):
		print(gene + "\t" + "\t".join(genotype[gene])) ;
