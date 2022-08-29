This tutorial is about how to create the reference sequences from VCF files, where the database represents the variations of a allele in VCF file. We will use the CYP2D6 from PharmVar as an example.

#### step 0: prerequisite files
You will need hg38 human reference genome (hg38.fa) and the gene annotation file such as from gencode (gencode.gtf).

#### step 1: download and process the VCF files

1.1 Click the "Download Complete Database" button at [https://www.pharmvar.org/download](https://www.pharmvar.org/download). 

1.2 Uncompress the pharmvar-XXX.zip file to the {T1K_PATH}/vcf_database/, and make {T1K_PATH}/vcf_database/ your current folder. You shall see the folder ./pharmvar-XXX/CYP2D6/ there. Put the hg38 VCF file names to the file by

	ls ./pharmvar-XXX/CY2D6/GRCh38/*.vcf > vcflist.out

1.3 Generate the combined VCF file
	
	perl ./CombinedVcf.pl CYP2D6.1 vcflist.out > cyp2d6_combined.vcf

	We need the first parameter "CYP2D6.1" because that the VCF files does not contain the primary allele "CYP2D6.1", and we need a place holder for it in the combined VCF file. 

#### step 2: create the reference files

2.1 Generate the EMBL-ENA format dat file by running:

	perl ./CombinedVcfToDat.pl genome.fa gencode.gtf cyp2d6_combined.vcf > cyp2d6.dat

	The dat file should look similar to the dat file from IPD-IMGT/HLA and IPD-KIR.

2.2 Generate the reference files

	perl {T1K_PATH}/t1k-build.pl -d cyp2d6 -g gencode.gtf -o cyp2d6_idx --prefix cyp2d6 

#### step 3: running T1K for genotyping

You can run T1K to genotype CYP2D6 with RNA-seq data using command:

	{T1K_PATH}/run-t1k -f cyp2d6_idx/cyp2d6_rna_seq.fa -1 read_1.fq -2 read_2.fq --alleleDigitUnits 1 --alleleDelimiter . -t 8

	The command will genotype CYP2D6 to the digit before "." symbol, i.e. "CYP2D6*1".
