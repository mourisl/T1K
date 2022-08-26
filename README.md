T1K
=======

Described in: 

Song, L., et al. T1K: efficient and accurate KIR and HLA genotyping with next-generation sequencing data

	Copyright (C) 2021- and GNU GPL by Li Song, Heng Li

Includes portions copyright from: 

	samtools - Copyright (C) 2008-, Genome Research Ltd, Heng Li
	

### What is T1K?

T1K (The ONE genotyper for Kir and HLA) is a computational tool to infer the alleles for the polymorphic genes such as KIR and HLA. T1K calculates the allele abundances based on the RNA-seq/WES/WGS read alignments on the provided allele reference sequences. The abundances are used to pick the true alleles for each gene. T1K provides the post analysis steps, including novel SNP detection and single-cell representation. T1K supports both single-end and paired-end sequencing data with any read length.   

### Install

1. Clone the [GitHub repo](https://github.com/liulab-dfci/T1K), e.g. with `git clone https://github.com/liulab-dfci/T1K.git`
2. Run `make` in the repo directory
3. Run `tar -xzf ipd_ref.tar.xzf` in the repo directory to decompress the pre-downloaded HLA and KIR sequences
4. Optional: If you want to update the allele reference sequences of IPD-IMGT/HLA and IPD-KIR from previous step, please run the

You will find the executable files in the downloaded directory. If you want to run T1K without specifying the directory, you can either add the directory of T1K to the environment variable PATH or create a soft link ("ln -s") of the file "run-t1k" to a directory in PATH.

T1K depends on [pthreads](http://en.wikipedia.org/wiki/POSIX_Threads) and samtools depends on [zlib](http://en.wikipedia.org/wiki/Zlib). 


### Usage

	Usage: ./run-t1k [OPTIONS]
	Required:
		-1 STRING -2 STRING: path to paired-end read files
			or
		-u STRING: path to single-end read file
			or
		-i STRING: path to interleaved read file
			or
		-b STRING: path to BAM file
		-f STRING: path to the reference sequence file
	Optional:
		-c STRING: path to the gene coordinate file (required when -b input)
		-o STRING: prefix of output files. (default: inferred from file prefix)
		--od STRING: the directory for output files. (default: ./)
		-t INT: number of threads (default: 1)
		-s FLOAT: minimum alignment similarity (default: 0.8)
		--frac FLOAT: filter if abundance is less than the frac of dominant allele (default: 0.15)
		--cov FLOAT: filter genes with average coverage less than the specified value (default: 1.0)
		--crossGeneRate FLOAT: the effect from other gene's expression (0.04)
		--alleleDigitUnits INT: the number of units in genotyping result. (default: automatic)
		--alleleDelimiter CHR: the delimiter character for digit unit. (default: automatic)
		--barcode STRING: if -b, BAM field for barcode; if -1 -2/-u, file containing barcodes (default: not used)
		--barcodeRange INT INT CHAR: start, end(-1 for length-1), strand in a barcode is the true barcode (default: 0 -1 +)
		--barcodeWhitelist STRING: path to the barcode whitelist (default: not used)
		--read1Range INT INT: start, end(-1 for length-1) in -1/-u files for genomic sequence (default: 0 -1)
		--read2Range INT INT: start, end(-1 for length-1) in -2 files for genomic sequence (default: 0 -1)
		--mateIdSuffixLen INT: the suffix length in read id for mate. (default: not used)
		--abnormalUnmapFlag: the flag in BAM for the unmapped read-pair is nonconcordant (default: not set)
		--relaxIntronAlign: allow one more mismatch in intronic alignment (default: false)
		--preset STRING: preset parameters for cases requiring non-default settings:
			hla: HLA genotyping
			kir-wgs: KIR genotyping on WGS data
			kir-wes: KIR genotyping on WES data
		--noExtraction: directly use the files from provided -1 -2/-u for genotyping (default: extraction first)
		--skipPostAnaysis: only conduct genotyping. (default: conduct the post analysis)
		--stage INT: start genotyping on specified stage (default: 0):
			0: start from beginning (candidate read extraction)
			1: start from genotype with candidate reads
			2: start from post analysis

### Input/Output

The primary input to T1K is the raw RNA-seq files in fasta/fastq format (-1/-2 for paired; -u for single-end; -i for interleaved), and the allele reference sequences (-f). The alternative input to T1K is the alignment BAM file (-b), which requires -f and the gene coordinate file (-b). For RNA-seq data, the user shall pick the "rna" reference file, e.g.: kiridx/kir_rna_seq.fa, for -f and -b option. For WES and WGS data, the user shall select the "dna" reference file for -f and -b.

T1K outputs several files: 

* t1k_genotype.tsv is the main output file holding the genotyping result, where the allele for each gene is on its own line with format

	gene_name num_alleles allele_1 abundance_1 quality_1 allele_2 abundance_2 quality_2

The other outputs files are: 

* t1k_candidate{_1/_2}.fq: the candidate reads extracted from raw data for genotyping
* t1k_aligned{_1/_2}.fq: the reads can be aligned to some alleles during genotyping 
* t1k_allele.tsv: the representative alleles with all digits and its quality score
* t1k_allele.vcf: the novel SNPs. Quality value "FAIL" represents ambiguous SNPs. The coordinates are with respected to the mRNA sequence (concatenation of the exons), even when genotyping WES/WGS data. If your reference does not contain

### Practical notes

* #### Update IPD-KIR and IPD-IMGT/HLA database


* #### Custom database based on VCF files, e.g. PharmVar
Please refer to the tutorial in the vcf_database folder.


* #### SMART-Seq data

We provide a wrapper "t1k-smartseq.pl" to process the files from platforms like SMART-seq. The user shall give the path to each file in a text file. An example command can be

	perl t1k-smartseq.pl -1 read1_list.txt -2 read2_list.txt -t 8 -f kiridx/kir_rna_seq.fa -o T1K
 
### Example

The directory './example' in this distribution contains two FASTQs as input for T1K. Run T1K with:

	./run-t1k -f kiridx/kir_rna_seq.fa -1 example/example_1.bam -2 example/example_2.fa -t 8 -o T1K_example

### Support

Create a [GitHub issue](https://github.com/mourisl/T1K/issues).
