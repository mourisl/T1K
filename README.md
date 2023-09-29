T1K
=======

Described in: 

Song, L., et al. Efficient and accurate KIR and HLA genotyping with massively parallel sequencing data. Genome Res. 2023 May 11;gr.277585.122. doi: 10.1101/gr.277585.122. 

	Copyright (c) 2021 Li Song, Bo Li, Heng Li

Includes portions copyright from: 

	samtools - Copyright (C) 2008-, Genome Research Ltd, Heng Li
	

### What is T1K?

T1K (The ONE genotyper for Kir and HLA) is a computational tool to infer the alleles for the polymorphic genes such as KIR and HLA. T1K calculates the allele abundances based on the RNA-seq/WES/WGS read alignments on the provided allele reference sequences. The abundances are used to pick the true alleles for each gene. T1K provides the post analysis steps, including novel SNP detection and single-cell representation. T1K supports both single-end and paired-end sequencing data with any read length.   

### Install

1. Clone the [GitHub repo](https://github.com/mourisl/T1K), e.g. with `git clone https://github.com/mourisl/T1K.git`
2. Run `make` in the repo directory
3. Generate the allele reference sequences of IPD-IMGT/HLA and IPD-KIR datbases. You can also find pre-built indices in the release page.
```
	perl t1k-build.pl -o hlaidx --download IPD-IMGT/HLA
	perl t1k-build.pl -o kiridx --download IPD-KIR
```
You will find the executable files in the downloaded directory. If you want to run T1K without specifying the directory, you can either add the directory of T1K to the environment variable PATH or create a soft link ("ln -s") of the file "run-t1k" to a directory in PATH.

T1K depends on [pthreads](http://en.wikipedia.org/wiki/POSIX_Threads) and samtools depends on [zlib](http://en.wikipedia.org/wiki/Zlib). 

T1K is also available from [Bioconda](https://bioconda.github.io/recipes/t1k/README.html?highlight=t1k). You can install T1K with `conda install -c bioconda t1k`.

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
		--alleleWhitelist STRING: only consider read aligned to the listed allele sereies. (default: not used)
		--barcode STRING: if -b, BAM field for barcode; if -1 -2/-u, file containing barcodes (default: not used)
		--barcodeRange INT INT CHAR: start, end(-1 for length-1), strand in a barcode is the true barcode (default: 0 -1 +)
		--barcodeWhitelist STRING: path to the barcode whitelist (default: not used)
		--read1Range INT INT: start, end(-1 for length-1) in -1/-u files for genomic sequence (default: 0 -1)
		--read2Range INT INT: start, end(-1 for length-1) in -2 files for genomic sequence (default: 0 -1)
		--mateIdSuffixLen INT: the suffix length in read id for mate. (default: not used)
		--abnormalUnmapFlag: the flag in BAM for the unmapped read-pair is nonconcordant (default: not set)
		--relaxIntronAlign: allow one more mismatch in intronic alignment (default: false)
		--preset STRING: preset parameters for cases requiring non-default settings:
			hla: HLA genotyping in general
			hla-wgs: HLA genotyping on WGS data
			kir-wgs: KIR genotyping on WGS data
			kir-wes: KIR genotyping on WES data
		--noExtraction: directly use the files from provided -1 -2/-u for genotyping (default: extraction first)
		--skipPostAnalysis: only conduct genotyping. (default: do the post analysis)
		--stage INT: start genotyping on specified stage (default: 0):
			0: start from beginning (candidate read extraction)
			1: start from genotype with candidate reads
			2: start from post analysis

* ##### User cases 
```	
	# Genotyping HLA on RNA-seq data
	./t1k -1 read_1.fq -2 read_2.fq --preset hla -f hlaidx/hlaidx_rna_seq.fa 
	# Genotyping KIR on whole genome sequencing data
	./t1k -1 read_1.fq -2 read_2.fq --preset kir-wgs -f kiridx/kiridx_dna_seq.fa 
```

### Input/Output

The primary input to T1K is the raw RNA-seq files in fasta/fastq format (-1/-2 for paired; -u for single-end; -i for interleaved), and the allele reference sequences (-f). For RNA-seq data, the user shall pick the "rna" reference file, e.g.: kiridx/kiridx_rna_seq.fa, for -f and -b option. For WES and WGS data, the user shall select the "dna" reference file for -f and -b.

The alternative input to T1K is the alignment BAM file (-b), which requires -f and the gene coordinate file (-c). To create the file for -c command, you can run command like "perl t1k-build.pl -o kiridx -d kiridx/kir.dat -g gencode.gtf" to create "_{dna,rna}_coord.fa" file.

T1K outputs several files. t1k_genotype.tsv is the main output file holding the genotyping result, where the allele for each gene is on its own line with format:

	gene_name num_diff_alleles allele_1 abundance_1 quality_1 allele_2 abundance_2 quality_2

In the case of missing alleles or homozygous alleles, the triple (allele, abundance, quality) will be ". 0 -1" as place holders. **We recommend to ignore alleles with quality less or equal to 0**. 

The other outputs files are: 

* t1k_candidate{_1/_2}.fq: the candidate reads extracted from raw data for genotyping
* t1k_aligned{_1/_2}.fq: the reads can be aligned to some alleles during genotyping 
* t1k_allele.tsv: the representative alleles with all digits and its quality score
* t1k_allele.vcf: the novel SNPs. Quality value "FAIL" represents ambiguous SNPs. The coordinates are with respected to the mRNA sequence (concatenation of the exons), even when genotyping WES/WGS data. 

### Practical notes

* #### Custom database based on VCF files such as PharmVar

Databases like PharmVar represent the variations of the alleles in the form of VCF file. T1K provides the scripts for generating the EMBL-ENA formatted dat file from VCF files. The dat file can then be used in "t1k-build" to create the reference files. Please refer to the tutorial in the vcf_database folder.

* #### Custom database with known sequences
If you have collected the linear sequences for the interested alleles, you can directly build the reference sequence. The allele name should be in the formation like "gene_name\*ABCDEFG". T1K reports the genotype at the allele series level, so by default it will report the allele "gene_name\*ABC" using three digits. If you have a special format for allele ids, you can feed the information to T1K throught the option "--alleleDigitUnits" and "--alleleDigitDelimiter". For example, the delimiter in HLA allele is ":", and if you want to genotype with first four digits (two digit units/groups), you can run T1K with options "--alleleDigitUnits 2 --alleleDigitDelimiter :".

* #### SMART-Seq data

We provide a wrapper "t1k-smartseq.pl" to process the files from platforms like SMART-seq. The user shall give the path to each file in a text file. An example command can be

	perl t1k-smartseq.pl -1 read1_list.txt -2 read2_list.txt -t 8 -f kiridx/kiridx_rna_seq.fa -o T1K

The final output file is "T1K_final_genotype.tsv". This file is a data matrix, where the rows are the cells and columns are the abundance for each identified allele.

* #### 10x Genomics data

For 10X Genomics data, usually the input is the BAM file from cell-ranger, and you can use "--barcode" to specify the field in the BAM file to specify the barcode: e.g. "--barcode CB".

If your input is raw FASTQ files, you can use "--barcode" to specify the barcode file and use "--barcodeRange" to tell T1K how to extract barcode information. If the barcode or UMI sequence is in the read sequence, you may use "--read1Range", "--read2Range" to tell T1K how to extract sequence information in the reads. T1K supports using wildcard in the -1 -2/-u option, so a typical way to run 10X Genomics single-end data is by:

	run-t1k -f kiridx/kiridx_rna_seq.fa -u path_to_10X_fastqs/*_R2_*.fastq.gz --barcode path_to_10X_fastqs/*_R1_*.fastq.gz --barcodeRange 0 15 + --barcodeWhitelist cellranger_folder/cellranger-cs/VERSION/lib/python/cellranger/barcodes/737K-august-2016.txt [other options]

The exact options depend on your 10x Genomics kit.

For barcoded file, T1K will generate the data matrix file "t1k_barcode_expr.tsv" file, where rows are the barcodes and columns are the allele abundances. Due to the shallow coverage in 10x Genomics data, the results need to be interpret with caution.

* #### Preset parameters
T1K's default parameter is based on RNA-seq data input, but it provides a series of preset parameters for different genotyping scenarios. The rationale for requiring different sets of parameters comes from two parts:

1. When the reference database is far from complete, we need a lower specificity alignment strategy to incorporate reads from a missing exact allele but may have a homologous allele in the same series (lower value for -s). Furthermore, the intron sequence in such database is even less representative, therefore we need to set "--relaxIntronAlign" to allow more variations in intron region during read alignment. 

2. When the sequencing data is WGS, the reads can come from more regions on the genomes. Therefore, we need a higher specificity alignment parameter to exclude reads from other regions (higher value for -s).

### Example

The directory './example' in this distribution contains two FASTQs as input for T1K. Run T1K with:

	./run-t1k -f kiridx/kiridx_rna_seq.fa -1 example/example_1.fq -2 example/example_2.fq -t 8 -o T1K_example

### Support

Create a [GitHub issue](https://github.com/mourisl/T1K/issues).
