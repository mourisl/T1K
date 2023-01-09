This tutorial is about how to create the reference sequences from Human Pangenome Reference Consortirum (HPRC) genomes. We will use the C4(C4A, C4B) genes from downloaded from https://zenodo.org/record/6617246#.Y1ih9-zMLzc (C4-96.agc) as an example.

#### step 0: prerequisite files
You will need the gene annotation file such as from gencode (gencode.gtf) and the program liftoff (https://github.com/agshumate/Liftoff). For this particular example, you need agc to decompress the phased genomes of HPRC.

#### step 1: organize HPRC genome seqeuences 

1.1 Download the C4-96.agc file from . 

1.2 Uncompress the C4-96.agc file by running: 

	agc getcol C4-96.agc > {T1K_PATH}/hprc_database/C4-96.fa

Make the {T1K_PATH}/hprc_database/ your current folder. 

Note that the first sequence in the fasta file (C4-96.fa here) has to be the geneome corresponding to the gene annotation GTF file. 

#### step 2: modify gene annotation GTF file  

2.1 Use manual regular pattern match find that the C4 region for the GRCh38 genome in C4-96.fa file starts at chr6:31932664 (1-based) 

2.2 Modify the GTF file to be compatible with GRCh38 genome in C4-96.fa . 

	grep "\"C4A\"\|\"C4B\"" gencode.gtf | awk -F'\t' 'BEGIN {OFS=FS} {$4-= 31932663; $5 -= 31932663; print $0}' > C4_shift.gtf

#### step 3: create the reference file

3.1 Generate the EMBL-ENA format dat file by running:

	perl ProcessMultipleGenomesToDat.pl -g C4-96.fa -a C4_shift.gtf > C4-96.dat 

This step will create a dummy allele id for C4A and C4B, respectively. For example, it assigns C4A*001 and C4B*001 for the first genome in C4-96.fa), and use 002, 003, so on so forth.

The origin phased genome id for each dummy allele is recorded in the dat file, and you can find the mappings by:
	
	grep source C4-96.dat | awk '{print $3,$4}'

, where the first column is the phased genome id and the second column is the gene name with allele id.  

3.2 Generate the reference files

	perl {T1K_PATH}/t1k-build.pl -d C4-96.dat -g genome.gtf -o C4_idx --prefix C4

#### step 4: running T1K for genotyping

You can run T1K to genotype C4 with RNA-seq data using command:

	{T1K_PATH}/run-t1k -f C4_idx/C4_rna_seq.fa -1 read_1.fq -2 read_2.fq -t 8
