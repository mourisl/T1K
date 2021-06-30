#!/usr/bin/env perl

use strict ;
use warnings ;

use Cwd 'cwd' ;
use Cwd 'abs_path' ;
use File::Basename ;

my $progName = "kir-genotype" ;

die "$progName usage: ./$progName [OPTIONS]:\n".
    "Required:\n".
    #"\t[Input]:\n".
    "\t-b STRING: path to bam file\n".
    "\t-1 STRING -2 STRING: path to paired-end read files\n".
    "\t-u STRING: path to single-end read file\n".
    "\t-f STRING: folder to the KIR annotation, BWA, kallisto index\n".
    "Optional:\n".
    "\t-o STRING: prefix of output files. (default: inferred from file prefix)\n".
    "\t--od STRING: the directory for output files. (default: ./)\n".
		"\t--mode STRING: \"rna\"-based or \"dna\"-based (default: rna)\n".
    "\t-t INT: number of threads (default: 1)\n".
    "\t--barcode STRING: if -b, bam field for barcode; if -1 -2/-u, file containing barcodes (default: not used)\n".
    "\t--barcodeRange INT INT CHAR: start, end(-1 for length-1), strand in a barcode is the true barcode (default: 0 -1 +)\n".
    "\t--barcodeWhitelist STRING: path to the barcode whitelist (default: not used)\n".
    "\t--read1Range INT INT: start, end(-1 for length-1) in -1/-u files for genomic sequence (default: 0 -1)\n".
    "\t--read2Range INT INT: start, end(-1 for length-1) in -2 files for genomic sequence (default: 0 -1)\n".
    #"\t--UMI STRING: if -b, bam field for UMI; if -1 -2/-u, file containing UMIs (default: not used)\n".
    #"\t--umiRange INT INT CHAR: start, end(-1 for lenght-1), strand in a umi is the true umi (default: 0 -1 +)\n".
    "\t--mateIdSuffixLen INT: the suffix length in read id for mate. (default: not used)\n".
    "\t--abnormalUnmapFlag: the flag in BAM for the unmapped read-pair is nonconcordant (default: not set)\n".
    "\t--noExtraction: directly use the files from provided -1 -2/-u to assemble (default: extraction first)\n".
    "\t--stage INT: start TRUST4 on specified stage (default: 0):\n".
    "\t\t0: start from beginning (candidate read extraction)\n".
    "\t\t1: start from bwa validation\n".
    "\t\t2: start from genotype with kallisto result\n".
    "" 
	if ( @ARGV == 0 ) ;

sub system_call
{
	print STDERR "[".localtime()."] SYSTEM CALL: ".join(" ",@_)."\n";
	system(@_) == 0
		or die "system @_ failed: $?";
	#print STDERR " finished\n";
} 

my $WD = dirname( abs_path( $0 ) ) ;
my $i ;
my $j ;


# process the options.
my @singleFiles ;
my @firstMateFiles ;
my @secondMateFiles ;
my $refFolder = "" ;
my @bamFiles ;
my @barcodeFiles ;
my $prefix = "" ;
my $bamExtractorArgs = "" ;
my $fastqExtractorArgs = "" ;
my $genotypeArgs = "" ;
my $threadCnt = 1 ;
my $stage = 0 ;
my $noExtraction = 0 ;
my $hasBarcode = 0 ;
my $hasUmi = 0 ;
my $outputDirectory = "" ;
my $mode = "rna" ;

print STDERR "[".localtime()."] $progName begins.\n" ;
for ( $i = 0 ; $i < @ARGV ; ++$i )
{
	if ( $ARGV[$i] eq "-1" )
	{
		for ($j = $i + 1; $j < @ARGV; ++$j )		
		{
			last if ($ARGV[$j] =~ /^-/) ;
			push @firstMateFiles, glob($ARGV[$j]) ;
		}
		$i = $j - 1 ;
	}
	elsif ( $ARGV[$i] eq "-2" )
	{	
		for ($j = $i + 1; $j < @ARGV; ++$j )		
		{
			last if ($ARGV[$j] =~ /^-/) ;
			push @secondMateFiles, glob($ARGV[$j]) ;
		}
		$i = $j - 1 ;
	}
	elsif ( $ARGV[ $i ] eq "-u" ) 
	{
		for ($j = $i + 1; $j < @ARGV; ++$j )		
		{
			last if ($ARGV[$j] =~ /^-/) ;
			push @singleFiles, glob($ARGV[$j]) ;
		}
		$i = $j - 1 ;
	}
	elsif ( $ARGV[$i] eq "-b" )
	{
		push @bamFiles, $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-f" )
	{	
		$refFolder = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-o" )
	{
		$prefix = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "--od" )
	{
		$outputDirectory = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-t" )
	{
		$threadCnt = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "--abnormalUnmapFlag" )
	{
		$bamExtractorArgs .= " -u" ;
	}
	elsif ( $ARGV[$i] eq "--mateIdSuffixLen" )
	{
		$bamExtractorArgs .= "--mateIdSuffixLen ".$ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "--noExtraction" )
	{
		$noExtraction = 1 ;
	}
	elsif ( $ARGV[$i] eq "--barcode" )
	{
		$hasBarcode = 1 ;
		$bamExtractorArgs .= " --barcode ".$ARGV[$i + 1] ;
		#$fastqExtractorArgs .= " --barcode ".$ARGV[$i + 1] ;
		for ($j = $i + 1; $j < @ARGV; ++$j )		
		{
			last if ($ARGV[$j] =~ /^-/) ;
			push @barcodeFiles, glob($ARGV[$j]) ;
		}
		$i = $j - 1 ;
	}
	elsif ( $ARGV[$i] eq "--barcodeRange" )
	{
		$fastqExtractorArgs .= " --barcodeStart ".$ARGV[$i + 1]." --barcodeEnd ".$ARGV[$i + 2] ;
		if ( $ARGV[$i + 3] eq "-" ) 
		{
			$fastqExtractorArgs .= " --barcodeRevComp" ;
		}

		$i += 3 ;
	}
	elsif ( $ARGV[$i] eq "--read1Range")
	{
		$fastqExtractorArgs .= " --read1Start ".$ARGV[$i + 1]." --read1End ".$ARGV[$i + 2] ;
		$i += 2 ;
	}
	elsif ( $ARGV[$i] eq "--read2Range")
	{
		$fastqExtractorArgs .= " --read2Start ".$ARGV[$i + 1]." --read2End ".$ARGV[$i + 2] ;
		$i += 2 ;
	}
	elsif ( $ARGV[$i] eq "--barcodeWhitelist" )
	{
		$fastqExtractorArgs .= " --barcodeWhitelist ".$ARGV[$i + 1] ;
		$i += 1 ;
	}
	elsif ( $ARGV[$i] eq "--UMI" )
	{
		$hasUmi = 1 ;
		$bamExtractorArgs .= " --UMI ".$ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "--mode")
	{
		$mode = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "--stage" )
	{
		$stage = $ARGV[$i + 1] ;
		++$i ;
	}
	else
	{
		die "Unknown parameter ".$ARGV[$i]."\n" ;
	}
}

if ( @bamFiles == 0 && @firstMateFiles == 0 && @singleFiles == 0 )
{
	die "Need to use -b/{-1,-2}/-u to specify input reads.\n" ;
}

if ( @bamFiles > 0 && $noExtraction == 1 )
{
	die "--noExtraction option can only be set when using -1 -2/-u as input.\n" ;
}

if ( $refFolder eq "" )
{
	die "Need to use -f to specify the folder for annotation/index files.\n" ;
}

my $kirCoordFasta = "$refFolder/kir_${mode}_coord.fa" ;
my $kirseqFasta = "$refFolder/kir_${mode}_seq.fa" ;
my $bwaIdx = "$refFolder/bwa_idx/bwa_${mode}" ;


# Infer the output prefix.
if ( $prefix eq "" )
{
	# infer the output prefix.
	if ( @bamFiles > 0 )
	{
		$prefix = "kir_".( split /\./, basename( $bamFiles[0] ) )[0] ;
	}
	elsif ( @firstMateFiles > 0 )
	{
		$prefix = "kir_".( split /\./, basename( $firstMateFiles[0] ) )[0] ;
	}
	elsif ( @singleFiles > 0 )
	{
		$prefix = "kir_".( split /\./, basename( $singleFiles[0] ) )[0] ;
	}
	else
	{
		$prefix = "kir" ;
	}
}

if ( $outputDirectory ne "" )
{
	mkdir $outputDirectory if ( !-d $outputDirectory ) ;
	$prefix = "$outputDirectory/$prefix" ;
}

# Extract the file
my $extractorPrefix = "${prefix}_possible" ;
my $possibleRd1 = "${extractorPrefix}_1.fq" ;
my $possibleRd2 = "${extractorPrefix}_2.fq" ;
my $possibleRd = "${extractorPrefix}.fq" ;
my $possibleFiles = "$possibleRd1 $possibleRd2" ;
if ( $stage <= 0 )
{
	if ( @bamFiles > 0 )
	{
		system_call( "$WD/bam-extractor -b ".$bamFiles[0]." -t $threadCnt -f $kirCoordFasta -o $extractorPrefix $bamExtractorArgs" ) ;
	}
	elsif ( @firstMateFiles > 0 && $noExtraction == 0 )
	{
		my $fname ; 
		foreach $fname (@firstMateFiles)
		{
			$fastqExtractorArgs .= " -1 ".$fname ;
		}
		foreach $fname (@secondMateFiles)
		{
			$fastqExtractorArgs .= " -2 ".$fname ;
		}
		foreach $fname (@barcodeFiles)
		{
			$fastqExtractorArgs .= " --barcode ".$fname ;
		}
		system_call( "$WD/fastq-extractor -t $threadCnt -f $kirseqFasta -o $extractorPrefix $fastqExtractorArgs" ) ;
	}
	elsif ( @singleFiles > 0 && $noExtraction == 0 )
	{
		my $fname ; 
		foreach $fname (@singleFiles)
		{
			$fastqExtractorArgs .= " -u ".$fname ;
		}
		foreach $fname (@barcodeFiles)
		{
			$fastqExtractorArgs .= " --barcode ".$fname ;
		}
		system_call( "$WD/fastq-extractor -t $threadCnt -f $kirseqFasta -o $extractorPrefix $fastqExtractorArgs" ) ;
		$possibleFiles = "$possibleRd" ;
	}
}

# determine paired-end or single-end
if ( $noExtraction == 0 )
{
	if ( -e $possibleRd1 )
	{
		;
	}
	elsif ( -e $possibleRd )
	{
		$possibleFiles = "$possibleRd" ;
	}
	elsif ( $stage <= 1 )
	{
		die "Could not find files like ${extractorPrefix}*.fq\n" ;
	}
}
else
{
	if ( @firstMateFiles > 0 )
	{
		$possibleFiles = $firstMateFiles[0]." ".$secondMateFiles[0] ;		
	}
	elsif ( @singleFiles > 0 )
	{
		$possibleFiles = $singleFiles[0] ;
	}
}

# Run bwa
my $bwaRd1 = "${prefix}_bwa_r1.fq" ;
my $bwaRd2 = "${prefix}_bwa_r2.fq" ;
my $bwaRd = "${prefix}_bwa.fq" ;
my $bwaReadFiles = "$bwaRd1 $bwaRd2" ;

if ( $stage <= 1 )
{
	system_call("bwa mem -t $threadCnt $bwaIdx $possibleFiles > ${prefix}_bwa_aligned.sam") ;
	my @cols = split /\s/, $possibleFiles ;
	if (scalar(@cols) > 1)
	{
		system_call("perl $WD/ExtractBamHits.pl ${prefix}_bwa_aligned.sam ".$cols[0]."> $bwaRd1") ;
		system_call("perl $WD/ExtractBamHits.pl ${prefix}_bwa_aligned.sam ".$cols[1]."> $bwaRd2") ;
	}
	else
	{
		system_call("perl $WD/ExtractBamHits.pl ${prefix}_bwa_aligned.sam ".$cols[0]."> $bwaRd") ;
		$bwaReadFiles = "$bwaRd" ;
	}
}

# Obtain the genotype
if ( $stage <= 2 )
{
	#system_call("python3 $WD/KirGenotype.py -a ${prefix}_kallisto/abundance.tsv > ${prefix}_genotype.tsv") ;
	system_call("$WD/genotyper -t $threadCnt -f $kirseqFasta -1 $bwaRd1 -2 $bwaRd2 > ${prefix}_genotype.tsv") ;
}

print STDERR "[".localtime()."] Finish.\n";
