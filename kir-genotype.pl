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
    "\t-b STRING: path to BAM file\n".
    "\t-1 STRING -2 STRING: path to paired-end read files\n".
    "\t-u STRING: path to single-end read file\n".
    "\t-f STRING: path to the KIR reference sequence file\n".
    "Optional:\n".
		"\t-c STRING: path to the gene coordinate file (required when -b input)\n".
    "\t-o STRING: prefix of output files. (default: inferred from file prefix)\n".
    "\t--od STRING: the directory for output files. (default: ./)\n".
    "\t-t INT: number of threads (default: 1)\n".
		"\t--frac FLOAT: filter if abundance is less than the frac of dominant allele (default: 0.15)\n".
		"\t--cov FLOAT: filter genes with average coverage less than the specified value (default: 1.0)\n".
		"\t--crossGeneRate FLOAT: the effect from other gene's expression (0.02)\n".
    "\t--barcode STRING: if -b, BAM field for barcode; if -1 -2/-u, file containing barcodes (default: not used)\n".
    "\t--barcodeRange INT INT CHAR: start, end(-1 for length-1), strand in a barcode is the true barcode (default: 0 -1 +)\n".
    "\t--barcodeWhitelist STRING: path to the barcode whitelist (default: not used)\n".
    "\t--read1Range INT INT: start, end(-1 for length-1) in -1/-u files for genomic sequence (default: 0 -1)\n".
    "\t--read2Range INT INT: start, end(-1 for length-1) in -2 files for genomic sequence (default: 0 -1)\n".
    #"\t--UMI STRING: if -b, bam field for UMI; if -1 -2/-u, file containing UMIs (default: not used)\n".
    #"\t--umiRange INT INT CHAR: start, end(-1 for lenght-1), strand in a umi is the true umi (default: 0 -1 +)\n".
    "\t--mateIdSuffixLen INT: the suffix length in read id for mate. (default: not used)\n".
    "\t--abnormalUnmapFlag: the flag in BAM for the unmapped read-pair is nonconcordant (default: not set)\n".
    "\t--noExtraction: directly use the files from provided -1 -2/-u for genotyping (default: extraction first)\n".
    "\t--stage INT: start genotyping on specified stage (default: 0):\n".
    "\t\t0: start from beginning (candidate read extraction)\n".
    "\t\t1: start from genotype with candidate reads\n".
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
my @firstMateFiles ;
my @secondMateFiles ;
my @bamFiles ;
my @barcodeFiles ;
my $prefix = "" ;
my $bamExtractorArgs = "" ;
my $fastqExtractorArgs = "" ;
my $genotyperArgs = "" ;
my $analyzerArgs = "" ;
my $threadCnt = 1 ;
my $stage = 0 ;
my $noExtraction = 0 ;
my $hasBarcode = 0 ;
my $hasUmi = 0 ;
my $outputDirectory = "" ;

my $refCoordFasta = "" ;
my $refSeqFasta = "" ;

#my $filterFrac = 0.15 ;
#my $filterCov = 1.0 ;
#my $crossGeneRate = 0.0005 ;
my %genotyperArgNames = ("--frac"=>0, "--cov"=>0, "--crossGeneRate"=>0, "-s"=>0) ;
my %analyzerArgNames = ("-s"=>0) ;
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
			push @firstMateFiles, glob($ARGV[$j]) ;
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
		$refSeqFasta = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-c" )
	{	
		$refCoordFasta = $ARGV[$i + 1] ;
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
	elsif ( $ARGV[$i] eq "--stage" )
	{
		$stage = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( defined( $genotyperArgNames{$ARGV[$i]} ) 
		|| defined ($analyzerArgNames{$ARGV[$i]}))
	{
		if ( defined( $genotyperArgNames{$ARGV[$i]} ) ) 
		{
			$genotyperArgs .= " ".$ARGV[$i]." ".$ARGV[$i + 1] ;
		}
		if ( defined( $analyzerArgNames{$ARGV[$i]} ) )
		{
			$analyzerArgs .= " ".$ARGV[$i]." ".$ARGV[$i + 1] ;
		}
		++$i ;
	}
	else
	{
		die "Unknown parameter ".$ARGV[$i]."\n" ;
	}
}

if ( @bamFiles == 0 && @firstMateFiles == 0 )
{
	die "Need to use -b/{-1,-2}/-u to specify input reads.\n" ;
}

if ( @bamFiles > 0 && $noExtraction == 1 )
{
	die "--noExtraction option can only be set when using -1 -2/-u as input.\n" ;
}

if ( $refSeqFasta eq "" )
{
	die "Need to use -f to specify the reference sequence file.\n" ;
}

if ( @bamFiles > 0 && $refCoordFasta eq "" )
{
	die "Need to use -c to specify gene coordinate file for BAM input.\n" ;
}



# Infer the output prefix.
if ( $prefix eq "" )
{
	# infer the output prefix.
	if ( @bamFiles > 0 )
	{
		$prefix = "genotype_".( split /\./, basename( $bamFiles[0] ) )[0] ;
	}
	elsif ( @firstMateFiles > 0 )
	{
		$prefix = "genotype_".( split /\./, basename( $firstMateFiles[0] ) )[0] ;
	}
	else
	{
		$prefix = "genotype" ;
	}
}

if ( $outputDirectory ne "" )
{
	mkdir $outputDirectory if ( !-d $outputDirectory ) ;
	$prefix = "$outputDirectory/$prefix" ;
}

# Extract the file
my $extractorPrefix = "${prefix}_candidate" ;
my $candidateRd1 = "${extractorPrefix}_1.fq" ;
my $candidateRd2 = "${extractorPrefix}_2.fq" ;
my $candidateRd = "${extractorPrefix}.fq" ;
my $candidateFiles = "$candidateRd1 $candidateRd2" ;
if ( $stage <= 0 )
{
	if ( @bamFiles > 0 )
	{
		system_call( "$WD/bam-extractor -b ".$bamFiles[0]." -t $threadCnt -f $refCoordFasta -o $extractorPrefix $bamExtractorArgs" ) ;
	}
	elsif ( @secondMateFiles > 0 && $noExtraction == 0 )
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
		system_call( "$WD/fastq-extractor -t $threadCnt -f $refSeqFasta -o $extractorPrefix $fastqExtractorArgs" ) ;
	}
	elsif ( @firstMateFiles > 0 && $noExtraction == 0 )
	{
		my $fname ; 
		foreach $fname (@firstMateFiles)
		{
			$fastqExtractorArgs .= " -u ".$fname ;
		}
		foreach $fname (@barcodeFiles)
		{
			$fastqExtractorArgs .= " --barcode ".$fname ;
		}
		system_call( "$WD/fastq-extractor -t $threadCnt -f $refSeqFasta -o $extractorPrefix $fastqExtractorArgs" ) ;
		$candidateFiles = "$candidateRd" ;
	}
}

# determine paired-end or single-end
if ( $noExtraction == 0 )
{
	if ( -e $candidateRd1 )
	{
		;
	}
	elsif ( -e $candidateRd )
	{
		$candidateFiles = "$candidateRd" ;
	}
	elsif ( $stage <= 1 )
	{
		die "Could not find files like ${extractorPrefix}*.fq\n" ;
	}
}
else
{
	if ( @secondMateFiles > 0 )
	{
		$candidateFiles = $firstMateFiles[0]." ".$secondMateFiles[0] ;		
	}
	elsif ( @firstMateFiles > 0 )
	{
		$candidateFiles = $firstMateFiles[0] ;
	}
}

# Run bwa
my $bwaRd1 = "${prefix}_bwa_r1.fq" ;
my $bwaRd2 = "${prefix}_bwa_r2.fq" ;
my $bwaRd = "${prefix}_bwa.fq" ;
my $bwaReadFiles = "$bwaRd1 $bwaRd2" ;

if ( 0 ) # no long needed #$stage <= 1 )
{
	#system_call("bwa mem -t $threadCnt $bwaIdx $candidateFiles > ${prefix}_bwa_aligned.sam") ;
	my @cols = split /\s/, $candidateFiles ;
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
if ( $stage <= 1 )
{
	#system_call("python3 $WD/KirGenotype.py -a ${prefix}_kallisto/abundance.tsv > ${prefix}_genotype.tsv") ;
	my @cols = split /\s/, $candidateFiles ;
	if (scalar(@cols) > 1)
	{
		system_call("$WD/genotyper $genotyperArgs -o $prefix -t $threadCnt -f $refSeqFasta -1 ".$cols[0]." -2 ".$cols[1]) ;
	}
	else
	{
		system_call("$WD/genotyper $genotyperArgs -o $prefix -t $threadCnt -f $refSeqFasta -u ".$cols[0]) ;
	}
}

if ($stage <= 2)
{
	my @cols = split /\s/, $candidateFiles ;
	if (scalar(@cols) > 1)
	{
		system_call("$WD/analyzer $analyzerArgs -o $prefix -t $threadCnt -f $refSeqFasta -a ${prefix}_allele.tsv -1 ${prefix}_aligned_1.fa -2 ${prefix}_aligned_2.fa") ;
	}
	else
	{
		system_call("$WD/analyzer $analyzerArgs -o $prefix -t $threadCnt -f $refSeqFasta -a ${prefix}_allele.tsv -u ${prefix}_aligned.fa") ;
	}
}

print STDERR "[".localtime()."] Finish.\n";
