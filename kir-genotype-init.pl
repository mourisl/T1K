#!/usr/bin/env perl

use strict ;
use warnings ;


use Cwd 'cwd' ;
use Cwd 'abs_path' ;
use File::Basename ;

my $progName = "kir-genotype-init" ;

die "$progName usage: ./$progName [OPTIONS]:\n".
    "Required:\n".
		"\t-f STRING: IPD KIR gene sequence file\n".
		"Optional:\n".
		"\t-o STRING: output folder (default: ./)\n".
		"\t-g STRING: genome annotation file (default: not used)\n" ;


sub system_call
{
	print STDERR "[".localtime()."] SYSTEM CALL: ".join(" ",@_)."\n";
	system(@_) == 0
		or die "system @_ failed: $?";
	#print STDERR " finished\n";
} 

my $i ;

my $ipdkirFasta = "" ;
my $outputDirectory = "./" ;
my $annotationFile = "" ;

for ($i = 0 ; $i < @ARGV ; ++$i) 
{
	if ($ARGV[$i] eq "-f")
	{
		$ipdkirFasta = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ($ARGV[$i] eq "-o")
	{
		$outputDirectory = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ($ARGV[$i] eq "-g") 
	{
		$annotationFile = $ARGV[$i + 1] ;
		++$i ;
	}
	else
	{
		die "Unknown parameter ".$ARGV[$i]."\n" ;
	}
}

if ($ipdkirFasta eq "")
{
	die "Need to use -f to specify IPDKIR_nuc file.\n" ;
}

if ( !-d $outputDirectory)
{
	mkdir $outputDirectory ;
}

my $kirSeqFile = "$outputDirectory/kir_seq.fa" ;

# Reheader the IPD KIR gene sequence file
open FP, $ARGV[0] ;
open FPout, $kirSeqFile ;
while (<FP>)
{
	if (!/^>/) 
	{
		print FPout $_ ;
		next ;
	}
	chomp ;
	my @cols = split /\s/, substr($_, 1) ;
	print FPout ">".$cols[1]."\n" ;
}
close FP ;
close FPout ;

# Build BWA index
if (!-d "$outputDirectory/bwa_idx")
{
	mkdir "$outputDirectory/bwa_idx" ;
}
system_call("bwa index -p $outputDirectory/bwa_idx/bwa $kirSeqFile") ;

# Build Kallisto index
if (!-d "$outputDirectory/kallisto_idx")
{
	mkdir "$outputDirectory/kallisto_idx" ;
}
system_call("kallisto index -i $outputDirectory/kallisto_idx/kallisto $kirSeqFile") ;

# Add the genome coordinate to fasta file.
if ($annotationFile ne "")
{
	system_call("perl AddKirCoord.pl $kirSeqFile $annotationFile > $outputDirectory/kir_coord.fa") ;
}
