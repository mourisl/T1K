#!/usr/bin/env perl

use strict ;
use warnings ;

use Cwd 'cwd' ;
use Cwd 'abs_path' ;
use File::Basename ;

die "T1K SMART-seq pipeline usage: perl pl [OPTIONS]:\n".
    "\t-1 STRING: file containing the list of read 1 (or single-end) files\n".
    "\t-2 STRING: file containing the list of read 2 files\n".
    "\t-f STRING: path to reference sequence file\n".
    "\t-o STRING: prefix of final output files. (default: T1K)\n".
    "\t-t INT: number of threads (default: 1)\n".
    "Other T1K parameters will be directly passed to T1K main program\n"
		#"\t--noclear: do not clear the intermediate results (default: clear)\n"
    if (@ARGV == 0) ;

sub system_call
{
	print STDERR "[".localtime()."] SYSTEM CALL: ".join(" ",@_)."\n";
	system(@_) == 0
		or die "system @_ failed: $?";
	#print STDERR " finished\n";
} 

my $WD = dirname( abs_path( $0 ) ) ;

my $i ;
my $readFile1 = "" ;
my $readFile2 = "" ;
my $outputPrefix = "T1K" ;
my $hasMate = 0 ;
my $t1kArgs = "" ;
my $referenceFile = "" ;

for ( $i = 0 ; $i < @ARGV ; ++$i )
{
	if ($ARGV[$i] eq "-1")	
	{
		$readFile1 = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ($ARGV[$i] eq "-2")
	{
		$readFile2 = $ARGV[$i + 1] ;
		++$i ;
		$hasMate = 1 ;
	}
	elsif ($ARGV[$i] eq "-o")
	{
		$outputPrefix = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ($ARGV[$i] eq "-f")
	{
		$referenceFile = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ($ARGV[$i] eq "-t")
	{
		$t1kArgs .= " ".$ARGV[$i]." ".$ARGV[$i + 1] ;
		++$i ;
	}
	elsif ($ARGV[$i] =~ /^-/)
	{
		$t1kArgs .= " ".$ARGV[$i] ;
		if (!($ARGV[$i + 1] =~ /^-/))
		{
			$t1kArgs .= " ".$ARGV[$i + 1] ;
			++$i ;
		}
	}
	else
	{
		die "Unknown parameter ".$ARGV[$i]."\n" ;
	}
}

die "Need to use -1 to specify the list of read 1 files.\n" if ($readFile1 eq "") ;

# Process the fastq file for each cell one by one
open FP1, $readFile1 ;
my $FP2 ;
open $FP2, $readFile2 if ($hasMate) ;
open FPfilelist, ">${outputPrefix}_genotype_list.out" ;
my @cells ;

while (<FP1>)
{
	chomp ;
	my $file1 = $_ ;
	my $file2 = "" ;
	
	my $fname = (fileparse($file1))[0] ;
	my $cellPrefix = (split /\./, $fname)[0] ; # use the content before the first "." as identifier.
	my $cellOutputFolder = "${outputPrefix}_${cellPrefix}" ;
	if ($hasMate)
	{
		$file2 = <$FP2> ;
		chomp $file2 ;
		system("$WD/run-t1k $t1kArgs -f $referenceFile -1 $file1 -2 $file2 --od $cellOutputFolder -o $cellPrefix") ;
	}
	else
	{
		system("$WD/run-t1k $t1kArgs -f $referenceFile -u $file1 --od $cellOutputFolder -o $cellPrefix") ;
	}

	print FPfilelist "${cellOutputFolder}/${cellPrefix}_genotype.tsv\n" ; 
	push @cells, $cellPrefix ;
}
close FPfilelist ;
close FP1 ;
close $FP2 if ($hasMate) ;

# Reduce the reference sequence 
my $cellCount = scalar(@cells) ;
my $qualityFilter = $cellCount / 2 ;
$qualityFilter = 30 if ($qualityFilter < 30) ;
my $mergedGenotypeFile = "${outputPrefix}_merged_genotype.tsv" ;
system("python3 $WD/t1k-merge.py --tq $qualityFilter -l ${outputPrefix}_genotype_list.out > $mergedGenotypeFile") ;
my $reducedReferenceFile = "${outputPrefix}_reduced_ref.fa" ;

open FP1, $mergedGenotypeFile ;
my $header = <FP1> ;
close FP1 ;

my %selectedAlleles ;
foreach my $allele (split /\s/, $header)
{
	next if ($allele eq "sample" || $allele eq "inconsistency") ;
	$allele =~ s/\*/\\\*/g ;
	$selectedAlleles{$allele} = 1 ;
}

open FPref, $referenceFile ;
open FPreducedref, ">$reducedReferenceFile" ;
while (<FPref>)
{
	my $header = $_ ;
	my $seq = <FPref> ;
	my $found = 0 ;
	foreach my $allele (keys %selectedAlleles)
	{
		if ($header =~ /$allele/)
		{
			$found = 1 ;
			last ;
		}
	}
	print FPreducedref "$header$seq" if ($found == 1) ;
}
close FPref ;
close $reducedReferenceFile ;

# Redo the genotyping
open FPfilelist, ">${outputPrefix}_reduced_genotype_list.out" ;
foreach my $cellPrefix (@cells)
{
	my $cellOutputFolder = "${outputPrefix}_${cellPrefix}" ;
	my $file1 = "${cellOutputFolder}/${cellPrefix}_candidate.fq" ;
	my $file2 = "" ;
	if ($hasMate)
	{
		$file1 = "${cellOutputFolder}/${cellPrefix}_candidate_1.fq" ;
		$file2 = "${cellOutputFolder}/${cellPrefix}_candidate_2.fq" ;
	}
	
	if ($hasMate)
	{
		system("$WD/run-t1k $t1kArgs -f $reducedReferenceFile -1 $file1 -2 $file2 --od $cellOutputFolder -o ${cellPrefix}_reduced --noExtraction") ;
	}
	else
	{
		system("$WD/run-t1k $t1kArgs -f $reducedReferenceFile -u $file1 --od $cellOutputFolder -o ${cellPrefix}_reduced --noExtraction") ;
	}
	print FPfilelist "${cellOutputFolder}/${cellPrefix}_reduced_genotype.tsv\n" ;
}
close FPfilelist ;

system("python3 $WD/t1k-merge.py --tq $qualityFilter -l ${outputPrefix}_reduced_genotype_list.out > ${outputPrefix}_final_genotype.tsv") ;
