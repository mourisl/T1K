#!/usr/bin/env perl

use strict ;
use warnings ;
use Getopt::Long qw(GetOptions) ;


sub usage {
	die "Usage: perl ProcessMultipleGenomesToDat.pl [OPTIONS] > output.dat\n".
		"\t-g genome.fa\n".
		"\t-a reference_annotation.gtf\n".
		"\t--tmp prefix for temporary file\n"
		;
}

sub system_call
{
	print STDERR "[".localtime()."] SYSTEM CALL: ".join(" ",@_)."\n";
	system(@_) == 0
		or die "system @_ failed: $?";
	#print STDERR " finished\n";
} 

usage() if (@ARGV <= 0) ;

my $genomeFile = "" ;
my $refAnnotation = "" ;
my $tmpPrefix = "tmp" ;

GetOptions(
	"g=s" => \$genomeFile,
	"a=s" => \$refAnnotation,
) or usage() ;


my %genomes ;
my $refGenome = "" ;
my @genomeList ;

open FP, $genomeFile ;
my $seq = "" ;
my $gname = "" ;
my $gcnt = 0 ;
while (<FP>)
{
	chomp ;
	my $line = $_ ;
	if (/^>/)
	{
		$genomes{$gname} = $seq if ($seq ne "") ;
		$seq = "" ;
		$gname = substr($line, 1) ;
		$refGenome = $gname if ($gcnt == 0) ;
		push @genomeList, $gname ;
		$gcnt += 1 ;
	}
	else
	{
		$seq .= $line ;
	}
}
$genomes{$gname} = $seq if ($seq ne "") ;
close FP ;

my $i ;
my $refTempFile = "${tmpPrefix}_ref.fa" ;
my $genomeTempFile = "${tmpPrefix}_genome.fa" ;
my $annotationTempFile = "${tmpPrefix}_genome.gtf" ;

open FPref, ">$refTempFile" ;
print FPref ">$refGenome\n".$genomes{$refGenome}."\n" ;
close FPref ;

# Output the dat information for each sample
for ($i = 0 ; $i < $gcnt ; ++$i)
{
	open FPgenome, ">$genomeTempFile"	 ;
	my $gname = $genomeList[$i] ;
	print FPgenome ">".$gname."\n".$genomes{$gname}."\n" ;
	close FPgenome ;

	system_call("liftoff -g $refAnnotation $genomeTempFile $refTempFile | awk '\$2==\"Liftoff\"'> $annotationTempFile") ;
	unlink "$genomeTempFile.mmi" ;
	#system_call("rm $genomeTempFile.fai > /dev/null") ;
	my $alleleId = sprintf("%03d", $i + 1) ;
	system_call("perl GtfToDat.pl $genomeTempFile $annotationTempFile $alleleId $gname") ;
}

unlink $refTempFile, $genomeTempFile, $annotationTempFile, "$refTempFile.fai", "$genomeTempFile.fai", "$genomeTempFile.mmi" ;
