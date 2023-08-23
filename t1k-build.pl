#!/usr/bin/env perl

use strict ;
use warnings ;


use Cwd 'cwd' ;
use Cwd 'abs_path' ;
use File::Basename ;
use File::Path 'make_path' ;

my $progName = "t1k-build.pl" ;

die "$progName usage: ./$progName [OPTIONS]:\n".
    "Required:\n".
		"\t-d STRING: EMBL-ENA dat file\n".
		"\t\tOr\n".
		"\t-f STRING: plain gene sequence file\n".
		"\t\tOr\n".
		"\t--download STRING: IPD-IMGT/HLA or IPD-KIR or user-specified dat file download link\n".
		"Optional:\n".
		"\t-o STRING: output folder (default: ./)\n".
		"\t-g STRING: genome annotation file (default: not used)\n".
		"\t--target STRING: gene name keyword (default: no filter)\n".
		"\t--prefix STRING: file prefix (default: based on --target or -o)\n"
		if (@ARGV == 0);


sub system_call
{
	print STDERR "[".localtime()."] SYSTEM CALL: ".join(" ",@_)."\n";
	system(@_) == 0
		or die "system @_ failed: $?";
	#print STDERR " finished\n";
} 

my $WD = dirname( abs_path( $0 ) ) ;

my $i ;

my $ipdFasta = "" ;
my $ipdDat = "" ;
my $outputDirectory = "./" ;
my $annotationFile = "" ;
my $downloadPath = "" ;
my $targetGene = "" ;
my $outputPrefix = "" ;

for ($i = 0 ; $i < @ARGV ; ++$i) 
{
	if ($ARGV[$i] eq "-f")
	{
		$ipdFasta = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ($ARGV[$i] eq "-d")
	{
		$ipdDat = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ($ARGV[$i] eq "--download")
	{
		$downloadPath = $ARGV[$i + 1] ;
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
	elsif ($ARGV[$i] eq "--target")
	{
		$targetGene = lc($ARGV[$i + 1]) ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "--prefix")
	{
		$outputPrefix = $ARGV[$i + 1] ;
		++$i ;
	}
	else
	{
		die "Unknown parameter ".$ARGV[$i]."\n" ;
	}
}

if ($ipdFasta eq "" && $ipdDat eq "" && $downloadPath eq "")
{
	die "Need to use -d/-f/--download to specify dat file, sequence file or dat file download link.\n" ;
}

if ( !-d $outputDirectory)
{
	make_path $outputDirectory ;
}

if ($ipdDat eq "" and $downloadPath ne "")
{
	if (uc($downloadPath) eq "IPD-IMGT/HLA")	
	{
		system_call("curl -o $outputDirectory/hla.dat http://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla.dat") ;
		$ipdDat = "$outputDirectory/hla.dat" ;
	}
	elsif (uc($downloadPath) eq "IPD-KIR")
	{
		system_call("curl -o $outputDirectory/kir.dat https://ftp.ebi.ac.uk/pub/databases/ipd/kir/kir.dat") ;
		$ipdDat = "$outputDirectory/kir.dat"
	}
	else
	{
		system_call("curl -o $outputDirectory/t1k_ref.dat $downloadPath") ;
		$ipdDat = "$outputDirectory/t1k_ref.dat" ;
	}
}

if ($outputPrefix eq "")
{
	if ($targetGene ne "")
	{
		$outputPrefix = $targetGene ;
	}
	elsif ($outputDirectory ne "./" )
	{
		$outputPrefix = (split /\//, $outputDirectory)[0] ;
	}
	else
	{
		$outputPrefix = "T1K_ref" ;
	}
}

my $rnaSeqFile = "$outputDirectory/${outputPrefix}_rna_seq.fa" ;
my $dnaSeqFile = "$outputDirectory/${outputPrefix}_dna_seq.fa" ;

if ($ipdDat ne "")
{
	my $options = "" ;
	$options .= " --gene $targetGene" if ($targetGene ne "") ;

	system_call("perl $WD/ParseDatFile.pl $ipdDat --mode dna $options > $dnaSeqFile") ;
	system_call("perl $WD/ParseDatFile.pl $ipdDat --mode rna $options > $rnaSeqFile") ;
}
else
{
	# Reheader the IPD gene sequence file
	open FP, $ipdFasta ;
	open FPout, ">$rnaSeqFile" ;
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
}

# Add the genome coordinate to fasta file.
if ($annotationFile ne "")
{
	system_call("perl $WD/AddGeneCoord.pl $rnaSeqFile $annotationFile > $outputDirectory/${outputPrefix}_rna_coord.fa") ;
	if ($ipdDat ne "")
	{
		system_call("perl $WD/AddGeneCoord.pl $dnaSeqFile $annotationFile > $outputDirectory/${outputPrefix}_dna_coord.fa") ;
	}
}
