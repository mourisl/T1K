#!/usr/bin/env perl

use strict ;
use warnings ;

die "Usage: perl GtfToDat.pl reference.fa annotation.gtf alleleid source > output.dat\n" if (@ARGV == 0) ;

my %genome ;
my %interestedGeneName ;

my $i ;

# Read in the reference genome.
open FP1, $ARGV[0] ;
my $chrom = "" ;
my $seq = "" ;

while ( <FP1> )
{
	if ( /^>/ )
	{
		$genome{ $chrom } = $seq if ( $chrom ne "" ) ;
		$seq = "" ;
		$chrom = substr( ( split )[0], 1 ) ;
	}
	else
	{
		chomp ;
		$seq .= $_ ;
	}
}
$genome{ $chrom } = $seq if ( $chrom ne "" ) ;
close FP1 ;

my $alleleId = "001" ;
if (defined $ARGV[2])
{
	$alleleId = $ARGV[2] ;
}

my $source = "" ;
if (defined $ARGV[3])
{
	$source = $ARGV[3] ;
}

# Obtain the exon coordinates from the gtf file
my %exons ;
open FP1, $ARGV[1] ;
my $prevTname = "-1" ;
my $strand = "." ;
my $gname = "-1" ;
my @range ;

sub GetExonsLength
{
	my @ranges = @_ ;
	my $len = 0 ;
	for (my $i = 0 ; $i < scalar(@range) ; $i += 3)
	{
		$len += $range[$i + 2] - $range[$i + 1] + 1 ;
	}
	return $len ;
}

while (<FP1>)
{
	next if ( /^#/ ) ;
	chomp ;
	my @cols = split /\t/ ;
	next if ( $cols[2] ne "exon" ) ;

	my $tname ;	
	if ( $cols[8] =~ /transcript_name \"(.*?)\"/ )
	{
		#print $1, "\n" ; 
		$tname = $1 ;
	}
	else
	{
		die "No transcript_name", $_, "\n" ;
	}

	if ( $tname ne $prevTname )
	{
		#print($interestedGeneName{"CYP2D6"}, " ", scalar(@range), "\n");
		if ( (defined $interestedGeneName{ $gname } ) && 
			 $interestedGeneName{$gname} eq "." && scalar(@range) > 0 )
		{
			$interestedGeneName{$gname}	= $strand ;
			# order the range in increasing order
			if (@range > 3 && $range[1] > $range[4])
			{
				{
					my $j ;
					for ($i = 0, $j = scalar(@range) - 3 ; $i < $j ; $i += 3, $j -= 3)
					{
						($range[$i + 1], $range[$j + 1]) = ($range[$j + 1], $range[$i + 1]) ;
						($range[$i + 2], $range[$j + 2]) = ($range[$j + 2], $range[$i + 2]) ;
					}
				}
			}
			
			if (!defined $exons{$gname} ||
				GetExonsLength(@range) > GetExonsLength(@{$exons{$gname}}))
			{
				@{$exons{$gname}} = @range ;
			}
		}

		$prevTname = $tname ;
		if ( $cols[8] =~ /gene_name \"(.*?)\"/ )
		{
			#print $1, "\n" ; 
			$gname = uc($1) ;
		}
		else
		{
			die "No gene_name: ", $_, "\n" ;
		}
		$interestedGeneName{$gname} = "." ;
		$strand = $cols[6] ;
		undef @range ;
	}
	
	push @range, $cols[0], $cols[3] - 1, $cols[4] - 1 ;			
}
close FP1 ;
 
foreach my $gname (keys %exons)
{
	my @alleleExon = @{$exons{$gname}} ;
	my $seq = "" ;
	my $padding = 500 ;
	
	my $chr = $alleleExon[0] ;
	my $start = $alleleExon[1] - $padding ;
	my $size = scalar(@alleleExon) ;
	$start = 0 if ($start < 0) ;
	my $end = $alleleExon[$size - 1] + $padding;
	$end = length($genome{$chr}) - 1 if ($end >= length($genome{$chr})) ;
	$seq = substr($genome{$chr}, $start, $end - $start + 1) ;
	
	# Move the exon coordinate according according to the padding.
	for ($i = 0 ; $i < scalar(@alleleExon) ; $i += 3)
	{
		$alleleExon[$i + 1] -= $start ;
		$alleleExon[$i + 2] -= $start ;
	}

	# Reverse
	$seq = uc($seq) ;
	my $len = length($seq) ;
	if ($interestedGeneName{$gname} eq '-')
	{
		$seq = reverse($seq) ;
		$seq =~ tr/ACGT/TGCA/ ;
		my $j ;
		for ($i = 0, $j = scalar(@alleleExon) - 3 ; $i < $j ; $i += 3, $j -= 3)
		{
			($alleleExon[$i + 1], $alleleExon[$j + 1]) = ($alleleExon[$j + 1], $alleleExon[$i + 1]) ;
			($alleleExon[$i + 2], $alleleExon[$j + 2]) = ($alleleExon[$j + 2], $alleleExon[$i + 2]) ;
		}
		
		for ($i = 0 ; $i < scalar(@alleleExon) ; $i += 3)
		{
			($alleleExon[$i + 1], $alleleExon[$i + 2]) = ($len - 1 - $alleleExon[$i + 2], 
								$len - 1 - $alleleExon[$i + 1]) ;
		}
	}

	# Output
	my $allele = $gname."*".$alleleId ;
	print("ID   $allele\n")	;
	print("DE   source $source $allele\n") if ($source ne "") ;
	print("FT   allele=\"$allele\"\n") ;
	if ($alleleExon[1] > 0)
	{
		print("FT   UTR            1..".$alleleExon[1]."\n") ;
	}
	for ($i = 0 ; $i < scalar(@alleleExon) ; $i += 3)
	{
		print("FT   exon          ".($alleleExon[$i + 1] + 1)."..".($alleleExon[$i + 2] + 1)."\n") ;
		if ($i + 3 < scalar(@alleleExon)) 
		{
			print("FT   intron        ".($alleleExon[$i + 2] + 2)."..".($alleleExon[$i + 4])."\n") ;
		}
	}
	if ($alleleExon[-1] < $len - 1)
	{
		print("FT   UTR            ".($alleleExon[-1] + 2)."..".$len."\n") ;
	}
	print("SQ  Sequence $len BP\n") ;
	print("$seq $len\n") ;
	print("//\n") ;
}
