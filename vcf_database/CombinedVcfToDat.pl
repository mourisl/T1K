#!/usr/bin/env perl

use strict ;
use warnings ;

die "Usage: perl CombinedVcfToDat.pl reference.fa annotation.gtf combined.vcf > output.dat\n" if (@ARGV == 0) ;

my %genome ;
my %interestedGeneName ;

my $i ;

# Read in the reference genome.
open FP1, $ARGV[0] ;
my $chrom = "" ;
my $seq = "" ;
my $hasChrPrefix = 0 ;
while ( <FP1> )
{
	if ( /^>/ )
	{
		$genome{ $chrom } = $seq if ( $chrom ne "" ) ;
		$seq = "" ;
		$chrom = substr( ( split )[0], 1 ) ;
		if ( $chrom =~ /^c/ )
		{
			$hasChrPrefix = 1 ;
		}
	}
	else
	{
		chomp ;
		$seq .= $_ ;
	}
}
$genome{ $chrom } = $seq if ( $chrom ne "" ) ;
close FP1 ;

# Read in the gene names from the combined vcf file
open FP1, $ARGV[2] ;
my %vcf ;
while (<FP1>)
{
	next if (/^#/) ;
	chomp ;
	my @cols = split ;
	my $gene = (split /\*/, $cols[0])[0] ;
	$interestedGeneName{$gene}	= "." ;
	push @{$vcf{$cols[0]}}, join("\t", @cols[1 .. $#cols ]) ;
}
close FP1 ;

# Obtain the exon coordinates from the gtf file
my %exons ;
open FP1, $ARGV[1] ;
my $prevTname = "-1" ;
my $strand = "." ;
my $gname = "-1" ;
my @range ;
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

			@{$exons{$gname}} = @range ;
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
		$strand = $cols[6] ;
		undef @range ;
	}
	
	if ( $hasChrPrefix == 1 && !( $cols[0] =~ /^c/) )
	{
		$cols[0] = "chr".$cols[0] ;
	}
	elsif ( $hasChrPrefix == 0 && $cols[0] =~ /^c/ )
	{
		$cols[0] = substr( $cols[0], 3 ) ;
	}
	push @range, $cols[0], $cols[3] - 1, $cols[4] - 1 ;			
}
close FP1 ;
 
foreach my $allele (keys %vcf)
{
	my @alleleVcf = @{$vcf{$allele}} ;
	my $gname = (split /\*/, $allele)[0] ;
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
	my $offset = $start ;
	my $firstOffset = $start ;
	# Change the sequence according to the vcf.
	my $tmp = $seq ;
	foreach my $v (@alleleVcf)
	{
		my @cols = split /\t/, $v ;
		my $pos = $cols[1] - 1 - $offset ;
		next if ($pos >= length($seq)) ;
		if ($cols[3] ne '.' && $cols[4] ne '.')
		{
			substr($seq, $pos, length($cols[3]), $cols[4]) ;
			$offset += (length($cols[3]) - length($cols[4])) ;
		}
		elsif ($cols[3] eq '.' && $cols[4] ne '.') # insertion
		{
			substr($seq, $pos, 0, $cols[4]) ;
			$offset -= length($cols[4]) ;
		}
		elsif ($cols[3] ne '.' && $cols[4] eq '.')
		{
			substr($seq, $pos, length($cols[3]), "") ;
			$offset += length($cols[3]) ; 
		}
	}
	# Move the exon coordinate according to the vcf file.
	for ($i = 0 ; $i < scalar(@alleleExon) ; $i += 3)
	{
		$alleleExon[$i + 1] -= $firstOffset ;
		$alleleExon[$i + 2] -= $firstOffset ;
	}

	foreach my $v (@alleleVcf)
	{
		my @cols = split /\t/, $v ;
		my $pos = $cols[1] - 1 ;
		my $shift = 0;
		if ($cols[3] ne '.' && $cols[4] ne '.')
		{
			$shift = (length($cols[3]) - length($cols[4])) ;
		}
		elsif ($cols[3] eq '.' && $cols[4] ne '.')
		{
			$shift = length($cols[4]) ;
		}
		elsif ($cols[3] ne '.' && $cols[4] eq '.')
		{
			$shift = -length($cols[3]) ; 
		}
		else
		{
			next ;
		}
		
		# exon coordinate should be reduced by 1 now.
		for ($i = 0 ; $i < scalar(@alleleExon) ; $i += 3)
		{
			$alleleExon[$i + 1] += $shift if ($alleleExon[$i + 1] >= $pos); 
			$alleleExon[$i + 2] += $shift if ($alleleExon[$i + 2] >= $pos);
		}
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
	print("ID   $allele\n")	;
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
