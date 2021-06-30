#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: a.pl xxx.dat [-f xxx_gene.fa --mode rna|dna] > yyy.fa\n" if (@ARGV == 0) ;

my %selectedAlleles ;
my $selectAlleleFile = "" ;
my $mode = "rna" ;

my $i ;
for ($i = 1 ; $i < scalar(@ARGV) ; ++$i)
{
	if ($ARGV[$i] eq "-f")
	{
		$selectAlleleFile = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ($ARGV[$i] eq "--mode")
	{
		$mode = $ARGV[$i + 1] ;
		++$i ;
	}
	else
	{
		die "Unknown option ".$ARGV[$i]."\n" ;
	}
}
	
if ($selectAlleleFile ne "" )
{
	open FP, $selectAlleleFile ;
	while (<FP>)
	{
		if (/^>/)
		{
			chomp ;
			#my $allele = (split /\s+/, $_)[1] ;#substr($_, 1) ;
			my $allele = substr($_, 1) ;
			my $seq = <FP> ;
			chomp $seq ;
			$selectedAlleles{$allele} = $seq ;
		}
	}
	close FP ;
}

open FP, $ARGV[0] ;
my @exons ;
my $seq = "" ;
my $allele = "" ;
my %partialAlleles ;
my $hasIntron ;
my %usedSeq ;
my @alleleOrder ;
my %alleleSeq ;
my %gene5UTRPadding ;
my %geneBestPossible5UTRPadding ;
my %gene3UTRPadding ;
my %geneBestPossible3UTRPadding ;
my %allelePaddingLength ;
my %alleleEffectiveLength ; # exon length + 2 * utrLength 

my $utrLength = 50 ;
my $intronPaddingLength = 200 ;

while (<FP>)
{
	if (/^FT/)
	{
		if (/allele=\"(.*?)\"/)
		{
			$allele = $1 ;
		}
		elsif (/\sCDS\s/)
		{
			undef @exons ;
			$hasIntron = 0 ;
			$seq = "" ;
		}
		elsif (/\sexon\s/)
		{
			chomp ;
			my @cols = split /\s+/, $_ ;
			my ($start, $end) = ($cols[2] =~ /(\d+)\.\.(\d+)/) ;
			push @exons, ($start - 1 ,$end - 1) ;
		}
		elsif (/pseudo$/)
		{
			pop @exons ; # psuedo exon
			pop @exons ;
		}
		elsif (/\sintron\s/)
		{
			$hasIntron = 1 ;
		}
		elsif (/partial$/)
		{
			$partialAlleles{$allele} = 1 ;
		}
	}
	elsif (/^SQ/)
	{
		# skip the header
		while (<FP>)
		{
			if (!/^\/\//)
			{
				chomp ;
				my @cols = split /\s+/, $_ ;
				foreach my $s (@cols[0..$#cols - 1])
				{
					$seq .= $s ;
				}
			}
			else
			{
				my $outputSeq = "" ;
				last if ($mode eq "dna" && $hasIntron == 0) ;

				# UTR before
				my $start = $exons[0] - $utrLength ;
				my $end = $exons[0] - 1 ;
				my $gene = (split /\*/, $allele)[0] ;
				@{$allelePaddingLength{$allele}} = (0, 0) ;
				if ($start < 0) 
				{
					$allelePaddingLength{$allele}[0] = -$start ;
					if (!defined $geneBestPossible5UTRPadding{$gene} ||
						$end > length($geneBestPossible5UTRPadding{$gene}))
					{
						$geneBestPossible5UTRPadding{$gene} = uc(substr($seq, 0, $end)) ;
					}
					$start = 0 ;
				}
				elsif (!defined $gene5UTRPadding{$gene})
				{
					$gene5UTRPadding{$gene} = uc(substr($seq, $start, $end - $start + 1)) ;
				}

				$outputSeq .= substr($seq, $start, $end - $start + 1) ;

				if ($mode eq "rna")
				{
					for (my $i = 0 ; $i < scalar(@exons) ; $i += 2)
					{
						$outputSeq .= substr($seq, $exons[$i], $exons[$i + 1] - $exons[$i] + 1) ;
					}
				}
				elsif ($mode eq "dna")
				{
					for (my $i = 0 ; $i < scalar(@exons) ; $i += 2)
					{
						$start = $exons[$i] ;
						$end = $exons[$i + 1] ;

						if ($i > 0)
						{
							$start = $exons[$i] - $intronPaddingLength ;
							$start = 0 if ($start < 0) ;
						}

						while ($i + 2 < scalar(@exons))
						{
							$end = $exons[$i + 1] + $intronPaddingLength ;
							$end = length($seq) - 1 if ($end >= length($seq)) ;
							if ($end >= $exons[$i + 2] - $intronPaddingLength)
							{
								# short intron.
								$i += 2 ;
							}
							else
							{
								last ;
							}
						}
						$outputSeq .= substr($seq, $start, $end - $start + 1) ;
					}
				}
				else
				{
					die "Unknown mode $mode\n" ;
				}
				# UTR after
				$start = $exons[scalar(@exons) - 1] + 1;
				$end = $start + $utrLength - 1 ;
				if ($end >= length($seq)) 
				{
					$allelePaddingLength{$allele}[1] = $end - length($seq) + 1 ;
					if (!defined $geneBestPossible3UTRPadding{$gene} ||
						length($seq) - $start > length($geneBestPossible3UTRPadding{$gene}))
					{
						$geneBestPossible3UTRPadding{$gene} = uc(substr($seq, $start)) ;
					}
					$end = length($seq) - 1 ;
				}
				elsif (!defined $gene3UTRPadding{$gene})
				{
					$gene3UTRPadding{$gene} = uc(substr($seq, $start, $end - $start + 1)) ;
				}
				$outputSeq .= substr($seq, $start, $end - $start + 1) ;
				
				$outputSeq = uc($outputSeq) ;				
				
				if (!defined $partialAlleles{$allele})
				{
					push @alleleOrder, $allele ;
					$alleleSeq{$allele} = $outputSeq ;
					#print(">$allele\n$outputSeq\n") ;
				}

				$alleleEffectiveLength{$allele} = 2 * $utrLength ;
				for (my $i = 0 ; $i < scalar(@exons) ; $i += 2)
				{
					$alleleEffectiveLength{$allele} += $exons[$i + 1] - $exons[$i] + 1 ;
				}
				last ;
			}
		}
	}
}
close FP ;

srand(17) ;
my @numToNuc = ("A", "C", "G", "T") ;
for my $allele (@alleleOrder)
{
	my $gene = (split /\*/, $allele)[0] ;
	if (!defined $gene5UTRPadding{$gene})
	{
		my $randomSeq = "" ;
		for ($i = 0 ; $i < $utrLength ; ++$i)
		{
			$randomSeq .= $numToNuc[int(rand(4))] ;
		}
		my $len = length($geneBestPossible5UTRPadding{$gene}) ;
		substr($randomSeq, -$len, $len, $geneBestPossible5UTRPadding{$gene}) ;
		$gene5UTRPadding{$gene} = $randomSeq ;
	}
	if (!defined $gene3UTRPadding{$gene})
	{
		my $randomSeq = "" ;
		for ($i = 0 ; $i < $utrLength ; ++$i)
		{
			$randomSeq .= $numToNuc[int(rand(4))] ;
		}
		my $len = length($geneBestPossible3UTRPadding{$gene}) ;
		substr($randomSeq, 0, $len, $geneBestPossible3UTRPadding{$gene}) ;
		$gene3UTRPadding{$gene} = $randomSeq ;
	}
}
for my $allele (@alleleOrder)
{
	my $outputSeq = $alleleSeq{$allele} ;
	my $gene = (split /\*/, $allele)[0] ;
	if ($allelePaddingLength{$allele}[0] > 0)
	{
		$outputSeq = substr($gene5UTRPadding{$gene}, 0, $allelePaddingLength{$allele}[0]).$outputSeq ;
	}
	if ($allelePaddingLength{$allele}[1] > 0)
	{
		$outputSeq = $outputSeq.substr($gene3UTRPadding{$gene}, -$allelePaddingLength{$allele}[1]) ;
	}

	next if (defined $usedSeq{$outputSeq}) ;
	$usedSeq{$outputSeq} = 1 ;
	print(">$allele/".$alleleEffectiveLength{$allele}."\n$outputSeq\n") ;
}
