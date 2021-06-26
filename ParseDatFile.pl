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
				# UTR before
				my $start = $exons[0] - 50 ;
				my $end = $exons[0] - 1 ;
				$start = 0 if ($start < 0) ;
				$outputSeq .= substr($seq, $start, $end - $start + 1) ;

				for (my $i = 0 ; $i < scalar(@exons) ; $i += 2)
				{
					$outputSeq .= substr($seq, $exons[$i], $exons[$i + 1] - $exons[$i] + 1) ;
				}
				# UTR after
				$start = $exons[scalar(@exons) - 1] + 1;
				$end = $start + 49 ;
				$end = length($seq) - 1 if ($end >= length($seq)) ;
				$outputSeq .= substr($seq, $start, $end - $start + 1) ;
				
				$outputSeq = uc($outputSeq) ;				
				last if (defined $usedSeq{$outputSeq}) ;

				$usedSeq{$outputSeq} = 1 ;
				if (!defined $partialAlleles{$allele})
				{
					print(">$allele\n$outputSeq\n") ;
				}

				last ;
			}
		}
	}
}
close FP ;
