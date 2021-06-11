#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: a.pl xxx.dat > yyy.fa\n" if (@ARGV == 0) ;

open FP, $ARGV[0] ;
my @selectedRegion ;
my $seq = "" ;
my $allele = "" ;

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
			chomp ;
			my $line = $_ ;
			my $catLine = "" ;
			$catLine = $line ;
			while (!($line =~ /\)/)) # the join might span multiple lines
			{
				$line = <FP> ;
				chomp $line ;
				my @cols = split /\s+/, $line ;
				$catLine .= $cols[1] ;
			}
			
			undef @selectedRegion ;
			if ($catLine =~ /join\((.*?)\)/)
			{
				my $tmp = $1 ;
				my @cols = split /,/, $tmp ;
				foreach my $c (@cols)
				{
					# note that the coorinate is 1-based in dat file.
					my ($start, $end) = ($c =~ /(\d+)\.\.(\d+)/) ;
					push @selectedRegion, ($start - 1 ,$end - 1) ;
				}
			}	
			else
			{
				die "Unknown format $catLine\n" ;
			}		
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
				my $start = $selectedRegion[0] - 50 ;
				my $end = $selectedRegion[0] - 1 ;
				$start = 0 if ($start < 0) ;
				$outputSeq .= substr($seq, $start, $end - $start + 1) ;

				for (my $i = 0 ; $i < scalar(@selectedRegion) ; $i += 2)
				{
					$outputSeq .= substr($seq, $selectedRegion[$i], $selectedRegion[$i + 1] - $selectedRegion[$i] + 1) ;
				}
				# UTR after
				$start = $selectedRegion[scalar(@selectedRegion) - 1] + 1;
				$end = $start + 49 ;
				$end = length($seq) - 1 if ($end >= length($seq)) ;
				$outputSeq .= substr($seq, $start, $end - $start + 1) ;
				
				$outputSeq = uc($outputSeq) ;				
				print(">$allele\n$outputSeq\n") ;

				last ;
			}
		}
	}
}
close FP ;
