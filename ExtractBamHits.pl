#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: a.pl bwa.sam candidate_read.fq > bwa_aligned_candidate_read.fq\n" if (@ARGV == 0) ;

open FP, $ARGV[0] ;
my %readid ;
while (<FP>)
{
	next if (/^@/) ;
	chomp ;
	my @cols = split ;
	next if ($cols[2] eq "*") ;
	$readid{$cols[0]} = 1 ;
}
close FP ;

open FP, $ARGV[1] ;
while (<FP>) 
{
	my $header = $_ ;
	my $seq = <FP> ;
	my $separator = "" ;
	my $qual = "" ;

	if ($header =~ /^@/)
	{
		$separator = <FP> ;
		$qual = <FP> ;
	}
	chomp $header ;
	my @cols = substr($header, 1) ;
	if (defined $readid{$cols[0]})
	{
		print("$header\n$seq$separator$qual") ;
	}
}
close FP ;
