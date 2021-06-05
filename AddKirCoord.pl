#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: a.pl ipdkir_seq.fa gencode.gtf > gene_coord.fa\n" if (@ARGV == 0) ;

my %geneCoord ;
my $hasChrPrefix = 1 ;
my $defaultChr = "chr19" ;
$defaultChr = "19" if ($hasChrPrefix == 0) ;

open FP, $ARGV[0] ;
while (<FP>) 
{
	chomp ;
	next if (!/^>/) ;
	
	my $gene = (split /\*/, substr($_, 1))[0] ;
	$geneCoord{$gene} = "$defaultChr -1 -1 +" ;
}
close FP ;

open FP, $ARGV[1] ;
while (<FP>)
{
	next if ( /^#/ ) ;
	chomp ;
	my @cols = split /\t/ ;
	next if ( $cols[2] ne "gene" ) ;

	my $gname ;	
	if ( $cols[8] =~ /gene_name \"(.*?)\"/ )
	{
		#print $1, "\n" ; 
		$gname = $1 ;
	}
	else
	{
		die "No gene_name", $_, "\n" ;
	}

	
	if ( $hasChrPrefix == 1 && !( $cols[0] =~ /^c/) )
	{
		$cols[0] = "chr".$cols[0] ;
	}
	elsif ( $hasChrPrefix == 0 && $cols[0] =~ /^c/ )
	{
		$cols[0] = substr( $cols[0], 3 ) ;
	}
	
	$geneCoord{$gname} = join(" ", ($cols[0], $cols[3], $cols[4], $cols[6])) ;
}
close FP ;

open FP, $ARGV[0] ;
my $seq = "" ;
while (<FP>)
{
	chomp ;
	if (!/^>/)
	{
		$seq .= $_ ;
	}	
	else 
	{
		print $seq, "\n" if ($seq ne "") ;
		my $header = $_ ;
		my $gene = (split /\*/, substr($header, 1))[0] ;
		print "$header ".$geneCoord{$gene}."\n" ;
		$seq = "" ;	
	}
}
close FP ;
print $seq, "\n" if ($seq ne "") ;



