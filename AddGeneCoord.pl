#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: a.pl ipdkir_seq.fa gencode.gtf [OPTIONS] > gene_coord.fa\n".
	"\tOPTIONS:\n".
	"\t--gtf-gene-name-mapping STRING: a comma-separated string for how to translate a gene name in GTF to IPD gene name. Each string is like \"HFE:HLA-HFE\". (default: HFE:HLA-HFE)\n" 
if (@ARGV == 0) ;

my $i ;
my %geneCoord ;
my $hasChrPrefix = 1 ; # whether the output include chr prefix
my $defaultChr = "chr19" ;
$defaultChr = "19" if ($hasChrPrefix == 0) ;
my $geneNameMappingString = "HFE:HLA-HFE" ;
my %geneNameMapping ;
for ($i = 2 ; $i < scalar(@ARGV) ; ++$i)
{
	if ($ARGV[$i] eq "--gtf-gene-name-mapping")
	{
		$geneNameMappingString = $ARGV[$i + 1] ;
		++$i ;
	}
	else
	{
		die "Unknown options ".$ARGV[$i]."\n" ;
	}
}

my @cols = split /,/,$geneNameMappingString ;
for my $c (@cols)
{
	my @subCols = split /:/,$c ;
	$geneNameMapping{$subCols[0]} = $subCols[1] ;
}

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
		if (defined $geneNameMapping{$gname})
		{
			$gname = $geneNameMapping{$gname} ;
		}
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

	if (defined $geneCoord{$gname} && (split /\s/, $geneCoord{$gname})[1] == -1)
	{
		$geneCoord{$gname} = join(" ", ($cols[0], $cols[3], $cols[4], $cols[6])) ;
	}
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
		my $header = (split /\s+/, $_)[0] ;
		my $gene = (split /\*/, substr($header, 1))[0] ;
		print "$header ".$geneCoord{$gene}."\n" ;
		$seq = "" ;	
	}
}
close FP ;
print $seq, "\n" if ($seq ne "") ;



