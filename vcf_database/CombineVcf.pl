#!/bin/perl

use strict ;
use warnings ;

use File::Basename ;

die "usage: a.pl default_allele_id vcf_file_list > output\n" if (@ARGV == 0) ;

open FPlist, $ARGV[1] ;
my $chr = "." ;

while (<FPlist>)
{
	chomp ;
	my $fname = $_ ;
	open FP, $fname ;

	while (<FP>)
	{
		next if (/^#/) ;
		chomp ;
		my $line = $_ ;
		my @cols = split ;
		$chr = $cols[0] ;
		
		$fname =~ s/.vcf// ;
		$fname =~ s/_/\*/ ;
		$fname = basename($fname) ;
		print(join("\t", ($fname, @cols[0..6])), "\n") ;
	}
	close FP ;
}

close FPlist ;

print(join("\t",($ARGV[0], $chr, 0, ".", ".", ".", ".", ".")), "\n")
