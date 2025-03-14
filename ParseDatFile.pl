#!/usr/bin/env perl

use strict ;
use warnings ;

#--partialInRnaMode INT: include explicitly annotated partial alleles if its length no less than the mode length by <int> for that gene in RNA mode. [0] 
die "usage: a.pl xxx.dat [-f xxx_gene.fa --mode rna|dna|genome --gene KIR|HLA|... --partialInRnaMode INT --partialIntronHasNoSeq --ignoreParital --intronPadding INT(200) --dedup] > yyy.fa\n" if (@ARGV == 0) ;

sub FindMode
{
	my %dist = %{$_[0]} ;
	my $ret = -1 ;
	my $max = -1 ;
	foreach my $k (keys %dist) 
	{
		if ($dist{$k} > $max) 
		{
			$max = $dist{$k} ;
			$ret = $k ;
		}
		elsif ($dist{$k} == $max && $k ge $ret)
		{
			$ret = $k ;
		}
	}
	return $ret ;
}

sub GetLastExonLength
{
	my @exons = @{$_[0]} ;
	return $exons[scalar(@exons) - 1] - $exons[scalar(@exons) - 2] + 1 ;
}

my %selectedAlleles ;
my $selectAlleleFile = "" ;
my $mode = "rna" ;
my $genePrefix = "" ;
my $fixGeneLength = 0 ;
my $ignorePartial = 0 ;
my $includePartialDiffLen = 0 ;
my $partialIntronHasNoSeq = 0 ; # the partial introns have no sequence in the .dat file. This is an issue from IPD-KIR 2.13.0.
my $intronPaddingLength = 200 ;
my $dedup = 0 ;

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
	elsif ($ARGV[$i] eq "--gene")
	{
		$genePrefix = uc($ARGV[$i + 1]) ;
		++$i ;
	}
  elsif ($ARGV[$i] eq "--ignorePartial")
  {
    $ignorePartial = 1 ;
    ++$i ;
  }
  elsif ($ARGV[$i] eq "--partialInRnaMode")
  {
    $includePartialDiffLen = $ARGV[$i + 1] ;
    ++$i ;
  }
  elsif ($ARGV[$i] eq "--partialIntronHasNoSeq")
  {
    $partialIntronHasNoSeq = 1 ;
  }
	elsif ($ARGV[$i] eq "--intronPadding")
	{
		$intronPaddingLength = $ARGV[$i + 1] ;
		++$i ;
	}
  elsif ($ARGV[$i] eq "--dedup")
  {
    $dedup = 1 ;
  }
	else
	{
		die "Unknown option ".$ARGV[$i]."\n" ;
	}
}

if ($mode eq "rna")
{
	$fixGeneLength = 1 ;
}
elsif ($mode eq "dna")
{
	$fixGeneLength = 1 ;
}
$includePartialDiffLen = -1 if ($mode eq "genome") ;

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
my @introns ;
my $seq = "" ;
my $allele = "" ;
my %partialAlleles ;
my $hasIntron ;
my $localIntronLen ;
my $descriptionState ; # 0: exon, 1: intron
my $partialIntronLen ;
my $isPartial ;
my $pseudoExonLen ;
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
my %alleleExonRegions ; # the actuall exon regions in the output sequence
my %alleleTrueExonRegions ; # the true exon regions in the underlying genome
my %geneLastExonLengthDist ;

if ($mode eq "genome")
{
	$utrLength = 0 ;
}

while (<FP>)
{
	if (/^ID/)
	{
			undef @exons ;
			$hasIntron = 0 ;
      $partialIntronLen = 0 ;
			$isPartial = 0 ;
			$seq = "" ;
			$allele = "-1" ;
      $pseudoExonLen = 0 ;
	}
	elsif (/^FT/)
	{
		if (/allele=\"(.*?)\"/)
		{
			$allele = $1 ; 
		}
		#elsif (/\sCDS\s/)
		#{
		#	undef @exons ;
		#	$hasIntron = 0 ;
		#	$seq = "" ;
		#}
		elsif (/\sexon\s/)
		{
			chomp ;
			my @cols = split /\s+/, $_ ;
			my ($start, $end) = ($cols[2] =~ /(\d+)\.\.(\d+)/) ;
			push @exons, ($start - 1 - $partialIntronLen, $end - 1 - $partialIntronLen) ;
      $descriptionState = 0 ;
      $pseudoExonLen = 0 ;
		}
		elsif (/pseudo$/)
		{
			my $end = pop @exons ; # psuedo exon
			my $start = pop @exons ;
      $pseudoExonLen = $end - $start + 1 ;
		}
		elsif (/\sintron\s/)
		{
      if ($partialIntronHasNoSeq == 1)
      {
        chomp ;
        my @cols = split /\s+/, $_ ;
        my ($start, $end) = ($cols[2] =~ /(\d+)\.\.(\d+)/) ;
        $localIntronLen = $end - $start + 1 ;
      }

			++$hasIntron ;
      $descriptionState = 1 ;
		}
		elsif (/partial$/)
		{
			#$partialAlleles{$allele} = 1 ;
      if ($descriptionState == 0 
        || $partialIntronHasNoSeq == 0)
      {
			  $isPartial = 1 ;
      }
      else
      {
        $partialIntronLen += $localIntronLen ;
        --$hasIntron ;
      }

      if ($pseudoExonLen > 0 && $partialIntronHasNoSeq == 1) # Assume partial comes after pseudo
      {
        $partialIntronLen += $pseudoExonLen ;
      }
		}
	}
	elsif (/^SQ/)
	{
		$partialAlleles{$allele} = 1 if ($isPartial) ;
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
				last if ($mode eq "genome" && $hasIntron == 0 && scalar(@exons) > 2) ;
				last if ($allele eq "-1") ;
				last if (scalar(@exons) == 0) ;
				# UTR before
				my $start = $exons[0] - $utrLength ;
				my $end = $exons[0] - 1 ;
				my $gene = (split /\*/, $allele)[0] ;
				@{$allelePaddingLength{$allele}} = (0, 0) ;
				my $exonOffset = 0 ;
				my @exonActualRegion ;
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

				$exonOffset = $utrLength ;
				if ($mode eq "rna")
				{
					for (my $i = 0 ; $i < scalar(@exons) ; $i += 2)
					{
						$outputSeq .= substr($seq, $exons[$i], $exons[$i + 1] - $exons[$i] + 1) ;
						push @exonActualRegion, $exonOffset ;
						push @exonActualRegion, $exonOffset + $exons[$i + 1] - $exons[$i] ;
						$exonOffset += ($exons[$i + 1] - $exons[$i] + 1) ;
					}
				}
				elsif ($mode eq "dna")
				{
					for (my $i = 2 ; $i < scalar(@exons) ; $i += 2)
					{
						if ($exons[$i] <= $exons[$i - 1] + 1)
						{
							$partialAlleles{$allele} = 1;						
						}
					}
					for (my $i = 0 ; $i < scalar(@exons) ; $i += 2)
					{
						$start = $exons[$i] ;
						$end = $exons[$i + 1] ;

						if ($i > 0)
						{
							$start = $exons[$i] - $intronPaddingLength ;
							$start = 0 if ($start < 0) ; # no need to worry about, if this happens, exons will merge
							$exonOffset += 1 + $intronPaddingLength ; # +1 here is for the 'N' separator
							$outputSeq .= 'N' ;
						}

						push @exonActualRegion, $exonOffset ; 
						push @exonActualRegion, $exonOffset + $exons[$i + 1] - $exons[$i] ; 

						my $k = $i ;
						while ($i + 2 < scalar(@exons))
						{
							$end = $exons[$i + 1] + $intronPaddingLength ;
							$end = length($seq) - 1 if ($end >= length($seq)) ;

							if ($end >= $exons[$i + 2] - $intronPaddingLength)
							{
								# short intron.
								$i += 2 ;
								$end = $exons[$i + 1] ;
								push @exonActualRegion, $exonOffset + $exons[$i] - $exons[$k] ; 
								push @exonActualRegion, $exonOffset + $exons[$i + 1] - $exons[$k] ; 
							}
							else
							{
								last ;
							}
						}

						$outputSeq .= substr($seq, $start, $end - $start + 1) ; 
						$exonOffset += ($exons[$i + 1] - $exons[$k] + 1) ;
						$exonOffset += $intronPaddingLength ;
					}
					@{$alleleTrueExonRegions{$allele}} = @exons ;
				}
				elsif ($mode eq "genome") 
				{
					for (my $i = 2 ; $i < scalar(@exons) ; $i += 2)
					{
						if ($exons[$i] <= $exons[$i - 1] + 1)
						{
							$partialAlleles{$allele} = 1;						
						}
					}
					$outputSeq = $seq ;
					@exonActualRegion = @exons ;
				}
				else
				{
					die "Unknown mode $mode\n" ;
				}
				++${$geneLastExonLengthDist{$gene}}{GetLastExonLength(\@exons)} ;

				# UTR after
				$start = $exons[scalar(@exons) - 1] + 1;
				if ($start > length($seq)) # partial allele that the exon ends after the sequence ends
				{
					$partialAlleles{$allele} = 1;						
				}
				else
				{
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
				}
				
				$outputSeq = uc($outputSeq) ;				
				
				if (!defined $partialAlleles{$allele})
				{
					push @alleleOrder, $allele ;
					#print(">$allele\n$outputSeq\n") ;
				}

        $alleleSeq{$allele} = $outputSeq ;
        @{$alleleExonRegions{$allele}} = @exonActualRegion ;
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

# Get statstiscs for exon count, exon length, and true (with respect to the underlying true genome)intron length 
my %geneEffectiveSeqLengthDist ;
my %geneLengthMode ;
my %geneExonCnt ;
my %geneExonCntMode ;
my %geneExonLengthDist ;
my %geneExonLengthMode ;
my %geneTrueIntronLengthDist ;
my %geneTrueIntronLengthMode ;
if ($mode eq "dna")
{
	foreach my $allele (@alleleOrder)
	{
		my $gene = (split /\*/, $allele)[0] ;
		++${$geneEffectiveSeqLengthDist{$gene}}{$alleleEffectiveLength{$allele}} ;

		my $exonCnt = scalar(@{$alleleExonRegions{$allele}}) / 2 ;
		++${$geneExonCnt{$gene}}{$exonCnt} ;
	}

	foreach my $gene (keys %geneEffectiveSeqLengthDist)
	{
		$geneLengthMode{$gene} = FindMode(\%{$geneEffectiveSeqLengthDist{$gene}}) ; 
	}

	foreach my $gene (keys %geneExonCnt)
	{
		$geneExonCntMode{$gene} = FindMode(\%{$geneExonCnt{$gene}}) ; 
	}

	# Get the representative exon length
	foreach my $allele (@alleleOrder)
	{
		my $gene = (split /\*/, $allele)[0] ;
		++${$geneEffectiveSeqLengthDist{$gene}}{$alleleEffectiveLength{$allele}} ;

		my @exons = @{$alleleExonRegions{$allele}} ;
		my @trueExons = @{$alleleTrueExonRegions{$allele}} ;
		my $exonCnt = scalar(@exons) / 2 ;
		next if ($exonCnt != $geneExonCntMode{$gene}) ;
		for (my $i = 0 ; $i < $exonCnt ; ++$i) 
		{
			++${${$geneExonLengthDist{$gene}}{$i}}{$exons[2*$i+1] - $exons[2*$i] + 1} ;
			if ($i < $exonCnt - 1)
			{
				++${${$geneTrueIntronLengthDist{$gene}}{$i}}{$trueExons[2*$i+2] - $trueExons[2*$i+1] - 1} ;
			}
		}
	}
	foreach my $gene (keys %geneExonLengthDist)
	{
		for my $i (keys %{$geneExonLengthDist{$gene}})
		{
			${$geneExonLengthMode{$gene}}{$i} = FindMode(\%{${$geneExonLengthDist{$gene}}{$i}}) ;
		}
		for my $i (keys %{$geneTrueIntronLengthDist{$gene}})
		{
			${$geneTrueIntronLengthMode{$gene}}{$i} = FindMode(\%{${$geneTrueIntronLengthDist{$gene}}{$i}}) ;
		}
	}
}

# Rescue partial alleles
if ($includePartialDiffLen >= 0 && $ignorePartial == 0)
{
	if (scalar(keys %geneLengthMode) == 0)
	{
		foreach my $allele (@alleleOrder)
		{
			my $gene = (split /\*/, $allele)[0] ;
			++${$geneEffectiveSeqLengthDist{$gene}}{$alleleEffectiveLength{$allele}} ;
		}

		foreach my $gene (keys %geneEffectiveSeqLengthDist)
		{
			$geneLengthMode{$gene} = FindMode(\%{$geneEffectiveSeqLengthDist{$gene}}) ; 
		}
	}

  # rna-mode rescue is easy, we just make sure the length is about right. 
  my @rescuedAlleles ;
  if ($mode eq "rna")
  {
    foreach my $allele (keys %partialAlleles)
    {
      my $gene = (split /\*/, $allele)[0] ;
      my $len = $alleleEffectiveLength{$allele} ;
      next if (!defined $geneLengthMode{$gene}) ;
      if ($len >= $geneLengthMode{$gene} - $includePartialDiffLen)
      {
        push @rescuedAlleles, $allele ;
      }
    }
  }
  elsif ($mode eq "dna")
  {
    # DNA mode we need to get the intronic anchor sequences 
    # The introns should be from the alleles with common exon count
        
    # Collect the consensus of introns
    my %geneIntronSeq ;
    my %geneIntronSeqMode ;
    foreach my $allele (@alleleOrder)
    {
      my $gene = (split /\*/, $allele)[0] ;
      my @exons = @{$alleleExonRegions{$allele}} ;
      my $exonCnt = scalar(@exons) / 2 ;
      next if ($exonCnt != $geneExonCntMode{$gene}) ;
      for ($i = 2 ; $i < 2 * $exonCnt ; $i += 2)
      {
        my $intronSeq = substr($alleleSeq{$allele}, $exons[$i - 1] + 1, 
          $exons[$i] - $exons[$i - 1] - 1) ;
        ++${${$geneIntronSeq{$gene}}{$i / 2 - 1}}{$intronSeq} ;
      }
    }
    
    foreach my $gene (keys %geneIntronSeq)
    {
      foreach $i (keys %{$geneIntronSeq{$gene}})
      {
        ${$geneIntronSeqMode{$gene}}{$i} = FindMode(\%{${$geneIntronSeq{$gene}}{$i}}) ;
        #printf("$gene $i ".${$geneIntronSeqMode{$gene}}{$i}."\n" ) ;
      }
    }
    
    # Adding introns to the partial alleles
    foreach my $allele (keys %partialAlleles)
    {
      my $gene = (split /\*/, $allele)[0] ;
      my $len = $alleleEffectiveLength{$allele} ;
      next if (!defined $geneLengthMode{$gene}) ;
      next if ($len < $geneLengthMode{$gene} - $includePartialDiffLen) ;
      my @exons = @{$alleleExonRegions{$allele}} ;
      my $exonCnt = scalar(@exons) / 2 ;
      next if ($exonCnt != $geneExonCntMode{$gene}) ;
    
      my $exonOffset = 0 ;
      my $outputSeq = $alleleSeq{$allele} ;
      
      # The exons holds the final exon offset, including the paddings 
      #   so we need to subtract those out first to get the 
      #   coordinate on real sequence at this step.
      my $extra5UtrLength = $allelePaddingLength{$allele}[0] ;
      for ($i = 0 ; $i < 2 * $exonCnt ; ++$i)
      {
        $exons[$i] -= $extra5UtrLength ;
      }

      for ($i = 2 ; $i < 2 * $exonCnt ; $i += 2)
      {
        if ($exons[$i] + $exonOffset == $exons[$i - 1] + 1)
        {
          my $intronSeq =  ${$geneIntronSeqMode{$gene}}{$i/2 - 1} ;
          #$outputSeq = substr($outputSeq, 0, $exons[$i - 1] + 1).
          #            $intronSeq.
          #            substr($outputSeq, $exons[$i]);
          substr($outputSeq, $exons[$i - 1] + 1, 0, $intronSeq) ;
          $exonOffset += length($intronSeq) ;
        }
        $exons[$i] += $exonOffset ;
        $exons[$i + 1] += $exonOffset ;
      }
      
      for ($i = 0 ; $i < 2 * $exonCnt ; ++$i)
      {
        $exons[$i] += $extra5UtrLength ;
      }
      @{$alleleExonRegions{$allele}} = @exons ;

      $alleleSeq{$allele} = $outputSeq ;
      push @rescuedAlleles, $allele ;
    }
  }

  push @alleleOrder, @rescuedAlleles ;
}

srand(17) ;
my @numToNuc = ("A", "C", "G", "T") ;
foreach my $allele (@alleleOrder)
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

foreach my $allele (@alleleOrder)
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
	$alleleSeq{$allele} = $outputSeq;
}

# Fix the exonization that may affect the anchor intron sequence
# Should be put after UTR padding so the exon regions matches
if ($mode eq "dna")
{
	foreach my $allele (@alleleOrder)
	{
		my $gene = (split /\*/, $allele)[0] ;
		my @exons = @{$alleleExonRegions{$allele}} ;
		my $exonCnt = scalar(@exons) / 2 ;
		next if ($exonCnt != $geneExonCntMode{$gene}) ;
		next if (!defined $alleleTrueExonRegions{$allele}) ;

		my $updated = 0 ;
		for (my $i = 0 ; $i < $exonCnt - 1 ; ++$i) # last exon has some special treatment
		{
			my $exonLength = $exons[2*$i+1] - $exons[2*$i] + 1 ;
			if (0)#$allele eq "HLA-C*07:02:01:17N")
			{
				print "$i ", $exons[2*$i]," ", $exons[2*$i + 1],
					" $exonLength ", ${$geneExonLengthMode{$gene}}{$i}, " ",
					$exons[2*$i + 1] + 1 + $intronPaddingLength, " ",
					substr($alleleSeq{$allele}, $exons[2*$i + 1] + 1 + $intronPaddingLength, 1), "\n" ;
				if ($i == 0)
				{
					print $alleleSeq{$allele}, "\n" ;
				}
			}
			# Handle the issue of exonization 
			if ($exonLength > ${$geneExonLengthMode{$gene}}{$i})
			{
				my $trim = $exonLength - ${$geneExonLengthMode{$gene}}{$i} ;
				my $trimSide = 0 ;
				my @trueExons = @{$alleleTrueExonRegions{$allele}} ;
				my $posN = 0 ;
				my $newSeq = "" ;
				if ($trueExons[2 * $i + 2] - $trueExons[2 * $i + 1] - 1 + $trim == ${$geneTrueIntronLengthMode{$gene}}{$i} 
					&& $exons[2*$i + 1] + 1 + $intronPaddingLength < length($alleleSeq{$allele})
					&& substr($alleleSeq{$allele}, $exons[2*$i + 1] + 1 + $intronPaddingLength, 1) eq "N") 
				{
					#print("Right exonization! $allele\n") ;
					$trimSide = 1 ; 
					$posN = $exons[2*$i + 1] + 1 + $intronPaddingLength ;
					$newSeq = substr($alleleSeq{$allele}, 0, $posN - $trim)
							.substr($alleleSeq{$allele}, $posN) ; #'N' is included in this portion
				}
				elsif ($i > 0 
					&& $trueExons[2 * $i] - $trueExons[2 * $i - 1] - 1 + $trim == ${$geneTrueIntronLengthMode{$gene}}{$i - 1} 
					&& $exons[2*$i] - 1 - $intronPaddingLength >= 0
					&& substr($alleleSeq{$allele}, $exons[2*$i - 1] - 1 - $intronPaddingLength, 1) eq "N") 
				{
					#print("Left exonization! $allele\n") ;
					$trimSide = -1 ;
					$posN = $exons[2*$i] - 1 - $intronPaddingLength ;
					$newSeq = substr($alleleSeq{$allele}, 0, $posN + 1) # 'N' is included in this portion
								.substr($alleleSeq{$allele}, $posN + $trim + 1) ;
				}
				
				next if ($trimSide == 0) ;

				$alleleSeq{$allele} = $newSeq ;
				
				#$alleleSeq{$allele} = substr($alleleSeq{$allele}, $posN-$trim, $trim, "") ;
				if ($trim > $intronPaddingLength) # the exon itself should be shrinked/trimmed
				{
					$exons[2*$i + 1] -= ($trim - $intronPaddingLength) if ($trimSide == 1) ;
					$exons[2*$i] += ($trim + $intronPaddingLength) if ($trimSide == -1) ;
				}

				if ($trimSide == -1)
				{
					$exons[2*$i] -= $trim ;
					$exons[2*$i+1] -= $trim ;
				}

				# Shift the reamining exons
				for (my $j = $i + 1 ; $j < $exonCnt ; ++$j)
				{
					$exons[2*$j] -= $trim ;
					$exons[2*$j+1] -= $trim ;
				}
				$updated = 1 ;
			}
		}
		@{$alleleExonRegions{$allele}} = @exons if ($updated == 1);
	}
}

my %geneSeqLengthDist ;
foreach my $allele (@alleleOrder)
{
	my $gene = (split /\*/, $allele)[0] ; 
	++${$geneSeqLengthDist{$gene}}{length($alleleSeq{$allele})} ;
}

my %geneSeqLength ;
my %geneLastExonLength ;
foreach my $gene (keys %geneSeqLengthDist)
{
	$geneSeqLength{$gene} = FindMode(\%{$geneSeqLengthDist{$gene}}) ;
	$geneLastExonLength{$gene} = FindMode(\%{$geneLastExonLengthDist{$gene}}) ;
}

if ($fixGeneLength == 1)
{
	# Trim the allele sequence if it is too long:
	#  1. Longer than other alleles
	#  2. Last exon is longer than other alleles
	foreach my $allele (@alleleOrder)	
	{
		my $outputSeq = $alleleSeq{$allele} ;
		my $gene = (split /\*/, $allele)[0] ;
		my $lastExonLength = GetLastExonLength(\@{$alleleExonRegions{$allele}}) ;
		my $trim = $lastExonLength - $geneLastExonLength{$gene} ;
		if (length($outputSeq) > $geneSeqLength{$gene} &&
			$trim > 0)
		{
			$outputSeq = substr($outputSeq, 0, length($outputSeq) - $trim) ;
		}
		$alleleSeq{$allele} = $outputSeq;
	}
}

foreach my $allele (@alleleOrder)
{
	my $outputSeq = $alleleSeq{$allele} ;
	next if ($outputSeq eq "") ;
	next if ($dedup == 1 && defined $usedSeq{$outputSeq}) ;
	next if (!($allele =~ /^$genePrefix/)) ;
	$usedSeq{$outputSeq} = 1 ;
	print(">$allele ".scalar(@{$alleleExonRegions{$allele}}) / 2 .
		" ".join(" ", @{$alleleExonRegions{$allele}}).
		"\n$outputSeq\n") ;
}
