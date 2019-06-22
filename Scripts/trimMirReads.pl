#!/usr/bin/perl -w
use strict;
use miRTRAP1_6;
use Math::CDF;
$| = 1;

my $USAGE = "USAGE:\n$0 <samFile> <output SamFile>\n";
my $MAXPVAL = 0.27;

my $samFile = $ARGV[0] or die $USAGE;
my $outputSamFile = $ARGV[1] or die $USAGE;

createTrimmedSamFile($samFile,$outputSamFile,$MAXPVAL);


sub createTrimmedSamFile {
    my($samFile,$outputSamFile,$maxPValue) = @_;
    open(SAM,"$samFile") or die "failed to open $samFile for reading\n";
    open(OUTSAM, ">$outputSamFile") or die "failed to open $outputSamFile for writing\n";
    while (<SAM>) {
	chomp;
	my $samLine = $_;
	if ( /^@/ ) {
	    print OUTSAM $samLine . "\n";
	} else {
	    my($id,$flags,$chrom,$pos,$mapQ,$extendedCigar,$MRNM,$mPos,$iSite,$seq,$qualityScores,@tags) = split("\t", $samLine);
	    if (isMapped($flags)) {
		my $cigar = miRTRAP1_6::getTagValue(\@tags,'MD','Z') or die "no cigar tag found";
		my $cigarArray = breakCigar($cigar);
		my $strand = getStrand($flags);
		my ($trimmedSeq,$trimmedQualityScores);
		if ($strand eq '+') {
		   ($trimmedSeq,$trimmedQualityScores) = trimSequence($seq,$qualityScores,$cigarArray,$maxPValue);
		} elsif ($strand eq '-') {
		    my @reversedCigarArray = reverse(@{$cigarArray});
		    my $revSeq = reverse($seq);
		    my $revQualityScores = reverse($qualityScores);
		    ($trimmedSeq,$trimmedQualityScores) = trimSequence($revSeq,$revQualityScores,\@reversedCigarArray,$maxPValue);
		    $trimmedSeq = reverse($trimmedSeq);
		    $trimmedQualityScores = reverse($trimmedQualityScores);
		    my $trimAmount = length($seq) - length($trimmedSeq);
		    $pos = $pos + $trimAmount;  #adjusting position since 3' end of negative strand was trimmed
		} else {
		    die "strand not identified.  Flags = $flags\n";
		}
#		print "sequence is trimmed:\nseqId = $id\nseq = $seq\ntrimmedSeq = $trimmedSeq\nstrand = $strand\ncigar=$cigarTag\n\n";
		my $newExtCigar = getNewExtendedCigarValue($trimmedSeq);
		my $newSamLine = makeSamLine($id,$flags,$chrom,$pos,$mapQ,$newExtCigar,$MRNM,$mPos,$iSite,$trimmedSeq,$trimmedQualityScores,\@tags);
		print OUTSAM $newSamLine . "\n";
	    } else {
		print OUTSAM $samLine . "\n";
	    }
	}
    }
}

sub isMapped {
    my($flags) = @_;
    return ($flags & 4) ? 0 : 1;
}

sub getStrand   {
    my($flags) = @_;
    return ($flags & 16) ? '-' : '+';
}

sub makeSamLine {
    my($id,$flags,$chrom,$pos,$mapQ,$CIGAR,$MRNM,$mPos,$iSite,$seq,$qualityScores,$tags) = @_;
    my $tagField = join("\t", @{$tags});
    my $samLine = join("\t",$id,$flags,$chrom,$pos,$mapQ,$CIGAR,$MRNM,$mPos,$iSite,$seq,$qualityScores,$tagField);
    return $samLine;
}

####################################################
#  Cigar Functions                                 #
####################################################

sub breakCigar {
    my($cigar) = @_;
    my @cigarArray;

    while (length($cigar) != 0) {
	if ($cigar =~ /^\d+/) {
	    my($matchLength) = $cigar =~ /^(\d+)/;
	    $cigar = substr($cigar,length($matchLength),length($cigar)-length($matchLength));
	    push(@cigarArray,$matchLength);
	}
	if ($cigar =~ /^\D/) {
	    my($mismatch) = $cigar =~ /^(\D)/;
	    $cigar = substr($cigar,length($mismatch),length($cigar)-length($mismatch));
	    push(@cigarArray,$mismatch);
	    die "gaps were unexpected in bam file"  if ($mismatch eq "^");
	}
    }

    return(\@cigarArray);
}

sub getMismatchesFromCigar {
    my($cigarArray) = @_;
    my @mismatches;
    my $pos = 0;
    for (my $i = 0; $i < @{$cigarArray}; $i++) {
	if ($cigarArray->[$i] =~ /^\d+$/) {
	    $pos += $cigarArray->[$i];
	}
	if ($cigarArray->[$i] =~ /^\D$/) {  #there is a mismatch
	    push(@mismatches,$pos);
	    $pos++;
	}
    }
    return \@mismatches;
}

sub getNewExtendedCigarValue {
    my($seq) = @_;
    my $extendedCigar = length($seq) . "M";
    return $extendedCigar;
}


####################################################
#   Functions involved in sequence trimming        #
####################################################

sub trimSequence {
    my($seq,$qualityScores,$cigarArray,$maxPValue) = @_;
    my $mismatches = getMismatchesFromCigar($cigarArray);
    my $matchesAtEachMismatch = getNumMatchesAtEachMismatch($mismatches,length($seq));
    my $pValues = getPValueAtEachMisMatch($matchesAtEachMismatch, length($seq));

    my $trimPosition = length($seq);
    for (my $i = 0; $i < @{$pValues}; $i++) {  #pvalues should be sorted from last mismatch to first
	my ($mismatchPos,$pValue) = @{$pValues->[$i]};
	$trimPosition = $mismatchPos if ($pValue > $maxPValue);
#	print "seq=$seq\nmismatchPos=$mismatchPos\ntrimPosition = $trimPosition\tpVal=$pValue\n";
    }
    my $trimmedSeq = substr($seq,0,$trimPosition);

    my $trimmedQualityScores = substr($qualityScores,0,$trimPosition);
    return ($trimmedSeq,$trimmedQualityScores);
} 

#returns an array of the number of matches from end of sequence at each mismatch location
sub getNumMatchesAtEachMismatch { 
    my ($mismatches,$sequenceLength) = @_;
    my @revSortedMismatches = sort {$b <=> $a} @{$mismatches};
    my @matches;

    my $totalMismatches = 0;
    foreach my $mismatchPosition (@revSortedMismatches) {
	$totalMismatches++;
	my $numMatches = $sequenceLength - $mismatchPosition - $totalMismatches;

	push(@matches, [$mismatchPosition,$numMatches]);
    }

    return \@matches;
}

sub getPValue {
    #gets the probability that the matches within $length from the end were due to chance.
    my($matches,$length,$p) = @_;
    my $probMatchesDueToChance = ($matches) ? 1 - Math::CDF::pbinom($matches-1,$length,$p) : 1;
    return $probMatchesDueToChance;
}

sub getPValueAtEachMisMatch {
    my($mismatches,$seqLength) = @_;
    my $gcProb = 0.41;
    my $atProb = 0.59;
    my $p = 2*($gcProb/2)**2 + 2*($atProb/2)**2;
    my @pValues;
    foreach my $mismatchData (@{$mismatches}) {
	my($mismatchPosition,$numMatches) = @{$mismatchData};
	my $lengthFromEnd = $seqLength - $mismatchPosition;
	my $pValue = getPValue($numMatches,$lengthFromEnd,$p);
	push(@pValues,[$mismatchPosition,$pValue]);
    }

    my @sortedPValues = sort {$b->[0] <=> $a->[0]} @pValues; 

    return \@sortedPValues;
}

sub testPrint {
    my($mismatch,$pValues,$trimmedSeq) = @_;
    print "position\tmatches\tP(M >= m)\n";
    foreach my $mismatchData (@{$mismatch}) {
	my ($pos,$matches) = @{$mismatchData};
	print "$pos\t$matches\n";
    }
    print "\n";
    foreach my $pValueData (@{$pValues}) {
	my ($pos, $pValue) = @{$pValueData};
	print "$pos\t$pValue\n";
    }
    print "\n$trimmedSeq\n";
}

