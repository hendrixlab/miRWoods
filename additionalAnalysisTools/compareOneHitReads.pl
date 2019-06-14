#!/usr/bin/perl
use strict;

my $usage = "USAGE:\n <readCount file list> <max counts>\n";

my $rcFileList = $ARGV[0] or die $usage;
my $maxCounts = $ARGV[1] or die $usage;

my $outputFile = $rcFileList . "_leq" . $maxCounts . ".txt";

creatLEQMaxCountMatrix($rcFileList,$maxCounts,$outputFile);

sub creatLEQMaxCountMatrix {
    my($rcFileList,$maxCounts,$outputFile) = @_;
    open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
    open(RCFL,$rcFileList) or die "failed to open $rcFileList\n";
    while (<RCFL>) {
	chomp;
	unless ( /^#/ ) {
	    my($sample,$miRWoodsCountFile,$mirdeepCountFile,$miReapCountFile) = split(/\s+/);
	    print OPTF $sample
		."\t".countReadsLEQMax($miRWoodsCountFile,$maxCounts) 
		."\t".countReadsLEQMax($mirdeepCountFile,$maxCounts)
		."\t".countReadsLEQMax($miReapCountFile,$maxCounts)."\n";
	}
    }
    close(RCFL);
    close(OPTF);
}

sub countReadsLEQMax {
    my($readsFile,$maxCounts) = @_;
    my $numLEQMax = 0;
    open(RFLE,$readsFile) or die "failed to open $readsFile\n";
    while (<RFLE>) {
	chomp;
	unless ( /^#/ ) {
	    my($tag,$annotationType,$product,@readCounts) = split(/\s+/);
	    my $total = 0;
	    foreach my $count (@readCounts) {
		$total += $count;
	    }
	    if ($total <= $maxCounts) {
		$numLEQMax += 1;
	    }
	    if ($total == 0) {
		print "Warning: total is equal to 0 for $product\n";
	    }
	}
    }
    close(RFLE);
    return $numLEQMax;
}
