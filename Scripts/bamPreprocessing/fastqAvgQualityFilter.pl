#!/usr/bin/perl -w
use strict;
$| = 1;

my $USAGE = "USAGE:\n$0 <fastq> <newFastq> <minQuality>\n";

my $fastq = $ARGV[0] or die $USAGE;
my $newFastq = $ARGV[1] or die $USAGE;
my $minQuality = $ARGV[2] or die $USAGE;

open(FQ,$fastq) or die "failed to open $fastq for reading\n";
open(NFQ,">$newFastq") or die "failed to open $newFastq for writing\n";
my @fastqEntry;
my $sizeGood = 0;
my $line = 0;
while (<FQ>) {
    chomp;
    if ( $line == 0 ) {
	@fastqEntry = ($_);
    } elsif ($line == 3) {
	my $qualLine = $_;
	if (checkQuality($qualLine,$minQuality)) {
	    print NFQ $fastqEntry[0] . "\n" . $fastqEntry[1] . "\n" . $fastqEntry[2] . "\n" . $qualLine . "\n";
	}
    } else {
	$fastqEntry[$line] = $_;
    }
    $line = ($line + 1) % 4;
}
close(FQ);
close(NFQ);

sub checkQuality {
    my($qualLine,$minQuality) = @_;
    my @qualScores = split('',$qualLine);
    my $qualTotal = 0;
    foreach my $qual (@qualScores) {
	$qualTotal += ord("$qual") - ord("!");
    }
    my $qualAvg = $qualTotal / (@qualScores);
    if ($qualAvg >= $minQuality) {
	return 1;
    } else {
	return 0;
    }
}
