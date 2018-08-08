#!/usr/bin/perl -w
use strict;
$| = 1;

my $USAGE = "\nUSAGE:\n$0 <sam file> <output sam file>\n";

my $samFile = $ARGV[0] or die $USAGE;
my $outputFile = $ARGV[1] or die $USAGE;
#my($fileBase) = $samFile =~ /([^\/]*)\.sam/;

open(SAM, $samFile) or die "failed to open $samFile\n";
open(SAMOUTPUT, ">$outputFile") or die "failed to open $outputFile for writing\n";

my @savedEntries;
my $prevId;
my $totalReads = 0;
while (<SAM>) {
    chomp;
    if ( /^@/ ) {
	print SAMOUTPUT "$_\n";
    } else {
	my $entry = $_;
	my ($id) = split(/\t/,$entry);
	$prevId = defined($prevId) ? $prevId : $id;  #if prevId is not initialized, initialize it
	if ($prevId ne $id) {
	    my $readCount;
	    if ($prevId =~ /.*_x(\d+)$/) {
		($readCount) = $prevId =~ /.*_x(\d+)$/;
	    } else {
		$readCount = 1;
	    }
	    my $hitCount = scalar @savedEntries;
	    $totalReads += $readCount;
	    while(@savedEntries) {
		my $savedEntry = shift(@savedEntries);
		print SAMOUTPUT "$savedEntry\tNH:i:$hitCount\n";
	    }
	    $prevId = $id;
	}
	push(@savedEntries, $entry);
    }
}
my $hitCount = scalar @savedEntries;
while(@savedEntries) {
    my $savedEntry = shift(@savedEntries);
    print SAMOUTPUT "$savedEntry\tNH:i:$hitCount\n";
}
my $countMinLocus = $totalReads / 1000000;
print "Total Reads = $totalReads\n";
print "Suggested setting for countMinLocus (-c) option in processReadRegions.pl and evaluateReadRegions.pl: $countMinLocus\n";
print "If more than one sam file is being used add countMinLocus values toghether and use their sum for the -c option.\n";

close(SAM);
close(SAMOUTPUT);

exit;
