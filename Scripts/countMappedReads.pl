#!/usr/bin/perl -w
use miRWoods;
use Bio::DB::Sam;
use strict;
$| = 1;

my $usage = "USAGE:\n$0 <bamListFile>\n";

my $bamListFile = $ARGV[0] or die $usage;

my $librarySizes = "readRegions_librarySizes.txt\n";

my $bamList = miRWoods::loadBamList($bamListFile);

open(LSZ,">$librarySizes") or die "failed to open $librarySizes for writing\n";
foreach my $bamData (@{$bamList}) {    
    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamData};
    my $mappedReads = `samtools view -F 0x904 $bamFile | awk '{a[\$1]=1} END {print length(a)}'`;
    chomp($mappedReads);
    print LSZ "$sample\t$mappedReads\n";
}
close(LSZ);


###################################
#  Functions not currently used ###
###################################

sub countReadIds {
    my($chromLengths,$bamData) = @_;
    my $mappedCount = 0;
    my $alignedCount = 0;
    foreach my $chrom (keys %{$chromLengths}) {
	my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamData};
	my($tid,$chromStart,$chromStop) = $bamHeader->parse_region("$chrom"); 
	my $callBack = sub {
	    my($alignment,$data) = @_;
	    my($mappedCount,$alignedCount) = @{$data};
	    my $id = $alignment->qname;
	    my $hitCount = $alignment->get_tag_values('NH');
	    my $count;
	    if ($id =~ /.*_x(\d+)$/) {
		($count) = $id =~ /.*_x(\d+)$/;
	    } else {
		$count = 1;
	    }
	    my $adjustedCount = $count / $hitCount;
	    $$mappedCount += $adjustedCount;
	    $$alignedCount += $count;
	};
	my $callBackData = [\$mappedCount,\$alignedCount];
	my $code = $bamIndex->fetch($bam,$tid,$chromStart,$chromStop,$callBack,$callBackData);
    }
    return($mappedCount,$alignedCount);
}
