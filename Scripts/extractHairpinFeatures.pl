#!/usr/bin/perl -w
use miRWoods;
use Bio::DB::Sam;
use Getopt::Long;
use strict;
$| = 1;

my $USAGE = "USAGE:\n$0 -L <config file>\n";

die $USAGE unless (@ARGV);

my $parameters = miRWoods::loadDefaultParameters();

Getopt::Long::Configure("no_ignore_case");

GetOptions ("lengthMin=i" => \$parameters->{lengthMin},
	    "totalLength=i" => \$parameters->{totalLength},
	    "hairpinShortLength=i" => \$parameters->{hairpinShortLength},
	    "distanceMin=i" => \$parameters->{distanceMin},
	    "reverseMax=i" => \$parameters->{reverseMax},
	    "countMinLocus=i" => \$parameters->{countMinLocus},
	    "fivePrimeHetMax=f" => \$parameters->{fivePrimeHetMax},
	    "shiftMin=i" => \$parameters->{shiftMin},
	    "OverlapMin=i" => \$parameters->{OverlapMin},
	    "InHairpinBuffer=i" => \$parameters->{InHairpinBuffer},
	    "OutHairpinBuffer=i" => \$parameters->{OutHairpinBuffer},
	    "RangeOfHairpin=i" => \$parameters->{RangeOfHairpin},
	    "outputPrefix=s" => \$parameters->{outputPrefix},
	    "bamListFile=s" => \$parameters->{bamListFile},
	    "FileRepeatRegions=s" => \$parameters->{FileRepeatRegions},
	    "genomeDir=s" => \$parameters->{genomeDir},
	    "SizesFile=s" => \$parameters->{SizesFile},
	    "LoadFromConfigFile=s" => \$parameters->{LoadFromConfigFile}
    );

unless ($parameters->{LoadFromConfigFile} eq "") {
    my $configFile = $parameters->{LoadFromConfigFile};
    miRWoods::readConfigFile($configFile,$parameters);
}
miRWoods::createOutputFileParameters($parameters);

my $chromSizesFile = $parameters->{SizesFile} or die "error: SizesFile parameter not set\n";
my $genomeDir = $parameters->{genomeDir} or die "error: genomeDir parameter not set\n";
my $readRegionsFile = $parameters->{allReadRegionsFile} or die print "allReadRegionsFile: parameter not loaded.\n";
my $hairpinsFile = $parameters->{hairpinsFile} or die "error: hairpinsFile parameter not set\n";
my $hairpinVectorFile = $parameters->{hairpinVectorFile} or die "error: hairpinVectorFile not set\n";
my $geneModelsFile = $parameters->{geneModels};
die "$chromSizesFile not found\n" unless (-e $chromSizesFile);
die "$hairpinsFile not found\n" unless (-e $hairpinsFile);
my $chromLengths = miRWoods::readChromLengths($chromSizesFile) or die "failed to load chrom.sizes file\n";
my($posClass,$negClass) = (1,-1);
my $hairpinList = readHairpinListFile($hairpinsFile);
my $readRegions = readReadRegions($readRegionsFile);
my $geneModels = ($geneModelsFile) ? readGeneModelFile($geneModelsFile) : 0;
my $introns = ($geneModels) ? getIntronsFromGeneModels($geneModels) : {};
my $neighborCounts = getNeighborCounts($hairpinList,$readRegions,$introns,$parameters);
#printing vector file
miRWoods::printHairpinVectorFile($hairpinsFile,$hairpinVectorFile,$genomeDir,$chromLengths,$negClass,$neighborCounts,$parameters);

##########################################
# Functions to Compute neighbor Counts   #
##########################################

#Note: We will NOT compute nonMiRNeighborCounts like was done in the original miRTRAP
#      because that will require having a set of predicted miR loci.  Instead we will
#      compute just the neighborCount, and let the random forest determine what influence
#      it will have given each hairpins other scores.  This feature helps reduce false
#      positives and yet the random forest still appears able to find miRs within clusters.
#      The decision boundary may be influenced by other scores in these cases.

sub getNeighborCounts {
    my($hairpinList,$readRegions,$introns,$parameters) = @_;
    my %neighborCounts;
    my $neighborWindow = $parameters->{neighborWindow} or die "neighborWindow: parameter not loaded.\n";
    foreach my $chrom (keys %{$hairpinList}) {
	for(my $i1=0;$i1<@{$hairpinList->{$chrom}};$i1++) {
	    my($tag1,$start1,$stop1,$strand1,$hStart1,$hStop1,$totalSense1,$totalAntiSense1,$mfe1,$seq1,$fold1) = @{$hairpinList->{$chrom}[$i1]};
	    my($leftBound,$rightBound) = getIntronBounds($introns,$chrom,$hairpinList->{$chrom}[$i1],$parameters);
	    my $negativeCount = 0;
	    if($readRegions->{$chrom}) {
		for(my $i2=0;$i2<@{$readRegions->{$chrom}};$i2++) {
		    my($tag2,$start2,$stop2,$strand2) = @{$readRegions->{$chrom}[$i2]};
		    if(($leftBound <= $start2)&&($stop2 <= $rightBound)) {
			$negativeCount++;
		    }
		}
	    }
	    my $adjNeighborCount = int($negativeCount * ($stop1-$start1+1+2*$neighborWindow) / ($rightBound-$leftBound+1) + 0.4999999999);
	    $neighborCounts{$tag1} = $adjNeighborCount;
	}
    }
    return(\%neighborCounts);
}

sub extractNegativeReadRegions {
    my($readRegions) = @_;
    my %negativeReadRegions;
    my $buffer = 50;
    foreach my $chrom (keys %{$readRegions}) {
	foreach my $region (@{$readRegions->{$chrom}}) {
	    my($rId,$rStart,$rStop,$rStrand) = @{$region};
	    push(@{$negativeReadRegions{$chrom}},$region);
	}
    }
    return \%negativeReadRegions;
}

sub getIntronBounds {
    my($introns,$chrom,$hairpin,$parameters) = @_;
    my $neighborWindow = $parameters->{neighborWindow} or die "neighborWindow: parameter not loaded.\n";
    my($tag,$start,$stop,$strand,$hStart,$hStop) = @{$hairpin};    
    my @intronList;
    my $thisStart = $start + 50;
    my $thisStop = $stop - 50;
    my $leftBound = $thisStart - $neighborWindow;
    my $rightBound = $thisStop + $neighborWindow;
    $leftBound = $leftBound > 0 ? $leftBound : 0;
    #print "starting with $leftBound $rightBound\n";
#    print "$chrom\t$strand\n";
    if($introns->{$chrom}{$strand}) {	
	foreach my $intron (@{$introns->{$chrom}{$strand}}) {
	    my($inStart,$inStop,$geneId,$type1,$type2) = @{$intron};
	    if(($inStart<=$thisStart+5)&&($stop-5<=$inStop)) {
		push(@intronList,[$inStart,$inStop,$geneId]);
	    }
	}
    }
    if(@intronList) {
	if(@intronList == 1) {
	    my $intron = shift(@intronList);
	    my($inStart,$inStop) = @{$intron};
	    #print "found: $inStart $inStop\n";
	    if($thisStart - $inStart <= $neighborWindow) {
		$leftBound = $inStart;
	    }
	    if($inStop - $stop <= $neighborWindow) {
		$rightBound = $inStop;
	    }
	    #print "loaded: $leftBound $rightBound\n";
	    return ($leftBound,$rightBound);
	} else {
	    my $minLength = $rightBound - $leftBound;
	    foreach my $intron (@intronList) {
		my($inStart,$inStop,$geneId) = @{$intron};
		if($inStop - $inStart < $minLength) {
		    $minLength = $inStop - $inStart;
		    $leftBound = $inStart;
		    $rightBound = $inStop;
		}		
	    }
	    return ($leftBound,$rightBound);
	}
    } else {
	return ($leftBound,$rightBound);
    }
}

sub getIntronsFromGeneModels {
    my($geneModels) = @_;
    my %introns;
    my %USED;
    foreach my $geneId (keys %{$geneModels}) {
	my @exons;
	foreach my $term (@{$geneModels->{$geneId}}) {
	    my($type,$chrom,$start,$stop,$strand) = @{$term};
	    if($type eq "exon") {
		push(@exons,$term);
	    }
	}
	my @sorted = sort {$a->[2] <=> $b->[2]} @exons;
	for(my $i=0;$i<@sorted-1;$i++) {
	    my($type1,$chrom1,$start1,$stop1,$strand1) = @{$exons[$i]};
	    my($type2,$chrom2,$start2,$stop2,$strand2) = @{$exons[$i+1]};
	    if($stop1+1 < $start2) {
		unless($USED{$chrom1.$stop1.$start2.$strand1}) {
		    push(@{$introns{$chrom1}{$strand1}},[$stop1+1,$start2-1,$geneId,$type1,$type2]);
		    $USED{$chrom1.$stop1.$start2.$strand1}++;
		}
	    }
	}
    }
    return \%introns;
}

sub readGeneModelFile {
    my($geneModelGffFile) = @_;
    my %geneModels;
    my %mRNA;
    open(FILEA,$geneModelGffFile) or die "could not open $geneModelGffFile in readGeneModelFile()\n";
    while(<FILEA>) {
        unless(/\#/) {
	    chomp;
            my($chrom,$source,$type,$start,$stop,$score,$strand,$phase,$info)=split(/\t/);
	    if($type eq "mRNA") {
		if($info =~ /ID=(.*?);/) {
		    my($transcriptId) = $info =~ /ID=(.*?);/;
		    $mRNA{$transcriptId}++;
		} elsif($info =~ /ID=(.*)/) {
		    my($transcriptId) = $info =~ /ID=(.*)/;
		    $mRNA{$transcriptId}++;
		}
	    }
	}
    }
    open(FILEA,$geneModelGffFile) or die "could not open $geneModelGffFile in readGeneModelFile()\n";
    while(<FILEA>) {
        unless(/\#/) {
            chomp;
            my($chrom,$source,$type,$start,$stop,$score,$strand,$phase,$info)=split(/\t/);
	    if($type eq "exon") {
                my($parent) = $info =~ /Parent=(.*);?/;
		my @parents = split(/\,/,$parent);		
		foreach my $geneId (@parents) {		    
		    if($mRNA{$geneId}) {
			push(@{$geneModels{$geneId}},[$type,$chrom,$start,$stop,$strand]);
		    }
		}
	    }
	}
    }
    return \%geneModels;
}

sub readHairpinListFile {
    my $hairpinFile = shift;
    my %hairpinList;
    open(MFILE,$hairpinFile) or die "could not open $hairpinFile in readHairpinListFile()\n";
    while(<MFILE>) {
        chomp;
        unless(/\#/) {
	    my($tag,$chrom,$start,$stop,$strand,$leftCenter,$rightCenter,$totalSense,$totalAntiSense,$mfe,$seq,$fold) = split(/\t/);
	    push(@{$hairpinList{$chrom}},[$tag,$start,$stop,$strand,$leftCenter,$rightCenter,$totalSense,$totalAntiSense,$mfe,$seq,$fold,]);
	}
    }
    return \%hairpinList;
}

sub readReadRegions {
    my $readRegionsFile = shift;
    my %readRegions;
    open(RRF,$readRegionsFile) or die "could not open $readRegionsFile in in readReadRegions()\n";
    while(<RRF>) {
	my($chrom,$start,$stop,$id,$rrMaxCount,$strand,$rejectionReason) = split();
	push(@{$readRegions{$chrom}},[$id,$start,$stop,$strand]);
    }    
    return \%readRegions;
}
