#!/usr/bin/perl -w
use miRWoods;
use Bio::DB::Sam;
use Getopt::Long;
use File::Copy;
use strict;
$| = 1;

my $USAGE = "USAGE:\n$0 -L <config file>";

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

my $productFile = $parameters->{productFile} or die "error: productFile parameter not set\n";
my $chromSizesFile = $parameters->{SizesFile} or die "error: SizesFile parameter not set\n";
my $genomeDir = $parameters->{genomeDir} or die "error: genomeDir parameter not set\n";
my $hairpinVectorFile = $parameters->{hairpinVectorFile} or die "error: hairpinVectorFile parameter not set\n";
my $hairpinVecTrainFile = $parameters->{hairpinTrainFile} or die "error: hairpinTrainFile parameter not set\n";
my $bamListFile = $parameters->{bamListFile} or die "error: bamListFile parameter not set\n";
my $bamList = miRWoods::loadBamList($bamListFile) or die "error: failed to load $bamListFile\n";
my $mirBaseGff = $parameters->{mirbaseGff} or die "error: mirBaseGFF parameter not set\n";
#my $balancedHairpinVecTrainFile = $parameters->{hairpinTrainFileSB};
#my $prodBalancedHPVecTrainFile = "readRegions_hairpinTrainFilePB.txt";
die "$bamListFile not found" unless (-e $bamListFile);
die "$mirBaseGff not found" unless (-e $mirBaseGff);
die "$hairpinVectorFile not found" unless (-e $hairpinVectorFile);
die "$productFile not found" unless (-e $productFile);
die "$chromSizesFile not found" unless (-e $chromSizesFile);
my $chromLengths = miRWoods::readChromLengths($chromSizesFile) or die "failed to load chrom.sizes file\n";
my($posClass,$negClass) = (1,-1);
my $annotations = miRWoods::readMirbaseGff3($mirBaseGff);
miRWoods::printHairpinTrainingData($hairpinVectorFile,$mirBaseGff,$genomeDir,$chromLengths,$posClass,$hairpinVecTrainFile,$productFile,$parameters);
#printHairpinTrainingData_usingMirbaseFolds($bamList,$hairpinVecTrainFile,$genomeDir,$chromLengths,$mirBaseGff,$parameters);
#copy($hairpinVecTrainFile,$balancedHairpinVecTrainFile) or die "Failed to copy $hairpinVecTrainFile to $balancedHairpinVecTrainFile\n";
#copy($hairpinVecTrainFile,$prodBalancedHPVecTrainFile) or die "Failed to copy $hairpinVecTrainFile to $prodBalancedHPVecTrainFile";

if ($parameters->{negGff}) {
    my $curatedNegGff = $parameters->{negGff} or die "error: negGff parameter not set\n";
    die "$curatedNegGff not found\n" unless (-e $curatedNegGff);
    miRWoods::printHairpinTrainingData($hairpinVectorFile,$curatedNegGff,$genomeDir,$chromLengths,$negClass,$hairpinVecTrainFile,$parameters);
} else {
    #negative products have not been curated so we will assume that most features from the productVectorFile are not true miRs and add from that
#    miRWoods::addRandomNegHPFeatures($hairpinVecTrainFile,$hairpinVectorFile);
    printNegTrainData($hairpinVectorFile,$mirBaseGff,$genomeDir,$chromLengths,$posClass,$hairpinVecTrainFile,$productFile,$parameters);
}

#miRWoods::addRandomNegHPFeatures_SizeBalanced($balancedHairpinVecTrainFile,$hairpinVectorFile,$mirBaseGff,$parameters);
#miRWoods::addRandomNegHPFeatures_ProdBalanced($prodBalancedHPVecTrainFile,$hairpinVectorFile,$mirBaseGff,$parameters);


sub printNegTrainData {
    my($hairpinVectorFile,$mirBaseGff,$genomeDir,$chromLengths,$posClass,$hairpinVecTrainFile,$productFile,$parameters) = @_;
    my($annoHairpins,$annoProducts) = miRWoods::readMirbaseGff3($mirBaseGff);
    open(HTF,">>$hairpinVecTrainFile") or die "failed to open $hairpinVecTrainFile for appending\n";
#    print "reading in $hairpinVectorFile\n";
    open(HVF,$hairpinVectorFile) or die "failed to open $hairpinVectorFile\n";
    while (<HVF>) {
	chomp;
	unless ( /^GeneId/ ) {
	    my($geneId,$class,@info) = split(/\t/);
	    my($tag,$location) = miRWoods::parseGeneId($geneId);
	    unless (checkAnnotationLocations($annoHairpins,$location,0)) {
		print HTF join("\t", $geneId,'-1',@info) . "\n";
	    }
	}
    }
    close(HVF);
    close(HTF);
}

sub checkAnnotationLocations {
    my($annoHairpins,$location,$strandDependent) = @_;
    #strandDependent variable determines whether to only return items that are only on the correct strand or not
    my($chrom,$start,$stop,$strand) = miRWoods::parseLocation($location);
    foreach my $hairpinInfo (@{$annoHairpins->{$chrom}}) {
	my($annoStart,$annoStop,$annoStrand,$id,$name) = @{$hairpinInfo};
	unless ($strandDependent && ($strand ne $annoStrand)) {
	    if (miRWoods::getOverlap($start,$stop,$annoStart,$annoStop)) {
		return $name;
	    }
	}
    }
    return 0;
}

##################################################
# Old Experimental Functions: not currently used #
##################################################

sub addRandomNegHPFeatures_ProdBalanced {
    #adds random negative features from feature vector file to the training data
    my($hairpinVecTrainFile,$hairpinScoresFile,$gff,$parameters) = @_;
    my $productFile = $parameters->{productFile};
    my $negVectorMultiplier = 1;  #determines ratio of positive to negative scores 
    my %vectorBuckets;
    #my %posNames;
    my %positiveHairpinEntries;
    my $posHairpinCount = 0;
    #my $minOverlapPercent = 0.7;
    my $minOverlapPercent = 0;   #we want to avoid any hairpins that overlap mirs for our negative set
    my($mirHairpins,$mirProducts) = miRWoods::readMirbaseGff3($gff);
    my $hairpins = getHairpinsFromGeneIds($hairpinScoresFile,1);
    #getHairpinFoldRegions is a temorary fix needed to product the correct mirOverlaps.  See function for more info.
    my $hairpinFoldRegions = getHairpinFoldRegions($parameters->{hairpinsFile});
    my($mirOverlaps, $overlapInfo) = getMirOverlaps($mirHairpins,$hairpins,$minOverlapPercent,$hairpinFoldRegions);
    foreach my $mirName (keys %{$mirOverlaps}) {
	#print "mirName = $mirName\n";
	foreach my $hairpinName (@{$mirOverlaps->{$mirName}}) {
	    my($start,$stop,$strand,$geneId) = @{$overlapInfo->{$hairpinName}};
	    $positiveHairpinEntries{$geneId} = 1;
	    #print "geneID = $geneId\n";
	}
    }
    open(HVF, $hairpinVecTrainFile) or die "failed to open $hairpinVecTrainFile\n";
    while(<HVF>) {
	chomp;
	my($hpName) = split("\t", $_);
	#$posNames{$hpName} = 1;
	$posHairpinCount++;
    }
    close(HVF);
    #print "posHairpinCount = $posHairpinCount\n";
    my $prodCounts = getHPProdCounts($productFile);
    my($prodTargetCounts,$maxSenseProds,$maxAntisenseProds) = getHPProdCountTargets($hairpinVecTrainFile,$prodCounts,$negVectorMultiplier);
    open(HSF, $hairpinScoresFile) or die "failed to open $hairpinScoresFile\n";
    my $FIRST = 1;
    while(<HSF>) {
	chomp;
	my $vectorLine = $_;
	unless ($FIRST) {
	    my($geneId,$class,$count) = split("\t", $vectorLine);
	    unless ($positiveHairpinEntries{$geneId}) {
		my $numSenseProds = $prodCounts->{$geneId}{'sense'};
		my $numAntisenseProds = $prodCounts->{$geneId}{'antisense'};
		push(@{$vectorBuckets{$numSenseProds}{$numAntisenseProds}}, $vectorLine);
	    }
	} else {
	    $FIRST = 0;
	}
    }
    close(HSF);
    open(HVFA, ">>$hairpinVecTrainFile") or die "failed to open $hairpinVecTrainFile for appending\n";
    my $carryOverAmount = 0;  #when not enough vector data is present in a higher bucket we will use more in a smaller bucket
    for (my $i = $maxAntisenseProds; $i > 0; $i--) {
	for (my $j = $maxSenseProds; $j > 0; $j--) {
	    my $targetAmount = $prodTargetCounts->{$j,$i} + $carryOverAmount;
	    my @shuffledVectorData = List::Util::shuffle(@{$vectorBuckets{$j}{$i}});
	    if ($targetAmount > @shuffledVectorData) {
		$carryOverAmount = $targetAmount - @shuffledVectorData;
		$targetAmount = @shuffledVectorData;
	    } else {
		$carryOverAmount = 0;
	    }
	    for (my $i = 0; $i < $targetAmount; $i++) {
		print HVFA $shuffledVectorData[$i] . "\n";
	    }
	}
    }
    close(HVFA);
}

sub addRandomNegHPFeatures_SizeBalanced {
    #adds random negative features from feature vector file to the training data
    my($hairpinVecTrainFile,$hairpinScoresFile,$gff,$parameters) = @_;
    my @bucketRanges = ([0,1],[1,100],[100,500],[500,1000],[1000,5000],[5000,10000],[10000,~0]);
    my %vectorBuckets;
    #my %posNames;
    my %positiveHairpinEntries;
    my $posHairpinCount = 0;
    #my $minOverlapPercent = 0.7;
    my $minOverlapPercent = 0;   #we want to avoid any hairpins that overlap mirs for our negative set
    my($mirHairpins,$mirProducts) = miRWoods::readMirbaseGff3($gff);
    my $hairpins = getHairpinsFromGeneIds($hairpinScoresFile,1);
    #getHairpinFoldRegions is a temorary fix needed to product the correct mirOverlaps.  See function for more info.
    my $hairpinFoldRegions = getHairpinFoldRegions($parameters->{hairpinsFile});
    my($mirOverlaps, $overlapInfo) = getMirOverlaps($mirHairpins,$hairpins,$minOverlapPercent,$hairpinFoldRegions);
    foreach my $mirName (keys %{$mirOverlaps}) {
	#print "mirName = $mirName\n";
	foreach my $hairpinName (@{$mirOverlaps->{$mirName}}) {
	    my($start,$stop,$strand,$geneId) = @{$overlapInfo->{$hairpinName}};
	    $positiveHairpinEntries{$geneId} = 1;
	    #print "geneID = $geneId\n";
	}
    }
    open(HVF, $hairpinVecTrainFile) or die "failed to open $hairpinVecTrainFile\n";
    while(<HVF>) {
	chomp;
	my($hpName) = split("\t", $_);
	#$posNames{$hpName} = 1;
	$posHairpinCount++;
    }
    close(HVF);
    #print "posHairpinCount = $posHairpinCount\n";
    open(HSF, $hairpinScoresFile) or die "failed to open $hairpinScoresFile\n";
    my $FIRST = 1;
    while(<HSF>) {
	chomp;
	my $vectorLine = $_;
	unless ($FIRST) {
	    my($geneId,$class,$count) = split("\t", $vectorLine);
	    unless ($positiveHairpinEntries{$geneId}) {
		for (my $i = 0; $i < @bucketRanges; $i++) {
		    my($startRange,$stopRange) = @{$bucketRanges[$i]};
		    if ($count > $startRange && $count <= $stopRange) {
			push(@{$vectorBuckets{$i}}, $vectorLine);
		    }
		}
	    }
	} else {
	    $FIRST = 0;
	}
    }
    close(HSF);
    open(HVFA, ">>$hairpinVecTrainFile") or die "failed to open $hairpinVecTrainFile for appending\n";
    my $carryOverAmount = 0;  #when not enough vector data is present in a higher bucket we will use more in a smaller bucket
    foreach my $bucketNum (sort {$b <=> $a} keys %vectorBuckets) {
	my @shuffledVectorData = List::Util::shuffle(@{$vectorBuckets{$bucketNum}});
	my $numRetrieve = $posHairpinCount / @bucketRanges + $carryOverAmount;
	if ($numRetrieve > @shuffledVectorData) {
	    $carryOverAmount = $numRetrieve - @shuffledVectorData;
	    $numRetrieve = @shuffledVectorData;
	}
	for (my $i = 0; $i < $numRetrieve; $i++) {
	    print HVFA $shuffledVectorData[$i] . "\n";
	}
    }
    close(HVFA);
}

sub getMirOverlaps {
    my($mirHairpins,$hairpins,$minOverlapPercent,$hairpinFoldRegions) = @_;
    my %mirOverlaps;
    my %overlapInfo;
    #print "minOverlapPercent = $minOverlapPercent\n";
    foreach my $chrom (keys %{$hairpins}) {
	foreach my $mirInfo (@{$mirHairpins->{$chrom}}) {
	    my($mirStart,$mirStop,$mirStrand,$mirId,$mirName) = @{$mirInfo};
	    #print "\noverlaps of $mirName:\n";
	    foreach my $hpInfo (@{$hairpins->{$chrom}}) {
		my($start,$stop,$strand,$geneId,$name) = @{$hpInfo};
		if ($strand eq $mirStrand) {
		    my($foldStart,$foldStop) = @{$hairpinFoldRegions->{$name}};
		    my $minOverlap = ($foldStop - $foldStart + 1) * $minOverlapPercent;
		    my $overlap = miRWoods::getOverlap($foldStart,$foldStop,$mirStart,$mirStop);
		    #print "\t$name\t$overlap\toverlapRequired=$minOverlap\n" if ($overlap > 0);
		    if ($overlap > $minOverlap) {
			#print "\tadded to mirOverlaps\n";
			push(@{$mirOverlaps{$mirName}}, $name);
			$overlapInfo{$name} = $hpInfo unless($overlapInfo{$name});
			#print "$name pushed on with $mirName\n";
		    }
		}
	    }
	}
    }
    return(\%mirOverlaps,\%overlapInfo);
}

sub getHairpinsFromGeneIds {
    #geneIdFile may be any tab delimeted file with a geneId in the first column
    my($geneIdFile,$header) = @_;
    my %hairpinLocations;  #a hash of hairpins similar to the hash created for miR's in readMirbaseGff3
    my $first = ($header) ? 1 : 0;
    open(GIDF, $geneIdFile) or die "failed to open $geneIdFile for reading\n";
    while (<GIDF>) {
	chomp;
	unless($first) {
	    my($geneId) = split("\t",$_);
	    my($name,$location) = miRWoods::parseGeneId($geneId);
	    my($chrom,$start,$stop,$strand) = miRWoods::parseLocation($location);
	    push(@{$hairpinLocations{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
	} else {
	    $first = 0;
	}
    }
    close(GIDF);
    return \%hairpinLocations
}

sub getHairpinFoldRegions {
    #this function was designed to take care of a bug in which the whole hairpin region folded and unfolded was being checked
    #to see if it was overlapping a true mir in getMirOverlaps.  This was done to avoid calling hairpins true mirs if their
    #unvolded region overlapped the annotation.  In the future, processReadRegions could add the start and stop of the fold region
    #to the hairpins file to avoid needing this function and to create cleaner code.
    my($hairpinsFile) = @_;
    my %foldRegions;
    open(HP,$hairpinsFile) or die "Failed to open $hairpinsFile for reading\n";
    while (<HP>) {
	chomp;
	unless (/^#/) {
	    my($name,$chrom,$start,$stop,$strand,$leftCenter,$rightCenter,$totalSense,$totalAntisense,$mfe,$sequence,$fold) = split("\t", $_);
	    my($middle,$left) = $fold =~ /^((\+*).*[^\+])/;
	    my $relStart = length($left);
	    my $relStop = length($middle);
	    my($foldStart,$foldStop);
	    if ($strand eq '+') {
		$foldStart = $start + $relStart;
		$foldStop = $start + $relStop;
	    } else {
		$foldStart = $stop - $relStop;
		$foldStop = $stop - $relStart;
	    }
	    @{$foldRegions{$name}} = ($foldStart,$foldStop);	    
	}
    }
    close(HP);
    return \%foldRegions;
}

sub getHPProdCounts {
    my($hairpinVecTrainFile,$productFile) = @_;
    my %prodCounts;
    my %geneIds;
    open(HVF,$hairpinVecTrainFile) or die "failed to open $hairpinVecTrainFile for reading\n";
    while (<HVF>) {
	chomp;
	my($geneId) = split("\t", $_);
	my($hpName,$hpLocation) = miRWoods::parseGeneId($geneId);
	$geneIds{$hpName} = $geneId;
    }
    close(HVF);
    open(PF,$productFile) or die "failed to open $productFile for reading\n";
    while (<PF>) {
	chomp;
	unless ( /^#/ ) {
	    my($tag,$side,$type,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$start,$stop,$strand) = split("\t");
	    my $geneId = $geneIds{$tag};
	    unless ($prodCounts{$geneId}) {
		#initializing
		$prodCounts{$geneId}{'sense'} = 0;
		$prodCounts{$geneId}{'antisense'} = 0;
	    }
	    if ($type =~ /as/) {
		$prodCounts{$geneId}{'antisense'} += 1;
	    } else {
		$prodCounts{$geneId}{'sense'} += 1;
	    }
	}
    }
    close(PF);
    return \%prodCounts;
}

sub getHPProdCountTargets {
    #takes the hairpins from the training set and loads a hash with the product counts
    #these counts are used in addRandomNegHPFeatures_ProdBalanced() to create a balanced negative set
    my($hairpinVecTrainFile,$prodCounts,$negVectorMultiplier) = @_;
    my %prodTargetCounts;
    my($maxSenseProds,$maxAntisenseProds) = (-1,-1);
    open(HVF,$hairpinVecTrainFile) or die "failed to open $hairpinVecTrainFile for reading\n";
    while (<HVF>) {
	chomp;
	my($hpName) = split("\t", $_);
	my $numSenseProds = $prodCounts->{$hpName}{'sense'};
	my $numAntisenseProds = $prodCounts->{$hpName}{'antisense'};
	$prodTargetCounts{$numSenseProds}{$numAntisenseProds} += 1 * $negVectorMultiplier;
	$maxSenseProds = $numSenseProds if ($numSenseProds > $maxSenseProds);
	$maxAntisenseProds = $numAntisenseProds if ($numAntisenseProds > $maxAntisenseProds);
    }
    close(HVF);
    return \%prodTargetCounts;
}

sub printHairpinTrainingData_usingMirbaseFolds {
    my($bamList,$trainFile,$genomeDir,$chromLengths,$gff,$parameters) = @_;
    my($mirHairpins,$mirProducts) = miRWoods::readMirbaseGff3($gff);
    my $hairpinVectorFile = $parameters->{outputPrefix} . "_hairpinScores";
    my $librarySizesFile = $parameters->{librarySizes} or die "no librarySizes entry found in config file\n";
    my $librarySizes = miRWoods::loadLibrarySizes($librarySizesFile);
    my $predProductScores = readPredProductClassFile("readRegions_predProductClasses.txt");
    my $testReadCounts = 0;
    my $readsCount = 0; #number of mirBase hairpins with reads
    my $noReadsCount = 0;  #number of mirBase hairpins with no reads
    my $readsLessThanCountMinLocus = 0; #number of hairpins with fewer than count min locus reads
    my $centersReadsCount = 0; #number of centers on mirBase hairpins with reads
    my $noCenterCount = 0;  #number of hairpins without centers (hopefully this will stay 0)
    #my $mu = Memory::Usage->new() if ($testMemory);
    #my $tk = TimeKeeper->new() if ($testTime);
    my @sampleList;
    my $minLength = $parameters->{lengthMin} or die "minLength: not loaded.\n";
    my $filteredCandidatesFile = $parameters->{filteredCandidatesFile};
    my $hairpinsFile = $parameters->{hairpinsFile};
    my $distinctReadsFile = $parameters->{distinctReadsFile};
    my $featureVectorFile = $parameters->{featureVectorFile};
    my $maxPredictedProductScores = (-e $parameters->{predProductClasses}) ? miRWoods::getMaxPredictedProductScores($parameters->{predProductClasses}) : 0;
    miRWoods::initializeRNAFold();
#    open(HL,">".$hairpinsFile) or die "failed to open $hairpinsFile for writing\n";
#    open(HDR,">".$distinctReadsFile) or die "failed to open $distinctReadsFile for writing\n";
#    open(HPL,">".$productFile) or die "failed to open $productFile for writing\n";
    open(HTF,">".$trainFile) or die "failed to open $hairpinVectorFile for writing\n";
    print HTF "GeneId\ty\ttotalSense\ttotalAntisense\tmfe\taapd\ttapd\turf\tahc\tafh\tpbp\tsameShift\tbothShift\tSPA\tOPA\t";
    print HTF "aa\tac\tag\tat\tca\tcc\tcg\tct\tga\tgc\tgg\tgt\tta\ttc\ttg\ttt\tduplexEnergy\tGCcontent\t";
    #print HTF "zScore\t";
    print HTF "foldDupCmp\texactFoldDupCmp\tdupPBP\tdupLoopLength";
    print HTF "\tAPV\twAPV";
    print HTF "\tARV\twARV";
    print HTF "\tmpCount\tdupCount\tdupOverlap\tmpLoopDistance\tdupLoopDistance\tloopSize\tmpOverlapScore\tdupOverlapScore\twghtMPOverlapIntensity\twghtDupOverlapIntensity\tOPA2\ttotalOverlapAmount\taverageOverlapAmount\ttotalRelativeOverlapAmount\taverageRelativeOverlapAmount\tmaxOverlap\tinnerLoopGapCount\ttotalSenseRPM\ttotalAntisenseRPM";
    print HTF "\tRFProductAvg\n";
#    print HL "#tag\tchrom\tstart\tstop\tstrand\tleftCenter\trightCenter\ttotalSense\ttotalAntisense\tmfe\tsequence\tfold\t";
#    print HL "aapd\ttapd\turf\tahc\tmpfh\tpbp\tsameShifted\tbothShift\tsplitProductAbundance\toverlapProductAbundance\t";
#    print HL "aa\tac\tag\tat\tca\tcc\tcg\tct\tga\tgc\tgg\tgt\tta\ttc\ttg\ttt\tduplexEnergy\tGCcontent\tfoldDupCmp\texactFoldDupCmp\tdupPBP\tdupLoopLength";
##    print HL "aa\tac\tag\tat\tca\tcc\tcg\tct\tga\tgc\tgg\tgt\tta\ttc\ttg\ttt\tduplexEnergy\tGCcontent\tzScore\tfoldDupCmp\texactFoldDupCmp\tdupPBP\tdupLoopLength";
#    print HL "\tAPV\twAPV";
#    print HL "\tARV\twARV";
#    print HL "\tmpCount\tdupCount\tdupOverlap\tmpLoopDistance\tdupLoopDistance\tloopSize\tmpOverlapScore\tdupOverlapScore\twghtMPOverlapIntensity\twghtDupOverlapIntensity\tOPA2\ttotalOverlapAmount\taverageOverlapAmount\ttotalRelativeOverlapAmount\taverageRelativeOverlapAmount\tmaxOverlap\tinnerLoopGapCount\tsenseRPM\tantiSenseRPM";
#    if ($maxPredictedProductScores) {
#	print HL "\tRFProductAvg\n";
#    } else {
#	print HL "\n";
#    }
#    print HDR "#tag\tstart\tstop\tstrand\ttotal reads\toffsetStart\trelative Start\trelative Stop\tsequence";
#    print HPL "#tag\tside type\ttype\ttotal reads\ttotal most abundant\tstart\tstop\tstrand\tsequence";
#    foreach my $bamListElement (@{$bamList}) {
#	my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamListElement};
#	print HDR "\t$sample";
#        print HPL "\t$sample";
#    }
#    print HDR "\n";
#    print HPL "\n";
#    #$tk->start("TotalTime") if ($testTime);
    foreach my $chrom (keys %{$mirHairpins}) {
	#print "$chrom\n";
	my $genome = miRWoods::loadGenome("$genomeDir/$chrom.fa");	
	foreach my $hairpinInfo (@{$mirHairpins->{$chrom}}) {
	    my($foldStart,$foldStop,$foldStrand,$id,$name) = @{$hairpinInfo};
	    my $foldLocation = "$chrom:$foldStart..$foldStop:$foldStrand";
	    my $foldSequence = miRWoods::getSequence($foldLocation,$genome);
	    my($fold,$mfe) = RNA::fold($foldSequence);
	    #print "$name\t$foldLocation\n";

	    #need to find scores from predicted scores file.
       
#	my $maxProductScore;
#	if ($maxPredictedProductScores) {
#	    $maxProductScore = ($maxPredictedProductScores->{$id}) ? $maxPredictedProductScores->{$id} : 0;
#	}

	    #the chroms in RRF should be sorted at this point so the next line of code will only be executed once per chrom


	    my $location = miRWoods::extendSequenceToFillBuffer($chrom,$foldStart,$foldStop,$foldStrand,$chromLengths->{$chrom},$parameters);
	    #print "extended location = $location\n";
	    my($chrom,$start,$stop,$strand) = miRWoods::parseLocation($location);
	    my $leftOutRegion = '+' x ($foldStart - $start);
	    my $rightOutRegion = '+' x ($stop - $foldStop);
	    my $newFold = $leftOutRegion . $fold . $rightOutRegion;

#	    print ' ' x ($foldStart - $start) . "$foldSequence" . ' ' x ($stop - $foldStop) . "\n";
#	    print ' ' x ($foldStart - $start) . "$fold" . ' ' x ($stop - $foldStop) . "\n";

	    
	    my $sequence = miRWoods::getSequence($location,$genome);
#	    print "$sequence\n$newFold\n";
	    my($distinctReads,$uniqueReadCount,$totalDistinctReads,$adjustedReadCount) = miRWoods::retrieveReadData($location,$bamList,$parameters);
	    if ($distinctReads->{$strand}) {  #neccessary if testing read reagions with known mirs and no real read exists in that read region
		my($senseProducts) = miRWoods::buildProducts($location,$distinctReads->{$strand},$strand,$parameters);
		my $annotProducts = getAnnotatedProducts($id,$senseProducts,$mirProducts,$location);
#		printProducts($annotProducts,$location,$sequence,$newFold,$id,$genome);
		if (checkAnnotatedProducts($annotProducts)) {
		    my $likelyFolds = miRWoods::getFold($chrom,$sequence,$annotProducts,$id,$parameters);
		    #bestDupFold and bestDupMFE will not be used,  we will use the fold and mfe previously calculated.
		    my $foldInfo = pop(@{$likelyFolds});
		    my($bestDupFold,$bestDupMFE,$maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure,$duplexEnergy) = @{$foldInfo};
#		    print "\nfoldInfo:\n";
#		    print "$bestDupFold\n$bestDupMFE\n";
#		    print "maxProductCoords: \(". $maxProductCoords->[0] . "," . $maxProductCoords->[1] . "\)\n";
#		    print "productDuplexCoords: \(". $productDuplexCoords->[0] . "," . $productDuplexCoords->[1] . "\)\n";
#		    print "seqDuplexCoords: \(". $seqDuplexCoords->[0] . "," . $seqDuplexCoords->[1] . "\)\n";
#		    print $dupStructure . "\n";
#		    print "dupEnergy:" . $duplexEnergy . "\n";
		    my($averageProductVariance,$weightedAverageProductVariance) = miRWoods::getReadDataProductVariance($senseProducts,$parameters);
		    my($averageReadVariance,$weightedAverageReadVariance) = miRWoods::getReadLocationVariance($senseProducts,$parameters);
		    #here productDuplexCoords is the relative coordinates of the product after duplexing.  seqDuplexCoords is the actual produts duplex.
		    my($foldDuplexCmp,$exactFoldDuplexCmp) = miRWoods::compareFoldAndDuplexStructs($newFold,$maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure,$id);
		    my($prodDupPBP,$mpDuplexLoopLength) = miRWoods::getProductDuplexPBP($newFold,$maxProductCoords);
		    #$tk->end("getFold") if ($testTime);
		    #$tk->start("getMergedHairpinCenters") if ($testTime);
		    my @centers = miRWoods::getMergedHairpinCenters($newFold,$parameters);
		    if (@centers == 1) {
			my $center = $centers[0];
			my $newId = $name;
#		    #$tk->end("getMergedHairpinCenters") if ($testTime);
#		    $centersReadsCount += scalar(@centers);
#		    unless (scalar(@centers)) {
#			$noCenterCount++;
#			print "$id has no centers\n" if ($testReadCounts);
#			print "$newFold\n" if ($testReadCounts);
#		    }
#		    my $COUNT;
#		    foreach my $center (@centers) {
#			my $asciiVal = ord('a') + $COUNT;
#			my $label = chr($asciiVal);
#			my $newId = $id . $label;
#			$COUNT++;


			my $basePairs = miRWoods::getBasePairs($center,$newFold);
			my $hairpinLength = miRWoods::getMaxHairpinLength($newFold,$center);
			my($hStart,$hStop) = miRWoods::getFullHairpinWindow($newFold,$center);
			if($hairpinLength >= $minLength) {
			    my $senseProducts=miRWoods::getProductInfo($center,$newFold,$basePairs,$location,$senseProducts,$parameters);
			    $senseProducts=miRWoods::rebuildProducts($center,$newFold,$basePairs,$senseProducts,$parameters);
			    $senseProducts=miRWoods::addProductSize($senseProducts,$parameters);
			    my $revStrand = miRWoods::revStrand($strand);
			    my($antiSenseProducts) = miRWoods::buildProducts($location,$distinctReads->{$revStrand},$revStrand,$parameters);
			    $antiSenseProducts = miRWoods::getProductInfo($center,$newFold,$basePairs,$location,$antiSenseProducts,$parameters);
			    $antiSenseProducts = miRWoods::rebuildProducts($center,$newFold,$basePairs,$antiSenseProducts,$parameters);
			    $antiSenseProducts=miRWoods::addProductSize($antiSenseProducts,$parameters);
			    my $adjTotalProductReads = miRWoods::getAdjTotalProductReads($senseProducts);
			    my $adjTotalRevProductReads = miRWoods::getAdjTotalProductReads($antiSenseProducts);
			    my($PASS,$REASON) = miRWoods::plausibleReads($adjTotalProductReads,$adjTotalRevProductReads,$senseProducts,$antiSenseProducts,
							       $parameters);
			    my($maxProduct,$mpDuplex,$maxProductCount,$mpDuplexCount,$duplexOverlap) = miRWoods::getMaxAndDuplexProducts($senseProducts,$newFold,$id);
			    my($mpLoopDistance,$dupLoopDistance,$loopSize) = miRWoods::getMaxProductLoopDistance($center,$maxProduct,$mpDuplex);
			    my($mpOverlapScore,$dupOverlapScore,$wghtMPOverlapIntensity,$wghtDupOverlapIntensity) = miRWoods::getMaxProductOverlapScores($senseProducts,$maxProduct,$mpDuplex);
			    my($tpd,$totalRP)=miRWoods::newGetReverseProductDisplacement($senseProducts,
									       $antiSenseProducts,
									       $adjTotalProductReads,
									       $parameters,$newId);
			    my $apd = $totalRP ? $tpd/$totalRP : 0.0;
			    my $urf = $adjustedReadCount->{$strand} ? $uniqueReadCount->{$strand} / $adjustedReadCount->{$strand} : 0;
			    my $ahc = miRWoods::computeMaxProdHitCount($senseProducts,$location,$distinctReads,$parameters);
			    my $afh = miRWoods::computeMaxProdFivePrimeHet($senseProducts,$parameters);
			    my $pbp = miRWoods::computeProductBasePairing($center,$senseProducts,$basePairs,$parameters);
			    my $sameShift = miRWoods::computeMaxSameShift($location,$distinctReads->{$strand},$senseProducts,$parameters);
			    my $bothShift = miRWoods::computeMaxBothShift($basePairs,$location,$distinctReads->{$strand},$senseProducts,$parameters,$newId);
			    my $innerLoopGapCount = miRWoods::getInnerLoopGapCount($center,$newFold);
			    my $splitProductAbundance = miRWoods::computeSplitProductAbundance($senseProducts,$parameters);
			    my $overlapProductAbundance = miRWoods::computeOverlapProductAbundance($senseProducts,$parameters);
			    my $opa2 = miRWoods::computeOverlapProductAbundance2($senseProducts,$parameters);
			    my ($totalOverlapAmount,$averageOverlapAmount,
				$totalRelativeOverlapAmount,$averageRelativeOverlapAmount) = miRWoods::computeOverlapAmounts($senseProducts,$parameters);
			    my $maxOverlap = miRWoods::getMaximumOverlap($senseProducts,$parameters);
			    my($aa,$ac,$ag,$at,$ca,$cc,$cg,$ct,$ga,$gc,$gg,$gt,$ta,$tc,$tg,$tt) = miRWoods::getDinucleotideFrequencies($sequence);
			    my $sequenceLength = length($sequence);
			    my $gcContent = miRWoods::getGCcontent($sequence,$sequenceLength);
			    # total is the read count for the whole hairpin
			    my $totalSense = 0;
			    foreach my $product (@{$senseProducts}) {
				my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
				   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
#				my $length = $relStop - $relStart + 1;
#				my $productSequence = substr($sequence,$relStart,$length);
				$totalSense += $adjProductCount;
#				my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($productList);
#				print HPL "$newId\t$side\t$newType\t$adjProductCount\t";
#				print HPL "$adjMaxProductCount\t$relStart\t$relStop\t$productStrand\t";
#				print HPL "$productSequence";
#				foreach my $bamElement (@{$bamList}) {
#				    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
#				    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
#				}
#				print HPL "\n";
#				foreach my $read (@{$productList}) {
#				    my($relStart,$relStop,$offset,$gStart,
#				       $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
#				    my $gStop = $gStart + ($relStop - $relStart);
#				    my $adjCount = $count * $adjSeqCount;
#				    print HDR "$newId\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
#				    print HDR "$offset\t$relStart\t$relStop\t$seq";
#				    foreach my $bamElement (@{$bamList}) {
#					my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
#					printf(HDR "\t%.3f",$adjLibraryCounts->{$sample});
#				    }
#				    print HDR "\n";
#				}				
			    }
#			    #############################
#			    # ANTISENSE PRODUCTS
			    my $totalAntisense = 0;
			    foreach my $product (@{$antiSenseProducts}) {
				my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
				   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
#				$newType = "as-".$newType;
				$totalAntisense += $adjProductCount;
#				my $length = $relStop - $relStart + 1;
#				my $productSequence = reverseComplement(substr($sequence,$relStart,$length));
#				my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($productList);
#				print HPL "$newId\t$side\t$newType\t$adjProductCount\t";
#				print HPL "$adjMaxProductCount\t$relStart\t$relStop\t$productStrand\t";
#				print HPL "$productSequence";
#				foreach my $bamElement (@{$bamList}) {
#				    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
#				    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
#				}
#				print HPL "\n";
#				foreach my $read (@{$productList}) {
#				    my($relStart,$relStop,$offset,$gStart,
#				       $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
#				    my $gStop = $gStart + ($relStop - $relStart);
#				    my $adjCount = $count * $adjSeqCount;
#				    print HDR "$newId\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
#				    print HDR "$offset\t$relStart\t$relStop\t$seq";
#				    foreach my $bamElement (@{$bamList}) {
#					my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
#					printf(HDR "\t%.3f",$adjLibraryCounts->{$sample});
#				    }
#				    print HDR "\n";
#				}				
			    }
#			    #############################
#			    #############################

			    my @maxAnnotProductInfo = @{$annotProducts->[0]};
			    my $maxAnnotRelStart = $maxAnnotProductInfo[6];
			    my $maxAnnotRelStop = $maxAnnotProductInfo[7];
			    my $maxAnnotStrand = $maxAnnotProductInfo[9];
#			    print "$maxAnnotStrand\n";
			    my($maxAnnotGStart,$maxAnnotGStop) = miRWoods::getGlobalStartStop($location,$maxAnnotRelStart,$maxAnnotRelStop,$maxAnnotStrand);
			    my $maxAnnotLocation = "$chrom:$maxAnnotGStart..$maxAnnotGStop:$maxAnnotStrand";
			    my $maxProductScore;
			    if ($predProductScores->{$chrom}{$maxAnnotStrand}{$maxAnnotGStart}) {
				if ($predProductScores->{$chrom}{$maxAnnotStrand}{$maxAnnotGStart} > 0.5) {
				    #this is a temporary solution to some products not being in the class file.
				    #if this function ends up being used as part of the final program then this problem will be fixed.
				    
				    $maxProductScore = $predProductScores->{$chrom}{$maxAnnotStrand}{$maxAnnotGStart};
#				print "Score = $maxProductScore\n";
				

				    my($leftCenter, $rightCenter) = @{$center};
#				print HL "$newId\t$chrom\t$start\t$stop\t$strand\t$leftCenter\t$rightCenter\t";
#				printf(HL "%.3f\t%.3f\t%.3f\t",$totalSense,$totalAntisense,$mfe);
#				print HL "$sequence\t$newFold\t";
#				printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\t",
#				       $apd,$tpd,$urf,$ahc,$afh,$pbp,$sameShift,$bothShift,$splitProductAbundance,$overlapProductAbundance);
#				printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t",
#				       $aa,$ac,$ag,$at,$ca,$cc,$cg,$ct,$ga,$gc,$gg,$gt,$ta,$tc,$tg,$tt,$duplexEnergy,$gcContent);
#				printf(HL "%.3f\t%.3f\t%.3f\t%d\t",
#				       $foldDuplexCmp,$exactFoldDuplexCmp,$prodDupPBP,$mpDuplexLoopLength);
##			printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%d\t",
##			       $zScore, $foldDuplexCmp,$exactFoldDuplexCmp,$prodDupPBP,$mpDuplexLoopLength);
#				printf(HL "%.3f\t%.3f",
#				   $averageProductVariance,$weightedAverageProductVariance);
#				printf(HL "\t%.3f\t%.3f",
#				       $averageReadVariance,$weightedAverageReadVariance);
#				printf(HL "\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
#				       $maxProductCount,$mpDuplexCount,$duplexOverlap,$mpLoopDistance,$dupLoopDistance,$loopSize,$mpOverlapScore,$dupOverlapScore,$wghtMPOverlapIntensity,$wghtDupOverlapIntensity);
#				printf(HL "\t%.3f",
#				       $opa2);
#				printf(HL "\t%.3f\t%.3f\t%.3f\t%.3f", $totalOverlapAmount,$averageOverlapAmount,$totalRelativeOverlapAmount,$averageRelativeOverlapAmount);
#				printf(HL "\t%.3f", $maxOverlap);
#				printf(HL "\t%.3f",  $innerLoopGapCount);
				    my $totalSenseRPM = miRWoods::getHairpinRPM($totalSense,$librarySizes);
				    my $totalAntisenseRPM = miRWoods::getHairpinRPM($totalAntisense,$librarySizes);
#				printf(HL "\t%.3f\t%.3f", $totalSenseRPM,$totalAntisenseRPM);
##			    if ($maxProductScore) {
#				printf(HL "\t%.3f\n", $maxProductScore);
##			    } else {
##				print HL "\n";
##			    }
				    print HTF "$newId\t1\t$totalSense\t$totalAntisense\t$mfe\t$apd\t$tpd\t$urf\t$ahc\t$afh\t$pbp\t$sameShift\t$bothShift\t$splitProductAbundance\t$overlapProductAbundance\t";
				    print HTF "$aa\t$ac\t$ag\t$at\t$ca\t$cc\t$cg\t$ct\t$ga\t$gc\t$gg\t$gt\t$ta\t$tc\t$tg\t$tt\t";
				    print HTF "$duplexEnergy\t$gcContent";
#		print HTF "\t$zScore";
				    print HTF "\t$foldDuplexCmp\t$exactFoldDuplexCmp\t$prodDupPBP\t$mpDuplexLoopLength";
				    print HTF "\t$averageProductVariance\t$weightedAverageProductVariance\t$averageReadVariance\t$weightedAverageReadVariance";
				    print HTF "\t$maxProductCount\t$mpDuplexCount\t$duplexOverlap\t$mpLoopDistance\t$dupLoopDistance\t$loopSize\t$mpOverlapScore\t$dupOverlapScore\t$wghtMPOverlapIntensity\t$wghtDupOverlapIntensity\t$opa2\t$totalOverlapAmount\t$averageOverlapAmount\t$totalRelativeOverlapAmount\t$averageRelativeOverlapAmount\t$maxOverlap\t$innerLoopGapCount\t$totalSenseRPM\t$totalAntisenseRPM";
				    print HTF "\t$maxProductScore";
				    print HTF "\n";
				}
			    }
			}
		    }
		}
	    }
	} 
    }
}

sub readPredProductClassFile {
    my($predProductClassFile) = @_;
    my %predProductScores;
    open(PPCF,$predProductClassFile) or die "failed to open $predProductClassFile for reading\n";
    while (<PPCF>) {
	chomp;
	my($geneId,$prediction,$score) = split("\t");
	my($rrId,$location) = miRWoods::parseGeneId($geneId);
	my($chrom,$start,$stop,$strand) = miRWoods::parseLocation($location);
	$predProductScores{$chrom}{$strand}{$start} = $score;
    }
    close(PPCF);
    return \%predProductScores;
}

sub getAnnotatedProducts {
    my($parentId,$senseProducts,$miRProducts,$location) = @_;
    my @annotatedProducts;
    for my $annotProductInfo (@{$miRProducts->{$parentId}}) {
	my($chrom,$start,$stop,$strand,$id,$name) = @{$annotProductInfo};
	my $productFound = 0;
	foreach my $productInfo (@{$senseProducts}) {
	    my($prodId,$productReads,$productCount,$maxProductCount,$adjustedProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$productInfo};
	    if ($productStrand eq $strand) {
		my($gStart,$gStop) = miRWoods::getGlobalStartStop($location,$maxRelStart,$maxRelStop,$productStrand);
		if ($gStart == $start) {
		    push(@annotatedProducts,[$name,$productReads,$productCount,$maxProductCount,$adjustedProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand]);
#		    print "pushing $name onto annotated Products\n";
#		    print "$name\t$productReads\t$adjMaxProductCount\t$adjustedProductCount\t$maxRelStart\t$maxRelStop\t$maxGStart\t$productStrand\n";
#		    print "$maxRelStart - $maxRelStop\n";
		    $productFound = 1;
		}
	    }
	}
	unless ($productFound) {
#	    print "$name not Found\n";
	    my($relStart,$relStop) = miRWoods::getRelStartStop($location,$start,$stop,$strand);
	    my $productReads = [];
	    push(@annotatedProducts,[$name,$productReads,0,0,0,0,$relStart,$relStop,$start,$strand]);
#	    print "$name\t$productReads\n";
#	    print "$relStart - $relStop\n";
	}
    }
    my @sortedAnnotProducts = sort {$b->[4] <=> $a->[4]} @annotatedProducts;
    return \@sortedAnnotProducts;
}

sub printProducts {
    #a test function to print out the products
    my($products,$location,$sequence,$fold,$hairpinId,$genome) = @_;
    my($chrom,$start,$stop,$strand) = miRWoods::parseLocation($location);
    my @sortedProducts = sort {$a->[6] <=> $b->[6]} @{$products};  #sorting by relative position
    print "$hairpinId\t$location\n";
    my $lastPosition = 0;
    foreach my $productInfo (@sortedProducts) {
	my($name,$productReads,$productCount,$maxProductCount,$adjustedProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$productInfo};
	my($gStart,$gStop) = miRWoods::getGlobalStartStop($location,$maxRelStart,$maxRelStop,$productStrand);
	my $productLocation = "$chrom:$gStart..$gStop:$productStrand";
	my $productSequence = miRWoods::getSequence($productLocation,$genome);
	print ' ' x ($maxRelStart - $lastPosition) . "$productSequence";
	$lastPosition = $maxRelStart + length($productSequence);
    }
    print "\n$sequence\n$fold\n";
}

sub checkAnnotatedProducts {
    my($annotatedProducts) = @_;
    for my $productInfo (@{$annotatedProducts}) {
	my($name,$productReads,$productCount,$maxProductCount,$adjustedProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$productInfo};
	if ($adjustedProductCount) {
	    #print "product $name found: $adjustedProductCount\n";
	    return 1;
	}
    }
    return 0;
}

