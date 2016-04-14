package miRWoods;
use Bio::DB::Sam;
use Memory::Usage;
use Statistics::R;
#use TimeKeeper;
use RNA;
use strict;
use warnings;
use threads;
use threads::shared;
use List::Util;


sub loadDefaultParameters {
    my $parameters = {
	"aapdLimit" => 3,
	"ahcLimit" => 5,
	"afhLimit" => 0.5,
	"minSameShift" => 7,
	"minBothShift" => 7,
	"NEIGHBOR" => 1,
	"CHECK_EXONS" => 0,
	"nonMirNeighborLimit" => 10,
	"neighborWindow" => 1000,
	"minGap" => 3,
	"minMirCount" => 1,  
	"MaxLength" => 160,
	"lengthMin" => 20,
	"totalLength" => 150,
	"distanceMin" => 10,
	"CountMax" => 5,
	"HitCountMax" => 50,
	"MaxMajor" => 22,
	"minMajor" => -44,
	"reverseMax" => 0.05,
	"countMinLocus" => 0,   
	"fivePrimeHetMax" => 0.5,
	"shiftMin" => 7,
	"OverlapMin" => 2,
	"BPdensityLimit" => 0.6,
	"hairpinShortLength" => 12,
	"InHairpinBuffer" => 3,
	"OutHairpinBuffer" => 3,
	"RangeOfHairpin" => 70,
	"outputPrefix" => "readRegions",
	"mirbaseGff" => "",
	"bamListFile" => "",
	"RepeatRegionsFile" => "",
	"genomeDir" => "",
	"SizesFile" => "",
	"geneModels" => "",
	"validList" => "",
	"knownMirsFile" => "",
	"MaxShortProductSize" => "14",
	"minLongProductSize" => "28",
	"LoadFromConfigFile" => "",
	"maxSplitProductAbundance" => 0.05,
	"maxOverlapProductAbundance" => 0.05,
	"minAntisenseProduct" => 0,   #the minimum amount of antisense product required to incorporate a sense antisense overlapping pair into the aapd
	"MLAlgorithm" => "RandomForest",
    };
    return $parameters;
}

sub loadDefaultMLProductParameters {
    #these parameters determine which variables will be used in machine learning for products
    my $prodMLParameters = {
	"adjMaxProductCount" => 0,
	"adjProductCount"=> 0,
	"fivePrimeHet"=> 1,
	"length"=> 1,
	"medianLength"=> 1,
	"gcContent"=> 1,
	"aa"=> 1,
	"ac"=> 1,
	"ag"=> 1,
	"at"=> 1,
	"ca"=> 1,
	"cc"=> 1,
	"cg"=> 1,
	"ct"=> 1,
	"ga"=> 1,
	"gc"=> 1,
	"gg"=> 1,
	"gt"=> 1,
	"ta"=> 1,
	"tc"=> 1,
	"tg"=> 1,
	"tt"=> 1,
	"r10"=> 0,
	"r9"=> 0,
	"r8"=> 0,
	"r7"=> 1,
	"r6"=> 1,
	"r5"=> 1,
	"r4"=> 1,
	"r3"=> 1,
	"r2"=> 1,
	"r1"=> 1,
	"s0"=> 1,
	"f1"=> 1,
	"f2"=> 1,
	"f3"=> 1,
	"f4"=> 1,
	"f5"=> 1,
	"f6"=> 1,
	"f7"=> 1,
	"f8"=> 0,
	"f9"=> 0,
	"f10"=> 0,
	"WFC"=> 1,
	"duplexEnergy"=> 1
    };
    return $prodMLParameters;
}

sub loadDefaultMLHairpinParameters {
    #these parameters determine which variables will be used in machine learning for products
    my $hairpinMLParameters = {
	"totalSense"=> 0,
	"totalAntisense"=> 0,
	"mfe"=> 1,
	"aapd"=> 1,
	"tapd"=> 1,
	"urf"=> 1,
	"ahc"=> 1,
	"afh"=> 1,
	"pbp"=> 1,
	"sameShift"=> 1,
	"bothShift"=> 1,
	"SPA"=> 0,
	"OPA"=> 0,
	"aa"=> 1,
	"ac"=> 1,
	"ag"=> 1,
	"at"=> 1,
	"ca"=> 1,
	"cc"=> 1,
	"cg"=> 1,
	"ct"=> 1,
	"ga"=> 1,
	"gc"=> 1,
	"gg"=> 1,
	"gt"=> 1,
	"ta"=> 1,
	"tc"=> 1,
	"tg"=> 1,
	"tt"=> 1,
	"duplexEnergy"=> 1,
	"GCcontent"=> 1,
	"foldDupCmp"=> 0,
	"exactFoldDupCmp"=> 0,
	"dupPBP"=> 0,
	"dupLoopLength"=> 0,
	"APV"=> 0,
	"wAPV"=> 0,
	"ARV"=> 0,
	"wARV"=> 0,
	"mpCount"=> 0,
	"dupCount"=> 0,
	"dupOverlap"=> 0,
	"mpLoopDistance"=> 0,
	"dupLoopDistance"=> 0,
	"loopSize"=> 0,
	"mpOverlapScore"=> 0,
	"dupOverlapScore"=> 0,
	"wghtMPOverlapIntensity"=> 0,
	"wghtDupOverlapIntensity"=> 0,
	"OPA2"=> 0,
	"totalOverlapAmount"=> 0,
	"averageOverlapAmount"=> 0,
	"totalRelativeOverlapAmount"=> 0,
	"averageRelativeOverlapAmount"=> 0,
	"maxOverlap"=> 0,
	"innerLoopGapCount"=> 0,
	"totalSenseRPM" => 1,
	"totalAntisenseRPM" => 1,
	"RFProductAvg"=> 1
    };
    return $hairpinMLParameters;
}

sub loadDefaultMLHairpinParameters_sizeAdj {
    #these parameters determine which variables will be used in machine learning for products
    my $hairpinMLParameters = {
	"totalSense"=> 0,
	"totalAntisense"=> 0,
	"mfe"=> 1,
	"aapd"=> 1,
	"tapd"=> 1,
	"urf"=> 1,
	"ahc"=> 1,
	"afh"=> 1,
	"pbp"=> 1,
	"sameShift"=> 1,
	"bothShift"=> 1,
	"SPA"=> 1,
	"OPA"=> 1,
	"aa"=> 1,
	"ac"=> 1,
	"ag"=> 1,
	"at"=> 1,
	"ca"=> 1,
	"cc"=> 1,
	"cg"=> 1,
	"ct"=> 1,
	"ga"=> 1,
	"gc"=> 1,
	"gg"=> 1,
	"gt"=> 1,
	"ta"=> 1,
	"tc"=> 1,
	"tg"=> 1,
	"tt"=> 1,
	"duplexEnergy"=> 1,
	"GCcontent"=> 1,
	"foldDupCmp"=> 1,
	"exactFoldDupCmp"=> 1,
	"dupPBP"=> 1,
	"dupLoopLength"=> 1,
	"APV"=> 1,
	"wAPV"=> 1,
	"ARV"=> 1,
	"wARV"=> 1,
	"mpCount"=> 0,
	"dupCount"=> 0,
	"dupOverlap"=> 1,
	"mpLoopDistance"=> 1,
	"dupLoopDistance"=> 1,
	"loopSize"=> 1,
	"mpOverlapScore"=> 0,
	"dupOverlapScore"=> 0,
	"wghtMPOverlapIntensity"=> 1,
	"wghtDupOverlapIntensity"=> 1,
	"OPA2"=> 0,
	"totalOverlapAmount"=> 1,
	"averageOverlapAmount"=> 1,
	"totalRelativeOverlapAmount"=> 1,
	"averageRelativeOverlapAmount"=> 1,
	"maxOverlap"=> 1,
	"innerLoopGapCount"=> 1,
	"totalSenseRPM" => 0,
	"totalAntisenseRPM" => 0,
	"RFProductAvg"=> 1
    };
    return $hairpinMLParameters;
}

sub readMLProdConfigFile {
    my($prodMLConfigFile,$prodMLParameters) = @_;
    open(CONFIGFILE,$prodMLConfigFile) or die "FAIL: could not open $prodMLConfigFile\n";
    while(<CONFIGFILE>) {
	chomp;
	my($key,$value) = split(/\s+=\s+/);
	$prodMLParameters->{$key} = $value;
    }
    return $prodMLParameters;
}

sub createOutputFileParameters {
    my($parameters) = shift;
    my $filePrefix = $parameters->{outputPrefix} or die "FAIL: filePrefix not loaded.\n";
    $parameters->{distancesFile} = $filePrefix."_distances.txt";
    $parameters->{readRegionsFile} = $filePrefix.".txt";
    $parameters->{longReadRegionsFile} = $filePrefix."_long.txt";
    $parameters->{allReadRegionsFile} = $filePrefix."_full.txt";
    $parameters->{readRegionsFastaFile} = $filePrefix . ".fasta";
    $parameters->{tRNAScanFasta} = $filePrefix . "_trnaScanFasta.fasta";
    $parameters->{tRNAScanOutputFile} = $filePrefix . ".trna";
    $parameters->{filteredCandidatesFile} = $filePrefix . "_filteredCandidates.txt";
    $parameters->{hairpinsFile} = $filePrefix . "_hairpins.txt";
    $parameters->{distinctReadsFile} = $filePrefix . "_distinctReads.txt";
    $parameters->{productFile} = $filePrefix . "_products.txt";    
    $parameters->{featureVectorFile} = $filePrefix . "_scores.txt";
    $parameters->{productVectorFile} = $filePrefix . "_productScores.txt";
    $parameters->{productTrainFile} = $filePrefix . "_productTrainFile.txt";
    $parameters->{predProductFile} = $filePrefix . "_predProductRF.txt";
    $parameters->{predProductClasses} = $filePrefix . "_predProductClasses.txt";
    $parameters->{hairpinVectorFile} = $filePrefix . "_hairpinScores.txt";
    $parameters->{newHairpinVectorFile} = $filePrefix . "_newHairpinScores.txt";
    $parameters->{hairpinTrainFile} = $filePrefix . "_hairpinTrainFile.txt";
    $parameters->{hairpinTrainFileSB} = $filePrefix . "_hairpinTrainFileSB.txt";
    $parameters->{predHairpinClasses} = $filePrefix . "_predHairpinClasses.txt";
    $parameters->{newPredHairpinClasses} = $filePrefix . "_newPredHairpinClasses.txt";
    return $parameters;
}

sub readConfigFile {
    my($configFile,$parameters) = @_;
    open(CONFIGFILE,$configFile) or die "FAIL: could not open $configFile\n";
    while(<CONFIGFILE>) {
	chomp;
	my($key,$value) = split(/\s+=\s+/);
	$parameters->{$key} = $value;
    }
    return $parameters;
}


######################################
# SEQUENCE / GENOME TOOLS            #
######################################

sub mapBamStrand {
    my $strand = shift;
    if($strand eq "-1") {
	return "-";
    }
    return "+";
}

sub revStrand {
    my $strand = shift;
    if($strand eq "-") {
	return "+";
    }
    return "-";
}


sub reverseComplement {
# Returns the reverse complement of the input sequence.
    my($seq)=@_;
    $seq =~ tr/acgtuACGTU/tgcaaTGCAA/;
    $seq=reverse($seq);
    return $seq;
}

sub parseGeneId {
    my($geneId) = @_;
    my($name,$location) = $geneId =~ /(^.*?)\_(.*)$/;
    return($name,$location);
}

sub parseLocation {
    my($location)=@_;
    my($chrom,$start,$end,$strand);    
    if($location =~ /(.*)\:(-?\d+)\-(-?\d+)\:(.*)/) {
        $chrom=$1;
        $start=$2;
        $end=$3;
        $strand=$4;
    } elsif($location =~ /(.*)\:(-?\d+)\-(-?\d+)/) {
        $chrom=$1;
        $start=$2;
        $end=$3;
        $strand="+";
    } elsif($location =~ /(.*)\:(-?\d+)\.\.(-?\d+)\:(.*)/) {
	$chrom=$1;
        $start=$2;
        $end=$3;
        $strand=$4;   	
    } elsif($location =~ /(.*)\:(-?\d+)\.\.(-?\d+)/) {
	$chrom=$1;
        $start=$2;
        $end=$3;
        $strand="+";   	
    }
    return ($chrom,$start,$end,$strand);
}    

sub loadGenome {
    my $genomeFile = shift;

    my %sequences;
    my $tag;

    open(GFILE,$genomeFile) or die "could not open $genomeFile\n";
    while(<GFILE>) {
	s/\r?\n|\r/\n/g; # in case of dos or mac
        chomp;
        if(/>/) {
            s/>//g;
            my @terms = split;
	    $tag = shift(@terms);
        } else {
            $sequences{$tag} .= $_;
        }
    }
    return \%sequences;
}

sub getSequence {
    my($location,$sequences)=@_;
    my($chrom,$start,$stop,$strand)=parseLocation($location);
    if($sequences->{$chrom}) {
	# converting from 1-based to 0-based below:
        my $string =  substr($sequences->{$chrom},$start-1,$stop-$start+1);
        if($strand eq "+") {
            return $string;
        } else {
	    return reverseComplement($string);
        }
    } else {
        die "could not find $chrom in function getSequence\nLoaded from $location\n";
    }
}


sub readChromLengths {
    my($chromSizesFile) = @_;
    my %chromLengths;
    open(CSF, "$chromSizesFile");
    while (<CSF>) {
	chomp;
	my($chrom,$size) = split(/\s+/);
	$chromLengths{$chrom} = $size;
    }
    close(CSF);
    return \%chromLengths;
}



############################
# BINARY SEARCH TOOLS      #
############################

sub significantlyInsideList {
    # must overlap at least 50% of max region.
    my($start,$stop,$list) = @_;
    if(@{$list}) {
	my($searchStart,$searchStop) = getSearchStart($start,$stop,$list);
	# for region in the list:
	for(my $j=$searchStart;$j<=$searchStop;$j++) {
	    my($thisStart,$thisStop,$foldInfo) = @{$list->[$j]};
	    # check to see if it overlaps the given coords.
	    if(significantOverlap($start,$stop,$thisStart,$thisStop)) {
		return 1;
	    }
 	}
    }
    return 0;
}

sub significantOverlap {
    my($start,$stop,$thisStart,$thisStop) = @_;
    if(lociContained($start,$stop,$thisStart,$thisStop)) {
	return 1;
    }
    my $overlap = getOverlap($start,$stop,$thisStart,$thisStop);
    if($overlap) {
	my $maxLen = max($stop-$start,$thisStop-$thisStart);
	if($overlap/$maxLen > 0.35) {
	    return 1;
	}
    }
    return 0;
}

sub insideList {
    my($start,$stop,$list) = @_;
    #print "testing $start $stop\n";
    #print "size of list: ", scalar(@{$list}), "\n";
    if($list) {
	if(@{$list}) {
	    my($searchStart,$searchStop) = getSearchStart($start,$stop,$list);
	    # for region in the list:
	    #print "searching from $searchStart to $searchStop\n";
	    for(my $j=$searchStart;$j<=$searchStop;$j++) {
		my($thisStart,$thisStop) = @{$list->[$j]};
		# check to see if it overlaps the given coords.
		if((($thisStart <= $start)&&($start <= $thisStop))||
		   (($thisStart <= $stop)&&($stop <= $thisStop))||
		   (($start <= $thisStart)&&($thisStart <= $stop))||
		   (($start <= $thisStop)&&($thisStop <= $stop))) {		
		    return 1;
		}
	    }
	}
    }
    return 0;
}

sub getSearchStart {
    my($start,$stop,$list) = @_;
    # searchStart and searchStop are INDICES for the repeat regions in this context
    my $searchStart = 0;
    my $searchStop = scalar(@{$list})-1;    
    return binarySearch($list,$searchStart,$searchStop,$start,$stop);
    #return ($searchStart,$searchStop);
}

sub binarySearch {
    my($list,$searchStart,$searchStop,$start,$stop) = @_;    
    # this 200 is arbitrary. I don't trust it...
    # if the separation is significantly large, try and shorten it.
    if($searchStop - $searchStart > 200) {
        my $searchMiddle = int(($searchStart+$searchStop)/2);
        my($fStart,$fStop) = @{$list->[$searchStart]};
        my($lStart,$lStop) = @{$list->[$searchStop]};
        my($mStart,$mStop) = @{$list->[$searchMiddle]};
        if(($fStart <= $start)&&($stop <= $mStop)) {
            return binarySearch($list,$searchStart,$searchMiddle,$start,$stop);
        } elsif(($mStart  <= $start)&&($stop <= $lStop)) {
            return binarySearch($list,$searchMiddle,$searchStop,$start,$stop);
        } else {
            return ($searchStart,$searchStop);
        }
    } else {
        return ($searchStart,$searchStop);
    }
}

sub insertEntry {
    my($list,$entry) = @_;
    if(@{$list}) {
	my $i = @{$list}-1;
	while(($i >= 0)&&($list->[$i]->[0] > $entry->[0])) {
	    $list->[$i+1] = $list->[$i];
	    $i--;
	}
	$list->[$i+1] = $entry;	
    } else {
	push(@{$list},$entry);
    }
}

##############
# MATH TOOLS #
##############

sub max {
    my @array = @_;
    my @sortedArray = sort {$a <=> $b} @array;
    my $max = pop(@sortedArray);
    return $max;
}

sub min {
    my @array = @_;
    my @sortedArray = sort {$a <=> $b} @array;
    my $min = shift(@sortedArray);
    return $min;
}

sub isInt {
    my $num = shift;
    if ($num =~ /^\d+$/) {
	return 1;
    }
    return 0;
}

sub logK {
    my($num,$base) = @_;
    return log($num) / log($base); 
}


#########################
# FOLD PROCESSING TOOLS #
#########################

sub extractUniqueHairpinList {
    my $hairpinList = shift;
    my %uniqueHairpinList;
    foreach my $chrom (keys %{$hairpinList}) {
	# sort by total reads, greatest to lowest.
	my @sortedHairpinList = sort {$b->[6] <=> $a->[6]} @{$hairpinList->{$chrom}};
	my @uniqueList;
	foreach my $foldInfo (@sortedHairpinList) {
	    my($tag,$start,$stop,$strand,$hStart,$hStop,$total,$mfe,
	       $aapd,$tapd,$urf,$ahc,$afh,$valid,$shifted) = @{$foldInfo};               
	    my $absStart = $start + $hStart;
	    my $absStop = $start + $hStop;
	    if($strand eq "-") {
		$absStart = $stop - $hStop;
		$absStop = $stop - $hStart;
	    }
	    unless(significantlyInsideList($absStart,$absStop,\@uniqueList)) {
		insertEntry(\@uniqueList,[$absStart,$absStop,$foldInfo]);	
	    }
	}
	foreach my $uniqueFold (@uniqueList) {
	    my($absStart,$absStop,$foldInfo) = @{$uniqueFold};
	    push(@{$uniqueHairpinList{$chrom}},$foldInfo);
	}
    }
    return \%uniqueHairpinList
}

sub getHairpinCenters {
    my($fold,$parameters) = @_;
    my $INRUN=0;
    my $lastLeft = 0;
    my @centers;
    for(my $i=0;$i<length($fold);$i++)  {
        if(substr($fold,$i,1) eq "(") {
            $lastLeft = $i;
            $INRUN=1;
        } elsif(substr($fold,$i,1) eq ")") {
            if($INRUN) {
		push(@centers,[$lastLeft,$i]);
	    }
            $INRUN=0;
        }
    }
    return @centers;
}

sub getMergedHairpinCenters {
    my($fold,$parameters) = @_;
    my @centers = getHairpinCenters($fold,$parameters);
    my($mergedCenters,$MERGED) = mergeCenters(\@centers,$fold,$parameters);
    while($MERGED) {
	($mergedCenters,$MERGED) = mergeCenters($mergedCenters,$fold,$parameters);
    }
    return @{$mergedCenters};
}

sub mergeCenters {
    my($centers,$fold,$parameters) = @_;
    my $minLength = $parameters->{lengthMin};
    my $shortLength = $parameters->{hairpinShortLength};
    my @mergedCenters;
    my $MERGED = 0;
    for(my $i=0;$i<@{$centers};$i++) {
	#print "checking $center\n";
	my($leftEnd,$rightEnd) = extendHairpinToFullEnds($fold,$centers->[$i]);
	my($leftCenter, $rightCenter) = @{$centers->[$i]};
	my $leftLength = $leftCenter - $leftEnd + 1;
	my $rightLength = $rightEnd - $rightCenter + 1;
	if(($leftLength >= $minLength)&&($rightLength >= $shortLength)) {
	    # sufficiently large of a loop, keep as is.
	    push(@mergedCenters,$centers->[$i]);
	} elsif(($leftLength >= $minLength)&&($rightLength < $shortLength)) {
	    # potentially merge the loop.
	    if($i<@{$centers}-1) {
		# if there is a neighboring center
		my($leftEnd2,$rightEnd2) = extendHairpinToFullEnds($fold,$centers->[$i+1]);
		my($leftCenter2, $rightCenter2) = @{$centers->[$i+1]};
		my $leftLength2 = $leftCenter2 - $leftEnd2 + 1;
		my $rightLength2 = $rightEnd2 - $rightCenter2 + 1;
		if($leftLength2 < $shortLength) {
		    #print "merging $center $leftCenter2,$rightCenter2\n";
		    my($tightLeftEnd,$tightRightEnd) = extendHairpinToEnds($fold,$centers->[$i]);
		    my($tightLeftEnd2,$tightRightEnd2) = extendHairpinToEnds($fold,$centers->[$i+1]);
		    my $newLeftCenter = getNextLeftPar($fold,$tightLeftEnd);
		    my $newRightCenter = getNextRightPar($fold,$tightRightEnd2);
		    if(($newLeftCenter != -1)&&($newRightCenter != -1)) {
			push(@mergedCenters,[$newLeftCenter,$newRightCenter]);
			$MERGED = 1;
		    }
		}
	    } else {
		push(@mergedCenters,$centers->[$i]);
	    }
	} elsif(($leftLength >= $shortLength)&&($rightLength >= $minLength)) {
	    # sufficiently large of a loop, keep as is.
	    push(@mergedCenters,$centers->[$i]);
	} elsif(($leftLength < $shortLength)&&($rightLength >= $minLength)) {
	    # potentially merge the loop.
	    if($i > 0) {
		# if there is a neighboring center
		my($leftEnd2,$rightEnd2) = extendHairpinToFullEnds($fold,$centers->[$i-1]);
		my($leftCenter2, $rightCenter2) = @{$centers->[$i-1]};
		my $leftLength2 = $leftCenter2 - $leftEnd2 + 1;
		my $rightLength2 = $rightEnd2 - $rightCenter2 + 1;
		if($rightLength2 < $shortLength) {
		    #print "merging $center $leftCenter2,$rightCenter2\n";
		    my($tightLeftEnd,$tightRightEnd) = extendHairpinToEnds($fold,$centers->[$i]);
		    my($tightLeftEnd2,$tightRightEnd2) = extendHairpinToEnds($fold,$centers->[$i-1]);
		    my $newLeftCenter = getNextLeftPar($fold,$tightLeftEnd2);
		    my $newRightCenter = getNextRightPar($fold,$tightRightEnd);
		    if(($newLeftCenter != -1)&&($newRightCenter != -1)) {
			push(@mergedCenters,[$newLeftCenter,$newRightCenter]);
			$MERGED = 1;
		    }
		}
	    } else {
		push(@mergedCenters,$centers->[$i]);
	    }
	}
    }
    my @filteredCenters;
    for(my $i=0;$i<@mergedCenters;$i++) {
	my($l1,$r1) = @{$mergedCenters[$i]};
	my $N = @filteredCenters;
	my $REMOVE=0;
	for(my $j=0;$j<$N;$j++) {
	    my($l2,$r2) = @{$filteredCenters[$j]};
	    if(($l1==$l2)&&($r1==$r2)) {
		$REMOVE=1;
	    }
	}
	unless($REMOVE) {
	    push(@filteredCenters,[$l1,$r1]);
	}
    }
    return(\@filteredCenters,$MERGED);
}

sub getFullHairpinWindow {
    # this method returns the window of paired bases, extends beyond minor hairpins    
    my($fold,$center) = @_;;
    return extendHairpinToFullEnds($fold,$center);
}

sub getBasePairs {
    my($center,$fold) = @_;
    my($leftCenter, $rightCenter) = @{$center};
    my @basePairs;
    my $leftCurrent = $leftCenter+1;
    my $rightCurrent = $rightCenter-1;
    my $STAY=1;
    my $COUNT=0;
    while($STAY) {
        my $left = getNextLeftPar($fold,$leftCurrent);
        my $right = getNextRightPar($fold,$rightCurrent);
        if(($left == -1)||($right == -1)) {
            $STAY=0;
        } else {
            push(@basePairs,[$left,$right]);
            $COUNT++;
            $leftCurrent = $left;
            $rightCurrent = $right;
        }
    }
    if(@basePairs) {
        return \@basePairs;
    } else {
        die "crapped out at $center:\n$fold\n";
    }
}

sub getMaxHairpinLength {
    my($fold,$center) = @_;
    my($leftCenter, $rightCenter) = @{$center};
    my($leftEnd,$rightEnd) = extendHairpinToFullEnds($fold,$center);
    my $leftLength = $leftCenter - $leftEnd + 1;
    my $rightLength = $rightEnd - $rightCenter + 1;
    return max($leftLength,$rightLength);
}

sub extendHairpinToEnds {
    my($fold,$center) = @_;
    my($leftCenter, $rightCenter) = @{$center};
    my @basePairs;
    # start left current one up, so the true left middle will be found first.
    my $leftCurrent = $leftCenter+1;
    # start the search for the right par on the next base.
    my $rightCurrent = $rightCenter-1;
    my $STAY=1;
    my $COUNT=0;
    while($STAY) {
        my $left = getNextLeftPar($fold,$leftCurrent);
        my $right = getNextRightPar($fold,$rightCurrent);
        if(($left == -1)||($right == -1)) {
            $STAY=0;
	} else {
            push(@basePairs,[$left,$right]);
	    $COUNT++;
            $leftCurrent = $left;
            $rightCurrent = $right;
        }
    }
    unless(@basePairs) {
	my $L = "L" x ($leftCenter+1);
	my $R = "R" x ($rightCenter+1);
	die "extendHairpinToEnds: FAIL:\n$fold\t$L\n$R\n";
    }
    my $endPairs = pop(@basePairs);
    return @{$endPairs};
}

sub extendHairpinToFullEnds {
    my($fold,$center) = @_;
    my($leftCenter, $rightCenter) = @{$center};
    my @basePairs;
    # start left current one up, so the true left middle will be found first.  
    my $leftCurrent = $leftCenter+1;
    # start the search for the right par on the next base.
    my $rightCurrent = $rightCenter-1;
    my $STAY=1;
    my $COUNT=0;
    while($STAY) {
        my $left = getNextLeftPar($fold,$leftCurrent);
        my $right = getNextRightPar($fold,$rightCurrent);
        if(($left == -1)||($right == -1)) {
            $STAY=0;
        } else {
            push(@basePairs,[$left,$right]);
	    $COUNT++;
            $leftCurrent = $left;
            $rightCurrent = $right;
        }
    }
    unless(@basePairs) {
	my $L = "L" x ($leftCenter+1);
	my $R = "R" x ($rightCenter+1);
        die "extendHairpinFullToEnds FAIL: ($center)\n$fold\n$L\n$R\n\n";	
    }
    my $endPairs = pop(@basePairs);
    my($leftEnd,$rightEnd) = @{$endPairs};
    while(($leftEnd > 0)&&(substr($fold,$leftEnd-1,1) ne ")")) {
	$leftEnd--;
    }
    while(($rightEnd < length($fold)-1)&&(substr($fold,$rightEnd+1,1) ne "(")) {
	$rightEnd++;
    }
    while(substr($fold,$leftEnd,1) eq ".") {
	$leftEnd++;
    }
    while(substr($fold,$rightEnd,1) eq ".") {
	$rightEnd--;
    }
    return ($leftEnd,$rightEnd);
}

sub getNextLeftPar {
    my($fold,$start) = @_;
    for(my $i=$start-1;$i>=0;$i--) {
        if(substr($fold,$i,1) eq "(") {
            return $i;
        } elsif(substr($fold,$i,1) eq ")") {
            return -1;
        }
    }
    return -1;
}

sub getNextRightPar {
    my($fold,$start)=@_;
    for(my $i=$start+1;$i<length($fold);$i++) {
	if(substr($fold,$i,1) eq ")") {
            return $i;
        } elsif(substr($fold,$i,1) eq "(") {
	    return -1
	}
    }
    return -1;
}


#############################
# BAM RETRIEVAL SUBROUTINES #
#############################

sub loadBamFile {
    my $bamFile = shift;
    my $bam = Bio::DB::Bam->open( $bamFile );
    return $bam;
}

sub loadBamIndex {
    my $bamFile = shift;
    my $reIndex;  #changed to 1 if the index file doesn't exist
    my $bamIndex =  Bio::DB::Bam->index($bamFile,$reIndex);
    die "failed to load index for $bamFile\n" if ($reIndex);
    return $bamIndex;
}

sub loadBamList {
    my $bamListFile = shift;
    open(BLF,$bamListFile) or die "could not open $bamListFile\n";
    my @bamList;
    while(<BLF>) {
	chomp;
	my($label,$bamFile) = split;
	my $bam = loadBamFile($bamFile);
	my $bamHeader = $bam->header;
	my $bamIndex = loadBamIndex($bamFile);
	my $totalMapped = getBamTotal($bamFile);
	push(@bamList,[$bam,$bamHeader,$bamIndex,$bamFile,$label,$totalMapped]);
    }
    return \@bamList;
}

sub getBamTotal {
    my $bamFile = shift;
    my $bamFlagstatFile = $bamFile . ".flagstat";
    unless(-e $bamFlagstatFile) {
	# if flagstat file hasn't been created, make one.
	system("samtools flagstat $bamFile > $bamFlagstatFile");
    }    
    open(BFS,$bamFlagstatFile) or die "Could not open $bamFlagstatFile.\n";
    while(<BFS>) {
	if(/(\d+) \+ \d+ mapped/) {
	    close(BFS);
	    return $1;
	}
    }
    close(BFS);
    # if we made it here, something went wrong with parsing the flagstat file
    die "Error. Could not parse the flagstat file here: $bamFlagstatFile. Older version of samtools perhaps?\n"
}

sub retrieveReadData {
    my($location,$bamList,$parameters) = @_;
    my $minLen = ($parameters->{minReadLength}) ? $parameters->{minReadLength} : 0;
    my $maxLen = ($parameters->{maxReadLength}) ? $parameters->{maxReadLength} : ~0;
    my($chrom,$start,$stop,$strand) = miRTRAP1_7::parseLocation($location);
    $location = "$chrom:$start-$stop";  
    my %distinctReadCount;  #normalized by dividing by the hitCount for each read before totalling the number of reads
    my %libraryCounts;  #normalized by dividing by the hitCount for each read before totalling the number of reads
    my %adjustedLibraryCounts;
    my %hitCounts;
    my %distinctReads;
    my %readTotal;
    my %adjustedSeqCount;
    my @sampleList;
    my %uniqueReadCount;
    my %readCount;
    my %adjustedReadCount;   #normalized by dividing by the hitCount
    foreach my $bamData (@{$bamList}) {
	my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamData};
	my($tid,$chromStart,$chromStop) = $bamHeader->parse_region("$chrom:$start-$stop"); #converts start and end to zero based coordinates
	push(@sampleList, $sample);
	my $callBack = sub {
	    my($alignment,$data) = @_;
	    my($chromStart,$chromStop,$distinctReadCount,$libraryCounts,$adjustedLibraryCounts,$hitCounts,
	       $adjustedSeqCount,$readTotal,$uniqueReadCount,$readCount,$adjustedReadCount) = @{$data};
	    my $id = $alignment->qname;
	    my $rStart = $alignment->start;  #returned in 1 based coordinates
	    my $rStop = $alignment->end;  #returned in 1 based coordinates
	    my $seq = $alignment->qseq;
	    my $hitCount = $alignment->get_tag_values('NH');
	    my $rStrand = miRTRAP1_7::mapBamStrand($alignment->strand);
	    my $count;
	    if ($id =~ /.*_x(\d+)$/) {
		($count) = $id =~ /.*_x(\d+)$/;
	    } else {
		$count = 1;
	    }
	    $seq = miRTRAP1_7::reverseComplement($seq) if ($strand eq '-');
	    my $length = length($seq);
	    if ($length >= $minLen && $length <= $maxLen) {
		if(($chromStart <= $rStart)&&($rStop <= $chromStop)) {
		    $distinctReadCount->{$rStrand}{$rStart}{$seq} += $count;
		    $libraryCounts->{$rStrand}{$rStart}{$seq}{$sample} += $count;
		    $adjustedLibraryCounts->{$rStrand}{$rStart}{$seq}{$sample} += $count / $hitCount;
		    $hitCounts->{$seq} += $count * $hitCount;
		    $adjustedSeqCount->{$seq} += $count/$hitCount;
		    $readTotal->{$seq} += $count;
		    $uniqueReadCount->{$rStrand} += $count if ($hitCount == 1);
		    $readCount->{$rStrand} += $count;
		    $adjustedReadCount->{$rStrand} += $count / $hitCount;
		}
	    }
	};
	my $callBackData = [$chromStart,$chromStop,\%distinctReadCount,\%libraryCounts,\%adjustedLibraryCounts,\%hitCounts,
			    \%adjustedSeqCount,\%readTotal,\%uniqueReadCount,\%readCount,\%adjustedReadCount];
	my $code = $bamIndex->fetch($bam,$tid,$chromStart,$chromStop,$callBack,$callBackData);
    }

    foreach my $rStrand (keys %distinctReadCount) {
	$uniqueReadCount{$rStrand} = 0 unless ($uniqueReadCount{$rStrand});
	$readCount{$rStrand} = 0 unless ($readCount{$rStrand});
	$adjustedReadCount{$rStrand} = 0 unless ($adjustedReadCount{$rStrand});
	foreach my $rStart (keys %{$distinctReadCount{$rStrand}}) {
	    foreach my $seq (keys %{$distinctReadCount{$rStrand}{$rStart}}) {
		foreach my $sample (@sampleList) {
		    unless ($libraryCounts{$rStrand}{$rStart}{$seq}{$sample}) {
			$libraryCounts{$rStrand}{$rStart}{$seq}{$sample} = 0;
			$adjustedLibraryCounts{$rStrand}{$rStart}{$seq}{$sample} = 0;
		    }
		}
		my $avgHitCount = $hitCounts{$seq} / $readTotal{$seq};
		my $adjustedSeqCount = $adjustedSeqCount{$seq} / $readTotal{$seq};
		push(@{$distinctReads{$rStrand}}, [$rStart, $rStart + length($seq) - 1, $distinctReadCount{$rStrand}{$rStart}{$seq}, $avgHitCount, 
						   $libraryCounts{$rStrand}{$rStart}{$seq}, $seq, $adjustedSeqCount, 
						   $adjustedLibraryCounts{$rStrand}{$rStart}{$seq}]);
	    }
	}
    }
    # sort distinct reads by abundance.  This is important in getProductInfo
    foreach my $strand (keys %distinctReads) {
	@{$distinctReads{$strand}} = sort {$b->[2] <=> $a->[2]} @{$distinctReads{$strand}};
    }
    return \%distinctReads, \%uniqueReadCount, \%readCount, \%adjustedReadCount;
}


sub getTagValue {
    my $tags = shift;
    my $tagName = shift;
    my $tagType = shift;
    foreach my $tag (@{$tags}) {
	my($tagValue) = $tag =~ /^$tagName:$tagType:(.*)/;
	if ($tagValue) {
	    return $tagValue;
	}
    }
    return 0;
}

sub getDistinctReadCounts {
    my $location = shift;
    my $bamList = shift;
    my %readCounts;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    foreach my $bamEntry (@{$bamList}) {
	my($label, $bam, $total, $bamFile) = @{$bamEntry};	
	my @reads = $bam->get_features_by_location(-seq_id => $chrom, -start  => $start, -end => $stop);
	foreach my $read (@reads) {
	    my $rStart = $read->start;
	    my $rStop = $read->end;
	    # 5' end methdod.
	    my $rStrand = mapBamStrand($read->strand);
	    # we could remove this restriction later for requiring the same strand
	    if($strand eq $rStrand) {
		$readCounts{$rStart."x".$rStop}{$label}++;
	    }
	}
    }
    my @readCounts;
    foreach my $key (keys %readCounts) {
	my($rStart,$rStop) = $key =~ /(\d+)x(\d+)/;
	my @counts;
	my $totalCounts = 0;
	foreach my $bamEntry (@{$bamList}) {
	    my($label, $bam, $total, $bamFile) = @{$bamEntry};
	    my $count = $readCounts{$key}{$label} ? $readCounts{$key}{$label} : 0;
	    $totalCounts += $count;
	    push(@counts,$count);
	}
	push(@readCounts,[$rStart,$rStop,$totalCounts,\@counts]);
    }
    return \@readCounts;
}

###############################
# tRNA PROCESSING SUBROUTINES #
###############################

sub printTRNAFastaFile {
    my($hairpinsFile,$tRNAScanFastaFile) = @_;
    open(HP,"$hairpinsFile") or die "failed to open $hairpinsFile for reading in printTRNAFastaFile()";
    open(TRNAFF, ">$tRNAScanFastaFile") or die "failed to open $tRNAScanFastaFile for writing in printTRNAFastaFile()";
    while(<HP>) {
	my($id,$chrom,$start,$stop,$strand,$leftCenter,$rightCenter,$adjTotalReads,$mfe,$seq,$fold) = split(/\t/);
	my $location = "$chrom:$start..$stop:$strand";
	print TRNAFF ">$id\t$location\n$seq\n";
    }
    close(HP);
    close(TRNAFF);
}  


###############################
# READ PROCESSING SUBROUTINES #
###############################

sub printReadRegions {
    my $bamList = shift;
    my $chromLengths = shift;
    my $repeatRegions = shift;
    my $parameters = shift;
    my $maxLength = $parameters->{MaxLength} or die "FAIL: maxLength not loaded (not found in parameters.)\n";
    my $readRegionsFile = $parameters->{readRegionsFile} or die "FAIL: readRegionsFile not loaded (not found in parameters.)\n";
    my $longReadRegionsFile = $parameters->{longReadRegionsFile} or die "FAIL: longReadRegionsFile not loaded (not found in parameters.)\n";
    my $allReadRegionsFile = $parameters->{allReadRegionsFile} or die "FAIL: allReadRegionsFile not loaded (not found in parameters.)(\n";
    my $initWindowLength = 200000;
    my $shortReadsCount = 1;     # for regions shorter than 100bp.
    my $longReadsCount = 1;     # for regions longer than 100bp.
    my $repeatReadsCount = 1;     # for repeat regions.

    my $mu = Memory::Usage->new();

    open(RF,">".$readRegionsFile) or die "Failed to load $readRegionsFile for writing\n";
    open(LRF,">".$longReadRegionsFile) or die "Failed to load $longReadRegionsFile for writing\n";
    open(ARF,">".$allReadRegionsFile) or die "Failed to load $allReadRegionsFile for writing\n";
    foreach my $chrom (keys %{$chromLengths}) {
	my $start = 0;
	while($start <= $chromLengths->{$chrom}) {
	    # push out the windowLength until no reads extend beyond the ends.
	    my $windowLength = $initWindowLength <= $chromLengths->{$chrom} ? $initWindowLength : $chromLengths->{$chrom};
	    # now the window is defined so that all contiguous read regions are properly contained within it.
	    my($plusArray,$minusArray) = getReadCountArrays($bamList,$chrom,$start,$start+$windowLength,$parameters);
	    # we now have plus and minus arrays for thie window.
 	    my $posReadRegions = extractReadRegionList($plusArray,$start,$windowLength,"+");
 	    my $negReadRegions = extractReadRegionList($minusArray,$start,$windowLength,"-");
	    # sort readRegions in descending order by rrStart position
	    my @readRegions = sort {$b->[0] <=> $a->[0]} (@{$posReadRegions}, @{$negReadRegions});
	    my $newStart = $start+$windowLength;

	    $mu->record("readRegions obtained for $chrom:$start-$newStart");

	    foreach my $readRegion (@readRegions) {
		my($rrStart,$rrStop,$rrStrand,$rrMaxCount) = @{$readRegion};
		my $prrStart = $rrStart-1; # printable 0-based start for bed file
		if($rrStop < $newStart) {
		    unless(insideList($rrStart,$rrStop,$repeatRegions->{$chrom})) {
			# print the info.
			my $length = $rrStop - $rrStart;
			if($length <= $maxLength) {
			    print ARF "$chrom\t$prrStart\t$rrStop\trr$shortReadsCount\t$rrMaxCount\t$rrStrand\n";
			    print RF "$chrom\t$prrStart\t$rrStop\trr$shortReadsCount\t$rrMaxCount\t$rrStrand\n";
			    $shortReadsCount++;
			} else {
			    print ARF "$chrom\t$prrStart\t$rrStop\tlrr$longReadsCount\t$rrMaxCount\t$rrStrand\tRejected:Region too long\n";
			    print LRF "$chrom\t$prrStart\t$rrStop\tlrr$longReadsCount\t$rrMaxCount\t$rrStrand\n";
			    $longReadsCount++;
			}
		    } else {
			print ARF "$chrom\t$prrStart\t$rrStop\trrr$repeatReadsCount\t$rrMaxCount\t$rrStrand\tRejected: region within repeat.\n";
			$repeatReadsCount++;
		    }
		} else {
		    $newStart = $rrStart;
		}
	    }
	    # increment the start to beyond this window.
	    die "Huge read region at $chrom:$start.." . ($start+$windowLength) . "\nPlease report to miRTRAP developers.\n" if($start == $newStart);
	    # put the next window at beginning of last read region
	    $start = $newStart;
	}
    }
    close(RF);
    close(LRF);
    close(ARF);
    $mu->dump();
}

sub getReadCountArrays {
    my($bamList,$chrom,$windowStart,$windowStop,$parameters) = @_;
    my $minLen = ($parameters->{minReadLength}) ? $parameters->{minReadLength} : 0;
    my $maxLen = ($parameters->{maxReadLength}) ? $parameters->{maxReadLength} : ~0;
    #print "minLen = $minLen\t maxLen = $maxLen\n";
    my $location = "$chrom:$windowStart-$windowStop";
    my(@plusArray,@minusArray);
    foreach my $bamElement (@{$bamList}) {
	my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
	my @bamARGS = ("samtools", "view", $bamFile, $location);
	open(BAM, "-|", @bamARGS) or die "Failed to open pipe to @bamARGS";
	while (<BAM>) {
	    chomp;
	    unless( /^@/ ) {
		my($id, $flag, $refSeqName, $rStart, $mapq, $cigar, $mrnm, $mpos, $iSize, $seq, $qualityScore, @tags) = split(/\t/);
		my $readCount;
		if ($id =~ /.*_x(\d+)$/) {
		    ($readCount) = $id =~ /.*_x(\d+)$/;
		} else {
		    $readCount = 1;
		}
		my $rStrand = $flag & 16 ? '-': '+';  #bitwise anding to see if sequence is reverse complemented
		my $length = length($seq);
		my $rStop = $rStart + $length - 1;
		if ($length >= $minLen and $length <= $maxLen) {
		    if($rStart >= $windowStart) {
			if($rStrand eq "+") {
			    # in plus strand
			    for(my $i=$rStart;$i<=$rStop;$i++) {
				$plusArray[$i-$windowStart]++;     #generating arrays relative to windowStart
			    }
			} else {
			    # in minus strand
			    for(my $i=$rStart;$i<=$rStop;$i++) {
				$minusArray[$i-$windowStart]++;
			    }
			}
		    } else {
			die "Found read at $chrom:$windowStart..$windowStop that extends beyond beginning: $chrom:$rStart..$rStop:$rStrand\n";
		    }
		}
	    }
	}
    }
    return(\@plusArray,\@minusArray);
}

sub extractReadRegionList {
    my($readCountArray,$start,$windowLength,$strand) = @_;
    my $regionStart;
    my $regionStop;
    my @readRegions;
    my $INREGION = 0;
    my $maxCount = 0;
    for (my $j=0; $j<@{$readCountArray}; $j++) {
	if($readCountArray->[$j]) {
	    if($INREGION) {
		# extend the regionStop
		$regionStop = $j;
	    } else {
		# assume the region is only 1bp long initially.
		$regionStart = $j;
		$regionStop = $j;
		$INREGION = 1;
	    }
	    $maxCount = $readCountArray->[$j] if($maxCount < $readCountArray->[$j]);
	} else {
	    # no reads here, store regionious region if just leaving one.
	    if($INREGION) {
		push(@readRegions, [$regionStart+$start,$regionStop+$start,$strand,$maxCount]);
		$INREGION = 0;
		$maxCount = 0;
	    }
	}
    }
    # if still in region at the end of the loop push start and stop into readRegions array
    if($INREGION) {
	push(@readRegions, [$regionStart+$start,$regionStop+$start,$strand,$maxCount]);
	$INREGION = 0;
    }
    return \@readRegions;
}

sub extendSequenceToFillBuffer {
    my($chrom, $start, $stop, $strand, $chromLength, $parameters) = @_;
    my $bufferLength = getBufferLength($start,$stop,$parameters);
    my $newStart = $start-$bufferLength > 1 ? $start - $bufferLength : 1;
    my $newStop = $stop+$bufferLength <= $chromLength ? $stop+$bufferLength : $chromLength;
    my $newLocation = $chrom . ":" . $newStart . ".." . $newStop  . ":" . $strand;
    return $newLocation;
}

sub getBufferLength {
    my($start,$stop,$parameters) = @_;
    my $totalLength = $parameters->{totalLength};
    if($stop - $start + 1 > 50 ) {
	if($stop - $start + 1 < $totalLength) {
	    return int(($totalLength - ($stop - $start + 1))/2);
	} else {
	    return 0;
	}
    } else {
        return int($totalLength/2);
    }
}

sub combineProcessRRFiles {
    my($parameters) = @_;
    my $numProcesses = $parameters->{numProcesses};
    my $outputPrefix = $parameters->{outputPrefix};
    my $hairpinFile = $outputPrefix . "_hairpins.txt";
    my $productFile = $outputPrefix . "_products.txt";
    my $distinctReadsFile = $outputPrefix . "_distinctReads.txt";
    my $filteredCandidatesFile = $outputPrefix . "_filteredCandidates.txt";
    open(HPF,">$hairpinFile") or die "failed to open $hairpinFile for writing\n";
    open(PF,">$productFile") or die "failed to open $productFile for writing\n";
    open(DRF,">$distinctReadsFile") or die "failed to open $distinctReadsFile for writing\n";
    open(FCF, ">$filteredCandidatesFile") or die "failed to open $filteredCandidatesFile for writing\n";
    for (my $i = 0; $i < $numProcesses; $i++) {
	my $fragHairpinsFile = $outputPrefix . "_" . $i . "_hairpins.txt";
	my $fragProductFile = $outputPrefix . "_" . $i . "_products.txt";
	my $fragDistinctReadsFile = $outputPrefix . "_" . $i . "_distinctReads.txt";
	my $fragFilteredCandidatesFile = $outputPrefix . "_". $i . "_filteredCandidates.txt";
	open(FHPF,"$fragHairpinsFile") or die "failed to open $fragHairpinsFile\n";
	open(FPF,"$fragProductFile") or die "failed to open $fragProductFile\n";
	open(FDRF,"$fragDistinctReadsFile") or die "failed to open $fragDistinctReadsFile\n";
	open(FFCF,"$fragFilteredCandidatesFile") or die "failed to open $fragFilteredCandidatesFile\n";
	while(<FHPF>) {
	    chomp;
	    if ( /^#/ ) {
		print HPF $_ . "\n" if ($i == 0);
	    } else {
		print HPF $_ . "\n";
	    }
	}
	while(<FPF>) {
	    chomp;
	    if ( /^#/ ) {
		print PF $_ . "\n" if ($i == 0);
	    } else {
		print PF $_ . "\n";
	    }
	}
	while(<FDRF>) {
	    chomp;
	    if ( /^#/ ) {
		print DRF $_ . "\n" if ($i == 0);
	    } else {
		print DRF $_ . "\n";
	    }
	}
	while(<FFCF>) {
	    chomp;
	    if ( /^#/ ) {
		print FCF $_ . "\n" if ($i == 0);
	    } else {
		print FCF $_ . "\n";
	    }
	}
	close(FHPF);
	close(FPF);
	close(FDRF);
	close(FFCF);
	unlink $fragHairpinsFile;
	unlink $fragProductFile;
	unlink $fragDistinctReadsFile;
	unlink $fragFilteredCandidatesFile;
    }
    close(HPF);
    close(PF);
    close(DRF);
    close(FCF);
}

sub breakupReadRegionsFile {
    my($readRegions,$numFiles) = @_;
    my $autoOverwrite = 1;  #if the file exists the program automatically writes over it.
    my %readRegionsFiles;
    my($fileBase) = $readRegions =~ /([^\/]*)\.txt/;
    open(RR,$readRegions) or die "failed to open $readRegions\n";
    my $header = <RR>;
    chomp($header);
    my $wcOut = `wc -l $readRegions` or die "failed to run command \"wc\" on command line";
    chomp($wcOut);
    my($numLines) = $wcOut =~ /(^\d+)/;
    my $curFileNum = 0;
    my $lastFileNum = $numFiles - 1;
    my $nextNum = 0;
    my $fairShare = $numLines / $numFiles;
    my $outCount = -1;
    #splititng up the read regions file
    while(<RR>) {
	chomp;
	my $nextLine = $_;
	if ($outCount == -1) {
	    my $nextFile = $fileBase . "_$nextNum.txt";
	    $nextNum++;
	    $outCount++;
	    unless ($autoOverwrite) {
		while (-e $nextFile) {
		    $nextFile = $fileBase . "_$nextNum.txt";
		    $nextNum++;
		}
	    }
	    $readRegionsFiles{$curFileNum} = $nextFile;
	    open(NRR, ">$nextFile") or die "cannot open $nextFile for writing\n";
	    print NRR "$header\n";
	}
	print NRR "$nextLine\n";
	$outCount++;
	if (($outCount >= $fairShare) && ($curFileNum != $lastFileNum)) {
	    close(NRR) or die "failed to close an output file\n";
	    $outCount = -1;
	    $curFileNum++;
	}
	die "more files were created than there were threads\n" if ($curFileNum > $lastFileNum);
    }
    close(NRR);
    close(RR);
    return \%readRegionsFiles;
}

#note: the following function causes a segfault
sub processReadRegionsMultiThreaded {
    my($bamList,$readRegionsFile,$genomeDir,$chromLengths,$parameters) = @_;
    my $numThreads = $parameters->{numThreads};
    my $readRegionsFiles = breakupReadRegionsFile($readRegionsFile,$numThreads);
    my @threads;
    #creating threads
    for (my $i = 0; $i < $numThreads; $i++) {
	my $threadReadRegionsFile = $readRegionsFiles->{$i};
	my($threadRRFileBase) = $threadReadRegionsFile =~ /([^\/]*)\.txt/;
	my $newParams = {"outputPrefix" => "$threadRRFileBase"};
	my $thr = threads->create('processReadRegionsWithNewParams',$bamList,$threadReadRegionsFile,$genomeDir,$chromLengths,$parameters,$newParams);
	push(@threads,$thr);
    }
    #joining threads
    my @running = threads->list(threads::running);
    while (@running != 0) {
	foreach my $thr (@threads) {
	   $thr->join if ($thr->is_joinable()); 
	}
	@running = threads->list(threads::running);
    }
}

sub processReadRegionsMultiProcess {
    my($bamList,$readRegionsFile,$genomeDir,$chromLengths,$parameters) = @_;
    my $numProcesses = $parameters->{numProcesses};
    my $readRegionsFiles = breakupReadRegionsFile($readRegionsFile,$numProcesses);
    my @children;
    #creating threads
    for (my $i = 0; $i < $numProcesses; $i++) {
	my $pid = fork();
	if( $pid )
	{
	    #If $pid is non zero, then the parent is running
	    print "PID $pid forked ($i)\n";
	    push(@children, $pid);
	}
	else
	{
	    # Else we are a child process ($pid == 0)
	    my $childReadRegionsFile = $readRegionsFiles->{$i};
	    my($childRRFileBase) = $childReadRegionsFile =~ /([^\/]*)\.txt/;
	    my $newParams = {"outputPrefix" => "$childRRFileBase"};
	    processReadRegionsWithNewParams($bamList,$childReadRegionsFile,$genomeDir,$chromLengths,$parameters,$newParams);
	}
    }
    #waiting for children
    foreach my $n (@children){
	my $pid = waitpid($n,0); # waitpid returns the pid that finished
	print $pid . " done\n";
    }
}

sub processReadRegionsWithNewParams {
    my($bamList,$readRegionsFile,$genomeDir,$chromLengths,$parameters,$newParams) = @_;
    foreach my $paramName (keys %{$parameters}) {
	$newParams->{$paramName} = $parameters->{$paramName} unless ($newParams->{$paramName});
    }
    if ($newParams->{outputPrefix}) {
	createOutputFileParameters($newParams);
    }
    processReadRegions($bamList,$readRegionsFile,$genomeDir,$chromLengths,$newParams);
}

sub getMaxPredictedProductScores {
    my($predProductFile) = @_;
    my %maxPredProductScores;
    open(PPF,$predProductFile) or die "failed to open $predProductFile for reading\n";
    while (<PPF>) {
	chomp;
	my($id_location, $predictedValue, $RFScore) = split("\t", $_);
	if ($predictedValue == 1) {
	    my($id,$location) = $id_location =~ /(.*?)\_(.*)/;
	    $maxPredProductScores{$id} = $RFScore;
	}
    }
    close(PPF);
    return \%maxPredProductScores;
}


sub getReadLocationVariance {
    my($products,$parameters) = @_;
    my(%productLocationVariance,%productReads);
    my($productReads,$totalReads) = (0,0);
    foreach my $productInfo (@{$products}) {
	my($id,$productReads,$adjMaxProductCount,$adjProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$productInfo};
	my($total,$sumOfSquares,$numCounts) = (0,0,0);
	foreach my $drInfo (@{$productReads}) {
	    my($relStart,$relStop,$dStart,$readCount,$hitCount,$libraryCounts,$seq,$adjustedSeqCount) = @{$drInfo};
	    for (my $i = 0; $i < $readCount; $i++) {
		$total += ($maxRelStart - $relStart);
		$sumOfSquares += ($maxRelStart - $relStart)**2;
		$numCounts += 1;
	    }
	    $productReads{$id} += $readCount;
	    $totalReads += $readCount;
	}
	my $vNumerator = $numCounts * $sumOfSquares - $total**2;
	my $vDenominator = $numCounts * ($numCounts - 1);
	$productLocationVariance{$id} = ($vDenominator) ? $vNumerator / $vDenominator : 0;
    }
    my($totalVariance, $weightedVariance) = (0,0);
    foreach my $id (keys %productLocationVariance) {
	my $variance = $productLocationVariance{$id};
	$totalVariance += $variance;
	$weightedVariance += $variance * $productReads{$id} / $totalReads;
    }
    my $avgProdVariance = $totalVariance / $totalReads;
    return($avgProdVariance,$weightedVariance);
}


sub getReadDataProductVariance {
    my($products,$parameters) = @_;
    my(%readCountHash,%productVariance,%productWeight);
    my($totalHPReadCount,$averageProductVariance,$weightedAverageProductVariance) = (0,0,0);
    #accumulating reads into a hash keyed by start position
    foreach my $productInfo (@{$products}) {
	my($id,$productReads,$adjMaxProductCount,$adjProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$productInfo};
	foreach my $drInfo (@{$productReads}) {
	    my($relStart,$relStop,$dStart,$readCount,$hitCount,$libraryCounts,$seq,$adjustedSeqCount) = @{$drInfo};
	    #print "$relStart x $relStop $adjustedSeqCount added\n";
	    $readCountHash{$id}{$relStart} += $adjustedSeqCount;
	    $totalHPReadCount += $adjustedSeqCount;
	    #print "totalHPReadCount now $totalHPReadCount\n";
	}
    }
    #geting variances of reads for each product
    foreach my $productInfo (@{$products}) {
	my($id,$productReads,$adjMaxProductCount,$adjProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$productInfo};
	my($total,$sumOfSquares,$numCounts) = (0,0,0);
	foreach my $relStart (keys %{$readCountHash{$id}}) {
	    my $count = $readCountHash{$id}{$relStart};
	    $total += $count;
	    $sumOfSquares += $count**2;
	    #print "total=$total\tsquared=$sumOfSquares\n";
	    $numCounts += 1;
	}
	my $vNumerator = $numCounts * $sumOfSquares - $total**2;
	#print "vNumerator = $vNumerator\n";
	my $vDenominator = $numCounts * ($numCounts - 1);
	#print "vDenominator = $vDenominator\n";
	$productVariance{$id} = ($vDenominator) ? $vNumerator / $vDenominator : 0;
	#print "productVariance = " . $productVariance{$id} . "\n";
	$productWeight{$id} = $total / $totalHPReadCount;
	#print "productWeight = " . $productWeight{$id} . "\n";
    }
    #computing the average and weighted average for the product variances
    my($sumOfVariances,$numProducts) = (0,0);
    foreach my $productInfo (@{$products}) {
	my($id,$productReads,$adjMaxProductCount,$adjProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$productInfo};
	$sumOfVariances += $productVariance{$id};
	#print "sumOfVariances = " . $sumOfVariances . "\n";
	$weightedAverageProductVariance += $productVariance{$id} * $productWeight{$id};
	#print "weightedAverageProductVariance = " . $weightedAverageProductVariance . "\n";
	$numProducts += 1;
    }
    $averageProductVariance = $sumOfVariances / $numProducts;
    #print "averageProductVariance = "  . $averageProductVariance .  "\n";
    return($averageProductVariance,$weightedAverageProductVariance);
}

sub getMaxAndDuplexProducts {
    my($senseProducts,$fold,$readId) = @_;
    my @sortedProducts = sort {$b->[3] <=> $a->[3]} @{$senseProducts};
    my $maxProduct = $sortedProducts[0];
    my($mpSide,$mpType,$mpProductList,$mpAdjProductCount,$mpAdjMaxProductCount,$mpRelStart,$mpRelStop,$mpOffset,$mpGStart,$mpProductStrand) = @{$maxProduct};  
    #print "$readId\n$fold\n";
    #print " " x $mpRelStart . substr($fold,$mpRelStart,$mpRelStop-$mpRelStart+1) . "\n";
    #print "total length = " . length($fold) . "\n"; 
    #print "$mpRelStart\t$mpRelStop\n";
    my $map = pairMap($fold);
    my($j,$k,$jFound) = (-1,-1,0);
    for(my $i = $mpRelStart; $i <= $mpRelStop; $i++) {
	#print "\$i = $i\n";
	unless ($jFound) {
	    my $duplexPair = ${$map}[$i];
	    #my $duplexPair = (defined(${$map}[$i])) ? ${$map}[$i] : '.';
	    #print "$i\t$duplexPair\n";
	    unless ($duplexPair eq '.') {
		$j = $duplexPair;
		$k = $j;
		$jFound = 1;
		#print "\$j = $j\t\$k = $k\n";
	    }
	} else {
	    #my $duplexPair = ${$map}[$i]; 
	    #sometimes the major product is not part of the readregion and when the sequence is extended in processReadRegions it does not cover the
	    #entire major product.  In the following line '.' is assigned to $duplexPair for parts of the major product not on the hairpin.  In the 
	    #future it may be better to extend the sequence further in order to fully cover the major product.
	    my $duplexPair = (defined(${$map}[$i])) ? ${$map}[$i] : '.';
	    #print "$i\t$duplexPair\n";
	    unless ($duplexPair eq '.') {
		$k = $duplexPair;
		#print "$i\t$k\n";
	    }
	}
    }
    #print "\n";
    my $duplexProduct = 0;
    my $duplexProductCount = 0;
    my $duplexProductOverlap = 0;
    unless ($j == -1) {  #j=-1 if major product was in a loop and there was no duplex
	my $duplexLeft = ($j < $k) ? $j : $k;
	my $duplexRight = ($k <= $j) ? $j : $k;
	#print "duplex of $mpRelStart $mpRelStop ($mpType)  mapped to $duplexLeft $duplexRight\n";
	foreach my $product (@sortedProducts) {
	    my($side,$type,$productList,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};  
	    unless ($product == $maxProduct) {
		my $overlap = getOverlap($relStart,$relStop,$duplexLeft,$duplexRight);
		if ($overlap > 0 && $overlap == $duplexProductOverlap && $duplexProductCount < $adjProductCount) {
		    #we will choose the more abundant duplex if two overlap the duplex region the same amount
		    $duplexProductCount = $adjProductCount;
		    $duplexProduct = $product;
		} elsif ($overlap > $duplexProductOverlap) {
		    $duplexProductOverlap = $overlap;
		    $duplexProduct = $product;
		    $duplexProductCount = $adjProductCount;
		}
	    }
	}
    }
    #duplex product may be null when returned
    return($maxProduct,$duplexProduct,$mpAdjProductCount,$duplexProductCount,$duplexProductOverlap);
}

sub getMaxProductLoopDistance {
    my($center,$maxProduct,$duplexProduct) = @_;
    my($leftCenter,$rightCenter) = @{$center};
    my $loopSize = $rightCenter - $leftCenter + 1;
    my($side,$type,$productList,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$maxProduct};
    my($loopDistance,$dupLoopDistance);  
    my $overlap = getOverlap($leftCenter,$rightCenter,$relStart,$relStop);
    unless($overlap) {
	if ($relStop <= $leftCenter) {
	    $loopDistance = $leftCenter - $relStop;
	} elsif ($relStart >= $rightCenter) {
	    $loopDistance = $relStart - $rightCenter;
	} else {
	    die "error: cannot get loop distance in getMaxProductLoopDistance ($side\t$type\t$relStart-$relStop  Loop:$leftCenter-$rightCenter  gStart=$gStart)\n";
	}
    } else {
	#products with negative distances from loops may be splits or loops and could help define a decision boundry in the model
	$loopDistance = -$overlap;
    }
    if ($duplexProduct) {  #if there is a duplex Product
	my($dupSide,$dupType,$dupProductList,$dupAdjProductCount,$dupAdjMaxProductCount,$dupRelStart,$dupRelStop,$dupOffset,$dupGStart,$dupProductStrand) = @{$duplexProduct};	
	my $dupOverlap = getOverlap($leftCenter,$rightCenter,$dupRelStart,$dupRelStop);
	unless($dupOverlap) {
	    if ($dupRelStop <= $leftCenter) {
		$dupLoopDistance = $leftCenter - $dupRelStop;
	    } elsif ($dupRelStart >= $rightCenter) {
		$dupLoopDistance = $dupRelStart - $rightCenter;
	    } else {
		die "error: cannot get loop distance in getMaxProductLoopDistance ($side\t$type\t$relStart-$relStop  Loop:$leftCenter-$rightCenter  gStart=$gStart)\n";
	    }
	} else {
	    #products with negative distances from loops may be splits and could help define a decision boundry in the model
	    $dupLoopDistance = -$dupOverlap;
	}
    } else {
	$dupLoopDistance = 0;
    }
    return($loopDistance,$dupLoopDistance,$loopSize);
}

sub getMaxProductOverlapScores {
    my($senseProducts,$maxProduct,$duplexProduct) = @_;
    my($mpSide,$mpType,$mpProductList,$mpAdjProductCount,$mpAdjMaxProductCount,$mpRelStart,$mpRelStop,$mpOffset,$mpGStart,$mpProductStrand) = @{$maxProduct};
    my($dupSide,$dupType,$dupProductList,$dupAdjProductCount,$dupAdjMaxProductCount,$dupRelStart,$dupRelStop,$dupOffset,$dupGStart,$dupProductStrand);
    if ($duplexProduct) {
	($dupSide,$dupType,$dupProductList,$dupAdjProductCount,$dupAdjMaxProductCount,$dupRelStart,$dupRelStop,$dupOffset,$dupGStart,$dupProductStrand) = @{$duplexProduct};
    }
    my($mpOverlapScore,$dupOverlapScore) = (0,0);
    my(@mpOverlaps,@dupOverlaps);
    my($totalMPOverlap,$totalDupOverlap) = (0,0);
    foreach my $product (@{$senseProducts}) {
	my($side,$type,$productList,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
	unless ($type =~ /long/) {  #long mirs tend to overlap so we exclude them
	    unless ($product == $maxProduct) {
		my $overlap = getOverlap($relStart,$relStop,$mpRelStart,$mpRelStop);
		if ($overlap) {
		    my $overlapScore = $overlap * $adjProductCount * $mpAdjProductCount;
		    $mpOverlapScore = $mpOverlapScore + $overlapScore;
		    push(@mpOverlaps,$product);
		    $totalMPOverlap = $totalMPOverlap + $overlap;
		}
	    }
	    if ($duplexProduct) {
		unless ($product == $duplexProduct) {
		    my $overlap = getOverlap($relStart,$relStop,$dupRelStart,$dupRelStop);
		    if ($overlap) {
			my $overlapScore = $overlap * $adjProductCount * $dupAdjProductCount;
			$dupOverlapScore = $dupOverlapScore + $overlapScore;
			push(@dupOverlaps,$product);
			$totalDupOverlap = $totalDupOverlap + $overlap;
		    }		
		}
	    }
	}
    }
    my $wghtMPOverlapIntensity = 0;
    #my $numMPOverlaps = @mpOverlaps;
    foreach my $product (@mpOverlaps) {
	my($side,$type,$productList,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
	my $overlap = getOverlap($relStart,$relStop,$mpRelStart,$mpRelStop);
	$wghtMPOverlapIntensity = $wghtMPOverlapIntensity + $adjProductCount * $overlap / $totalMPOverlap;
    }
    my $wghtDupOverlapIntensity = 0;
    #my $numDupOverlaps = @dupOverlaps;
    foreach my $product (@dupOverlaps) {
	my($side,$type,$productList,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
	my $overlap = getOverlap($relStart,$relStop,$dupRelStart,$dupRelStop);
	$wghtDupOverlapIntensity = $wghtDupOverlapIntensity + $adjProductCount * $overlap / $totalDupOverlap;
    }
    return($mpOverlapScore,$dupOverlapScore,$wghtMPOverlapIntensity,$wghtDupOverlapIntensity);
}


sub getMaxProductOverlapScores_old_ForTesting {
    #this function was created to score overlaps in a lot of different ways in an attempt to figure out which scoring method worked best
    my($senseProducts,$maxProduct,$duplexProduct) = @_;
    my($mpSide,$mpType,$mpProductList,$mpAdjProductCount,$mpAdjMaxProductCount,$mpRelStart,$mpRelStop,$mpOffset,$mpGStart,$mpProductStrand) = @{$maxProduct};
    my($dupSide,$dupType,$dupProductList,$dupAdjProductCount,$dupAdjMaxProductCount,$dupRelStart,$dupRelStop,$dupOffset,$dupGStart,$dupProductStrand);
    if ($duplexProduct) {
	($dupSide,$dupType,$dupProductList,$dupAdjProductCount,$dupAdjMaxProductCount,$dupRelStart,$dupRelStop,$dupOffset,$dupGStart,$dupProductStrand) = @{$duplexProduct};
    }
    #at most only two of the following scores will be kept after the best scoring method is found
    my($mpOverlapScore,$dupOverlapScore,$mpOverlapScoreV2,$dupOverlapScoreV2,$mpOverlapScoreV3,$dupOverlapScoreV3,$mpOverlapScoreV4,$dupOverlapScoreV4) = (0,0,0,0,0,0,0,0);
    my(@mpOverlaps,@dupOverlaps);
    my($totalMPOverlap,$totalDupOverlap,$totalMPOverlapCount,$totalDupOverlapCount) = (0,0,0,0);
    foreach my $product (@{$senseProducts}) {
	my($side,$type,$productList,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
	unless ($type =~ /long/) {  #long mirs tend to overlap so we exclude them
	    unless ($product == $maxProduct) {
		my $overlap = getOverlap($relStart,$relStop,$mpRelStart,$mpRelStop);
		if ($overlap) {
		    my $overlapScore = $overlap * $adjProductCount / ($mpAdjProductCount + $adjProductCount);
		    my $overlapScoreV2 = $overlap * $adjProductCount;
		    my $overlapScoreV3 = $overlap * $adjProductCount / $mpAdjProductCount;
		    my $overlapScoreV4 = $overlap * $adjProductCount * $mpAdjProductCount;
		    $mpOverlapScore = $mpOverlapScore + $overlapScore;
		    $mpOverlapScoreV2 = $mpOverlapScoreV2 + $overlapScoreV2;
		    $mpOverlapScoreV3 = $mpOverlapScoreV3 + $overlapScoreV3;
		    $mpOverlapScoreV4 = $mpOverlapScoreV4 + $overlapScoreV4;
		    push(@mpOverlaps,$product);
		    $totalMPOverlap = $totalMPOverlap + $overlap;
		    $totalMPOverlapCount = $totalMPOverlapCount + $adjProductCount;
		}
	    }
	    if ($duplexProduct) {
		unless ($product == $duplexProduct) {
		    my $overlap = getOverlap($relStart,$relStop,$dupRelStart,$dupRelStop);
		    if ($overlap) {
			my $overlapScore = $overlap * $adjProductCount / ($dupAdjProductCount + $adjProductCount);
			$dupOverlapScore = $dupOverlapScore + $overlapScore;
			if ($dupAdjProductCount >= $adjProductCount) {
			    my $overlapScoreV2 = $overlap * $adjProductCount;
			    my $overlapScoreV3 = $overlap * $adjProductCount / $dupAdjProductCount;
			    $dupOverlapScoreV2 = $dupOverlapScoreV2 + $overlapScoreV2;
			    $dupOverlapScoreV3 = $dupOverlapScoreV3 + $overlapScoreV3;
			} else {
			    my $overlapScoreV2 = $overlap * $dupAdjProductCount;
			    my $overlapScoreV3 = $overlap * $dupAdjProductCount / $adjProductCount;
			    $dupOverlapScoreV2 = $dupOverlapScoreV2 + $overlapScoreV2;
			    $dupOverlapScoreV3 = $dupOverlapScoreV3 + $overlapScoreV3;
			}
			my $overlapScoreV4 = $overlap * $adjProductCount * $dupAdjProductCount;
			$dupOverlapScoreV4 = $dupOverlapScoreV4 + $overlapScoreV4;
			push(@dupOverlaps,$product);
			$totalDupOverlap = $totalDupOverlap + $overlap;
			$totalDupOverlapCount = $totalDupOverlapCount + $adjProductCount;
		    }		
		}
	    }
	}
    }
    my($wghtMPOverlap,$wghtMPOverlapSize) = (0,0);
    my $numMPOverlaps = @mpOverlaps;
    foreach my $product (@mpOverlaps) {
	my($side,$type,$productList,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
	my $overlap = getOverlap($relStart,$relStop,$mpRelStart,$mpRelStop);
	$wghtMPOverlap = $wghtMPOverlap + $overlap * $adjProductCount / $totalMPOverlapCount;
	$wghtMPOverlapSize = $wghtMPOverlapSize + $adjProductCount * $overlap / $totalMPOverlap;
    }
    my($wghtDupOverlap,$wghtDupOverlapSize) = (0,0);
    my $numDupOverlaps = @dupOverlaps;
    foreach my $product (@dupOverlaps) {
	my($side,$type,$productList,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
	my $overlap = getOverlap($relStart,$relStop,$dupRelStart,$dupRelStop);
	$wghtDupOverlap = $wghtDupOverlap + $overlap * $adjProductCount / $totalDupOverlapCount;
	$wghtDupOverlapSize = $wghtDupOverlapSize + $adjProductCount * $overlap / $totalDupOverlap;
    }
    return($mpOverlapScore,$dupOverlapScore,$mpOverlapScoreV2,$dupOverlapScoreV2,$mpOverlapScoreV3,$dupOverlapScoreV3,$mpOverlapScoreV4,$dupOverlapScoreV4,$numMPOverlaps,$numDupOverlaps,$totalMPOverlap,$totalDupOverlap,$totalMPOverlapCount,$totalDupOverlapCount,$wghtMPOverlap,$wghtDupOverlap,$wghtMPOverlapSize,$wghtDupOverlapSize);
}

sub getHairpinRPM {
    my($readCount,$librarySizes) = @_;
    my $libSizeTotal = 0;
    foreach my $sample (keys %{$librarySizes}) {
	$libSizeTotal += $librarySizes->{$sample};
    }
    my $RPM = $readCount / ($libSizeTotal / 1e6);
    return $RPM;
}

sub processReadRegions {
    my($bamList,$readRegionsFile,$genomeDir,$chromLengths,$parameters) = @_;
    my $hairpinVectorFile = $parameters->{outputPrefix} . "_hairpinScores";
    my $librarySizesFile = $parameters->{librarySizes} or die "no librarySizes entry found in config file\n";
    my $librarySizes = loadLibrarySizes($librarySizesFile);
    my $testTime = 0;
    my $testMemory = 0;
    my $testReadCounts = 0;
    my $readsCount = 0; #number of mirBase hairpins with reads
    my $noReadsCount = 0;  #number of mirBase hairpins with no reads
    my $readsLessThanCountMinLocus = 0; #number of hairpins with fewer than count min locus reads
    my $centersReadsCount = 0; #number of centers on mirBase hairpins with reads
    my $noCenterCount = 0;  #number of hairpins without centers (hopefully this will stay 0)
    my $mu = Memory::Usage->new() if ($testMemory);
    #my $tk = TimeKeeper->new() if ($testTime);
    my @sampleList;
    my $minLength = $parameters->{lengthMin} or die "minLength: not loaded.\n";
    my $filteredCandidatesFile = $parameters->{filteredCandidatesFile};
    my $hairpinsFile = $parameters->{hairpinsFile};
    my $distinctReadsFile = $parameters->{distinctReadsFile};
    my $productFile = $parameters->{productFile};
    my $featureVectorFile = $parameters->{featureVectorFile};
    my $maxPredictedProductScores = (-e $parameters->{predProductClasses}) ? getMaxPredictedProductScores($parameters->{predProductClasses}) : 0;
    initializeRNAFold();
    open(HL,">".$hairpinsFile) or die "failed to open $hairpinsFile for writing\n";
    open(HDR,">".$distinctReadsFile) or die "failed to open $distinctReadsFile for writing\n";
    open(HPL,">".$productFile) or die "failed to open $productFile for writing\n";
    print HL "#tag\tchrom\tstart\tstop\tstrand\tleftCenter\trightCenter\ttotalSense\ttotalAntisense\tmfe\tsequence\tfold\t";
    print HL "aapd\ttapd\turf\tahc\tmpfh\tpbp\tsameShifted\tbothShift\tsplitProductAbundance\toverlapProductAbundance\t";
    print HL "aa\tac\tag\tat\tca\tcc\tcg\tct\tga\tgc\tgg\tgt\tta\ttc\ttg\ttt\tduplexEnergy\tGCcontent\tfoldDupCmp\texactFoldDupCmp\tdupPBP\tdupLoopLength";
#    print HL "aa\tac\tag\tat\tca\tcc\tcg\tct\tga\tgc\tgg\tgt\tta\ttc\ttg\ttt\tduplexEnergy\tGCcontent\tzScore\tfoldDupCmp\texactFoldDupCmp\tdupPBP\tdupLoopLength";
    print HL "\tAPV\twAPV";
    print HL "\tARV\twARV";
    print HL "\tmpCount\tdupCount\tdupOverlap\tmpLoopDistance\tdupLoopDistance\tloopSize\tmpOverlapScore\tdupOverlapScore\twghtMPOverlapIntensity\twghtDupOverlapIntensity\tOPA2\ttotalOverlapAmount\taverageOverlapAmount\ttotalRelativeOverlapAmount\taverageRelativeOverlapAmount\tmaxOverlap\tinnerLoopGapCount\tsenseRPM\tantiSenseRPM";
    if ($maxPredictedProductScores) {
	print HL "\tRFProductAvg\n";
    } else {
	print HL "\n";
    }
    print HDR "#tag\tstart\tstop\tstrand\ttotal reads\toffsetStart\trelative Start\trelative Stop\tsequence";
    print HPL "#tag\tside type\ttype\ttotal reads\ttotal most abundant\tstart\tstop\tstrand\tsequence";
    foreach my $bamListElement (@{$bamList}) {
	my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamListElement};
	print HDR "\t$sample";
        print HPL "\t$sample";
    }
    print HDR "\n";
    print HPL "\n";
    
    #$tk->start("TotalTime") if ($testTime);
    
    open(RRF, $readRegionsFile) or die "failed to open $readRegionsFile\n";
    my $prevChrom;
    my $genome;
    while (<RRF>) {
        chomp;
	unless(/\#/) {	  
	    my($chrom,$start,$stop,$id,$score,$strand) = split(/\t/,$_);
	    #print "rrlocation = $chrom:$start..$stop..$strand\n";
	    $start++; #converting to one based
	    #the chroms in RRF should be sorted at this point so the next line of code will only be executed once per chrom
	    my $maxProductScore;
	    if ($maxPredictedProductScores) {
		$maxProductScore = ($maxPredictedProductScores->{$id}) ? $maxPredictedProductScores->{$id} : 0;
	    }
	    if ($prevChrom) {
		$genome = loadGenome("$genomeDir/$chrom.fa") unless ($prevChrom eq $chrom);	
	    } else {
		$genome = loadGenome("$genomeDir/$chrom.fa");	
	    }
	    $prevChrom = $chrom;
	    my $location = extendSequenceToFillBuffer($chrom,$start,$stop,$strand,$chromLengths->{$chrom},$parameters);
	    #print "extended location = $location\n";
	    ($chrom,$start,$stop,$strand) = parseLocation($location);
	    my $sequence = getSequence($location,$genome);
    	    #$tk->start("retrieveReadData") if ($testTime);
	    my($distinctReads,$uniqueReadCount,$totalDistinctReads,$adjustedReadCount) = retrieveReadData($location,$bamList,$parameters);
	    #$tk->end("retrieveReadData") if ($testTime);
	    if ($distinctReads->{$strand}) {  #neccessary if testing read reagions with known mirs and no real read exists in that read region
		if ($adjustedReadCount->{$strand} < $parameters->{countMinLocus}) {
		    $readsLessThanCountMinLocus++;
		    print "$id has fewer reads than count min locus\n" if ($testReadCounts);
		}
		$readsCount++;
		#$tk->start("buildProducts") if ($testTime);
		my($senseProducts) = buildProducts($location,$distinctReads->{$strand},$strand,$parameters);
		my($averageProductVariance,$weightedAverageProductVariance) = getReadDataProductVariance($senseProducts,$parameters);
		my($averageReadVariance,$weightedAverageReadVariance) = getReadLocationVariance($senseProducts,$parameters);
		#print "$averageProductVariance\t$weightedAverageProductVariance\n";
		#$tk->end("buildProducts") if ($testTime);
		#$tk->start("getFold") if ($testTime);
		#here productDuplexCoords is the relative coordinates of the product after duplexing.  seqDuplexCoords is the actual produts duplex.
		my($fold,$mfe,$maxProductCoords,$productDuplexCoords,
		   $seqDuplexCoords,$dupStructure,$duplexEnergy) = getFold($chrom,$sequence,$senseProducts,$id);
		my($foldDuplexCmp,$exactFoldDuplexCmp) = compareFoldAndDuplexStructs($fold,$maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure,$id);
		my($prodDupPBP,$mpDuplexLoopLength) = getProductDuplexPBP($fold,$maxProductCoords);
		#$tk->end("getFold") if ($testTime);
		#$tk->start("getMergedHairpinCenters") if ($testTime);
		my @centers = getMergedHairpinCenters($fold,$parameters);
		#$tk->end("getMergedHairpinCenters") if ($testTime);
		$centersReadsCount += scalar(@centers);
		unless (scalar(@centers)) {
		    $noCenterCount++;
		    print "$id has no centers\n" if ($testReadCounts);
		    print "$fold\n" if ($testReadCounts);
		}
		my $COUNT= 0;
		foreach my $center (@centers) {
		    my $asciiVal = ord('a') + $COUNT;
		    my $label = chr($asciiVal);
		    my $newId = $id . $label;
		    $COUNT++;
		    #$tk->start("getBasePairs") if ($testTime);
		    my $basePairs = getBasePairs($center,$fold);
		    #$tk->end("getBasePairs") if ($testTime);
		    #$tk->start("getMaxHairpinLength") if ($testTime);
		    my $hairpinLength = getMaxHairpinLength($fold,$center);
		    #$tk->end("getMaxHairpinLength") if ($testTime);
		    #$tk->start("getFullHairpinWindow") if ($testTime);
		    my($hStart,$hStop) = getFullHairpinWindow($fold,$center);
		    #$tk->end("getFullHairpinWindow") if ($testTime);
		    if($hairpinLength >= $minLength) {
			#$tk->start("getProductInfo") if ($testTime);
			my $senseProducts=getProductInfo($center,$fold,$basePairs,$location,$senseProducts,$parameters);
			#$tk->end("getProductInfo") if ($testTime);
			#$tk->start("rebuildProducts") if ($testTime);
			$senseProducts=rebuildProducts($center,$fold,$basePairs,$senseProducts,$parameters);
			#$tk->end("rebuildProducts") if ($testTime);
			$senseProducts=addProductSize($senseProducts,$parameters);
#			if(goodProducts($senseProducts,$parameters)) {
			my $revStrand = revStrand($strand);
			#$tk->start("buildProducts") if ($testTime);
			my($antiSenseProducts) = buildProducts($location,$distinctReads->{$revStrand},$revStrand,$parameters);
			#$tk->end("buildProducts") if ($testTime);
			#$tk->start("getProductInfo") if ($testTime);
			$antiSenseProducts = getProductInfo($center,$fold,$basePairs,$location,$antiSenseProducts,$parameters);
			#$tk->end("getProductInfo") if ($testTime);
			#$tk->start("rebuildProducts") if ($testTime);
			$antiSenseProducts = rebuildProducts($center,$fold,$basePairs,$antiSenseProducts,$parameters);
			#$tk->end("rebuildProducts") if ($testTime);
			$antiSenseProducts=addProductSize($antiSenseProducts,$parameters);
			#$tk->start("getAdjTotalProductReads") if ($testTime);
			my $adjTotalProductReads = getAdjTotalProductReads($senseProducts);
			#$tk->end("getAdjTotalProductReads") if ($testTime);
			#$tk->start("getAdjTotalProductReads") if ($testTime);
			my $adjTotalRevProductReads = getAdjTotalProductReads($antiSenseProducts);
			#$tk->end("getAdjTotalProductReads") if ($testTime);
			my($PASS,$REASON) = plausibleReads($adjTotalProductReads,$adjTotalRevProductReads,$senseProducts,$antiSenseProducts,
							   $parameters);
#			    if($PASS) {
			my($maxProduct,$mpDuplex,$maxProductCount,$mpDuplexCount,$duplexOverlap) = getMaxAndDuplexProducts($senseProducts,$fold,$id);
			#both mpLoopDiatance and dupLoopDistance can be negative if overlapping the loop
			my($mpLoopDistance,$dupLoopDistance,$loopSize) = getMaxProductLoopDistance($center,$maxProduct,$mpDuplex);
			my($mpOverlapScore,$dupOverlapScore,$wghtMPOverlapIntensity,$wghtDupOverlapIntensity) = getMaxProductOverlapScores($senseProducts,$maxProduct,$mpDuplex);
			#$tk->start("getReverseProductDisplacement") if ($testTime);
			my($tpd,$totalRP)=newGetReverseProductDisplacement($senseProducts,
									   $antiSenseProducts,
									   $adjTotalProductReads,
									   $parameters,$newId);
			#$tk->end("getReverseProductDisplacement") if ($testTime);
			my $apd = $totalRP ? $tpd/$totalRP : 0.0;
			my $urf = $adjustedReadCount->{$strand} ? $uniqueReadCount->{$strand} / $adjustedReadCount->{$strand} : 0;
			#$tk->start("computeMaxProdHitCount") if ($testTime);
			my $ahc = computeMaxProdHitCount($senseProducts,$location,$distinctReads,$parameters);
			#$tk->end("computeMaxProdHitCount") if ($testTime);
			#$tk->start("computeMaxProdFivePrimeHet") if ($testTime);
			my $afh = computeMaxProdFivePrimeHet($senseProducts,$parameters);
			#$tk->end("computeMaxProdFivePrimeHet") if ($testTime);
			#$tk->start("computeProductBasePairing") if ($testTime);
			my $pbp = computeProductBasePairing($center,$senseProducts,$basePairs,$parameters);
			#$tk->end("computeProductBasePairing") if ($testTime);
			#$tk->start("computeMaxSameShift") if ($testTime);
			my $sameShift = computeMaxSameShift($location,$distinctReads->{$strand},$senseProducts,$parameters);
			#$tk->end("computeMaxSameShift") if ($testTime);
			#$tk->start("computeMaxBothShift") if ($testTime);
			my $bothShift = computeMaxBothShift($basePairs,$location,$distinctReads->{$strand},$senseProducts,$parameters,$newId);
			#$tk->end("computeMaxBothShift") if ($testTime);
#			my($zScore,$consensusMFE) = getRNAzData($sequence);
			my $innerLoopGapCount = getInnerLoopGapCount($center,$fold);
			my $splitProductAbundance = computeSplitProductAbundance($senseProducts,$parameters);
			my $overlapProductAbundance = computeOverlapProductAbundance($senseProducts,$parameters);
			my $opa2 = computeOverlapProductAbundance2($senseProducts,$parameters);
			my ($totalOverlapAmount,$averageOverlapAmount,
			 $totalRelativeOverlapAmount,$averageRelativeOverlapAmount) = computeOverlapAmounts($senseProducts,$parameters);
			my $maxOverlap = getMaximumOverlap($senseProducts,$parameters);
			my($aa,$ac,$ag,$at,$ca,$cc,$cg,$ct,$ga,$gc,$gg,$gt,$ta,$tc,$tg,$tt) = getDinucleotideFrequencies($sequence);
			my $sequenceLength = length($sequence);
			my $gcContent = getGCcontent($sequence,$sequenceLength);
			# total is the read count for the whole hairpin
			my $totalSense = 0;
			foreach my $product (@{$senseProducts}) {
			    my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
			       $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
			    my $length = $relStop - $relStart + 1;
			    my $productSequence = substr($sequence,$relStart,$length);
			    $totalSense += $adjProdCount;
			    #$tk->start("getAdjTotalLibraryCounts") if ($testTime);
			    my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($prodList);
			    #$tk->end("getAdjTotalLibraryCounts") if ($testTime);
			    print HPL "$newId\t$side\t$newType\t$adjProdCount\t";
			    print HPL "$adjMaxProdCount\t$relStart\t$relStop\t$productStrand\t";
			    print HPL "$productSequence";
			    foreach my $bamElement (@{$bamList}) {
				my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
				printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
			    }
			    print HPL "\n";
			    foreach my $read (@{$prodList}) {
				my($relStart,$relStop,$offset,$gStart,
				   $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
				my $gStop = $gStart + ($relStop - $relStart);
				my $adjCount = $count * $adjSeqCount;
				print HDR "$newId\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
				print HDR "$offset\t$relStart\t$relStop\t$seq";
				foreach my $bamElement (@{$bamList}) {
					    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
					    printf(HDR "\t%.3f",$adjLibraryCounts->{$sample});
				}
				print HDR "\n";
			    }				
			}
			#############################
			# ANTISENSE PRODUCTS
			my $totalAntisense = 0;
			foreach my $product (@{$antiSenseProducts}) {
			    my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
			       $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
			    $newType = "as-".$newType;
			    $totalAntisense += $adjProdCount;
			    my $length = $relStop - $relStart + 1;
			    my $productSequence = reverseComplement(substr($sequence,$relStart,$length));
			    #$tk->start("getAdjTotalLibraryCounts") if ($testTime);
			    my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($prodList);
			    #$tk->end("getAdjTotalLibraryCounts") if ($testTime);
			    print HPL "$newId\t$side\t$newType\t$adjProdCount\t";
			    print HPL "$adjMaxProdCount\t$relStart\t$relStop\t$productStrand\t";
			    print HPL "$productSequence";
			    foreach my $bamElement (@{$bamList}) {
				my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
					printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
			    }
			    print HPL "\n";
			    foreach my $read (@{$prodList}) {
				my($relStart,$relStop,$offset,$gStart,
				   $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
				my $gStop = $gStart + ($relStop - $relStart);
				my $adjCount = $count * $adjSeqCount;
				print HDR "$newId\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
				print HDR "$offset\t$relStart\t$relStop\t$seq";
				foreach my $bamElement (@{$bamList}) {
				    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
				    printf(HDR "\t%.3f",$adjLibraryCounts->{$sample});
				}
				print HDR "\n";
			    }				
			}
			#############################
			#############################
			my($leftCenter, $rightCenter) = @{$center};
			print HL "$newId\t$chrom\t$start\t$stop\t$strand\t$leftCenter\t$rightCenter\t";
			printf(HL "%.3f\t%.3f\t%.3f\t",$totalSense,$totalAntisense,$mfe);
			print HL "$sequence\t$fold\t";
			printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\t",
			       $apd,$tpd,$urf,$ahc,$afh,$pbp,$sameShift,$bothShift,$splitProductAbundance,$overlapProductAbundance);
			printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t",
			       $aa,$ac,$ag,$at,$ca,$cc,$cg,$ct,$ga,$gc,$gg,$gt,$ta,$tc,$tg,$tt,$duplexEnergy,$gcContent);
			printf(HL "%.3f\t%.3f\t%.3f\t%d\t",
			       $foldDuplexCmp,$exactFoldDuplexCmp,$prodDupPBP,$mpDuplexLoopLength);
#			printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%d\t",
#			       $zScore, $foldDuplexCmp,$exactFoldDuplexCmp,$prodDupPBP,$mpDuplexLoopLength);
			printf(HL "%.3f\t%.3f",
			       $averageProductVariance,$weightedAverageProductVariance);
			printf(HL "\t%.3f\t%.3f",
			       $averageReadVariance,$weightedAverageReadVariance);
			printf(HL "\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
			       $maxProductCount,$mpDuplexCount,$duplexOverlap,$mpLoopDistance,$dupLoopDistance,$loopSize,$mpOverlapScore,$dupOverlapScore,$wghtMPOverlapIntensity,$wghtDupOverlapIntensity);
			printf(HL "\t%.3f",
			       $opa2);
			printf(HL "\t%.3f\t%.3f\t%.3f\t%.3f", $totalOverlapAmount,$averageOverlapAmount,$totalRelativeOverlapAmount,$averageRelativeOverlapAmount);
			printf(HL "\t%.3f", $maxOverlap);
			printf(HL "\t%.3f",  $innerLoopGapCount);
			my $totalSenseRPM = getHairpinRPM($totalSense,$librarySizes);
			my $totalAntisenseRPM = getHairpinRPM($totalAntisense,$librarySizes);
			printf(HL "\t%.3f\t%.3f", $totalSenseRPM,$totalAntisenseRPM);
			if ($maxProductScore) {
			    printf(HL "\t%.3f\n", $maxProductScore);
			} else {
			    print HL "\n";
			}
#			    } else {
#				open(FRC,">>$parameters->{filteredCandidatesFile}");
#				print FRC "$newId\t$location\t$REASON\n";
#				close(FRC);
#			    }
#			} else {
#			    open(FRC,">>$parameters->{filteredCandidatesFile}");
#			    print FRC "$newId\t$location\trejected: no good products.\n";
#			    close(FRC);
#			}
		    } 
		}
	    } else {
		$noReadsCount++;
		print "no reads found at $id\t$location\n" if ($testReadCounts);
	    }
	}
    }
    $mu->record("after finishing") if ($testMemory);
    #$tk->end("TotalTime") if ($testTime);
    #$tk->dumpTimes() if ($testTime);

    print "hairpins with reads = $readsCount\n" if ($testReadCounts);
    print "hairpins without reads = $noReadsCount\n" if ($testReadCounts);
    print "hairpin centers = $centersReadsCount\n" if ($testReadCounts);
    print "hairpins without centers (not including those with no reads) = $noCenterCount\n" if ($testReadCounts);
    print "hairpins with reads less than countMinLocus (not including those with no reads) = $readsLessThanCountMinLocus\n" if ($testReadCounts);

    close(RRF);
    close(HPL);
    close(HDR);
    close(HL);
#    close(FVF);
    $mu->dump() if ($testMemory);
}

sub oldProcessReadRegions {
    my($bamList,$readRegionsFile,$genomeDir,$chromLengths,$parameters) = @_;
    my $hairpinVectorFile = $parameters->{outputPrefix} . "_hairpinScores";
    my $testTime = 0;
    my $testMemory = 0;
    my $testReadCounts = 0;
    my $readsCount = 0; #number of mirBase hairpins with reads
    my $noReadsCount = 0;  #number of mirBase hairpins with no reads
    my $readsLessThanCountMinLocus = 0; #number of hairpins with fewer than count min locus reads
    my $centersReadsCount = 0; #number of centers on mirBase hairpins with reads
    my $noCenterCount = 0;  #number of hairpins without centers (hopefully this will stay 0)
    my $mu = Memory::Usage->new() if ($testMemory);
    my #$tk = TimeKeeper->new() if ($testTime);
    my @sampleList;
    my $minLength = $parameters->{lengthMin} or die "minLength: not loaded.\n";
    my $filteredCandidatesFile = $parameters->{filteredCandidatesFile};
    my $hairpinsFile = $parameters->{hairpinsFile};
    my $distinctReadsFile = $parameters->{distinctReadsFile};
    my $productFile = $parameters->{productFile};
    my $featureVectorFile = $parameters->{featureVectorFile};
    my $maxPredictedProductScores = (-e $parameters->{predProductClasses}) ? getMaxPredictedProductScores($parameters->{predProductClasses}) : 0;
    open(HL,">".$hairpinsFile) or die "failed to open $hairpinsFile for writing\n";
    open(HDR,">".$distinctReadsFile) or die "failed to open $distinctReadsFile for writing\n";
    open(HPL,">".$productFile) or die "failed to open $productFile for writing\n";
    print HL "#tag\tchrom\tstart\tstop\tstrand\tleftCenter\trightCenter\ttotalSense\ttotalAntisense\tmfe\tsequence\tfold\t";
    print HL "aapd\ttapd\turf\tahc\tmpfh\tpbp\tsameShifted\tbothShift\tsplitProductAbundance\toverlapProductAbundance\t";
    print HL "aa\tac\tag\tat\tca\tcc\tcg\tct\tga\tgc\tgg\tgt\tta\ttc\ttg\ttt\tduplexEnergy\tGCcontent\tzScore";
    if ($maxPredictedProductScores) {
	print HL "\tRFProductAvg\n";
    } else {
	print HL "\n";
    }
    print HDR "#tag\tstart\tstop\tstrand\ttotal reads\toffsetStart\trelative Start\trelative Stop\tsequence";
    print HPL "#tag\tside type\ttype\ttotal reads\ttotal most abundant\tstart\tstop\tstrand\tsequence";
    foreach my $bamListElement (@{$bamList}) {
	my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamListElement};
	print HDR "\t$sample";
        print HPL "\t$sample";
    }
    print HDR "\n";
    print HPL "\n";

    #$tk->start("TotalTime") if ($testTime);

    open(RRF, $readRegionsFile) or die "failed to open $readRegionsFile\n";
    my $prevChrom;
    my $genome;
    while (<RRF>) {
        chomp;
	unless(/\#/) {	  
	    my($chrom,$start,$stop,$id,$score,$strand) = split(/\t/,$_);
	    $start++; #converting to one based
	    #the chroms in RRF should be sorted at this point so the next line of code will only be executed once per chrom
	    my $maxProductScore = ($maxPredictedProductScores->{$id}) ? $maxPredictedProductScores->{$id} : 0;
	    if ($prevChrom) {
		$genome = loadGenome("$genomeDir/$chrom.fa") unless ($prevChrom eq $chrom);	
	    } else {
		$genome = loadGenome("$genomeDir/$chrom.fa");	
	    }
	    $prevChrom = $chrom;
	    my $location = extendSequenceToFillBuffer($chrom,$start,$stop,$strand,$chromLengths->{$chrom},$parameters);
	    ($chrom,$start,$stop,$strand) = parseLocation($location);
	    my $sequence = getSequence($location,$genome);
    	    #$tk->start("retrieveReadData") if ($testTime);
	    my($distinctReads,$uniqueReadCount,$totalDistinctReads,$adjustedReadCount) = retrieveReadData($location,$bamList,$parameters);
	    #$tk->end("retrieveReadData") if ($testTime);
	    if ($distinctReads->{$strand}) {  #neccessary if testing read reagions with known mirs and no real read exists in that read region
		if ($adjustedReadCount->{$strand} < $parameters->{countMinLocus}) {
		    $readsLessThanCountMinLocus++;
		    print "$id has fewer reads than count min locus\n" if ($testReadCounts);
		}
		$readsCount++;
		#$tk->start("buildProducts") if ($testTime);
		my($senseProducts) = buildProducts($location,$distinctReads->{$strand},$strand,$parameters);
		#$tk->end("buildProducts") if ($testTime);
		#$tk->start("getFold") if ($testTime);
		my($fold,$mfe,$maxProductCoords,$productDuplexCoords,
		   $seqDuplexCoords,$dupStructure,$duplexEnergy) = getFold($chrom,$sequence,$senseProducts,$id);
		my $foldDuplexComparison = compareFoldAndDuplexStructs($fold,$maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure);
		#$tk->end("getFold") if ($testTime);
		#$tk->start("getMergedHairpinCenters") if ($testTime);
		my @centers = getMergedHairpinCenters($fold,$parameters);
		#$tk->end("getMergedHairpinCenters") if ($testTime);
		$centersReadsCount += scalar(@centers);
		unless (scalar(@centers)) {
		    $noCenterCount++;
		    print "$id has no centers\n" if ($testReadCounts);
		    print "$fold\n" if ($testReadCounts);
		}
		my $COUNT= 0;
		foreach my $center (@centers) {
		    my $asciiVal = ord('a') + $COUNT;
		    my $label = chr($asciiVal);
		    my $newId = $id . $label;
		    $COUNT++;
		    #$tk->start("getBasePairs") if ($testTime);
		    my $basePairs = getBasePairs($center,$fold);
		    #$tk->end("getBasePairs") if ($testTime);
		    #$tk->start("getMaxHairpinLength") if ($testTime);
		    my $hairpinLength = getMaxHairpinLength($fold,$center);
		    #$tk->end("getMaxHairpinLength") if ($testTime);
		    #$tk->start("getFullHairpinWindow") if ($testTime);
		    my($hStart,$hStop) = getFullHairpinWindow($fold,$center);
		    #$tk->end("getFullHairpinWindow") if ($testTime);
		    if($hairpinLength >= $minLength) {
			#$tk->start("getProductInfo") if ($testTime);
			my $senseProducts=getProductInfo($center,$fold,$basePairs,$location,$senseProducts,$parameters);
			#$tk->end("getProductInfo") if ($testTime);
			#$tk->start("rebuildProducts") if ($testTime);
			$senseProducts=rebuildProducts($center,$fold,$basePairs,$senseProducts,$parameters);
			#$tk->end("rebuildProducts") if ($testTime);
			$senseProducts=addProductSize($senseProducts,$parameters);
			if(goodProducts($senseProducts,$parameters)) {
			    my $revStrand = revStrand($strand);
			    #$tk->start("buildProducts") if ($testTime);
			    my($antiSenseProducts) = buildProducts($location,$distinctReads->{$revStrand},$revStrand,$parameters);
			    #$tk->end("buildProducts") if ($testTime);
			    #$tk->start("getProductInfo") if ($testTime);
			    $antiSenseProducts = getProductInfo($center,$fold,$basePairs,$location,$antiSenseProducts,$parameters);
			    #$tk->end("getProductInfo") if ($testTime);
			    #$tk->start("rebuildProducts") if ($testTime);
			    $antiSenseProducts = rebuildProducts($center,$fold,$basePairs,$antiSenseProducts,$parameters);
			    #$tk->end("rebuildProducts") if ($testTime);
			    $antiSenseProducts=addProductSize($antiSenseProducts,$parameters);
			    #$tk->start("getAdjTotalProductReads") if ($testTime);
			    my $adjTotalProductReads = getAdjTotalProductReads($senseProducts);
			    #$tk->end("getAdjTotalProductReads") if ($testTime);
			    #$tk->start("getAdjTotalProductReads") if ($testTime);
			    my $adjTotalRevProductReads = getAdjTotalProductReads($antiSenseProducts);
			    #$tk->end("getAdjTotalProductReads") if ($testTime);
			    my($PASS,$REASON) = plausibleReads($adjTotalProductReads,$adjTotalRevProductReads,$senseProducts,$antiSenseProducts,
							       $parameters);
			    if($PASS) {
				#$tk->start("getReverseProductDisplacement") if ($testTime);
				my($tpd,$totalRP)=newGetReverseProductDisplacement($senseProducts,
										$antiSenseProducts,
										$adjTotalProductReads,
										$parameters,$newId);
				#$tk->end("getReverseProductDisplacement") if ($testTime);
				my $apd = $totalRP ? $tpd/$totalRP : 0.0;
				my $urf = $adjustedReadCount->{$strand} ? $uniqueReadCount->{$strand} / $adjustedReadCount->{$strand} : 0;
				#$tk->start("computeMaxProdHitCount") if ($testTime);
				my $ahc = computeMaxProdHitCount($senseProducts,$location,$distinctReads,$parameters);
				#$tk->end("computeMaxProdHitCount") if ($testTime);
				#$tk->start("computeMaxProdFivePrimeHet") if ($testTime);
				my $afh = computeMaxProdFivePrimeHet($senseProducts,$parameters);
				#$tk->end("computeMaxProdFivePrimeHet") if ($testTime);
				#$tk->start("computeProductBasePairing") if ($testTime);
				my $pbp = computeProductBasePairing($center,$senseProducts,$basePairs,$parameters);
				#$tk->end("computeProductBasePairing") if ($testTime);
				#$tk->start("computeMaxSameShift") if ($testTime);
				my $sameShift = computeMaxSameShift($location,$distinctReads->{$strand},$senseProducts,$parameters);
				#$tk->end("computeMaxSameShift") if ($testTime);
				#$tk->start("computeMaxBothShift") if ($testTime);
				my $bothShift = computeMaxBothShift($basePairs,$location,$distinctReads->{$strand},$senseProducts,$parameters,$newId);
				#$tk->end("computeMaxBothShift") if ($testTime);
				my($zScore,$consensusMFE) = getRNAzData($sequence);
				my $splitProductAbundance = computeSplitProductAbundance($senseProducts,$parameters);
				my $overlapProductAbundance = computeOverlapProductAbundance($senseProducts,$parameters);
				my($aa,$ac,$ag,$at,$ca,$cc,$cg,$ct,$ga,$gc,$gg,$gt,$ta,$tc,$tg,$tt) = getDinucleotideFrequencies($sequence);
				my $sequenceLength = length($sequence);
				my $gcContent = getGCcontent($sequence,$sequenceLength);
				# total is the read count for the whole hairpin
				my $totalSense = 0;
				foreach my $product (@{$senseProducts}) {
				    my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
				       $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
				    my $length = $relStop - $relStart + 1;
				    my $productSequence = substr($sequence,$relStart,$length);
				    $totalSense += $adjProdCount;
				    #$tk->start("getAdjTotalLibraryCounts") if ($testTime);
				    my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($prodList);
				    #$tk->end("getAdjTotalLibraryCounts") if ($testTime);
				    print HPL "$newId\t$side\t$newType\t$adjProdCount\t";
				    print HPL "$adjMaxProdCount\t$relStart\t$relStop\t$productStrand\t";
				    print HPL "$productSequence";
				    foreach my $bamElement (@{$bamList}) {
					my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
					printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
				    }
				    print HPL "\n";
				    foreach my $read (@{$prodList}) {
					my($relStart,$relStop,$offset,$gStart,
					   $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
					my $gStop = $gStart + ($relStop - $relStart);
					my $adjCount = $count * $adjSeqCount;
					print HDR "$newId\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
					print HDR "$offset\t$relStart\t$relStop\t$seq";
					foreach my $bamElement (@{$bamList}) {
					    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
					    printf(HDR "\t%.3f",$adjLibraryCounts->{$sample});
					}
					print HDR "\n";
				    }				
				}
				#############################
				# ANTISENSE PRODUCTS
				my $totalAntisense = 0;
				foreach my $product (@{$antiSenseProducts}) {
				    my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
				       $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
				    $newType = "as-".$newType;
				    $totalAntisense += $adjProdCount;
				    my $length = $relStop - $relStart + 1;
				    my $productSequence = reverseComplement(substr($sequence,$relStart,$length));
				    #$tk->start("getAdjTotalLibraryCounts") if ($testTime);
				    my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($prodList);
				    #$tk->end("getAdjTotalLibraryCounts") if ($testTime);
				    print HPL "$newId\t$side\t$newType\t$adjProdCount\t";
				    print HPL "$adjMaxProdCount\t$relStart\t$relStop\t$productStrand\t";
				    print HPL "$productSequence";
				    foreach my $bamElement (@{$bamList}) {
					my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
					printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
				    }
				    print HPL "\n";
				    foreach my $read (@{$prodList}) {
					my($relStart,$relStop,$offset,$gStart,
					   $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
					my $gStop = $gStart + ($relStop - $relStart);
					my $adjCount = $count * $adjSeqCount;
					print HDR "$newId\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
					print HDR "$offset\t$relStart\t$relStop\t$seq";
					foreach my $bamElement (@{$bamList}) {
					    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
					    printf(HDR "\t%.3f",$adjLibraryCounts->{$sample});
					}
					print HDR "\n";
				    }				
				}
				#############################
				#############################
				my($leftCenter, $rightCenter) = @{$center};
				print HL "$newId\t$chrom\t$start\t$stop\t$strand\t$leftCenter\t$rightCenter\t";
				printf(HL "%.3f\t%.3f\t%.3f\t",$totalSense,$totalAntisense,$mfe);
				print HL "$sequence\t$fold\t";
				printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\t",
				       $apd,$tpd,$urf,$ahc,$afh,$pbp,$sameShift,$bothShift,$splitProductAbundance,$overlapProductAbundance);
				printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t",
				       $aa,$ac,$ag,$at,$ca,$cc,$cg,$ct,$ga,$gc,$gg,$gt,$ta,$tc,$tg,$tt,$duplexEnergy,$gcContent);
				printf(HL "%.3f",
				       $zScore);
				if ($maxProductScore) {
				    printf(HL "\t%.3f\n", $maxProductScore);
				} else {
				    print HL "\n";
				}
			    } else {
				open(FRC,">>$parameters->{filteredCandidatesFile}");
				print FRC "$newId\t$location\t$REASON\n";
				close(FRC);
			    }
			} else {
			    open(FRC,">>$parameters->{filteredCandidatesFile}");
			    print FRC "$newId\t$location\trejected: no good products.\n";
			    close(FRC);
			}
		    } 
		}
	    } else {
		$noReadsCount++;
		print "no reads found at $id\t$location\n" if ($testReadCounts);
	    }
	}
    }
    $mu->record("after finishing") if ($testMemory);
    #$tk->end("TotalTime") if ($testTime);
    #$tk->dumpTimes() if ($testTime);

    print "hairpins with reads = $readsCount\n" if ($testReadCounts);
    print "hairpins without reads = $noReadsCount\n" if ($testReadCounts);
    print "hairpin centers = $centersReadsCount\n" if ($testReadCounts);
    print "hairpins without centers (not including those with no reads) = $noCenterCount\n" if ($testReadCounts);
    print "hairpins with reads less than countMinLocus (not including those with no reads) = $readsLessThanCountMinLocus\n" if ($testReadCounts);

    close(RRF);
    close(HPL);
    close(HDR);
    close(HL);
#    close(FVF);
    $mu->dump() if ($testMemory);
}



sub mirPreprocess {
    my($bamList,$hairpins,$genomeDir,$chromLengths,$parameters) = @_;
    my @sampleList;
    my $minLength = $parameters->{lengthMin} or die "minLength: not loaded.\n";
    my $filteredCandidatesFile = $parameters->{filteredCandidatesFile};
    my $hairpinsFile = $parameters->{hairpinsFile};
    my $distinctReadsFile = $parameters->{distinctReadsFile};
    my $productFile = $parameters->{productFile};
    my $readsCount = 0; #number of mirBase hairpins with reads
    my $centersReadsCount = 0; #number of centers on mirBase hairpins with reads
    my $noReadsCount = 0;  #number of mirBase hairpins with no reads
    open(FRC,">".$filteredCandidatesFile) or die "failed to open $filteredCandidatesFile for writing\n";
    close(FRC);
    open(HL,">".$hairpinsFile) or die "failed to open $hairpinsFile for writing\n";
    open(HDR,">".$distinctReadsFile) or die "failed to open $distinctReadsFile for writing\n";
    open(HPL,">".$productFile) or die "failed to open $productFile for writing";
    print HL "#tag\tchrom\tstart\tstop\tstrand\tleftCenter\trightCenter\ttotalSense\ttotalAntisense\tmfe\tsequence\tfold\t";
    print HL "aapd\ttapd\turf\tahc\tafh\tvalid\tsameShifted\tbothShift\tsplitProductAbundance\toverlapProductAbundance\n";
    print HDR "#tag\tstart\tstop\tstrand\ttotal reads\toffsetStart\trelative Start\trelative Stop\tsequence";
    print HPL "#tag\tside type\ttype\ttotal reads\ttotal most abundant\tstart\tstop\tstrand\tsequence";
    foreach my $bamListElement (@{$bamList}) {
	my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamListElement};
	print HDR "\t$sample";
        print HPL "\t$sample";
    }
    print HDR "\n";
    print HPL "\n";

    my $prevChrom;
    my $genome;
    foreach my $chrom (keys %{$hairpins}) {
	if ($prevChrom) {
	    $genome = loadGenome("$genomeDir/$chrom.fa") unless ($prevChrom eq $chrom);	
	} else {
	    $genome = loadGenome("$genomeDir/$chrom.fa");	
	}
	$prevChrom = $chrom;
	foreach my $hairpin (@{$hairpins->{$chrom}}) {
	    my ($start, $stop, $strand, $id, $mirName) = @$hairpin;
	    my $location = "$chrom:$start..$stop:$strand";
	    my $sequence = getSequence($location,$genome);
	    my($distinctReads,$uniqueReadCount,$totalDistinctReads,$adjustedReadCount) = retrieveReadData($location,$bamList,$parameters);
	    if ($distinctReads->{$strand}) {  #neccessary if testing read reagions with known mirs and no real read exists in that read region
		$readsCount++;
		my($senseProducts) = buildProducts($location,$distinctReads->{$strand},$strand,$parameters);
		my($fold,$mfe,$maxProductCoords,$productDuplexCoords,
		   $seqDuplexCoords,$dupStructure,$duplexEnergy) = getFold($chrom,$sequence,$senseProducts,$mirName);
		my @centers = getMergedHairpinCenters($fold,$parameters);
		$centersReadsCount += scalar(@centers);
		my $COUNT= 0;
		foreach my $center (@centers) {
		    my $asciiVal = ord('a') + $COUNT;
		    my $label = chr($asciiVal);
		    my $newId = $mirName . $label;
		    print ("$mirName has more than one label\n") if ($label eq 'b');
		    $COUNT++;
		    my $basePairs = getBasePairs($center,$fold);
		    my $hairpinLength = getMaxHairpinLength($fold,$center);
		    my($hStart,$hStop) = getFullHairpinWindow($fold,$center);
		    my $senseProducts=getProductInfo($center,$fold,$basePairs,$location,$senseProducts,$parameters);
		    $senseProducts=rebuildProducts($center,$fold,$basePairs,$senseProducts,$parameters);
		    my $revStrand = revStrand($strand);
		    my($antiSenseProducts) = buildProducts($location,$distinctReads->{$revStrand},$revStrand,$parameters);
		    $antiSenseProducts = getProductInfo($center,$fold,$basePairs,$location,$antiSenseProducts,$parameters);
		    $antiSenseProducts = rebuildProducts($center,$fold,$basePairs,$antiSenseProducts,$parameters);
		    my $adjTotalProductReads = getAdjTotalProductReads($senseProducts);
		    my $adjTotalRevProductReads = getAdjTotalProductReads($antiSenseProducts);
		    my($tpd,$totalRP)=getReverseProductDisplacement($senseProducts,
								    $antiSenseProducts,
								    $adjTotalProductReads,
								    $parameters);
		    my $apd = $totalRP ? $tpd/$totalRP : 0.0;
		    my $urf = $adjustedReadCount->{$strand} ? $uniqueReadCount->{$strand} / $adjustedReadCount->{$strand} : 0;
		    my $ahc = computeMaxProdHitCount($senseProducts,$location,$distinctReads,$parameters);
		    my $afh = computeMaxProdFivePrimeHet($senseProducts,$parameters);
		    my $pbp = computeProductBasePairing($center,$senseProducts,$basePairs,$parameters);
		    my $sameShift = computeMaxSameShift($location,$distinctReads->{$strand},$senseProducts,$parameters);
		    my $bothShift = computeMaxBothShift($basePairs,$location,$distinctReads->{$strand},$senseProducts,$parameters);
		    my $splitProductAbundance = computeSplitProductAbundance($senseProducts,$parameters);
		    my $overlapProductAbundance = computeOverlapProductAbundance($senseProducts,$parameters);
		    # total is the read count for the whole hairpin
		    my $totalSense = 0;
		    foreach my $product (@{$senseProducts}) {
			my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
			   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
			    my $length = $relStop - $relStart + 1;
			my $productSequence = substr($sequence,$relStart,$length);
			$totalSense += $adjProdCount;
			my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($prodList);
			print HPL "$newId\t$side\t$newType\t$adjProdCount\t";
			print HPL "$adjMaxProdCount\t$relStart\t$relStop\t$productStrand\t";
			print HPL "$productSequence";
			foreach my $bamElement (@{$bamList}) {
			    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
			    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
			}
			print HPL "\n";
			    foreach my $read (@{$prodList}) {
				my($relStart,$relStop,$offset,$gStart,
				   $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
				my $gStop = $gStart + ($relStop - $relStart);
				my $adjCount = $count * $adjSeqCount;
				print HDR "$newId\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
				print HDR "$offset\t$relStart\t$relStop\t$seq";
				foreach my $bamElement (@{$bamList}) {
				    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
				    printf(HDR "\t%.3f",$adjLibraryCounts->{$sample});
					}
				print HDR "\n";
			    }				
		    }
		    #############################
		    # ANTISENSE PRODUCTS
		    my $totalAntisense = 0;
		    foreach my $product (@{$antiSenseProducts}) {
			my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
			   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
			$newType = "as-".$newType;
			$totalAntisense += $adjProdCount;
			my $length = $relStop - $relStart + 1;
			my $productSequence = reverseComplement(substr($sequence,$relStart,$length));
			my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($prodList);
			print HPL "$newId\t$side\t$newType\t$adjProdCount\t";
			print HPL "$adjMaxProdCount\t$relStart\t$relStop\t$productStrand\t";
			print HPL "$productSequence";
			foreach my $bamElement (@{$bamList}) {
			    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
			    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
			}
			print HPL "\n";
			foreach my $read (@{$prodList}) {
			    my($relStart,$relStop,$offset,$gStart,
			       $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
			    my $gStop = $gStart + ($relStop - $relStart);
			    my $adjCount = $count * $adjSeqCount;
			    print HDR "$newId\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
			    print HDR "$offset\t$relStart\t$relStop\t$seq";
			    foreach my $bamElement (@{$bamList}) {
				my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
				printf(HDR "\t%.3f",$adjLibraryCounts->{$sample});
			    }
			    print HDR "\n";
			}				
		    }
		    #############################
		    #############################
		    my($leftCenter, $rightCenter) = @{$center};
		    print HL "$newId\t$chrom\t$start\t$stop\t$strand\t$leftCenter\t$rightCenter\t";
		    printf(HL "%.3f\t%.3f\t%.3f\t",$totalSense,$totalAntisense,$mfe);
		    print HL "$sequence\t$fold\t";
		    printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\n",
			   $apd,$tpd,$urf,$ahc,$afh,$pbp,$sameShift,$bothShift,$splitProductAbundance,$overlapProductAbundance);
		}
	    } else {
		$noReadsCount++;
		print "no reads found at $mirName\t$location\n";
	    }
	}
    }
    print "hairpins with reads = $readsCount\n";
    print "hairpin centers with reads = $centersReadsCount\n";
    print "hairpins without reads = $noReadsCount\n";

    close(RRF);
    close(HPL);
    close(HDR);
}

#pairMap subroutine from Michelle Wiley
#(altered to create map for 0-based coords)
sub pairMap {  
    my ($structure) = @_;
    my @map;
    my $len = length($structure);
    my @left;
    for (my $i=0; $i<$len; $i++) {
	if (substr($structure,$i,1) eq "("){
	    push(@left,$i);
	} elsif (substr($structure,$i,1) eq ")"){
	    $map[$left[scalar(@left)-1]] = $i; #or $map[$left[$#left]]
	    $map[$i] = pop(@left);
	} else {
	    $map[$i] = ".";
	}
    }
    return(\@map);
}

sub duplexMap {
    my($maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure) = @_;
    my %duplexMap;
    my($prodDupStruct,$seqDupStruct) = $dupStructure =~ /^(.*)&(.*)/;
    my $len = length($prodDupStruct);
    #print "$prodDupStruct & $seqDupStruct    len=$len\n";
    #adding the first coordinant returned by RNAduplex for the product to the products coordinants
    #to get a starting location relative to the sequence.
    my $productStart = $maxProductCoords->[0] + $productDuplexCoords->[0];
    my $j = length($seqDupStruct) - 1;
    #print "productStart = $productStart\t\$j = $j\n";
    for (my $i=0; $i < length($prodDupStruct); $i++) {
	#print "\$i = $i\n";
	if (substr($prodDupStruct,$i,1) eq "(") {
	    #print "prodDupStruct equals \(\n";
	    while ($j >= 0 && substr($seqDupStruct,$j,1) ne ")") {
		$duplexMap{$j + $seqDuplexCoords->[0]} = ".";
		#print "\$j = $j\t" . ($j + $seqDuplexCoords->[0]) . " = .\n";
		$j--;
	    }
	    $duplexMap{$productStart + $i} = $j + $seqDuplexCoords->[0];
	    $duplexMap{$j + $seqDuplexCoords->[0]} = $productStart + $i;
	    $j--;
	} else {
	    $duplexMap{$i} = ".";
	}
    }
    while ($j >= 0) {
	$duplexMap{$j + $seqDuplexCoords->[0]} = ".";
	$j--;
    }
    return \%duplexMap;
}

sub compareFoldAndDuplexStructs {
    my($fold,$maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure,$readId) = @_;
    my($foldDuplexCmp,$exactFoldDuplexCmp) = (0,0);
    #print "$readId\n$fold\n";
    #print "length = " . length($fold) . "\n";
    #print "maxProductCoords = \(". $maxProductCoords->[0] . "," . $maxProductCoords->[1] . ")\n";
    #print "prodDuplexCoords = \(". $productDuplexCoords->[0] . "," . $productDuplexCoords->[1] . ")\n";
    #print "seqDuplexCoords = \(". $seqDuplexCoords->[0] . "," . $seqDuplexCoords->[1] . ")\n";
    my $map = pairMap($fold);
    my $duplexMap = duplexMap($maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure);
    my($prodDupStruct,$seqDupStruct) = $dupStructure =~ /^(.*)&(.*)/;  
    my $prodStart = $maxProductCoords->[0] + $productDuplexCoords->[0];
    my $total = length($prodDupStruct) + length($seqDupStruct);
    my $totalSame = 0;
    my $totalSimilar = 0;
    foreach my $i (sort {$a <=> $b} keys %{$duplexMap}) {
	#print "\$i = $i\n";
	#print $duplexMap->{$i} . "\t";
	#print $map->[$i] . "\n";
	if ($duplexMap->{$i} eq $map->[$i]) {
	    $totalSame++;
	}
	if (($duplexMap->{$i} ne "." && $map->[$i] ne ".") ||
	    ($duplexMap->{$i} eq "." && $map->[$i] eq ".")) {
	    $totalSimilar++;
	}
    }
    $foldDuplexCmp = $totalSimilar / $total;
    $exactFoldDuplexCmp = $totalSame / $total;
    return($foldDuplexCmp,$exactFoldDuplexCmp);
}

sub getProductDuplexPBP {
    my($fold,$maxProductCoords) = @_;
    #note: this is for the duplex created by RNAfold rather than RNAduplex
    my($prodDupPBP,$mpDuplexLoopLength) = (0,0);
    my $map = pairMap($fold);
    #print "\n$fold\n";
    #print "maxProductCoords: " . $maxProductCoords->[0] . "-" . $maxProductCoords->[1] . "\n";
    #finding duplex from RNA structure
    my $i = $maxProductCoords->[0];
    while ($map->[$i] eq "." || $map->[$i] eq "+") {
	#print "$i ---> " . $map->[$i];
	if ($i > $maxProductCoords->[1]) {
	    #print "no mapping found returning 0 and $i\n";
	    return(0,$i-1);
	}
	$i++;
    }
    my @prodDuplex = ($map->[$i],$map->[$i]);
    while ($i <= $maxProductCoords->[1]) {
	last if ($i >= length($fold));  #some of the max products extend over the hairpin
	#print $i . "\n";
	unless ($map->[$i] eq "." || $map->[$i] eq "+") {
	    if ($map->[$i] < $prodDuplex[0]) {
		$prodDuplex[0] = $map->[$i];
	    } elsif ($map->[$i] > $prodDuplex[1]) {
		$prodDuplex[1] = $map->[$i];
	    }
	}
	$i++;
	#print "prodDuplex: " . $prodDuplex[0] . "\t" . $prodDuplex[1] . "\n";
    }
    #getting the longest loop count and pbp of the duplex
    my($newLoopCount,$bindCount,$total) = (0,0,0);
    for (my $j = $prodDuplex[0]; $j <= $prodDuplex[1]; $j++) {
	my $mapped = ($map->[$j] eq ".") ? 0 : 1;
	if ($mapped) {
	    $bindCount++;
	    if ($newLoopCount > $mpDuplexLoopLength) {
		$mpDuplexLoopLength = $newLoopCount;
	    }
	    $newLoopCount = 0;
	} else {
	    $newLoopCount++;
	}
    }
    $prodDupPBP = $bindCount / ($prodDuplex[1] - $prodDuplex[0] + 1);
    #print "prodDupPBP = $prodDupPBP\n";
    return($prodDupPBP,$mpDuplexLoopLength);
}

sub initializeRNAFold {
    $RNA::noLonelyPairs = 1;
};

sub getFold {
    my($chrom,$sequence,$senseProducts,$hairpinId) = @_;
    my $test = 0;
    my $majorProduct = $senseProducts->[0];  #products must be sorted by abundance at this point
    my($id,$productReads,$adjMaxProductCount,$adjustedProductCount,$maxRelStart,$maxRelStop,$maxGStart,$strand) = @{$majorProduct};
    $maxRelStart = 0 if ($maxRelStart < 0);  #sometimes the major product begins before the sequence if the original read region did not cover it before extending the sequence
    print "$hairpinId\n" if ($test);
    my $productSequence = substr($sequence,$maxRelStart,$maxRelStop-$maxRelStart+1);
    my @maxProductCoords = ($maxRelStart,$maxRelStop);
    print "Most Abundant Product Sequence:$productSequence\n" if ($test);
    my $productSeqLength = length($productSequence);
    my $maskedSequence = maskSequence($sequence,$maxRelStart,$maxRelStop);
    print "sequence:\t\t$sequence\n" if ($test);
    print "masked sequence:\t$maskedSequence\n" if ($test);
    #$productSequence = reverse($productSequence);  #testing out reverseing the product before getting the duplex
    my ($productDuplexCoords,$seqDuplexCoords,$dupStructure,$duplexEnergy) = getDuplex($productSequence,$maskedSequence);
    print "Most Abundant Product Coords: $maxRelStart $maxRelStop\n" if ($test);
    print "productDuplexCoords: $productDuplexCoords->[0] $productDuplexCoords->[1]\n" if ($test);
    print "sequenceDuplexCoords: $seqDuplexCoords->[0] $seqDuplexCoords->[1]\n" if ($test);
    print "duplexSequence: " . substr($maskedSequence,$seqDuplexCoords->[0],$seqDuplexCoords->[1]-$seqDuplexCoords->[0]+1) . "\n" if ($test);

    my $duplexSide = ($maxRelStart > $seqDuplexCoords->[0]) ? '5p' : '3p'; #in reality this is only true for real mirs but it works for our purpose here
    print "duplex side = $duplexSide\n" if ($test);

    
    #mirStarShift is the length of the mir* not covered by the duplex
    #in rare cases the duplex region may overlap a mor but at the moment there is no way to check if a product is a mor until the hairpin is folded
    my $mirStarShift = getMirStarShift($seqDuplexCoords,$duplexSide,$senseProducts);
    print "mirStarShift = $mirStarShift\n" if ($test);
    #mirShift is the buffer to add due to the duplexFold not finding the duplex energy for the full product sequence
    my $mirShift = ($duplexSide eq '5p') ? $productSeqLength - $productDuplexCoords->[1] - 1 : $productDuplexCoords->[0];
    print "mirShift = $mirShift\n" if ($test);
    my $maxMirShift = ($mirStarShift > $mirShift) ? $mirStarShift: $mirShift;
    my($foldStart,$foldStop) = extendHairpinToCoverMirs($majorProduct,$seqDuplexCoords,$maxMirShift,$duplexSide);
    print "adding $maxMirShift to both sides and 2 to the $duplexSide side\nfoldStart = $foldStart foldStop = $foldStop\n" if ($test);
    $foldStart = 0 if ($foldStart < 0);
    $foldStop = length($sequence) - 1 if ($foldStop > (length($sequence)-1));

    my $foldSequence = substr($sequence,$foldStart,$foldStop-$foldStart+1);
    print "foldSequence=\t$foldSequence\n" if ($test);
    my($fold,$mfe) = RNA::fold($foldSequence);
    my $leftOutRegion = '+' x $foldStart;
    my $rightOutRegion = '+' x (length($sequence)-$foldStop-1);
    my $newFold = $leftOutRegion . $fold . $rightOutRegion;
    print "sequence=\t$sequence\n" if ($test);
    print "fold=    \t$newFold\n" if ($test);
    print "\n\n" if ($test);
    return $newFold,$mfe,\@maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure,$duplexEnergy;
}

sub getDuplex {
    my($seq1,$seq2) = @_;
    my $dup = RNA::duplexfold($seq1,$seq2);
    my $duplexEnergy = $dup->{energy};
    my $dupStructure = $dup->{structure};
    my($seq1DupStruct,$seq2DupStruct) = $dupStructure =~ /^(.*)&(.*)/;  
    my $seq1DupLength = length($seq1DupStruct);             
    my $seq2DupLength = length($seq2DupStruct); 
    my @seq1DupCoords = ($dup->{i}-$seq1DupLength, $dup->{i} - 1);   
    my @seq2DupCoords = ($dup->{j} - 1,$dup->{j}+$seq2DupLength-2); 
    return (\@seq1DupCoords,\@seq2DupCoords,$dupStructure,$duplexEnergy);
}

sub oldGetDuplex {
    my($seq1,$seq2) = @_;
    my $dup = RNA::duplexfold($seq1,$seq2);
    my $duplexEnergy = $dup->{energy};
    my $dupStructure = $dup->{structure};
    my($seq1DupStruct,$seq2DupStruct) = $dupStructure =~ /^(.*)&(.*)/;  #needed to get the lengths of the duplex coords.
    my $seq1Length = length($seq1);
    my $seq1DupStructLength = length($seq1DupStruct);             
    my $seq2DupStructLength = length($seq2DupStruct);                          
    my @seq1DuplexCoords = ($dup->{i}-$seq1DupStructLength, $dup->{i} - 1);   
    my @seq2DuplexCoords = ($dup->{j} - 1,$dup->{j}+$seq2DupStructLength-2); 
    return (\@seq1DuplexCoords,\@seq2DuplexCoords,$duplexEnergy);
}

sub getMirStarShift {
    my($duplexCoords,$duplexSide,$products) = @_;
    my($leftDuplexCoord,$rightDuplexCoord) = @{$duplexCoords};
    if ($duplexSide eq '5p') {
	foreach my $product (@{$products}) {
	    my($id,$productReads,$adjMaxProductCount,$adjustedProductCount,$maxRelStart,$maxRelStop,$maxGStart,$strand) = @{$product};
	    if (getOverlap($maxRelStart,$maxRelStop,$leftDuplexCoord,$rightDuplexCoord)) {
		if ($maxRelStart < $leftDuplexCoord) {
		    return $leftDuplexCoord - $maxRelStart;
		}
	    }
	}
    } elsif ($duplexSide eq '3p') {
	foreach my $product (@{$products}) {
	    my($id,$productReads,$adjMaxProductCount,$adjustedProductCount,$maxRelStart,$maxRelStop,$maxGStart,$strand) = @{$product};
	    if (getOverlap($maxRelStart,$maxRelStop,$leftDuplexCoord,$rightDuplexCoord)) {
		if ($maxRelStop > $rightDuplexCoord) {
		    return $maxRelStop - $rightDuplexCoord;
		}
	    }
	}
    } else {
	die "no duplexSide given\n";
    }
    return 0;
}

sub extendHairpinToCoverMirs {
    my($majorProduct,$seqDuplexCoords,$mirBuffer,$mirStarSide) = @_;
    my($id,$productReads,$adjMaxProductCount,$adjustedProductCount,$maxRelStart,$maxRelStop,$maxGStart,$strand) = @{$majorProduct};
    my($foldStart,$foldStop);
    if ($mirStarSide eq '5p') {
	$foldStart = $seqDuplexCoords->[0] - $mirBuffer;
	$foldStop = $maxRelStop + $mirBuffer;
	$foldStart = $foldStart - 2;  #since there are two unpaired nucleotides on the hairpin
    } else {
	$foldStart = $maxRelStart - $mirBuffer;
	$foldStop = $seqDuplexCoords->[1] + $mirBuffer;
	$foldStop = $foldStop + 2;   #since there are two unpaired nucleotides on the hairpin
    }

    return ($foldStart,$foldStop);
}

sub oldGetFold {
    my($chrom,$sequence,$senseProducts,$hairpinId) = @_;
    my $test = 1;
    my $majorProduct = $senseProducts->[0];  #products must be sorted by this point
    my($id,$productReads,$adjMaxProductCount,$adjustedProductCount,$maxRelStart,$maxRelStop,$maxGStart,$strand) = @{$majorProduct};
    print "$hairpinId\n" if ($test);
    print "maxRelStart = $maxRelStart\tmaxRelStop = $maxRelStop\n" if ($test);
    my $productSequence = substr($sequence,$maxRelStart,$maxRelStop-$maxRelStart+1); #converted to zero based to get substring
    print "product Sequence:$productSequence\n" if ($test);
    my $maskedSequence = maskSequence($sequence,$maxRelStart,$maxRelStop);
    print "sequence:\t\t$sequence\n" if ($test);
    print "masked sequence:\t$maskedSequence\n" if ($test);
    my $dup = RNA::duplexfold($productSequence,$maskedSequence);
    my $dupStructure = $dup->{structure};
    my($productSeqDupStruct,$seqDupStruct) = $dupStructure =~ /^(.*)&(.*)/;  #needed to get the lengths.
    my $productSeqLength = $maxRelStop - $maxRelStart + 1;
    my $productSeqDupLength = length($productSeqDupStruct);             
    my $seqDupLength = length($seqDupStruct);                           #a better way is to use the lengths from the structure to find coords
    my @productDuplexCoords = ($dup->{i}-$productSeqDupLength, $dup->{i} - 1);   #duplex coords of main product in 0 based
    print "ProductCoords: $maxRelStart $maxRelStop\n" if ($test);
    print "productDuplexCoords: $productDuplexCoords[0] $productDuplexCoords[1]\n" if ($test);
    my @seqDuplexCoords = ($dup->{j} - 1,$dup->{j}+$seqDupLength-2); #duplex coords of sequence in 0 based
    print "sequenceDuplexCoords: $seqDuplexCoords[0] $seqDuplexCoords[1]\n" if ($test);
    print "duplexSequence: " . substr($maskedSequence,$seqDuplexCoords[0],$seqDuplexCoords[1]-$seqDuplexCoords[0]+1) . "\n" if ($test);
    my $pairedDuplexSide = ($maxRelStart > $seqDuplexCoords[0]) ? 'left' : 'right';
    print "duplex side = $pairedDuplexSide\n" if ($test);
    #overlappingProductBuffer is the buffer to add to the region that is due to overlapping products
    my $overlappingProductBuffer = getOverlappingProductBuffer(\@seqDuplexCoords,$pairedDuplexSide,$senseProducts);
    print "overlappingProductBuffer = $overlappingProductBuffer\n" if ($test);
    #maxProductBuffer is the buffer to add due to the duplexFold not finding the duplex energy for the full product sequence
    my $maxProductBuffer = ($pairedDuplexSide eq 'left') ? $productSeqLength - $productDuplexCoords[1] - 1 : $productDuplexCoords[0];
    print "maxProductBuffer = $maxProductBuffer\n" if ($test);

    my($foldStart,$foldStop,$foldSequence);    
    if ($pairedDuplexSide eq 'left') {
	if ($overlappingProductBuffer) {
	    $foldStart = $seqDuplexCoords[0] - $overlappingProductBuffer - 2;
	    print "subtracting overlapping product buffer: FoldStart = $seqDuplexCoords[0] - $overlappingProductBuffer - 2 = $foldStart\n" if ($test);
	} else {
	    $foldStart = $seqDuplexCoords[0] - $maxProductBuffer - 2;
	    print "no Overlapping product: subtracting maxProductBuffer: FoldStart = $seqDuplexCoords[0] - $maxProductBuffer - 2 = $foldStart\n" if ($test);
	}
	$foldStop = $maxRelStop + $maxProductBuffer;
	print "adding maxProductBuffer: FoldStop = $maxRelStop + $maxProductBuffer = $foldStop\n" if ($test);
    } elsif ($pairedDuplexSide eq 'right') {
	if ($overlappingProductBuffer) {
	    $foldStop = $seqDuplexCoords[1] + $overlappingProductBuffer + 2;
	    print "adding overlapping product buffer: FoldStop = $seqDuplexCoords[1] + $overlappingProductBuffer + 2 = $foldStop\n" if ($test);
	} else {
	    $foldStop = $seqDuplexCoords[1] + $maxProductBuffer + 2;
	    print "no Overlapping product: adding maxProductBuffer: FoldStop = $seqDuplexCoords[1] + $maxProductBuffer + 2 = $foldStop\n" if ($test);
	}
	$foldStart = $maxRelStart - $maxProductBuffer;
	print "subtracting maxProductBuffer: FoldStart = $maxRelStart - $maxProductBuffer = $foldStart\n" if ($test);
    }
    $foldStart = 0 if ($foldStart < 0);
    $foldStop = length($sequence) - 1 if ($foldStop > (length($sequence)-1));
    $foldSequence = substr($sequence,$foldStart,$foldStop-$foldStart+1);
    print "foldSequence=\t$foldSequence\n" if ($test);
    my($fold,$mfe) = RNA::fold($foldSequence);
    my $leftOutRegion = '+' x $foldStart;
    my $rightOutRegion = '+' x (length($sequence)-$foldStop-1);
    my $newFold = $leftOutRegion . $fold . $rightOutRegion;
    print "sequence=\t$sequence\n" if ($test);
    print "fold=    \t$newFold\n" if ($test);
    print "\n\n" if ($test);
    return $newFold,$mfe,\@seqDuplexCoords;
}

sub getRNAzData {
    my($sequence) = @_;
    my $output = `echo -n "CLUSTAL W(1.81) multiple sequence alignment\n\n\nseq1 $sequence\nseq2 $sequence" | RNAz`;
    my @lines = split("\n",$output);
    my ($zScore,$consensusMFE);
    foreach my $line (@lines) {
	if ($line =~ /Mean\sz\-score\:/) {
	    ($zScore) = $line =~ /Mean\sz\-score\:\s+(\-?\d+\.\d+)$/;
	}
	if ($line =~ /Consensus\sMFE:/) {
	    ($consensusMFE) = $line =~ /Consensus\sMFE:\s+(\-?\d+\.\d+)$/;
	}
    }
    return($zScore,$consensusMFE);
}

sub computeOverlapRatio {
    my($productAbundance1,$productAbundance2) = @_;
    my $overlapRatio;
    if ($productAbundance1 > $productAbundance2) {
	$overlapRatio = $productAbundance2 / $productAbundance1;
    } else {
	$overlapRatio = $productAbundance1 / $productAbundance2;
    }
    return $overlapRatio;
}

sub getInnerLoopGapCount {
    my($center,$fold) = @_;
    my $minInnerLoopGapSize = 3;
    my $inLoopCount = 0;
    my $totalGaps = 0;
    my($leftCenter,$rightCenter) = @{$center};
    my $foldSubStr = substr($fold,$leftCenter,$rightCenter - $leftCenter + 1);
    my @loopFold = split('',$foldSubStr);
    for (my $i = 0; $i < @loopFold; $i++) {
	if ($loopFold[$i] eq ".") {
	    $inLoopCount = $inLoopCount + 1;
	} else {
	    if ($inLoopCount >= $minInnerLoopGapSize) {
		$totalGaps = $totalGaps + 1;
	    }
	    $inLoopCount = 0;
	}
    }
    if ($inLoopCount >= $minInnerLoopGapSize) {
	$totalGaps = $totalGaps + 1;
    }
#    if ($totalGaps > 1) {
#	print "$fold\nleftCenter, $rightCenter\n";
#	print "$foldSubStr\n";
#    }
    #print "total gaps = $totalGaps\n";
    return $totalGaps;
}

sub getMaximumOverlap {
    my($productInfo,$parameters) = @_;
    my $maxOverlap = 0;
    foreach my $product (@{$productInfo}) {
	foreach my $product2(@{$productInfo}) {
	    my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop) = @{$product};
	    my($side2,$newType2,$prodList2,$adjProdCount2,$adjMaxProdCount2,$relStart2,$relStop2) = @{$product2};
	    my $overlap = getOverlap($relStart,$relStop,$relStart2,$relStop2);
	    unless ($product == $product2) {
		if ($overlap > $maxOverlap) {
		    $maxOverlap = $overlap;
		}
	    }
	}
    } 
    return $maxOverlap;
}

sub computeOverlapAmounts {
    my($productInfo,$parameters) = @_;
    my($totalOverlapAmount,$averageOverlapAmount,$totalRelativeOverlapAmount,$averageRelativeOverlapAmount,$productTotal) = (0,0,0,0,0);
    foreach my $product(@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop) = @{$product};
	$productTotal += $adjProdCount;
    }
    my @sortedProducts = sort {$a->[3] <=> $b->[3]} @{$productInfo};  #sorted from least to most abundant
    for (my $i = 0; $i < @sortedProducts; $i++) {
	my $product1 = $sortedProducts[$i];
	my($side1,$newType1,$prodList1,$adjProdCount1,$adjMaxProdCount1,$relStart1,$relStop1) = @{$product1};
	unless(($newType1 =~ /long/) || ($newType1 =~ /short/)) {
	    for (my $j = $i + 1; $j < @sortedProducts; $j++) {
		my $product2 = $sortedProducts[$j];
		my($side2,$newType2,$prodList2,$adjProdCount2,$adjMaxProdCount2,$relStart2,$relStop2) = @{$product2};
		unless (($newType2 =~ /long/) || ($newType2 =~ /short/)) {
		    if (my $overlap = getOverlap($relStart1,$relStop1,$relStart2,$relStop2)) {
			my $smallProductCount = ($adjProdCount1 < $adjProdCount2) ? $adjProdCount1 : $adjProdCount2;
			my $bigProductCount = ($adjProdCount1 < $adjProdCount2) ? $adjProdCount2 : $adjProdCount1;
			$totalOverlapAmount = $totalOverlapAmount + $overlap * $smallProductCount;
			$averageOverlapAmount = $averageOverlapAmount + $overlap * $smallProductCount / $productTotal;
			$totalRelativeOverlapAmount = $totalRelativeOverlapAmount + $overlap * $smallProductCount / $bigProductCount;
			$averageRelativeOverlapAmount = $averageRelativeOverlapAmount + ($overlap * $smallProductCount / $bigProductCount) / $productTotal; 
		    }
		}
	    }
	}	
    }
    return ($totalOverlapAmount,$averageOverlapAmount,$totalRelativeOverlapAmount,$averageRelativeOverlapAmount);
}


sub computeOverlapProductAbundance {
    my($productInfo,$parameters) = @_;
    my $minOverlap = $parameters->{OverlapMin};
    my $productTotal = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop) = @{$product};
	$productTotal += $adjProdCount;
    }
    my @sortedProducts = sort {$a->[3] <=> $b->[3]} @{$productInfo};  #sorted from least to most abundant
    my $maxOverlapProductAbundance = 0;
    for (my $i = 0; $i < @sortedProducts; $i++) {
	my $product1 = $sortedProducts[$i];
	my($side1,$newType1,$prodList1,$adjProdCount1,$adjMaxProdCount1,$relStart1,$relStop1) = @{$product1};
	unless(($newType1 =~ /long/) || ($newType1 =~ /short/)) {
	    for (my $j = $i + 1; $j < @sortedProducts; $j++) {
		my $product2 = $sortedProducts[$j];
		my($side2,$newType2,$prodList2,$adjProdCount2,$adjMaxProdCount2,$relStart2,$relStop2) = @{$product2};
		unless (($newType2 =~ /long/) || ($newType2 =~ /short/)) {
		    if (my $overlap = getOverlap($relStart1,$relStop1,$relStart2,$relStop2)) {
			if ($overlap > $minOverlap) {
			    my $overlapProductAbundance = $adjProdCount1 / $productTotal;
			    $maxOverlapProductAbundance = $overlapProductAbundance if ($overlapProductAbundance > $maxOverlapProductAbundance);
			}
		    }
		}
	    }
	}	
    }
    return $maxOverlapProductAbundance;
}

sub computeOverlapProductAbundance2 {
    my($productInfo,$parameters) = @_;
    my $minOverlap = $parameters->{OverlapMin};
    my $productTotal = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop) = @{$product};
	$productTotal += $adjProdCount;
    }
    my @sortedProducts = sort {$a->[3] <=> $b->[3]} @{$productInfo};  #sorted from least to most abundant
    my $overlapProductAbundance = 0;
    for (my $i = 0; $i < @sortedProducts; $i++) {
	my $product1 = $sortedProducts[$i];
	my($side1,$newType1,$prodList1,$adjProdCount1,$adjMaxProdCount1,$relStart1,$relStop1) = @{$product1};
	unless(($newType1 =~ /long/) || ($newType1 =~ /short/)) {
	    for (my $j = $i + 1; $j < @sortedProducts; $j++) {
		my $product2 = $sortedProducts[$j];
		my($side2,$newType2,$prodList2,$adjProdCount2,$adjMaxProdCount2,$relStart2,$relStop2) = @{$product2};
		unless (($newType2 =~ /long/) || ($newType2 =~ /short/)) {
		    if (my $overlap = getOverlap($relStart1,$relStop1,$relStart2,$relStop2)) {
			if ($overlap > $minOverlap) {
			    $overlapProductAbundance = $overlapProductAbundance + $overlap * $adjProdCount1 * $adjProdCount2;
			}
		    }
		}
	    }
	}	
    }
    return $overlapProductAbundance;
}

#this version of computeOverlapProductAbundance does not take into account product size
sub oldComputeOverlapProductAbundance {
    my($productInfo,$parameters) = @_;
    my $productTotal = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop) = @{$product};
	$productTotal += $adjProdCount;
    }
    my @sortedProducts = sort {$a->[3] <=> $b->[3]} @{$productInfo};  #sorted from least to most abundant
    my $maxOverlapProductAbundance = 0;
    for (my $i = 0; $i < @sortedProducts; $i++) {
	my $product1 = $sortedProducts[$i];
	my($side1,$newType1,$prodList1,$adjProdCount1,$adjMaxProdCount1,$relStart1,$relStop1) = @{$product1};
	for (my $j = $i + 1; $j < @sortedProducts; $j++) {
	    my $product2 = $sortedProducts[$j];
	    my($side2,$newType2,$prodList2,$adjProdCount2,$adjMaxProdCount2,$relStart2,$relStop2) = @{$product2};
	    if (getOverlap($relStart1,$relStop1,$relStart2,$relStop2)) {
		my $overlapProductAbundance = $adjProdCount1 / $productTotal;
		$maxOverlapProductAbundance = $overlapProductAbundance if ($overlapProductAbundance > $maxOverlapProductAbundance);
	    }
	}	
    }
    return $maxOverlapProductAbundance;
}

sub computeSplitProductAbundance {
    my($productInfo,$parameters) = @_;
    my $productTotal = 0;
    my $splitProductTotal = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop) = @{$product};
	unless(($newType =~ /long/) || ($newType =~ /short/)) {
	    $splitProductTotal += $adjProdCount if ($newType eq 'split');
	    $productTotal += $adjProdCount;
	}
    }
    return ($productTotal == 0) ? 0 : $splitProductTotal / $productTotal;
}

#this version of computeSplitProductAbundance does not take into account product size
sub oldComputeSplitProductAbundance {
    my($productInfo,$parameters) = @_;
    my $productTotal = 0;
    my $splitProductTotal = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop) = @{$product};
	$splitProductTotal += $adjProdCount if ($newType eq 'split');
	$productTotal += $adjProdCount;
    }
    return ($productTotal == 0) ? 0 : $splitProductTotal / $productTotal;
}

sub buildProducts {
    my($location,$distinctReads,$productStrand,$parameters) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $PRODUCTCOUNT = 1;
    my %productHash;
    my %relStartHash;
    foreach my $dRead (@{$distinctReads}) {
	# dStart is genomic position.
	my($dStart,$dStop,$total,$hitCount,$libraryCounts,$seq,$adjustedSeqCount,$adjustedLibraryCounts) = @{$dRead};
	my $relStart = ($strand eq "+") ? $dStart - $start : $stop - $dStop;
	my $relStop = ($strand eq "+") ? $dStop - $start : $stop - $dStart;
	#print "absolute: $dStart..$dStop:$productStrand\n";
	#print "relative: $relStart..$relStop\n";
	#offsets will be determined later after the fold
	my $parsedRead = [$relStart,$relStop,$dStart,$total,$hitCount,$libraryCounts,$seq,$adjustedSeqCount,$adjustedLibraryCounts,$productStrand];
	if (my $id = overlapCurrent(\%productHash,$parsedRead,$parameters)) {
	    push(@{$productHash{$id}},$parsedRead);
	} else {
	    push(@{$productHash{$PRODUCTCOUNT}},$parsedRead);
	    $PRODUCTCOUNT++;
	}
    }
    my @productInfo;
    foreach my $id (keys %productHash) {
	my $adjTotal=0;
	my($maxRelStart,$maxRelStop,$maxGStart);
	my $FIRST=1;
	my $adjProductCount = 0;
	my $adjMaxProductCount = 0;
	my %dReadCounts;
	foreach my $storedRead (@{$productHash{$id}}) {
	    my($relStart,$relStop,$dStart,$readCount,$hitCount,$libraryCounts,$seq,$adjustedSeqCount) = @{$storedRead};
	    $dReadCounts{$relStart."x".$relStop} += $readCount * $adjustedSeqCount;
	    $relStartHash{$relStart."x".$relStop} = $relStart;
	}
	my @keys = sort {$dReadCounts{$b} <=> $dReadCounts{$a} || $relStartHash{$a} <=> $relStartHash{$b}} keys %dReadCounts;
#	my @keys = sort {$dReadCounts{$b} <=> $dReadCounts{$a}} keys %dReadCounts;
	my $maxKey = shift(@keys);
	($maxRelStart,$maxRelStop) = split(/x/,$maxKey); 
	foreach my $storedRead (@{$productHash{$id}}) {
	    my($relStart,$relStop,$dStart,$readCount,$hitCount,$libraryCounts,$seq,$adjustedSeqCount) = @{$storedRead};
	    $adjProductCount += $readCount * $adjustedSeqCount;
	    if(($relStart == $maxRelStart)&&($relStop == $maxRelStop)) {
		$adjMaxProductCount += $readCount * $adjustedSeqCount;
		$maxGStart = $dStart;
	    }
	}
	push(@productInfo,[$id,$productHash{$id},$adjMaxProductCount,$adjProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand]);
    }
    my @sortedProducts = sort {$b->[3] <=> $a->[3]} @productInfo;
    return \@sortedProducts;
}

sub maskSequence {
    my($sequence,$start,$stop) = @_;
    my $seqLength = length($sequence);
    my $maskLength = $stop-$start+1;
    my $maskedSequence;
    #note the product sequence is in one base but the sequence is in 0 base
    if ($stop == $seqLength) {
	return substr($sequence,0,$start) . 'N' x ($maskLength - 1);
    }
    return substr($sequence,0,$start) . 'N' x ($maskLength) . substr($sequence, $stop + 1, $seqLength - $stop - 1);
}

sub extractProducts {
    my($center,$fold,$basePairs,$location,$distinctReads,$drStrand,$parameters) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $productInfo = getProductInfo($center,$fold,$basePairs,
				     $start,$stop,$strand,
				     $distinctReads,$drStrand,$parameters);
    $productInfo = rebuildProducts($center,$fold,$basePairs,
				    $productInfo,$parameters);    
    return $productInfo;
}

sub computeMaxSameShift {
    my($location,$strandDistinctReads,$productInfo, $parameters) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $newParameters = {};
    foreach my $key (keys %{$parameters}) {
	$newParameters->{$key} = $parameters->{$key};
    }
    $newParameters->{distanceMin} = $parameters->{shiftMin};
    my $adjTotalReads = getAdjTotalReads($strandDistinctReads,$parameters);
    return getSameArmShift($productInfo,$adjTotalReads,$newParameters);
}

sub computeMaxBothShift {
    my($basePairs,$location,$strandDistinctReads,$productInfo,$parameters,$mirId) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $newParameters = {};
    foreach my $key (keys %{$parameters}) {
	$newParameters->{$key} = $parameters->{$key};
    }
    $newParameters->{distanceMin} = $parameters->{shiftMin};
    my $adjTotalReads = getAdjTotalReads($strandDistinctReads,$parameters);
    if(bothArmProducts($basePairs,$productInfo,$newParameters)) {
	#print "possible both arm shift\n";
	return newGetBothArmShift($basePairs,$productInfo,$adjTotalReads,$newParameters,$location,$mirId);
    }
    return 0;
}

sub bothArmProducts {
    my($basePairs,$productInfo,$parameters) = @_;
    my %sideHash;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop,$offset,$gStart) = @{$product};
	unless(($newType =~ /long/) || ($newType =~ /short/)) {
	     if($side eq "3p") {
		 if(overlapsBasePairs($basePairs,$relStart,$relStop)) {
		     $sideHash{$side}++;
		 }
	     } else {
		 $sideHash{$side}++;
	     }
	}
    }
    if(($sideHash{"5p"})&&($sideHash{"3p"})) {
	return 1;
    }
    return 0;
}

sub oldBothArmProducts {
    my($basePairs,$productInfo,$parameters) = @_;
    my %sideHash;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$product};	
	if($side eq "3p") {
	    if(overlapsBasePairs($basePairs,$relStart,$relStop)) {
		$sideHash{$side}++;
	    }
	} else {
	    $sideHash{$side}++;
	}
    }
    if(($sideHash{"5p"})&&($sideHash{"3p"})) {
	return 1;
    }
    return 0;
}


sub newGetBothArmShift {
    my($basePairs,$productInfo,$adjTotalReads,$parameters,$location,$newId) = @_;
    my $test = 0;
    my $minReadFraction = 0.01;
    my $countThreshold = $adjTotalReads*$minReadFraction;
    my $minOverlap = $parameters->{OverlapMin} or die "FAIL: no minOverlap loaded.\n";
    my %FIVE_PRIME_USED;
    my %THREE_PRIME_USED;
    my @productOverlaps;
    my $maxShift = 0;
    print "\n\n$newId\n" if ($test);
    print "minReadFraction = $minReadFraction\tcountThreshold = $countThreshold\tminOverlap = $minOverlap\n" if ($test);
     for(my $i=0;$i<@{$productInfo};$i++) {
	 my($side1,$newType1,$prodList1,$adjProdCount1,$adjMaxProdCount1,$relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	 unless(($newType1 =~ /long/) || ($newType1 =~ /short/)) {
	     if($adjProdCount1 >= $countThreshold) {	
		 if($side1 eq "5p") {
		     print "adjProdCount1 = $adjProdCount1 >= $countThreshold = countThreshold for $side1-$newType1\n" if ($test); 
		     for(my $j=0;$j<@{$productInfo};$j++) {
			 if($i != $j) {
			     my($side2,$newType2,$prodList2,$adjProdCount2,$adjMaxProdCount2,$relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
			     if (computeOverlapRatio($adjProdCount1,$adjProdCount2) > $parameters->{maxOverlapProductAbundance}) {
				 unless(($newType2 =~ /long/) || ($newType2 =~ /short/)) {
				     if($adjProdCount2 >= $countThreshold) {
					 if($side2 eq "3p") {				
					     print "\tadjProdCount2 = $adjProdCount2 >= $countThreshold = countThreshold for $side2-$newType2\n" if ($test);
					     if(overlapsBasePairs($basePairs,$relStart2,$relStop2)) {
						 my($relStart2bp,$relStop2bp) = getOuterBasePairsFivePrime($basePairs,
													   $relStart2,
													   $relStop2,
													   $location);
						 my $overlap = getOverlap($relStart1,$relStop1,
									  $relStart2bp,$relStop2bp);
						 my $shift = getShift($relStart1,$relStop1,
								      $relStart2bp,$relStop2bp);
						 push(@productOverlaps, [$overlap,$shift,$i,$j]);
						 print "\tBase pairs for $side1-$newType1 and $side2-$newType2 overlap by $overlap" if ($test);
						 print "(shift = $shift)\n" if ($test);
					     }	
					 }
				     }
				 }
			     }
			 }
		     }
		 }
	     }
	 }
     }
    #finding the maximum shift for products that are sorted from most abundant to least abundant overlap
    @productOverlaps = sort {$b->[0] <=> $a->[0]} @productOverlaps;
    foreach my $item (@productOverlaps) {
	my($overlap,$shift,$i,$j) = @{$item};
	unless(($FIVE_PRIME_USED{$i})||($THREE_PRIME_USED{$j})) {
	    $maxShift = $shift if ($shift > $maxShift);
	    $FIVE_PRIME_USED{$i}++;
	     $THREE_PRIME_USED{$j}++;  
	}
    }
    print "maximum shift =  $maxShift\n" if ($test);
    
    return $maxShift;
}



sub getBothArmShift {
    my($basePairs,$productInfo,$adjTotalReads,$parameters,$location,$newId) = @_;
    my $test = 0;
    my $minReadFraction = 0.01;
    my $countThreshold = $adjTotalReads*$minReadFraction;
    my $minOverlap = $parameters->{OverlapMin} or die "FAIL: no minOverlap loaded.\n";
    my $maxShift = 0;
    print "\n\n$newId\n" if ($test);
    print "minReadFraction = $minReadFraction\tcountThreshold = $countThreshold\tminOverlap = $minOverlap\tmaxShfit = $maxShift\n" if ($test);
    for(my $i=0;$i<@{$productInfo};$i++) {
	my($side1,$newType1,$prodList1,$adjProdCount1,$adjMaxProdCount1,$relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	unless(($newType1 =~ /long/) || ($newType1 =~ /short/)) {
	    if($adjProdCount1 >= $countThreshold) {	
		if($side1 eq "5p") {
		    print "adjProdCount1 = $adjProdCount1 >= $countThreshold = countThreshold for $side1-$newType1\n" if ($test); 
		    for(my $j=0;$j<@{$productInfo};$j++) {
			if($i != $j) {
			    my($side2,$newType2,$prodList2,$adjProdCount2,$adjMaxProdCount2,$relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
			    unless(($newType2 =~ /long/) || ($newType2 =~ /short/)) {
				if($adjProdCount2 >= $countThreshold) {
				    if($side2 eq "3p") {				
					print "\tadjProdCount2 = $adjProdCount2 >= $countThreshold = countThreshold for $side2-$newType2\n" if ($test);
					if(overlapsBasePairs($basePairs,$relStart2,$relStop2)) {
					    my($relStart2bp,$relStop2bp) = getOuterBasePairsFivePrime($basePairs,
												      $relStart2,
												      $relStop2,
												      $location);
					    my $overlap = getOverlap($relStart1,$relStop1,
								     $relStart2bp,$relStop2bp);
					     print "\tBase pairs for $side1-$newType1 and $side2-$newType2 overlap by $overlap\n" if ($test);
					     my $shift = getShift($relStart1,$relStop1,
								  $relStart2bp,$relStop2bp);
					    print "\tshift = $shift\n" if ($test);
					    if($overlap > $minOverlap) {
						if(abs($shift) > abs($maxShift)) {				
						    $maxShift = $shift;
						    print "\tmaxShift = $maxShift\n" if ($test);
						}
					    }
					}	
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }    
    return $maxShift;
}

#this version of getBothArmShift does not take into account product size
sub oldGetBothArmShift {
    my($basePairs,$productInfo,$adjTotalReads,$parameters,$location) = @_;
    my $minReadFraction = 0.01;
    my $countThreshold = $adjTotalReads*$minReadFraction;
    my $minOverlap = $parameters->{OverlapMin} or die "FAIL: no minOverlap loaded.\n";
    my $maxShift = 0;
    for(my $i=0;$i<@{$productInfo};$i++) {
	 my($side1,$newType1,$prodList1,$adjProdCount1,$adjMaxProdCount1,
	    $relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	 if($adjProdCount1 >= $countThreshold) {	
	     if($side1 eq "5p") {
		 for(my $j=0;$j<@{$productInfo};$j++) {
		     if($i != $j) {
			 my($side2,$newType2,$prodList2,$adjProdCount2,$adjMaxProdCount2,
			    $relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
			 if($adjProdCount2 >= $countThreshold) {
			     if($side2 eq "3p") {				
				 if(overlapsBasePairs($basePairs,$relStart2,$relStop2)) {
				     my($relStart2bp,$relStop2bp) = getOuterBasePairsFivePrime($basePairs,
											       $relStart2,
											       $relStop2,
											       $location);
				     my $overlap = getOverlap($relStart1,$relStop1,
							      $relStart2bp,$relStop2bp);
				     my $shift = getShift($relStart1,$relStop1,
							  $relStart2bp,$relStop2bp);
				     if($overlap > $minOverlap) {
					 if(abs($shift) > abs($maxShift)) {				
					     $maxShift = $shift;
					 }
				     }
				 }	
			     }
			 }
		     }
		 }
	     }
	 }
    }
    return $maxShift;
}

sub overlapsBasePairs {
    my($basePairs,$relStart,$relStop) = @_;
    my @pairs;
    foreach my $pair (@{$basePairs}) {
	my($bp1,$bp2) = @{$pair};
	if(($relStart <= $bp2)&&($bp2 <= $relStop)) {
	     push(@pairs,[$bp1,$bp2]);
	}
    }
    if(@pairs > 1) {
	 return 1;
    } else {
	return 0;
    }
}

sub computeBasePairDensity {
    my($basePairs,$relStart,$relStop) = @_;
    my $basePairCount = 0;
     foreach my $pair (@{$basePairs}) {
	 my($bp1,$bp2) = @{$pair};
	 if(($relStart <= $bp2)&&($bp2 <= $relStop)||
	    ($relStart <= $bp1)&&($bp1 <= $relStop)) {	    
	     $basePairCount++;
	 }
     }
    return sprintf("%.2f",$basePairCount / ($relStop - $relStart + 1));
}

sub getOuterBasePairs {
    # NOTE: This assumes that relStart and relStop are on the 3p arm!
     my($basePairs,$relStart,$relStop,$location) = @_;
     my @pairs;
     foreach my $pair (@{$basePairs}) {
	 my($bp1,$bp2) = @{$pair};
	 # here assumes bp2 is basepair on 3p arm
	 if(($relStart <= $bp2)&&($bp2 <= $relStop)) {
	     push(@pairs,[$bp1,$bp2]);
	 }
     }
     if(@pairs > 1) {
	 # sort by bp position on 5' arm
	 @pairs = sort {$a->[0] <=> $b->[0]} @pairs;
	 my @runList;
	 my $runStart = 0;
	 my $runStop = 0;
	 my $runStartBp = 0;
	 my $runStopBp = 0;
	 my $INRUN = 0;
	 for(my $i=0;$i<@pairs-1;$i++) {
	     my($bp1,$bp2) = @{$pairs[$i]};
	     my($nbp1,$nbp2) = @{$pairs[$i+1]};	    
	     # if base pairs on both arms are adjacent nucleotides.
	     if(($bp1+1==$nbp1)&&($bp2-1==$nbp2)) {
		 # current base pairs are neighbors.
		 if($INRUN) {
		     # in a run of neighboring base pairs, extend the run
		     $runStop = $nbp1;
		     $runStopBp = $nbp2;
		 } else {
		     # not in a run, start a new run
		     $runStart = $bp1;
		     $runStop = $nbp1;
		     $runStartBp = $bp2;
		     $runStopBp = $nbp2;
		     $INRUN = 1;
		 }
	     } else {
		 if($INRUN) {
		     push(@runList,[$runStart,$runStop,$runStartBp,$runStopBp]);
		 }
		 $INRUN = 0;
	     }
	 }
	 if($INRUN) {
	     push(@runList,[$runStart,$runStop,$runStartBp,$runStopBp]);
	 }
	 my $maxLength = 0;
	 my $leftPair;
	 my $rightPair;
	 # here runStart is the most 5' basepair that pairs with a base that overlaps
	 # the relStart and relStop.
	 foreach my $run (@runList) {
	     my($runStart,$runStop,$runStartBp,$runStopBp) = @{$run};
	     if($runStop - $runStart > $maxLength) {
		 $maxLength = $runStop - $runStart;
		 $leftPair = [$runStart,$runStartBp];
		 $rightPair = [$runStop,$runStopBp];
	     }
	 }
	 if($maxLength) {
	     unless($leftPair) {
		 my $L = "L" x ($relStart+1);
		 my $R = "R" x ($relStop+1);
		 die "could not get leftPair:\n$L\n$R\nat $location\n";
	     }
	     # lbp1 is the most 5' base overlapping with relStart,relStop
	     my($lbp1,$lbp2) = @{$leftPair};
	     my $deltaL = $relStop - $lbp2;
	     unless($rightPair) {
		 my $L = "L" x ($relStart+1);
		 my $R = "R" x ($relStop+1);
		 die "could not get rightPair:\n$L\n$R\nat $location\n";
	     }
	     # lbp1 is the most 3' base overlapping with relStart,relStop
	     my($rbp1,$rbp2) = @{$rightPair};
	     my $deltaR = $rbp2 - $relStart;
	     if(($deltaL < 0)||($deltaR < 0)) {
		 die "left: $lbp1 $lbp2\nright:$rbp1 $rbp2\nwindow: $relStart $relStop\n";
	     }
	     #print "longest run: $lbp1,$lbp2 .. $rbp1,$rbp2\n";
	     return ($lbp1-$deltaL,$rbp1+$deltaR);
	 }
     }
     return (0,0);
}

sub getOuterBasePairsFivePrime {
    # NOTE: This assumes that relStart and relStop are on the 3p arm!
    my($basePairs,$relStart,$relStop,$location) = @_;
    my @pairs;
    foreach my $pair (@{$basePairs}) {
	my($bp1,$bp2) = @{$pair};
	# here assumes bp2 is basepair on 3p arm
	if(($relStart <= $bp2)&&($bp2 <= $relStop)) {
	    push(@pairs,[$bp1,$bp2]);
	}
    }
    # sort by bp position on 5' arm
    if(@pairs > 1) {
	@pairs = sort {$a->[0] <=> $b->[0]} @pairs;
	# lbp1 is the most 5' base overlapping with relStart,relStop
	my $leftPair = shift(@pairs);
	my $rightPair = pop(@pairs);
	my($lbp1,$lbp2) = @{$leftPair};
	my $deltaL = $relStop - $lbp2;
	unless($rightPair) {
	    my $L = "L" x ($relStart+1);
	    my $R = "R" x ($relStop+1);
	    die "could not get rightPair:\n$L\n$R\nat $location\n";
	}
	# lbp1 is the most 3' base overlapping with relStart,relStop
	my($rbp1,$rbp2) = @{$rightPair};
	my $deltaR = $rbp2 - $relStart;
	if(($deltaL < 0)||($deltaR < 0)) {
	    die "left: $lbp1 $lbp2\nright:$rbp1 $rbp2\nwindow: $relStart $relStop\n";
	}
	#print "longest run: $lbp1,$lbp2 .. $rbp1,$rbp2\n";
	return ($lbp1-$deltaL,$rbp1+$deltaR);
    }
    return (0,0);
}

sub getSameArmShift {
    my($productInfo,$adjTotalReads,$parameters) = @_;
    my $test = 0;
     my $minOverlap = $parameters->{OverlapMin} or die "FAIL: no minOverlap loaded.\n";
    my $minReadFraction = 0.01;
    my $countThreshold = $adjTotalReads*$minReadFraction;
    my $maxShift = 0;
    print "minOverlap = $minOverlap\tminReadFraction = $minReadFraction\ncountThreshold = $countThreshold\tmaxShift = $maxShift\n" if ($test);
    for(my $i=0;$i<@{$productInfo};$i++) {
	my($side1,$newType1,$prodList1,$adjProdCount1,$adjMaxProdCount1,$relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	 unless (($newType1 =~ /long/) || ($newType1 =~ /short/)) {
	     #print "$i: $side1 $newType1 $relStart1 $relStop1 $prodCount1 vs $countThreshold\n";
	     if($adjProdCount1 >= $countThreshold) {
		 print "relStart1 = $relStart1\trelStop1 = $relStop1\n" if ($test);
		 print "adjProdCount1 = $adjProdCount1 >= $countThreshold = countthreshold\n" if ($test);
		 for(my $j=0;$j<@{$productInfo};$j++) {
		     if($i != $j) {
			 my($side2,$newType2,$prodList2,$adjProdCount2,$adjMaxProdCount2,$relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
			 if (computeOverlapRatio($adjProdCount1,$adjProdCount2) > $parameters->{maxOverlapProductAbundance}) {
			     unless (($newType2 =~ /long/) || ($newType2 =~ /short/)) {
				 if($adjProdCount2 >= $countThreshold) {
				     print "relStart2 = $relStart2\trelStop2 = $relStop2\n" if ($test);
				     print "adjProdCount2 = $adjProdCount2 >= $countThreshold = countThreshold\n" if ($test);
				     my $overlap = getOverlap($relStart1,$relStop1,$relStart2,$relStop2);
				     my $shift = getShift($relStart1,$relStop1,$relStart2,$relStop2);
				     print "overlap = $overlap\tshift = $shift\n" if ($test);
				     if($overlap > $minOverlap) {
					 print "overlap = $overlap > $minOverlap = minOverlap\n" if ($test);
					 if($shift > $maxShift) {
					     print "shift = $shift > $maxShift = maxShift\n" if ($test);
					     $maxShift = $shift;
					     print "now maxShift = $shift\n" if ($test);
					 }
				     }
				 }
			     }
			 }
		     }
		 }
	     }
	 }
    }
    print "maxShift = $maxShift\n" if ($test);
    return $maxShift;
}

#this version of getSameArmShift does not take into account product size
sub oldGetSameArmShift {
    my($productInfo,$adjTotalReads,$parameters) = @_;
    my $minOverlap = $parameters->{OverlapMin} or die "FAIL: no minOverlap loaded.\n";
    my $minReadFraction = 0.01;
    my $countThreshold = $adjTotalReads*$minReadFraction;
    my $maxShift = 0;
    for(my $i=0;$i<@{$productInfo};$i++) {
	my($side1,$newType1,$prodList1,$adjProdCount1,$adjMaxProdCount1,$relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	#print "$i: $side1 $newType1 $relStart1 $relStop1 $prodCount1 vs $countThreshold\n";
	if($adjProdCount1 >= $countThreshold) {	    
	    for(my $j=0;$j<@{$productInfo};$j++) {
		if($i != $j) {
		    my($side2,$newType2,$prodList2,$adjProdCount2,$adjMaxProdCount2,
		       $relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
		     if($adjProdCount2 >= $countThreshold) {
			 my $overlap = getOverlap($relStart1,$relStop1,$relStart2,$relStop2);
			 my $shift = getShift($relStart1,$relStop1,$relStart2,$relStop2);
			 if($overlap > $minOverlap) {
			     if($shift > $maxShift) {
				 $maxShift = $shift;				
			     }
			 }
		     }
		}
	    }
	}
    }
    return $maxShift;
}

sub lociContained {
    my($relStart1,$relStop1,$relStart2,$relStop2) = @_;    
    if(($relStart1 <= $relStart2)&&($relStop2 <= $relStop1)) {
	return 1;
    } elsif(($relStart2 <= $relStart1)&&($relStop1 <= $relStop2)) {
	return 1;
    }
}

sub getOverlap {
    my($relStart1,$relStop1,$relStart2,$relStop2) = @_;    
    if(($relStart1 <= $relStart2)&&($relStart2 <= $relStop1)) {
	# relStart2 within first window
	return $relStop1 - $relStart2 + 1;
    } elsif(($relStart1 <= $relStop2)&&($relStop2 <= $relStop1)) {
	# if relStop2 within first window
	return $relStop2 - $relStart1 + 1;	
    } elsif(($relStart2 <= $relStart1)&&($relStart1 <= $relStop2)) {
	return $relStop2 - $relStart1 + 1;
    } elsif(($relStart2 <= $relStop1)&&($relStop1<=$relStop2)) {
	return $relStop2-$relStop1 + 1;
    }
    return 0;
}

sub getShift {
    my($relStart1,$relStop1,$relStart2,$relStop2) = @_;    
    return $relStart2-$relStart1;
}

sub doubleStrandedProducts {
    my($readCount,$revReadCount,$parameters) = @_;
    my $maxReverse = $parameters->{reverseMax};
    my $totalReads = $readCount + $revReadCount;
    my $fraction = $totalReads ? $revReadCount/$totalReads : 0;
    #print "doubleStrandedProducts $totalReads $readCount $revReadCount $fraction\n";
    if($fraction <= $maxReverse) {
	return 0;
    } else {
	return 1;
    }
}

sub overlappingProducts {
    my($productInfo,$revProductInfo,$parameters) = @_;
    for(my $i1=0;$i1<@{$productInfo};$i1++) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$productInfo->[$i1]};
	my $start = $relStart;
	my $stop = $relStop;
	for(my $i2=0;$i2<@{$revProductInfo};$i2++) {
	    my($aType,$aNewType,$aProdList,$aAdjProdCount,$aAdjMaxProdCount,
	       $aRelStart,$aRelStop,$aOffset,$aGStart) = @{$revProductInfo->[$i2]};
	    my $aStart = $aRelStart;
	    my $aStop = $aRelStop;
	    # if the products overlap..
	    if((($start <= $aStart)&&($aStart <= $stop))||
	       (($start <= $aStop)&&($aStop <= $stop))||
	       (($aStart <= $start)&&($start <= $aStop))||
	       (($aStart <= $stop)&&($stop <= $aStop))) {
		return 1;
	    }
	}
    }
    # if made it here, then no overlaps.
    return 0;
}


sub newGetReverseProductDisplacement {
    my($productInfo,$revProductInfo,$adjTotalReads,$parameters,$hairpinId) = @_;
    my $test = 0;
    print "\n\n$hairpinId\n" if ($test);
    my $minDispInit = 100;
    # make global variable VVV
    my $minFrac = 0.0005;
    my $sum = 0;
    my $total = 0;
    my %USED;
    my %AS_USED;
    my @distances;
    for(my $i1=0;$i1<@{$productInfo};$i1++) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop,$offset,$gStart) = @{$productInfo->[$i1]};
	unless (($newType =~ /long/) || ($newType =~ /short/)) {
	    my $start = $relStart;
	    my $stop = $relStop;
	    for(my $i2=0;$i2<@{$revProductInfo};$i2++) {
		my($aType,$aNewType,$aProdList,$aAdjProdCount,$aAdjMaxProdCount,
		   $aRelStart,$aRelStop,$aOffset,$aGStart) = @{$revProductInfo->[$i2]};
		unless (($aNewType =~ /long/) || ($aNewType =~ /short/)) {
		    if($aAdjProdCount >= $minFrac*$adjTotalReads) {
			#print "aAdjProdCount = $aAdjProdCount >= $minFrac * $adjTotalReads = minFrac x adjTotalReads\n" if ($test);
			my $aStart = $aRelStart;
			my $aStop = $aRelStop;
			# if the products overlap..
			if((($start <= $aStart)&&($aStart <= $stop))||
			   (($start <= $aStop)&&($aStop <= $stop))||
			   (($aStart <= $start)&&($start <= $aStop))||
			   (($aStart <= $stop)&&($stop <= $aStop))) {
			    if ((computeOverlapRatio($adjProdCount,$aAdjProdCount) > $parameters->{maxOverlapProductAbundance}) &&
				($aAdjProdCount > $parameters->{minAntisenseProduct})) {
				my $disp = abs($aStart - $start);
				push(@distances,[$disp,$i1,$i2]);
				print "\t$side-$newType overlaps as-$aNewType. Displacement = $disp\n" if ($test);
			    }
			}		    
		    }
		}
	    }
	}
    }
    print "\t" if ($test);
    @distances = sort {$a->[0] <=> $b->[0]} @distances;
    foreach my $item (@distances) {
	my($disp,$i1,$i2) = @{$item};
	unless(($USED{$i1})||($AS_USED{$i2})) {
	    $sum += $disp;
	    print "$disp + " if ($test);
	    $total++;
	    $USED{$i1}++;
	    $AS_USED{$i2}++;  
	}
    }
    print " = $sum\n" if ($test);
    if ($total) {
	print "\taapd = $sum / $total = " . $sum / $total . "\n" if ($test);
    }
    return($sum,$total);
}



sub getReverseProductDisplacement {
    my($productInfo,$revProductInfo,$adjTotalReads) = @_;
    my $test = 0;
    my $minDispInit = 100;
    my $minFrac = 0.0005;
    my $sum = 0;
    my $total = 0;
    my %USED;
    for(my $i1=0;$i1<@{$productInfo};$i1++) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop,$offset,$gStart) = @{$productInfo->[$i1]};
	unless (($newType =~ /long/) || ($newType =~ /short/)) {
	    print "side=$side\ttype=$newType\tadjProdCount=$adjProdCount\tadjMaxProdCount=$adjMaxProdCount\n" if ($test);
	    print "relStart=$relStart\trelStop=$relStop\toffset=$offset\tgStart=$gStart\n" if ($test);
	    my $start = $relStart;
	    my $stop = $relStop;
	    my $minDisp = $minDispInit;
	    print "minDisp=$minDisp\n" if ($test);
	    my $minI2;
	    for(my $i2=0;$i2<@{$revProductInfo};$i2++) {
		my($aType,$aNewType,$aProdList,$aAdjProdCount,$aAdjMaxProdCount,
		   $aRelStart,$aRelStop,$aOffset,$aGStart) = @{$revProductInfo->[$i2]};
		unless (($aNewType =~ /long/) || ($aNewType =~ /short/)) {
		    print "aSide=$aType\taType=$aNewType\taAdjProdCount=$aAdjProdCount\taAdjMaxProdCount=$aAdjMaxProdCount\n" if ($test);
		    print "aRelStart=$aRelStart\taRelStop=$aRelStop\taOffset=$aOffset\taGStart=$aGStart\n" if ($test);
		    if($aAdjProdCount >= $minFrac*$adjTotalReads) {
			print "aAdjProdCount = $aAdjProdCount >= $minFrac * $adjTotalReads = minFrac x adjTotalReads\n" if ($test);
			my $aStart = $aRelStart;
			my $aStop = $aRelStop;
			# if the products overlap..
			if((($start <= $aStart)&&($aStart <= $stop))||
			   (($start <= $aStop)&&($aStop <= $stop))||
			   (($aStart <= $start)&&($start <= $aStop))||
			   (($aStart <= $stop)&&($stop <= $aStop))) {
			    my $disp = abs($aStart - $start);
			    print "products overlap: disp = $disp\n" if ($test);
			    unless($USED{$i2}) {
				if($disp < $minDisp) {
				    $minDisp = $disp;
				    $minI2 = $i2;
				    print "disp = $disp < $minDisp = minDisp: now minDisp = $disp\n" if ($test);
				    print "minI2 = $i2\n" if ($test);
				}
			    }
			}		    
		    }
		}
	    }
	    if($minDisp < $minDispInit) {
		$sum += $minDisp;
		$total++;
		$USED{$minI2}++;	    
		print "minDisp = $minDisp < $minDispInit = minDispInit: now sum = sum + $minDisp = $sum\ntotal=$total\n" if ($test);
		print "minI2 = $minI2\n" if ($test);
	    }
	}
    }
    return($sum,$total);
}

#this version of getReverseProductDisplacement does not take into account product size
sub oldGetReverseProductDisplacement {
    my($productInfo,$revProductInfo,$adjTotalReads) = @_;
    my $test = 1;
    my $minDispInit = 100;
    my $minFrac = 0.0005;
    my $sum = 0;
    my $total = 0;
    my %USED;
    for(my $i1=0;$i1<@{$productInfo};$i1++) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$productInfo->[$i1]};
	print "side=$side\ttype=$newType\tadjProdCount=$adjProdCount\tadjMaxProdCount=$adjMaxProdCount\n" if ($test);
	print "relStart=$relStart\trelStop=$relStop\toffset=$offset\tgStart=$gStart\n" if ($test);
	my $start = $relStart;
	my $stop = $relStop;
	my $minDisp = $minDispInit;
	print "minDisp=$minDisp\n" if ($test);
	my $minI2;
	for(my $i2=0;$i2<@{$revProductInfo};$i2++) {
	    my($aType,$aNewType,$aProdList,$aAdjProdCount,$aAdjMaxProdCount,
	       $aRelStart,$aRelStop,$aOffset,$aGStart) = @{$revProductInfo->[$i2]};
	    print "aSide=$aType\taType=$aNewType\taAdjProdCount=$aAdjProdCount\taAdjMaxProdCount=$aAdjMaxProdCount\n" if ($test);
	    print "aRelStart=$aRelStart\taRelStop=$aRelStop\taOffset=$aOffset\taGStart=$aGStart\n" if ($test);
	    if($aAdjProdCount >= $minFrac*$adjTotalReads) {
		print "aAdjProdCount = $aAdjProdCount >= $minFrac * $adjTotalReads = minFrac x adjTotalReads\n" if ($test);
		my $aStart = $aRelStart;
		my $aStop = $aRelStop;
		# if the products overlap..
		if((($start <= $aStart)&&($aStart <= $stop))||
		   (($start <= $aStop)&&($aStop <= $stop))||
		   (($aStart <= $start)&&($start <= $aStop))||
		   (($aStart <= $stop)&&($stop <= $aStop))) {
		    my $disp = abs($aStart - $start);
		    print "products overlap: disp = $disp\n" if ($test);
		    unless($USED{$i2}) {
			if($disp < $minDisp) {
			    $minDisp = $disp;
			    $minI2 = $i2;
			    print "disp = $disp < $minDisp = minDisp: now minDisp = $disp\n" if ($test);
			    print "minI2 = $i2\n" if ($test);
			}
		    }
		}
	    }
	}
	if($minDisp < $minDispInit) {
	    $sum += $minDisp;
	    $total++;
	    $USED{$minI2}++;	    
	    print "minDisp = $minDisp < $minDispInit = minDispInit: now sum = sum + $minDisp = $sum\ntotal=$total\n" if ($test);
	    print "minI2 = $minI2\n" if ($test);
	}
    }
    return($sum,$total);
}

sub plausibleReads {
    my($adjTotalProductReads,$adjTotalRevProductReads,$productInfo,$revProductInfo,$parameters) = @_;
    #    #$tk->start("plausibleReads");
    if(doubleStrandedProducts($adjTotalProductReads, $adjTotalRevProductReads, $parameters)) {
	if(overlappingProducts($productInfo,$revProductInfo,$parameters)) {
	    #	    #$tk->end("plausibleReads");
	    return (1,"");
	} else {
 #    #$tk->end("plausibleReads");
	    return (0,"rejected: has abundant reads on each strand, but none overlap.");
	}
    } else {
	##$tk->end("plausibleReads");
	return (1,"");
    }
}

sub computeProductBasePairing {
    my($center,$productInfo,$basePairs,$parameters) = @_;
    my $minBasePairDensity = 1.0;
    my $foundMir = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop,$offset,$gStart) = @{$product};
	unless (($newType =~ /long/) || ($newType =~ /short/)) {
	    if($newType eq "miR") {	    
		$foundMir = 1;
		my $basePairDensity = computeBasePairDensity($basePairs,$relStart,$relStop);
		if($basePairDensity < $minBasePairDensity) {
		    $minBasePairDensity = $basePairDensity;
		}
	    }
	}
    }
    $minBasePairDensity = 0.0 unless ($foundMir);
    return $minBasePairDensity;    
}

#this version of computeProductBasePairing does not take into account product size
sub oldComputeProductBasePairing {
    my($center,$productInfo,$basePairs,$parameters) = @_;
    my $minBasePairDensity = 1.0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$product};
	if($newType eq "miR") {	    
	    my $basePairDensity = computeBasePairDensity($basePairs,$relStart,$relStop);
	    if($basePairDensity < $minBasePairDensity) {
		$minBasePairDensity = $basePairDensity;
	    }
	}
    }
    return $minBasePairDensity;    
}

sub goodProducts {
    my($productInfo,$parameters) = @_;
    #    #$tk->start("goodProducts");
    my $maxFivePrimeHet = $parameters->{fivePrimeHetMax} or die "FAIL: no maxFivePrimeHet loaded.\n";
    my $minLocusCount = $parameters->{countMinLocus} or die "FAIL: no minLocusCount loaded.\n";
    my $adjReadCount = getAdjTotalProductReads($productInfo);
    #print "readCount = $readCount\n";
    if($adjReadCount < $minLocusCount) {
	#	#$tk->end("goodProducts");
	return 0;
    }
     foreach my $product (@{$productInfo}) {
	 my($side,$newType,$prodList,$prodCount,$maxProdCount) = @{$product};
	 my $fivePrimeHet = computeFivePrimeHet($prodList);
	 #print "prodCount: $prodCount, 5' het: $fivePrimeHet\n";
	 if(($fivePrimeHet < $maxFivePrimeHet)&&($prodCount > 1)) {
	     #	    $tk->end("goodProducts");
	     return 1;
	 }
     }
 #    $tk->end("goodProducts");
    return 0;
}

sub computeMaxProdFivePrimeHet {
    my($productInfo,$parameters) = @_;
    my $maxCount = 0;
    my $maxProdFivePrimeHet = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount) = @{$product};
	unless (($newType =~ /long/) || ($newType =~ /short/)) {
	    if($newType eq "miR") {
		my $fivePrimeHet = computeFivePrimeHet($prodList);
		if($adjProdCount > $maxCount) {
		    $maxProdFivePrimeHet = $fivePrimeHet;
		    $maxCount = $adjMaxProdCount;
		}
	    }
	}
    }
    return $maxProdFivePrimeHet;
}

#this version of computeMaxProdFivePrimeHet does not take into account product size
sub oldComputeMaxProdFivePrimeHet {
    my($productInfo,$parameters) = @_;
    my $maxCount = 0;
    my $maxProdFivePrimeHet = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount) = @{$product};
	if($newType eq "miR") {
	    my $fivePrimeHet = computeFivePrimeHet($prodList);
	    if($adjProdCount > $maxCount) {
		$maxProdFivePrimeHet = $fivePrimeHet;
		$maxCount = $adjMaxProdCount;
	    }
	}
    }
    return $maxProdFivePrimeHet;
}

sub computeFivePrimeHet {
    my($productList) = @_;
    my $FIRST=1;
    my $fivePrimeMaxPos;
    my $adjFivePrimeMaxCount = 0;
    my $adjFivePrimeTotalCount = 0;
    my %startCount;
    foreach my $read (@{$productList}) {
	my($relStart,$relStop,$offset,$gStart,$count,$hitCount,$libraryCounts,$seq,$adjustedSeqCount,$adjustedLibraryCounts) = @{$read};
	$startCount{$relStart} += $count * $adjustedSeqCount;
    }
    my @starts = sort {$startCount{$b} <=> $startCount{$a}} keys(%startCount);
    my $topStart = shift(@starts);
    foreach my $read (@{$productList}) {
	my($relStart,$relStop,$offset,$gStart,$count,$hitCount,$libraryCounts,$seq,$adjustedSeqCount,$adjustedLibraryCounts) = @{$read};
	if($relStart == $topStart) {
	    $adjFivePrimeMaxCount += $count * $adjustedSeqCount;
	}
	$adjFivePrimeTotalCount += $count * $adjustedSeqCount;
    }
    my $fivePrimeHet = ($adjFivePrimeTotalCount-$adjFivePrimeMaxCount)/$adjFivePrimeTotalCount;
    return $fivePrimeHet;
}

sub otherSide {
    my $side = shift;
    if($side eq "5p") {
	return "3p";
    } elsif($side eq "3p") {
	return "5p";
    } else {
	die "unexpected side in otherSide: $side\n";
    }
}

sub overlapsMir {
    my($center,$basePairs,$side,$relStart,$relStop,$newProductInfo,$parameters) = @_;
    my $minShift = $parameters->{shiftMin} or die "FAIL: no minShift loaded.\n";
    foreach my $oProduct (@{$newProductInfo}) {
	my($oSide,$oNewType,$oProductlist,$oAdjustedProductCount,$oAdjustedMaxProdCount,$oRelStart,$oRelStop) = @{$oProduct};
	if($oNewType eq "miR") {
	    if(($side eq "5p")&&($oSide eq "3p")) {
		my($oRelStartBp,$oRelStopBp) = getOuterBasePairs($basePairs,$oRelStart,$oRelStop);		
		my $overlap = getOverlap($relStart,$relStop,$oRelStartBp,$oRelStopBp);
		my $shift = getShift($relStart,$relStop,$oRelStartBp,$oRelStopBp);
		#print "1:$side ($relStart,$relStop) ($oRelStartBp,$oRelStopBp) $shift\n";
		if($overlap > 0) { 
		    if(abs($shift) <= $minShift) {
			return 1;
		     }
		}
	    } elsif(($side eq "3p")&&($oSide eq "5p")) {
		my($relStartBp,$relStopBp) = getOuterBasePairs($basePairs,$relStart,$relStop);
		my $shift = getShift($oRelStart,$oRelStop,$relStartBp,$relStopBp);
		 my $overlap = getOverlap($oRelStart,$oRelStop,$relStartBp,$relStopBp);
		#print "2:$side ($relStartBp,$relStopBp) ($oRelStart,$oRelStop) $shift\n";
		if($overlap > 0) {
		    if(abs($shift) <= $minShift) {
			return 1;
		    }
		}
	    }
	 }
    }
    return 0;
}

sub withinHairpin {
    my($leftEnd,$rightEnd,$relStart,$relStop,$parameters) = @_;
    my $inBuffer = $parameters->{InHairpinBuffer};
    my $outBuffer = $parameters->{OutHairpinBuffer};   
    if(($leftEnd-$outBuffer<=$relStart)&&
       ($relStop<=$rightEnd+$outBuffer)) {	    
	return 1;
    }
    return 0;
}

sub onHairpinArm {
    my($center,$leftEnd,$rightEnd,$relStart,$relStop,$parameters) = @_;
    my($leftCenter, $rightCenter) = @{$center};
    my $inBuffer = $parameters->{InHairpinBuffer};
    my $outBuffer = $parameters->{OutHairpinBuffer};   
    if((($leftEnd-$outBuffer<=$relStart)&&($relStop<=$leftCenter+$inBuffer))||
       (($rightCenter-$inBuffer<=$relStart)&&($relStop<=$rightEnd+$outBuffer))) {	    
	return 1;
    }
    return 0;
}

sub closeToHairpin {
    my($center,$relStart,$relStop,$parameters) = @_;
    my($leftCenter, $rightCenter) = @{$center};
    my $hairpinRange = $parameters->{RangeOfHairpin} or die "no hairpinRange loaded\n";
    if(($leftCenter-$hairpinRange <= $relStart)&&
       ($relStop <= $rightCenter+$hairpinRange)) {
	return 1;
    }
    return 0;
}

sub computeMaxProdHitCount {
    my($productInfo,$location,$distinctReads,$parameters) = @_;
    my $minDist = $parameters->{distanceMin} or die "failed to load parameter distanceMin\n";
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $maxCount = 0;
    my $maxProdHitCount = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop,$offset,$gStart) = @{$product};
	unless (($newType =~ /long/) || ($newType =~ /short/)) {
	    if($newType eq "miR") {
		my $sum = 0;
		my $total = 0;
		foreach my $distinctRead (@{$distinctReads->{$strand}}) {
		    # for each read...
		    my($rStart,$rStop,$readCounts,$hitCount,$libraryCounts) = @{$distinctRead};
		    # if it is associated with a product...
		    if(abs($gStart-$rStart) < $minDist) {
			# compute the average number of hits to the genome
			$sum += $hitCount * $readCounts;
			$total += $readCounts;
			#print "incrementing: ( $rStart $rStop ) $rID $hitCount->{$sample.$rID}\n";
			#print "--> $sum $total\n";
		    }
		}
		my $avg = $total ? $sum/$total : 0.0;
		#print " = $avg\n";
		if($adjProdCount > $maxCount) {
		    $maxCount = $adjProdCount;
		    $maxProdHitCount = $avg;		
		}
	    }
	 }
    }
    return $maxProdHitCount;
}

#this version of computeMaxProdHitCount does not take into account product size
sub oldComputeMaxProdHitCount {
    my($productInfo,$location,$distinctReads,$parameters) = @_;
    my $minDist = $parameters->{distanceMin} or die "failed to load parameter distanceMin\n";
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $maxCount = 0;
    my $maxProdHitCount = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$product};
	 if($newType eq "miR") {
	     my $sum = 0;
	     my $total = 0;
	     foreach my $distinctRead (@{$distinctReads->{$strand}}) {
		 # for each read...
		 my($rStart,$rStop,$readCounts,$hitCount,$libraryCounts) = @{$distinctRead};
		 # if it is associated with a product...
		 if(abs($gStart-$rStart) < $minDist) {
		     # compute the average number of hits to the genome
		     $sum += $hitCount * $readCounts;
		     $total += $readCounts;
		     #print "incrementing: ( $rStart $rStop ) $rID $hitCount->{$sample.$rID}\n";
		     #print "--> $sum $total\n";
		 }
	     }
	     my $avg = $total ? $sum/$total : 0.0;
	     #print " = $avg\n";
	     if($adjProdCount > $maxCount) {
		 $maxCount = $adjProdCount;
		 $maxProdHitCount = $avg;		
	     }
	 }
    }
    return $maxProdHitCount;
}

sub getAdjTotalReads {
    my($distinctReads,$parameters) = @_;
    my $adjTotalReads = 0;
    my $shortSize = $parameters->{"MaxShortProductSize"};
    my $longSize = $parameters->{"minLongProductSize"};
    foreach my $dRead (@{$distinctReads}) {
	# dStart is genomic position.
	my($dStart,$dStop,$total,$hitCount,$libraryCounts) = @{$dRead};
	my $length = $dStop - $dStart + 1;
	if (($length > $shortSize) && ($length < $longSize)) {
	    $adjTotalReads += $total / $hitCount;
	}
    }
    return $adjTotalReads;
}

#this version of getAdjTotalReads does not take into account product size
sub oldGetAdjTotalReads {
    my($distinctReads) = @_;
    my $adjTotalReads = 0;
    foreach my $dRead (@{$distinctReads}) {
	# dStart is genomic position.
	my($dStart,$dStop,$total,$hitCount,$libraryCounts) = @{$dRead};
	$adjTotalReads += $total / $hitCount;
    }
    return $adjTotalReads;
}


sub getProductInfo {
    my($center,$fold,$basePairs,$location,$senseProducts,$parameters) = @_;
    my $test = 0;
    my ($chrom,$start,$stop,$strand) = parseLocation($location);
    my($leftCenter, $rightCenter) = @{$center};
    my @productInfo;
    my($leftEnd,$rightEnd) = getFullHairpinWindow($fold,$center);
    foreach my $product (@{$senseProducts}) {
	my($id,$productReads,$adjMaxProductCount,$adjProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$product};
	my @newProductReadsData;
	my $FIRST=1;
	 my $maxOffset;
	foreach my $read (@{$productReads}) {
	    my($dRelStart,$dRelStop,$dStart,$total,$hitCount,$libraryCounts,$seq,$adjustedSeqCount,$adjustedLibraryCounts,$rStrand) = @{$read};
	    print "dRelStart=$dRelStart\tdRelStop=$dRelStop\n" if ($test);
	    my $offset = $dRelStop - $leftCenter;
	     print "offset = $offset = $dRelStop - $leftCenter = dRelStaop - leftCenter\n" if ($test);
	    if(($offset > 0)&&($dRelStop < $rightCenter)) {
		$offset = 0;
		print "offset = $offset > 0 and drelStop = $dRelStop < $rightCenter = rightCenter\noffset set to 0" if ($test);
	    } elsif($offset > 0) {
		$offset = $dRelStart - $rightCenter;
		print "offset = $offset > 0\noffset now equals $offset = $dRelStart - $rightCenter = dRelStart - rightCenter\n" if ($test);
	    }	    
	    if (($dRelStart == $maxRelStart) && ($dRelStop == $maxRelStop)) {
		$maxOffset = $offset;
	    }
	    push(@newProductReadsData,[$dRelStart,$dRelStop,$offset,$dStart,$total,$hitCount,
				       $libraryCounts,$seq,$adjustedSeqCount,$adjustedLibraryCounts,$rStrand]);
	}
	my $side = getProductSide($center,$leftEnd,$rightEnd,
				  $maxRelStart,$maxRelStop,$parameters);
	push(@productInfo,[$side,\@newProductReadsData,$adjProductCount,$adjMaxProductCount,
			   $maxRelStart,$maxRelStop,$maxOffset,$maxGStart,$productStrand]);	
    }
    return \@productInfo;
}

sub rebuildProducts {
    my($center,$fold,$basePairs,$productInfo,$parameters) = @_;
    my($leftEnd,$rightEnd) = getFullHairpinWindow($fold,$center);
    my %usedMirSide;
    my %usedMorSide;
    my %newTypes;
    my @newProductInfo1;
    my @newProductInfo2;
    # sort by distance to loop, ascending order.
    my $totalReads = 0;
    foreach my $products (@{$productInfo}) {
	my($side,$productList,$adjProductCount,$adjMaxProductCount) = @{$products};
	$totalReads += $adjProductCount;
    }   
    my @sortedProductInfo = sort {abs($a->[6]) <=> abs($b->[6])} @{$productInfo};
    foreach my $products (@sortedProductInfo) {
	my($side,$productList,$adjProductCount,$adjMaxProductCount,
	   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$products};
	my $newType;
	#print "$side $relStart $relStop $offset\n";
	if($side eq "loop") {
	    $newType = "loop";
	} elsif($side eq "split") {
	    $newType = "split";
	} elsif($side eq "out") {
	    $newType = "out";
	} elsif(withinHairpin($leftEnd,$rightEnd,$relStart,$relStop,$parameters)) {
	    if(onHairpinArm($center,$leftEnd,$rightEnd,
			    $relStart,$relStop,$parameters)) {		
		if($usedMirSide{$side}) {
		    if($usedMorSide{$side}) {
			$newType = "out";
		    } else {
			$newType = "moR";
			$usedMorSide{$side}++;
		    }
		} elsif($usedMirSide{otherSide($side)}) {
		    if(overlapsMir($center,$basePairs,$side,$relStart,$relStop,\@newProductInfo1,$parameters)) {
			$newType = "miR";
			$usedMirSide{$side}++;
		    } else {
			if($usedMorSide{$side}) {
			    $newType = "out";
			} else {
			    $newType = "moR";
			    $usedMorSide{$side}++;
			}
		    }
		} else {
		     if($adjProductCount > 0.05*$totalReads) {
			 $newType = "miR";
			 $usedMirSide{$side}++;
		     } else {
			 $newType = "loop";
		     }
		}
	    } else {
		# this should be redundant, but here just in case.
		$newType = "loop";
	     }
	} elsif(closeToHairpin($center,$relStart,$relStop,$parameters)) {
	    if($usedMirSide{$side}) {
		if($usedMorSide{$side}) {
		    $newType = "out";
		} else {
		    $newType = "moR";
		    $usedMorSide{$side}++
		}
	    } elsif($usedMirSide{otherSide($side)}) {
		if(overlapsMir($center,$basePairs,
			       $side,$relStart,$relStop,\@newProductInfo1,$parameters)) {
		    $newType = "miR";
		    $usedMirSide{$side}++;
		} else {
		    if($usedMorSide{$side}) {
			$newType = "out";
		    } else {
			$newType = "moR";
			$usedMorSide{$side}++;
		    }
		 }
	    } else {
		$newType = "out";
	    }
	} else {
	    $newType = "out";
	}
	$newTypes{$relStart."x".$relStop} = $newType;
	push(@newProductInfo1,[$side,$newType,$productList,
			       $adjProductCount,$adjMaxProductCount,
			       $relStart,$relStop,$offset,$gStart,$productStrand]);
    }
    foreach my $products (@sortedProductInfo) {
	my($side,$productList,$adjProductCount,$adjMaxProductCount,
	   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$products};
	if(($side eq "5p")||($side eq "3p")) {
	    if($newTypes{$relStart."x".$relStop} eq "loop") {
		if(overlapsMir($center,$basePairs,$side,
			       $relStart,$relStop,\@newProductInfo1,$parameters)) {
		    $newTypes{$relStart."x".$relStop} = "miR";
		}
	    }
	}
	push(@newProductInfo2,[$side,$newTypes{$relStart."x".$relStop},$productList,
			       $adjProductCount,$adjMaxProductCount,
			       $relStart,$relStop,$offset,$gStart,$productStrand]);
    }
    return \@newProductInfo2;
}

sub addProductSize {
    my($products,$parameters) = @_;
    my $shortSize = $parameters->{"MaxShortProductSize"};
    my $longSize = $parameters->{"minLongProductSize"};
    my @reNamedProducts;
    
    foreach my $product (@{$products}) {
	my($side,$type,$productList,$adjProductcount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
	my $length = $relStop - $relStart + 1;
	if ($length >= $longSize) {
	    $product->[1] = "long-" . $product->[1]; #adding long to the begginging of product type
	} elsif ($length <= $shortSize) {
	    $product->[1] = "short-" . $product->[1]; #adding short to the begginging of product type
	}
	push(@reNamedProducts,$product);
    }
    return \@reNamedProducts;
}

sub overlapCurrent {
    my($productHash,$parsedRead,$parameters)=@_;
    my $minDist = $parameters->{distanceMin};
    my($readRelStart,$readRelStop,$gStart,$count,$hitCount,$libraryCounts) = @{$parsedRead};
    # sort by the first entries total count/hit count ( the most abundant entry by construction, since
    # distinctReads is already sorted by abundance. ) highest to lowest...
    my @sortedIdList = sort {$productHash->{$b}[0][3]/$productHash->{$b}[0][4] <=> $productHash->{$a}[0][3]/$productHash->{$a}[0][4]} keys(%{$productHash});
    my @idList;
    foreach my $id (@sortedIdList) {
	my($sRelStart,$sRelStop,$sGStart,$sCount,$hitCount)=@{$productHash->{$id}[0]};
	#print "comparing($id): ($readRelStart,$readRelStop) to ($sRelStart,$sRelStop) : ",abs($sRelStart-$readRelStart),"\n";
	if(abs($sRelStart-$readRelStart) <= $minDist) {
	    push(@idList,$id);
	    #print "ACCPETED -> $id\n";
	}
    }
    if(scalar(@idList) == 0) {
	# no overlaps, return FALSE.
        return "";
    } else {
	# more than one overlap, return the most abundant, first one by construction..
        return shift(@idList);
    }
}

sub getProductSide {
    my($center,$leftEnd,$rightEnd,$relStart,$relStop,$parameters) = @_;    
    my($leftCenter, $rightCenter) = @{$center};
    my $inBuffer = $parameters->{InHairpinBuffer} or die "no inHairpinBuffer loaded\n";
    my $outBuffer = $parameters->{OutHairpinBuffer} or die "no outHairpinBuffer loaded\n";   
    my $hairpinRange = $parameters->{RangeOfHairpin} or die "no hairpinRange loaded";
    if(withinHairpin($leftEnd,$rightEnd,$relStart,$relStop,$parameters)) {
	if(onHairpinArm($center,$leftEnd,$rightEnd,$relStart,$relStop,$parameters)) {
	    if($relStop <= $leftCenter+$inBuffer) {
		return "5p";
	    } elsif($rightCenter-$inBuffer<=$relStart) {	    
		return "3p";	
	    } else {
		die "unknow situation reached in getProductSide!!! ouch.\n";
	    }
	} elsif(($relStart<=$leftCenter+1)&&($rightCenter-1<=$relStop)) {
#	    print "loop named here: relStart = $relStart leftCenter = $leftCenter rightCenter = $rightCenter relStop = $relStop\n";
	    return "loop";
	} elsif(($leftCenter<=$relStart)&&($relStop<=$rightCenter)) {
	    return "loop";
	} else {
	    return "split";
	}
    } elsif(($leftCenter-$hairpinRange <= $relStart)&&
	    ($relStop<=$rightCenter+$hairpinRange)) {
	if(($leftCenter-$hairpinRange <= $relStart)&&
	   ($relStop <= $leftCenter+$inBuffer)) {
	    return "5p";
	} elsif(($rightCenter-$inBuffer<=$relStart)&&
		($relStop<=$rightCenter+$hairpinRange)) {	    
	    return "3p";
	}	
    } else {
	return "out";
    }
    return "out";
}

sub getAdjTotalProductReads {
    my $productInfo = shift;
    my $adjTotalReads = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,$relStart,$relStop,$offset,$gStart) = @{$product};
	unless (($newType =~ /long/) || ($newType =~ /short/)) { 
	    $adjTotalReads += $adjProdCount;
	}
    }
    return $adjTotalReads;
}

#this version of getAdjTotalProductReads does not take into account product size
sub oldGetAdjTotalProductReads {
    my $productInfo = shift;
    my $adjTotalReads = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$adjProdCount,$adjMaxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$product};
	$adjTotalReads += $adjProdCount;
    }
    return $adjTotalReads;
}

sub getMaxProductInfo {
    my $productList = shift;
    my %dReadCounts;
    my $adjProductCount = 0;
    my $adjMaxProductCount = 0;
    my($maxRelStart,$maxRelStop,$maxOffset,$maxGStart);
    foreach my $read (@{$productList}) {
        my($rRelStart,$rRelStop,$rOffset,$rGStart,$rCount,$rHitCount,$libraryCounts,$seq,$adjustedSeqCount)=@{$read};
	$dReadCounts{$rRelStart."x".$rRelStop} += $rCount * $adjustedSeqCount;
    }
    my @keys = sort {$dReadCounts{$b} <=> $dReadCounts{$a}} keys %dReadCounts;
    my $maxKey = shift(@keys);
    ($maxRelStart,$maxRelStop) = split(/x/,$maxKey); 
    foreach my $read (@{$productList}) {
        my($rRelStart,$rRelStop,$rOffset,$rGStart,$rCount,$rHitCount,$libraryCounts,$seq,$adjustedSeqCount)=@{$read};
	$adjProductCount += $rCount * $adjustedSeqCount;
	if(($rRelStart == $maxRelStart)&&($rRelStop == $maxRelStop)) {
	    $adjMaxProductCount += $rCount * $adjustedSeqCount;
	    $maxOffset = $rOffset;
	    $maxGStart = $rGStart;	    
	}
    }
    return($maxRelStart,$maxRelStop,$maxOffset,$maxGStart,$adjMaxProductCount,$adjProductCount);
}

sub getAdjTotalLibraryCounts {
    my $productList = shift;
    my %adjTotalLibraryCounts;
    foreach my $read (@{$productList}) {
	my($rRelStart,$rRelStop,$rOffset,$rGStart,$rCount,$rHitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts)=@{$read};
	foreach my $sample (keys %{$adjLibraryCounts}) {
	    $adjTotalLibraryCounts{$sample} += $adjLibraryCounts->{$sample};
	}
    }
    return \%adjTotalLibraryCounts;
}

########################################
#  General Machine Learning Functions  #
########################################

sub evaluateFeatureVectorFile_new {
    my($R,$trainFile,$featureVectorFile,$outputFile,$scoresToRemove,$parameters) = @_;
    my $out = $R->set('featureVecTrainFile', $trainFile);
    $out .= $R->set('featureVectorFile', $featureVectorFile);
    $out .= $R->set('outputFile', $outputFile);
    $out .= $R->run(q`unlink(outputFile)`);
    $out .= $R->run(q`print(paste("reading in ",featureVecTrainFile, " for the training data", sep=""))`);
    $out .= $R->run(q`print(paste("reading in ",featureVectorFile, " for the testing data", sep=""))`);
    $out .= $R->run(q`print(paste("reading in ",outputFile, " for the output of the training data", sep=""))`);
    # Loading Libraries and Parameters 
    $out .= $R->run(q`library(ROCR)`);
    $out .= $R->run(q`library(randomForest)`);
    $out .= $R->run(q`classificationtype <-'C-classification'`);
    $out .= $R->run(q`crossvalidation <- 10`);
    $out .= $R->run(q`options(scipen=999)`);
    #loading train data
    $out .= $R->run(q`trainData <- subset(read.table(featureVecTrainFile, header=TRUE, sep="\t"), select=c(-GeneId))`);
    foreach my $score (@{$scoresToRemove}) {
	my $command = "if (\"". $score . "\" %in% colnames(trainData)) trainData\$" . $score. " <- NULL";
	$out .= $R->run($command);
	$R->run(q`columnNames <- colnames(trainData)`);
#	my $names = $R->get( 'columnNames' );
#	foreach my $name (@{$names}) {
#	    print "$name\t";
#	}
#	print "\n";
    }
    $out .= $R->run(q`trainClass <- as.matrix(subset(trainData, select = y))`);
    $out .= $R->run(q`trainValues <- as.matrix(subset(trainData, select = -y))`);
    #loading test data
    $out .= $R->run(q`testRows <- read.table(featureVectorFile,header=TRUE)`);
    $out .= $R->run(q`num <- nrow(testRows)`);
    $out .= $R->run(q`trainIndicies <- c(1:num)`);
    $out .= $R->run(q`testData <- subset(testRows, select = c(-GeneId))`);
    foreach my $score (@{$scoresToRemove}) {
	my $command = "if (\"". $score . "\" %in% colnames(testData)) testData\$" . $score. " <- NULL";
	$out .= $R->run($command);
	$R->run(q`columnNames <- colnames(testData)`);
#	my $names = $R->get( 'columnNames' );
#	foreach my $name (@{$names}) {
#	    print "$name\t";
#	}
#	print "\n";
    }
    $out .= $R->run(q`testClass <- as.matrix(subset(testData, select = y))`);
    $out .= $R->run(q`testValues<- as.matrix(subset(testData, select = c(-y)))`);
    #creating random forest from train data
    $out .= $R->run(q`model <-rfcv(trainy=trainClass,trainx=trainValues,scale = "log" ,cv.fold=crossvalidation)`);
    $out .= $R->run(q`iMIN = which.min(model$error.cv)`);
    $out .= $R->run(q`predicted <-model$predicted[[iMIN]]`);
    $out .= $R->run(q`predict.rf <-randomForest(as.factor(trainClass) ~ .,data=trainValues, mtry=iMIN, ntree=1000,keep.forest=TRUE,importance=TRUE,test=data$val)`);
    $out .= $R->run(q`predict.rf$predicted = predicted[,"y"]`);
    $out .= $R->run(q`predict.rf.pr = predict(predict.rf,type="prob",newdata=testValues)[,2]`);
    $out .= $R->run(q`classPredictions <- as.matrix(predict(predict.rf, testValues, decision.values=TRUE))`);
    $out .= $R->run(q`testedGeneIds <- subset(testRows[attributes(predict.rf.pr)$names,], select = GeneId)`);
    $out .= $R->run(q`testedProbs <- as.matrix(predict.rf.pr)`);
#    $out .= $R->run(q`testAdjMaxCount <- subset(testRows[attributes(predict.rf.pr)$names,], select = adjProductCount)`);
#    $out .= $R->run(q`testOutput <- data.frame(testedGeneIds,classPredictions,testedProbs,testAdjMaxCount)`);
#    $out .= $R->run(q`write(t(testOutput),outputFile,ncolumns=4,sep="\t")`);
    $out .= $R->run(q`testOutput <- data.frame(testedGeneIds,classPredictions,testedProbs)`);
    $out .= $R->run(q`write(t(testOutput),outputFile,ncolumns=3,sep="\t")`);
    print "$out\n";
}


sub evaluateFeatureVectorFile {
    my($trainFile,$featureVectorFile,$outputFile,$parameters) = @_;
    my $R = Statistics::R->new();
    $R->start();
    my $out = $R->set('featureVecTrainFile', $trainFile);
    $out .= $R->set('featureVectorFile', $featureVectorFile);
    $out .= $R->set('outputFile', $outputFile);
    $out .= $R->run(q`unlink(outputFile)`);
    $out .= $R->run(q`print(paste("reading in ",featureVecTrainFile, " for the training data", sep=""))`);
    $out .= $R->run(q`print(paste("reading in ",featureVectorFile, " for the testing data", sep=""))`);
    $out .= $R->run(q`print(paste("reading in ",outputFile, " for the output of the training data", sep=""))`);
    # Loading Libraries and Parameters 
    $out .= $R->run(q`library(ROCR)`);
    $out .= $R->run(q`library(randomForest)`);
    $out .= $R->run(q`classificationtype <-'C-classification'`);
    $out .= $R->run(q`crossvalidation <- 10`);
    $out .= $R->run(q`options(scipen=999)`);
    #loading train data
    $out .= $R->run(q`trainData <- subset(read.table(featureVecTrainFile, header=TRUE, sep="\t"), select=c(-GeneId))`);
    $out .= $R->run(q`trainClass <- as.matrix(subset(trainData, select = y))`);
    $out .= $R->run(q`trainValues <- as.matrix(subset(trainData, select = -y))`);
    #loading test data
    $out .= $R->run(q`testRows <- read.table(featureVectorFile,header=TRUE)`);
    $out .= $R->run(q`num <- nrow(testRows)`);
    $out .= $R->run(q`trainIndicies <- c(1:num)`);
    $out .= $R->run(q`testData <- subset(testRows, select = c(-GeneId))`);
    $out .= $R->run(q`testClass <- as.matrix(subset(testData, select = y))`);
    $out .= $R->run(q`testValues<- as.matrix(subset(testData, select = c(-y)))`);
    #creating random forest from train data
    $out .= $R->run(q`model <-rfcv(trainy=trainClass,trainx=trainValues,scale = "log" ,cv.fold=crossvalidation)`);
    $out .= $R->run(q`iMIN = which.min(model$error.cv)`);
    $out .= $R->run(q`predicted <-model$predicted[[iMIN]]`);
    $out .= $R->run(q`predict.rf <-randomForest(as.factor(trainClass) ~ .,data=trainValues, mtry=iMIN, ntree=1000,keep.forest=TRUE,importance=TRUE,test=data$val)`);
    $out .= $R->run(q`predict.rf$predicted = predicted[,"y"]`);
    $out .= $R->run(q`predict.rf.pr = predict(predict.rf,type="prob",newdata=testValues)[,2]`);
    $out .= $R->run(q`classPredictions <- as.matrix(predict(predict.rf, testValues, decision.values=TRUE))`);
    $out .= $R->run(q`testedGeneIds <- subset(testRows[attributes(predict.rf.pr)$names,], select = GeneId)`);
    $out .= $R->run(q`testedProbs <- as.matrix(predict.rf.pr)`);
#    $out .= $R->run(q`testAdjMaxCount <- subset(testRows[attributes(predict.rf.pr)$names,], select = adjProductCount)`);
#    $out .= $R->run(q`testOutput <- data.frame(testedGeneIds,classPredictions,testedProbs,testAdjMaxCount)`);
#    $out .= $R->run(q`write(t(testOutput),outputFile,ncolumns=4,sep="\t")`);
    $out .= $R->run(q`testOutput <- data.frame(testedGeneIds,classPredictions,testedProbs)`);
    $out .= $R->run(q`write(t(testOutput),outputFile,ncolumns=3,sep="\t")`);
    $R->stop();
    print "$out\n";
}


sub testTrainDataRF {
    my($featureVectorFile,$outputFile,$parameters) = @_;
    my $R = Statistics::R->new();
    my $out = $R->set('featureVecTrainFile', $featureVectorFile);
    $out = $R->set('outputFile', $outputFile);
    # Loading Libraries and Parameters 
    $out .= $R->run(q`library(ROCR)`);
    $out .= $R->run(q`library(randomForest)`);
    $out .= $R->run(q`classificationtype <-'C-classification'`);
    $out .= $R->run(q`crossvalidation <- 10`);
    $out .= $R->run(q`options(scipen=999)`);
    $out .= $R->run(q`unlink(outputFile)`);
    #reading in training file
    $out .= $R->run(q`allRows <- read.table(featureVecTrainFile,header=TRUE)`);
    $out .= $R->run(q`num <- nrow(allRows)`);
    $out .= $R->run(q`allIndicies <- c(1:num)`);
    #testing predictions of a model trained from one half of data (testing on other half)
    $out .= $R->run(q`trainIndicies <- sort(sample(allIndicies, size=num/2,replace=FALSE))`);
    $out .= $R->run(q`testIndicies <- setdiff(allIndicies,trainIndicies)`);
    $out .= $R->run(q`trainRows <- allRows[trainIndicies,]`);
    $out .= $R->run(q`trainData <- subset(trainRows, select = c(-GeneId))`);
    $out .= $R->run(q`trainClass <- as.matrix(subset(trainData, select = y))`);
    $out .= $R->run(q`trainValues <- as.matrix(subset(trainData, select = c(-y)))`);
    $out .= $R->run(q`testRows <- allRows[testIndicies,]`);
    $out .= $R->run(q`testData <- subset(testRows, select = c(-GeneId))`);
    $out .= $R->run(q`testClass<- as.matrix(subset(testData, select = y))`);
    $out .= $R->run(q`testValues<- as.matrix(subset(testData, select = c(-y)))`);
    $out .= $R->run(q`model <-rfcv(trainy=trainClass,trainx=trainValues,scale = "log" ,cv.fold=crossvalidation)`);
    $out .= $R->run(q`iMIN = which.min(model$error.cv)`);
    $out .= $R->run(q`predicted <-model$predicted[[iMIN]]`);
    $out .= $R->run(q`predict.rf <-randomForest(as.factor(trainClass) ~ .,data=trainValues, mtry=iMIN, ntree=1000,keep.forest=TRUE,importance=TRUE,test=data$val)`);
    $out .= $R->run(q`predict.rf$predicted = predicted[,"y"]`);
    $out .= $R->run(q`predict.rf.pr = predict(predict.rf,type="prob",newdata=testValues)[,2]`);
    $out .= $R->run(q`predict.rf.pred = prediction(predict.rf.pr, testClass)`);
    $out .= $R->run(q`predict.rf.perf = performance(predict.rf.pred,"tpr","fpr")`);
    $out .= $R->run(q`aucValue1<-round(performance(predict.rf.pred, 'auc')@y.values[[1]],4)`);
    $out .= $R->run(q`classPredictions <- as.matrix(predict(predict.rf, testValues, decision.values=TRUE))`);
    $out .= $R->run(q`testedGeneIds <- subset(allRows[attributes(predict.rf.pr)$names,], select = GeneId)`);
    $out .= $R->run(q`testedGeneClass <- subset(allRows[attributes(predict.rf.pr)$names,], select = y)`);
    $out .= $R->run(q`testedProbs <- as.matrix(predict.rf.pr)`);
    #$out .= $R->run(q`testedAdjMaxProdCount <- subset(allRows[attributes(predict.rf.pr)$names,], select = adjMaxProductCount)`);
    #$out .= $R->run(q`testedAdjMaxCount <- subset(allRows[attributes(predict.rf.pr)$names,], select = adjProductCount)`);
    #$out .= $R->run(q`test1Output <- data.frame(testedGeneIds,testedGeneClass,classPredictions,testedProbs,testedAdjMaxProdCount,testedAdjMaxCount)`);
    $out .= $R->run(q`test1Output <- data.frame(testedGeneIds,testedGeneClass,classPredictions,testedProbs)`);
    $out .= $R->run(q`write(t(test1Output),outputFile,ncolumns=4,sep="\t")`);
    #testing predictions of a model trained on the other half of the sample
    $out .= $R->run(q`temp <- trainIndicies`);
    $out .= $R->run(q`trainIndicies <- testIndicies`);
    $out .= $R->run(q`testIndicies <- temp`);
    $out .= $R->run(q`trainRows <- allRows[trainIndicies,]`);
    $out .= $R->run(q`trainData <- subset(trainRows, select = c(-GeneId))`);
    $out .= $R->run(q`trainClass <- as.matrix(subset(trainData, select = y))`);
    $out .= $R->run(q`trainValues <- as.matrix(subset(trainData, select = c(-y)))`);
    $out .= $R->run(q`testRows <- allRows[testIndicies,]`);
    $out .= $R->run(q`testData <- subset(testRows, select = c(-GeneId))`);
    $out .= $R->run(q`testClass<- as.matrix(subset(testData, select = y))`);
    $out .= $R->run(q`trainValues <- as.matrix(subset(trainData, select = c(-y)))`);
    $out .= $R->run(q`testRows <- allRows[testIndicies,]`);
    $out .= $R->run(q`testData <- subset(testRows, select = c(-GeneId))`);
    $out .= $R->run(q`testClass<- as.matrix(subset(testData, select = y))`);
    $out .= $R->run(q`testValues<- as.matrix(subset(testData, select = c(-y)))`);
    $out .= $R->run(q`model <-rfcv(trainy=trainClass,trainx=trainValues,scale = "log" ,cv.fold=crossvalidation)`);
    $out .= $R->run(q`iMIN = which.min(model$error.cv)`);
    $out .= $R->run(q`predicted <-model$predicted[[iMIN]]`);
    $out .= $R->run(q`predict.rf <-randomForest(as.factor(trainClass) ~ .,data=trainValues, mtry=iMIN, ntree=1000,keep.forest=TRUE,importance=TRUE,test=data$val)`);
    $out .= $R->run(q`predict.rf$predicted=predicted[,"y"]`);
    $out .= $R->run(q`predict.rf.pr = predict(predict.rf,type="prob",newdata=testValues)[,2]`);
    $out .= $R->run(q`predict.rf.pred = prediction(predict.rf.pr, testClass)`);
    $out .= $R->run(q`predict.rf.perf = performance(predict.rf.pred,"tpr","fpr")`);
    $out .= $R->run(q`aucValue2<-round(performance(predict.rf.pred, 'auc')@y.values[[1]],4)`);
    $out .= $R->run(q`classPredictions <- as.matrix(predict(predict.rf, testValues, decision.values=TRUE))`);
    $out .= $R->run(q`testedGeneIds <- subset(allRows[attributes(predict.rf.pr)$names,], select = GeneId)`);
    $out .= $R->run(q`testedGeneClass <- subset(allRows[attributes(predict.rf.pr)$names,], select = y)`);
    $out .= $R->run(q`testedProbs <- as.matrix(predict.rf.pr)`);
    #$out .= $R->run(q`testedAdjMaxProdCount <- subset(allRows[attributes(predict.rf.pr)$names,], select = adjMaxProductCount)`);
    #$out .= $R->run(q`testedAdjMaxCount <- subset(allRows[attributes(predict.rf.pr)$names,], select = adjProductCount)`);
    #$out .= $R->run(q`test2Output <- data.frame(testedGeneIds,testedGeneClass,classPredictions,testedProbs,testedAdjMaxProdCount,testedAdjMaxCount)`);
    $out .= $R->run(q`test2Output <- data.frame(testedGeneIds,testedGeneClass,classPredictions,testedProbs)`);
    $out .= $R->run(q`write(t(test2Output),outputFile,ncolumns=4,sep="\t",append=TRUE)`);
    $R->stop();
    print "$out\n";
}

sub generateRandForestAucFile {
    my($featureVecTrainingFile,$aucFile,$trials) = @_;
    my $R = Statistics::R->new();
    my $out = $R->set('featureFile',$featureVecTrainingFile);
    $out = $R->set('aucFile',$aucFile);
    $out = $R->run(q`unlink(aucFile)`);
    $out = $R->run(q`library(ROCR)`);
    $out = $R->run(q`library(randomForest)`);
    $out = $R->run(q`classificationtype <-'C-classification'`);
    $out = $R->run(q`crossvalidation <- 10`);
    $out = $R->run(q`trials <- 50`);
    $out = $R->run(q`geneIds <- subset(read.table(featureFile, header=TRUE, sep="\t"), select=c(GeneId))`);
    $out = $R->run(q`num = nrow(geneIds)`);
    for (my $i = 0; $i < $trials; $i++) {
	$out = $R->run(q`trainIndices<-sort(sample(c(1:num),size=2*num/3,replace=FALSE))`);
	$out = $R->run(q`testIndices<-setdiff(c(1:num),trainIndices)`);
	$out = $R->run(q`auc<-NULL`);
	$out = $R->run(q`data<-subset(read.table(featureFile, header=TRUE, sep="\t"), select=c(-GeneId))`);
	$out = $R->run(q`trainData<-data[trainIndices,]`);
	$out = $R->run(q`trainClass<-as.matrix(subset(trainData, select = y))`);
	$out = $R->run(q`trainValues<- as.matrix(subset(trainData, select = c(-y)))`);
	$out = $R->run(q`testData<-data[testIndices,]`);
	$out = $R->run(q`testClass<-as.matrix(subset(testData, select = y))`);
	$out = $R->run(q`testValues<- as.matrix(subset(testData, select = c(-y)))`);
	$out = $R->run(q`model <-rfcv(trainy=trainClass,trainx=trainValues,scale = "log" ,cv.fold=crossvalidation)`);
	$out = $R->run(q`iMIN = which.min(model$error.cv)`);
	$out = $R->run(q`predicted <-model$predicted[[iMIN]]`);
	$out = $R->run(q`predict.rf <-randomForest(as.factor(trainClass) ~ .,data=trainValues, mtry=iMIN, ntree=1000,keep.forest=TRUE,importance=TRUE,test=data$val)`);
	$out = $R->run(q`predict.rf$predicted=predicted[,"y"]`);
	$out = $R->run(q`predict.rf.pr = predict(predict.rf,type="prob",newdata=testValues)[,2]`);
	$out = $R->run(q`predict.rf.pred = prediction(predict.rf.pr, testClass)`);
	$out = $R->run(q`predict.rf.perf = performance(predict.rf.pred,"tpr","fpr")`);
	$out = $R->run(q`aucValue<-round(performance(predict.rf.pred, 'auc')@y.values[[1]],4)`);    
	$out = $R->run(q`write(aucValue, aucFile, ncolumns=length(aucValue), sep="\t", append=TRUE)`);
    }
    $R->stop();
    print "$out\n";
}

sub addRandomNegFeatures {
    #adds random negative features from feature vector file to the training data
    my($trainingData,$featureFile) = @_;
    my $R = Statistics::R->new();
    my $out = $R->set('posVectorFile',$trainingData);
    $out = $R->set('vectorFile',$featureFile);
    $out = $R->run(q`print(paste("reading in ",posVectorFile, " for the poitive Vector File", sep=""))`);
    $out = $R->run(q`print(paste("reading in ",vectorFile, " for the negative Vector File", sep=""))`);
    $out = $R->run(q`options(scipen=999)`);
    $out = $R->run(q`posTrainData <- read.table(posVectorFile, header=TRUE, sep="\t")`);
    $out = $R->run(q`Data <- read.table(vectorFile, header=TRUE, sep="\t")`);
    $out = $R->run(q`colNames <- as.matrix(names(posTrainData))`);
    $out = $R->run(q`numPos <- nrow(posTrainData)`);
    $out = $R->run(q`numNeg <- nrow(Data)`);
    $out = $R->run(q`numCols <- nrow(colNames)`);
    $out = $R->run(q`posIndicies <- c(1:numPos)`);
    $out = $R->run(q`sampledNegativeIndicies <- sort(sample(c(1:numNeg),size=numPos,replace=FALSE),decreasing=FALSE)`);
    $out = $R->run(q`negTrainData<-as.matrix(Data[sampledNegativeIndicies,])`);
    $out = $R->run(q`write(t(negTrainData),posVectorFile,ncolumns=numCols,sep="\t",append=TRUE)`);
    $R->stop();
}

#############################################
# MACHINE LEARNING PRODUCT DATA SUBROUTINES #
#############################################

sub getClassValues {
    my($parameters) = @_;
    my $MLType = $parameters->{MLAlgorithm};
    my($posClass,$negClass);
    if ($MLType eq 'RandomForest') {
	($posClass,$negClass) = (1,-1);
    } elsif ($MLType eq 'LogisticRegression') {
	($posClass,$negClass) = (1,0);	
    } else {
	die "unknown machine learning algorithm";
    }
    return($posClass,$negClass);
}

sub oldnearAnnotatedMir {
    my($location,$locationList,$cutoffDistance) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    foreach my $mirLocation (keys %{$locationList}) {
	my($mirChrom,$mirStart,$mirStop,$mirStrand) = parseLocation($mirLocation);
	my $distance = abs($start - $mirStart);
	if ($distance <= $cutoffDistance) {
	    return 1;
	}
    }
    return 0;
}

sub nearAnnotatedMir {
    my($location,$locationList,$cutoffDistance) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    foreach my $mirLocation (@{$locationList->{$chrom}}) {
	my($mirChrom,$mirStart,$mirStop,$mirStrand) = parseLocation($mirLocation);
	my $distance = abs($start - $mirStart);    #intended to be fast rather than absolutelty correct
	if ($distance <= $cutoffDistance) {
	    return 1;
	}
    }
    return 0;
}

sub getRandomProduct {
    my($products) = @_;
    my $numProducts = @{$products};
    my $randProductNum = int(rand($numProducts));
    my $randProduct = $products->[$randProductNum];
    return $randProduct;
}

sub getGreatestOverlappingProduct {
    #gets the greatest overlapping product with the same 5' end
    my($products,$location,$mirbaseId) = @_;
    my $printWarning = 1;
    my $minOverlap = 0;
    my $mirProductCount = 0;
    my $mirGStart;
    my @greatestOverlappingProduct;
    my($chrom,$mirbaseStart,$mirbaseStop,$strand) = parseLocation($location);
    my $length = $mirbaseStop - $mirbaseStart + 1;
    my $mirbaseRelStart = 0;
    my $mirbaseRelStop = $length - 1;
    foreach my $productInfo (@{$products}) {
	my($id,$productReads,$adjMaxProductCount,$adjProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$productInfo};
	my $overlap = getOverlap($mirbaseRelStart,$mirbaseRelStop,$maxRelStart,$maxRelStop);
	if ($overlap > $minOverlap) {
	    $minOverlap = $overlap;
	    @greatestOverlappingProduct = @{$productInfo};
	    $mirProductCount = $adjProductCount;
	    $mirGStart = $maxGStart;
	    $printWarning  = ($mirbaseRelStart == $maxRelStart) ? 0 : 1;
	}
    }
    if ($printWarning) {
	print "mir $mirbaseId at $location has reads that overlap the mir but don't line up on the 5' end (count = $mirProductCount)\n";
	print "start = $mirGStart\tannotatedStart = $mirbaseStart\n";
    }

    my $percentOverlap = $minOverlap / $length;

    return(\@greatestOverlappingProduct,$percentOverlap);
}

sub combineProductLocations {
    my($posProductLocations,$negProductLocations) = @_;
    my %trainingProductLocations;
    foreach my $chrom (keys %{$posProductLocations}) {
	foreach my $location (%{$posProductLocations->{$chrom}}) {
	    $trainingProductLocations{$chrom}{$location} = 1;
	}
    }
    foreach my $chrom (keys %{$negProductLocations}) {
	foreach my $location (%{$negProductLocations->{$chrom}}) {
	    $trainingProductLocations{$chrom}{$location} = 1;
	}
    }
    return \%trainingProductLocations;
}

sub getGenomicLocation {  #returns the genomic start and stop positions of a product
    my($product,$chrom) = @_;
    my($id,$reads,$adjMaxProductCount,$adjProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$product};
    my $location = "$chrom:$maxGStart.." . ($maxGStart + $maxRelStop - $maxRelStart) . ":$productStrand";
    return $location;
}

sub getProductMedianLength {
    my($reads) = @_;
    my @lengths;
    my $total;
    foreach my $read (@{$reads}) {
	my($relStart,$relStop,$dStart,$readCount,$hitCount,$libraryCounts,$seq,$adjustedSeqCount) = @{$read};
	my $length = length($seq);
	my $numReads = $readCount * $adjustedSeqCount;
	push(@lengths, [$length, $numReads]);
	$total += $numReads;
    }
    my @sortedLengths = sort {$a->[0] <=> $b->[0]} @lengths;
    my $medianLocation = $total / 2;
    my $numChecked = 0;
    my $i = 0;
    do {
	my($length,$count) = @{$sortedLengths[$i]};
	$numChecked += $count;
	$i++ unless ($i == @sortedLengths - 1);  #if we are already at the last index incrimenting it will cause mediumLength to be uninitialized
    } while ($numChecked < $medianLocation);
    my $medianLength = $sortedLengths[$i][0];
    return $medianLength;
}

sub factorial {
    my($num) = @_;
    if ($num == 0) {
	return 1;
    } 	
    my $fact = 1;
    for (my $i = 1; $i <= $num; $i++) {
	$fact *= $i;
    }
    return $fact;
}

sub getWFC {
    #determines WFComplexity
    my($seq) = @_;
    my $length = length($seq);
    my $aCount = $seq =~ s/(a|A)/$1/g;
    my $cCount = $seq =~ s/(c|C)/$1/g;
    my $gCount = $seq =~ s/(g|G)/$1/g;
    my $tCount = $seq =~ s/(t|T)/$1/g;
    my $nFactorial = factorial($length);
    my $aFactorial = factorial($aCount);
    my $cFactorial = factorial($cCount);
    my $gFactorial = factorial($gCount);
    my $tFactorial = factorial($tCount);
    my $mulNucFactorials = $aFactorial * $cFactorial * $gFactorial * $tFactorial;
    my $WFCquotient = $nFactorial / $mulNucFactorials;
    my $WFClog = logK($WFCquotient,4);
    my $WFC = $WFClog / $length;
    return $WFC;
}

sub getWFC_Stirling {
    #determines WFComplexity for longer sequences using the sterling approximation
    my($seq) = @_;
    my $length = length($seq);
    my $aCount = $seq =~ s/(a|A)/$1/g;
    my $cCount = $seq =~ s/(c|C)/$1/g;
    my $gCount = $seq =~ s/(g|G)/$1/g;
    my $tCount = $seq =~ s/(t|T)/$1/g;
    my $WFCnum = ($length * log($length) - $aCount * log($aCount) - $cCount * log($cCount) - $gCount * log($gCount) - $tCount * log($tCount) - $length + $aCount + $cCount + $gCount + $tCount) * logK(exp(1),4);
    my $WFC = $WFCnum / $length;
    return $WFC;
}

sub getWFC_Dinuc {
    #determines dinucleotide WFComplexity for longer sequences using the sterling approximation
    my($seq) = @_;
    my $length = length($seq);
    my $aa = $seq =~ s/(aa|aA|Aa|AA)/$1/g;
    my $ac = $seq =~ s/(ac|aC|Ac|AC)/$1/g;
    my $ag = $seq =~ s/(ag|aG|Ag|AG)/$1/g;
    my $at = $seq =~ s/(at|aT|At|AT)/$1/g;
    my $ca = $seq =~ s/(ca|cA|Ca|CA)/$1/g;
    my $cc = $seq =~ s/(cc|cC|Cc|CC)/$1/g;
    my $cg = $seq =~ s/(cg|cG|Cg|CG)/$1/g;
    my $ct = $seq =~ s/(ct|cT|Ct|CT)/$1/g;
    my $ga = $seq =~ s/(ga|gA|Ga|GA)/$1/g;
    my $gc = $seq =~ s/(gc|gC|Gc|GC)/$1/g;
    my $gg = $seq =~ s/(gg|gG|Gg|GG)/$1/g;
    my $gt = $seq =~ s/(gt|gT|Gt|GT)/$1/g;
    my $ta = $seq =~ s/(ta|tA|Ta|TA)/$1/g;
    my $tc = $seq =~ s/(tc|tC|Tc|TC)/$1/g;
    my $tg = $seq =~ s/(tg|tG|Tg|TG)/$1/g;
    my $tt = $seq =~ s/(tt|tT|Tt|TT)/$1/g;
    my $WFCnum = ($length * log($length) - $aa * log($aa) - $ac * log($ac) - $ag * log($ag) - $at * log($at) - $ca * log($ca) - $cc * log($cc) - $cg * log($cg) - $ct * log($ct) - $ga * log($ga) - $gc * log($gc) - $gg * log($gg) - $gt * log($gt) - $ta * log($ta) - $tc * log($tc) - $tg * log($tg) - $tt * log($tt) - $length + $aa + $ac + $ag + $at + $ca + $cc + $cg + $ct + $ga + $gc + $gg + $gt + $ta + $tc + $tg + $tt) * logK(exp(1),16);
    my $WFC = $WFCnum / $length;
    return $WFC;
}



sub getProductDuplexEnergy {
    my ($genome,$sequence,$rrLocation,$chromLengths,$parameters) = @_;
    my ($chrom,$start,$stop,$strand) = parseLocation($rrLocation);
    my $bufferedLocation = extendSequenceToFillBuffer($chrom,$start,$stop,$strand,$chromLengths->{$chrom},$parameters);
    my $bufferedSeq = getSequence($bufferedLocation,$genome);
    my($seq1Coords,$seq2Coords,$dupStructure,$duplexEnergy) = getDuplex($sequence,$bufferedSeq);
    return $duplexEnergy;
}

sub getGCcontent {
    my($seq,$length) = @_;
    my $gc = $seq =~ tr/[gcGC]//;
    my $gcContent = $gc / $length;
    return $gcContent;
}

sub getProductRelStartCounts {
    my($product) = @_;
    my($r10,$r9,$r8,$r7,$r6,$r5,$r4,$r3,$r2,$r1,$s0,$f1,$f2,$f3,$f4,$f5,$f6,$f7,$f8,$f9,$f10) = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
    my($id,$productReads,$adjMaxProductCount,$adjProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$product};
    foreach my $dRead (@{$productReads}) {
	my($relStart,$relStop,$dStart,$readCount,$hitCount,$libraryCounts,$seq,$adjustedSeqCount) = @{$dRead};
	$relStart = $relStart - $maxRelStart;  #adjusting to the products rel start
	my $adjReadCount = $readCount * $adjustedSeqCount;
	if ($relStart == -10) {
	    $r10 += $adjReadCount;
	} elsif ($relStart == -9) {
	    $r9 += $adjReadCount;
	} elsif ($relStart == -8) {
	    $r8 += $adjReadCount;
	} elsif ($relStart == -7) {
	    $r7 += $adjReadCount;
	} elsif ($relStart == -6) {
	    $r6 += $adjReadCount;
	} elsif ($relStart == -5) {
	    $r5 += $adjReadCount;
	} elsif ($relStart == -4) {
	    $r4 += $adjReadCount;
	} elsif ($relStart == -3) {
	    $r3 += $adjReadCount;
	} elsif ($relStart == -2) {
	    $r2 += $adjReadCount;
	} elsif ($relStart == -1) {
	    $r1 += $adjReadCount;
	} elsif ($relStart == 0) {
	    $s0 += $adjReadCount;
	} elsif ($relStart == 1) {
	    $f1 += $adjReadCount;
	} elsif ($relStart == 2) {
	    $f2 += $adjReadCount;
	} elsif ($relStart == 3) {
	    $f3 += $adjReadCount;
	} elsif ($relStart == 4) {
	    $f4 += $adjReadCount;
	} elsif ($relStart == 5) {
	    $f5 += $adjReadCount;
	} elsif ($relStart == 6) {
	    $f6 += $adjReadCount;
	} elsif ($relStart == 7) {
	    $f7 += $adjReadCount;
	} elsif ($relStart == 8) {
	    $f8 += $adjReadCount;
	} elsif ($relStart == 9) {
	    $f9 += $adjReadCount;
	} elsif ($relStart == 10) {
	    $f10 += $adjReadCount;
	}
    }
    return ($r10,$r9,$r8,$r7,$r6,$r5,$r4,$r3,$r2,$r1,$s0,$f1,$f2,$f3,$f4,$f5,$f6,$f7,$f8,$f9,$f10);
}

sub getDinucleotideFrequencies {
    my($seq) = @_;
    my $length = length($seq);
    my $aa = $seq =~ s/(aa|aA|Aa|AA)/$1/g;
    $aa = $aa/$length;
    my $ac = $seq =~ s/(ac|aC|Ac|AC)/$1/g;
    $ac = $ac/$length;
    my $ag = $seq =~ s/(ag|aG|Ag|AG)/$1/g;
    $ag = $ag/$length;
    my $at = $seq =~ s/(at|aT|At|AT)/$1/g;
    $at = $at/$length;
    my $ca = $seq =~ s/(ca|cA|Ca|CA)/$1/g;
    $ca = $ca/$length;
    my $cc = $seq =~ s/(cc|cC|Cc|CC)/$1/g;
    $cc = $cc/$length;
    my $cg = $seq =~ s/(cg|cG|Cg|CG)/$1/g;
    $cg = $cg/$length;
    my $ct = $seq =~ s/(ct|cT|Ct|CT)/$1/g;
    $ct = $ct/$length;
    my $ga = $seq =~ s/(ga|gA|Ga|GA)/$1/g;
    $ga = $ga/$length;
    my $gc = $seq =~ s/(gc|gC|Gc|GC)/$1/g;
    $gc = $gc/$length;
    my $gg = $seq =~ s/(gg|gG|Gg|GG)/$1/g;
    $gg = $gg/$length;
    my $gt = $seq =~ s/(gt|gT|Gt|GT)/$1/g;
    $gt = $gt/$length;
    my $ta = $seq =~ s/(ta|tA|Ta|TA)/$1/g;
    $ta = $ta/$length;
    my $tc = $seq =~ s/(tc|tC|Tc|TC)/$1/g;
    $tc = $tc/$length;
    my $tg = $seq =~ s/(tg|tG|Tg|TG)/$1/g;
    $tg = $tg/$length;
    my $tt = $seq =~ s/(tt|tT|Tt|TT)/$1/g;
    $tt = $tt/$length;
    return ($aa,$ac,$ag,$at,$ca,$cc,$cg,$ct,$ga,$gc,$gg,$gt,$ta,$tc,$tg,$tt);
}

sub getProductVectorInfo {
    my($product,$seq,$class,$rrLocation,$genome,$chromLengths,$prodMLParameters,$parameters) = @_;
    #calculating scores for products
    my($id,$productReads,$adjMaxProductCount,$adjProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$product};
    my $fivePrimeHet = computeFivePrimeHet($productReads) if ($prodMLParameters->{fivePrimeHet});
    my $length = length($seq) if ($prodMLParameters->{length});
    my $medianLength = getProductMedianLength($productReads) if ($prodMLParameters->{medianLength});
    my $gcContent = getGCcontent($seq,$length) if ($prodMLParameters->{gcContent});
    my $WFC = getWFC($seq) if ($prodMLParameters->{WFC});
    my($aa,$ac,$ag,$at,$ca,$cc,$cg,$ct,$ga,$gc,$gg,$gt,$ta,$tc,$tg,$tt) = getDinucleotideFrequencies($seq) unless ($prodMLParameters->{noDNF});
    my($r10,$r9,$r8,$r7,$r6,$r5,$r4,$r3,$r2,$r1,$s0,
       $f1,$f2,$f3,$f4,$f5,$f6,$f7,$f8,$f9,$f10) = getProductRelStartCounts($product) unless ($prodMLParameters->{noPositionalData});
    my $duplexEnergy = getProductDuplexEnergy($genome,$seq,$rrLocation,$chromLengths,$parameters) if ($prodMLParameters->{duplexEnergy});
    my @productVector = ($class);
    #pushing product scores into product vector array
    push(@productVector,$adjMaxProductCount) if ($prodMLParameters->{adjMaxProductCount});
    push(@productVector,$adjProductCount) if ($prodMLParameters->{adjProductCount});
    push(@productVector,$fivePrimeHet) if ($prodMLParameters->{fivePrimeHet});
    push(@productVector,$length) if ($prodMLParameters->{length});
    push(@productVector,$medianLength) if ($prodMLParameters->{medianLength});
    push(@productVector,$gcContent) if ($prodMLParameters->{gcContent});
    push(@productVector,$aa) if ($prodMLParameters->{aa});
    push(@productVector,$ac) if ($prodMLParameters->{ac});
    push(@productVector,$ag) if ($prodMLParameters->{ag});
    push(@productVector,$at) if ($prodMLParameters->{at});
    push(@productVector,$ca) if ($prodMLParameters->{ca});
    push(@productVector,$cc) if ($prodMLParameters->{cc});
    push(@productVector,$cg) if ($prodMLParameters->{cg});
    push(@productVector,$ct) if ($prodMLParameters->{ct});
    push(@productVector,$ga) if ($prodMLParameters->{ga});
    push(@productVector,$gc) if ($prodMLParameters->{gc});
    push(@productVector,$gg) if ($prodMLParameters->{gg});
    push(@productVector,$gt) if ($prodMLParameters->{gt});
    push(@productVector,$ta) if ($prodMLParameters->{ta});
    push(@productVector,$tc) if ($prodMLParameters->{tc});
    push(@productVector,$tg) if ($prodMLParameters->{tg});
    push(@productVector,$tt) if ($prodMLParameters->{tt});
    push(@productVector,$r10) if ($prodMLParameters->{r10});
    push(@productVector,$r9) if ($prodMLParameters->{r9});
    push(@productVector,$r8) if ($prodMLParameters->{r8});
    push(@productVector,$r7) if ($prodMLParameters->{r7});
    push(@productVector,$r6) if ($prodMLParameters->{r6});
    push(@productVector,$r5) if ($prodMLParameters->{r5});
    push(@productVector,$r4) if ($prodMLParameters->{r4});
    push(@productVector,$r3) if ($prodMLParameters->{r3});
    push(@productVector,$r2) if ($prodMLParameters->{r2});
    push(@productVector,$r1) if ($prodMLParameters->{r1});
    push(@productVector,$s0) if ($prodMLParameters->{s0});
    push(@productVector,$f1) if ($prodMLParameters->{f1});
    push(@productVector,$f2) if ($prodMLParameters->{f2});
    push(@productVector,$f3) if ($prodMLParameters->{f3});
    push(@productVector,$f4) if ($prodMLParameters->{f4});
    push(@productVector,$f5) if ($prodMLParameters->{f5});
    push(@productVector,$f6) if ($prodMLParameters->{f6});
    push(@productVector,$f7) if ($prodMLParameters->{f7});
    push(@productVector,$f8) if ($prodMLParameters->{f8});
    push(@productVector,$f9) if ($prodMLParameters->{f9});
    push(@productVector,$f10) if ($prodMLParameters->{f10});
    push(@productVector,$WFC) if ($prodMLParameters->{WFC});
    push(@productVector,$duplexEnergy) if ($prodMLParameters->{duplexEnergy});
    return \@productVector;
}

sub printProdTrainingData {
    my($bamList,$genomeDir,$chromLengths,$gff,$class,$prodTrainFile,$productFileHeader,$prodMLParameters,$parameters) = @_;
    my($hairpins,$mirBaseProducts) = readMirbaseGff3($gff);
    my $minMirOverlap = 0.60;  #require at least a 60% overlap of mirbase annotation for positive cases
    my $extraScanDistance = 15;  #distance added to the scanRegion so that retrieveReadData won't ignore reads that extend beyond region
    open(PVF,">>".$prodTrainFile) or die "failed to open $prodTrainFile for appending\n";
    my($posClass,$negClass) = getClassValues($parameters);
    foreach my $hairpinChrom (keys %{$hairpins}) {
	my $genome = loadGenome("$genomeDir/$hairpinChrom.fa");	
	foreach my $hairpinInfo (@{$hairpins->{$hairpinChrom}}) {
	    my($hairpinStart,$hairpinStop,$strand,$hairpinId,$hairpinName) = @{$hairpinInfo};
	    foreach my $product (@{$mirBaseProducts->{$hairpinId}}) {
		my($chrom,$start,$stop,$strand,$id,$name) = @{$product};
		my $location = "$chrom:$start..$stop:$strand";
		my $scanRegion = "$chrom:". ($start - $extraScanDistance) . ".." . ($stop + $extraScanDistance) . ":$strand";
		my($distinctReads,$uniqueReadCount,$totalDistinctReads,$adjustedReadCount) = retrieveReadData($scanRegion,$bamList,$parameters);
		if ($distinctReads->{$strand}) {  
		    my($products) = buildProducts($location,$distinctReads->{$strand},$strand,$parameters);
		    my ($mirProduct,$percentOverlap) = getGreatestOverlappingProduct($products,$location,$name);
		    my $mirProductLocation = getGenomicLocation($mirProduct,$hairpinChrom);
		    my $sequence = getSequence($mirProductLocation,$genome);
		    if ($percentOverlap > $minMirOverlap) {
			my $productName = $name . "_" . $mirProductLocation;
			my $productVector = getProductVectorInfo($mirProduct,$sequence,$class,$location,$genome,$chromLengths,$prodMLParameters,$parameters);
			my $productVectorLine = join("\t", @{$productVector});
			print PVF "$productName\t$productVectorLine\n";
		    }
		}
	    }
	}
    }
    close(PVF);
}

sub printProductVectorFile {
    my($readRegionsFile,$bamList,$genomeDir,$chromLengths,$classNum,$productFileHeader,$prodMLParameters,$parameters) = @_;
    my $productVectorFile = $parameters->{productVectorFile};
    open(PVF,">".$productVectorFile) or die "failed to open $productVectorFile for writing\n";
    open(RRF, $readRegionsFile) or die "failed to open $readRegionsFile\n";
    print PVF "$productFileHeader\n";
    my $prevChrom;
    my $genome;
    while (<RRF>) {
        chomp;
	unless(/\#/) {	  
	    my($chrom,$start,$stop,$rrId,$score,$strand) = split(/\t/,$_);
	    $start++; #converting to one based
	    #the chroms in RRF should be sorted at this point so the next line of code will only be executed once per chrom
	    if ($prevChrom) {
		$genome = loadGenome("$genomeDir/$chrom.fa") unless ($prevChrom eq $chrom);	
	    } else {
		$genome = loadGenome("$genomeDir/$chrom.fa");	
	    }
	    $prevChrom = $chrom;
	    my $location = "$chrom:$start..$stop:$strand";
	    my($distinctReads,$uniqueReadCount,$totalDistinctReads,$adjustedReadCount) = retrieveReadData($location,$bamList,$parameters);
	    if ($distinctReads->{$strand}) {  
		my($products) = buildProducts($location,$distinctReads->{$strand},$strand,$parameters);
		foreach my $product (@{$products}) {
		    my $productLocation = getGenomicLocation($product,$chrom);
		    my($productId) = @{$product};
			my $productName = $rrId . "-". $productId . "_" . $productLocation;
		    my $sequence = getSequence($productLocation,$genome);
		    my $productVector = getProductVectorInfo($product,$sequence,$classNum,$location,$genome,$chromLengths,$prodMLParameters,$parameters);
		    my $productVectorLine = join("\t", @{$productVector});
		    print PVF "$productName\t$productVectorLine\n";
		}
	    }
	}
    }
    close(PVF);
    close(RRF);
}

sub getProdFileHeader {
    #returns a header with the variables that were either defaults or specified by the user
    my($prodMLParameters) = @_;
    my $header = "GeneId\ty";
    $header .= "\tadjMaxProductCount" if $prodMLParameters->{adjMaxProductCount};
    $header .= "\tadjProductCount" if $prodMLParameters->{adjProductCount};
    $header .= "\tfivePrimeHet" if $prodMLParameters->{fivePrimeHet};
    $header .= "\tlength" if $prodMLParameters->{length};
    $header .= "\tmedianLength" if $prodMLParameters->{medianLength};
    $header .= "\tgcContent" if $prodMLParameters->{gcContent};
    $header .= "\taa" if $prodMLParameters->{aa};
    $header .= "\tac" if $prodMLParameters->{ac};
    $header .= "\tag" if $prodMLParameters->{ag};
    $header .= "\tat" if $prodMLParameters->{at};
    $header .= "\tca" if $prodMLParameters->{ca};
    $header .= "\tcc" if $prodMLParameters->{cc};
    $header .= "\tcg" if $prodMLParameters->{cg};
    $header .= "\tct" if $prodMLParameters->{ct};
    $header .= "\tga" if $prodMLParameters->{ga};
    $header .= "\tgc" if $prodMLParameters->{gc};
    $header .= "\tgg" if $prodMLParameters->{gg};
    $header .= "\tgt" if $prodMLParameters->{gt};
    $header .= "\tta" if $prodMLParameters->{ta};
    $header .= "\ttc" if $prodMLParameters->{tc};
    $header .= "\ttg" if $prodMLParameters->{tg};
    $header .= "\ttt" if $prodMLParameters->{tt};
    $header .= "\tr10" if $prodMLParameters->{r10};
    $header .= "\tr9" if $prodMLParameters->{r9};
    $header .= "\tr8" if $prodMLParameters->{r8};
    $header .= "\tr7" if $prodMLParameters->{r7};
    $header .= "\tr6" if $prodMLParameters->{r6};
    $header .= "\tr5" if $prodMLParameters->{r5};
    $header .= "\tr4" if $prodMLParameters->{r4};
    $header .= "\tr3" if $prodMLParameters->{r3};
    $header .= "\tr2" if $prodMLParameters->{r2};
    $header .= "\tr1" if $prodMLParameters->{r1};
    $header .= "\ts0" if $prodMLParameters->{s0};
    $header .= "\tf1" if $prodMLParameters->{f1};
    $header .= "\tf2" if $prodMLParameters->{f2};
    $header .= "\tf3" if $prodMLParameters->{f3};
    $header .= "\tf4" if $prodMLParameters->{f4};
    $header .= "\tf5" if $prodMLParameters->{f5};
    $header .= "\tf6" if $prodMLParameters->{f6};
    $header .= "\tf7" if $prodMLParameters->{f7};
    $header .= "\tf8" if $prodMLParameters->{f8};
    $header .= "\tf9" if $prodMLParameters->{f9};
    $header .= "\tf10" if $prodMLParameters->{f10};
    $header .= "\tWFC" if $prodMLParameters->{WFC};
    $header .= "\tduplexEnergy" if ($prodMLParameters->{duplexEnergy});
    return $header;
}

sub createProductBedFile {
    my($predProductClassFile,$parameters) = @_;
    my $productFile = $parameters->{predProductRF};
    open(PPCF,$predProductClassFile) or die "failed to open $predProductClassFile\n";
    open(PF,">$productFile") or die "failed to open $productFile for writing\n";
    while (<PPCF>) {
	chomp;
	my($product,$classNum,$RFavg) = split(/\t/, $_);
	my $count = 1000;  #temporary dummy val
	if ($classNum == 1) {
	    my($productName,$productLocation) = $product =~ /(.*?)\_(.*)/;
	    my($chrom,$start,$stop,$strand) = parseLocation($productLocation);
	    print PF "$chrom\t$start\t$stop\t$productName\t$count\t$strand\n";
	}
    }
    close(PPCF);
    close(PF);
}

#############################################
# MACHINE LEARNING Hairpin DATA SUBROUTINES #
#############################################

sub getFoldInfo {
    my($hairpinFile) = @_;
    my %foldInfo;
    open(HP,$hairpinFile) or die "failed to open $hairpinFile for reading\n";
    while(<HP>) {
	chomp;
	unless ( /^#/ ) {
	    my($name,$chrom,$start,$stop,$strand,$leftCenter,$rightCenter,$totalSense,$totalAntisense,$mfe,$sequence,$fold) = split("\t",$_);
	    my($unfoldedStart,$trueFold,$unfoldedEnd) = $fold =~ /(^\+*)(.*?)(\+*$)/;
	    my $foldStart = -1;  #dummy var
	    if ($strand eq "+") {
		$foldStart = $start + length($unfoldedStart);
	    } else {
		$foldStart = $start + length($unfoldedEnd);		
	    }
	    my $foldLength = length($trueFold);
	    my $foldStop = $foldStart + $foldLength - 1;
	    $foldInfo{$chrom}{$name} = [$start,$stop,$strand,$foldStart,$foldStop,$trueFold];
	}
    }
    close(HP);
    return \%foldInfo;
}

sub getOverlappingFolds {
    my($foldInfo) = @_;
    my %overlappingFolds;
    #collecting hairpins with the same folds toghether into a hash
    foreach my $chrom (keys %{$foldInfo}) {
	foreach my $hpName (keys %{$foldInfo->{$chrom}}) {
	    my($start,$stop,$strand,$foldStart,$foldStop,$fold) = @{$foldInfo->{$chrom}{$hpName}};
	    push(@{$overlappingFolds{"$chrom:$foldStart..$foldStop:$strand"}{$fold}}, [$hpName,$start,$stop]);
	}
    }
    return \%overlappingFolds;
}

sub dropPartialOverlappingFolds {
    my($overlappingFolds) = @_;
    my %chromPositions;
    my %newOverlappingFolds;
    foreach my $location (keys %{$overlappingFolds}) {
	my ($chrom,$start,$stop,$strand) = parseLocation($location);
	push(@{$chromPositions{$chrom}{$strand}}, [$start,$stop]);
    }
    foreach my $chrom (keys %chromPositions) {
	foreach my $strand (keys %{$chromPositions{$chrom}}) {
	    my($lastStart,$lastStop) = (-1,-1);
	    #the command below should sort positions so that the longest sequence at a particular start position is given before a shorter one
	    foreach my $position (sort {$a->[0] <=> $b->[0] || $b->[1] <=> $a->[1]} @{$chromPositions{$chrom}{$strand}}) {
		my($start,$stop) = @{$position};
		unless (getOverlap($start,$stop,$lastStart,$lastStop) == $stop - $start + 1) {
		    $newOverlappingFolds{"$chrom:$start..$stop:$strand"} = $overlappingFolds->{"$chrom:$start..$stop:$strand"};
		    ($lastStart,$lastStop) = ($start,$stop);
		}  #else - these positions are already part of a larger fold and not helpful to our anaylysis
#		else {
#		    #this part is for testing
#		    foreach my $fold (keys %{$overlappingFolds->{"$chrom:$start..$stop:$strand"}}) {
#		        foreach my $hpInfo (@{$overlappingFolds->{"$chrom:$start..$stop:$strand"}{$fold}}) {
#			    my($hpName,$hpStart,$hpStop) = @{$hpInfo};
#			    print "$hpName dropped.  partially overlaps another fold \($start,$stop\) overlaps \($lastStart,$lastStop\)\n";
#			}
#		    }
#		}
	    }
	}
    }
    return \%newOverlappingFolds;
}

sub dropDuplicateFolds {
    my($hairpinFile) = @_;
    my %nonDuplicateHairpins;
    my $foldInfo = getFoldInfo($hairpinFile);
    my $overlappingFolds = getOverlappingFolds($foldInfo);
    $overlappingFolds = dropPartialOverlappingFolds($overlappingFolds);
    #now determining the most centered fold to keep
    foreach my $location (keys %{$overlappingFolds}) {
	my($chrom,$foldStart,$foldStop,$strand) = parseLocation($location);
	foreach my $fold (keys %{$overlappingFolds->{$location}}) {
	    #print "\n$location\n$fold\n";
	    my $minCenterProximity = ~0;
	    my $mostCentered;
	    foreach my $info (@{$overlappingFolds->{$location}{$fold}}) {
		my($hpName,$start,$stop) = @{$info};
		#print "$hpName\n";
		my $centerProximity = abs(($foldStart - $start) - ($stop - $foldStop));
		if($centerProximity < $minCenterProximity) {
		    $minCenterProximity = $centerProximity;
		    $mostCentered = $hpName;
		    #print "most centered now $mostCentered\n";
		}
	    }
	    $nonDuplicateHairpins{$mostCentered} = 1;
	}
    }
    return \%nonDuplicateHairpins;
}

sub getHairpinFoldRegions {
    #this function was designed to take care of a bug in which the whole hairpin region folded and unfolded was being checked
    #to see if it was overlapping a true mir in getMirOverlaps.  This was done to avoid calling hairpins true mirs if their
    #unvolded region overlapped the annotation.  In the future, processReadRegions could add the start and stop of the fold region
    #to the hairpins file to avoid needing this function and to create cleaner code.
    my($hairpinsFile) = @_;
    my %foldRegions;
    open(HP,$hairpinsFile) or die "Failed to open Hairpins File for reading\n";
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

sub getHPReadCounts {
    my($hairpinsFile) = @_;
    my %hpCounts;
    open(HP, $hairpinsFile) or die "Failed to open hairpins file for reading\n";
    while (<HP>) {
	chomp;
	unless ( /^#/ ) {
	    my($name,$chrom,$start,$stop,$strand,$leftCenter,$rightCenter,$totalSense) = split("\t",$_);
	    $hpCounts{$name} = $totalSense;
	}
    }
    close(HP);
    return \%hpCounts;
}

sub extractGffReadRegions {
    my($gffPrefix,$hairpinsFile,$productsFile,$distinctReadsFile) = @_;
    my $gffFile = $gffPrefix . ".gff";
    my $newHairpinsFile = $gffPrefix . "_hairpins.txt";
    my $newProductFile = $gffPrefix . "_products.txt";
    my $newDistinctReadsFile = $gffPrefix . "_distinctReads.txt";
    my %gffEntry;
    my %RFscore;
    open(GFF,$gffFile) or die "failed to open $gffFile for reading\n";
    while (<GFF>) {
	chomp;
	unless ( /^#/ ) {
	    my($chrom,$source,$type,$start,$stop,$score,$strand,$phase,$info) = split(/\t/);
	    my %info;
	    my @terms = split(/;/,$info);
	    foreach my $term (@terms) {
		my($key,$value) = $term =~ /(.*)=(.*)/;
		$info{$key} = $value;
	    }
	    my $name = $info{Name} or die "No ID found for the line:\n$_\n"; 
	    $gffEntry{$name} = 1;
	    $RFscore{$name} = $score;
	}
    }
    close(GFF);
    open(HP,$hairpinsFile) or die "failed to open $hairpinsFile for reading\n";
    open(NHP,">$newHairpinsFile") or die "failed to open $newHairpinsFile for writing\n";
    while (<HP>) {
	chomp;
	my $line = $_;
	unless ( $line =~ /^#/ ) {
	    my($name) = split("\t",$line);
	    if ($gffEntry{$name}) {
		my $score = $RFscore{$name};
		print NHP "$line\t$score\n";
	    }
	} else {
	    print NHP "$line\tHPscore\n";
	}
    }
    close(NHP);
    close(HP);
    open(PF,$productsFile) or die "failed to open $productsFile for reading\n";
    open(NPF,">$newProductFile") or die "failed to open $newProductFile for writing\n";
    while (<PF>) {
	chomp;
	my $line = $_;
	unless ( $line =~ /^#/ ) {
	    my($name) = split("\t",$line);
	    if ($gffEntry{$name}) {
		print NPF "$line\n";
	    }
	} else {
	    print NPF "$line\n";
	}
    }
    close(NPF);
    close(PF);
    open(DR,$distinctReadsFile) or die "failed to open $distinctReadsFile for reading\n";
    open(NDR,">$newDistinctReadsFile") or die "failed to open $newDistinctReadsFile for writing\n";
    while (<DR>) {
	chomp;
	my $line = $_;
	unless ( $line =~ /^#/ ) {
	    my($name) = split("\t",$line);
	    if ($gffEntry{$name}) {
		print NDR "$line\n";
	    }
	} else {
	    print NDR "$line\n";
	}
    }
    close(NDR);
    close(DR);    
}

#sub overlapsPositivePrediction {
#    my($hpScoreInfo,$usedScores) = @_;
#    foreach my $scoreInfo (@{$usedScores}) {
#	my($start1,$stop1) = @{$hpScoreInfo};
#	my($start2,$stop2) = @{$scoreInfo};
#	if (getOverlap($start1,$stop1,$start2,$stop2)) {
#	    return 1;
#	}
#    }
#    return 0;
#}

#sub dropPosPredictedOverlaps {
#    my($positivePredictions,$RFPredictions) = @_;
#    my %newPositivePredictions;
#    my %posOverlaps;
#    foreach my $chrom (keys %{$positivePredictions}) {
#	my %hpScores;
#	foreach my $predictionInfo (@{$positivePredictions->{$chrom}}) {
#	    my($start,$stop,$strand,$geneId,$name) = @{$predictionInfo};
#	    push(@{$hpScores{$strand}}, [$start,$stop,$RFPredictions->{$name},$geneId,$name]);
#	}
#	foreach my $strand (keys %hpScores) {
#	    my @usedScores;
#	    my @sortedHPScores = sort {$b->[2] <=> $a->[2]} @{$hpScores{$strand}};
#	    foreach my $hpScoreInfo (@sortedHPScores) {
#		my($start,$stop,$hpScore,$geneId,$name) = @{$hpScoreInfo};
#		if (overlapsPositivePrediction($hpScoreInfo,\@usedScores)) {
#		    push(@{$posOverlaps{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
#		} else {
#		    push(@{$newPositivePredictions{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
#		    push(@usedScores, [$start,$stop,$hpScore,$geneId,$name]);
#		}
#	    }
#	}
#    }
#    return (\%newPositivePredictions, \%posOverlaps);
#}

sub overlapsPrediction {
    my($hpScoreInfo,$usedScores,$overlapCutoff,$hairpinFoldRegions) = @_;
    my $greatestOverlapAmt = 0;
    my $greatestOverlappingHP = 0;
    my($start1,$stop1,$hpScoreInfo1,$geneId1,$name1) = @{$hpScoreInfo};
    my($foldStart1,$foldStop1) = @{$hairpinFoldRegions->{$name1}};
    foreach my $scoreInfo (@{$usedScores}) {
	my($start2,$stop2,$hpScoreInfo2,$geneId2,$name2) = @{$scoreInfo};
	my($foldStart2,$foldStop2) = @{$hairpinFoldRegions->{$name2}};
	my $overlapAmt = getOverlap($foldStart1,$foldStop1,$foldStart2,$foldStop2);
	if ($overlapAmt > $overlapCutoff && $overlapAmt > $greatestOverlapAmt) {
	    $greatestOverlapAmt = $overlapAmt;
	    $greatestOverlappingHP = $name2;
	}
    }
    return $greatestOverlappingHP;
}

sub classifyPredictions {
    my($predictions,$RFPredictions,$RFCutoff,$dropOverlaps,$hairpinFoldRegions) = @_;
    my $overlapCutoff = 1;  #require at least this many nucleotides to overlap before dropping
                             #this may be changed later to drop based on overlapping major products
    my %positivePredictions;
    my %negativePredictions;
    my %overlaps;
    foreach my $chrom (keys %{$predictions}) {
	my %hpScores;
	foreach my $predictionInfo (@{$predictions->{$chrom}}) {
	    my($start,$stop,$strand,$geneId,$name) = @{$predictionInfo};
	    push(@{$hpScores{$strand}}, [$start,$stop,$RFPredictions->{$name},$geneId,$name]);
	}
	foreach my $strand (keys %hpScores) {
	    my @usedScores;
	    my @sortedHPScores = sort {$b->[2] <=> $a->[2]} @{$hpScores{$strand}};
	    foreach my $hpScoreInfo (@sortedHPScores) {
		my($start,$stop,$hpScore,$geneId,$name) = @{$hpScoreInfo};
		my $predictedOverlap = overlapsPrediction($hpScoreInfo,\@usedScores,$overlapCutoff,$hairpinFoldRegions);
		if ($dropOverlaps && $predictedOverlap) {
		    push(@{$overlaps{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
		    print "dropping $name because it overlaps $predictedOverlap\n";
		} elsif ($hpScore >= $RFCutoff) {
		    push(@{$positivePredictions{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
		} else {
		    push(@{$negativePredictions{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
		}
		push(@usedScores, [$start,$stop,$hpScore,$geneId,$name]);
	    }
	}
    }
    return (\%positivePredictions,\%negativePredictions,\%overlaps);
}

sub getAnnotatedAntisense {
    my($annotatedMirs) = @_;
}

sub addAntiSenseToAnnotations {
    my($annotatedMirs) = @_;
    my $knowAntisense = getAnnotatedAntisense($annotatedMirs);
    foreach my $chrom (keys %{$annotatedMirs}) {
	my %antiSenseAnnotationExists;
	foreach my $hairpinInfo1 (@{$annotatedMirs->{$chrom}}) {
	    foreach my $hairpinInfo2 (@{$annotatedMirs->{$chrom}}) {
		my($start1,$stop1,$strand1,$id1,$name1) = @{$hairpinInfo1};
		my($start2,$stop2,$strand2,$id2,$name2) = @{$hairpinInfo2};
		if ($strand1 ne $strand2 && getOverlap($start1,$stop1,$start2,$stop2)) {
		    $antiSenseAnnotationExists{$name1} = 1;
		    $antiSenseAnnotationExists{$name2} = 1;
		}
	    }
	}
	foreach my $hairpinInfo (@{$annotatedMirs->{$chrom}}) {
	    my($start,$stop,$strand,$id,$name) = @{$hairpinInfo};
	    unless ($antiSenseAnnotationExists{$name}) {
		my $asName = "as-$name";
		my $asId = "as-$id";
		my $newStrand = ($strand eq '+') ? '-' : '+';
		push(@{$annotatedMirs->{$chrom}}, [$start,$stop,$newStrand,$asId,$asName]);
	    }
	}
    }
    return $annotatedMirs;
}

sub printPredictedMirGff {
    my($predHairpinClassFile,$gff,$productFile,$hairpinsFile,$parameters) = @_; 
    print "in predicted mir gff\n";
    print "predHairpinClassFile = $predHairpinClassFile\n";
    print "gff = $gff\n";
    print "productFile = $productFile\n";
    print "hairpinsFile = $hairpinsFile\n";
    my($posGff,$negGff,$overlapGff) = ("positivePredictions.gff","negativePredictions.gff","droppedOverlaps.gff");
    my $useLargestProduct = 1;
    my $dropOverlaps = 1;
    my $minOverlapPercent = 0.7;
    my $RFCutoff = 0.5;
    my %predictions;
    my %RFPredictions;
    my $hpCounts = getHPReadCounts($hairpinsFile);
    open(PHC,$predHairpinClassFile) or die "failed to open $predHairpinClassFile\n";
    while (<PHC>) {
	chomp;
	my($geneId,$predClass,$RFavg) = split("\t", $_);
	my($name,$location) = parseGeneId($geneId);
	my($chrom,$start,$stop,$strand) = parseLocation($location);
	push(@{$predictions{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
	$RFPredictions{$name} = $RFavg;
    }
    close(PHC);
    #updating with annotated names if one exists
    my $annotatedMirNames;
    print "classes retrieved\n";
    my $hairpinFoldRegions = getHairpinFoldRegions($parameters->{hairpinsFile});
    if ($gff) {
	my($annotatedMirs) = readMirbaseGff3($gff);
	print "adding antiSense to annotations\n";
	
	#too slow.  I'll have to change this
#	$annotatedMirs = addAntiSenseToAnnotations($annotatedMirs);


	#getHairpinFoldRegions is a temorary fix needed to product the correct mirOverlaps.  See function for more info.
	print "getting hairpin fold regions\n";
	$annotatedMirNames = renameMirOverlaps(\%predictions,$annotatedMirs,$productFile,$useLargestProduct,$minOverlapPercent,$hairpinFoldRegions);
    }
    #printing out gffs
    my($posPredictions,$negPredictions,$overlaps) = classifyPredictions(\%predictions,\%RFPredictions,$RFCutoff,$dropOverlaps,$hairpinFoldRegions);
    print "printing gff\'s\n";
    open(PGFF,">$posGff") or die "failed to open $posGff for writing\n";
    foreach my $chrom (%{$posPredictions}) {
	foreach my $hpInfo (@{$posPredictions->{$chrom}}) {
	    my($start,$stop,$strand,$geneId,$name) = @{$hpInfo};
	    my $RFavg = $RFPredictions{$name};
	    my $readCount = $hpCounts->{$name};
	    my $hpInfo = "ID=$name;Alias=$name;Name=$name";
	    if($annotatedMirNames->{$name}) {
		my $mirName = $annotatedMirNames->{$name};
		$hpInfo = "ID=$name;Alias=$mirName;Name=$name";
	    }
	    print PGFF "$chrom\tmiRWOODS\tmiRNA_primary_transcript\t$start\t$stop\t$RFavg\t$strand\t.\t$hpInfo\n";
	}
    }
    close(PGFF);
    open(NGFF,">$negGff") or die "failed to open $negGff for writing\n";
    foreach my $chrom (%{$negPredictions}) {
	foreach my $hpInfo (@{$negPredictions->{$chrom}}) {
	    my($start,$stop,$strand,$geneId,$name) = @{$hpInfo};
	    my $RFavg = $RFPredictions{$name};
	    my $readCount = $hpCounts->{$name};
	    my $hpInfo = "ID=$name;Alias=$name;Name=$name";
	    if($annotatedMirNames->{$name}) {
		my $mirName = $annotatedMirNames->{$name};
		$hpInfo = "ID=$name;Alias=$mirName;Name=$name";
	    }
	    print NGFF "$chrom\tmiRWOODS\tmiRNA_primary_transcript\t$start\t$stop\t$RFavg\t$strand\t.\t$hpInfo\n";
	}
    }
    close(NGFF);
    #the following is for testing purposes and will eventually be included in the negGFF predicitons
    if ($dropOverlaps) {
	open(OGFF,">$overlapGff") or die "failed to open overlappingPosPredictions.gff for writing\n";
	foreach my $chrom (%{$overlaps}) {
	    foreach my $hpInfo (@{$overlaps->{$chrom}}) {
		my($start,$stop,$strand,$geneId,$name) = @{$hpInfo};
		my $RFavg = $RFPredictions{$name};
		my $readCount = $hpCounts->{$name};
		my $hpInfo = "ID=$name;Alias=$name;Name=$name";
		if($annotatedMirNames->{$name}) {
		    my $mirName = $annotatedMirNames->{$name};
		    $hpInfo = "ID=$name;Alias=$mirName;Name=$name";
		}
		print OGFF "$chrom\tmiRWOODS\tmiRNA_primary_transcript\t$start\t$stop\t$RFavg\t$strand\t.\t$hpInfo\n";
	    }
	}
    }
    close(OGFF);
}

sub overlapsMirHairpin {
    my($location,$mirHairpins,$minOverlapPercent) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $minOverlap = $minOverlapPercent * ($stop - $start + 1);
    foreach my $mirHairpin (@{$mirHairpins->{$chrom}}) {
	my($mirStart,$mirStop,$strand) = @{$mirHairpin};
	my $overlap = getOverlap($mirStart,$mirStop,$start,$stop);
	if ($overlap > $minOverlap) {
	    return 1;
	}
    }
    return 0;
}

sub renameMirOverlaps {
    my($hairpins,$annotatedMirs,$productFile,$useLargestProduct,$minOverlapPercent,$hairpinFoldRegions) = @_;
    my %mirNames;
    #the follwoing code retrives hairpins that fit the criteria to be labeled as the real mir
    my($overlaps, $overlapInfo) = getMirOverlaps($annotatedMirs, $hairpins, $minOverlapPercent,$hairpinFoldRegions);
    my $overlapProducts = getOverlapProducts($productFile, $overlapInfo);
    if ($useLargestProduct) {
	$overlaps = improveMirOverlaps($annotatedMirs, $overlaps, $overlapProducts);
    }
    foreach my $mirName (keys %{$overlaps}) {
	my $i = 0;
	foreach my $overlapName (@{$overlaps->{$mirName}}) {
	    $mirNames{$overlapName} = ($i) ? $mirName . "_MW" . chr(ord('a') + $i) : $mirName;
	    $i++;
	}
    }
    #the following code retrieves hairpins that overlap a mir but do not fit the criteria
    my($allOverlaps, $allOverlapInfo) = getMirOverlaps($annotatedMirs, $hairpins, 0,$hairpinFoldRegions);
    foreach my $mirName (keys %{$allOverlaps}) {
	my $i = 0;
	foreach my $overlapName (@{$allOverlaps->{$mirName}}) {
	    unless ($mirNames{$overlapName}) {
		$mirNames{$overlapName} = ($i) ? $mirName . "overlap" . chr(ord('a') + $i) : $mirName . "overlap";
		$i++;
	    }
	}
    }
    return \%mirNames;
}

sub oldRenameMirOverlaps {
    my($mirOverlaps) = @_;
    my %mirNames;
    foreach my $mirName (keys %{$mirOverlaps}) {
	my $i = 0;
	foreach my $overlapName (@{$mirOverlaps->{$mirName}}) {
	    $mirNames{$overlapName} = ($i) ? $mirName + "_MW" + chr("$i") : $mirName;
	    $i++;
	}
    }
    return \%mirNames;
}


sub getMirOverlaps {
    my($mirHairpins,$hairpins,$minOverlapPercent,$hairpinFoldRegions) = @_;
    my %mirOverlaps;
    my %overlapInfo;
    #print "minOverlapPercent = $minOverlapPercent\n";
    foreach my $chrom (keys %{$hairpins}) {
	foreach my $mirInfo (@{$mirHairpins->{$chrom}}) {
	    my($mirStart,$mirStop,$mirStrand,$mirId,$mirName) = @{$mirInfo};
	    foreach my $hpInfo (@{$hairpins->{$chrom}}) {
		my($start,$stop,$strand,$geneId,$name) = @{$hpInfo};
		if ($strand eq $mirStrand) {
		    my($foldStart,$foldStop) = @{$hairpinFoldRegions->{$name}};
		    my $minOverlap = ($foldStop - $foldStart + 1) * $minOverlapPercent;
		    my $overlap = getOverlap($foldStart,$foldStop,$mirStart,$mirStop);
		    if ($overlap > $minOverlap) {
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

sub getProductsFromFile {
    my($hairpins,$productFile) = @_;
    my %products;
    open(PF,"$productFile") or die "failed to open $productFile\n";
    while(<PF>) {
	chomp;
	my($name,$sidetype,$type,$total,$totalMostAbundant,$relStart,$relStop,$strand) = split(/\t/,$_);
	foreach my $chrom (keys %{$hairpins}) {
	    foreach my $hairpinName (keys %{$hairpins->{$chrom}}) {
		my($hairpinInfo) = @{$hairpins->{$chrom}{$hairpinName}};
		my($hpStart,$hpStop,$hpStrand) = @{$hairpinInfo};
		my($prodStart,$prodStop) = (-2,-1);  #initialized with dummy variables for testing
		if ($strand eq '+') {
		    $prodStart = $hpStart + $relStart;
		    $prodStop = $hpStart + $relStop;
		} elsif ($strand eq '-') {
		    $prodStop = $hpStop - $relStart;
		    $prodStart = $hpStop - $relStop;
		}
		if ($prodStart > $prodStop) {
		    die "prodStart is greater than prodStop\n";
		} 
		if ($hairpinName eq $name) {
		    push(@{$products{$chrom}{$hairpinName}}, [$sidetype,$type,$total,$totalMostAbundant,$prodStart,$prodStop,$strand]);
		}
	    }
	}
    }
    close(PF);
    return (\%products);
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
	    my($name,$location) = parseGeneId($geneId);
	    my($chrom,$start,$stop,$strand) = parseLocation($location);
	    push(@{$hairpinLocations{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
	} else {
	    $first = 0;
	}
    }
    close(GIDF);
    return \%hairpinLocations
}

sub getOverlapProducts {
    my($productFile, $overlapInfo) = @_;
    my %overlapProducts;
    open(PF,"$productFile") or die "failed to open $productFile\n";
    while(<PF>) {
	chomp;
	my($name,$sidetype,$type,$total,$totalMostAbundant,$relStart,$relStop,$strand) = split(/\t/,$_);
	if ($overlapInfo->{$name}) {
	    my($hpStart,$hpStop) = @{$overlapInfo->{$name}};
	    my($prodStart,$prodStop) = (-1,-1);
	    if ($strand eq '+') {
		$prodStart = $hpStart + $relStart;
		$prodStop = $hpStart + $relStop;
	    } elsif ($strand eq '-') {
		$prodStop = $hpStop - $relStart;
		$prodStart = $hpStop - $relStop;
	    }
	    push(@{$overlapProducts{$name}}, [$sidetype,$type,$total,$totalMostAbundant,$prodStart,$prodStop,$strand]);
	}
    }
    close(PF);
    return \%overlapProducts;
}

sub improveMirOverlaps {
    my($mirHairpins,$mirOverlaps,$overlapProducts) = @_;
    my $mirsOnly = 1;  #use mirs as largest products (unless there are no mirs)
    my(%improvedOverlaps);
    #gathering all products that overlap a mir
    foreach my $chrom (keys %{$mirHairpins}) {
	foreach my $mirInfo (@{$mirHairpins->{$chrom}}) {
	    my($mirStart,$mirStop,$mirStrand,$mirId,$mirName) = @{$mirInfo};
	    my($largestProd,$largestMir) = (0,0);
	    foreach my $hpName (@{$mirOverlaps->{$mirName}}) {
		foreach my $overlapProducts ($overlapProducts->{$hpName}) {
		    foreach my $productInfo (@{$overlapProducts}) {
			my($sidetype,$type,$total,$totalMostAbundant,$prodStart,$prodStop,$strand) = @{$productInfo};
			if (($strand eq $mirStrand) && getOverlap($mirStart,$mirStop,$prodStart,$prodStop)) {
			    if ($mirsOnly && ($type =~ /miR/)) {
				#a mir is present.  Checking to see if the total is greater than the largest mir
				if ($total > $largestMir) {
				    @{$improvedOverlaps{$mirName}} = ($hpName); #resetting array
				} elsif ($total == $largestMir) {
				    push(@{$improvedOverlaps{$mirName}}, $hpName)
				}
			    } elsif ($largestMir == 0) {
				#a mir is not present yet.  Checking to see if the total is greater than the largest product
				if ($total > $largestProd) {
				    @{$improvedOverlaps{$mirName}} = ($hpName); #resetting array
				} elsif ($total == $largestProd) {
				    push(@{$improvedOverlaps{$mirName}}, $hpName)
				}
			    }
			}
		    }
		}
	    }
	}
    }
    return \%improvedOverlaps;
}

sub renameProducts {
    my($productFile,$annotatedMirNames) = @_;
    my $newProductFile = "newProducts.txt";
    open(NPF, ">$newProductFile") or die "failed to open $newProductFile for writing\n";
    open(PF, $productFile) or die "failed to open product file\n";
    while (<PF>) {
	chomp;
	unless ( /^#/ ) {
	    my($name,@other) = split("\t",$_);
	    my $newName = ($annotatedMirNames->{$name}) ? $annotatedMirNames->{$name} : $name;
	    print NPF $newName . "\t" . join("\t", @other) . "\n";
	} else {
	    print NPF $_ . "\n";
	}
    }
    close(PF);
    close(NPF);
}

sub printHairpinTrainingData {
    my($scoresFile,$gff,$genomeDir,$chromLengths,$class,$trainFile,$productFile,$parameters) = @_;
    my $useLargestProduct = 1;
    my $minOverlapPercent = 0.7;
    my($mirHairpins,$mirProducts) = readMirbaseGff3($gff);
    my $hairpins = getHairpinsFromGeneIds($scoresFile,1);
    #getHairpinFoldRegions is a temorary fix needed to product the correct mirOverlaps.  See function for more info.
    my $hairpinFoldRegions = getHairpinFoldRegions($parameters->{hairpinsFile});
    my($mirOverlaps, $overlapInfo) = getMirOverlaps($mirHairpins,$hairpins,$minOverlapPercent,$hairpinFoldRegions);
    my $overlapProducts = getOverlapProducts($productFile,$overlapInfo);
    if ($useLargestProduct) {
	#we wish to identify hairpins by a positive class only if their largest product overlaps with annotated mirs
	$mirOverlaps = improveMirOverlaps($mirHairpins,$mirOverlaps,$overlapProducts);
    }
    my %positiveHairpinEntries;
    foreach my $mirName (keys %{$mirOverlaps}) {
	foreach my $hairpinName (@{$mirOverlaps->{$mirName}}) {
	    my($start,$stop,$strand,$geneId) = @{$overlapInfo->{$hairpinName}};
	    $positiveHairpinEntries{$geneId} = 1;
	}
    }
    open(HVF,">".$trainFile) or die "failed to open $trainFile for writing\n";
    open(HPS, $scoresFile) or die "failed to open $scoresFile for reading\n";
    my $first = 1;
    while (<HPS>) {
	chomp;
	my $line = $_;
	unless($first) {
	    my($geneId,$class,@other) = split(/\t/,$line);
	    if ($positiveHairpinEntries{$geneId}) {
		print HVF "$geneId\t1\t" . join("\t", @other)  . "\n";
	    }
	} else {
	    $first = 0;
	    print HVF $line . "\n";
	}	
    }
    close(HPS);
    close(HVF);
}

sub getHPProdCounts {
    my($hairpinVecTrainFile,$productFile) = @_;
    my %prodCounts;
    my %geneIds;
    open(HVF,$hairpinVecTrainFile) or die "failed to open $hairpinVecTrainFile for reading\n";
    while (<HVF>) {
	chomp;
	my($geneId) = split("\t", $_);
	my($hpName,$hpLocation) = parseGeneId($geneId);
	$geneIds{$hpName} = $geneId;
    }
    close(HVF);
    open(PF,$productFile) or die "failed to open $productFile for reading\n";
    while (<PF>) {
	chomp;
	unless ( /^#/ ) {
	    my($tag,$side,$type,$total,$totalMostAbundant,$start,$stop,$strand) = split("\t");
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
    my($mirHairpins,$mirProducts) = readMirbaseGff3($gff);
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
    my($mirHairpins,$mirProducts) = readMirbaseGff3($gff);
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

sub addRandomNegHPFeatures {
    #adds random negative features from feature vector file to the training data
    my($hairpinVecTrainFile,$hairpinScoresFile) = @_;
    my %posNames;
    my $posHairpinCount = 0;
    my @negVectorData;
    open(HVF, $hairpinVecTrainFile) or die "failed to open $hairpinVecTrainFile\n";
    while(<HVF>) {
	chomp;
	my($hpName) = split("\t", $_);
	$posNames{$hpName} = 1;
	$posHairpinCount++;
    }
    close(HVF);
    open(HSF, $hairpinScoresFile) or die "failed to open $hairpinScoresFile\n";
    while(<HSF>) {
	chomp;
	my $vectorLine = $_;
	my($hpName) = split("\t", $vectorLine);
	unless ($posNames{$hpName}) {
	    push(@negVectorData, $vectorLine);
	}
    }
    close(HSF);
    my @shuffledNegVectorData = List::Util::shuffle(@negVectorData);
    open(HVFA, ">>$hairpinVecTrainFile") or die "failed to open $hairpinVecTrainFile for appending\n";
    my $negHairpinCount = ($posHairpinCount < @shuffledNegVectorData) ? $posHairpinCount : $posHairpinCount - @shuffledNegVectorData;
    for (my $i = 0; $i < $negHairpinCount; $i++) {
	print HVFA pop(@shuffledNegVectorData) . "\n";
    }
    close(HVFA);
}

sub printHairpinVectorFile {
    my($hairpinsFile,$hairpinVectorFile,$genomeDir,$chromLengths,$class,$parameters) = @_;
    open(HVF,">".$hairpinVectorFile) or die "failed to open $hairpinVectorFile for writing\n";
    print HVF "GeneId\ty\ttotalSense\ttotalAntisense\tmfe\taapd\ttapd\turf\tahc\tafh\tpbp\tsameShift\tbothShift\tSPA\tOPA\t";
    print HVF "aa\tac\tag\tat\tca\tcc\tcg\tct\tga\tgc\tgg\tgt\tta\ttc\ttg\ttt\tduplexEnergy\tGCcontent\t";
    #print HVF "zScore\t";
    print HVF "foldDupCmp\texactFoldDupCmp\tdupPBP\tdupLoopLength";
    print HVF "\tAPV\twAPV";
    print HVF "\tARV\twARV";
    print HVF "\tmpCount\tdupCount\tdupOverlap\tmpLoopDistance\tdupLoopDistance\tloopSize\tmpOverlapScore\tdupOverlapScore\twghtMPOverlapIntensity\twghtDupOverlapIntensity\tOPA2\ttotalOverlapAmount\taverageOverlapAmount\ttotalRelativeOverlapAmount\taverageRelativeOverlapAmount\tmaxOverlap\tinnerLoopGapCount\ttotalSenseRPM\ttotalAntisenseRPM";
    print HVF "\tRFProductAvg\n";
    my $nonDuplicateFolds = dropDuplicateFolds($hairpinsFile);
    open(HPF, $hairpinsFile) or die "failed to open $hairpinsFile\n";
    while (<HPF>) {
        chomp;
	unless(/\#/) {	  
	    my($name,$chrom,$start,$stop,$strand,$leftCenter,$rightCenter,$totalSense,$totalAntisense,$mfe,$sequence,$fold,$aapd,$tapd,$urf,$ahc,$afh,$pbp,$sameShift,$bothShift,$spa,$opa,$aa,$ac,$ag,$at,$ca,$cc,$cg,$ct,$ga,$gc,$gg,$gt,$ta,$tc,$tg,$tt,$duplexEnergy,$GCcontent,$foldDupCmp,$exactFoldDupCmp,$dupPBP,$dupLoopLength,$apv,$wApv,$arv,$wArv,$maxProductCount,$mpDuplexCount,$duplexOverlap,$mpLoopDistance,$dupLoopDistance,$loopSize,$mpOverlapScore,$dupOverlapScore,$wghtMPOverlapIntensity,$wghtDupOverlapIntensity,$opa2,$totalOverlapAmount,$averageOverlapAmount,$totalRelativeOverlapAmount,$averageRelativeOverlapAmount,$maxOverlap,$innerLoopGapCount,$totalSenseRPM,$totalAntisenseRPM,$RFProductAvg) = split(/\t/,$_);
	    if($nonDuplicateFolds->{$name}) {
		my $location = "$chrom:$start..$stop:$strand";
		print HVF "$name" . "_" . "$location\t$class\t$totalSense\t$totalAntisense\t$mfe\t$aapd\t$tapd\t$urf\t$ahc\t$afh\t$pbp\t$sameShift\t$bothShift\t$spa\t$opa\t";
		print HVF "$aa\t$ac\t$ag\t$at\t$ca\t$cc\t$cg\t$ct\t$ga\t$gc\t$gg\t$gt\t$ta\t$tc\t$tg\t$tt\t";
		print HVF "$duplexEnergy\t$GCcontent";
#		print HVF "\t$zScore";
		print HVF "\t$foldDupCmp\t$exactFoldDupCmp\t$dupPBP\t$dupLoopLength";
		print HVF "\t$apv\t$wApv\t$arv\t$wArv";
		print HVF "\t$maxProductCount\t$mpDuplexCount\t$duplexOverlap\t$mpLoopDistance\t$dupLoopDistance\t$loopSize\t$mpOverlapScore\t$dupOverlapScore\t$wghtMPOverlapIntensity\t$wghtDupOverlapIntensity\t$opa2\t$totalOverlapAmount\t$averageOverlapAmount\t$totalRelativeOverlapAmount\t$averageRelativeOverlapAmount\t$maxOverlap\t$innerLoopGapCount\t$totalSenseRPM\t$totalAntisenseRPM";
		print HVF "\t$RFProductAvg";
		print HVF "\n";
	    }
	}
    }
    close(HPF);
    close(HVF);
}


##################################################
# MACHINE LEARNING Gene Regions DATA SUBROUTINES #
##################################################

sub getAnnotatedMirReadRegions {
    my($mirHairpins,$readRegionsFile,$parameters) = @_;
    my %mirRRHash;
    my %annotatedMirBuckets;
    my $bucketSize = $parameters->{bucketSize} or die "bucket size not defined\n";
    foreach my $chrom (keys %{$mirHairpins}) {
	foreach my $hairpinInfo (@{$mirHairpins->{$chrom}}) {
	    my($start,$stop,$strand,$id,$name) = @{$hairpinInfo};
	    my $bucket = int($start / $bucketSize);
	    my $nextBucket = int($stop / $bucketSize);
	    push(@{$annotatedMirBuckets{$chrom}{$bucket}}, [$start,$stop,$strand,$id,$name]);
	    #print "pushing $name into bucket $bucket\n";
	    if ($nextBucket > $bucket) {
		#we will have some sequences in two buckets but it doesn't really matter for this function
		push(@{$annotatedMirBuckets{$chrom}{$nextBucket}}, [$start,$stop,$strand,$id,$name]);
	    }
	}
    }
    open(RR,$readRegionsFile) or die "failed to open $readRegionsFile\n";
    while (<RR>) {
	chomp;
	my($chrom,$rrStart,$rrStop,$rrName,$size,$rrStrand) = split("\t",$_);
	my $bucket = int($rrStart / $bucketSize);
	if ($annotatedMirBuckets{$chrom}{$bucket}) {
	    foreach my $hairpinInfo (@{$annotatedMirBuckets{$chrom}{$bucket}}) {
		my($start,$stop,$strand,$id,$name) = @{$hairpinInfo};
		if (getOverlap($rrStart,$rrStop,$start,$stop)) {
		    #print "pushing $rrName onto $name\n";
		    $mirRRHash{$rrName} = 1;
		}		    
	    }
	}
    }
    close(RR);
    return \%mirRRHash;
}


sub getGeneRegionRRCounts {
    my($grLocations,$readRegionsFile,$mirRRHash,$parameters) = @_;
    my $bucketSize = $parameters->{bucketSize} or die "bucket size not defined\n";
    #print "$bucketSize\n";
    my %grSequenceBuckets;
    my(%grTotalRRLoci,%grTotalMirLoci,%grTotalNonMirLoci,%grTotalRRSize,%grTotalMirSize,%grTotalNonMirSize,%hpTotalRRLoci,%hpTotalRRSize);
    my(%grLeftLoci,%grRightLoci,%grLeftSize,%grRightSize,%symmetryByLoci,%symmetryBySize,%grASLeftLoci,%grASRightLoci,%grASLeftSize,%grASRightSize,%crossSymmetryByLoci,%crossSymmetryBySize);
    foreach my $chrom (keys %{$grLocations}) {
	foreach my $geneId (keys %{$grLocations->{$chrom}}) {
	    my $location = $grLocations->{$chrom}{$geneId};
	    my($chrom,$grStart,$grStop,$strand) = parseLocation($location);
	    my $bucket = int($grStart / $bucketSize);
	    my $nextBucket = int($grStop / $bucketSize);
	    push(@{$grSequenceBuckets{$chrom}{$bucket}}, [$geneId,$grStart,$grStop]);
	    #print "pushing $geneId into bucket $bucket\n";
	    if ($nextBucket > $bucket) {
		#we will have some sequences in two buckets but it doesn't really matter for this function
		push(@{$grSequenceBuckets{$chrom}{$nextBucket}}, [$geneId,$grStart,$grStop]);
	    }
	    $grTotalRRLoci{$geneId} = 0;
	    $grTotalRRSize{$geneId} = 0;
	    $grTotalMirLoci{$geneId} = 0;
	    $grTotalMirSize{$geneId} = 0;
	    $grTotalNonMirLoci{$geneId} = 0;
	    $grTotalNonMirSize{$geneId} = 0;
	    $grLeftLoci{$geneId} = 0;
	    $grLeftSize{$geneId} = 0;
	    $grRightLoci{$geneId} = 0;
	    $grRightSize{$geneId} = 0;
	    $symmetryByLoci{$geneId} = 0;
	    $symmetryBySize{$geneId} = 0;
	    $grASLeftLoci{$geneId} = 0;
	    $grASLeftSize{$geneId} = 0;
	    $grASRightLoci{$geneId} = 0;
	    $grASRightSize{$geneId} = 0;
	    $crossSymmetryByLoci{$geneId} = 0;
	    $crossSymmetryBySize{$geneId} = 0;
	}
    }
    open(RR,$readRegionsFile) or die "failed to open $readRegionsFile\n";
    while (<RR>) {
	chomp;
	my($chrom,$rrStart,$rrStop,$rrName,$size,$rrStrand) = split("\t",$_);
	my $bucket = int($rrStart / $bucketSize);
	if ($grSequenceBuckets{$chrom}{$bucket}) {
	    #print "in bucket\n";
	    foreach my $geneRegionInfo (@{$grSequenceBuckets{$chrom}{$bucket}}) {
		my($geneId,$grStart,$grStop) = @{$geneRegionInfo};
		my($hpName,$hpLocation) = parseGeneId($geneId);
		my($hpChrom,$hpStart,$hpStop,$hpStrand) = parseLocation($hpLocation);
		if (getOverlap($rrStart,$rrStop,$hpStart,$hpStop)) {
		    #print "$rrName overlaps hairpin\n";
		    $hpTotalRRLoci{$geneId} += 1;
		    $hpTotalRRSize{$geneId} += $size;
		} elsif (getOverlap($rrStart,$rrStop,$grStart,$grStop)) {
		    #print "$rrName overlaps gene region\n";
		    $grTotalRRLoci{$geneId} += 1;
		    $grTotalRRSize{$geneId} += $size;
		    if ($mirRRHash->{$rrName}) {
			#print "$rrName is a mir region\n";
			$grTotalMirLoci{$geneId} += 1;
			$grTotalMirSize{$geneId} += $size;
		    } else {
			#print "$rrName is in gene region\n";
			$grTotalNonMirLoci{$geneId} += 1;
			$grTotalNonMirSize{$geneId} += $size;
		    }
		    #collecting left and right RR info in order to determine the degree of symmetry of the reads later
		    if ($rrStop < $hpStart && $rrStrand eq $hpStrand) {
			$grLeftLoci{$geneId} += 1;
			$grLeftSize{$geneId} += $size;
		    } elsif ($rrStart > $hpStop && $rrStrand eq $hpStrand) {
			$grRightLoci{$geneId} += 1;
			$grRightSize{$geneId} += $size;
		    } elsif ($rrStop < $hpStart && $rrStrand ne $hpStrand) {
			$grASLeftLoci{$geneId} += 1;
			$grASLeftSize{$geneId} += $size;
		    } elsif ($rrStart > $hpStop && $rrStrand ne $hpStrand) {
			$grASRightLoci{$geneId} += 1;
			$grASRightSize{$geneId} += $size;
		    }
		}
	    }
	}
    }
    close(RR);    
    foreach my $geneId (keys %hpTotalRRLoci) {
	if ($grLeftLoci{$geneId} == 0 && $grRightLoci{$geneId} == 0) {
	    $symmetryByLoci{$geneId} = 1;
	} else {
	    $symmetryByLoci{$geneId} = ($grLeftLoci{$geneId} < $grRightLoci{$geneId}) ? $grLeftLoci{$geneId} / $grRightLoci{$geneId} : $grRightLoci{$geneId} / $grLeftLoci{$geneId};
	}
	if ($grLeftSize{$geneId} == 0 && $grRightSize{$geneId} == 0) {
	    $symmetryBySize{$geneId} = 1;
	} else {
	    $symmetryBySize{$geneId} = ($grLeftSize{$geneId} < $grRightSize{$geneId}) ? $grLeftSize{$geneId} / $grRightSize{$geneId} : $grRightSize{$geneId} / $grLeftSize{$geneId};
	}
	if ($grLeftLoci{$geneId} == 0 && $grASRightLoci{$geneId} == 0 && $grASLeftLoci{$geneId} == 0 && $grRightLoci{$geneId} == 0) {
	    $crossSymmetryByLoci{$geneId} = 1;
	} else {
	    my $forwardCross = $grLeftLoci{$geneId} + $grASRightLoci{$geneId};
	    my $backCross = $grASLeftLoci{$geneId} + $grRightLoci{$geneId};
	    $crossSymmetryByLoci{$geneId} = ($forwardCross < $backCross) ? $forwardCross / $backCross : $backCross / $forwardCross;
	}
	if ($grLeftSize{$geneId} == 0 && $grASRightSize{$geneId} == 0 && $grASLeftSize{$geneId} == 0 && $grRightSize{$geneId} == 0) {
	    $crossSymmetryBySize{$geneId} = 1;
	} else {
	    my $forwardCross = $grLeftSize{$geneId} + $grASRightSize{$geneId};
	    my $backCross = $grASLeftSize{$geneId} + $grRightSize{$geneId};
	    $crossSymmetryBySize{$geneId} = ($forwardCross < $backCross) ? $forwardCross / $backCross : $backCross / $forwardCross;
	}

    }
    return(\%grTotalRRLoci,\%grTotalMirLoci,\%grTotalNonMirLoci,\%grTotalRRSize,\%grTotalMirSize,\%grTotalNonMirSize,\%hpTotalRRLoci,\%hpTotalRRSize,\%symmetryByLoci,\%symmetryBySize,\%crossSymmetryByLoci,\%crossSymmetryBySize);
}

sub getGeneRegionRRData {
    my($hairpins,$grSpan,$readRegionsFile,$mirHairpins,$chromLengths,$parameters) = @_;
    my $mirRRHash = getAnnotatedMirReadRegions($mirHairpins,$readRegionsFile,$parameters);
    my(%grRRInfo);
    my(%totalLoci,%totalMirLoci,%totalNonMirLoci,%totalSize,%totalMirSize,%totalNonMirSize,%hpTotalLoci,%hpTotalSize,%symByLoci,%symBySize);
    for (my $i = 0; $i < @{$grSpan}; $i++) {
	my $span = $grSpan->[$i];
	$parameters->{totalLength} = $span;
	my %grLocations;	
	#print "collecting for $grSpan->[$i]\n";
	foreach my $chrom (keys %{$hairpins}) {
	    #print "chrom = $chrom\n";
	    my $chromLength = $chromLengths->{$chrom};
	    #my $genome = loadGenome("$genomeDir/$chrom.fa");
	    foreach my $hairpinInfo (@{$hairpins->{$chrom}}) {
		my($start,$stop,$strand,$geneId) = @{$hairpinInfo};
		#print "geneId = $geneId\n";
		my $grLocation = extendSequenceToFillBuffer($chrom,$start,$stop,$strand,$chromLength,$parameters);
		$grLocations{$chrom}{$geneId} = $grLocation;
		#print "grLocation = $grLocation\n";
	    }
	}
	my($grTotalRRLoci,$grTotalMirLoci,$grTotalNonMirLoci,$grTotalRRSize,$grTotalMirSize,
	   $grTotalNonMirSize,$hpTotalRRLoci,$hpTotalRRSize,$symmetryByLoci,$symmetryBySize,
	   $crossSymmetryByLoci,$crossSymmetryBySize) = getGeneRegionRRCounts(\%grLocations,$readRegionsFile,$mirRRHash,$parameters);
	foreach my $geneId (keys %{$hpTotalRRLoci}) {
	    @{$grRRInfo{$geneId}{$span}} = ($grTotalRRLoci->{$geneId}, $grTotalMirLoci->{$geneId}, $grTotalNonMirLoci->{$geneId}, $grTotalRRSize->{$geneId}, $grTotalMirSize->{$geneId}, $grTotalNonMirSize->{$geneId}, $hpTotalRRLoci->{$geneId}, $hpTotalRRSize->{$geneId}, $symmetryByLoci->{$geneId}, $symmetryBySize->{$geneId}, $crossSymmetryByLoci->{$geneId}, $crossSymmetryBySize->{$geneId});
	}
    }
    return \%grRRInfo;
}

sub generateGeneRegionsTrainData {
    my($hpTrainFile,$productRFFile,$productFile,$readRegionsFile,$chromLengths,$mirBaseGff,$genomeDir,$parameters) = @_;
    my @grSpan = (1000,2000,3000,4000);  #must be sorted from least to greatest
    $parameters->{bucketSize} = 100000;
    my $geneRegionTrainFile = "geneTrainFile.txt";
    #my $hpTrainFile = "trainTest.txt";
    my %hairpins;
    my %hairpinTrainFileData;
    my %grSequenceData;
    my %grWFC;
    my %grInfo;
    my($mirHairpins) = readMirbaseGff3($mirBaseGff);
    open(TF,$hpTrainFile) or die "failed to open $hpTrainFile\n";
    <TF>; #skipping first line
    while (<TF>) {
	chomp;
	my($geneId,$class,$totalSense,$totalAntisense,$mfe,$aapd,$tapd,$urf,$ahc,$afh,$pbp,$sameShift,$bothShift,$SPA,$OPA,$aa,$ac,$ag,$at,
	   $ca,$cc,$cg,$ct,$ga,$gc,$gg,$gt,$ta,$tc,$tg,$tt,$duplexEnergy,$gcContent,$zScore,$RFProductAvg) = split("\t",$_);
	my($name,$location) = parseGeneId($geneId);
	my($chrom,$start,$stop,$strand) = parseLocation($location);
	#push(@{$hairpins{$chrom}}, [$start,$stop,$strand,$geneId]);
	push(@{$hairpins{$class}{$chrom}}, [$start,$stop,$strand,$geneId]);
    }
    close(TF);
    open(GRTF,">$geneRegionTrainFile") or die "failed to open $geneRegionTrainFile for writing\n";
    #creating file header
    print GRTF "GeneId\ty";
    foreach my $span (@grSpan) {
	print GRTF "\tnonMirLoci$span\tneigToHPSize$span";
	print GRTF "\tsymLoci$span\tsymSize$span\tcrossSymLoci$span\tcrossSymSize$span";
    }
    for (my $i = 0; $i < @grSpan; $i++) {
	for (my $j = $i + 1; $j < @grSpan; $j++) {
	    my $span1 = $grSpan[$i];
	    my $span2 = $grSpan[$j];
	    print GRTF "\tnonMirLoci$span1\_$span2\tnonMirSize$span1\_$span2";
	}
    }
    foreach my $span (@grSpan) {
	print GRTF "\taa\_$span\tac\_$span\tag\_$span\tat\_$span\tca\_$span\tcc\_$span\tcg\_$span\tct\_$span\tga\_$span\tgc\_$span\tgg\_$span\tgt\_$span\tta\_$span\ttc\_$span\ttg\_$span\ttt\_$span";
	print GRTF "\tWFC\_$span";
    }
    print GRTF "\n";
    #adding vector info
    foreach my $class (keys %hairpins) {
	my $grRRInfo = getGeneRegionRRData($hairpins{$class},\@grSpan,$readRegionsFile,$mirHairpins,$chromLengths,$parameters);
	foreach my $chrom (keys %{$hairpins{$class}}) {
	    my $genome = loadGenome("$genomeDir/$chrom.fa");
	    foreach my $geneInfo (@{$hairpins{$class}{$chrom}}) {
		my($start,$stop,$strand,$geneId) = @{$geneInfo};
		print GRTF "$geneId\t$class";
		foreach my $span (@grSpan) {
		    my($totalLoci,$totalMirLoci,$totalNonMirLoci,$totalSize,$totalMirSize,$totalNonMirSize,$hpTotalLoci,$hpTotalSize,$symByLoci,$symBySize,$crossSymByLoci,$crossSymBySize) = @{$grRRInfo->{$geneId}{$span}};
		    my $neighborToHPSizeRatio = $totalNonMirSize / $hpTotalSize;
		    #I'm holding off on analyzing information for neighboring miRs for now since this program will eventually be used for genomes with no annotated mirs
		    print GRTF "\t$totalNonMirLoci\t$neighborToHPSizeRatio";
		    print GRTF "\t$symByLoci\t$symBySize\t$crossSymByLoci\t$crossSymBySize";
		}
		for (my $i = 0; $i < @grSpan; $i++) {
		    for (my $j = $i + 1; $j < @grSpan; $j++) {
			my $span1 = $grSpan[$i];
			my $span2 = $grSpan[$j];
			my($totalLoci1,$totalMirLoci1,$totalNonMirLoci1,$totalSize1,$totalMirSize1,$totalNonMirSize1,
			   $hpTotalLoci1,$hpTotalSize1,$symByLoci1,$symBySize1,$crossSymByLoci1,$crossSymBySize1) = @{$grRRInfo->{$geneId}{$span1}};
			my($totalLoci2,$totalMirLoci2,$totalNonMirLoci2,$totalSize2,$totalMirSize2,$totalNonMirSize2,
			   $hpTotalLoci2,$hpTotalSize2,$symByLoci2,$symBySize2,$crossSymByLoci2,$crossSymBySize2) = @{$grRRInfo->{$geneId}{$span2}};
			#the reason for using the smaller span over the larger span is to avoid division by zero
			my $nonMirLociRatio = ($totalNonMirLoci2 == 0 && $totalNonMirLoci1 == 0) ? 1 : $totalNonMirLoci1 / $totalNonMirLoci2;
			my $nonMirSizeRatio = ($totalNonMirSize2 == 0 && $totalNonMirSize1 == 0) ? 1 : $totalNonMirSize1 / $totalNonMirSize2;
			print GRTF "\t$nonMirLociRatio\t$nonMirSizeRatio";
		    }
		}
		foreach my $span (@grSpan) {
		    $parameters->{totalLength} = $span;
		    my $grLocation = extendSequenceToFillBuffer($chrom,$start,$stop,$strand,$chromLengths->{$chrom},$parameters);
		    my $grSequence = getSequence($grLocation,$genome);
		    my($aa,$ac,$ag,$at,$ca,$cc,$cg,$ct,$ga,$gc,$gg,$gt,$ta,$tc,$tg,$tt) = getDinucleotideFrequencies($grSequence);
		    my $WFC = getWFC_Stirling($grSequence);
		    print GRTF "\t$aa\t$ac\t$ag\t$at\t$ca\t$cc\t$cg\t$ct\t$ga\t$gc\t$gg\t$gt\t$ta\t$tc\t$tg\t$tt";
		    print GRTF "\t$WFC";
		}
		print GRTF "\n";
	    }
	}
    }
    close(GRTF);
}


#####################
# INPUT SUBROUTINES #
#####################

sub loadLibrarySizes {
    my($librarySizesFile) = @_;
    my %librarySizes;
    open(LS,$librarySizesFile) or die "failed to open $librarySizesFile for reading\n";
    while (<LS>) {
	chomp;
	my($sample,$readCount) = split("\t");
	$librarySizes{$sample} = $readCount;
    }
    close(LS);
    return \%librarySizes;
}

sub loadRepeatRegions {
    my $repeatRegionsFile = shift;
    my $chromSizes = shift;
    my $repeatRegions = {};
    my $fileType = checkFileFormat($repeatRegionsFile, $chromSizes);
    if ($fileType eq "fileList") {
	open(RRFL,$repeatRegionsFile) or die "could not open $repeatRegionsFile\n";
	while(<RRFL>) {
	    chomp;
	    my $repeatRegionsFile = $_;
	    my $listedFileType = checkFileFormat($repeatRegionsFile, $chromSizes);
	    $repeatRegions = readRepeatRegionsFile($repeatRegionsFile, $repeatRegions,  $listedFileType);
	}
	close(RRFL);
    } else {
	$repeatRegions = readRepeatRegionsFile($repeatRegionsFile, $repeatRegions, $fileType);
    }
    return $repeatRegions;
}

sub readRepeatRegionsFile {
    my $repeatRegionsFile = shift;
    my $repeatRegions = shift;
    my $fileType = shift;
    open(RRF, $repeatRegionsFile) or die "Could not read repeat regions file:$repeatRegionsFile\n";
    while (<RRF>) {
	if ($fileType eq 'txt' || $fileType eq 'bed') {
	    while (<RRF>) {
		chomp;
		my($chrom,$start,$stop) = split(/\t/,$_);
		push(@{$repeatRegions->{$chrom}},[$start,$stop]);
	    }
	} elsif ($fileType eq 'gff') {
	    while (<RRF>) {
		chomp;
		unless (/\#/) {
		    my ($chrom, $source, $type, $start, $stop, $score, $strand, $frame, $attribute) = split(/\t/);
		    push(@{$repeatRegions->{$chrom}},[$start,$stop]);
		}
	    }
	} else {
	    die print "unknown format of $repeatRegionsFile.  mirTRAP readRepeatRegions takes in .txt and .gff files.\n";
	}
    }
    close(RRF);
    return $repeatRegions;
}

sub readMirbaseGff3 {
    my($mirbaseGff3) = @_;
    my %hairpins;
    my %products;
    open(MBGFF3,$mirbaseGff3) or die "could not open $mirbaseGff3\n";
    while(<MBGFF3>) {
	unless(/^\#/) {
	    chomp;
	    my($chrom,$source,$type,$start,$stop,$score,$strand,$phase,$info) = split(/\t/);
	    my %info;
	    my @terms = split(/;/,$info);
	    foreach my $term (@terms) {
		my($key,$value) = $term =~ /(.*)=(.*)/;
		$info{$key} = $value;
	    }
	    if($type eq "miRNA_primary_transcript") {
		# hairpin region from the gff file
		my $id = $info{ID} or die "No ID found for the line:\n$_\n"; 
		my $name = $info{Name} or die "No Name found for the line:\n$_\n"; 
		push(@{$hairpins{$chrom}},[$start,$stop,$strand,$id,$name]);
	    }
	    if($type eq "miRNA") {
		# mature product from the hairpin file
		my $id = $info{ID} or die "No ID found for the line:\n$_\n";
                my $name = $info{Name} or die "No Name found for the line:\n$_\n";
		my $parentId = $info{Derives_from} or die "No Derives_from found for the line:\n$_\n";
		push(@{$products{$parentId}},[$chrom,$start,$stop,$strand,$id,$name]);
	    }
	}
    }
    close(MBGFF3);
    return(\%hairpins,\%products);
}

sub getLoopsFromMirbaseGff {
    my($mirBaseGff) = @_;
    my($hairpins,$products) = readMirbaseGff3($mirBaseGff);
    my %positiveList;
    foreach my $hairpinChrom (keys %{$hairpins}) {
	foreach my $hairpinInfo (@{$hairpins->{$hairpinChrom}}) {
	    my($hairpinStart,$hairpinStop,$strand,$hairpinId,$hairpinName) = @{$hairpinInfo};
	    if ($products->{$hairpinId}) {
		if (@{$products->{$hairpinId}} == 2) {
		    my($loopStart, $loopStop, $i) = (0,0,0);
		    foreach my $mirProduct (@{$products->{$hairpinId}}) {
			my($chrom,$start,$stop,$strand,$id,$name) = @{$mirProduct};
			if ($i == 1) {
			    $loopStop = $start;
			} else {
			    $loopStart = $stop;
			}
			$i++;
		    }
		    my $annotatedLocation = "$hairpinChrom:$loopStart..$loopStop:$strand";
		    push(@{$positiveList{$hairpinChrom}}, ["$hairpinName-loop",$annotatedLocation]);
		} elsif (@{$products->{$hairpinId}} > 2) {
		    print "$hairpinName has more than 2 products\n";
		}
	    }
	}
    }
    return (\%positiveList);    
}

sub checkFileFormat {
    my $fileName = shift;
    my $chromSizes = shift;
    my @chroms = keys %{$chromSizes};
    my $fileType;

    if (($fileType) = $fileName =~ /\.(bed|gff)$/) {
	return $fileType;
    }
    open(FP, $fileName) or die "failed to open $fileName in miRTRAP1_7::checkFileFormat";
    while (<FP>) {
	unless (/\#/) {
	    my @column = split;
	    if ($chromSizes->{$column[0]} && isInt($column[1]) && isInt($column[2])) {
		close(FP);
		return "txt";
	    } elsif (@column == 1) {
		close(FP);
		return "fileList";
	    } else {
		close(FP);
		die print "unable to determine file type for $fileName fileType\n";
	    }
	}
    }
    die print "unable to determine file type for $fileName\n";
}

sub readtRNAScanFile {
    my $tRNAScanFile = shift;
    my %tRNAs;
    my $RECORD = 0;
    if($tRNAScanFile) {
        open(TRNA,$tRNAScanFile) or die "failed to open $tRNAScanFile for reading in readtRNAScanFile()";
        while(<TRNA>) {
            chomp;
            my($id,$count,$tStart,$tStop,$type,$antiCodon,$seqId,$length,$score) = split(/\s+/);
            if($RECORD) {
                $tRNAs{$id} = 1;
            }
            if($id =~ /\-+/) {
                $RECORD = 1;
            }
        }
    }
    return \%tRNAs;
}


1;
__END__


=head1 NAME

miRTRAP (miR Tests for Read Analysis and Prediction) - Perl module (and script) for discovering microRNAs from highthroughput sequencing data.

=head1 SYNOPSIS

#the following code is the contents of the associated "printReadRegions.pl" script.
use miRTRAP;
my $usage = "Usage:\n\n$0 <configFile>\n";

my $parameters = {
    "maxLength" => 160,
    "maxCount" => 5,       
    "maxHitCount" => 50,
    "filePrefix" => "readRegions",
    "readListFile" => "",
    "repeatRegionsFile" => ""   
};

my $configFile = $ARGV[0] or die $usage;
miRTRAP::readConfigFile($configFile,$parameters);
my $readGffListFile = $parameters->{readListFile} or die "FAIL: readListFile not loaded in configFile.\n";
my $repeatRegionsFile = $parameters->{repeatRegionsList};
my($dataSet,$sampleTotal,$sampleList,$hitCount) = miRTRAP::loadReadGffList($readGffListFile,$parameters);
my $repeatRegions = {};
if($repeatRegionsFile) {
    $repeatRegions = miRTRAP::readRepeatRegions($repeatRegionsFile);
}
miRTRAP::processOverlaps($dataSet,$repeatRegions,$hitCount,$parameters);

=head1 DESCRIPTION

MicroRNAs (miRs) have been broadly implicated in animal development and disease. However, the systematic, whole-genome identification of miRs is complicated by their small size.  Current identification methods rely mainly on the projected stability of putative stem-loop structures (pre-miRs) and on sequence comparisons with known miRs. 

This method, miRTRAP, incorporates the mechanisms of miR biogenesis, and includes additional criteria regarding the prevalence and quality of small RNAs arising from the antisense strand and neighboring loci.

Please see associate scripts and README file for more information and examples.

=head2 EXPORT

None by default.

=head1 SEE ALSO

http://flybuzz.berkeley.edu/miRTRAP.html

=head1 AUTHOR

David Hendrix, davidhendrix@berkeley.edu

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010 by David Hendrix

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

=cut

#  LocalWords:  printProdTrainingData returnLargestProductOverlap
#  LocalWords:  posOverlapProducts
