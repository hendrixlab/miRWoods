package miRBamTools;
use Bio::DB::Bam::Alignment;
use Math::CDF;
use strict; 
return 1;

sub readConfigFile {
    my($configFile,$parameters) = @_;
    open(CONFIGFILE,$configFile) or die "FAIL: could not open $configFile\n";
    while(<CONFIGFILE>) {
	chomp;
	my($key,$value) = split(/\s+=\s+/);
	$parameters->{$key} = $value;
    }
    my $filePrefix = $parameters->{filePrefix} or die "FAIL: filePrefix not loaded.\n";
    $parameters->{distancesFile} = $filePrefix."_distances.txt";
    $parameters->{readRegionsFile} = $filePrefix.".txt";
    $parameters->{longReadRegionsFile} = $filePrefix."_long.txt";
    $parameters->{allReadRegionsFile} = $filePrefix."_full.txt";
    $parameters->{readRegionsFastaFile} = $filePrefix . ".fasta";
    $parameters->{rnaFoldOutput} = $filePrefix . ".mfe";
    $parameters->{filteredCandidatesFile} = $filePrefix . "_filteredCandidates.txt";
    $parameters->{hairpinsFile} = $filePrefix . "_hairpins.txt";
    $parameters->{distinctReadsFile} = $filePrefix . "_distinctReads.txt";
    $parameters->{productFile} = $filePrefix . "_products.txt";    
    $parameters->{configFile} = $configFile;
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
    $seq =~ tr/acgtACGT/tgcaTGCA/;
    $seq=reverse($seq);
    return $seq;
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
    my $chromRegExp = shift;
    
    my($chromRegExp1,$chromRegExp2);
    if($chromRegExp) {
	($chromRegExp1,$chromRegExp2) = $chromRegExp =~ /s\/(.*)\/(.*)\// or
	    die "regular expression must be of the form \"s/a/b/\"\nto change from 'a' to 'b'\n";;
    }
    my %sequences;
    my $tag;
    open(GFILE,$genomeFile) or die "could not open $genomeFile\n";
    while(<GFILE>) {
        chomp;
        if(/>/) {
            s/>//g;
            my @terms = split;
	    $tag = shift(@terms);
	    if($chromRegExp) {		
		$tag =~ s/$chromRegExp1/$chromRegExp2/;
	    }
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
	# gff is 1-based as is bam
	# assuming the coord is 1-based so subtract 1 to get a chromosome sequence substring
	# for a given genomic coordinate.
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

sub binomialPValue {
    my($expCounts,$ctlCounts) = @_;
    if(my $cdf = Math::CDF::pbinom($expCounts-1,int($expCounts+$ctlCounts),0.5)) {
	return 1 - $cdf;
    } else {
        print STDERR "...Could not get cdf for $expCounts $ctlCounts\n";
    }
    return 1.0;
}

sub readLinearRegressionFile {
    my($linearRegressionFile) = @_;
    open(SF,$linearRegressionFile) or die "Could not open $linearRegressionFile\n.";
    my $A;
    my $B;
    while(<SF>) {
        chomp;
        if(/a=/) {
            ($A) = $_ =~ /a=(.*)/;
        } elsif(/b=/) {
            ($B) = $_ =~ /b=(.*)/;
        }
    }
    unless((defined $A)&&(defined $B)) {die "No values for a and b? $A $B\n";}
    return ($A,$B);
}

#########################
# FOLD PROCESSING TOOLS #
#########################

sub computeProductLinearRegression {
    my($bamExpCtl,$chromLengths,$parameters) = @_;    
    my $minDist = $parameters->{minDist} or die "minDist: parameter not loaded.";
    my $minCount = $parameters->{minCount} or die "minCount: parameter not loaded.";
    my $minHet = $parameters->{minHet} or die "minHet: parameter not loaded.";
    my $maxLength = $parameters->{maxLength} or die "maxLength: parameter not loaded.";
    my($expLabel,$expBam,$expTotal,$expFile) = @{$bamExpCtl->{"exp"}};
    my($ctlLabel,$ctlBam,$ctlTotal,$ctlFile) = @{$bamExpCtl->{"ctl"}};
    my @strands = ("+","-");
    my @chroms = sort keys(%{$chromLengths});
    my $sumX = 0;
    my $sumY = 0;
    my $sumX2 = 0;
    my $sumY2 = 0;
    my $sumXY = 0;
    my $N = 0;
    my $windowSize = 10000;
    foreach my $chrom (@chroms) {
	foreach my $strand (@strands) {
	    # for the chrom/strand identify all distinct start/stop positions.
	    for(my $startPos = 0; $startPos <= $chromLengths->{$chrom}-$windowSize+1; $startPos += $windowSize) {
		my %distinctReadCount;
		my @reads = $expBam->get_features_by_location(-seq_id => $chrom,
							      -start => $startPos,
							      -end => $startPos + $windowSize);		
		foreach my $read (@reads) {
		    my $rStart = $read->start;
		    my $rStop = $read->end;
		    my $rStrand = mapBamStrand($read->strand);
		    if($strand eq $rStrand) {
			$distinctReadCount{$rStart."x".$rStop}++;
		    }
		}
		# repeat for the ctl reads
		@reads = $ctlBam->get_features_by_location(-seq_id => $chrom,
							      -start => $startPos,
							      -end => $startPos + $windowSize);		
		foreach my $read (@reads) {
		    my $rStart = $read->start;
		    my $rStop = $read->end;
		    my $rStrand = mapBamStrand($read->strand);
		    if($strand eq $rStrand) {
			$distinctReadCount{$rStart."x".$rStop}++;
		    }
		}
		# for each start/stop position, sort by count.
		my @products;
		foreach my $key (sort {$distinctReadCount{$b} <=> $distinctReadCount{$a}} keys %distinctReadCount) {
		    my($drStart,$drStop) = split(/x/,$key);
		    unless(overlapsProducts(\@products,$drStart,$drStop,$strand,$parameters)) {
			push(@products,[$drStart,$drStop,$distinctReadCount{$key}]);
		    }
		}
		%distinctReadCount = ();
		foreach my $product (@products) {
		    my($pStart,$pStop) = @{$product};		
		    if($pStop - $pStart <= $maxLength) {
			my($expCount,$expFiveHet,$expThreeHet) = getBamProductCount($expBam,$chrom,$pStart,$pStop,$strand,$parameters);
			my($ctlCount,$ctlFiveHet,$ctlThreeHet) = getBamProductCount($ctlBam,$chrom,$pStart,$pStop,$strand,$parameters);
			$sumX += $ctlCount;
			$sumY += $expCount;
			$sumX2 += $ctlCount*$ctlCount;
			$sumY2 += $expCount*$expCount;	    
			$sumXY += $ctlCount*$expCount;
			$N++;
		    }
		}
		@products = ();
	    }
	}
    }
    my $xAvg = $sumX / $N;
    my $yAvg = $sumY / $N;
    my $B = ($N*$sumXY - $sumX*$sumY)/($N*$sumX2 - $sumX*$sumX);
    my $A = $yAvg - $B * $xAvg;
    my $scaleFactorFile = $expLabel."_vs_".$ctlLabel."_linReg.txt";
    open(SFF,">$scaleFactorFile") or die "could not open $scaleFactorFile for writing.\n";
    print SFF "a=$A\nb=$B";
    close(SFF);
    return($A,$B);
}

sub getReplicateProductCounts {
    my($products,$bamFile,$parameters) = @_;
    my($fileBase) = $bamFile =~ /([^\/]*)\.bam/;
    my $outputFile = $fileBase . "_productReplicate.txt";
    my $bam = loadBamFile($bamFile);
    open(OUTPUT,">$outputFile");
    foreach my $chrom (keys %{$products}) {
	foreach my $product (@{$products->{$chrom}}) {
	    my($id,$location,$count,$fiveHet) = @{$product};
	    my($chr,$start,$stop,$strand) = parseLocation($location);
	    my($expCount,$expFiveHet,$expThreeHet) = getBamProductCount($bam,$chrom,$start,$stop,$strand,$parameters);
	    print OUTPUT "$id\t$location\t$expCount\t$expFiveHet\t$expThreeHet\n";
	}
    }
    close(OUTPUT);
}

sub getReplicateMatureExpressionCounts {
    my($products,$bamFile,$parameters) = @_;
    my($fileBase) = $bamFile =~ /([^\/]*)\.bam/;
    my $outputFile = $fileBase . "_matureExpressionReplicate.txt";
    my $bam = loadBamFile($bamFile);
    open(OUTPUT,">$outputFile");
    foreach my $chrom (keys %{$products}) {
	foreach my $product (@{$products->{$chrom}}) {
	    my($id,$location,$ctlCount,$count,$fiveHet) = @{$product};
	    my($chr,$start,$stop,$strand) = parseLocation($location);
	    my($expCount,$expFiveHet,$expThreeHet) = getBamProductCount($bam,$chrom,$start,$stop,$strand,$parameters);
	    print OUTPUT "$id\t$location\t$expCount\t$expFiveHet\t$expThreeHet\n";
	}
    }
    close(OUTPUT);
}

sub extractAllProducts {
    my($bamFile,$chromLengths,$parameters) = @_;    
    my $bam = loadBamFile($bamFile);
    my $minDist = $parameters->{minDist} or die "minDist: parameter not loaded.";
    my $minCount = $parameters->{minCount} or die "minCount: parameter not loaded.";
    my $minHet = $parameters->{minHet} or die "minHet: parameter not loaded.";
    my @strands = ("+","-");
    my @chroms = sort keys(%{$chromLengths});
    my($bamBase) = $bamFile =~ /(.*).bam/;
    my $productsFile = $bamBase."_products.txt";
    open(PRODUCTS,">$productsFile");
    print PRODUCTS "#ID\tLocation\t$bamBase\t5'Het\n";
    my $longProductsFile = $bamBase."_longProducts.txt";
    open(LONGPRODUCTS,">$longProductsFile");
    print LONGPRODUCTS "#ID\tLocation\t$bamBase\t5'Het\n";
    my $windowLength = 1000000;
    my $COUNT = 1;
    foreach my $chrom (@chroms) {
	foreach my $strand (@strands) {
	    for(my $start=0;$start<=$chromLengths->{$chrom};$start+=$windowLength) {
		# for the chrom/strand identify all distinct start/stop positions.
		my %distinctReadCount;
		my @reads = $bam->get_features_by_location(-seq_id => $chrom,
							   -start => $start,
							   -end => $start+$windowLength);
		foreach my $read (@reads) {
		    my $rStart = $read->start;
		    my $rStop = $read->end;
		    my $rStrand = mapBamStrand($read->strand);
		    if($strand eq $rStrand) {
			$distinctReadCount{$rStart."x".$rStop}++;
		    }
		}
		# for each start/stop position, sort by count.
		my @products;
		foreach my $key (sort {$distinctReadCount{$b} <=> $distinctReadCount{$a}} keys %distinctReadCount) {
		    my($drStart,$drStop) = split(/x/,$key);
		    unless(overlapsProducts(\@products,$drStart,$drStop,$strand,$parameters)) {
			push(@products,[$drStart,$drStop,$distinctReadCount{$key}]);
		    }
		}
		foreach my $product (@products) {
		    my($pStart,$pStop) = @{$product};		
		    my($expCount,$expFiveHet,$expThreeHet) = getBamProductCount($bam,$chrom,$pStart,$pStop,$strand,$parameters);
		    my $id = "srp".$COUNT;		    
		    if($expCount >= $minCount) {
			if($pStop - $pStart + 1 <= 27) {
			    # products with a low heterogeneity.
			    print PRODUCTS "$id\t$chrom:$pStart..$pStop:$strand\t$expCount\t$expFiveHet\t$expThreeHet\n";
			} else {
			    # products with a low heterogeneity.
			    print LONGPRODUCTS "$id\t$chrom:$pStart..$pStop:$strand\t$expCount\t$expFiveHet\t$expThreeHet\n";
			}
			$COUNT++;
		    }
		}
	    }
	}
    }
    close(PRODUCTS);
    close(LONGPRODUCTS);
    return $productsFile;
}

sub extractAllProductsForChrom {
    my($bamFile,$chromLength,$chrom,$parameters) = @_;    
    my $bam = loadBamFile($bamFile);
    my $minDist = $parameters->{minDist} or die "minDist: parameter not loaded.";
    my $minCount = $parameters->{minCount} or die "minCount: parameter not loaded.";
    my $minHet = $parameters->{minHet} or die "minHet: parameter not loaded.";
    my @strands = ("+","-");
    my($bamBase) = $bamFile =~ /(.*).bam/;
    my $allProductsFile = $bamBase."_allProducts.".$chrom.".txt";
    open(ALLPRODUCTS,">$allProductsFile"); 
    print ALLPRODUCTS "#ID\tLocation\t$bamBase\t5'Het\n";
    my $productsFile = $bamBase."_products.".$chrom.".txt";
    open(PRODUCTS,">$productsFile");
    print PRODUCTS "#ID\tLocation\t$bamBase\t5'Het\n";
    my $longProductsFile = $bamBase."_longProducts.".$chrom.".txt";
    open(LONGPRODUCTS,">$longProductsFile");
    print LONGPRODUCTS "#ID\tLocation\t$bamBase\t5'Het\n";
    my $windowLength = 1000000;
    my $COUNT = 1;
    foreach my $strand (@strands) {
	for(my $start=0;$start<=$chromLength;$start+=$windowLength) {
	    # for the chrom/strand identify all distinct start/stop positions.
	    my %distinctReadCount;
	    my @reads = $bam->get_features_by_location(-seq_id => $chrom,
						       -start => $start,
						       -end => $start+$windowLength);
	    foreach my $read (@reads) {
		my $rStart = $read->start;
		my $rStop = $read->end;
		my $rStrand = mapBamStrand($read->strand);
		if($strand eq $rStrand) {
		    $distinctReadCount{$rStart."x".$rStop}++;
		}
	    }
	    # for each start/stop position, sort by count.
	    my @products;
	    foreach my $key (sort {$distinctReadCount{$b} <=> $distinctReadCount{$a}} keys %distinctReadCount) {
		my($drStart,$drStop) = split(/x/,$key);
		unless(overlapsProducts(\@products,$drStart,$drStop,$strand,$parameters)) {
		    push(@products,[$drStart,$drStop,$distinctReadCount{$key}]);
		}
	    }
	    foreach my $product (@products) {
		my($pStart,$pStop) = @{$product};		
		my($expCount,$expFiveHet,$expThreeHet) = getBamProductCount($bam,$chrom,$pStart,$pStop,$strand,$parameters);
		my $id = "srp".$COUNT;		    
		if($expCount >= $minCount) {
		    # the full product list
		    print ALLPRODUCTS "$id\t$chrom:$pStart..$pStop:$strand\t$expCount\t$expFiveHet\t$expThreeHet\n";
		    if($expFiveHet <= $minHet) {
			if($pStop - $pStart + 1 <= 27) {
			    # products with a low heterogeneity.
			    print PRODUCTS "$id\t$chrom:$pStart..$pStop:$strand\t$expCount\t$expFiveHet\t$expThreeHet\n";
			} else {
			    # products with a low heterogeneity.
			    print LONGPRODUCTS "$id\t$chrom:$pStart..$pStop:$strand\t$expCount\t$expFiveHet\t$expThreeHet\n";
			}
		    }
		    $COUNT++;
		}
	    }
	}
    }
    close(ALLPRODUCTS);
    close(PRODUCTS);
    close(LONGPRODUCTS);
    return $productsFile;
}

sub extractAllProductsExpCtl {
    my($bamExpCtl,$chromLengths,$parameters) = @_;    
    my $minDist = $parameters->{minDist} or die "minDist: parameter not loaded.";
    my $minCount = $parameters->{minCount} or die "minCount: parameter not loaded.";
    my $minHet = $parameters->{minHet} or die "minHet: parameter not loaded.";
    my($expLabel,$expBam,$expTotal,$expFile) = @{$bamExpCtl->{"exp"}};
    my($ctlLabel,$ctlBam,$ctlTotal,$ctlFile) = @{$bamExpCtl->{"ctl"}};
    my @strands = ("+","-");
    my @chroms = sort keys(%{$chromLengths});
    my $productsFile = $expLabel."_vs_".$ctlLabel."_products.txt";
    open(PRODUCTS,">$productsFile");
    print PRODUCTS "#ID\tLocation\t$ctlLabel\t$expLabel\t5'Het\n";
    my $longProductsFile = $expLabel."_vs_".$ctlLabel."_longProducts.txt";
    open(LONGPRODUCTS,">$longProductsFile");
    print LONGPRODUCTS "#ID\tLocation\t$ctlLabel\t$expLabel\t5'Het\n";
    my $initWindowLength = 1000;
    my $COUNT = 1;
    foreach my $chrom (@chroms) {
	foreach my $strand (@strands) {
	    # for the chrom/strand identify all distinct start/stop positions.
	    # first loop through all start positions.
	    my $startPos = 0;
	    while($startPos <= $chromLengths->{$chrom}) {
		# push out the windowLength until no reads extend beyond the ends.
		my $windowLength = $initWindowLength;
		# check if any reads overlap the endpoint.
		my @expReads = $expBam->get_features_by_location(-seq_id => $chrom, -start => $startPos+$windowLength, -end => $startPos+$windowLength);
		my @ctlReads = $ctlBam->get_features_by_location(-seq_id => $chrom, -start => $startPos+$windowLength, -end => $startPos+$windowLength);
		while(@expReads|@ctlReads) {
		    @expReads = $expBam->get_features_by_location(-seq_id => $chrom, -start => $startPos+$windowLength, -end => $startPos+$windowLength);
		    @ctlReads = $ctlBam->get_features_by_location(-seq_id => $chrom, -start => $startPos+$windowLength, -end => $startPos+$windowLength);
		    $windowLength++;
		}
		@expReads = ();
		@ctlReads = ();
		# now we have a good endpoint, proceed to extract the read loci.
		print "$chrom $strand $startPos ", $startPos + $windowLength, "\n";
		@expReads = $expBam->get_features_by_location(-seq_id => $chrom,
								 -start => $startPos,
								 -end => $startPos + $windowLength);
		my %distinctReadCount;
		foreach my $read (@expReads) {
		    my $rStart = $read->start;
		    my $rStop = $read->end;
		    my $rStrand = mapBamStrand($read->strand);
		    if($strand eq $rStrand) {
			$distinctReadCount{$rStart."x".$rStop}++;
		    }
		}
		@expReads = ();
		@ctlReads = $ctlBam->get_features_by_location(-seq_id => $chrom,
								 -start => $startPos,
								 -end => $startPos + $windowLength);
		foreach my $read (@ctlReads) {
		    my $rStart = $read->start;
		    my $rStop = $read->end;
		    my $rStrand = mapBamStrand($read->strand);
		    if($strand eq $rStrand) {
			$distinctReadCount{$rStart."x".$rStop}++;
		    }
		}
		# for each start/stop position, sort by count.
		my @products;
		foreach my $key (sort {$distinctReadCount{$b} <=> $distinctReadCount{$a}} keys %distinctReadCount) {
		    my($drStart,$drStop) = split(/x/,$key);
		    unless(overlapsProducts(\@products,$drStart,$drStop,$strand,$parameters)) {
			push(@products,[$drStart,$drStop,$distinctReadCount{$key}]);
		    }
		}
		%distinctReadCount = ();
		foreach my $product (@products) {
		    my($pStart,$pStop) = @{$product};
		    my($expCount,$expFiveHet,$expThreeHet) = getBamProductCount($expBam,$chrom,$pStart,$pStop,$strand,$parameters);
		    my($ctlCount,$ctlFiveHet,$ctlThreeHet) = getBamProductCount($ctlBam,$chrom,$pStart,$pStop,$strand,$parameters);
		    my $id = "srp".$COUNT;
		    if($expCount >= $minCount) {
			if($pStop - $pStart + 1 <= 27) {
			    # products with a low heterogeneity.
			    print PRODUCTS "$id\t$chrom:$pStart..$pStop:$strand\t$ctlCount\t$expCount\t$expFiveHet\t$expThreeHet\n";
			} else {
			    # print the longer smallRNAs
			    print LONGPRODUCTS "$id\t$chrom:$pStart..$pStop:$strand\t$ctlCount\t$expCount\t$expFiveHet\t$expThreeHet\n";
			}
			$COUNT++;
		    }
		}
		@products = ();
		$startPos += $windowLength;
	    }
	}
    }
    close(ALLPRODUCTS);
    close(PRODUCTS);
    close(LONGPRODUCTS);
    return $productsFile;
}

sub extractAllProductsExpCtlForChrom {
    my($bamExpCtl,$chromLengths,$chrom,$parameters) = @_;    
    my $minDist = $parameters->{minDist} or die "minDist: parameter not loaded.";
    my $minCount = $parameters->{minCount} or die "minCount: parameter not loaded.";
    my $minHet = $parameters->{minHet} or die "minHet: parameter not loaded.";
    my($expLabel,$expBam,$expTotal,$expFile) = @{$bamExpCtl->{"exp"}};
    my($ctlLabel,$ctlBam,$ctlTotal,$ctlFile) = @{$bamExpCtl->{"ctl"}};
    my @strands = ("+","-");
    my $productsFile = $expLabel."_vs_".$ctlLabel."_".$chrom."_products.txt";
    open(PRODUCTS,">$productsFile");
    print PRODUCTS "#ID\tLocation\t$ctlLabel\t$expLabel\t5'Het\n";
    my $longProductsFile = $expLabel."_vs_".$ctlLabel."_".$chrom."_longProducts.txt";
    open(LONGPRODUCTS,">$longProductsFile");
    print LONGPRODUCTS "#ID\tLocation\t$ctlLabel\t$expLabel\t5'Het\n";
    my $initWindowLength = 1000;
    my $COUNT = 1;
    foreach my $strand (@strands) {
	# for the chrom/strand identify all distinct start/stop positions.
	# first loop through all start positions.
	my $startPos = 0;
	while($startPos <= $chromLengths->{$chrom}) {
	    # push out the windowLength until no reads extend beyond the ends.
	    my $windowLength = $initWindowLength;
	    # check if any reads overlap the endpoint.
	    my @expReads = $expBam->get_features_by_location(-seq_id => $chrom, -start => $startPos+$windowLength, -end => $startPos+$windowLength);
	    my @ctlReads = $ctlBam->get_features_by_location(-seq_id => $chrom, -start => $startPos+$windowLength, -end => $startPos+$windowLength);
	    while(@expReads|@ctlReads) {
		@expReads = $expBam->get_features_by_location(-seq_id => $chrom, -start => $startPos+$windowLength, -end => $startPos+$windowLength);
		@ctlReads = $ctlBam->get_features_by_location(-seq_id => $chrom, -start => $startPos+$windowLength, -end => $startPos+$windowLength);
		$windowLength++;
	    }
	    @expReads = ();
	    @ctlReads = ();
	    # now we have a good endpoint, proceed to extract the read loci.
	    print "$chrom $strand $startPos ", $startPos + $windowLength, "\n";
	    @expReads = $expBam->get_features_by_location(-seq_id => $chrom,
							  -start => $startPos,
							  -end => $startPos + $windowLength);
	    my %distinctReadCount;
	    foreach my $read (@expReads) {
		my $rStart = $read->start;
		my $rStop = $read->end;
		my $rStrand = mapBamStrand($read->strand);
		if($strand eq $rStrand) {
		    $distinctReadCount{$rStart."x".$rStop}++;
		}
		}
	    @expReads = ();
	    @ctlReads = $ctlBam->get_features_by_location(-seq_id => $chrom,
							  -start => $startPos,
							  -end => $startPos + $windowLength);
	    foreach my $read (@ctlReads) {
		my $rStart = $read->start;
		my $rStop = $read->end;
		my $rStrand = mapBamStrand($read->strand);
		if($strand eq $rStrand) {
		    $distinctReadCount{$rStart."x".$rStop}++;
		}
	    }
	    # for each start/stop position, sort by count.
	    my @products;
	    foreach my $key (sort {$distinctReadCount{$b} <=> $distinctReadCount{$a}} keys %distinctReadCount) {
		my($drStart,$drStop) = split(/x/,$key);
		unless(overlapsProducts(\@products,$drStart,$drStop,$strand,$parameters)) {
		    push(@products,[$drStart,$drStop,$distinctReadCount{$key}]);
		}
	    }
	    %distinctReadCount = ();
	    foreach my $product (@products) {
		my($pStart,$pStop) = @{$product};
		my($expCount,$expFiveHet,$expThreeHet) = getBamProductCount($expBam,$chrom,$pStart,$pStop,$strand,$parameters);
		my($ctlCount,$ctlFiveHet,$ctlThreeHet) = getBamProductCount($ctlBam,$chrom,$pStart,$pStop,$strand,$parameters);
		my $id = $chrom.".".$COUNT;
		if($expCount >= $minCount) {
		    if($pStop - $pStart + 1 <= 27) {
			# products with a low heterogeneity.
			print PRODUCTS "$id\t$chrom:$pStart..$pStop:$strand\t$ctlCount\t$expCount\t$expFiveHet\t$expThreeHet\n";
		    } else {
			# print the longer smallRNAs
			print LONGPRODUCTS "$id\t$chrom:$pStart..$pStop:$strand\t$ctlCount\t$expCount\t$expFiveHet\t$expThreeHet\n";
		    }
		    $COUNT++;
		}
	    }
	    @products = ();
	    $startPos += $windowLength;
	}
    }
    close(ALLPRODUCTS);
    close(PRODUCTS);
    close(LONGPRODUCTS);
    return $productsFile;
}

sub extractSignificantProducts {
    my($productsFile,$A,$B,$parameters) = @_;
    # not the products could be renamed at this point, hence the regexp below.
    my($fileBase) = $productsFile =~ /([^\/]*)_(products|allProducts)/ or die "Could not parse the productsFile: $productsFile\n";
    my($expLabel,$ctlLabel) = $fileBase =~ /(.*)_vs_(.*)/ or die "Could not parse the file base: $fileBase\n";
    my $FDR = $parameters->{FDR} or die "FDR: parameter not loaded.\n";
    open(PF,$productsFile) or die "could not open $productsFile\n";
    my %products;
    while(<PF>) {
	chomp;
	unless(/^\#/) {
	    my($id,$location,$ctlCount,$expCount,$expFiveHet,$expThreeHet) = split(/\t/);
	    $products{$id}{$location} = [$expCount,$ctlCount,$expFiveHet,$expThreeHet];
	}
    }
    my @allLoci;
    foreach my $id (keys %products) {
	my $expSum = 0;
	my $ctlSum = 0;
	my $tot = 0;
	my $firstLocation = "";
	my $firstExpFiveHet;
	my $firstExpThreeHet;
	# must loop through all distinct locations, because here duplicated products (by sequence) have the same id, but different genomic locations
	foreach my $location (keys %{$products{$id}}) {
	    unless($products{$id}{$location}) {die "no information for $id $location\n";}
	    my($expCount,$ctlCount,$expFiveHet,$expThreeHet) = @{$products{$id}{$location}};
	    unless($firstLocation) {
		$firstLocation = $location;
		$firstExpFiveHet = $expFiveHet;
		$firstExpThreeHet = $expThreeHet;
	    }
	    $expSum += $expCount;
	    $ctlSum += $ctlCount;
	    $tot++;	    
	}
	my $expAvg = $expSum / $tot;
	my $ctlAvg = $ctlSum / $tot;
	$ctlAvg = $A + $B*$ctlAvg;
	$ctlAvg = $ctlAvg > 0 ? $ctlAvg : 1;
	$ctlAvg = sprintf("%.1f",$ctlAvg);
 	my $pValue = binomialPValue($expAvg,$ctlAvg);
	push(@allLoci,[$pValue,$id,$firstLocation,$expAvg,$ctlAvg,$firstExpFiveHet,$firstExpThreeHet]);
    }
    @allLoci = sort {$a->[0] <=> $b->[0]} @allLoci;
    my $total = scalar(@allLoci);
    my $rank = 0;
    my $sigProductsFile = $fileBase."_significantProducts.txt";
    open(SIGPRODUCTS,">$sigProductsFile") or die "could not open $sigProductsFile for writing.\n";
    print SIGPRODUCTS "#ID\tLocation\t$ctlLabel\t$expLabel\tLogFold\tPValue\tQValue\tRank\t5'Het\t3'Het\n";
    my $allProductsFile = $fileBase . "_allProductSignificance.txt";
    open(ALLPRODUCTS,">$allProductsFile") or die "Could not open $allProductsFile for writing.\n";
    print ALLPRODUCTS "#ID\tLocation\t$ctlLabel\t$expLabel\tLogFold\tPValue\tQValue\tRank\t5'Het\t3'Het\n";
    foreach my $locus (@allLoci) {
	$rank++;
	my($pValue,$id,$location,$expCount,$ctlCount,$expFiveHet,$expThreeHet) = @{$locus};
	my $qValue = $pValue*$total/$rank;
	my $logFold = log(($expCount+1)/($ctlCount+1));
	if($qValue <= $FDR) {
	    print SIGPRODUCTS "$id\t$location\t$ctlCount\t$expCount\t$logFold\t$pValue\t$qValue\t$rank\t$expFiveHet\t$expThreeHet\n";
	}
	print ALLPRODUCTS "$id\t$location\t$ctlCount\t$expCount\t$logFold\t$pValue\t$qValue\t$rank\t$expFiveHet\t$expThreeHet\n";
    }
    close(SIGPRODUCTS);
    close(ALLPRODUCTS);
}

sub processSignificantProducts {
    my($products,$genomeDir,$parameters) = @_;
    my $minDist = $parameters->{minDist} or die "minDist: parameter not loaded.";
    my $ctlLabel = $parameters->{ctlLabel} or die "ctlLabel: parameter not loaded.";
    my $expLabel = $parameters->{expLabel} or die "expLabel: parameter not loaded.";
    # productInfo is the core list of products and their sequences.
    my $source = $expLabel."_vs_".$ctlLabel;
    my $type = "match";
    my $label = $parameters->{label} ? $parameters->{label} : "product";
    my $productInfo = $source."_".$label."_Expression.txt";
    open(INFO,">$productInfo") or die "Count not open $productInfo for writing.\n";
    print INFO "#ID\tLocation\t5'Het\t3'Het\tSequence\t$ctlLabel\t$expLabel\n";
    my $productGff = $source."_products.gff";
    open(GFF,">$productGff") or die "Could not open $productGff for writing.\n";
    foreach my $chrom (keys %{$products}) {
	# the file for the chrom
	my $chromFile = $genomeDir . "/" . $chrom . ".fa";
        # actually just loading the chrom sequence here.
        my $genome = loadGenome($chromFile);
	foreach my $products (@{$products->{$chrom}}) {
	    my($id,$location,$ctlCount,$expCount,$pValue,$qValue,$rank,$expFiveHet,$expThreeHet) = @{$products};
	    # by construction, this should be the sequence of the most abundant product.
	    my $sequence = getSequence($location,$genome);
	    my($chrom,$start,$stop,$strand) = parseLocation($location);
	    $ctlCount = sprintf("%.1f",$ctlCount);
	    print INFO "$id\t$location\t$expFiveHet\t$expThreeHet\t$sequence\t$ctlCount\t$expCount\n";
	    print GFF "$chrom\t$source\t$type\t$start\t$stop\t$qValue\t$strand\t.\tID=$id\n";
	}
    }
    close(INFO);
    close(GFF);
}

sub overlapsProducts {
    my($products,$drStart,$drStop,$strand,$parameters)=@_;
    my $minDist = $parameters->{minDist} or die "minDist: parameter not loaded.";
    foreach my $product (@{$products}) {
	my($start,$stop)=@{$product};
	if($strand eq "+") {
	    # use start for 5' end.
	    if(abs($start-$drStart) <= $minDist) {
		return 1;
	    }
	} else {
	    # use stop for 5' end when on "-" strand.
	    if(abs($stop-$drStop) <= $minDist) {
		return 1;
	    }
	}
    }
    return 0;
}

sub getBamProductCount {
    my($bam,$chrom,$start,$stop,$strand,$parameters) = @_;
    my $minDist = $parameters->{minDist} or die "minDist: parameter not loaded.";
    my @alignments = $bam->get_features_by_location(-seq_id => $chrom,
						    -start  => $start - $minDist,
						    -end    => $stop + 2*$minDist);
    my $count = 0;
    my %fivePosHash;
    my %threePosHash;
    foreach my $read (@alignments) {
	my $rStart = $read->start;
	my $rStop = $read->end;
	# 5' end methdod.
	my $rStrand = mapBamStrand($read->strand);
	if($strand eq $rStrand) {
	    if($strand eq "+") {
		# for forward strand, 5' end is "start"
		if(($start - $minDist <= $rStart)&&($rStart <= $start + $minDist)) {
		    $count++;
		    $fivePosHash{$rStart}++;
		    $threePosHash{$rStop}++;
		}
	    } else {
		# for reverse strand, 5' end is "stop"
		if(($stop - $minDist <= $rStop)&&($rStop <= $stop + $minDist)) {
		    $count++;
		    $fivePosHash{$rStop}++;
		    $threePosHash{$rStart}++;
		}
	    }
	}
    }
    if($count) {
	# sort the positions, highest to lowest
	# compute the 5' heterogeneity
	my @fivePos = sort {$fivePosHash{$b} <=> $fivePosHash{$a}} keys(%fivePosHash);
	my $maxFivePos = shift(@fivePos);
	my $maxFivePosCount = $fivePosHash{$maxFivePos};
	my $fiveHet = sprintf("%.3f",1 - $maxFivePosCount/$count);
	# compute the 3' heterogeneity
	my @threePos = sort {$threePosHash{$b} <=> $threePosHash{$a}} keys(%threePosHash);
	my $maxThreePos = shift(@threePos);
	my $maxThreePosCount = $threePosHash{$maxThreePos};
	my $threeHet = sprintf("%.3f",1 - $maxThreePosCount/$count);	
	return($count,$fiveHet,$threeHet);
    } else {
	return (0,1,1);
    }
}

sub getBamOverlap {
    my($locusFile,$bamFile,$parameters)=@_;   
    my $buffer = $parameters->{buffer} or die "buffer: parameter not loaded.\n";
    my $loci = readLocusGff($locusFile);
    my($locusBase) = $locusFile =~ /\/?([^\/]*).gff/;
    my($bamBase) = $bamFile =~ /\/?([^\/]*).bam/;
    my $outFile = $locusBase . "_vs_" . $bamBase . ".overlap"; 
    my $bam = Bio::DB::Sam->new( -bam  => $bamFile );
    open(OUTFILE,">".$outFile) or die "could not open $outFile for writing.\n";
    foreach my $chrom (keys %{$loci}) {
	foreach my $locus (@{$loci->{$chrom}}) {
	    my($start,$stop,$strand,$id) = @{$locus};
	    # bam is 1-based
	    # gff is 1-based ( as is most input files here ) 
            # assuming the coordinates are 1-based as in gff
	    my @alignments = $bam->get_features_by_location(-seq_id => $chrom,
							    -start  => $start - $buffer,
							    -end    => $stop + $buffer);
	}
    }
    close(OUTFILE);
}

sub getBamLocusCounts {
    my($location,$bam,$paramters) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my @alignments = $bam->get_features_by_location(-seq_id => $chrom,
						    -start  => $start,
						    -end    => $stop);
    my $count = 0;
    foreach my $read (@alignments) {
	my $rStart = $read->start;
	my $rStop = $read->end;
	if(($start <= $rStart)&&($rStop <= $stop)) {
	    my $rStrand = mapBamStrand($read->strand);	    
	    if($strand eq $rStrand) {
		$count++;
	    }
	}
    }
    return $count;
}

sub processUnknownRegions {
    my($folds,$bamList,$parameters)=@_;   
    my $minLength = $parameters->{minLength};
    open(HL,">locusList.txt");
    open(HDR,">locusDistinctReads.txt");
    open(HPL,">locusProducts.txt");
    print HL "#tag\tchrom\tstart\tstop\tstrand\thStart\thStop\t";
    print HL "totalReads\tmfe\taapd\tsequence\tfold";
    print HDR "#tag\tstart\tstop\tstrand\ttotalReads\toffsetStart\trelativeStart"; 
    print HPL "#tag\tside\ttype\ttotalReads\ttotalMostAbundant\tstart\tstop\t5'Het\t3'Het\tsequence";
    foreach my $entry (@{$bamList}) {
	my($label,$bam,$total) = @{$entry};
	print HL "\t$label";
	print HPL "\t$label";
	print HDR "\t$label";
    }
    print HL "\n";
    print HDR "\n";
    print HPL "\n";
    foreach my $chrom (keys %{$folds}) {
        foreach my $strand (keys %{$folds->{$chrom}}) {	    
	    foreach my $foldInfo (@{$folds->{$chrom}{$strand}}) {
		my($tag,$start,$stop,$leftCenter,$rightCenter,$mfe,$sequence,$fold) = @{$foldInfo};				
		my $location = "$chrom:$start..$stop:$strand";		
		my $basePairs = getBasePairs($leftCenter,$rightCenter,
					     $fold,$parameters);
		my($hStart,$hStop) = getFullHairpinWindow($fold,$leftCenter,$rightCenter);
		my $distinctReads=getDistinctReadInfo($chrom,$start,$stop,$strand,$bamList);
		my $products = getUnknownProductList($leftCenter,$rightCenter,
						     $fold,$basePairs,$location,
						     $distinctReads,$parameters);
		my $revStrand = revStrand($strand);
		my $revDistinctReads=getDistinctReadInfo($chrom,$start,$stop,$revStrand,$bamList);
		my $revProducts=getUnknownProductList($leftCenter,$rightCenter,
						      $fold,$basePairs,$location,
						      $revDistinctReads,$parameters);
		my($tpd,$totalRP)=getReverseProductDisplacement($products,$revProducts,$parameters);
		my $aapd = $totalRP ? $tpd/$totalRP : 0.0;
		my $COUNT=1;
		my $total = 0;
		my %mirLibraryCounts;
		foreach my $product (@{$products}) {
		    # process the product information and print...
		    my($side,$fakeType,$prodList,$prodCount,$maxProdCount,$relStart,$relStop,$offset,$gStart,$productLibraryCounts) = @{$product};
		    $total += $prodCount;
		    my $type = $COUNT;
		    my $length = $relStop - $relStart + 1;
		    my $productSequence = substr($sequence,$relStart,$length);
		    my $fivePrimeHet = computeFivePrimeHet($prodList);
		    my $threePrimeHet = computeThreePrimeHet($prodList);
		    print HPL "${tag}\t$side\t$type\t$prodCount\t";
		    print HPL "$maxProdCount\t$relStart\t$relStop\t";
		    print HPL "$fivePrimeHet\t$threePrimeHet\t";
		    print HPL "$productSequence";
		    foreach my $entry (@{$bamList}) {
			my($label,$bam,$total) = @{$entry};
			print HPL "\t", $productLibraryCounts->{$label};
			$mirLibraryCounts{$label} += $productLibraryCounts->{$label};
		    }
		    print HPL "\n";
		    $COUNT++;		    
		    # print the distinct reads...
		    foreach my $read (@{$prodList}) {
			my($relStart,$relStop,$offset,$gStart,
			   $count,$libraryCounts) = @{$read};
			my $gStop = $gStart + ($relStop - $relStart);
			print HDR "$tag\t$gStart\t$gStop\t$strand\t$count\t";
			print HDR "$offset\t$relStart";
			foreach my $entry (@{$bamList}) {
			    my($label,$bam) = @{$entry};
			    printf(HDR "\t%.3f",$libraryCounts->{$label});
			}
			print HDR "\n";
		    }
		}
		foreach my $product (@{$revProducts}) {
		    # process the product information and print...
		    my($side,$fakeType,$prodList,$prodCount,$maxProdCount,$relStart,$relStop,$offset,$gStart,$productLibraryCounts) = @{$product};
		    $total += $prodCount;
		    my $type = $COUNT . "-as";
		    my $length = $relStop - $relStart + 1;
		    my $productSequence = substr($sequence,$relStart,$length);
		    my $fivePrimeHet = computeFivePrimeHet($prodList);
		    my $threePrimeHet = computeThreePrimeHet($prodList);
		    print HPL "${tag}\t$side\t$type\t$prodCount\t";
		    print HPL "$maxProdCount\t$relStart\t$relStop\t";
		    print HPL "$fivePrimeHet\t$threePrimeHet\t";
		    print HPL "$productSequence";
		    foreach my $entry (@{$bamList}) {
			my($label,$bam,$total) = @{$entry};
			print HPL "\t", $productLibraryCounts->{$label};
			$mirLibraryCounts{$label} += $productLibraryCounts->{$label};
		    }
		    print HPL "\n";
		    $COUNT++;		    
		    # print the distinct reads...
		    foreach my $read (@{$prodList}) {
			my($relStart,$relStop,$offset,$gStart,
			   $count,$libraryCounts) = @{$read};
			my $gStop = $gStart + ($relStop - $relStart);
			print HDR "$tag\t$gStart\t$gStop\t$revStrand\t$count\t";
			print HDR "$offset\t$relStart";
			foreach my $entry (@{$bamList}) {
			    my($label,$bam) = @{$entry};
			    printf(HDR "\t%.3f",$libraryCounts->{$label});
			}
			print HDR "\n";
		    }
		}
		print HL "$tag\t$chrom\t$start\t$stop\t$strand\t$hStart\t$hStop\t";
		printf(HL "%.3f\t%.3f\t%.3f\t",$total,$mfe,$aapd);
		print HL "$sequence\t$fold";
		foreach my $entry (@{$bamList}) {
		    my($label,$bam) = @{$entry};
		    if($mirLibraryCounts{$label}) {
			printf(HL "\t%.3f",$mirLibraryCounts{$label});
		    } else {
			printf(HL "\t0");
		    }
		}
		print HL "\n";
	    }
	}
    }
}

sub processMirRegions {
    my($folds,$bamList,$parameters)=@_;   
    my $minLength = $parameters->{minLength};
    open(HL,">miRList.txt");
    open(HDR,">miRDistinctReads.txt");
    open(HPL,">miRProductList.txt");
    print HL "#tag\tchrom\tstart\tstop\tstrand\thStart\thStop\t";
    print HL "totalReads\tmfe\taapd\tsequence\tfold";
    print HDR "#tag\tstart\tstop\tstrand\ttotalReads\toffsetStart\trelative Start";
    print HPL "#tag\tside\ttype\ttotalReads\ttotalMostAbundant\tstart\tstop\t5'Het\t3'Het\tsequence";
    foreach my $entry (@{$bamList}) {
	my($label,$bam) = @{$entry};
	print HL "\t$label";
	print HPL "\t$label";
 	print HDR "\t$label";
    }
    print HL "\n";
    print HDR "\n";
    print HPL "\n";
    foreach my $chrom (keys %{$folds}) {
        foreach my $strand (keys %{$folds->{$chrom}}) {	    
	    foreach my $foldInfo (@{$folds->{$chrom}{$strand}}) {
		my($tag,$start,$stop,$leftCenter,$rightCenter,$mfe,$sequence,$fold) = @{$foldInfo};
		#print "checking $tag ";
		my $location = "$chrom:$start..$stop:$strand";		
		my $basePairs = getBasePairs($leftCenter,$rightCenter,
					     $fold,$parameters);
		my $hairpinLength = getMaxHairpinLength($fold,$leftCenter,$rightCenter);
		my($hStart,$hStop) = getFullHairpinWindow($fold,$leftCenter,$rightCenter);
		my $distinctReads=getDistinctReadInfo($chrom,$start,$stop,$strand,$bamList);
		my $productInfo=extractProducts($leftCenter,$rightCenter,
						$fold,$basePairs,$location,
						$distinctReads,$parameters);
		my $revStrand = revStrand($strand);
		my $revDistinctReads=getDistinctReadInfo($chrom,$start,$stop,$revStrand,$bamList);
		my $revProductInfo=extractProducts($leftCenter,$rightCenter,
						   $fold,$basePairs,$location,
						   $revDistinctReads,$parameters);
		my($tpd,$totalRP)=getReverseProductDisplacement($productInfo,$revProductInfo,$parameters);
		my $aapd = $totalRP ? $tpd/$totalRP : 0.0;		
		my $total = 0;
		my %mirLibraryCounts;
		foreach my $product (@{$productInfo}) {
		    my($side,$newType,$prodList,$prodCount,$maxProdCount,
		       $relStart,$relStop,$offset,$gStart) = @{$product};
		    my $length = $relStop - $relStart + 1;
		    my $productSequence = substr($sequence,$relStart,$length);
		    my $fivePrimeHet = computeFivePrimeHet($prodList);
		    my $threePrimeHet = computeThreePrimeHet($prodList);
		    unless($newType eq "out") {
			$total += $prodCount;
		    }
		    my $totalLibraryCounts = getTotalLibraryCounts($prodList);
		    print HPL "$tag\t$side\t$newType\t$prodCount\t";
		    print HPL "$maxProdCount\t$relStart\t$relStop\t";
		    print HPL "$fivePrimeHet\t$threePrimeHet\t";
		    print HPL "$productSequence";				
		    foreach my $entry (@{$bamList}) {
			my($label,$bam) = @{$entry};
			printf(HPL "\t%.3f",$totalLibraryCounts->{$label});
			$mirLibraryCounts{$label} += $totalLibraryCounts->{$label};
		    }
		    print HPL "\n";
		    foreach my $read (@{$prodList}) {
			my($relStart,$relStop,$offset,$gStart,
			   $count,$libraryCounts) = @{$read};
			my $gStop = $gStart + ($relStop - $relStart);
			print HDR "$tag\t$gStart\t$gStop\t$strand\t$count\t";
			print HDR "$offset\t$relStart";
			foreach my $entry (@{$bamList}) {
			    my($label,$bam) = @{$entry};			    
			    printf(HDR "\t%.3f",$libraryCounts->{$label});
			}
			print HDR "\n";
		    }
		}
		foreach my $product (@{$revProductInfo}) {
		    my($side,$newType,$prodList,$prodCount,$maxProdCount,
		       $relStart,$relStop,$offset,$gStart) = @{$product};
		    $newType .= "roR";
		    my $length = $relStop - $relStart + 1;
		    my $productSequence = substr($sequence,$relStart,$length);
		    my $totalLibraryCounts = getTotalLibraryCounts($prodList);
		    print HPL "$tag\t$side\t$newType\t$prodCount\t";
		    print HPL "$maxProdCount\t$relStart\t$relStop\t";
		    print HPL "$productSequence";				
		    foreach my $entry (@{$bamList}) {
			my($label,$bam) = @{$entry};
			printf(HPL "\t%.3f",$totalLibraryCounts->{$label});
			$mirLibraryCounts{$label} += $totalLibraryCounts->{$label};
		    }
		    print HPL "\n";
		    foreach my $read (@{$prodList}) {
			my($relStart,$relStop,$offset,$gStart,
			   $count,$libraryCounts) = @{$read};
			my $gStop = $gStart + ($relStop - $relStart);
			print HDR "$tag\t$gStart\t$gStop\t$revStrand\t$count\t";
			print HDR "$offset\t$relStart";
			foreach my $entry (@{$bamList}) {
			    my($label,$bam) = @{$entry};
			    printf(HDR "\t%.3f",$libraryCounts->{$label});
			}
			print HDR "\n";
		    }
		}
		print HL "$tag\t$chrom\t$start\t$stop\t$strand\t$hStart\t$hStop\t";
		printf(HL "%.3f\t%.3f\t%.3f\t",$total,$mfe,$aapd);
		print HL "$sequence\t$fold";
		foreach my $entry (@{$bamList}) {
		    my($label,$bam) = @{$entry};
		    if($mirLibraryCounts{$label}) {
			printf(HL "\t%.3f",$mirLibraryCounts{$label});
		    } else {
			print HL "\t0";
		    }
		}
		print HL "\n";
	    }
	}
    }		    
}

sub getDistinctReadInfo {
    my($chrom,$start,$stop,$strand,$bamList) = @_;
    my %distinctReadCount;
    my %libraryCounts;
    foreach my $entry (@{$bamList}) {
	my($label,$bam) = @{$entry};
	# bam uses 1-based
	# gff is 1-based.
	# assuming that input start and stop coordinates are 1-based as with gff format.
	my @alignments = $bam->get_features_by_location(-seq_id => $chrom,
							-start  => $start,
							-end    => $stop);
	foreach my $read (@alignments) {
	    my $rStart = $read->start;
	    my $rStop = $read->end;	    
	    my $rStrand = mapBamStrand($read->strand);
	    # make sure it is completely within the locus.
	    # note that bam uses 1-based.
	    if(($start <= $rStart)&&($rStop <= $stop)) {
		if($rStrand eq $strand) {
		    $distinctReadCount{$rStart."x".$rStop}++;
		    $libraryCounts{$rStart."x".$rStop}{$label}++;
		}
	    }
	}
    }
    my @distinctReads;
    foreach my $key (keys %distinctReadCount) {
	my($drStart,$drStop) = split(/x/,$key);
	foreach my $entry (@{$bamList}) {
	    my($label,$bam) = @{$entry};	    
	    unless($libraryCounts{$key}{$label}) {
		$libraryCounts{$key}{$label} = 0;
	    }
	}
	push(@distinctReads,[$drStart,$drStop,$distinctReadCount{$key},$libraryCounts{$key}]);
    }
    # sort distinct reads by abundance.  This is important in getProductInfo
    @distinctReads = sort {$b->[2] <=> $a->[2]} @distinctReads;
    return \@distinctReads;
}

sub getUnknownProductList {
    my($leftCenter,$rightCenter,$fold,$basePairs,$location,$distinctReads,$parameters) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    # first build a product hash:
    my $PRODUCTCOUNT = 1;
    my %productHash;
    foreach my $dRead (@{$distinctReads}) {
	# dStart is genomic position.
	my($dStart,$dStop,$total,$libraryCounts) = @{$dRead};
	my $relStart = $dStart - $start;
	my $relStop = $dStop - $start;
	# relStart is the position relative to the 5' end of readRegion. numbers start from zero.
	if($strand eq "-") {
	    $relStart = $stop - $dStop;
	    $relStop = $stop - $dStart;
	}
	# offset is relative to the leftMiddlePos (last 5' base pair in hairpin)
	my $offset = $relStop - $leftCenter;	
	if(($offset > 0)&&($relStop < $rightCenter)) {
	    $offset = 0;
	} elsif($offset > 0) {
	    $offset = $relStart - $rightCenter;
	}
	my $parsedRead = [$relStart,$relStop,$offset,$dStart,
			  $total,$libraryCounts];
	if(my $id = overlapCurrent(\%productHash,$parsedRead,$parameters)) {
	    push(@{$productHash{$id}},$parsedRead);	    
	} else {
	    push(@{$productHash{$PRODUCTCOUNT}},$parsedRead);
	    $PRODUCTCOUNT++;
	}
    }
    my @productList;
    # for each distinct product...
    foreach my $id (keys %productHash) {
	my $total=0;
	my $productStart;
	my $productStop;
	my $FIRST=1;
	my %productLibraryCounts;
	# foreach distinct read assocaited with that product:
	foreach my $storedRead (@{$productHash{$id}}) {
	    my($relStart,$relStop,$offset,$dStart,$readCount,$libraryCounts) = @{$storedRead};
	    $total += $readCount;
	    foreach my $label (keys %{$libraryCounts}) {
		$productLibraryCounts{$label} += $libraryCounts->{$label};
	    }
	    if($FIRST) {
		$FIRST=0;  
		$productStart = $relStart;
		$productStop = $relStop;
	    }
	}
	push(@productList,[$id,$total,$productStart,$productStop,\%productLibraryCounts]);
    }
    # now sort the products, by total, highest to lowest.
    if(@productList) {
	my @sortedProductList = sort {$b->[1] <=> $a->[1]} @productList;
	my($leftEnd,$rightEnd) = getFullHairpinWindow($fold,$leftCenter,$rightCenter);
	my @finalProductList;
	foreach my $idInfo (@sortedProductList) {
	    my($id,$totalCount,$productStart,$productStop,$productLibraryCounts) = @{$idInfo};
	    # get the relative position that is most abundant:
	    my($maxRelStart,$maxRelStop,$offset,$gStart) = getMaxProductPos($productHash{$id});
	    my($maxProductCount,$productCount) = getProductCounts($productHash{$id});
	    my $side = getProductSide($leftCenter,$rightCenter,$leftEnd,$rightEnd,$maxRelStart,$maxRelStop,$parameters);
	    push(@finalProductList,[$side,"",$productHash{$id},$productCount,$maxProductCount,$maxRelStart,$maxRelStop,$offset,$gStart,$productLibraryCounts]);
	}
	return \@finalProductList;
    } else {
	return [];
    }
}

sub extractProducts {
    my($leftCenter,$rightCenter,$fold,$basePairs,$location,$distinctReads,$parameters) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $productInfo = getProductInfo($leftCenter,$rightCenter,$fold,$basePairs,
				     $start,$stop,$strand,
				     $distinctReads,$parameters);
    $productInfo = rebuildProducts($leftCenter,$rightCenter,$fold,$basePairs,
				    $productInfo,$parameters);    
    return $productInfo;
}

sub getProductInfo {
    my($leftCenter,$rightCenter,$fold,$basePairs,
       $start,$stop,$strand,$distinctReads,$parameters) = @_;
    # first build a product hash:
    my $PRODUCTCOUNT = 1;
    my %productHash;
    foreach my $dRead (@{$distinctReads}) {
	# dStart is genomic position.
	my($dStart,$dStop,$total,$libraryCounts) = @{$dRead};
	my $relStart = $dStart - $start;
	my $relStop = $dStop - $start;
	# relStart is the position relative to the 5' end of readRegion. numbers start from zero.
	if($strand eq "-") {
	    $relStart = $stop - $dStop;
	    $relStop = $stop - $dStart;
	}
	# offset is relative to the leftMiddlePos (last 5' base pair in hairpin)
	my $offset = $relStop - $leftCenter;	
	if(($offset > 0)&&($relStop < $rightCenter)) {
	    $offset = 0;
	} elsif($offset > 0) {
	    $offset = $relStart - $rightCenter;
	}
	my $parsedRead = [$relStart,$relStop,$offset,$dStart,
			  $total,$libraryCounts];
	if(my $id = overlapCurrent(\%productHash,$parsedRead,$parameters)) {
	    push(@{$productHash{$id}},$parsedRead);	    
	} else {
	    push(@{$productHash{$PRODUCTCOUNT}},$parsedRead);
	    $PRODUCTCOUNT++;
	}
    }
    my @productList;
    foreach my $id (keys %productHash) {
	my $total=0;
	my $productStart;
	my $productStop;
	my $FIRST=1;
	foreach my $storedRead (@{$productHash{$id}}) {
	    my($relStart,$relStop,$offset,$dStart,$readCount,$libraryCounts) = @{$storedRead};
	    $total += $readCount;
	    if($FIRST) {
		$FIRST=0;  
		$productStart = $relStart;
		$productStop = $relStop;
	    }
	}
	push(@productList,[$id,$total,$productStart,$productStop]);
    }
    # now sort the products, by total, highest to lowset.                    
    my @sortedProductList = sort {$b->[1] <=> $a->[1]} @productList;
    my @productInfo;
    my($leftEnd,$rightEnd) = getFullHairpinWindow($fold,$leftCenter,$rightCenter);
    foreach my $idInfo (@sortedProductList) {
	my($id,$totalCount,$productStart,$productStop) = @{$idInfo};
	# get the relative position that is most abundant:
	my($maxRelStart,$maxRelStop,$offset,$gStart) = getMaxProductPos($productHash{$id});
	my($maxProductCount,$productCount) = getProductCounts($productHash{$id});
	my $side = getProductSide($leftCenter,$rightCenter,$leftEnd,$rightEnd,
				  $maxRelStart,$maxRelStop,$parameters);
	push(@productInfo,[$side,$productHash{$id},$productCount,$maxProductCount,$maxRelStart,$maxRelStop,$offset,$gStart]);
    }
    return \@productInfo;
}

sub rebuildProducts {
    my($leftCenter,$rightCenter,$fold,$basePairs,$productInfo,$parameters) = @_;
    my($leftEnd,$rightEnd) = getFullHairpinWindow($fold,$leftCenter,$rightCenter);
    my %usedMirSide;
    my %usedMorSide;
    my %newTypes;
    my @newProductInfo1;
    my @newProductInfo2;
    # sort by distance to loop, ascending order.
    my $totalReads = 0;
    foreach my $products (@{$productInfo}) {
	my($side,$productList,$productCount,$maxProductCount) = @{$products};
	$totalReads += $productCount;
    }   
    my @sortedProductInfo = sort {abs($a->[6]) <=> abs($b->[6])} @{$productInfo};
    foreach my $products (@sortedProductInfo) {
	my($side,$productList,$productCount,$maxProductCount,
	   $relStart,$relStop,$offset,$gStart) = @{$products};
	my $newType;
	#print "$side $relStart $relStop $offset\n";
	if($side eq "loop") {
	    $newType = "loop";
	} elsif($side eq "split") {
	    $newType = "split";
	} elsif($side eq "out") {
	    $newType = "out";
	} elsif(withinHairpin($leftEnd,$rightEnd,$relStart,$relStop,$parameters)) {
	    if(onHairpinArm($leftCenter,$rightCenter,$leftEnd,$rightEnd,
			    $relStart,$relStop,$parameters)) {		
		if($usedMirSide{$side}) {
		    if($usedMorSide{$side}) {
			$newType = "out";
		    } else {
			$newType = "moR";
			$usedMorSide{$side}++;
		    }
		} elsif($usedMirSide{otherSide($side)}) {
		    if(overlapsMir($leftCenter,$rightCenter,$basePairs,$side,$relStart,$relStop,\@newProductInfo1,$parameters)) {
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
		    if($productCount > 0.01*$totalReads) {
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
	} elsif(closeToHairpin($leftCenter,$rightCenter,$relStart,$relStop,$parameters)) {
	    if($usedMirSide{$side}) {
		if($usedMorSide{$side}) {
		    $newType = "out";
		} else {
		    $newType = "moR";
		    $usedMorSide{$side}++
		}
	    } elsif($usedMirSide{otherSide($side)}) {
		if(overlapsMir($leftCenter,$rightCenter,$basePairs,
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
			      $productCount,$maxProductCount,
			      $relStart,$relStop,$offset,$gStart]);
    }
    foreach my $products (@sortedProductInfo) {
	my($side,$productList,$productCount,$maxProductCount,
	   $relStart,$relStop,$offset,$gStart) = @{$products};
	if(($side eq "5p")||($side eq "3p")) {
	    if($newTypes{$relStart."x".$relStop} eq "loop") {
		if(overlapsMir($leftCenter,$rightCenter,$basePairs,$side,
			       $relStart,$relStop,\@newProductInfo1,$parameters)) {
		    $newTypes{$relStart."x".$relStop} = "miR";
		}
	    }
	}
	push(@newProductInfo2,[$side,$newTypes{$relStart."x".$relStop},$productList,
			       $productCount,$maxProductCount,
			       $relStart,$relStop,$offset,$gStart]);
    }
    return \@newProductInfo2;
}

sub getSimpleProductList {
    # this is a simpler version of getUnknownProductList for use when you don't know/need the fold information.
    my($location,$distinctReads,$parameters) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    # first build a product hash:
    my $PRODUCTCOUNT = 1;
    my %productHash;
    foreach my $dRead (@{$distinctReads}) {
	# dStart is genomic position.
	my($dStart,$dStop,$total,$libraryCounts) = @{$dRead};
	my $relStart = $dStart - $start;
	my $relStop = $dStop - $start;
	# relStart is the position relative to the 5' end of readRegion. numbers start from zero.
	if($strand eq "-") {
	    $relStart = $stop - $dStop;
	    $relStop = $stop - $dStart;
	}
	# offset is a placeholder for a more complex version of this function.
	my $offset = 0;
	my $parsedRead = [$relStart,$relStop,$offset,$dStart,$total,$libraryCounts];
	if(my $id = overlapCurrent(\%productHash,$parsedRead,$parameters)) {
	    push(@{$productHash{$id}},$parsedRead);	    
	} else {
	    push(@{$productHash{$PRODUCTCOUNT}},$parsedRead);
	    $PRODUCTCOUNT++;
	}
    }
    my @productList;
    # for each distinct product...
    foreach my $id (keys %productHash) {
	my $total=0;
	my $productStart;
	my $productStop;
	my $FIRST=1;
	my %productLibraryCounts;
	# foreach distinct read assocaited with that product:
	foreach my $storedRead (@{$productHash{$id}}) {
	    my($relStart,$relStop,$offset,$dStart,$readCount,$libraryCounts) = @{$storedRead};
	    $total += $readCount;
	    foreach my $label (keys %{$libraryCounts}) {
		$productLibraryCounts{$label} += $libraryCounts->{$label};
	    }
	    if($FIRST) {
		$FIRST=0;  
		$productStart = $relStart;
		$productStop = $relStop;
	    }
	}
	push(@productList,[$id,$total,$productStart,$productStop,\%productLibraryCounts]);
    }
    # now sort the products, by total, highest to lowest.
    if(@productList) {
	my @sortedProductList = sort {$b->[1] <=> $a->[1]} @productList;
	my @finalProductList;
	foreach my $idInfo (@sortedProductList) {
	    my($id,$totalCount,$productStart,$productStop,$productLibraryCounts) = @{$idInfo};
	    # get the relative position that is most abundant:
	    my($maxRelStart,$maxRelStop,$offset,$gStart) = getMaxProductPos($productHash{$id});
	    my($maxProductCount,$productCount) = getProductCounts($productHash{$id});
	    # note: the first two entries here a blank as they are not used or computable with the inputs to this simpler version
	    # of this function. In the getUnknownProducts defined above, it is the "side" and "type" of the products.
	    push(@finalProductList,["","",$productHash{$id},$productCount,$maxProductCount,$maxRelStart,$maxRelStop,$offset,$gStart,$productLibraryCounts]);
	}
	return \@finalProductList;
    } else {
	return [];
    }
}

sub getReverseProductDisplacement {
    my($productInfo,$revProductInfo) = @_;
    # define some parameters for what to consider.
    my $minDispInit = 20;
    my $minFrac = 0.0005;
    # get the total number of sense reads.
    my $totalReads = getTotalProductReads($productInfo);
    # first build a list of all pair-wise distances.
    my @dispList;
    for(my $i1=0;$i1<@{$productInfo};$i1++) {
	my($side,$newType,$prodList,$prodCount,$maxProdCount,$relStart,$relStop,$offset,$gStart) = @{$productInfo->[$i1]};
	#print "$side $newType @ $relStart:\n";
	for(my $i2=0;$i2<@{$revProductInfo};$i2++) {
	    my($aSide,$aNewType,$aProdList,$aProdCount,$aMaxProdCount,$aRelStart,$aRelStop,$aOffset,$aGStart) = @{$revProductInfo->[$i2]};
	    #print "\tcompared to: $aSide $aNewType @ $aRelStart\n";
	    if($aProdCount >= $minFrac*$totalReads) {
		# if the products overlap..
		if((($relStart <= $aRelStart)&&($aRelStart <= $relStop))||
		   (($relStart <= $aRelStop)&&($aRelStop <= $relStop))||
		   (($aRelStart <= $relStart)&&($relStart <= $aRelStop))||
		   (($aRelStart <= $relStop)&&($relStop <= $aRelStop))) {
		    my $disp = abs($aRelStart - $relStart);
		    #print "---> $disp\n";
		    push(@dispList,[$disp,$i1,$i2]);
		}
	    }
	}
    }
    @dispList = sort {$a->[0] <=> $b->[0]} @dispList;
    my %USED;
    my $sum = 0;
    my $total = 0;
    foreach my $overlap (@dispList) {
	my($disp,$i1,$i2) = @{$overlap};
	#print "checking: $disp $i1 $i2\n";
	unless(($USED{$i1})||($USED{$i2})) {
	    #print "\t--->adding $i1 $i2: $disp\n";
	    $sum += $disp;
	    $total++;
	    $USED{$i1}++;
	    $USED{$i2}++;
	}
    }
    return($sum,$total);
}

sub oldGetReverseProductDisplacement {
    my($productInfo,$revProductInfo) = @_;
    my $minDispInit = 20;
    my $minFrac = 0.0005;
    my $sum = 0;
    my $total = 0;
    my %USED;
    my $totalReads = getTotalProductReads($productInfo);
    for(my $i1=0;$i1<@{$productInfo};$i1++) {
	my($side,$newType,$prodList,$prodCount,$maxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$productInfo->[$i1]};
	my $start = $relStart;
	my $stop = $relStop;
	my $minDisp = $minDispInit;
	my $minI2;
	for(my $i2=0;$i2<@{$revProductInfo};$i2++) {
	    my($aType,$aNewType,$aProdList,$aProdCount,$aMaxProdCount,
	       $aRelStart,$aRelStop,$aOffset,$aGStart) = @{$revProductInfo->[$i2]};
	    if($aProdCount >= $minFrac*$totalReads) {
		my $aStart = $aRelStart;
		my $aStop = $aRelStop;
		# if the products overlap..
		if((($start <= $aStart)&&($aStart <= $stop))||
		   (($start <= $aStop)&&($aStop <= $stop))||
		   (($aStart <= $start)&&($start <= $aStop))||
		   (($aStart <= $stop)&&($stop <= $aStop))) {
		    my $disp = abs($aStart - $start);
		    unless($USED{$i2}) {
			if($disp < $minDisp) {
			    $minDisp = $disp;
			    $minI2 = $i2;
			}
		    }
		}
	    }
	}
	if($minDisp < $minDispInit) {
	    $sum += $minDisp;
	    $total++;
	    $USED{$minI2}++;
	}
    }
    return($sum,$total);
}

sub overlapCurrent {
    my($productHash,$parsedRead,$parameters)=@_;
    my $minDist = $parameters->{minDist} or die "minDist: parameter not loaded.";
    my($readRelStart,$readRelStop,$offset,$gStart,$count,$libraryCounts) = @{$parsedRead};
    # sort by the first entries total count ( the most abundant entry by construction, since
    # distinctReads is already sorted by abundance. ) highest to lowest...
    my @sortedIdList = sort {$productHash->{$b}[0][4] <=> $productHash->{$a}[0][4]} keys(%{$productHash});
    my @idList;
    foreach my $id (@sortedIdList) {
	my($sRelStart,$sRelStop,$sOffset,$sGStart,$sCount)=@{$productHash->{$id}[0]};
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
    my($leftCenter,$rightCenter,$leftEnd,$rightEnd,$relStart,$relStop,$parameters) = @_;    
    my $inBuffer = $parameters->{inHairpinBuffer} or die "no inHairpinBuffer loaded\n";
    my $outBuffer = $parameters->{outHairpinBuffer} or die "no outHairpinBuffer loaded\n";   
    my $hairpinRange = $parameters->{hairpinRange} or die "no hairpinRange loaded";
    if(withinHairpin($leftEnd,$rightEnd,$relStart,$relStop,$parameters)) {
	if(onHairpinArm($leftCenter,$rightCenter,$leftEnd,$rightEnd,$relStart,$relStop,$parameters)) {
	    if($relStop <= $leftCenter+$inBuffer) {
		return "5p";
	    } elsif($rightCenter-$inBuffer<=$relStart) {	    
		return "3p";	
	    } else {
		die "unknow situation reached in getProductSide!!! ouch.\n";
	    }
	} elsif(($relStart<=$leftCenter+1)&&($rightCenter-1<=$relStop)) {
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

sub getMaxProductPos {
    my $productList = shift;
    my($maxRelStart,$maxRelStop,$maxOffset,$maxGStart);
    my $FIRST = 1;
    # list is assumed to be sorted by abundance. must be sorted in getProductInfo()
    foreach my $read (@{$productList}) {
        my($sRelStart,$sRelStop,$sOffset,$sGStart,$sCount)=@{$read};
        if($FIRST) {
            $FIRST=0;
            $maxRelStart = $sRelStart;
	    $maxRelStop = $sRelStop;
	    $maxOffset = $sOffset;
	    $maxGStart = $sGStart;
	}
    }
    return($maxRelStart,$maxRelStop,$maxOffset,$maxGStart);
}

sub computeThreePrimeHet {
    my($productList) = @_;
    my $FIRST=1;
    my $threePrimeMaxPos;
    my $threePrimeMaxCount = 0;
    my $threePrimeTotalCount = 0;
    my %stopCount;
    foreach my $read (@{$productList}) {
	my($relStart,$relStop,$offset,$gStart,$count,$libraryCounts) = @{$read};
	$stopCount{$relStop} += $count;
    }
    # sort stops, greatest abundance to lowest.
    my @stops = sort {$stopCount{$b} <=> $stopCount{$a}} keys(%stopCount);
    my $topStop = shift(@stops);
    foreach my $read (@{$productList}) {
	my($relStart,$relStop,$offset,$gStart,$count,$libraryCounts) = @{$read};
	#print "3': comparing $relStop to $topStop\n";
	if($relStop == $topStop) {
	    $threePrimeMaxCount += $count;
	    #print "3': adding $count to max\n";
	}
        $threePrimeTotalCount += $count;
    }
    my $threePrimeHet = ($threePrimeTotalCount-$threePrimeMaxCount)/$threePrimeTotalCount;
    return $threePrimeHet;
}

sub computeFivePrimeHet {
    my($productList) = @_;
    my $FIRST=1;
    my $fivePrimeMaxPos;
    my $fivePrimeMaxCount = 0;
    my $fivePrimeTotalCount = 0;
    my %startCount;
    foreach my $read (@{$productList}) {
	my($relStart,$relStop,$offset,$gStart,$count,$libraryCounts) = @{$read};
	$startCount{$relStart} += $count;
    }
    my @starts = sort {$startCount{$b} <=> $startCount{$a}} keys(%startCount);
    my $topStart = shift(@starts);
    foreach my $read (@{$productList}) {
	my($relStart,$relStop,$offset,$gStart,$count,$libraryCounts) = @{$read};
	#print "5': comparing $relStart to $topStart\n";
	if($relStart == $topStart) {
	    $fivePrimeMaxCount += $count;
	    #print "5': adding $count to max\n";
	}
        $fivePrimeTotalCount += $count;
    }
    my $fivePrimeHet = ($fivePrimeTotalCount-$fivePrimeMaxCount)/$fivePrimeTotalCount;
    return $fivePrimeHet;
}

sub getTotalProductReads {
    my $productInfo = shift;
    my $totalReads = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$prodList,$prodCount,$maxProdCount,
	   $relStart,$relStop,$offset,$gStart) = @{$product};
	$totalReads += $prodCount;
    }
    return $totalReads;
}

sub getProductCounts {
    my $productList = shift;
    my $productCounts = 0;
    my $maxProductCount = 0;
    my $FIRST = 1;
    foreach my $read (@{$productList}) {
        my($sRelStart,$sRelStop,$sOffset,$sGStart,$sCount)=@{$read};
        $productCounts += $sCount;
        if($FIRST) {
            $FIRST=0;
            $maxProductCount = $sCount;
        }
    }
    return($maxProductCount,$productCounts);
}

sub getTotalLibraryCounts {
    my $productList = shift;
    my %totalLibraryCounts;
    foreach my $read (@{$productList}) {
	my($sRelStart,$sRelStop,$sOffset,$sGStart,$sCount,$libraryCounts)=@{$read};
	foreach my $sample (keys %{$libraryCounts}) {
	    $totalLibraryCounts{$sample} += $libraryCounts->{$sample};
	}
    }
    return \%totalLibraryCounts;
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
    my($leftCenter,$rightCenter,$basePairs,$side,$relStart,$relStop,$newProductInfo,$parameters) = @_;
    my $minShift = $parameters->{minShift} or die "FAIL: no minShift loaded.\n";
    foreach my $oProduct (@{$newProductInfo}) {
	my($oSide,$oNewType,$oProductlist,$oProductCount,$oMaxProdCount,$oRelStart,$oRelStop) = @{$oProduct};
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
    my($productInfo,$totalReads,$parameters) = @_;
    my $minOverlap = $parameters->{minOverlap} or die "FAIL: no minOverlap loaded.\n";
    my $minReadFraction = 0.01;
    my $countThreshold = $totalReads*$minReadFraction;
    my $maxShift = 0;
    for(my $i=0;$i<@{$productInfo};$i++) {
	my($side1,$newType1,$prodList1,$prodCount1,$maxProdCount1,
	   $relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	#print "$i: $side1 $newType1 $relStart1 $relStop1 $prodCount1 vs $countThreshold\n";
	if($prodCount1 >= $countThreshold) {	    
	    for(my $j=0;$j<@{$productInfo};$j++) {
		if($i != $j) {
		    my($side2,$newType2,$prodList2,$prodCount2,$maxProdCount2,
		       $relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
		    if($prodCount2 >= $countThreshold) {
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

sub getShift {
    my($relStart1,$relStop1,$relStart2,$relStop2) = @_;    
    return $relStart2-$relStart1;
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

sub withinHairpin {
    my($leftEnd,$rightEnd,$relStart,$relStop,$parameters) = @_;
    my $inBuffer = $parameters->{inHairpinBuffer};
    my $outBuffer = $parameters->{outHairpinBuffer};   
    if(($leftEnd-$outBuffer<=$relStart)&&
       ($relStop<=$rightEnd+$outBuffer)) {	    
	return 1;
    }
    return 0;
}

sub onHairpinArm {
    my($leftCenter,$rightCenter,$leftEnd,$rightEnd,$relStart,$relStop,$parameters) = @_;
    my $inBuffer = $parameters->{inHairpinBuffer};
    my $outBuffer = $parameters->{outHairpinBuffer};   
    if((($leftEnd-$outBuffer<=$relStart)&&($relStop<=$leftCenter+$inBuffer))||
       (($rightCenter-$inBuffer<=$relStart)&&($relStop<=$rightEnd+$outBuffer))) {	    
	return 1;
    }
    return 0;
}

sub closeToHairpin {
    my($leftCenter,$rightCenter,$relStart,$relStop,$parameters) = @_;
    my $hairpinRange = $parameters->{hairpinRange} or die "no hairpinRange loaded\n";
    if(($leftCenter-$hairpinRange <= $relStart)&&
       ($relStop <= $rightCenter+$hairpinRange)) {
	return 1;
    }
    return 0;
}

sub doubleStranded {
    my($start,$stop,$strand,$chromDataSet,$hitCount,$parameters) = @_;
    my $maxReverse = $parameters->{maxReverse};
    my $senseReadCount = readCount($start,$stop,$chromDataSet->{$strand},$hitCount);
    my $revStrand = revStrand($strand);
    my $reverseReadCount = readCount($start,$stop,$chromDataSet->{$revStrand},$hitCount);
    my $readCount = $senseReadCount + $reverseReadCount;
    # if reverse is less than 5%, it is singleStranded
    my $fraction = $readCount ? $reverseReadCount/$readCount : 0;
    #print "doubleStranded $readCount $senseReadCount $reverseReadCount $fraction\n";
    if($fraction <= $maxReverse) {
	return 0;
    } else {
	return 1;
    }
}

sub getBasePairs {
    my($leftCenter,$rightCenter,$fold) = @_;
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
    return \@basePairs;
}

sub getMaxHairpinLength {
    my($fold,$leftCenter,$rightCenter) = @_;
    my($leftEnd,$rightEnd) = extendHairpinToFullEnds($fold,$leftCenter,$rightCenter);
    my $leftLength = $leftCenter - $leftEnd + 1;
    my $rightLength = $rightEnd - $rightCenter + 1;
    return max($leftLength,$rightLength);
}

sub getFullHairpinWindow {
    # this method returns the window of paired bases, extends beyond minor hairpins    
    my($fold,$leftCenter,$rightCenter) = @_;;
    return extendHairpinToFullEnds($fold,$leftCenter,$rightCenter);
}

sub extendHairpinToFullEnds {
    my($fold,$leftCenter,$rightCenter) = @_;
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
        print STDERR "extendHairpinFullToEnds FAIL: ($leftCenter,$rightCenter)\n$fold\n$L\n$R\n\n";	
	return($leftCenter,$rightCenter);
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

sub getMaxLeftCenter {
    my $fold = shift;
    my $INRUN=0;
    my $maxLength = 0;
    my $maxRunPos=0;
    my $lastLeft = 0;
    for(my $i=0;$i<length($fold);$i++)  {
        if(substr($fold,$i,1) eq "(") {
	    $lastLeft = $i;
            $INRUN=1;
        } elsif(substr($fold,$i,1) eq ")") {
            if($INRUN) {
                my($leftEnd,$rightEnd) = extendHairpinToFullEnds($fold,$lastLeft,$i);
                my $length = $rightEnd - $leftEnd + 1;
                if($length > $maxLength) {
                    $maxLength = $length;
                    $maxRunPos = $lastLeft;
                }
            }
            $INRUN=0;
	}
    }
    return $maxRunPos;
}

sub getMaxRightCenter {
    my $fold = shift;
    my $leftCenter = getMaxLeftCenter($fold);
    for(my $i=$leftCenter;$i<length($fold);$i++)  {
        if(substr($fold,$i,1) eq ")") {
            return $i;
        }
    }
    return length($fold)-1;
}


#####################
# INPUT SUBROUTINES #
#####################

sub readMirRegionsFile {
    my $mirRegionsFile = shift;
    my %folds;
    my($id,$sequence,$tag,$location);
    open(MFILE,$mirRegionsFile) or die "could not open $mirRegionsFile\n";
    while(<MFILE>) {
        chomp;
	my($tag,$location,$leftCenter,$rightCenter,$mfe,$sequence,$fold) = split;
	my($chrom,$start,$stop,$strand) = parseLocation($location);
	push(@{$folds{$chrom}{$strand}},[$tag,$start,$stop,$leftCenter,$rightCenter,$mfe,$sequence,$fold]);
	#print "attempting to load: $tag $chrom $start $stop $strand\n";
    }
    close(MFILE);
    return \%folds;
}

sub loadBamFile {
    my $bamFile = shift;
    my $bam = Bio::DB::Sam->new( -bam  => $bamFile );
    return $bam;
}

sub loadBamList {
    my $bamListFile = shift;
    open(BLF,$bamListFile) or die "could not open $bamListFile\n";
    my @bamList;
    while(<BLF>) {
	chomp;
	my($label,$bamFile,$total) = split;
	my $bam = Bio::DB::Sam->new( -bam  => $bamFile );
	push(@bamList,[$label,$bam,$total,$bamFile]);
    }
    return \@bamList;
}

sub loadBamExpCtl {
    # it is a bamList where the first entry is a experiment, and the second is control.
    my $bamListFile = shift;
    open(BLF,$bamListFile) or die "could not open $bamListFile\n";
    my %bamExpCtl;
    while(<BLF>) {
	chomp;
	my($label,$bamFile,$total) = split;
	my $bam = Bio::DB::Sam->new( -bam  => $bamFile );
	if($bamExpCtl{"exp"}) {
	    $bamExpCtl{"ctl"} = [$label,$bam,$total,$bamFile];
	} else {
	    $bamExpCtl{"exp"} = [$label,$bam,$total,$bamFile];
	}
    }
    return \%bamExpCtl;
}

sub readLocusGff {
    my($locusGff) = @_;
    my %loci;
    open(LF,$locusGff) or die "could not open $locusGff";
    while(<LF>) {
	unless(/\#/) {
	    chomp;
	    my($chrom,$source,$type,$start,$stop,$score,$strand,$dot,$id)=split(/\t/,$_);	
	    my $uniqueId = $id;
	    if(/ID=\"(.*)\";/) {
		($uniqueId) = $id =~ /ID=\"(.*)\";/ or die "error1: could not parse the id $id\n";
	    } elsif(/ID=\"(.*)\"/) {
		($uniqueId) = $id =~ /ID=\"(.*)\"/ or die "error2: could not parse the id $id\n";
	    } elsif(/ID=(.*);/) {
		($uniqueId) = $id =~ /ID=(.*?);/ or die "error3: could not parse the id $id\n";
	    } elsif(/ID=(.*)/) {
		($uniqueId) = $id =~ /ID=(.*)/ or die "error4: could not parse the id $id\n";
	    } elsif(/Name=(.*);/) {
		($uniqueId) = $id =~ /Name=(.*?);/ or die "error5: could not parse the id $id\n";
	    } elsif(/Match .*/) {
		($uniqueId) = $id =~ /Match (.*)/;
	    } else {
		die "error6: could not parse the id:\n$_\n";
	    }
	    $chrom =~ s/chrMT/chrM/;
	    push(@{$loci{$chrom}},[$start,$stop,$strand,$uniqueId]);
	}
    }
    close(LF);
    return \%loci;
}

sub readChromLengths {
    my($genomeFile) = @_;
    my %chromLengths;
    my $tag;
    open(GFILE,$genomeFile) or die "could not open $genomeFile\n";
    while(<GFILE>) {
        chomp;
        if(/>/) {
            s/>//g;
            my @terms = split;
	    $tag = shift(@terms);
        } else {
            $chromLengths{$tag} += length($_);
        }
    }
    return \%chromLengths;
}

sub readSignificantProducts {
    my($significantProductFile,$parameters) = @_;
    my($expLabel,$ctlLabel);
    if($significantProductFile =~ /allProductSignificance/) {
	($expLabel,$ctlLabel) = $significantProductFile =~ /(.*)_vs_(.*)_allProductSignificance.txt/ or die "Could not parse significant products file.";
    } else {
	($expLabel,$ctlLabel) = $significantProductFile =~ /(.*)_vs_(.*)_(in|)(significant|remaining)Products.txt/ or die "Could not parse significant products file.";
    }
    $parameters->{expLabel} = $expLabel;
    $parameters->{ctlLabel} = $ctlLabel;
    open(SIGPRODUCTS,$significantProductFile) or die "could not open $significantProductFile\n";
    my %products;
    while(<SIGPRODUCTS>) {
	chomp;
	unless(/^#/) {
	       my($id,$location,$ctlCount,$expCount,$logFold,$pValue,$qValue,$rank,$expFiveHet,$expThreeHet) = split(/\t/);
	       my($chrom,$start,$stop,$strand) = parseLocation($location);
	       push(@{$products{$chrom}},[$id,$location,$ctlCount,$expCount,$pValue,$qValue,$rank,$expFiveHet,$expThreeHet]);
	   }
    }
    return \%products;
}

sub readSignificantProductsFile {
    my($productFile,$parameters) = @_;
    my($expLabel) = $productFile =~ /([^\/]*)_significantProducts.txt/ or die "Could not parse signficant products file.";
    $parameters->{expLabel} = $expLabel;
    open(PRODUCTS,$productFile) or die "could not open $productFile\n";
    my %products;
    while(<PRODUCTS>) {
	chomp;
	unless(/^#/) {
	       my($id,$location,$ctlCount,$expCount,$pValue,$qValue,$rank,$expFiveHet,$expThreeHet) = split(/\t/);
	       my($chrom,$start,$stop,$strand) = parseLocation($location);
	       push(@{$products{$chrom}},[$id,$location,$ctlCount,$expCount,$expFiveHet,$expThreeHet]);
	   }
    }
    return \%products;
}

sub readProductsFile {
    my($productFile,$parameters) = @_;
    my($expLabel) = $productFile =~ /([^\/]*)_products.(rename\.|)txt/ or die "Could not parse products file.";
    $parameters->{expLabel} = $expLabel;
    open(PRODUCTS,$productFile) or die "could not open $productFile\n";
    my %products;
    while(<PRODUCTS>) {
	chomp;
	unless(/^#/) {
	       my($id,$location,$expCount,$expFiveHet) = split(/\t/);
	       my($chrom,$start,$stop,$strand) = parseLocation($location);
	       push(@{$products{$chrom}},[$id,$location,$expCount,$expFiveHet]);
	   }
    }
    return \%products;
}

sub readProductsExpCtlFile {
    my($productFile,$parameters) = @_;
    # here the product file could be renamed, hence the regexp below:
    my($expLabel) = $productFile =~ /([^\/]*)_(products|allProducts)/ or die "Could not parse products file.";
    $parameters->{expLabel} = $expLabel;
    open(PRODUCTS,$productFile) or die "could not open $productFile\n";
    my %products;
    while(<PRODUCTS>) {
	chomp;
	unless(/^#/) {
	       my($id,$location,$ctlCount,$expCount,$expFiveHet,$expThreeHet) = split(/\t/);
	       my($chrom,$start,$stop,$strand) = parseLocation($location);
	       push(@{$products{$chrom}},[$id,$location,$ctlCount,$expCount,$expFiveHet,$expThreeHet]);
	   }
    }
    return \%products;
}

sub readMatureExpressionFile {
    my($productFile,$parameters) = @_;
    my($expLabel) = $productFile =~ /([^\/]*)_Expression\.txt/ or die "Could not parse mature expression file.";
    $parameters->{expLabel} = $expLabel;
    open(PRODUCTS,$productFile) or die "could not open $productFile\n";
    my %products;
    while(<PRODUCTS>) {
	chomp;
	unless(/^#/) {
	       my($id,$location,$expFiveHet,$expThreeHet,$sequence,$ctlCount,$expCount) = split(/\t/);
	       my($chrom,$start,$stop,$strand) = parseLocation($location);
	       push(@{$products{$chrom}},[$id,$location,$ctlCount,$expCount,$expFiveHet]);
	   }
    }
    return \%products;
}
