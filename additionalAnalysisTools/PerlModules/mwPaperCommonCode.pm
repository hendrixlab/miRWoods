package mwPaperCommonCode;
use miRWoods;
use strict;


sub loadDefaultParameters {
    my $parameters = { 
	"judgeOverlap" => "Novel",  
	"judgeAntisense" => "Novel",
	#"mirDeepCutoff" => 0,
	"minSignalToNoise" => 10,
	"minDist" => 10,
	"minCount" => 10,
	"minHet" => 0.1,
	"MSEProducts" => "One",  #may be set to "One" to get MSE for Major product or "All" for all products
	"minWTCount" => 5,
	"pseudoCount" => 0.015
    };
    return $parameters;
}

####################################################################
#         Euler plot creation functions                            #
####################################################################

sub getEulerIntersectStats {
    my($predictions1,$predictions2) = @_;
    my @both;
    my @pred1Only;
    my @pred2Only;
    my %pred2OverlapHash;
    foreach my $chrom (keys %{$predictions1}) {
	foreach my $hpInfo (@{$predictions1->{$chrom}}) {
	    my($start,$stop,$strand,$id,$annotation,$annotationType) = @{$hpInfo};
	    my $location = "$chrom:$start..$stop:$strand";
	    my $pred2OverlapId = findEulerOverlap($predictions2,$location);
	    if ($pred2OverlapId) {
		if ($pred2OverlapHash{$pred2OverlapId}) {
		    #only counting one overlap for the Euler plot
		    push(@pred1Only,$id);
		    print "pred2OverlapId overlapping $pred2OverlapHash{$pred2OverlapId} and $id\n";
		} else {
		    push(@both,[$id,$pred2OverlapId]);
		    $pred2OverlapHash{$pred2OverlapId} = $id;
		}
	    } else {
		push(@pred1Only,$id);
	    }
	}
    }
    foreach my $chrom (keys %{$predictions2}) {
	foreach my $hpInfo (@{$predictions2->{$chrom}}) {
	    my($start,$stop,$strand,$id,$annotation,$annotationType) = @{$hpInfo};
	    unless ($pred2OverlapHash{$id}) {
		push(@pred2Only,$id);
	    }
	}
    }
    return(\@pred1Only,\@pred2Only,\@both);
}

sub findEulerOverlap {
    my($miRDeepPredictions,$location) = @_;
    my($chrom,$start,$stop,$strand) = miRWoods::parseLocation($location);
    foreach my $hpInfo (@{$miRDeepPredictions->{$chrom}}) {
	my($start2,$stop2,$strand2,$id,$annotation,$annotationType) = @{$hpInfo};
	if (miRWoods::getOverlap($start,$stop,$start2,$stop2) && $strand eq $strand2) {
	    return $id;
	}
    }
    return 0;
}

#########################################################
#    Stats Functions                                   ##
#########################################################

sub printInfoFile {
    my($infoFile,$annotPredictions,$warning) = @_;
    open(INFLE,">$infoFile") or die "failed to open $infoFile for writing\n";
    my @annotated;
    my @antisense;
    my @overlaps;
    my @homologs;
    my @novel;
    foreach my $chrom (keys %{$annotPredictions}) {
	foreach my $predData (@{$annotPredictions->{$chrom}}) {
	    my($start,$stop,$strand,$id,$annotation,$annotationType) = @{$predData};
	    my $location = "$chrom:$start..$stop:$strand";
	    if ($annotationType eq "Annotated") {
		push(@annotated,[$id,$location,$annotation]);
	    } elsif ($annotationType eq "Antisense") {
		push(@antisense,[$id,$location,$annotation]);
	    } elsif ($annotationType eq "Overlap") {
		push(@overlaps,[$id,$location,$annotation]);
	    } elsif ($annotationType eq "Novel") {
		push(@novel,[$id,$location]);
	    } elsif ($annotationType eq "Homolog") {
		push(@homologs,[$id,$location,$annotation]);
	    } else {
		die "annotationType $annotationType not known in function printInfoFile";
	    }
	}
    }
    print INFLE "$warning\n\n";
    print INFLE "Annotated:\n";
    foreach my $annotInfo (@annotated) {
	my($id,$location,$annotation) = @{$annotInfo};
	print INFLE "$id\t$location\t$annotation\n";
    }
    print INFLE "\n\n\nHomologous:\n";
    foreach my $annotInfo (@homologs) {
	my($id,$location,$annotation) = @{$annotInfo};
	print INFLE "$id\t$location\t$annotation\n";
    }
    print INFLE "\n\n\nAntisense:\n";
    foreach my $annotInfo (@antisense) {
	my($id,$location,$annotation) = @{$annotInfo};
	print INFLE "$id\t$location\t$annotation\n";
    }
    print INFLE "\n\n\nAnnotation Overlaps:\n";
    foreach my $annotInfo (@overlaps) {
	my($id,$location,$annotation) = @{$annotInfo};
	print INFLE "$id\t$location\t$annotation\n";
    }
    print INFLE "\n\n\nNovel:\n";
    foreach my $annotInfo (@novel) {
	my($id,$location) = @{$annotInfo};
	print INFLE "$id\t$location\n";
    }
    close(INFLE);
}

sub calcPrecisionAndRecall {
    my($numAnnotated,$numNovel,$numAntisense,$numOverlaps,$numHomologs,$totalExpressed,$parameters) = @_;
    my($annoCount,$novelCount,$totalPos) = (0,0,0);
    $annoCount += $numAnnotated;
    $novelCount += $numNovel;
    $novelCount += $numHomologs;
    if ($parameters->{"judgeOverlap"} eq "Annotated") {
	$annoCount += $numOverlaps;
    } elsif ($parameters->{"judgeOverlap"} eq "Novel") {
	$novelCount += $numOverlaps;
    } elsif ($parameters->{"judgeOverlap"} eq "None") {
	#do nothing
    } else {
	die "error: judgeOverlap parameter not properly set";
    }
    if ($parameters->{"judgeAntisense"} eq "Annotated") {
	$annoCount += $numAntisense;
    } elsif ($parameters->{"judgeAntisense"} eq "Novel") {
	$novelCount += $numAntisense;
    } elsif ($parameters->{"judgeAntisense"} eq "None") {
	#do nothing
    } else {
	die "error: judgeAntisense parameter not properly set";
    }
    my $numPred = $annoCount + $novelCount;
    my $recall = $annoCount / $totalExpressed;
    my $precision = $annoCount / $numPred;
    my $falseNegatives = $totalExpressed - $annoCount;
    my $f1Score = calculateF1($annoCount,$novelCount,$falseNegatives);
    return($annoCount,$novelCount,$totalExpressed,$precision,$recall,$f1Score);
}

sub countExpressedAnnotations {
    my($hpExpressionInfo) = @_;
    my $totalExpressed = 0;
    foreach my $name (keys %{$hpExpressionInfo}) {
   	my($hpId,$hpLocation,$hpExpCount,$hpNormCount) = @{$hpExpressionInfo->{$name}};
	if ($hpExpCount) {
	    $totalExpressed++;
	}
    }
    return $totalExpressed;
}

sub countAnnotationTypes {
    my($annotPredictions) = @_; 
    my($numAnnotated,$numNovel,$numAntisense,$numOverlaps,$numHomologs) = (0,0,0,0,0);
    foreach my $chrom (keys %{$annotPredictions}) {
	foreach my $predData (@{$annotPredictions->{$chrom}}) {
	    my($start,$stop,$strand,$id,$annotation,$annotationType) = @{$predData};
	    if ($annotationType eq "Annotated") {
		$numAnnotated++;
	    } elsif ($annotationType eq "Antisense") {
		$numAntisense++;
	    } elsif ($annotationType eq "Overlap") {
		$numOverlaps++;
	    } elsif ($annotationType eq "Novel") {
		$numNovel++;
	    } elsif ($annotationType eq "Homolog") {
		$numHomologs++;
	    } else {
		die "error: annotation $annotationType not expected in countAnnotationTypes()";
	    }
	}
    }
    return($numAnnotated,$numNovel,$numAntisense,$numOverlaps,$numHomologs);
}

sub calculateF1 {
    my($TP,$FP,$FN) = @_;
    my $posPredictions = $TP + $FP;
    my $totalPositives = $TP + $FN;
    my $recall = $TP / $totalPositives;
    if ($posPredictions == 0) {
	return 0;
    }
    my $precision = $TP / $posPredictions;
    my $f1num = $precision * $recall;
    my $f1denom = $precision + $recall;
    my $f1score;
    if ($f1num == 0) {
	return 0;
    } else {
	my $f1score = $f1num / $f1denom;
	return $f1score;
    }
}

###########################################################
#  Annotation Functions                                   #
###########################################################

sub getKnownAndNovelPredictions {
    my($annotPredictions,$parameters) = @_;
    my %knownHairpins;
    my %novelHairpins;
#    my $knownCount = 0;
#    my $novelCount = 0;
    foreach my $chrom (keys %{$annotPredictions}) {
	foreach my $annotData (@{$annotPredictions->{$chrom}}) {
	    my($start,$stop,$strand,$id,$annotation,$annotationType) = @{$annotData};
	    if ($annotationType eq "Overlap") {
		if ($parameters->{"judgeOverlap"} eq "Annotated") {
		    push(@{$knownHairpins{$chrom}},[$start,$stop,$strand,$id,$annotation,$annotationType]);
#		    $knownCount++;
		} elsif ($parameters->{"judgeOverlap"} eq "Novel") {
		    push(@{$novelHairpins{$chrom}},[$start,$stop,$strand,$id,$annotation,$annotationType]);
#		    $novelCount++;
		} elsif ($parameters->{"judgeOverlap"} eq "None") {
			#do nothing
		} else {
		    die "error: judgeOverlap parameter not properly set";
		}
	    } elsif ($annotationType eq "Antisense") {
		if ($parameters->{"judgeAntisense"} eq "Annotated") {
		    push(@{$knownHairpins{$chrom}},[$start,$stop,$strand,$id,$annotation,$annotationType]);
#		    $knownCount++;
		} elsif ($parameters->{"judgeAntisense"} eq "Novel") {
		    push(@{$novelHairpins{$chrom}},[$start,$stop,$strand,$id,$annotation,$annotationType]);
#		    $novelCount++;
		} elsif ($parameters->{"judgeAntisense"} eq "None") {
		    #do nothing
		} else {
		    die "error: judgeAntisense parameter not properly set";
		}
	    } elsif ($annotationType eq "Annotated") {
		push(@{$knownHairpins{$chrom}},[$start,$stop,$strand,$id,$annotation,$annotationType]);
#		    $knownCount++;
	    } elsif ($annotationType eq "Novel") {
		push(@{$novelHairpins{$chrom}},[$start,$stop,$strand,$id,$annotation,$annotationType]);
#		    $novelCount++;
	    } elsif ($annotationType eq "Homolog") {
		push(@{$novelHairpins{$chrom}},[$start,$stop,$strand,$id,$annotation,$annotationType]);
	    } else {
		die "Error: $annotationType not a known annotation in getKnownAndNovelPredicitons()\n";
	    }
	}
    }
#   print "known = $knownCount\tnovel = $novelCount\n";
    return(\%knownHairpins,\%novelHairpins);
}

sub checkForDuplicatePredictions {
    my($annotPredictions,$warnings) = @_;
    my %annotCount;
    my %newAnnotPredictions;
    my %duplicatePredictions;
    my $duplicatedPredictionCount = 0;
    #searching for duplicate predicitons overlapping an annotation
    foreach my $chrom (keys %{$annotPredictions}) {
	foreach my $annotData (@{$annotPredictions->{$chrom}}) {
	    my($start,$stop,$strand,$id,$annotation,$annotationType) = @{$annotData};
	    if ($annotCount{$annotation}) {
		$annotCount{$annotation} += 1;
		$duplicatedPredictionCount++;
	    } elsif ($annotation) {
		$annotCount{$annotation} = 1;
	    }
	}
    }
    #fixing annotations where more than one $annotationType for an $annotation shows up as "Annotated"
    #in these cases they will be changed to "Overlap".
    #additionally the number of duplicate predictions summed in $duplicatePredictionCount
    foreach my $chrom (keys %{$annotPredictions}) {
	my %correctFoldFound;
	foreach my $annotData (@{$annotPredictions->{$chrom}}) {
	    my($start,$stop,$strand,$id,$annotation,$annotationType) = @{$annotData};
	    if ($correctFoldFound{$annotation}) {
		if ($annotationType eq "Annotated") {
		    $$warnings .= "Warning [Annot_Exists]: the annotation for $annotation already exist.  Changing $id to overlap\n";
		    $annotationType = "Overlap";
		    @{$annotData} = ($start,$stop,$strand,$id,$annotation,$annotationType);
		}
	    } elsif ($annotationType eq "Annotated") {
		$correctFoldFound{$annotation} = 1;
	    }
	    push(@{$newAnnotPredictions{$chrom}},$annotData);
	    #adding to duplicatePredictions so that warnings can be created using ids of all duplicates
	    if (($annotCount{$annotation}) && $annotCount{$annotation} > 1) {
		my $location = "$chrom:$start..$stop:$strand";
		push(@{$duplicatePredictions{$annotation}},[$id,$annotationType,$location]); 
	    }
	}
    }
    #creating warnings
    foreach my $annotation (keys %duplicatePredictions) {
	$$warnings .= "warning [SAME_ANNOTATION]: the following ids are annotated as $annotation - ";
	foreach my $dupData (@{$duplicatePredictions{$annotation}}) {
	    my($id,$annotationType,$location) = @{$dupData};
	    $$warnings .= "$id($annotationType) ";
	}
	$$warnings .= "\n";
    }
    return(\%newAnnotPredictions,\%duplicatePredictions,$duplicatedPredictionCount);
}

sub getOverlappingAnnotation {
    my($location,$hairpins) = @_;
    my($chrom,$start,$stop,$strand) = miRWoods::parseLocation($location);
    #checking sense strand
    my $bestAnnotation = 0;
    my $bestOverlap = 0;
    foreach my $hairpinInfo (@{$hairpins->{$chrom}}) {
	my($annoStart,$annoStop,$annoStrand,$id,$name) = @{$hairpinInfo};
	my $annotationOverlap = miRWoods::getOverlap($start,$stop,$annoStart,$annoStop);
	if ($annotationOverlap > $bestOverlap && $strand eq $annoStrand) {
	    $bestOverlap = $annotationOverlap;
	    $bestAnnotation = $name;
	}
    }
    if ($bestAnnotation) {
	return $bestAnnotation;
    }
    #checking antisense strand
    foreach my $hairpinInfo (@{$hairpins->{$chrom}}) {
	my($annoStart,$annoStop,$annoStrand,$id,$name) = @{$hairpinInfo};
	my $annotationOverlap = miRWoods::getOverlap($start,$stop,$annoStart,$annoStop);
	if ($annotationOverlap > $bestOverlap && $strand ne $annoStrand) {
	    $bestOverlap = $annotationOverlap;
	    $bestAnnotation = "as-$name";
	}
    }
    if ($bestAnnotation) {
	return $bestAnnotation;
    }
    return 0;
}

sub provideCommonAnnotation {
    #this function is used to make certain that annotations are being fairly applied to all pipelines
    #and will be shared by miRWoods, miRDeep, mireap, etc.
    my($location,$fold,$annotHairpins,$annotProducts) = @_;
    my $minMatCoverage = 0.9; #$parameters->{aliasMinMatCoverage};
    my $annotationType;
    my $annotation = getOverlappingAnnotation($location,$annotHairpins);
    if ($annotation) {
	my($hpAnnotationId,$hpAnnotationName,$hairpinCoverage,
	 $annotCoverage,$annotMatCoverage,$numAnnotProducts,$annotProdsOverlapped) = miRWoods::annotateHairpin($location,$fold,
													       $annotHairpins,$annotProducts);
	if (($numAnnotProducts == 1 && $annotCoverage >= 0.50) ||
	    ($numAnnotProducts == 1 && $hairpinCoverage >= 1) || 
	    ($numAnnotProducts > 1 && $annotProdsOverlapped == $numAnnotProducts)) {
	    if ( $hpAnnotationName =~ /^as\-/ ) {
		$annotationType = "Antisense";
	    } else {
		$annotationType = "Annotated";
#		print "$hpAnnotationName\t$location\t$hairpinCoverage\t$annotCoverage\t$annotMatCoverage\t$numAnnotProducts\t$annotProdsOverlapped\n";
	    }
	} else {
	    if ( $hpAnnotationName =~ /^as\-/ ) {
		$annotationType = "Antisense";
#		print "$hpAnnotationName\t$location\t$hairpinCoverage\t$annotCoverage\t$annotMatCoverage\t$numAnnotProducts\t$annotProdsOverlapped\n";
	    } else {
		$annotationType = "Overlap";
	    }
	}
    } else {
	$annotationType = "Novel";
    }
    return($annotation,$annotationType);
}

####################################################
#       General Homology Functions                 #
####################################################

sub addHomology {
    my($annotPredictions,$homologyInfo) = @_;
    my %newAnnotPredictions;
    foreach my $chrom (keys %{$annotPredictions}) {
	foreach my $annotData (@{$annotPredictions->{$chrom}}) {
	    my($start,$stop,$strand,$id,$annotation,$annotationType) = @{$annotData};
	    if ($homologyInfo->{$id} && $annotationType eq "Novel") {
		$annotData->[4] = $homologyInfo->{$id};
		$annotData->[5] = "Homolog";
	    }
	    push(@{$newAnnotPredictions{$chrom}},$annotData);
	}
    }
    return \%newAnnotPredictions;
}

sub readHomologyFile {
    my($homologyFile) = @_;
    my %homology;
    open(HF,$homologyFile) or die "failed to open $homologyFile\n";
    while (<HF>) {
	chomp;
	my($hpName,$familyName,$eVal,$closestHomolog) = split(/\t/);
	if ($closestHomolog ne "-") {
	    if ($familyName eq "na" || $familyName eq "Novel") {
		$homology{$hpName} = $closestHomolog;
	    } else {
		$homology{$hpName} = $familyName;
	    }
	}
    }
    close(HF);
    return \%homology;
}

###########################################
#   Annotated Expression Functions        #
###########################################


sub readExpressionFile {
    my($expressionFile) = @_;
    my %matureExpressionInfo;
    my %hpExpressionInfo;
    my %hpIds;
    my %hpLocations;
    open(EPF,$expressionFile) or die "failed to open $expressionFile\n";
    while (<EPF>) {
	chomp;
	unless ( /^#/ ) {
	    my($hairpinName,$hairpinId,$hairpinLocation,$mirname,$mirid,$location,$expcount,$normcount) = split(/\t/);
	    push(@{$matureExpressionInfo{$hairpinName}},[$mirname,$mirid,$location,$expcount,$normcount]);
	    $hpIds{$hairpinName} = $hairpinId; 
	    $hpLocations{$hairpinName} = $hairpinLocation;
	}
    }
    foreach my $hairpinName (keys %matureExpressionInfo) {
	my($hpExpCount,$hpNormCount) = (0,0);
	foreach my $matureInfo (@{$matureExpressionInfo{$hairpinName}}) {
	    my($mirname,$mirid,$location,$expcount,$normcount) = @{$matureInfo};
	    $hpExpCount += $expcount;
	    $hpNormCount += $normcount;
	}
	my $hpId = $hpIds{$hairpinName};
	my $hpLocation = $hpLocations{$hairpinName};
	if ($hpExpressionInfo{$hairpinName}) {
	    print "Warning: $hairpinName already found in %hpExpressionInfo hash\n";
	}
	@{$hpExpressionInfo{$hairpinName}} = ($hpId,$hpLocation,$hpExpCount,$hpNormCount);
    }
    close(EPF);
    return(\%hpExpressionInfo,\%matureExpressionInfo);
}

sub countExpressed {
    my($hpExpressionInfo) = @_;
    my $count = 0;
    foreach my $name (keys %{$hpExpressionInfo}) {
   	my($hpId,$hpLocation,$hpExpCount,$hpNormCount) = @{$hpExpressionInfo->{$name}};
	if ($hpExpCount) {
	    $count++;
	}
    }
    return $count;
}

##################################################
#    miRWoods File Parsing Functions             #
##################################################

sub addMiRWoodsHomology {
    my($mwAnnotPredictions,$homologyFile) = @_;
    my $homologyInfo = readHomologyFile($homologyFile);
    $mwAnnotPredictions = addHomology($mwAnnotPredictions,$homologyInfo);
    return $mwAnnotPredictions;
}

sub annotateMirWoods {
    #Note: annotatateMirWoods retrieves annotations through provideCommonAnnotation in order to
    #guarantee that the same annotation functions are used to annotate predictions in all pipelines.
    my($mwPredictions,$mwFolds,$annoHairpins,$annoProducts,$warning) = @_;
    my %mwAnnotPredictions;
    my %annotationHash;
    foreach my $chrom (keys %{$mwPredictions}) {
	foreach my $predictionInfo (@{$mwPredictions->{$chrom}}) {
	    my($start,$stop,$strand,$id,$name) = @{$predictionInfo};
	    my $location = "$chrom:$start..$stop:$strand";
	    my $foldLocation = miRWoods::getFoldGenomicLocation($location,$mwFolds->{$id});
	    my($locChrom,$locStart,$locStop,$locStrand) = miRWoods::parseLocation($foldLocation);
	    my $fold = '.' x ($locStop - $locStart + 1);  #just to keep the annotation process as similar to mirdeep as possible
	    my($annotation,$annotationType) = mwPaperCommonCode::provideCommonAnnotation($foldLocation,$fold,$annoHairpins,$annoProducts);
	    my $mwPredictionClass = "Novel";
	    if ( $name =~ /verlap/ ) {
		$mwPredictionClass = "Overlap";
	    } elsif ( $name =~ /as-/ ) {
		$mwPredictionClass = "Antisense";
	    } elsif ( $name =~ /mir-/ ) {
		$mwPredictionClass = "Annotated";		
	    } elsif ( $name =~ /let-/ ) {
		$mwPredictionClass = "Annotated";		
	    }
	    #the same set of warning should probably be applied to other pipelines being tested
	    if ($annotationType ne $mwPredictionClass) {
		$$warning .= "warning [DIFF_ANNO_TYPE]: $id $name predicted as $mwPredictionClass but annotated as $annotationType\n";
	    }
	    #we are just using the fold region, not the surounding area
	    push(@{$mwAnnotPredictions{$chrom}},[$locStart,$locStop,$locStrand,$id,$annotation,$annotationType]);
	}
    }
    return \%mwAnnotPredictions;
}

sub getFolds {
    my($hairpinsFile) = @_;
    my %mwFolds;
    open(HP,$hairpinsFile) or die "failed to open $hairpinsFile\n";
    while (<HP>) {
	chomp;
	unless ( /^#/ ) {
	    my($tag,$chrom,$start,$stop,$strand,$lc,$rc,$totalSense,$totalAntisense,$mfe,$seq,$fold) = split(/\t/);
	    $mwFolds{$tag} = $fold;
	}
    }
    close(HP);
    return \%mwFolds;
}

sub getMiRWoodsMatProducts {
    my($mwPredictions,$productFile) = @_;
    my $products = miRWoods::readProductFile($productFile);
    my %mwProducts;
#    my $genomeDir = '/nfs0/Hendrix_Lab/Human/Genome/hg19';
    foreach my $chrom (keys %{$mwPredictions}) {
#	my $genome = miRWoods::loadGenome("$genomeDir/$chrom.fa");	
	foreach my $hairpinInfo (@{$mwPredictions->{$chrom}}) {
	    my($start,$stop,$strand,$id) = @{$hairpinInfo};
	    my $location = "$chrom:$start..$stop:$strand";
	    foreach my $productInfo (sort {$b->[3] <=> $a->[3]} @{$products->{$id}}) {
		my($tag,$side,$type,$total,$totalMostAbundant,$adjTotal,$adjTotalMostAbundant,$relStart,$relStop,$prodStrand,$prodSeq) = @{$productInfo};
		if ($strand eq $prodStrand) {
		    if ($type =~ /miR/) {
			my $prodName = "$tag-$side-$type";
			my($chrom,$start,$stop,$strand) = miRWoods::parseLocation($location);
			my($prodStart,$prodStop);
			if ($strand eq '+') {
			    $prodStart = $start + $relStart;
			    $prodStop = $start + $relStop;
			} elsif ($strand eq '-') {
			    $prodStop = $stop - $relStart;
			    $prodStart = $stop - $relStop;
			}
			my $matProdLocation = "$chrom:$prodStart..$prodStop:$strand";
			push(@{$mwProducts{$id}},[$chrom,$prodStart,$prodStop,$prodStrand,$prodName,$prodName,$prodSeq]);	      
		    }
		}
	    }
	}
    }
    return \%mwProducts;
}

sub parseMiRWoodsResults {
    my($posPredGFF) = @_;
    my %mwPredictions;
    open(PPGFF,$posPredGFF) or die "could not open $posPredGFF\n";
    while(<PPGFF>) {
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
		my $alias = $info{Alias} or die "No Alias found for the line:\n$_\n";
		my $name = $info{Name} or die "No Name found for the line:\n$_\n";
		push(@{$mwPredictions{$chrom}},[$start,$stop,$strand,$id,$alias]); 
	    }
	}
    }
    close(PPGFF);
    return \%mwPredictions;
}

sub testMiRWoodsProducts {
    my($mwPredictions,$mwProducts,$genomeDir) = @_;
    foreach my $chrom (keys %{$mwPredictions}) {
	print "chrom = $chrom\n";
	my $genome = miRWoods::loadGenome("$genomeDir/$chrom.fa");
	foreach my $mdInfo (@{$mwPredictions->{$chrom}}) {
	    my($start,$stop,$strand,$id) = @{$mdInfo};
	    foreach my $mdProdInfo (@{$mwProducts->{$id}}) {
		my($matChrom,$matStart,$matStop,$matStrand,$matId,$matName,$sequence) = @{$mdProdInfo};
		my $matLocation = "$matChrom:$matStart..$matStop:$matStrand";
		my $testSeq = miRWoods::getSequence($matLocation,$genome);
		$testSeq =~ tr/acgtACGT/acguacgu/;
		$sequence =~ tr/acgtACGT/acguacgu/;
		if ($sequence ne $testSeq) {
		    print "$testSeq ne $sequence for $matId $matLocation\n";
		}
	    }
	}
    }    
}

##################################################
#    miRDeep  File Parsing Functions             #
##################################################

sub addMiReapHomology {
    my($mrAnnotPredictions,$homologyFile) = @_;
    my $homologyInfo = readHomologyFile($homologyFile);
    $mrAnnotPredictions = addHomology($mrAnnotPredictions,$homologyInfo);
    return $mrAnnotPredictions;
}

sub getHomologyFromMDPredictions {
    my($mdPredictions) = @_;
    my %mdPredHomology;
    foreach my $chrom (keys %{$mdPredictions}) {
	foreach my $predictionInfo (@{$mdPredictions->{$chrom}}) {
	    my($id,$score,$TPProb,$rfamAlert,$total,$matureCount,$loopCount,$starCount,$randfoldSignificant,$miRNA,$otherSpecMiRNA,$UCSCbrowser,$NCBIblastn,$mature,$star,$precursor,$location,$mdPredictionClass) = @{$predictionInfo};
	    unless ($miRNA eq "-" && $otherSpecMiRNA eq "-") {
		$mdPredHomology{$id} = ($miRNA eq "-") ? $otherSpecMiRNA : $miRNA;
	    }
	}
    }
    return \%mdPredHomology;
}

sub addMiRDeepHomology {
    my($mdPredictions,$mdAnnotPredictions,$homologyFile) = @_;
    my $homologyInfo = readHomologyFile($homologyFile);
    $mdAnnotPredictions = addHomology($mdAnnotPredictions,$homologyInfo);
#    my $mdPredHomology = getHomologyFromMDPredictions($mdPredictions);
#    $mdAnnotPredictions = addHomology($mdAnnotPredictions,$mdPredHomology);
    return $mdAnnotPredictions;
}


sub annotateMirDeep {
    #Note: annotatateMirDeep retrieves annotations through provideCommonAnnotation in order to
    #guarantee that the same annotation functions are used to annotate predictions in all pipelines.
    my($mdPredictions,$annoHairpins,$annoProducts,$warning) = @_;
    my %mdAnnotPredictions;
    my %annotationHash;
    foreach my $chrom (keys %{$mdPredictions}) {
	foreach my $predictionInfo (@{$mdPredictions->{$chrom}}) {
	    my($id,$score,$TPProb,$rfamAlert,$total,$matureCount,$loopCount,$starCount,$randfoldSignificant,$miRNA,$otherSpecMiRNA,$UCSCbrowser,$NCBIblastn,$mature,$star,$precursor,$location,$mdPredictionClass) = @{$predictionInfo};
	    my($locChrom,$start,$stop,$strand) = miRWoods::parseLocation($location);
	    my $fold = '.' x ($stop - $start + 1);
	    my($annotation,$annotationType) = mwPaperCommonCode::provideCommonAnnotation($location,$fold,$annoHairpins,$annoProducts);
	    #the same set of warning should probably be applied to other pipelines being tested
	    if ($annotationType ne $mdPredictionClass) {
		$$warning .= "warning [DIFF_ANNO_TYPE]: $id predicted as $mdPredictionClass but annotated as $annotationType\n";
	    }
	    push(@{$mdAnnotPredictions{$chrom}},[$start,$stop,$strand,$id,$annotation,$annotationType]);
	}
    }
    return \%mdAnnotPredictions;
}

sub determineMirdeepCutoff {
    my($mirDeepAccData,$parameters) = @_;
    my $minSignalToNoise = $parameters->{minSignalToNoise} or die "minSignalToNoise parameter not set\n";
    my $cutoffScore = ~0;
    my $minNoise = ~0;
    for my $score (keys %{$mirDeepAccData}) {
	my($novelReported,$estimatedFalsePositives,$estimatedTruePositives,
	   $miRsInSpecies,$miRsInData,$knownDetected,$signalToNoise,$exisionGearing) = @{$mirDeepAccData->{$score}};
	if ($minSignalToNoise <= $signalToNoise && $signalToNoise <= $minNoise) {
	    $cutoffScore = $score;
	    $minNoise = $signalToNoise;
	}
    }
    print "mirdeep cutoff Score = $cutoffScore\n";
    return $cutoffScore;
}

sub applyMirdeepCutoff {
    my($mirDeepAccData,$mdPredictions,$parameters) = @_;
    #my $mirDeepCutoffScore = $parameters->{mirDeepCutoff};
    my $mirDeepCutoffScore = determineMirdeepCutoff($mirDeepAccData,$parameters);
    my %newMDPredictions;
    foreach my $chrom (%{$mdPredictions}) {
	foreach my $mdInfo (@{$mdPredictions->{$chrom}}) {
	    my($id,$score) = @{$mdInfo};
	    if ($score >= $mirDeepCutoffScore) {
		push(@{$newMDPredictions{$chrom}},$mdInfo);
	    }
	}
    }
    return \%newMDPredictions;
}

sub getMDProductLocation {
    my($location,$precursor,$productSeq) = @_;
    my($chrom,$start,$stop,$strand) = miRWoods::parseLocation($location);
    $start++;  #converting to 1-based
    #Getting the start and stop positions relative to the hairpin
    my $matLen = length($productSeq);
    my $precLen = length($precursor);
    my $lenDiff = $precLen - $matLen;
    my $seq1 = substr($precursor,0,$matLen);
    my $seq2 = substr($precursor,$lenDiff,$matLen);
    my($relStart,$relStop) = (-1,-1);
    my $side;
    if ($productSeq eq $seq1) {
	$side = "5p";
	($relStart,$relStop) = (0,($matLen - 1));
    } elsif ($productSeq eq $seq2) {
	$side = "3p";
	($relStart,$relStop) = ($lenDiff,($precLen - 1));
    } else {
	die "Error: $productSeq\t$precursor\t$seq1\t$seq2\n";
    }
    #Getting the product locations within the genome
    my($prodStart,$prodStop);
    if ($strand eq '+') {
	$prodStart = $start + $relStart;
	$prodStop = $start + $relStop;
    } elsif ($strand eq '-') {
	$prodStop = $stop - $relStart;
	$prodStart = $stop - $relStop;
    }
    my $prodLocation = "$chrom:$prodStart..$prodStop:$strand";
    return($prodLocation,$side);
}

sub parseMirDeepResults {
    my($mirDeepResults) = @_;
    my %mirDeepAccData;  #accuracy data keyed by Score
    my %mirDeepPredictions;
    my %mirDeepProducts;
    open(MDR,$mirDeepResults) or die "failed to open $mirDeepResults\n";
    my $parsingPredictedAccuracy = 0;
    my $parsingNovelPredictions = 0;
    my $parsingKnownPredictions = 0;
    my $parsingNonDetectedMirs = 0;
    while (<MDR>) {
	chomp;
	my $line = $_;
	#determining the type of line to parse
	if ($line =~ /^miRDeep2\ score/) {
	    $parsingPredictedAccuracy = 1;
	    $parsingNovelPredictions = 0;
	    $parsingKnownPredictions = 0;
	    $parsingNonDetectedMirs = 0;
	} elsif ($line =~ /^novel\ miRNAs\ predicted\ by\ miRDeep2/) {
	    $line = <MDR>;  #removing extra headers
	    $parsingPredictedAccuracy = 0;
	    $parsingNovelPredictions = 1;
	    $parsingKnownPredictions = 0;
	    $parsingNonDetectedMirs = 0;
	} elsif ($line =~ /^mature\ miRBase\ miRNAs\ detected\ by\ miRDeep2/) {
	    $line = <MDR>;  #removing extra headers
	    $parsingPredictedAccuracy = 0;
	    $parsingNovelPredictions = 0;
	    $parsingKnownPredictions = 1;
	    $parsingNonDetectedMirs = 0;
	} elsif ($line =~ /^#miRBase\ miRNAs\ not\ detected\ by\ miRDeep2/) {
	    last;  #currently the non-detected mirs are not output properly by miRDeep
#	    $parsingPredictedAccuracy = 0;
#	    $parsingNovelPredictions = 0;
#	    $parsingKnownPredictions = 0;
#	    $parsingNonDetectedMirs = 1;
	} elsif ($line ne "") {
	    #parsing the lines based on their type
	    if ($parsingPredictedAccuracy) {
		my($score,$novelReported,$estimatedFalsePositives,$estimatedTruePositives,$miRsInSpecies,$miRsInData,$knownDetected,$signalToNoise,$exisionGearing) = split(/\t/,$line);
		@{$mirDeepAccData{$score}} = ($novelReported,$estimatedFalsePositives,$estimatedTruePositives,$miRsInSpecies,$miRsInData,$knownDetected,$signalToNoise,$exisionGearing);
	    } elsif ($parsingNovelPredictions || $parsingKnownPredictions) {
		my($id,$score,$TPProb,$rfamAlert,$total,$matureCount,$loopCount,$starCount,$randfoldSignificant,$miRNA,$otherSpecMiRNA,$UCSCbrowser,$NCBIblastn,$mature,$star,$precursor,$location) = split(/\t/,$line);
		my($chrom) = miRWoods::parseLocation($location);
		my $mdPredictionClass = ($parsingNovelPredictions) ? "Novel" : "Annotated";
		push(@{$mirDeepPredictions{$chrom}},[$id,$score,$TPProb,$rfamAlert,$total,$matureCount,$loopCount,$starCount,$randfoldSignificant,$miRNA,$otherSpecMiRNA,$UCSCbrowser,$NCBIblastn,$mature,$star,$precursor,$location,$mdPredictionClass]);
		my($matProdLocation,$matProdSide) = getMDProductLocation($location,$precursor,$mature);
		my($matChrom,$matStart,$matStop,$matStrand) = miRWoods::parseLocation($matProdLocation);
		push(@{$mirDeepProducts{$id}},[$matChrom,$matStart,$matStop,$matStrand,"$id-mat-$matProdSide","$id-mat-$matProdSide",$mature]);
		my($starProdLocation,$starProdSide) = getMDProductLocation($location,$precursor,$star);
		my($starChrom,$starStart,$starStop,$starStrand) = miRWoods::parseLocation($starProdLocation);
		push(@{$mirDeepProducts{$id}},[$starChrom,$starStart,$starStop,$starStrand,"$id-star-$starProdSide","$id-star-$starProdSide",$star]);	      
	    } else {
		last;
	    }
	}
    }
    close(MDR);
    return(\%mirDeepAccData,\%mirDeepPredictions,\%mirDeepProducts);
}

sub testMiRDeepProducts {
    my($mirDeepPredictions,$mirDeepProducts,$genomeDir) = @_;
    foreach my $chrom (keys %{$mirDeepPredictions}) {
	print "chrom = $chrom\n";
	my $genome = miRWoods::loadGenome("$genomeDir/$chrom.fa");
	foreach my $mdInfo (@{$mirDeepPredictions->{$chrom}}) {
	    my($id) = @{$mdInfo};
	    foreach my $mdProdInfo (@{$mirDeepProducts->{$id}}) {
		my($matChrom,$matStart,$matStop,$matStrand,$matId,$matName,$sequence) = @{$mdProdInfo};
		my $matLocation = "$matChrom:$matStart..$matStop:$matStrand";
		my $testSeq = miRWoods::getSequence($matLocation,$genome);
		$testSeq =~ tr/acgtACGT/acguacgu/;
		if ($sequence ne $testSeq) {
		    print "$testSeq ne $sequence for $matId $matLocation\n";
		}
	    }
	}
    }
}

##################################################
#    miReap  File Parsing Functions              #
##################################################



sub annotateMiReap {
    #Note: annotatateMiReap retrieves annotations through provideCommon Annotation in order to
    #guarantee that the same annotation functions are used to annotate predictions in all pipelines.
    my($mrPredictions,$annoHairpins,$annoProducts,$warning) = @_;
    my %mrAnnotPredictions;
    my %annotationHash;
    foreach my $chrom (keys %{$mrPredictions}) {
	foreach my $predictionInfo (@{$mrPredictions->{$chrom}}) {
	    my($start,$stop,$strand,$id,$name) = @{$predictionInfo};
	    my $location = "$chrom:$start..$stop:$strand";
	    my $fold = '.' x ($stop - $start + 1);  #just to keep the annotation process as similar to mirdeep as possible
	    my($annotation,$annotationType) = mwPaperCommonCode::provideCommonAnnotation($location,$fold,$annoHairpins,$annoProducts);
	    #we are just using the fold region, not the surounding area
	    push(@{$mrAnnotPredictions{$chrom}},[$start,$stop,$strand,$id,$annotation,$annotationType]);
	}
    }
    return \%mrAnnotPredictions;
}

sub parseMiReapResults {
    my($posPredGFF) = @_;
    my %mrPredictions;
    my %mrProducts;
    open(PPGFF,$posPredGFF) or die "could not open $posPredGFF\n";
    while(<PPGFF>) {
	unless(/^\#/) {
	    chomp;
	    my($chrom,$source,$type,$start,$stop,$score,$strand,$phase,$info) = split(/\t/);
	    my %info;
	    my @terms = split(/;/,$info);
	    foreach my $term (@terms) {
		my($key,$value) = $term =~ /(.*)=(.*)/;
		$info{$key} = $value;
	    }
	    if($type eq "precursor") {
		# hairpin region from the gff file
		my $id = $info{ID} or die "No ID found for the line:\n$_\n"; 
		my $alias = $info{ID} or die "No Alias found for the line:\n$_\n";
		my $name = $info{ID} or die "No Name found for the line:\n$_\n";
		push(@{$mrPredictions{$chrom}},[$start,$stop,$strand,$id,$alias]); 
	    } elsif ($type =~ /^mature/) {
		my $id = $info{ID} or die "No ID found for the line:\n$_\n"; 
		my $alias = $info{ID} or die "No Alias found for the line:\n$_\n";
		my $name = $info{ID} or die "No Name found for the line:\n$_\n";
		my $parentId = $info{Parent} or die "No Parent found for the line:\n$_\n";
		my $sequence = $info{Seq} or die "No Seq found for the line:\n$_\n";
		push(@{$mrProducts{$parentId}},[$chrom,$start,$stop,$strand,$id,$name,$sequence]);
	    }
	}
    }
    close(PPGFF);
    return(\%mrPredictions,\%mrProducts);
}

#############################################################
#   Bam Processing and Dicer Knockdown comparison Functions #
#############################################################

sub compareOtherDkrd {
    my($majorProductRegions,$dkrdOutputFile,$bam,$bamTotal,$dkrdBam,$dkrdBamTotal,$parameters) = @_;
    my $pseudoCount = $parameters->{pseudoCount};
    my $minWTCount = $parameters->{minWTCount};
    open(DKRD,">$dkrdOutputFile") or die "failed to open $dkrdOutputFile for writing\n";
    print DKRD "Tag\twtExpCount\tdkrdExpCount\tlfc\n";
    foreach my $chrom (%{$majorProductRegions}) {
	foreach my $productInfo (@{$majorProductRegions->{$chrom}}) {
	    my($prodName,$prodLocation) = @{$productInfo};
	    my($chrom,$start,$stop,$strand) = miRWoods::parseLocation($prodLocation);
	    my($wtExpCount,$wtAdjCount) = getAdjBamProductCount($bam,$chrom,$start,$stop,$strand,$parameters);
	    if ($wtExpCount >= $minWTCount) {
		my($dkrdExpCount,$dkrdAdjCount) = getAdjBamProductCount($dkrdBam,$chrom,$start,$stop,$strand,$parameters);
		my $wtAdjARPM = $wtAdjCount / ($bamTotal/1e6);
		my $dkrdAdjARPM = $dkrdAdjCount / ($dkrdBamTotal/1e6);
		my $change = ($dkrdAdjARPM + $pseudoCount) / ($wtAdjARPM + $pseudoCount);
		my $lfc = log2($change);
		print DKRD "$prodName\t$wtAdjARPM\t$dkrdAdjARPM\t$lfc\n";
	    }
	}
    }
    close(DKRD);
}

sub comparePredicitonsDkrd {
    my($hairpins,$products,$dkrdOutputFile,$bam,$bamTotal,$dkrdBam,$dkrdBamTotal,$parameters) = @_;
    my $pseudoCount = $parameters->{pseudoCount};
    my $minWTCount = $parameters->{minWTCount};
    open(DKRD,">$dkrdOutputFile") or die "failed to open $dkrdOutputFile for writing\n";
    print DKRD "Tag\twtExpCount\tdkrdExpCount\tlfc\n";
    foreach my $chrom (keys %{$hairpins}) {
	foreach my $hpInfo (@{$hairpins->{$chrom}}) {
	    my($start,$stop,$strand,$id) = @{$hpInfo};
	    my @lfcData;
	    foreach my $prodInfo (@{$products->{$id}}) {
		my($prodChrom,$prodStart,$prodStop,$prodStrand,$prodId,$prodName) = @{$prodInfo};
		my($wtExpCount,$wtAdjCount) = getAdjBamProductCount($bam,$prodChrom,$prodStart,$prodStop,$prodStrand,$parameters);
		my($dkrdExpCount,$dkrdAdjCount) = getAdjBamProductCount($dkrdBam,$prodChrom,$prodStart,$prodStop,$prodStrand,$parameters);
		my $wtAdjARPM = $wtAdjCount / ($bamTotal/1e6);
		my $dkrdAdjARPM = $dkrdAdjCount / ($dkrdBamTotal/1e6);
		my $change = ($dkrdAdjARPM + $pseudoCount) / ($wtAdjARPM + $pseudoCount);
		my $lfc = log2($change);
		push(@lfcData,[$prodId,$wtExpCount,$dkrdExpCount,$wtAdjCount,$dkrdAdjCount,$wtAdjARPM,$dkrdAdjARPM,$lfc]);
	    }
	    my @sortedLFCData = sort {$b->[1] <=> $a->[1]} @lfcData;
	    if ($parameters->{MSEProducts} eq "One") {
		#printing MSE for major product only
		my $lfcInfo = $sortedLFCData[0];
		my($prodId,$wtExpCount,$dkrdExpCount,$wtAdjCount,$dkrdAdjCount,$wtAdjARPM,$dkrdAdjARPM,$lfc) = @{$lfcInfo};
		if ($wtExpCount >= $minWTCount) {
		    print DKRD "$prodId\t$wtAdjARPM\t$dkrdAdjARPM\t$lfc\n";
		}
	    } elsif ($parameters->{MSEProducts} eq "All") {
		#printing MSE for all products
		foreach my $lfcInfo (@sortedLFCData) {
		    my($prodId,$wtExpCount,$dkrdExpCount,$wtAdjCount,$dkrdAdjCount,$wtAdjARPM,$dkrdAdjARPM,$lfc) = @{$lfcInfo};
		    if ($wtExpCount >= $minWTCount) {
			print DKRD "$prodId\t$wtAdjARPM\t$dkrdAdjARPM\t$lfc\n";
		    }
		}
	    } else {
		#Error
		die "Error: MSEProducts parameter not set\n";
	    }
	}
    }
}

sub getAdjBamProductCount {
    my($bam,$chrom,$start,$stop,$strand,$parameters) = @_;
    my $minDist = $parameters->{minDist} or die "minDist: parameter not loaded.";
    my @alignments = $bam->get_features_by_location(-seq_id => $chrom,
						    -start  => $start - $minDist,
						    -end    => $stop + 2*$minDist);
    my $count = 0;
    my $adjCount = 0;
    my %fivePosHash;
    my %threePosHash;
    foreach my $read (@alignments) {
	my $rStart = $read->start;
	my $rStop = $read->end;
	my $hitCount = $read->get_tag_values('NH');
	# 5' end methdod.
	my $rStrand = miRBamTools::mapBamStrand($read->strand);
	if($strand eq $rStrand) {
	    if($strand eq "+") {
		# for forward strand, 5' end is "start"
		if(($start - $minDist <= $rStart)&&($rStart <= $start + $minDist)) {
		    $count++;
		    $adjCount += 1/$hitCount;
		    $fivePosHash{$rStart}++;
		    $threePosHash{$rStop}++;
		}
	    } else {
		# for reverse strand, 5' end is "stop"
		if(($stop - $minDist <= $rStop)&&($rStop <= $stop + $minDist)) {
		    $count++;
		    $adjCount += 1/$hitCount;
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
	return($count,$adjCount,$fiveHet,$threeHet);
    } else {
	return (0,0,1,1);
    }
}

sub getBamTotal {
    my($bamFile) = @_;
    my $mappedReads = `samtools view -F 0x904 $bamFile | awk '{a[\$1]=1} END {print length(a)}'`;
    chomp($mappedReads);
    return $mappedReads;
}

sub readProductScoresFile {
    my($productScoresFile) = @_;
    my %products;
    open(PS,$productScoresFile) or die "failed to open $productScoresFile for reading\n";
    my $header = <PS>;
    chomp($header);
    while (<PS>) {
	chomp;
	my($geneId) = split(/\t/,$_);
	my($name,$location) = miRWoods::parseGeneId($geneId);
	my($chrom,$start,$stop,$strand) = miRWoods::parseLocation($location);
	push(@{$products{$chrom}},[$name,$location]);
    }
    close(PS);
    return \%products;
}

sub getMajorProductRegions {
    my($bedFile,$bamFile,$totalMapped) = @_;
    my %majorProductRegions;
    my $miRWoodsParameters = miRWoods::loadDefaultParameters();
    #i'm doing the following since im using some miRWoods functions and the way
    #miRWoods loads bam files is slightly different from miRBamTools.pm
    my @bamList;
    my $bam = miRWoods::loadBamFile($bamFile);
    my $bamHeader = $bam->header;
    my $bamIndex = miRWoods::loadBamIndex($bamFile);
    push(@bamList,[$bam,$bamHeader,$bamIndex,$bamFile,"otherRNA",$totalMapped]);
    #now gathering major product regions for each region in the RNA bed file
    open(BED, $bedFile) or die "failed to open $bedFile\n";
    while (<BED>) {
        chomp;
	unless(/\#/) {	  
	    my($chrom,$start,$stop,$rrId,$score,$strand) = split(/\t/,$_);
	    $start++; #converting to one based
	    my $location = "$chrom:$start..$stop:$strand";
#	    print "$chrom\t$start\t$stop\t$rrId\t$score\t$strand\n";
	    my($distinctReads,$uniqueReadCount,$totalDistinctReads,$adjustedReadCount) = miRWoods::retrieveReadData($location,\@bamList,$miRWoodsParameters);
	    if ($distinctReads->{$strand}) {  
		my($products) = miRWoods::buildProducts($location,$distinctReads->{$strand},$strand,$miRWoodsParameters);
		my @sortedProducts = sort {$b->[4] <=> $a->[4]} @{$products};
		#my($id,$productReads,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = $sortedProducts[0];
		my $prodLocation = miRWoods::getProductGenomicLocation($sortedProducts[0],$chrom);
		push(@{$majorProductRegions{$chrom}},[$rrId,$prodLocation]);
#		print $prodLocation . "\n";
	    }
	}
    }
    close(BED);
    return \%majorProductRegions;
}

sub log2 {
    my $num = shift;
    return log($num) / log(2); 
}


1;


