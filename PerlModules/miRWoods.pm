package miRWoods;
use Bio::DB::Sam;
#use Memory::Usage;
use Statistics::R;
use Storable;
use RNA;
use strict;
use warnings;
#use threads;
#use threads::shared;
use List::Util;


sub loadDefaultParameters {
    my $parameters = {
	"outputPrefix" => "readRegions",
	"mirbaseGff" => "",
	"bamListFile" => "",
	"RepeatRegionsFile" => "",
	"genomeDir" => "",
	"SizesFile" => "",
	"geneModels" => "",
	"validList" => "",
	"knownMirsFile" => "",
	"otherAnnotFileList" => "",
	"bucketSize" => 100000,  #the span of a genome that makes up one bucket.  Used for long lists of annotations to save time.
	"minReadLength" => 0,  #min size of a read to retrieve from a bam file
	"maxReadLength" => ~0, #max size of a read to retrieve from a bam file
	"prrWindowLength" => 200000, #the window size for retrieval of reads in printReadRegions()
	"prrMaxLength" => 160,  #lengths greater than prrMaxLength are placed in the "*_long.txt" read regions file.
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
	"bamMaxMismatces" => ~0,   #the maximum number of mismatches to accept within a seed region (see bamMismatchSeedLength) for a bam file
	"bamMismatchSeedLength" => ~0,  #seperate seed length to check for mismatches in bam entry (note: this is different from seed used in bowtie)
	"MaxShortProductSize" => "14",  #the maximum size of a product which should be appended with the word short (ie. short-miR, short-loop, etc.)
	"minLongProductSize" => "28",   #the minimum size of a product which should be appended with the word long (ie. long-miR, long-loop, etc.)
	"LoadFromConfigFile" => "",
	"maxSplitProductAbundance" => 0.05,
	"maxOverlapProductAbundance" => 0.05, 
	"minAntisenseProduct" => 0,   #the minimum amount of antisense product required to incorporate a sense antisense overlapping pair into the aapd
	"HPDVCutoff" => 0.5,      #the cutoff for the hairpin decision value
	"ARPMCutoff" => 0,        #the cutoff for the minimum ARPM to consider for positive predictions
	"LetThroughCutoff" => 2,  #wont have an effect unless it is set in the config file to a value between [0..1]
	"minMajProdCount" => 1,   #minimum count of major product to be considered a positive prediction (see classifyPredictions())
	"minAdjMajProdCount" => 1, #minimum count of adj major product to be considered a positive prediction (see classifyPredictions())
	"noAbundCheckMiRCutoff" => 2, #abundance of major product will not be checked if number of mirs >= noAbundCheckMiRCutoff (see classifyPredictions())
	"verbose" => 0,           #for some functions this will print out extra lines.
	"aliasMinHPCoverage" => 0,   #the amount an annotated precursor must overlap the fold of the hairpin
	"aliasMinPrecCoverage" => 0, #the amount a fold must overlap an annotated precursror
	"aliasMinMatCoverage" => 0.9,  #the amount the hairpin must overlap the mature annotations for a precursor
	"aliasMinHPCoverageSingle" => 0.7,  #same as aliasMinHPCoverage, applied in the case where a precurser only has one annotated mature mir
	"aliasMinPrecCoverageSingle" => 0, #same as aliasMinPrecCoverage, applied in the case where a precurser only has one annotated mature mir
	"aliasMinMatCoverageSingle" => 0.9, #same as aliasMinMatCoverage, applied in the case where a precurser only has one annotated mature mir
	"trainModels" => 0,              #will make miRWoods train new product and hairpin models when set to 1
	"hairpinRF" => "hairpinRF.model",  #the file containting the random forest model for the hairpin phase
	"productRF" => "productRF.model",  #the file containting the random forest model for the product phase
	"MLAlgorithm" => "RandomForest",  #like an appendix this parameter will one day be removed
    };
    return $parameters;
}

sub loadDefaultMLProductParameters {
    #these parameters determine which variables will be used in machine learning for products
    my $prodMLParameters = {
	"adjMaxProductCount" => 0,
	"adjProductCount"=> 0,
	"fivePrimeHet"=> 1,
	"length"=> 0,
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

sub loadHairpinReadScoreParameters {
    my $readScoreParameters = {
	"aapd"=> 1,
	"tapd"=> 1,
	"urf"=> 1,
	"ahc"=> 1,
	"afh"=> 1,
	"sameShift"=> 1,
	"bothShift"=> 1,
	"SPA"=> 1,
	"OPA"=> 1,
	"APV"=> 1, 
	"wAPV"=> 1,
	"ARV"=> 1,
	"wARV"=> 1,
	"totalOverlapAmount"=> 1,
	"averageOverlapAmount"=> 1,
	"totalRelativeOverlapAmount"=> 1,
	"averageRelativeOverlapAmount"=> 1,
	"maxOverlap"=> 1,
	"totalSenseRPM" => 1,
	"totalAntisenseRPM" => 1,
	"RFProductAvg"=> 1
    };
    return $readScoreParameters;
}

sub loadHairpinStructureParameters {
    my $hpStructureParameters = {
	"mfe"=> 1,
	"pbp"=> 1,
	"duplexEnergy"=> 1,
	"foldDupCmp"=> 1, 
	"exactFoldDupCmp"=> 1,
	"dupPBP"=> 1,
	"dupLoopLength"=> 1,
	"dupOverlap"=> 1,
	"mpLoopDistance"=> 1,
	"dupLoopDistance"=> 1,
	"loopSize"=> 1,
	"innerLoopGapCount"=> 1
    };
    return $hpStructureParameters;
}

sub loadHPSequenceParameters {
    my $hpSequenceParameters = {
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
	"GCcontent"=> 1,
    };
    return $hpSequenceParameters;
}

sub loadHPStructSeqParameters {
    my $hpStructSeqParameters = {
	"mfe"=> 1,
	"pbp"=> 1,
	"duplexEnergy"=> 1,
	"foldDupCmp"=> 1, 
	"exactFoldDupCmp"=> 1,
	"dupPBP"=> 1,
	"dupLoopLength"=> 1,
	"dupOverlap"=> 1,
	"mpLoopDistance"=> 1,
	"dupLoopDistance"=> 1,
	"loopSize"=> 1,
	"innerLoopGapCount"=> 1,
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
	"GCcontent"=> 1,
    };
    return $hpStructSeqParameters;
}

sub loadDefaultMLHairpinParameters_prodOverlaps {
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
	"exactFoldDupCmp"=> 0,
	"dupPBP"=> 1,
	"dupLoopLength"=> 1,
	"APV"=> 1, 
	"wAPV"=> 1,
	"ARV"=> 1,
	"wARV"=> 1,
	"mpCount"=> 0,
	"dupCount"=> 0,
	"dupOverlap"=> 0,
	"mpLoopDistance"=> 1,  
	"dupLoopDistance"=> 1,  
	"loopSize"=> 1,
	"mpOverlapScore"=> 0,
	"dupOverlapScore"=> 0,
	"wghtMPOverlapIntensity"=> 0, 
	"wghtDupOverlapIntensity"=> 0,
	"OPA2"=> 0,
	"totalOverlapAmount"=> 0,
	"averageOverlapAmount"=> 1, 
	"totalRelativeOverlapAmount"=> 1, 
	"averageRelativeOverlapAmount"=> 0,
	"maxOverlap"=> 0, 
	"innerLoopGapCount"=> 1,
	"totalSenseRPM" => 1,
	"totalAntisenseRPM" => 1,
	"zScore" => 0,
	"overlapSizeRatioSum" => 0,
	"totalOverlap" => 1,
	"maxBulge" => 1,
	"maxInteriorLoop" => 1,
	"intLoopSideDiff" => 1,
	"maxUnboundOverhang" => 1,
	"numOffshoots" => 1,
	"dupSize" => 1,
	"neighborCount" => 1,
	"RFProductAvg" => 1,
	"rel5pOutCount" => 1, 
	"rel5pMorCount" => 1, 
	"rel5pMirCount" => 1, 
	"relLoopCount" => 1, 
	"relSplitCount" => 1, 
	"rel3pMirCount" => 1, 
	"rel3pMorCount" => 1, 
	"rel3pOutCount" => 1, 
	"extraOutCount" => 0,
#	"miRmoR5pOverlap" => 0,  
#	"miRmoR3pOverlap" => 0, 
#	"miRLoop5pOverlap" => 0, 
#	"miRLoop3pOverlap" => 0, 
#	"loopLoopOverlap" => 0, 
#	"out5pOverlap" => 0, 
#	"out3pOverlap" => 0, 
#	"outOut5pOverlap" => 0, 
#	"outOut3pOverlap" => 0, 
#	"inProdOverlap" => 0, 
#	"miRSplit5pOverlap" => 0, 
#	"miRSplit3pOverlap" => 0, 
#	"moRSplit5pOverlap" => 0, 
#	"moRSplit3pOverlap" => 0, 
#	"splitLoopOverlap" => 0, 
#	"loopSplitOverlap" => 0
	"miRmoR5pOverlap" => 1,  
	"miRmoR3pOverlap" => 1, 
	"miRLoop5pOverlap" => 1, 
	"miRLoop3pOverlap" => 1, 
	"loopLoopOverlap" => 1, 
	"out5pOverlap" => 1, 
	"out3pOverlap" => 0, 
	"outOut5pOverlap" => 1, 
	"outOut3pOverlap" => 1, 
	"inProdOverlap" => 1, 
	"miRSplit5pOverlap" => 1, 
	"miRSplit3pOverlap" => 1, 
	"moRSplit5pOverlap" => 0, 
	"moRSplit3pOverlap" => 0, 
	"splitLoopOverlap" => 0, 
	"loopSplitOverlap" => 0

    };
    return $hairpinMLParameters;
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
	"mpLoopDistance"=> 0,
	"dupLoopDistance"=> 0,
	"loopSize"=> 1,
	"mpOverlapScore"=> 0,
	"dupOverlapScore"=> 0,
	"wghtMPOverlapIntensity"=> 1, 
	"wghtDupOverlapIntensity"=> 1, #10 real vs 58
	"OPA2"=> 0,
	"totalOverlapAmount"=> 1,
	"averageOverlapAmount"=> 1,
	"totalRelativeOverlapAmount"=> 1,
	"averageRelativeOverlapAmount"=> 1,
	"maxOverlap"=> 1,
	"innerLoopGapCount"=> 1,
	"totalSenseRPM" => 1,
	"totalAntisenseRPM" => 1,
	"zScore" => 1,
	"overlapSizeRatioSum" => 1,
	"totalOverlap" => 1,
	"maxBulge" => 1,
	"maxInteriorLoop" => 1,
	"intLoopSideDiff" => 1,
	"maxUnboundOverhang" => 1,
	"numOffshoots" => 0,
	"dupSize" => 1,
	"neighborCount" => 1,
	"RFProductAvg"=> 1
    };
    return $hairpinMLParameters;
}

sub loadDefaultMLHairpinParameters_test_old {
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
	"dupPBP"=> 0,
	"dupLoopLength"=> 1,
	"APV"=> 0, 
	"wAPV"=> 0,
	"ARV"=> 1,
	"wARV"=> 1,
	"mpCount"=> 0,
	"dupCount"=> 0,
	"dupOverlap"=> 0,
	"mpLoopDistance"=> 0,  
	"dupLoopDistance"=> 0,  
	"loopSize"=> 1,
	"mpOverlapScore"=> 0,
	"dupOverlapScore"=> 0,
	"wghtMPOverlapIntensity"=> 0, 
	"wghtDupOverlapIntensity"=> 0,
	"OPA2"=> 0,
	"totalOverlapAmount"=> 0,
	"averageOverlapAmount"=> 1,
	"totalRelativeOverlapAmount"=> 1,
	"averageRelativeOverlapAmount"=> 0,
	"maxOverlap"=> 1,
	"innerLoopGapCount"=> 1,
	"totalSenseRPM" => 1,
	"totalAntisenseRPM" => 0,
	"zScore" => 0,
	"RFProductAvg"=> 1
    };
    return $hairpinMLParameters;
}

sub loadDefaultMLHairpinParameters_test {
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
	"exactFoldDupCmp"=> 0,
	"dupPBP"=> 1,
	"dupLoopLength"=> 1,
	"APV"=> 1, 
	"wAPV"=> 1,
	"ARV"=> 1,
	"wARV"=> 1,
	"mpCount"=> 0,
	"dupCount"=> 0,
	"dupOverlap"=> 0,
	"mpLoopDistance"=> 0,  
	"dupLoopDistance"=> 0,  
	"loopSize"=> 1,
	"mpOverlapScore"=> 0,
	"dupOverlapScore"=> 0,
	"wghtMPOverlapIntensity"=> 0, 
	"wghtDupOverlapIntensity"=> 0,
	"OPA2"=> 0,
	"totalOverlapAmount"=> 0,
	"averageOverlapAmount"=> 1, 
	"totalRelativeOverlapAmount"=> 1, 
	"averageRelativeOverlapAmount"=> 0,
	"maxOverlap"=> 0, 
	"innerLoopGapCount"=> 1,
	"totalSenseRPM" => 1,
	"totalAntisenseRPM" => 0,
	"zScore" => 0,
	"overlapSizeRatioSum" => 0,
	"totalOverlap" => 1,
	"maxBulge" => 1,
	"maxInteriorLoop" => 1,
	"intLoopSideDiff" => 1,
	"maxUnboundOverhang" => 1,
	"numOffshoots" => 0,
	"dupSize" => 1,
	"neighborCount" => 1,
	"RFProductAvg" => 1
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
	"foldDupCmp"=> 1, #10 lost vs 2
	"exactFoldDupCmp"=> 1,
	"dupPBP"=> 1,
	"dupLoopLength"=> 1,
	"APV"=> 1, #
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
	"wghtDupOverlapIntensity"=> 1, #10 real vs 58
	"OPA2"=> 0,
	"totalOverlapAmount"=> 1,
	"averageOverlapAmount"=> 1,
	"totalRelativeOverlapAmount"=> 1,
	"averageRelativeOverlapAmount"=> 1,
	"maxOverlap"=> 1,
	"innerLoopGapCount"=> 1,
#	"totalSenseRPM" => 1,
#	"totalAntisenseRPM" => 1,
	"totalSenseRPM" => 0,
	"totalAntisenseRPM" => 0,
	"zScore" => 1,
	"RFProductAvg"=> 1
    };
    return $hairpinMLParameters;
}

sub loadDefaultMLHairpinParameters_otherRF {
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
	"totalSenseRPM" => 1,
	"totalAntisenseRPM" => 1,
#	"totalSenseRPM" => 0,
#	"totalAntisenseRPM" => 0,
	"zScore" => 1,
	"piRNA" => 0,
	"piRNA_vs_mir" => 0,
	"rRNA" => 0,
	"rRNA_vs_mir" => 0,
	"snoRNA" => 0,
	"snoRNA_vs_mir" => 0,
	"tRNA" => 0,
	"tRNA_vs_mir" => 0,
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
    $parameters->{hpAnnotFile} = $filePrefix . "_hairpinAnnotations.txt";
    $parameters->{prodAnnotFile} = $filePrefix . "_productAnnotations.txt";
    $parameters->{otherAnnotFile} = $filePrefix . "_otherAnnotations.txt";
    $parameters->{featureVectorFile} = $filePrefix . "_scores.txt";
    $parameters->{productFeatVectorFile} = $filePrefix . "_productScores.txt";
    $parameters->{productTrainFile} = $filePrefix . "_productTrainFile.txt";
    $parameters->{predProductFile} = $filePrefix . "_predProductRegions.txt";
    $parameters->{predProductClasses} = $filePrefix . "_predProductClasses.txt";
    $parameters->{hairpinVectorFile} = $filePrefix . "_hairpinScores.txt";
    $parameters->{newHairpinVectorFile} = $filePrefix . "_newHairpinScores.txt";
    $parameters->{hairpinTrainFile} = $filePrefix . "_hairpinTrainFile.txt";
    $parameters->{hairpinTrainFileSB} = $filePrefix . "_hairpinTrainFileSB.txt";
    $parameters->{predHairpinClasses} = $filePrefix . "_predHairpinClasses.txt";
    $parameters->{newPredHairpinClasses} = $filePrefix . "_newPredHairpinClasses.txt";
    $parameters->{librarySizes} = $filePrefix . "_librarySizes.txt";
    return $parameters;
}

sub readConfigFile {
    my($configFile,$parameters) = @_;
    open(CONFIGFILE,$configFile) or die "FAIL: could not open $configFile\n";
    while(<CONFIGFILE>) {
	$_ =~ s/^\s+//;
	chomp;
	unless ( $_ eq "" ) {
	    my(@lineData) = split(/\s+/);
	    unless ( @lineData == 3 && $lineData[1] eq "=") {
		die "Error: could not read following line:\n$_\n\nConfig file should be of the form:\nParameter\t=\tValue\n";
	    }
	    my($key,$eqsign,$value) = @lineData;
	    $parameters->{$key} = $value;
	}
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

############################
# Bucket Processing Tools  #
############################

sub readBedFileIntoBuckets {
    my($buckets,$bedFile,$bucketSize) = @_;
    open(BED,$bedFile) or die "failed to open $bedFile\n";
    while (<BED>) {
	chomp;
	my $line = $_;
	unless ( /^#/ ) {
#	    print "$line\n";
	    my($chrom,$start,$stop,$id,$score,$strand) = split(/\t/,$line);
	    my $bucketNum = int($start/$bucketSize);
	    push(@{$buckets->{$chrom}{$bucketNum}}, [$start,$stop,$id,$score,$strand]);
#	    print "$chrom\t$bucketNum\t$start\t$stop\t$id\t$score\t$strand\n";
	}
    }
    close(BED);
    return $buckets;
}

sub getBucketData {
    my($buckets,$location,$bucketSize) = @_;
    my($chrom,$start,$stop) = parseLocation($location);
    my $startBucket = int($start/$bucketSize);
    my $stopBucket = int($stop/$bucketSize);
    my @bucketEntries;
    for (my $bucket = $startBucket; $bucket <= $stopBucket; $bucket++) {
	if ($buckets->{$chrom}{$bucket}) {
	    push(@bucketEntries,@{$buckets->{$chrom}{$bucket}});
	}
    }
    return \@bucketEntries;
}

############################
# FOLD PROCESSING TOOLS    #
############################

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
    my($bamFile) = @_;
    my $bam = Bio::DB::Bam->open( $bamFile );
    return $bam;
}

sub loadBamIndex {
    my($bamFile) = @_;
    my $reIndex;  #changed to 1 if the index file doesn't exist
    my $bamIndex =  Bio::DB::Bam->index($bamFile,$reIndex);
    die "failed to load index for $bamFile\n" if ($reIndex);
    return $bamIndex;
}

sub loadBamList {
    my($bamListFile) = @_;
    open(BLF,$bamListFile) or die "could not open $bamListFile\n";
    my @bamList;
    while(<BLF>) {
	$_ =~ s/^\s+//;
	chomp;
	unless ( $_ eq "" ) {
	    my(@lineData) = split(/\s+/);
	    unless ( @lineData == 2 ) {
		die "Error: could not read following line:\n$_\n\nBam file should be of the form:\nSample\tBamFile.bam\n";
	    }
	    my($label,$bamFile) = @lineData;
	    my $bam = loadBamFile($bamFile);
	    my $bamHeader = $bam->header;
	    my $bamIndex = loadBamIndex($bamFile);
	    my $totalMapped = getBamTotal($bamFile);
	    push(@bamList,[$bam,$bamHeader,$bamIndex,$bamFile,$label,$totalMapped]);
	}
    }
    return \@bamList;
}

sub getBamTotal {
    my($bamFile) = @_;
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
    my($chrom,$start,$stop,$strand) = parseLocation($location);
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
	    my $rStrand = mapBamStrand($alignment->strand);
	    my $count;
	    if ($id =~ /.*_x(\d+)$/) {
		($count) = $id =~ /.*_x(\d+)$/;
	    } else {
		$count = 1;
	    }
	    $seq = reverseComplement($seq) if ($strand eq '-');
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
    return(\%distinctReads,\%uniqueReadCount,\%readCount,\%adjustedReadCount);
}

#sub getTagValue {
#    my($tags,$tagName,$tagType) = @_;
#    foreach my $tag (@{$tags}) {
#	my($tagValue) = $tag =~ /^$tagName:$tagType:(.*)/;
#	if ($tagValue) {
#	    return $tagValue;
#	}
#    }
#    return 0;
#}

sub getDistinctReadCounts {
    my($location,$bamList);
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
    my($bamList,$chromLengths,$repeatRegions,$parameters) = @_;
    my $maxLength = $parameters->{prrMaxLength} or die "Error: prrMaxLength parameter not set\n";
    my $readRegionsFile = $parameters->{readRegionsFile} or die "Error: readRegionsFile parameter not set\n";
    my $longReadRegionsFile = $parameters->{longReadRegionsFile} or die "Error: longReadRegionsFile parameter not set\n";
    my $allReadRegionsFile = $parameters->{allReadRegionsFile} or die "Error: allReadRegionsFile parameter not set\n";
    my $maxMismatches = $parameters->{bamMaxMismatces} or die "Error: bamMaxMisMatches parameter not set\n";
    print "setting bamMaxMismatces for reads in bamfile to $maxMismatches\n" unless($maxMismatches == ~0);
    my $seedLength = $parameters->{bamMismatchSeedLength} or die "Error: bamMismatchSeedLength parameter not set\n";
    print "warning: bamMaxMismatces set but not bamMismatchSeedLength\n" if ($maxMismatches < ~0 && $seedLength == ~0);
    my $initWindowLength = $parameters->{prrWindowLength} or die "Error: prrWindowLength parameter not set\n";
    my $shortReadsCount = 1;     # for regions shorter than 100bp.
    my $longReadsCount = 1;     # for regions longer than 100bp.
    my $repeatReadsCount = 1;     # for repeat regions.
#    my $mu = Memory::Usage->new();
    open(RF,">".$readRegionsFile) or die "Failed to load $readRegionsFile for writing\n";
    open(LRF,">".$longReadRegionsFile) or die "Failed to load $longReadRegionsFile for writing\n";
    open(ARF,">".$allReadRegionsFile) or die "Failed to load $allReadRegionsFile for writing\n";
    foreach my $chrom (keys %{$chromLengths}) {
	my $start = 0;
	while($start <= $chromLengths->{$chrom}) {
	    # push out the windowLength until no reads extend beyond the ends.
	    my $windowLength = $initWindowLength <= $chromLengths->{$chrom} ? $initWindowLength : $chromLengths->{$chrom};
	    # now the window is defined so that all contiguous read regions are properly contained within it.
	    my($plusArray,$minusArray) = getReadCountArrays($bamList,$chrom,$start,$start+$windowLength,$maxMismatches,$seedLength,$parameters);
	    # we now have plus and minus arrays for thie window.
 	    my $posReadRegions = extractReadRegionList($plusArray,$start,$windowLength,"+");
 	    my $negReadRegions = extractReadRegionList($minusArray,$start,$windowLength,"-");
	    # sort readRegions in descending order by rrStart position
	    my @readRegions = sort {$b->[0] <=> $a->[0]} (@{$posReadRegions}, @{$negReadRegions});
	    my $newStart = $start+$windowLength;
#	    $mu->record("readRegions obtained for $chrom:$start-$newStart");
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
#    $mu->dump();
}

sub getReadCountArrays {
    my($bamList,$chrom,$windowStart,$windowStop,$maxMismatches,$seedLength,$parameters) = @_;
    my $minLen = ($parameters->{minReadLength}) ? $parameters->{minReadLength} : 0;
    my $maxLen = ($parameters->{maxReadLength}) ? $parameters->{maxReadLength} : ~0;
    #print "minLen = $minLen\t maxLen = $maxLen\n";
    my $location = "$chrom:$windowStart-$windowStop";
    my(@plusArray,@minusArray);
    foreach my $bamData (@{$bamList}) {
	my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamData};
	my($tid,$chromStart,$chromStop) = $bamHeader->parse_region("$location"); #chromStart and chromStop returned as zero based coordinates
	# callback function to be used with the $bamIndex->fetch() command from the Bio::DB::Sam module
	my $callBack = sub {
	    my($alignment,$data) = @_;
	    my($windowStart,$windowStop,$plusArr,$minusArr) = @{$data};
	    my $rStart = $alignment->start;  #returned in 1 based coordinates
	    my $rStop = $alignment->end;  #returned in 1 based coordinates
	    my $rStrand = mapBamStrand($alignment->strand);
	    my $seq = $alignment->qseq;
	    my $length = length($seq);
	    my $numMismatchesInSeed = 0;
	    if ($maxMismatches < ~0) {
		my $cigar = $alignment->get_tag_values('MD');
		my @cigarArray = ($rStrand eq '-') ? reverse(@{breakCigar($cigar)}) : @{breakCigar($cigar)};
		my $mismatches = getMismatchesFromCigar(\@cigarArray);
		$numMismatchesInSeed = getNumMismatchesInSeedRegion($mismatches,$seedLength);
	    }
	    if ($rStrand eq '-') {
		$seq = reverseComplement($seq);
	    }
	    if ($length >= $minLen and $length <= $maxLen and $numMismatchesInSeed <= $maxMismatches) {
		if($rStart >= $windowStart) {
		    if($rStrand eq "+") {
			# in plus strand
			for(my $i=$rStart;$i<=$rStop;$i++) {
			    $plusArr->[$i-$windowStart]++;     #generating arrays relative to windowStart
			}
		    } else {
			# in minus strand
			for(my $i=$rStart;$i<=$rStop;$i++) {
			    $minusArr->[$i-$windowStart]++;
			}
		    }
		} else {
		    die "Found read at $chrom:$windowStart..$windowStop that extends beyond beginning: $chrom:$rStart..$rStop:$rStrand\n";
		}
	    }  
	};
	my $callBackData = [$windowStart,$windowStop,\@plusArray,\@minusArray];
	my $code = $bamIndex->fetch($bam,$tid,$chromStart,$chromStop,$callBack,$callBackData);
    }
    return(\@plusArray,\@minusArray);
}

sub getNumMismatchesInSeedRegion {
    my($mismatches,$seedLength) = @_;
    my $numMismatches = 0;
    for (my $i = 0; $i < @{$mismatches}; $i++) {
	if ($mismatches->[$i] < $seedLength) {
	    $numMismatches++;
	} else {
	    last;
	}
    }
    return $numMismatches;
}

sub breakCigar {
    #breakCigar returns an array of mismatches found in the cigar of the bamfile
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
#sub processReadRegionsMultiThreaded {
#    my($bamList,$readRegionsFile,$genomeDir,$chromLengths,$parameters) = @_;
#    my $numThreads = $parameters->{numThreads};
#    my $readRegionsFiles = breakupReadRegionsFile($readRegionsFile,$numThreads);
#    my @threads;
#    #creating threads
#    for (my $i = 0; $i < $numThreads; $i++) {
#	my $threadReadRegionsFile = $readRegionsFiles->{$i};
#	my($threadRRFileBase) = $threadReadRegionsFile =~ /([^\/]*)\.txt/;
#	my $newParams = {"outputPrefix" => "$threadRRFileBase"};
#	my $thr = threads->create('processReadRegionsWithNewParams',$bamList,$threadReadRegionsFile,$genomeDir,$chromLengths,$parameters,$newParams);
#	push(@threads,$thr);
#    }
#    #joining threads
#    my @running = threads->list(threads::running);
#    while (@running != 0) {
#	foreach my $thr (@threads) {
#	   $thr->join if ($thr->is_joinable()); 
#	}
#	@running = threads->list(threads::running);
#    }
#}

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
    #my($productReads,$totalReads) = (0,0);
    my $totalReads = 0;
    foreach my $productInfo (@{$products}) {
	my($id,$productReads,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$productInfo};
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
	my($id,$productReads,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$productInfo};
	foreach my $drInfo (@{$productReads}) {
	    my($relStart,$relStop,$dStart,$readCount,$hitCount,$libraryCounts,$seq,$adjustedSeqCount) = @{$drInfo};
	    #print "$relStart x $relStop $adjustedSeqCount added\n";
	    $readCountHash{$id}{$relStart} += $readCount;
	    $totalHPReadCount += $readCount;
	    #print "totalHPReadCount now $totalHPReadCount\n";
	}
    }
    #geting variances of reads for each product
    foreach my $productInfo (@{$products}) {
	my($id,$productReads,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$productInfo};
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
	my($id,$productReads,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$productInfo};
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

sub getReadDataProductVariance_old {
    my($products,$parameters) = @_;
    my(%readCountHash,%productVariance,%productWeight);
    my($totalHPReadCount,$averageProductVariance,$weightedAverageProductVariance) = (0,0,0);
    #accumulating reads into a hash keyed by start position
    foreach my $productInfo (@{$products}) {
	my($id,$productReads,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$productInfo};
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
	my($id,$productReads,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$productInfo};
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
	my($id,$productReads,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$productInfo};
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
    my @sortedProducts = sort {$b->[5] <=> $a->[5]} @{$senseProducts};
    my $maxProduct = $sortedProducts[0];
    my($mpSide,$mpType,$mpProductList,$mpProductCount,$mpMaxProductCount,$mpAdjProductCount,$mpAdjMaxProductCount,$mpRelStart,$mpRelStop,$mpOffset,$mpGStart,$mpProductStrand) = @{$maxProduct};  
#    print "$readId\n$fold\n";
#    print " " x $mpRelStart . substr($fold,$mpRelStart,$mpRelStop-$mpRelStart+1) . "\n";
#    print "total length = " . length($fold) . "\n"; 
#    print "$mpRelStart\t$mpRelStop\n";
    my $map = pairMap($fold);
    my($j,$k,$jFound) = (-1,-1,0);
    for(my $i = $mpRelStart; $i <= $mpRelStop; $i++) {
#	print "\$i = $i\n";
	unless ($jFound) {
	    my $duplexPair = (defined(${$map}[$i])) ? ${$map}[$i] : '.';
#	    print "$i\t$duplexPair\n";
	    unless ($duplexPair eq '.') {
		$j = $duplexPair;
		$k = $j;
		$jFound = 1;
#		print "\$j = $j\t\$k = $k\n";
	    }
	} else {
	    #my $duplexPair = ${$map}[$i]; 
	    #sometimes the major product is not part of the readregion and when the sequence is extended in processReadRegions it does not cover the
	    #entire major product.  In the following line '.' is assigned to $duplexPair for parts of the major product not on the hairpin.  In the 
	    #future it may be better to extend the sequence further in order to fully cover the major product.
	    my $duplexPair = (defined(${$map}[$i])) ? ${$map}[$i] : '.';
#	    print "$i\t$duplexPair\n";
	    unless ($duplexPair eq '.') {
		$k = $duplexPair;
#		print "$i\t$k\n";
	    }
	}
    }
#    print "\n";
    my $duplexProduct = 0;
    my $duplexProductCount = 0;
    my $duplexProductOverlap = 0;
    unless ($j == -1) {  #j=-1 if major product was in a loop and there was no duplex
	my $duplexLeft = ($j < $k) ? $j : $k;
	my $duplexRight = ($k <= $j) ? $j : $k;
	#print "duplex of $mpRelStart $mpRelStop ($mpType)  mapped to $duplexLeft $duplexRight\n";
	foreach my $product (@sortedProducts) {
	    my($side,$type,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};  
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
    my($side,$type,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$maxProduct};
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
	my($dupSide,$dupType,$dupProductList,$dupProductCount,$dupMaxProductCount,$dupAdjProductCount,$dupAdjMaxProductCount,$dupRelStart,$dupRelStop,$dupOffset,$dupGStart,$dupProductStrand) = @{$duplexProduct};	
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
    my($mpSide,$mpType,$mpProductList,$mpProductCount,$mpMaxProductCount,$mpAdjProductCount,$mpAdjMaxProductCount,$mpRelStart,$mpRelStop,$mpOffset,$mpGStart,$mpProductStrand) = @{$maxProduct};
    my($dupSide,$dupType,$dupProductList,$dupProductCount,$dupMaxProductCount,$dupAdjProductCount,$dupAdjMaxProductCount,$dupRelStart,$dupRelStop,$dupOffset,$dupGStart,$dupProductStrand);
    if ($duplexProduct) {
	($dupSide,$dupType,$dupProductList,$dupProductCount,$dupMaxProductCount,$dupAdjProductCount,$dupAdjMaxProductCount,$dupRelStart,$dupRelStop,$dupOffset,$dupGStart,$dupProductStrand) = @{$duplexProduct};
    }
    my($mpOverlapScore,$dupOverlapScore) = (0,0);
    my(@mpOverlaps,@dupOverlaps);
    my($totalMPOverlap,$totalDupOverlap) = (0,0);
    foreach my $product (@{$senseProducts}) {
	my($side,$type,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
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
	my($side,$type,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
	my $overlap = getOverlap($relStart,$relStop,$mpRelStart,$mpRelStop);
	$wghtMPOverlapIntensity = $wghtMPOverlapIntensity + $adjProductCount * $overlap / $totalMPOverlap;
    }
    my $wghtDupOverlapIntensity = 0;
    #my $numDupOverlaps = @dupOverlaps;
    foreach my $product (@dupOverlaps) {
	my($side,$type,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
	my $overlap = getOverlap($relStart,$relStop,$dupRelStart,$dupRelStop);
	$wghtDupOverlapIntensity = $wghtDupOverlapIntensity + $adjProductCount * $overlap / $totalDupOverlap;
    }
    return($mpOverlapScore,$dupOverlapScore,$wghtMPOverlapIntensity,$wghtDupOverlapIntensity);
}

sub correctMorNames {
    #a temporary function used to correct a bug with the moR naming.  This function may be replaced eventually.
    my($products,$fold,$tagName,$parameters) = @_;
    my $minMorDistance = 4;
    my $miRProducts = retrieveMiRProducts($products);
    #setting mir products into seperate variables
    my($miR3p,$miR5p) = (0,0);
    foreach my $miR (@{$miRProducts}) {
	my($side,$type) = @{$miR};
	if ($side eq '5p') {
	    $miR5p = $miR;
	} elsif ($side eq '3p') {
	    $miR3p = $miR;
	} else {
	    die "unidentified side in correctMorNames() side = $side";
	}
	unless ($type =~ /miR/) {
	    die "retrieveMiRProducts returned a $type";
	}
    }
    #checking for to see if moRs exist
    my($hasMoR5p,$hasMoR3p) = (0,0);
    foreach my $product (@{$products}) {
	my($side,$type) = @{$product};
	if ($type =~ /moR/) {
	    unless ($type =~ /long/) {
		if ($side eq '5p') {
#		    print "has 5p mor\n";
		    $hasMoR5p = 1;
		} elsif ($side eq '3p') {
#		    print "has 3p mor\n";
		    $hasMoR3p = 1;
		}
	    } else {
#		print "\n$tagName $side $type changed:\n";
#		printFoldAndProducts($products,$fold,$tagName);
		$product->[1] = "long-out";
#		print "new...\n";
#		printFoldAndProducts($products,$fold,$tagName);
	    }
	}
    }
    my @sortedProducts = sort {$a->[7] <=> $b->[7]} @{$products};
    my @newProducts;
    foreach my $product (@sortedProducts) {
	my($side,$type,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
	if ($type eq "out" && $side eq '5p' && not $hasMoR5p) {
	    if ($miR5p) {
		my($miRSide,$miRType,$miRProductList,$miRProductCount,$miRMaxProductCount,$miRAdjProductCount,$miRAdjMaxProductCount,$miRRelStart,$miRRelStop,$miROffset,$miRGStart,$miRProductStrand) = @{$miR5p};
		if (($miRRelStart - $relStop) <= $minMorDistance && $relStop < $miRRelStop) {
#		    print "\n$tagName $side $type changed:\n";
#		    printFoldAndProducts($products,$fold,$tagName);
		    $product->[1] = "moR";
		    $hasMoR5p = 1;
#		    print "new...\n";
#		    printFoldAndProducts($products,$fold,$tagName);
		}
	    }
	} elsif ($type eq "out" && $side eq '3p' && not $hasMoR3p) {
	    if ($miR3p) {
		my($miRSide,$miRType,$miRProductList,$miRProductCount,$miRMaxProductCount,$miRAdjProductCount,$miRAdjMaxProductCount,$miRRelStart,$miRRelStop,$miROffset,$miRGStart,$miRProductStrand) = @{$miR3p};
		if (($relStart - $miRRelStop) <= $minMorDistance && $relStart > $miRRelStart) {
#		    print "\n$tagName $side $type changed:\n";
#		    printFoldAndProducts($products,$fold,$tagName);
		    $product->[1] = "moR";
		    $hasMoR3p = 1;
#		    print "new...\n";
#		    printFoldAndProducts($products,$fold,$tagName);
		}
	    }
	} elsif ($type eq "moR" && $side eq '5p') {
	    if ($miR5p) {
		my($miRSide,$miRType,$miRProductList,$miRProductCount,$miRMaxProductCount,$miRAdjProductCount,$miRAdjMaxProductCount,$miRRelStart,$miRRelStop,$miROffset,$miRGStart,$miRProductStrand) = @{$miR5p};
		if (($miRRelStart - $relStop) > $minMorDistance) {
#		    print "\n$tagName $side $type changed:\n";
#		    printFoldAndProducts($products,$fold,$tagName);
		    $product->[1] = "out";
#		    print "new...\n";
#		    printFoldAndProducts($products,$fold,$tagName);
		}
	    }
	} elsif ($type eq "moR" && $side eq '3p') {
	    if ($miR3p) {
		my($miRSide,$miRType,$miRProductList,$miRProductCount,$miRMaxProductCount,$miRAdjProductCount,$miRAdjMaxProductCount,$miRRelStart,$miRRelStop,$miROffset,$miRGStart,$miRProductStrand) = @{$miR3p};
		if (($relStart - $miRRelStop) > $minMorDistance) {
#		    print "\n$tagName $side $type changed:\n";
#		    printFoldAndProducts($products,$fold,$tagName);
		    $product->[1] = "out";
#		    print "new...\n";
#		    printFoldAndProducts($products,$fold,$tagName);
		}
	    }
	}
	push(@newProducts,$product);
    }
    my @newSortedProducts = sort {$b->[5] <=> $a->[5]} @newProducts;
    return \@newSortedProducts;
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

sub getFoldRelStartStop {
    my($fold) = @_;
    my($middle,$left) = $fold =~ /^((\+*).*[^\+])/;
    my $relStart = length($left);
    my $relStop = length($middle);
    return($relStart,$relStop);
}

sub getFoldGenomicLocation {
    my($location,$fold) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my($relStart,$relStop) = getFoldRelStartStop($fold);
    my($foldStart,$foldStop);
    if ($strand eq '+') {
	$foldStart = $start + $relStart;
	$foldStop = $start + $relStop;
    } else {
	$foldStart = $stop - $relStop;
	$foldStop = $stop - $relStart;
    }
    my $foldLocation = "$chrom:$foldStart..$foldStop:$strand";
    return $foldLocation;
}

sub processReadRegions {
    my($bamList,$predProductRR,$genomeDir,$chromLengths,$librarySizes,$parameters) = @_;
    my $testTime = 0;
    my $testMemory = 0;
    my $testReadCounts = 0;
    my $readsCount = 0; #number of mirBase hairpins with reads
    my $noReadsCount = 0;  #number of mirBase hairpins with no reads
    my $readsLessThanCountMinLocus = 0; #number of hairpins with fewer than count min locus reads
    my $centersReadsCount = 0; #number of centers on mirBase hairpins with reads
    my $noCenterCount = 0;  #number of hairpins without centers (hopefully this will stay 0)
    my $noMiRCount = 0;  #number of hairpins without any miR's being labeled
    #my $mu = Memory::Usage->new() if ($testMemory);
    my @sampleList;
    my $minLength = $parameters->{lengthMin} or die "minLength: not loaded.\n";
    my $filteredCandidatesFile = $parameters->{filteredCandidatesFile};
    my $hairpinsFile = $parameters->{hairpinsFile};
    my $distinctReadsFile = $parameters->{distinctReadsFile};
    my $productFile = $parameters->{productFile};
    my $featureVectorFile = $parameters->{featureVectorFile};
    my $maxPredictedProductScores = getMaxPredictedProductScores($parameters->{predProductClasses});
    #loading annotation data
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
    print HL "\tmpCount\tdupCount\tdupOverlap\tmpLoopDistance\tdupLoopDistance\tloopSize\tmpOverlapScore\tdupOverlapScore\twghtMPOverlapIntensity\twghtDupOverlapIntensity\tOPA2\ttotalOverlapAmount\taverageOverlapAmount\ttotalRelativeOverlapAmount\taverageRelativeOverlapAmount\tmaxOverlap\tinnerLoopGapCount\tsenseRPM\tantiSenseRPM\toverlapSizeRatioSum\ttotalOverlap";
    print HL "\tmaxBulge\tmaxInteriorLoop\tintLoopSideDiff\tmaxUnboundOverhang\tnumOffshoots\tdupSize";
    print HL "\trel5pOutCount\trel5pMorCount\trel5pMirCount\trelLoopCount\trelSplitCount\trel3pMirCount\trel3pMorCount\trel3pOutCount\textraOutCount";
    print HL "\tmiRmoR5pOverlap\tmiRmoR3pOverlap\tmiRLoop5pOverlap\tmiRLoop3pOverlap\tloopLoopOverlap\tout5pOverlap\tout3pOverlap\toutOut5pOverlap\toutOut3pOverlap\tinProdOverlap\tmiRSplit5pOverlap\tmiRSplit3pOverlap\tmoRSplit5pOverlap\tmoRSplit3pOverlap\tsplitLoopOverlap\tloopSplitOverlap";
    if ($maxPredictedProductScores) {
	print HL "\tRFProductAvg\n";
    } else {
	print HL "\n";
    }
    print HDR "#tag\tstart\tstop\tstrand\ttotal reads\toffsetStart\trelative Start\trelative Stop\tsequence";
    print HPL "#tag\tside type\ttype\ttotal reads\ttotal most abundant\tadjusted reads\tadj total most abundant\tstart\tstop\tstrand\tsequence";
    foreach my $bamListElement (@{$bamList}) {
	my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamListElement};
	print HDR "\t$sample";
        print HPL "\t$sample";
    }
    print HDR "\n";
    print HPL "\n";
    open(RRF, $predProductRR) or die "failed to open $predProductRR\n";
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
	    ($chrom,$start,$stop,$strand) = parseLocation($location);
	    my $sequence = getSequence($location,$genome);
	    my($distinctReads,$uniqueReadCount,$totalDistinctReads,$adjustedReadCount) = retrieveReadData($location,$bamList,$parameters);
	    if ($distinctReads->{$strand}) {  #neccessary if testing read reagions with known mirs and no real read exists in that read region
		my $COUNT= 0;  #count of the number of folds for each read region.
		if ($adjustedReadCount->{$strand} < $parameters->{countMinLocus}) {
		    $readsLessThanCountMinLocus++;
		    print "$id has fewer reads than count min locus\n" if ($testReadCounts);
		}
		$readsCount++;
		my($senseProducts) = buildProducts($location,$distinctReads->{$strand},$strand,$parameters);
		my($averageProductVariance,$weightedAverageProductVariance) = getReadDataProductVariance($senseProducts,$parameters);
		my($averageReadVariance,$weightedAverageReadVariance) = getReadLocationVariance($senseProducts,$parameters);
		my $likelyFolds = getFold($chrom,$sequence,$senseProducts,$id,$parameters);

		foreach my $foldInfo (@$likelyFolds) {
		    my($fold,$mfe,$maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure,$duplexEnergy) = @{$foldInfo};
		    my($foldDuplexCmp,$exactFoldDuplexCmp) = compareFoldAndDuplexStructs($fold,$maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure,$id);
		    my($prodDupPBP_old,$mpDuplexLoopLength_old) = getProductDuplexPBP($fold,$maxProductCoords);
		    my($maxBulge_old,$maxInteriorLoop_old,$intLoopSideDiff_old,$maxUnboundOverhang_old,$numOffshoots_old,$dupSize_old) = getMaxProductFoldInfo($fold,$maxProductCoords);
		    my @centers = getMergedHairpinCenters($fold,$parameters);
		    $centersReadsCount += scalar(@centers);
		    unless (scalar(@centers)) {
			$noCenterCount++;
			print "$id has no centers\n" if ($testReadCounts);
			print "$fold\n" if ($testReadCounts);
		    }
		    foreach my $center (@centers) {
			my $asciiVal = ord('a') + $COUNT;
			my $label = chr($asciiVal);
			my $newId = $id . $label;
			$COUNT++;
			my $basePairs = getBasePairs($center,$fold);
			my $hairpinLength = getMaxHairpinLength($fold,$center);
			my($hStart,$hStop) = getFullHairpinWindow($fold,$center);
			if($hairpinLength >= $minLength) {
			    my $senseProducts= getProductInfo($center,$fold,$basePairs,$location,$senseProducts,$parameters);
			    $senseProducts= rebuildProducts($center,$fold,$basePairs,$senseProducts,$parameters);
			    $senseProducts= addProductSize($senseProducts,$parameters);

			    #a function being used temporarily to correct a bug in the naming of mors
			    #later versions will have the naming function corrected
			    $senseProducts = correctMorNames($senseProducts,$fold,$newId);

			    my($viableHairpin,$reason) = checkHairpinViability($senseProducts,$parameters);
			    if ($viableHairpin) {
				my($rel5pOutCount,$rel5pMorCount,$rel5pMirCount,$relLoopCount,$relSplitCount,$rel3pMirCount,$rel3pMorCount,$rel3pOutCount,$extraOutCount) = getRelativeCounts($senseProducts,$center);
				my($prodDupPBP,$mpDuplexLoopLength) = getMiRProductDuplexPBP($fold,$senseProducts);
				my($maxBulge,$maxInteriorLoop,$intLoopSideDiff,$maxUnboundOverhang,$numOffshoots,$dupSize) = getMiRProductFoldInfo($fold,$senseProducts);
				my $revStrand = revStrand($strand);
				my($antiSenseProducts) = buildProducts($location,$distinctReads->{$revStrand},$revStrand,$parameters);
				$antiSenseProducts = getProductInfo($center,$fold,$basePairs,$location,$antiSenseProducts,$parameters);
				$antiSenseProducts = rebuildProducts($center,$fold,$basePairs,$antiSenseProducts,$parameters);
				$antiSenseProducts= addProductSize($antiSenseProducts,$parameters);
				my $adjTotalProductReads = getAdjTotalProductReads($senseProducts);
				my $adjTotalRevProductReads = getAdjTotalProductReads($antiSenseProducts);
				my($PASS,$REASON) = plausibleReads($adjTotalProductReads,$adjTotalRevProductReads,
									     $senseProducts,$antiSenseProducts,$parameters);
				my($maxProduct,$mpDuplex,$maxProductCount,
				   $mpDuplexCount,$duplexOverlap) = getMaxAndDuplexProducts($senseProducts,$fold,$id);
				my($mpLoopDistance_old,$dupLoopDistance_old,$loopSize_old) = getMaxProductLoopDistance($center,$maxProduct,$mpDuplex);
				my($mpOverlapScore_old,$dupOverlapScore_old,
				   $wghtMPOverlapIntensity_old,$wghtDupOverlapIntensity_old) = getMaxProductOverlapScores($senseProducts,$maxProduct,
															    $mpDuplex);
				my($mpLoopDistance,$dupLoopDistance,$loopSize) = getMiRProductLoopDistance($center,$senseProducts);
				my($mpOverlapScore,$dupOverlapScore,
				   $wghtMPOverlapIntensity,$wghtDupOverlapIntensity) = getMiRProductOverlapScores($senseProducts);
				my ($miRmoR5pOverlap,$miRmoR3pOverlap,
				    $miRLoop5pOverlap,$miRLoop3pOverlap,
				    $loopLoopOverlap,$out5pOverlap,$out3pOverlap,
				    $outOut5pOverlap,$outOut3pOverlap,$inProdOverlap,
				    $miRSplit5pOverlap,$miRSplit3pOverlap,$moRSplit5pOverlap,
				    $moRSplit3pOverlap,$splitLoopOverlap,$loopSplitOverlap) = getProductOverlaps($senseProducts,$fold,$newId,$parameters);
				my($tpd,$totalRP) = newGetReverseProductDisplacement($senseProducts,
											       $antiSenseProducts,
											       $adjTotalProductReads,
											       $parameters,$newId);
				my $apd = $totalRP ? $tpd/$totalRP : 0.0;
				my $urf = $adjustedReadCount->{$strand} ? $uniqueReadCount->{$strand} / $adjustedReadCount->{$strand} : 0;
				my $ahc = computeMaxProdHitCount($senseProducts,$location,$distinctReads,$parameters);
				my $afh = computeMaxProdFivePrimeHet($senseProducts,$parameters);
				my $pbp = computeProductBasePairing($center,$senseProducts,$basePairs,$parameters);
				my $sameShift = computeMaxSameShift($location,$distinctReads->{$strand},$senseProducts,$parameters);
				my $bothShift = computeMaxBothShift($basePairs,$location,$distinctReads->{$strand},
									      $senseProducts,$parameters,$newId);
#			my($zScore,$consensusMFE) = getRNAzData($sequence);
				my $innerLoopGapCount = getInnerLoopGapCount($center,$fold);
				my $splitProductAbundance = computeSplitProductAbundance($senseProducts,$parameters);
				my $overlapProductAbundance = computeOverlapProductAbundance($senseProducts,$parameters);
				my $opa2 = computeOverlapProductAbundance2($senseProducts,$parameters);
				my ($totalOverlapAmount,$averageOverlapAmount,
				    $totalRelativeOverlapAmount,$averageRelativeOverlapAmount) = computeOverlapAmounts($senseProducts,$parameters);
				my($overlapSizeRatioSum,$totalOverlap) = getOverlapRatio($senseProducts,$parameters,$overlapProductAbundance);
				my $maxOverlap = getMaximumOverlap($senseProducts,$parameters);
				my($aa,$ac,$ag,$at,$ca,$cc,$cg,$ct,$ga,$gc,$gg,$gt,$ta,$tc,$tg,$tt) = getDinucleotideFrequencies($sequence);
				my $sequenceLength = length($sequence);
				my $gcContent = getGCcontent($sequence,$sequenceLength);
				# total is the read count for the whole hairpin
				my $totalSense = 0;
				foreach my $product (@{$senseProducts}) {
				    my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
				       $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
				    my $length = $relStop - $relStart + 1;
				    my $productSequence = substr($sequence,$relStart,$length);
				    $totalSense += $adjProductCount;
				    my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($productList);
				    print HPL "$newId\t$side\t$newType\t$productCount\t$maxProductCount\t$adjProductCount\t";
				    print HPL "$adjMaxProductCount\t$relStart\t$relStop\t$productStrand\t";
				    print HPL "$productSequence";
				    foreach my $bamElement (@{$bamList}) {
					my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
					printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
				    }
				    print HPL "\n";
				    foreach my $read (@{$productList}) {
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
				    my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
				       $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
				    $newType = "as-".$newType;
				    $totalAntisense += $adjProductCount;
				    my $length = $relStop - $relStart + 1;
				    my $productSequence = reverseComplement(substr($sequence,$relStart,$length));
				    my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($productList);
				    print HPL "$newId\t$side\t$newType\t$productCount\t$maxProductCount\t$adjProductCount\t";
				    print HPL "$adjMaxProductCount\t$relStart\t$relStop\t$productStrand\t";
				    print HPL "$productSequence";
				    foreach my $bamElement (@{$bamList}) {
					my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
					printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
				    }
				    print HPL "\n";
				    foreach my $read (@{$productList}) {
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
				       $maxProductCount,$mpDuplexCount,$duplexOverlap,$mpLoopDistance,$dupLoopDistance,
				       $loopSize,$mpOverlapScore,$dupOverlapScore,$wghtMPOverlapIntensity,$wghtDupOverlapIntensity);
				printf(HL "\t%.3f",
				       $opa2);
				printf(HL "\t%.3f\t%.3f\t%.9f\t%.3f", $totalOverlapAmount,$averageOverlapAmount,
				       $totalRelativeOverlapAmount,$averageRelativeOverlapAmount);
				printf(HL "\t%.3f", $maxOverlap);
				printf(HL "\t%.3f",  $innerLoopGapCount);
				my $totalSenseRPM = getHairpinRPM($totalSense,$librarySizes);
				my $totalAntisenseRPM = getHairpinRPM($totalAntisense,$librarySizes);
				printf(HL "\t%.3f\t%.3f\t%.9f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f", $totalSenseRPM,
				       $totalAntisenseRPM,$overlapSizeRatioSum,$totalOverlap,$maxBulge,
				       $maxInteriorLoop,$intLoopSideDiff,$maxUnboundOverhang,$numOffshoots,$dupSize);
				printf(HL "\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f",$rel5pOutCount,$rel5pMorCount,
				       $rel5pMirCount,$relLoopCount,$relSplitCount,$rel3pMirCount,$rel3pMorCount,$rel3pOutCount,$extraOutCount);
				print HL "\t$miRmoR5pOverlap\t$miRmoR3pOverlap\t$miRLoop5pOverlap\t$miRLoop3pOverlap\t$loopLoopOverlap\t$out5pOverlap\t$out3pOverlap";
				print HL "\t$outOut5pOverlap\t$outOut3pOverlap\t$inProdOverlap\t$miRSplit5pOverlap\t$miRSplit3pOverlap\t$moRSplit5pOverlap\t$moRSplit3pOverlap\t$splitLoopOverlap\t$loopSplitOverlap";
				if ($maxProductScore) {
				    printf(HL "\t%.3f\n", $maxProductScore);
				} else {
				    print HL "\n";
				}
				#############################
				#############################
			    } else {
#				print "Dropped $newId: $reason\n";
				$noMiRCount++;
			    }
			} 
		    }
		}
	    } else {
		$noReadsCount++;
		print "no reads found at $id\t$location\n" if ($testReadCounts);
	    }
	}
    }
    #$mu->record("after finishing") if ($testMemory);
    print "hairpins with reads = $readsCount\n" if ($testReadCounts);
    print "hairpins without reads = $noReadsCount\n" if ($testReadCounts);
    print "hairpin centers = $centersReadsCount\n" if ($testReadCounts);
    print "hairpins without centers (not including those with no reads) = $noCenterCount\n" if ($testReadCounts);
    print "hairpins with reads less than countMinLocus (not including those with no reads) = $readsLessThanCountMinLocus\n" if ($testReadCounts);
    print "hairpins without any miRs being labeled = $noMiRCount (out of $centersReadsCount)\n";

    close(RRF);
    close(HPL);
    close(HDR);
    close(HL);
    #$mu->dump() if ($testMemory);
}

sub processReadRegions_old {
    my($bamList,$predProductRR,$genomeDir,$chromLengths,$librarySizes,$parameters) = @_;
    my $testTime = 0;
    my $testMemory = 0;
    my $testReadCounts = 0;
    my $readsCount = 0; #number of mirBase hairpins with reads
    my $noReadsCount = 0;  #number of mirBase hairpins with no reads
    my $readsLessThanCountMinLocus = 0; #number of hairpins with fewer than count min locus reads
    my $centersReadsCount = 0; #number of centers on mirBase hairpins with reads
    my $noCenterCount = 0;  #number of hairpins without centers (hopefully this will stay 0)
    #my $mu = Memory::Usage->new() if ($testMemory);
    my @sampleList;
    my $minLength = $parameters->{lengthMin} or die "minLength: not loaded.\n";
    my $filteredCandidatesFile = $parameters->{filteredCandidatesFile};
    my $hairpinsFile = $parameters->{hairpinsFile};
    my $distinctReadsFile = $parameters->{distinctReadsFile};
    my $productFile = $parameters->{productFile};
    my $featureVectorFile = $parameters->{featureVectorFile};
    my $maxPredictedProductScores = getMaxPredictedProductScores($parameters->{predProductClasses});
    #loading annotation data
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
    print HL "\tmpCount\tdupCount\tdupOverlap\tmpLoopDistance\tdupLoopDistance\tloopSize\tmpOverlapScore\tdupOverlapScore\twghtMPOverlapIntensity\twghtDupOverlapIntensity\tOPA2\ttotalOverlapAmount\taverageOverlapAmount\ttotalRelativeOverlapAmount\taverageRelativeOverlapAmount\tmaxOverlap\tinnerLoopGapCount\tsenseRPM\tantiSenseRPM\toverlapSizeRatioSum\ttotalOverlap";
    print HL "\tmaxBulge\tmaxInteriorLoop\tintLoopSideDiff\tmaxUnboundOverhang\tnumOffshoots\tdupSize";
    if ($maxPredictedProductScores) {
	print HL "\tRFProductAvg\n";
    } else {
	print HL "\n";
    }
    print HDR "#tag\tstart\tstop\tstrand\ttotal reads\toffsetStart\trelative Start\trelative Stop\tsequence";
    print HPL "#tag\tside type\ttype\ttotal reads\ttotal most abundant\tadjusted reads\tadj total most abundant\tstart\tstop\tstrand\tsequence";
    foreach my $bamListElement (@{$bamList}) {
	my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamListElement};
	print HDR "\t$sample";
        print HPL "\t$sample";
    }
    print HDR "\n";
    print HPL "\n";
    open(RRF, $predProductRR) or die "failed to open $predProductRR\n";
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
	    ($chrom,$start,$stop,$strand) = parseLocation($location);
	    my $sequence = getSequence($location,$genome);
	    my($distinctReads,$uniqueReadCount,$totalDistinctReads,$adjustedReadCount) = retrieveReadData($location,$bamList,$parameters);
	    if ($distinctReads->{$strand}) {  #neccessary if testing read reagions with known mirs and no real read exists in that read region
		my $COUNT= 0;  #count of the number of folds for each read region.
		if ($adjustedReadCount->{$strand} < $parameters->{countMinLocus}) {
		    $readsLessThanCountMinLocus++;
		    print "$id has fewer reads than count min locus\n" if ($testReadCounts);
		}
		$readsCount++;
		my($senseProducts) = buildProducts($location,$distinctReads->{$strand},$strand,$parameters);
		my($averageProductVariance,$weightedAverageProductVariance) = getReadDataProductVariance($senseProducts,$parameters);
		my($averageReadVariance,$weightedAverageReadVariance) = getReadLocationVariance($senseProducts,$parameters);
		my $likelyFolds = getFold($chrom,$sequence,$senseProducts,$id,$parameters);
#		if (@{$likelyFolds} > 1) {
#		    print "$id:\n";
#		    foreach my $foldInfo (@{$likelyFolds}) {
#			my($fold,$mfe,$maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure,$duplexEnergy) = @{$foldInfo};
#			print "$fold\n";
#		    }
#		    print "\n";
#		}
		foreach my $foldInfo (@$likelyFolds) {
		    my($fold,$mfe,$maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure,$duplexEnergy) = @{$foldInfo};
		    my($foldDuplexCmp,$exactFoldDuplexCmp) = compareFoldAndDuplexStructs($fold,$maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure,$id);
		    my($prodDupPBP,$mpDuplexLoopLength) = getProductDuplexPBP($fold,$maxProductCoords);
		    my($maxBulge,$maxInteriorLoop,$intLoopSideDiff,$maxUnboundOverhang,$numOffshoots,$dupSize) = getMaxProductFoldInfo($fold,$maxProductCoords);
		    my @centers = getMergedHairpinCenters($fold,$parameters);
		    $centersReadsCount += scalar(@centers);
		    unless (scalar(@centers)) {
			$noCenterCount++;
			print "$id has no centers\n" if ($testReadCounts);
			print "$fold\n" if ($testReadCounts);
		    }
		    foreach my $center (@centers) {
			my $asciiVal = ord('a') + $COUNT;
			my $label = chr($asciiVal);
			my $newId = $id . $label;
			$COUNT++;
			my $basePairs = getBasePairs($center,$fold);
			my $hairpinLength = getMaxHairpinLength($fold,$center);
			my($hStart,$hStop) = getFullHairpinWindow($fold,$center);
			if($hairpinLength >= $minLength) {
			    my $senseProducts=getProductInfo($center,$fold,$basePairs,$location,$senseProducts,$parameters);
			    $senseProducts=rebuildProducts($center,$fold,$basePairs,$senseProducts,$parameters);
			    $senseProducts=addProductSize($senseProducts,$parameters);
			    my $revStrand = revStrand($strand);
			    my($antiSenseProducts) = buildProducts($location,$distinctReads->{$revStrand},$revStrand,$parameters);
			    $antiSenseProducts = getProductInfo($center,$fold,$basePairs,$location,$antiSenseProducts,$parameters);
			    $antiSenseProducts = rebuildProducts($center,$fold,$basePairs,$antiSenseProducts,$parameters);
			    $antiSenseProducts=addProductSize($antiSenseProducts,$parameters);
			    my $adjTotalProductReads = getAdjTotalProductReads($senseProducts);
			    my $adjTotalRevProductReads = getAdjTotalProductReads($antiSenseProducts);
			    my($PASS,$REASON) = plausibleReads($adjTotalProductReads,$adjTotalRevProductReads,$senseProducts,$antiSenseProducts,
							       $parameters);
			    my($maxProduct,$mpDuplex,$maxProductCount,$mpDuplexCount,$duplexOverlap) = getMaxAndDuplexProducts($senseProducts,$fold,$id);
			    my($mpLoopDistance,$dupLoopDistance,$loopSize) = getMaxProductLoopDistance($center,$maxProduct,$mpDuplex);
			    my($mpOverlapScore,$dupOverlapScore,$wghtMPOverlapIntensity,$wghtDupOverlapIntensity) = getMaxProductOverlapScores($senseProducts,$maxProduct,$mpDuplex);
			    my($tpd,$totalRP)=newGetReverseProductDisplacement($senseProducts,
									       $antiSenseProducts,
									       $adjTotalProductReads,
									       $parameters,$newId);
			    my $apd = $totalRP ? $tpd/$totalRP : 0.0;
			    my $urf = $adjustedReadCount->{$strand} ? $uniqueReadCount->{$strand} / $adjustedReadCount->{$strand} : 0;
			    my $ahc = computeMaxProdHitCount($senseProducts,$location,$distinctReads,$parameters);
			    my $afh = computeMaxProdFivePrimeHet($senseProducts,$parameters);
			    my $pbp = computeProductBasePairing($center,$senseProducts,$basePairs,$parameters);
			    my $sameShift = computeMaxSameShift($location,$distinctReads->{$strand},$senseProducts,$parameters);
			    my $bothShift = computeMaxBothShift($basePairs,$location,$distinctReads->{$strand},$senseProducts,$parameters,$newId);
#			my($zScore,$consensusMFE) = getRNAzData($sequence);
			    my $innerLoopGapCount = getInnerLoopGapCount($center,$fold);
			    my $splitProductAbundance = computeSplitProductAbundance($senseProducts,$parameters);
			    my $overlapProductAbundance = computeOverlapProductAbundance($senseProducts,$parameters);
			    my $opa2 = computeOverlapProductAbundance2($senseProducts,$parameters);
			    my ($totalOverlapAmount,$averageOverlapAmount,
				$totalRelativeOverlapAmount,$averageRelativeOverlapAmount) = computeOverlapAmounts($senseProducts,$parameters);
			    my($overlapSizeRatioSum,$totalOverlap) = getOverlapRatio($senseProducts,$parameters,$overlapProductAbundance);
			    my $maxOverlap = getMaximumOverlap($senseProducts,$parameters);
			    my($aa,$ac,$ag,$at,$ca,$cc,$cg,$ct,$ga,$gc,$gg,$gt,$ta,$tc,$tg,$tt) = getDinucleotideFrequencies($sequence);
			    my $sequenceLength = length($sequence);
			    my $gcContent = getGCcontent($sequence,$sequenceLength);
			    # total is the read count for the whole hairpin
			    my $totalSense = 0;
			    foreach my $product (@{$senseProducts}) {
				my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
				   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
				my $length = $relStop - $relStart + 1;
				my $productSequence = substr($sequence,$relStart,$length);
				$totalSense += $adjProductCount;
				my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($productList);
				print HPL "$newId\t$side\t$newType\t$productCount\t$maxProductCount\t$adjProductCount\t";
				print HPL "$adjMaxProductCount\t$relStart\t$relStop\t$productStrand\t";
				print HPL "$productSequence";
				foreach my $bamElement (@{$bamList}) {
				    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
				    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
				}
				print HPL "\n";
				foreach my $read (@{$productList}) {
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
				my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
				   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
				$newType = "as-".$newType;
				$totalAntisense += $adjProductCount;
				my $length = $relStop - $relStart + 1;
				my $productSequence = reverseComplement(substr($sequence,$relStart,$length));
				my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($productList);
				print HPL "$newId\t$side\t$newType\t$productCount\t$maxProductCount\t$adjProductCount\t";
				print HPL "$adjMaxProductCount\t$relStart\t$relStop\t$productStrand\t";
				print HPL "$productSequence";
				foreach my $bamElement (@{$bamList}) {
				    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
				    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
				}
				print HPL "\n";
				foreach my $read (@{$productList}) {
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
			    printf(HL "\t%.3f\t%.3f\t%.9f\t%.3f", $totalOverlapAmount,$averageOverlapAmount,$totalRelativeOverlapAmount,$averageRelativeOverlapAmount);
			    printf(HL "\t%.3f", $maxOverlap);
			    printf(HL "\t%.3f",  $innerLoopGapCount);
			    my $totalSenseRPM = getHairpinRPM($totalSense,$librarySizes);
			    my $totalAntisenseRPM = getHairpinRPM($totalAntisense,$librarySizes);
			    printf(HL "\t%.3f\t%.3f\t%.9f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f", $totalSenseRPM,$totalAntisenseRPM,$overlapSizeRatioSum,
				   $totalOverlap,$maxBulge,$maxInteriorLoop,$intLoopSideDiff,$maxUnboundOverhang,$numOffshoots,$dupSize);
			    if ($maxProductScore) {
				printf(HL "\t%.3f\n", $maxProductScore);
			    } else {
				print HL "\n";
			    }
			    #############################
			    #############################
			} 
		    }
		}
	    } else {
		$noReadsCount++;
		print "no reads found at $id\t$location\n" if ($testReadCounts);
	    }
	}
    }
    #$mu->record("after finishing") if ($testMemory);

    print "hairpins with reads = $readsCount\n" if ($testReadCounts);
    print "hairpins without reads = $noReadsCount\n" if ($testReadCounts);
    print "hairpin centers = $centersReadsCount\n" if ($testReadCounts);
    print "hairpins without centers (not including those with no reads) = $noCenterCount\n" if ($testReadCounts);
    print "hairpins with reads less than countMinLocus (not including those with no reads) = $readsLessThanCountMinLocus\n" if ($testReadCounts);

    close(RRF);
    close(HPL);
    close(HDR);
    close(HL);
    close(HPA);
    close(PA);
    #$mu->dump() if ($testMemory);
}

sub mirPreprocess {
    my($bamList,$hairpins,$genomeDir,$chromLengths,$librarySizes,$parameters) = @_;
    my @sampleList;
    my $foldFull = $parameters->{foldFullSequence};
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
    print HL "aapd\ttapd\turf\tahc\tmpfh\tpbp\tsameShifted\tbothShift\tsplitProductAbundance\toverlapProductAbundance\t";
    print HL "aa\tac\tag\tat\tca\tcc\tcg\tct\tga\tgc\tgg\tgt\tta\ttc\ttg\ttt\tduplexEnergy\tRFDV\n";
    print HDR "#tag\tstart\tstop\tstrand\ttotal reads\toffsetStart\trelative Start\trelative Stop\tsequence";
    print HPL "#tag\tside type\ttype\ttotal reads\ttotal most abundant\tadjusted reads\tadj total most abundant\tstart\tstop\tstrand\tsequence";
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
#		my($fold,$mfe,$maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure,$duplexEnergy) = getFoldFromBestDuplex($chrom,$sequence,$senseProducts,$hairpinId);
		my($fold,$mfe) = ($foldFull) ? RNA::fold($sequence) : getFoldFromBestDuplex($chrom,$sequence,$senseProducts,$id);
		my @centers = getMergedHairpinCenters($fold,$parameters);
		my($center) = @centers;
		my $basePairs;
		if ($center) {
		    $basePairs = getBasePairs($center,$fold);
		    my $hairpinLength = getMaxHairpinLength($fold,$center);
		    my($hStart,$hStop) = getFullHairpinWindow($fold,$center);
		    $senseProducts=getProductInfo($center,$fold,$basePairs,$location,$senseProducts,$parameters);
		    $senseProducts=rebuildProducts($center,$fold,$basePairs,$senseProducts,$parameters);
		    my $revStrand = revStrand($strand);
		    my($antiSenseProducts) = buildProducts($location,$distinctReads->{$revStrand},$revStrand,$parameters);
		    $antiSenseProducts = getProductInfo($center,$fold,$basePairs,$location,$antiSenseProducts,$parameters);
		    $antiSenseProducts = rebuildProducts($center,$fold,$basePairs,$antiSenseProducts,$parameters);
		    my $urf = $adjustedReadCount->{$strand} ? $uniqueReadCount->{$strand} / $adjustedReadCount->{$strand} : 0;
		    # total is the read count for the whole hairpin
		    #####################################
		    #   sense Products
		    my $totalSense = 0;
		    foreach my $product (@{$senseProducts}) {
			my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
			   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
			my $length = $relStop - $relStart + 1;
			my $productSequence = substr($sequence,$relStart,$length);
			$totalSense += $adjProductCount;
			my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($productList);
			print HPL "$mirName\t$side\t$newType\t$productCount\t$maxProductCount\t$adjProductCount\t";
			print HPL "$adjMaxProductCount\t$relStart\t$relStop\t$productStrand\t";
			print HPL "$productSequence";
			foreach my $bamElement (@{$bamList}) {
			    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
			    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
			}
			print HPL "\n";
			foreach my $read (@{$productList}) {
			    my($relStart,$relStop,$offset,$gStart,
			       $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
			    my $gStop = $gStart + ($relStop - $relStart);
			    my $adjCount = $count * $adjSeqCount;
			    print HDR "$mirName\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
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
			my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
			   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
			$newType = "as-".$newType;
			$totalAntisense += $adjProductCount;
			my $length = $relStop - $relStart + 1;
			my $productSequence = reverseComplement(substr($sequence,$relStart,$length));
			my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($productList);
			print HPL "$mirName\t$side\t$newType\t$productCount\t$maxProductCount\t$adjProductCount\t";
			print HPL "$adjMaxProductCount\t$relStart\t$relStop\t$productStrand\t";
			print HPL "$productSequence";
			foreach my $bamElement (@{$bamList}) {
			    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
			    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
			}
			print HPL "\n";
			foreach my $read (@{$productList}) {
			    my($relStart,$relStop,$offset,$gStart,
			       $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
			    my $gStop = $gStart + ($relStop - $relStart);
			    my $adjCount = $count * $adjSeqCount;
			print HDR "$mirName\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
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
		    print HL "$mirName\t$chrom\t$start\t$stop\t$strand\t$leftCenter\t$rightCenter\t";
		    printf(HL "%.3f\t%.3f\t%.3f\t",$totalSense,$totalAntisense,$mfe);
		    print HL "$sequence\t$fold\t";   
		    printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\t",
			   0,0,$urf,0,0,0,0,0,0,0);
		    printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
			   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1);
		} else {
		    my $revStrand = revStrand($strand);
		    my($antiSenseProducts) = buildProducts($location,$distinctReads->{$revStrand},$revStrand,$parameters);
		    #####################################
		    #   sense Products
		    my $totalSense = 0;
		    #in this case the products aren't named because no center was found
		    my $offset = 0;
		    my $side = 'loop';
		    my $newType = 'loop';
		    foreach my $product (@{$senseProducts}) {
			#these are the variable of the product returned from buildProducts
			#as it hasn't gone through the full processing
			my($id,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
			   $relStart,$relStop,$gStart,$productStrand) = @{$product};
			my $length = $relStop - $relStart + 1;
			my $productSequence = substr($sequence,$relStart,$length);
			$totalSense += $adjProductCount;
			my $adjTotalLibraryCounts = getAdjTotalLibraryCounts_noCenter($productList);
			print HPL "$mirName\t$side\t$newType\t$productCount\t$maxProductCount\t$adjProductCount\t";
			print HPL "$adjMaxProductCount\t$relStart\t$relStop\t$productStrand\t";
			print HPL "$productSequence";
			foreach my $bamElement (@{$bamList}) {
			    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
			    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
			}
			print HPL "\n";
			foreach my $read (@{$productList}) {
			    my($relStart,$relStop,$gStart,
			       $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
			    my $gStop = $gStart + ($relStop - $relStart);
			    my $adjCount = $count * $adjSeqCount;
			    print HDR "$mirName\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
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
			#these are the variable of the product returned from buildProducts
			#as it hasn't gone through the full processing
			my($id,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
			   $relStart,$relStop,$gStart,$productStrand) = @{$product};
			$newType = "as-".$newType;
			$totalAntisense += $adjProductCount;
			my $length = $relStop - $relStart + 1;
			my $productSequence = reverseComplement(substr($sequence,$relStart,$length));
			my $adjTotalLibraryCounts = getAdjTotalLibraryCounts_noCenter($productList);
			print HPL "$mirName\t$side\t$newType\t$productCount\t$maxProductCount\t$adjProductCount\t";
			print HPL "$adjMaxProductCount\t$relStart\t$relStop\t$productStrand\t";
			print HPL "$productSequence";
			foreach my $bamElement (@{$bamList}) {
			    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
			    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
			}
			print HPL "\n";
			foreach my $read (@{$productList}) {
			    my($relStart,$relStop,$gStart,
			       $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
			    my $gStop = $gStart + ($relStop - $relStart);
			    my $adjCount = $count * $adjSeqCount;
			print HDR "$mirName\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
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
		    print "no center was found for $mirName\: all products will be named as loops\n";
		    my($leftCenter, $rightCenter) = (0,0);
		    print HL "$mirName\t$chrom\t$start\t$stop\t$strand\t$leftCenter\t$rightCenter\t";
		    printf(HL "%.3f\t%.3f\t%.3f\t",$totalSense,$totalAntisense,$mfe);
		    print HL "$sequence\t$fold\t";   
		    printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\t",
			   0,0,0,0,0,0,0,0,0,0);
		    printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
			   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1);
		}
	    } else {
		my($fold,$mfe) = RNA::fold($sequence);
		my @centers = getMergedHairpinCenters($fold,$parameters);
		my($center) = @centers;
		my($leftCenter, $rightCenter) = (0,0);
		if ($center) {
		    ($leftCenter, $rightCenter) = @{$center};
		}
		print HL "$mirName\t$chrom\t$start\t$stop\t$strand\t$leftCenter\t$rightCenter\t";
		printf(HL "%.3f\t%.3f\t%.3f\t",0,0,$mfe);
		print HL "$sequence\t$fold\t";   
		printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\t",
		       0,0,0,0,0,0,0,0,0,0);
		printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
		       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1);
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

sub mirPreprocess_old {
    my($bamList,$hairpins,$genomeDir,$chromLengths,$librarySizes,$parameters) = @_;
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
    print HL "aapd\ttapd\turf\tahc\tmpfh\tpbp\tsameShifted\tbothShift\tsplitProductAbundance\toverlapProductAbundance\t";
    print HL "aa\tac\tag\tat\tca\tcc\tcg\tct\tga\tgc\tgg\tgt\tta\ttc\ttg\ttt\tduplexEnergy\tRFDV\n";
    print HDR "#tag\tstart\tstop\tstrand\ttotal reads\toffsetStart\trelative Start\trelative Stop\tsequence";
    print HPL "#tag\tside type\ttype\ttotal reads\ttotal most abundant\tadjusted reads\tadj total most abundant\tstart\tstop\tstrand\tsequence";
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
	    my($fold,$mfe) = RNA::fold($sequence);
	    my @centers = getMergedHairpinCenters($fold,$parameters);
	    my($center) = @centers;
	    my($distinctReads,$uniqueReadCount,$totalDistinctReads,$adjustedReadCount) = retrieveReadData($location,$bamList,$parameters);
	    if ($distinctReads->{$strand}) {  #neccessary if testing read reagions with known mirs and no real read exists in that read region
		$readsCount++;
		my($senseProducts) = buildProducts($location,$distinctReads->{$strand},$strand,$parameters); 
		my $basePairs;
		if ($center) {
		    $basePairs = getBasePairs($center,$fold);
		    my $hairpinLength = getMaxHairpinLength($fold,$center);
		    my($hStart,$hStop) = getFullHairpinWindow($fold,$center);
		    $senseProducts=getProductInfo($center,$fold,$basePairs,$location,$senseProducts,$parameters);
		    $senseProducts=rebuildProducts($center,$fold,$basePairs,$senseProducts,$parameters);
		    my $revStrand = revStrand($strand);
		    my($antiSenseProducts) = buildProducts($location,$distinctReads->{$revStrand},$revStrand,$parameters);
		    $antiSenseProducts = getProductInfo($center,$fold,$basePairs,$location,$antiSenseProducts,$parameters);
		    $antiSenseProducts = rebuildProducts($center,$fold,$basePairs,$antiSenseProducts,$parameters);
		    # total is the read count for the whole hairpin
		    #####################################
		    #   sense Products
		    my $totalSense = 0;
		    foreach my $product (@{$senseProducts}) {
			my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
			   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
			my $length = $relStop - $relStart + 1;
			my $productSequence = substr($sequence,$relStart,$length);
			$totalSense += $adjProductCount;
			my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($productList);
			print HPL "$mirName\t$side\t$newType\t$productCount\t$maxProductCount\t$adjProductCount\t";
			print HPL "$adjMaxProductCount\t$relStart\t$relStop\t$productStrand\t";
			print HPL "$productSequence";
			foreach my $bamElement (@{$bamList}) {
			    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
			    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
			}
			print HPL "\n";
			foreach my $read (@{$productList}) {
			    my($relStart,$relStop,$offset,$gStart,
			       $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
			    my $gStop = $gStart + ($relStop - $relStart);
			    my $adjCount = $count * $adjSeqCount;
			    print HDR "$mirName\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
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
			my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
			   $relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
			$newType = "as-".$newType;
			$totalAntisense += $adjProductCount;
			my $length = $relStop - $relStart + 1;
			my $productSequence = reverseComplement(substr($sequence,$relStart,$length));
			my $adjTotalLibraryCounts = getAdjTotalLibraryCounts($productList);
			print HPL "$mirName\t$side\t$newType\t$productCount\t$maxProductCount\t$adjProductCount\t";
			print HPL "$adjMaxProductCount\t$relStart\t$relStop\t$productStrand\t";
			print HPL "$productSequence";
			foreach my $bamElement (@{$bamList}) {
			    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
			    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
			}
			print HPL "\n";
			foreach my $read (@{$productList}) {
			    my($relStart,$relStop,$offset,$gStart,
			       $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
			    my $gStop = $gStart + ($relStop - $relStart);
			    my $adjCount = $count * $adjSeqCount;
			print HDR "$mirName\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
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
		    print HL "$mirName\t$chrom\t$start\t$stop\t$strand\t$leftCenter\t$rightCenter\t";
		    printf(HL "%.3f\t%.3f\t%.3f\t",$totalSense,$totalAntisense,$mfe);
		    print HL "$sequence\t$fold\t";   
		    printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\t",
			   0,0,0,0,0,0,0,0,0,0);
		    printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
			   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1);
		} else {
		    my $revStrand = revStrand($strand);
		    my($antiSenseProducts) = buildProducts($location,$distinctReads->{$revStrand},$revStrand,$parameters);
		    #####################################
		    #   sense Products
		    my $totalSense = 0;
		    #in this case the products aren't named because no center was found
		    my $offset = 0;
		    my $side = 'loop';
		    my $newType = 'loop';
		    foreach my $product (@{$senseProducts}) {
			#these are the variable of the product returned from buildProducts
			#as it hasn't gone through the full processing
			my($id,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
			   $relStart,$relStop,$gStart,$productStrand) = @{$product};
			my $length = $relStop - $relStart + 1;
			my $productSequence = substr($sequence,$relStart,$length);
			$totalSense += $adjProductCount;
			my $adjTotalLibraryCounts = getAdjTotalLibraryCounts_noCenter($productList);
			print HPL "$mirName\t$side\t$newType\t$productCount\t$maxProductCount\t$adjProductCount\t";
			print HPL "$adjMaxProductCount\t$relStart\t$relStop\t$productStrand\t";
			print HPL "$productSequence";
			foreach my $bamElement (@{$bamList}) {
			    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
			    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
			}
			print HPL "\n";
			foreach my $read (@{$productList}) {
			    my($relStart,$relStop,$gStart,
			       $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
			    my $gStop = $gStart + ($relStop - $relStart);
			    my $adjCount = $count * $adjSeqCount;
			    print HDR "$mirName\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
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
			#these are the variable of the product returned from buildProducts
			#as it hasn't gone through the full processing
			my($id,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
			   $relStart,$relStop,$gStart,$productStrand) = @{$product};
			$newType = "as-".$newType;
			$totalAntisense += $adjProductCount;
			my $length = $relStop - $relStart + 1;
			my $productSequence = reverseComplement(substr($sequence,$relStart,$length));
			my $adjTotalLibraryCounts = getAdjTotalLibraryCounts_noCenter($productList);
			print HPL "$mirName\t$side\t$newType\t$productCount\t$maxProductCount\t$adjProductCount\t";
			print HPL "$adjMaxProductCount\t$relStart\t$relStop\t$productStrand\t";
			print HPL "$productSequence";
			foreach my $bamElement (@{$bamList}) {
			    my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamElement};
			    printf(HPL "\t%.3f",$adjTotalLibraryCounts->{$sample});
			}
			print HPL "\n";
			foreach my $read (@{$productList}) {
			    my($relStart,$relStop,$gStart,
			       $count,$hitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts) = @{$read};
			    my $gStop = $gStart + ($relStop - $relStart);
			    my $adjCount = $count * $adjSeqCount;
			print HDR "$mirName\t$gStart\t$gStop\t$productStrand\t$adjCount\t";
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
		    print "no center was found for $mirName\: all products will be named as loops\n";
		    my($leftCenter, $rightCenter) = (0,0);
		    print HL "$mirName\t$chrom\t$start\t$stop\t$strand\t$leftCenter\t$rightCenter\t";
		    printf(HL "%.3f\t%.3f\t%.3f\t",$totalSense,$totalAntisense,$mfe);
		    print HL "$sequence\t$fold\t";   
		    printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\t",
			   0,0,0,0,0,0,0,0,0,0);
		    printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
			   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1);
		}
	    } else {
		my($leftCenter, $rightCenter) = (0,0);
		if ($center) {
		    ($leftCenter, $rightCenter) = @{$center};
		}
		print HL "$mirName\t$chrom\t$start\t$stop\t$strand\t$leftCenter\t$rightCenter\t";
		printf(HL "%.3f\t%.3f\t%.3f\t",0,0,$mfe);
		print HL "$sequence\t$fold\t";   
		printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\t",
		       0,0,0,0,0,0,0,0,0,0);
		printf(HL "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
		       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1);
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

sub checkHairpinViability {
    my($senseProducts,$parameters) = @_;
    my $miRCount = 0;
    foreach my $product (@{$senseProducts}) {
	my($side,$type,$productList,$productCount,$maxProductCount,$adjProductcount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
	#long miRs could cause problems with some features and so are excluded
	if ($type eq "miR" || $type eq "short-miR") {
	    $miRCount++;
	}
	####
	#  more tests could be included here to limit to a specific product count if needed
	###############
    }
    if ($miRCount) {
	return(1,"");
    }
    return(0,"hairpin contains no miRs");
}

sub getRelativeCounts {
    my($products,$center) = @_;
    my($leftCenter,$rightCenter) = @{$center};
    my($out5pCount,$mor5pCount,$mir5pCount,$loopCount,$splitCount,$mir3pCount,$mor3pCount,$out3pCount,$extraOutCount) = (0,0,0,0,0,0,0,0,0);
    my $totalCount = 0;
    foreach my $product (@{$products}) {
	my($side,$type,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
	unless ($type =~ /long/) {
	    $totalCount += $productCount;
	    if ($side eq '5p') {
		if ($type =~ /miR/) {
		    $mir5pCount += $productCount;
		} elsif ($type =~ /moR/) {
		    $mor5pCount += $productCount;
		} elsif ($type =~ /out/) {
		    $out5pCount += $productCount;
		} elsif ($type =~ /loop/) {
		    $loopCount += $productCount;
		} else {
		    die "no condition found for 5p type = $type in getRelativeCounts\n";
		}
	    } elsif ($side eq '3p') {
		if ($type =~ /miR/) {
		    $mir3pCount += $productCount;
		} elsif ($type =~ /moR/) {
		    $mor3pCount += $productCount;
		} elsif ($type =~ /out/) {
		    $out3pCount += $productCount;
		} elsif ($type =~ /loop/) {
		    $loopCount += $productCount;
		} else {
		    die "no condition found for 3p type = $type in getRelativeCounts\n";
		}
	    } elsif ($side eq 'loop') {
		if ($type =~ /loop/) {
		    $loopCount += $productCount;
		} else {
		    die "no condition found for loop type = $type in getRelativeCounts\n";
		}
	    } elsif ($side eq 'out') {
		if ($relStart > $rightCenter) {
		    $out3pCount += $productCount;
		} elsif ($relStop < $leftCenter) {
		    $out5pCount += $productCount;
		} else {
		    $extraOutCount += $productCount;
		}
	    } elsif ($side eq 'split') {
		if ($type =~ /split/) {
		    $splitCount += $productCount;
		} else {
		    die "no condition found for split type = $type in getRelativeCounts\n";
		}
	    } else {
		die "no condition found for side = $side in getRelativeCounts\n";
	    }
	} #else {
#	    print "long product found";
	#}
    }
    if ($totalCount == 0) {
	print "only long products found\n";
	return(0,0,0,0,0,0,0,0,0);
    }
    my $rel5pOutCount = $out5pCount / $totalCount;
    my $rel5pMorCount = $mor5pCount / $totalCount;
    my $rel5pMirCount = $mir5pCount / $totalCount;
    my $relLoopCount = $loopCount / $totalCount;
    my $relSplitCount = $splitCount / $totalCount;
    my $rel3pMirCount = $mir3pCount / $totalCount;
    my $rel3pMorCount = $mor3pCount / $totalCount;
    my $rel3pOutCount = $out3pCount / $totalCount;
    my $relExtraOutCount = $extraOutCount / $totalCount;
    return($rel5pOutCount,$rel5pMorCount,$rel5pMirCount,$relLoopCount,$relSplitCount,$rel3pMirCount,$rel3pMorCount,$rel3pOutCount,$relExtraOutCount);
}

sub getProductOverlaps {
    #distances are negative if products overlap
    #in cases where there are two choices for distance the least distance wins
    my($products,$fold,$tagName,$parameters) = @_;
    my @sortedProducts = sort {$a->[7] <=> $b->[7]} @{$products};
    my($miRmoR5pOverlap,$miRmoR3pOverlap,$miRLoop5pOverlap,$miRLoop3pOverlap,$loopLoopOverlap,$out5pOverlap,$out3pOverlap,$outOut5pOverlap,$outOut3pOverlap,$inProdOverlap) = (-~0,-~0,-~0,-~0,-~0,-~0,-~0,-~0,-~0,-~0);
    my($miRSplit5pOverlap,$miRSplit3pOverlap,$moRSplit5pOverlap,$moRSplit3pOverlap,$splitLoopOverlap,$loopSplitOverlap) = (-~0,-~0,-~0,-~0,-~0,-~0);
    my $startItr = 0;
    my $prevProduct = $sortedProducts[$startItr];
    my($initSide,$initType) = @{$prevProduct};
#    print "\n$tagName:\n";
#    printFoldAndProducts($products,$fold,$tagName);
    while ($initType =~ /long/) {
#	print "start product is long in getProdProdOverlaps\n";
	$startItr++;
	if ($startItr >= @sortedProducts) {
	    last;
	} else {
	    $prevProduct = $sortedProducts[$startItr];
	    ($initSide,$initType) = @{$prevProduct};
	}
    }
    for (my $itr = $startItr + 1; $itr < @sortedProducts; $itr++) {
	my $currProduct = $sortedProducts[$itr];
	my($side,$type,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$prevProduct};
	my($side2,$type2,$productList2,$productCount2,$maxProductCount2,$adjProductCount2,$adjMaxProductCount2,$relStart2,$relStop2,$offset2,$gStart2,$productStrand2) = @{$currProduct};
	unless ($type2 =~ /long/) {
	    my $overlap = -($relStart2 - $relStop - 1);
	    if ($side eq '5p' && $type eq 'miR' && $type2 eq 'split') {
		$miRSplit5pOverlap = $overlap if ($miRSplit5pOverlap < $overlap);
#		print "$side-$type overlaps with $side2-$type2 changing miRSplit5pOverlap to $overlap\n";
	    } elsif ($side eq '5p' && $type eq 'moR' && $type2 eq 'split') {
		$moRSplit5pOverlap = $overlap if ($moRSplit5pOverlap < $overlap);
#		print "$side-$type overlaps with $side2-$type2 changing moRSplit5pOverlap to $overlap\n";
	    } elsif ($type eq 'split' && $type2 eq 'loop') {
		$splitLoopOverlap = $overlap if ($splitLoopOverlap < $overlap);
#		print "$side-$type overlaps with $side2-$type2 changing splitLoop5pOverlap to $overlap\n";
	    } elsif ($type eq 'loop' && $type2 eq 'split') {
		$loopSplitOverlap = $overlap if ($loopSplitOverlap < $overlap);
#		print "$side-$type overlaps with $side2-$type2 changing loopSplit3pOverlap to $overlap\n";
	    } elsif ($type eq 'split' && $side2 eq '3p' && $type2 eq 'miR') {
		$miRSplit3pOverlap = $overlap if ($miRSplit3pOverlap < $overlap);
#		print "$side-$type overlaps with $side2-$type2 changing miRSplit3pOverlap to $overlap\n";
	    } elsif ($type eq 'split' && $side2 eq '3p' && $type2 eq 'moR') {
		$moRSplit3pOverlap = $overlap if ($moRSplit3pOverlap < $overlap);		
#		print "$side-$type overlaps with $side2-$type2 changing moRSplit3pOverlap to $overlap\n";
	    } elsif ($side eq '5p' && $type eq 'out' && $side2 eq '5p' && $type2 eq 'out') {
		$outOut5pOverlap = $overlap if ($outOut5pOverlap < $overlap);
#		print "$side-$type overlaps with $side2-$type2 changing outOut5pOverlap to $overlap\n";
	    } elsif ($side eq '5p' && $type eq 'out') {
		$out5pOverlap = $overlap if ($out5pOverlap < $overlap);
#		print "$side-$type overlaps with $side2-$type2 changing Out5pOverlap to $overlap\n";
	    } elsif ($side eq '5p' && $type eq 'moR' && $side2 eq '5p' && $type2 eq 'miR') {
		$miRmoR5pOverlap = $overlap if ($miRmoR5pOverlap < $overlap);
#		print "$side-$type overlaps with $side2-$type2 changing miRmoR5pOverlap to $overlap\n";
	    } elsif ($side eq '5p' && $type eq 'moR' && $side2 eq '5p' && $type2 eq 'out') {
		$inProdOverlap = $overlap if ($inProdOverlap < $overlap);
#		print "$side-$type overlaps with $side2-$type2 changing inProdOverlap to $overlap\n";
	    } elsif ($side eq '5p' && $type eq 'moR' && $overlap > 0) {
		#overlaps with products that are named out but are inside the main part of the hairpin will be called in product overlaps
		$inProdOverlap = $overlap if ($inProdOverlap < $overlap);
#		print "$side-$type overlaps with $side2-$type2 changing inProdOverlap to $overlap\n";
	    } elsif ($side eq '5p' && $type eq 'miR' && $type2 eq 'loop') {
		$miRLoop5pOverlap = $overlap if ($miRLoop5pOverlap < $overlap);
#		print "$side-$type overlaps with $side2-$type2 changing miRLoop5pOverlap to $overlap\n";
	    } elsif ($side eq '5p' && $type eq 'miR' && $side2 eq '5p' && $type2 eq 'out') {
		$inProdOverlap = $overlap if ($inProdOverlap < $overlap);
#		print "$side-$type overlaps with $side2-$type2 changing inProdOverlap to $overlap\n";
	    } elsif ($side eq '5p' && $type eq 'miR' && $overlap > 0) {
		#overlaps with products that are named out but are inside the main part of the hairpin will be called in product overlaps
		$inProdOverlap = $overlap if ($inProdOverlap < $overlap);
#		print "$side-$type overlaps with $side2-$type2 changing inProdOverlap to $overlap\n";
	    } elsif ($type eq 'loop' && $type2 eq 'loop') {
		$loopLoopOverlap = $overlap if ($loopLoopOverlap < $overlap);
#		print "$side-$type overlaps with $side2-$type2 changing loopLoopOverlap to $overlap\n";
	    } elsif ($type eq 'loop' && $side2 eq '3p' && $type2 eq 'miR') {
		$miRLoop3pOverlap = $overlap if ($miRLoop3pOverlap < $overlap);
#		print "$side-$type overlaps with $side2-$type2 changing miRLoop3pOverlap to $overlap\n";
	    } elsif ($type eq 'loop' && $overlap > 0) {
		#overlaps with products that are named out but are inside the main part of the hairpin will be called in product overlaps
		$inProdOverlap = $overlap if ($inProdOverlap < $overlap);
#		print "$side-$type overlaps with $side2-$type2 changing inProdOverlap to $overlap\n";
	    } elsif ($side eq '3p' && $type eq 'miR' && $side2 eq '3p' && $type2 eq 'moR') {
		$miRmoR3pOverlap = $overlap if ($miRmoR3pOverlap < $overlap);		
#		print "$side-$type overlaps with $side2-$type2 changing miRmoR3pOverlap to $overlap\n";
	    } elsif ($side eq '3p' && $type eq 'miR' && $side2 eq '3p' && $type2 eq 'out') {
		$out3pOverlap = $overlap if ($out3pOverlap < $overlap);		
#		print "$side-$type overlaps with $side2-$type2 changing out3pOverlap to $overlap\n";
	    }  elsif ($side eq '3p' && $type eq 'moR' && $side2 eq '3p' && $type2 eq 'out') {
		$out3pOverlap = $overlap if ($out3pOverlap < $overlap);		
#		print "$side-$type overlaps with $side2-$type2 changing out3pOverlap to $overlap\n";
	    } elsif ($side eq '3p' && $type eq 'out') {
		$outOut3pOverlap = $overlap if ($outOut3pOverlap < $overlap);		
#		print "$side-$type overlaps with $side2-$type2 changing outOut3pOverlap to $overlap\n";
	    } elsif ($side2 eq '3p' && $type2 eq 'out') {
		$out3pOverlap = $overlap if ($out3pOverlap < $overlap);
#		print "$side-$type overlaps with $side2-$type2 changing out3pOverlap to $overlap\n";
	    } elsif (($side eq '5p' && $type eq 'miR' && $side2 eq '3p' && $type2 eq 'miR') ||
		     ($side eq '5p' && $type eq 'miR' && $side2 eq '3p' && $type2 eq 'moR') ||
		     ($side eq '5p' && $type eq 'moR' && $side2 eq '3p' && $type2 eq 'miR') ||
		     ($side eq '5p' && $type eq 'moR' && $side2 eq '3p' && $type2 eq 'moR')) {
		#do nothing
	    } else {
#		print "not sure what to do with $side-$type $side2-$type2\n";
	    }
	    $prevProduct = $currProduct;
	}
	unless ($productStrand eq $productStrand2) {
	    die "product strands not equal in getProdProdOverlaps()\n";
	}
    }
    $miRmoR5pOverlap = 0 if ($miRmoR5pOverlap == -~0);
    $miRmoR3pOverlap = 0 if ($miRmoR3pOverlap == -~0);
    $miRLoop5pOverlap = 0 if ($miRLoop5pOverlap == -~0);
    $miRLoop3pOverlap = 0 if ($miRLoop3pOverlap == -~0);
    $loopLoopOverlap = 0 if ($loopLoopOverlap == -~0);
    $out5pOverlap = 0 if ($out5pOverlap == -~0);
    $out3pOverlap = 0 if ($out3pOverlap == -~0);
    $outOut5pOverlap = 0 if ($outOut5pOverlap == -~0);
    $outOut3pOverlap = 0 if ($outOut3pOverlap == -~0);
    $inProdOverlap = 0 if ($inProdOverlap == -~0);
    $miRSplit5pOverlap = 0 if ($miRSplit5pOverlap == -~0);
    $miRSplit3pOverlap = 0 if ($miRSplit3pOverlap == -~0);
    $moRSplit5pOverlap = 0 if ($moRSplit5pOverlap == -~0);
    $moRSplit3pOverlap = 0 if ($moRSplit3pOverlap == -~0);
    $splitLoopOverlap = 0 if ($splitLoopOverlap == -~0);
    $loopSplitOverlap = 0 if ($loopSplitOverlap == -~0);

#    print "miRmoR5pOverlap  = $miRmoR5pOverlap\n";
#    print "miRmoR3pOverlap  = $miRmoR3pOverlap\n";
#    print "miRLoop5pOverlap  = $miRLoop5pOverlap\n";
#    print "miRLoop3pOverlap  = $miRLoop3pOverlap\n";
#	print "loopLoopOverlap  = $loopLoopOverlap\n";
#	print "out5pOverlap  = $out5pOverlap\n";
#	print "out3pOverlap  = $out3pOverlap\n";
#	print "outOut5pOverlap  = $outOut5pOverlap\n";
#	print "outOut3pOverlap  = $outOut3pOverlap\n";
#	print "inProdOverlap  = $inProdOverlap\n";
#	print "miRSplit5pOverlap = $miRSplit5pOverlap\n";
#	print "miRSplit3pOverlap = $miRSplit3pOverlap\n";
#	print "moRSplit5pOverlap = $moRSplit5pOverlap\n";
#	print "moRSplit3pOverlap = $moRSplit3pOverlap\n";
#	print "splitLoop5pOverlap = $splitLoopOverlap\n";
#    print "loopSplit3pOverlap = $loopSplitOverlap\n";

    return($miRmoR5pOverlap,$miRmoR3pOverlap,$miRLoop5pOverlap,$miRLoop3pOverlap,$loopLoopOverlap,$out5pOverlap,$out3pOverlap,$outOut5pOverlap,$outOut3pOverlap,$inProdOverlap,
	   $miRSplit5pOverlap,$miRSplit3pOverlap,$moRSplit5pOverlap,$moRSplit3pOverlap,$splitLoopOverlap,$loopSplitOverlap);
}

sub retrieveMiRProducts {
    my($products) = @_;
    my @miRProducts;
    foreach my $product (@{$products}) {
	my($side,$type,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
	#long miRs could cause problems with some features and so are excluded
	if ($type eq "miR" || $type eq "short-miR") {
	    push(@miRProducts,[$side,$type,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand]);
	}
#	print "$side-$type: $relStart .. $relStop\tcount = $productCount\tadjProductCount = $adjProductCount\n";
    }
    my @sortedMiRProducts = sort {$b->[5] <=> $a->[5]} @miRProducts;
    return \@sortedMiRProducts;
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

sub checkUnfolded {
    my($fold,$coords) = @_;
    my($start,$stop) = @{$coords};
    #sometimes a product is found on the very edge and goes over the folded region
    if ($stop >= length($fold)) {
	$stop = length($fold) - 1;
    }
#    print "checkUnfolded $start..$stop\n";
    my $coordFold = substr($fold,$start,$stop-$start+1);
#    print "$fold\n";
#    print ' ' x $start . 'X' x ($stop - $start + 1) . "\n";
    my $foldCount = ($coordFold =~ tr/(\(|\))//);
#    print "foldCount = $foldCount\n";
    unless ($foldCount) {
	return 1;
    }
    return 0;
}

sub getDuplexCoords {
    my($map,$maxProductCoords,$offshoots) = @_;
    my @duplexCoords = (-1,-1);
    my $FIRST = 1;
    for (my $i = $maxProductCoords->[0]; $i <= $maxProductCoords->[1]; $i++) {
	unless ($map->[$i] eq "." || $map->[$i] eq "+" || insideOffshoot($map->[$i],$offshoots)) {
	    unless ($FIRST) {
		if ($map->[$i] > $duplexCoords[1]) {
		    $duplexCoords[1] = $map->[$i];
		} elsif ($map->[$i] < $duplexCoords[0]) {
		    $duplexCoords[0] = $map->[$i];
		}
	    } else {
		$FIRST = 0;
		$duplexCoords[0] = $map->[$i];
		$duplexCoords[1] = $map->[$i];
	    }
	}
    }
    return \@duplexCoords;
}

sub subtractOffshoots {
    my($unboundRegion,$offshoots) = @_;
    my $ubLen = $unboundRegion->[1] - $unboundRegion->[0] + 1;
    $ubLen -= 2;  #dropping the ends which are bound
    foreach my $offshoot (@{$offshoots}) {
	if ($offshoot->[0] >= $unboundRegion->[0] && $offshoot->[1] <= $unboundRegion->[1]) {
	    my $offLen = $offshoot->[1] - $offshoot->[0] + 1;
	    $ubLen = $ubLen - $offLen + 2;
	}
    }
    return $ubLen;
}

sub getOffshoots {
    my($map,$productCoords) = @_;
    my @offshoots;
    for (my $i = $productCoords->[0]; $i <= $productCoords->[1]; $i++) {
	unless ($map->[$i] eq "." || $map->[$i] eq "+") {
	    if ($map->[$i] >= $productCoords->[0] && $map->[$i] <= $productCoords->[1]) {
		#the nucleotide bonds within the coordinants of the product and we've reached a hairpin 
		#that's an offshoot of the one being analyzed
		push(@offshoots,[$i,$map->[$i]]);
		$i = $map->[$i];
	    }
	}
    }
    return \@offshoots;
}

sub combineOffshoots {
    my($mpOffshoots,$dupOffshoots) = @_;
    my @offshoots = (@{$mpOffshoots},@{$dupOffshoots});
    my @sortedOffshoots = sort {$a->[0] <=> $b->[0] && $b->[1] <=> $a->[1]} @offshoots;
    my @filteredOffshoots;
    my $lastStart = -1;
    for (my $i = 0; $i < @sortedOffshoots; $i++) {
	if ($sortedOffshoots[$i][0] != $lastStart) {
	    push(@filteredOffshoots,$sortedOffshoots[$i]);
	}
    }
    return \@filteredOffshoots;
}

sub getMaxProductFoldInfo {
    my($fold,$maxProductCoords) = @_;
    my ($maxBulge,$maxInteriorLoop,$intLoopSideDiff,$maxUnboundOverhang,$numOffshoots,$dupSize) = (0,0,0,0,0,0);
    my $map = pairMap($fold);
    my $length = $maxProductCoords->[1] - $maxProductCoords->[0] + 1;
    my $itr = $maxProductCoords->[0];
    if ($maxProductCoords->[1] >= length($fold)) {
	$maxProductCoords->[1] = length($fold) - 1;
    }
    #if the product is unfolded the the maxUnboundOverhang will be the size of the product
    if (checkUnfolded($fold,$maxProductCoords)) {
	return(0,0,0,$length,0,0);
    }
    #finding the unbound overhang at the start of the product
    while ($itr <= $maxProductCoords->[1]) {
	if ($map->[$itr] ne "." && $map->[$itr] ne "+") {
	    if ($map->[$itr] >= $maxProductCoords->[0] && $map->[$itr] <= $maxProductCoords->[1]) {
		#we've reached an offshoot hairpin
		$itr = $map->[$itr];
		$maxUnboundOverhang += 2; 
	    } else  {
		last;
	    }
	} else {
	    $maxUnboundOverhang++;
	}
	$itr++;
    }
    #as a special case loops will return the maxUnbound overhang and have a dupSize of 0
    if ($itr > $maxProductCoords->[1]) {
	return(0,0,0,$maxUnboundOverhang,0,0);
    }
    #determining number of offshoots 
    my $mpOffshoots = getOffshoots($map,$maxProductCoords);
    my $duplexCoords = getDuplexCoords($map,$maxProductCoords,$mpOffshoots);
    my $dupOffshoots = getOffshoots($map,$duplexCoords);
    my $offshoots = combineOffshoots($mpOffshoots,$dupOffshoots);
    $numOffshoots = @{$offshoots};
    #determining duplex size
    $dupSize = $duplexCoords->[1] - $duplexCoords->[0] + 1;
    #determining maxBulge, maxInteriorLoop, and intLoopSideDiff
    my $currUB = 0;
    my $lastBoundPos = $itr;
    my $lastDupPos = $map->[$itr];
    while ($itr <= $maxProductCoords->[1]) {
	if ($map->[$itr] eq "." || $map->[$itr] eq "+") {
	    #we're within an unbound region
	    $currUB++;
	} elsif ($currUB > 0) {
	    #we're now on a nucleotide next to an unbound nucleotide on the major product
	    if (abs($map->[$itr] - $lastDupPos) == 1) {
		#we're next to a bulge on the major product
		if ($currUB > $maxBulge) {
		    $maxBulge = $currUB;
		}
	    } else {
		#we're next to a loop
		my @dupRegion = ($map->[$itr] > $lastDupPos) ? ($lastDupPos,$map->[$itr]) : ($map->[$itr],$lastDupPos);
		my $dupLoopUB = subtractOffshoots(\@dupRegion,$dupOffshoots);
		my $loopSize = $dupLoopUB + $currUB;
		my $sideDiff = abs($dupLoopUB - $currUB);
		if (($loopSize > $maxInteriorLoop) || 
		    ($loopSize == $maxInteriorLoop && $sideDiff > $intLoopSideDiff)) {
		    $maxInteriorLoop = $loopSize;
		    $intLoopSideDiff = $sideDiff;
		}
	    }
	    $currUB = 0;
	    $lastBoundPos = $itr;
	    $lastDupPos = $map->[$itr];
	} elsif ($map->[$itr] >= $maxProductCoords->[0] && $map->[$itr] <= $maxProductCoords->[1]) {
	    #we've reached an offshoot hairpin
	    $itr = $map->[$itr];
	    $currUB += 2;
	} else {
	    #we're at a bound nucleotide
	    if (abs($map->[$itr] - $lastDupPos) != 1) {
		#bulge on the duplex
		my @dupRegion = ($map->[$itr] > $lastDupPos) ? ($lastDupPos,$map->[$itr]) : ($map->[$itr],$lastDupPos); 
		my $bulge = subtractOffshoots(\@dupRegion,$dupOffshoots);
		if ($bulge > $maxBulge) {
		    $maxBulge = $bulge;
		}
	    }
	    $lastBoundPos = $itr;
	    $lastDupPos = $map->[$itr];
	}
	$itr++;
    }
    #checking the overhang at the other side and updating maxUnboundOverhang.
    if ($currUB > $maxUnboundOverhang) {
	$maxUnboundOverhang = $currUB;
    }
#    printMPAndDup($fold,$maxProductCoords);
#    print "maxBulge = $maxBulge\nmaxInteriorLoop = $maxInteriorLoop\nintLoopSideDiff = $intLoopSideDiff\nmaxUnboundOverhand = $maxUnboundOverhang\n";
#    print "numOffshoots = $numOffshoots\ndupSize = $dupSize\n";
    return($maxBulge,$maxInteriorLoop,$intLoopSideDiff,$maxUnboundOverhang,$numOffshoots,$dupSize);
}

sub getProductDuplexPBP {
    my($fold,$maxProductCoords) = @_;
    #note: this is for the duplex created by RNAfold rather than RNAduplex
    my($prodDupPBP,$mpDuplexLoopLength) = (0,0);
    my $map = pairMap($fold);
#    print "maxProductCoords: " . $maxProductCoords->[0] . "-" . $maxProductCoords->[1] . "\n";
    my $prodLength = $maxProductCoords->[1] - $maxProductCoords->[0] + 1;
    #finding duplex from RNA structure
    my $i = $maxProductCoords->[0];
#    print "i = $i\n";
    if (checkUnfolded($fold,$maxProductCoords)) {
	return(0,$prodLength);
    }
    while ($map->[$i] eq "." || $map->[$i] eq "+") {
#	print "$i ---> " . $map->[$i] . "\n";
	if ($i >= $maxProductCoords->[1] && ($map->[$i] eq "." || $map->[$i] eq "+")) {
	    #print "no mapping found returning 0 and $i\n";
	    return(0,$prodLength);
	}
	$i++;
    }
#    print "$fold\n";
#    print ' ' x $maxProductCoords->[0] . 'X' x ($maxProductCoords->[1] - $maxProductCoords->[0] + 1) . "\n";
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
#    print "mpDupLoopLength = $mpDuplexLoopLength\n";
    return($prodDupPBP,$mpDuplexLoopLength);
}

sub getMiRProductDuplexPBP {
    my($fold,$products) = @_;
    my $sortedMiRProducts = retrieveMiRProducts($products);
    my($side,$type,$productList,$productCount,$maxProductCount,$adjProductcount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$sortedMiRProducts->[0]};
    $relStart = 0 if ($relStart < 0);
    my $maxProductCoords = [$relStart,$relStop];
    #note: this is for the duplex created by RNAfold rather than RNAduplex
    my($prodDupPBP,$mpDuplexLoopLength) = (0,0);
    my $map = pairMap($fold);
#    print "maxProductCoords: " . $relStart . "-" . $relStop . "\n";
    my $prodLength = $relStop - $relStart + 1;
    #finding duplex from RNA structure
    my $i = $relStart;
#    print "i = $i\n";
    if (checkUnfolded($fold,$maxProductCoords)) {
	return(0,$prodLength);
    }
    while ($map->[$i] eq "." || $map->[$i] eq "+") {
#	print "$i ---> " . $map->[$i] . "\n";
	if ($i >= $relStop && ($map->[$i] eq "." || $map->[$i] eq "+")) {
	    #print "no mapping found returning 0 and $i\n";
	    return(0,$prodLength);
	}
	$i++;
    }
#    print "$fold\n";
#    print ' ' x $relStart . 'X' x ($relStop - $relStart + 1) . "\n";
    my @prodDuplex = ($map->[$i],$map->[$i]);
    while ($i <= $relStop) {
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
#    print "mpDupLoopLength = $mpDuplexLoopLength\n";
    return($prodDupPBP,$mpDuplexLoopLength);
}

sub getMiRProductFoldInfo {
    my($fold,$products) = @_;
    my $sortedMiRProducts = retrieveMiRProducts($products);
    my($side,$type,$productList,$productCount,$maxProductCount,$adjProductcount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$sortedMiRProducts->[0]};
    $relStart = 0 if ($relStart < 0);
    my $maxProductCoords = [$relStart,$relStop];
    my ($maxBulge,$maxInteriorLoop,$intLoopSideDiff,$maxUnboundOverhang,$numOffshoots,$dupSize) = (0,0,0,0,0,0);
    my $map = pairMap($fold);
    my $length = $maxProductCoords->[1] - $maxProductCoords->[0] + 1;
    my $itr = $maxProductCoords->[0];
    if ($maxProductCoords->[1] >= length($fold)) {
	$maxProductCoords->[1] = length($fold) - 1;
    }
    #if the product is unfolded the the maxUnboundOverhang will be the size of the product
    if (checkUnfolded($fold,$maxProductCoords)) {
	return(0,0,0,$length,0,0);
    }
    #finding the unbound overhang at the start of the product
    while ($itr <= $maxProductCoords->[1]) {
	if ($map->[$itr] ne "." && $map->[$itr] ne "+") {
	    if ($map->[$itr] >= $maxProductCoords->[0] && $map->[$itr] <= $maxProductCoords->[1]) {
		#we've reached an offshoot hairpin
		$itr = $map->[$itr];
		$maxUnboundOverhang += 2; 
	    } else  {
		last;
	    }
	} else {
	    $maxUnboundOverhang++;
	}
	$itr++;
    }
    #as a special case loops will return the maxUnbound overhang and have a dupSize of 0
    if ($itr > $maxProductCoords->[1]) {
	return(0,0,0,$maxUnboundOverhang,0,0);
    }
    #determining number of offshoots 
    my $mpOffshoots = getOffshoots($map,$maxProductCoords);
    my $duplexCoords = getDuplexCoords($map,$maxProductCoords,$mpOffshoots);
    my $dupOffshoots = getOffshoots($map,$duplexCoords);
    my $offshoots = combineOffshoots($mpOffshoots,$dupOffshoots);
    $numOffshoots = @{$offshoots};
    #determining duplex size
    $dupSize = $duplexCoords->[1] - $duplexCoords->[0] + 1;
    #determining maxBulge, maxInteriorLoop, and intLoopSideDiff
    my $currUB = 0;
    my $lastBoundPos = $itr;
    my $lastDupPos = $map->[$itr];
    while ($itr <= $maxProductCoords->[1]) {
	if ($map->[$itr] eq "." || $map->[$itr] eq "+") {
	    #we're within an unbound region
	    $currUB++;
	} elsif ($currUB > 0) {
	    #we're now on a nucleotide next to an unbound nucleotide on the major product
	    if (abs($map->[$itr] - $lastDupPos) == 1) {
		#we're next to a bulge on the major product
		if ($currUB > $maxBulge) {
		    $maxBulge = $currUB;
		}
	    } else {
		#we're next to a loop
		my @dupRegion = ($map->[$itr] > $lastDupPos) ? ($lastDupPos,$map->[$itr]) : ($map->[$itr],$lastDupPos);
		my $dupLoopUB = subtractOffshoots(\@dupRegion,$dupOffshoots);
		my $loopSize = $dupLoopUB + $currUB;
		my $sideDiff = abs($dupLoopUB - $currUB);
		if (($loopSize > $maxInteriorLoop) || 
		    ($loopSize == $maxInteriorLoop && $sideDiff > $intLoopSideDiff)) {
		    $maxInteriorLoop = $loopSize;
		    $intLoopSideDiff = $sideDiff;
		}
	    }
	    $currUB = 0;
	    $lastBoundPos = $itr;
	    $lastDupPos = $map->[$itr];
	} elsif ($map->[$itr] >= $maxProductCoords->[0] && $map->[$itr] <= $maxProductCoords->[1]) {
	    #we've reached an offshoot hairpin
	    $itr = $map->[$itr];
	    $currUB += 2;
	} else {
	    #we're at a bound nucleotide
	    if (abs($map->[$itr] - $lastDupPos) != 1) {
		#bulge on the duplex
		my @dupRegion = ($map->[$itr] > $lastDupPos) ? ($lastDupPos,$map->[$itr]) : ($map->[$itr],$lastDupPos); 
		my $bulge = subtractOffshoots(\@dupRegion,$dupOffshoots);
		if ($bulge > $maxBulge) {
		    $maxBulge = $bulge;
		}
	    }
	    $lastBoundPos = $itr;
	    $lastDupPos = $map->[$itr];
	}
	$itr++;
    }
    #checking the overhang at the other side and updating maxUnboundOverhang.
    if ($currUB > $maxUnboundOverhang) {
	$maxUnboundOverhang = $currUB;
    }
#    printMPAndDup($fold,$maxProductCoords);
#    print "maxBulge = $maxBulge\nmaxInteriorLoop = $maxInteriorLoop\nintLoopSideDiff = $intLoopSideDiff\nmaxUnboundOverhand = $maxUnboundOverhang\n";
#    print "numOffshoots = $numOffshoots\ndupSize = $dupSize\n";
    return($maxBulge,$maxInteriorLoop,$intLoopSideDiff,$maxUnboundOverhang,$numOffshoots,$dupSize);
}

sub getMiRProductLoopDistance {
    my($center,$products) = @_;
    my $sortedMiRProducts = retrieveMiRProducts($products);
    my $maxProduct = $sortedMiRProducts->[0];
    my $duplexProduct = (@{$sortedMiRProducts} == 2) ? $sortedMiRProducts->[1] : 0;
    my($leftCenter,$rightCenter) = @{$center};
    my $loopSize = $rightCenter - $leftCenter + 1;
    my($side,$type,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$maxProduct};
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
	my($dupSide,$dupType,$dupProductList,$dupProductCount,$dupMaxProductCount,$dupAdjProductCount,$dupAdjMaxProductCount,$dupRelStart,$dupRelStop,$dupOffset,$dupGStart,$dupProductStrand) = @{$duplexProduct};	
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

sub getMiRProductOverlapScores {
    my($senseProducts) = @_;
    my $sortedMiRProducts = retrieveMiRProducts($senseProducts);
    my $maxProduct = $sortedMiRProducts->[0];
    my $duplexProduct = (@{$sortedMiRProducts} == 2) ? $sortedMiRProducts->[1] : 0;
    my($mpSide,$mpType,$mpProductList,$mpProductCount,$mpMaxProductCount,$mpAdjProductCount,$mpAdjMaxProductCount,$mpRelStart,$mpRelStop,$mpOffset,$mpGStart,$mpProductStrand) = @{$maxProduct};
    my($dupSide,$dupType,$dupProductList,$dupProductCount,$dupMaxProductCount,$dupAdjProductCount,$dupAdjMaxProductCount,$dupRelStart,$dupRelStop,$dupOffset,$dupGStart,$dupProductStrand);
    if ($duplexProduct) {
	($dupSide,$dupType,$dupProductList,$dupProductCount,$dupMaxProductCount,$dupAdjProductCount,$dupAdjMaxProductCount,$dupRelStart,$dupRelStop,$dupOffset,$dupGStart,$dupProductStrand) = @{$duplexProduct};
    }
    my($mpOverlapScore,$dupOverlapScore) = (0,0);
    my(@mpOverlaps,@dupOverlaps);
    my($totalMPOverlap,$totalDupOverlap) = (0,0);
    foreach my $product (@{$senseProducts}) {
	my($side,$type,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
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
	my($side,$type,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
	my $overlap = getOverlap($relStart,$relStop,$mpRelStart,$mpRelStop);
	$wghtMPOverlapIntensity = $wghtMPOverlapIntensity + $adjProductCount * $overlap / $totalMPOverlap;
    }
    my $wghtDupOverlapIntensity = 0;
    #my $numDupOverlaps = @dupOverlaps;
    foreach my $product (@dupOverlaps) {
	my($side,$type,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
	my $overlap = getOverlap($relStart,$relStop,$dupRelStart,$dupRelStop);
	$wghtDupOverlapIntensity = $wghtDupOverlapIntensity + $adjProductCount * $overlap / $totalDupOverlap;
    }
    return($mpOverlapScore,$dupOverlapScore,$wghtMPOverlapIntensity,$wghtDupOverlapIntensity);
}

sub initializeRNAFold {
    $RNA::noLonelyPairs = 1;
};

sub getFoldFromBestDuplex {
    #this function finds the fold based on the greatest duplex of the major product.
    my($chrom,$sequence,$senseProducts,$hairpinId) = @_;
    my $test = 0;
    my $majorProduct = $senseProducts->[0];  #products must be sorted by abundance at this point
    my($id,$productReads,$productCount,$maxProductCount,$adjustedProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$strand) = @{$majorProduct};
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
    return ($newFold,$mfe,\@maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure,$duplexEnergy);
}

sub getFold {
    my($chrom,$sequence,$senseProducts,$hairpinId,$parameters) = @_;
    my @likelyFolds;
    my $minProdFoldDistance = 5;
    my($bestDupFold,$mfe,$maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure,$duplexEnergy) = getFoldFromBestDuplex($chrom,$sequence,$senseProducts,$hairpinId);
    push(@likelyFolds, [$bestDupFold,$mfe,$maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure,$duplexEnergy]);
    unless ($parameters->{noMultipleFolds}) {
	my($maxProdRelStart,$maxProdRelStop) = @{$maxProductCoords};
	my $maxProductSequence = substr($sequence,$maxProdRelStart,$maxProdRelStop-$maxProdRelStart+1);
	for (my $i = 1; $i < @{$senseProducts}; $i++) {
	    my($id,$productReads,$productCount,$maxProductCount,$adjustedProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$strand) = @{$senseProducts->[$i]};
	    $maxRelStart = 0 if ($maxRelStart < 0);  #sometimes the product starts before the hairpin
	    if ($maxRelStart - $maxProdRelStop >= $minProdFoldDistance) {
		my $productSequence = substr($sequence,$maxRelStart,$maxRelStop-$maxRelStart+1);
		($productDuplexCoords,$seqDuplexCoords,$dupStructure,$duplexEnergy) = getDuplex($maxProductSequence,$productSequence);
		$seqDuplexCoords->[0] += $maxRelStart;
		$seqDuplexCoords->[1] += $maxRelStart;
		my $foldSequence = substr($sequence,$maxProdRelStart,$maxRelStop-$maxProdRelStart+1);
		my($fold,$mfe) = RNA::fold($foldSequence);
		my $leftOutRegion = '+' x $maxProdRelStart;
		my $rightOutRegion = '+' x (length($sequence)-$maxRelStop-1);
		my $newFold = $leftOutRegion . $fold . $rightOutRegion;
		push(@likelyFolds, [$newFold,$mfe,$maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure,$duplexEnergy]) unless ($newFold eq $bestDupFold);
	    } elsif ($maxProdRelStart - $maxRelStop >= $minProdFoldDistance) {
		my $productSequence = substr($sequence,$maxRelStart,$maxRelStop-$maxRelStart+1);
		($productDuplexCoords,$seqDuplexCoords,$dupStructure,$duplexEnergy) = getDuplex($maxProductSequence,$productSequence);
		$seqDuplexCoords->[0] += $maxRelStart;
		$seqDuplexCoords->[1] += $maxRelStart;
		my $foldSequence = substr($sequence,$maxRelStart,$maxProdRelStop-$maxRelStart+1);
		my($fold,$mfe) = RNA::fold($foldSequence);
		my $leftOutRegion = '+' x $maxRelStart;
		my $rightOutRegion = '+' x (length($sequence)-$maxProdRelStop-1);
		my $newFold = $leftOutRegion . $fold . $rightOutRegion;
		push(@likelyFolds, [$newFold,$mfe,$maxProductCoords,$productDuplexCoords,$seqDuplexCoords,$dupStructure,$duplexEnergy]) unless ($newFold eq $bestDupFold);
	    }
	}
    }
    return \@likelyFolds;
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

sub getMirStarShift {
    my($duplexCoords,$duplexSide,$products) = @_;
    my($leftDuplexCoord,$rightDuplexCoord) = @{$duplexCoords};
    if ($duplexSide eq '5p') {
	foreach my $product (@{$products}) {
	    my($id,$productReads,$productCount,$maxProductCount,$adjustedProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$strand) = @{$product};
	    if (getOverlap($maxRelStart,$maxRelStop,$leftDuplexCoord,$rightDuplexCoord)) {
		if ($maxRelStart < $leftDuplexCoord) {
		    return $leftDuplexCoord - $maxRelStart;
		}
	    }
	}
    } elsif ($duplexSide eq '3p') {
	foreach my $product (@{$products}) {
	    my($id,$productReads,$productCount,$maxProductCount,$adjustedProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$strand) = @{$product};
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
    my($id,$productReads,$productCount,$maxProductCount,$adjustedProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$strand) = @{$majorProduct};
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
	    my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop) = @{$product};
	    my($side2,$newType2,$productList2,$productCount2,$maxProductCount2,$adjProductCount2,$adjMaxProductCount2,$relStart2,$relStop2) = @{$product2};
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

sub getOverlapRatio {
    my($productInfo,$parameters,$opa) = @_;
    my($overlapSizeRatioSum,$totalOverlap,$totalLength) = (0,0,0);
    for (my $i = 0; $i < @{$productInfo}; $i++) {
	my $product1 = $productInfo->[$i];
	my($side1,$newType1,$productList1,$productCount1,$maxProductCount1,$adjProductCount1,$adjMaxProductCount1,$relStart1,$relStop1) = @{$product1};
	my $length1 = $relStop1 - $relStart1 + 1;
	$totalLength += $length1;
	unless(($newType1 =~ /long/) || ($newType1 =~ /short/)) {
	    for (my $j = $i + 1; $j < @{$productInfo}; $j++) {
		my $product2 = $productInfo->[$j];
		my($side2,$newType2,$productList2,$productCount2,$maxProductCount2,$adjProductCount2,$adjMaxProductCount2,$relStart2,$relStop2) = @{$product2};
		my $length2 = $relStop2 - $relStart2 + 1;
		unless (($newType2 =~ /long/) || ($newType2 =~ /short/)) {
		    if (my $overlap = getOverlap($relStart1,$relStop1,$relStart2,$relStop2)) {
			my $smallProductCount = ($adjProductCount1 < $adjProductCount2) ? $adjProductCount1 : $adjProductCount2;
			my $bigProductCount = ($adjProductCount1 < $adjProductCount2) ? $adjProductCount2 : $adjProductCount1;
			$overlapSizeRatioSum = $overlapSizeRatioSum + $smallProductCount / $bigProductCount;
			$totalOverlap += $overlap;
		    }
		}
	    }
	}	
    }
    #printProducts($productInfo) if ($totalOverlap);
    #print "totalRelativeOverlap = $totalOverlap\ntotalRelativeOverlapAmount = $overlapSizeRatioSum\n$opa\n" if ($totalOverlap);
    return ($overlapSizeRatioSum,$totalOverlap);
}

sub computeOverlapAmounts {
    my($productInfo,$parameters) = @_;
    my($totalOverlapAmount,$averageOverlapAmount,$totalRelativeOverlapAmount,$averageRelativeOverlapAmount,$productTotal) = (0,0,0,0,0);
    foreach my $product(@{$productInfo}) {
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop) = @{$product};
	$productTotal += $adjProductCount;
    }
    my @sortedProducts = sort {$a->[5] <=> $b->[5]} @{$productInfo};  #sorted from least to most abundant
    for (my $i = 0; $i < @sortedProducts; $i++) {
	my $product1 = $sortedProducts[$i];
	my($side1,$newType1,$productList1,$productCount1,$maxProductCount1,$adjProductCount1,$adjMaxProductCount1,$relStart1,$relStop1) = @{$product1};
	unless(($newType1 =~ /long/) || ($newType1 =~ /short/)) {
	    for (my $j = $i + 1; $j < @sortedProducts; $j++) {
		my $product2 = $sortedProducts[$j];
		my($side2,$newType2,$productList2,$productCount2,$maxProductCount2,$adjProductCount2,$adjMaxProductCount2,$relStart2,$relStop2) = @{$product2};
		unless (($newType2 =~ /long/) || ($newType2 =~ /short/)) {
		    if (my $overlap = getOverlap($relStart1,$relStop1,$relStart2,$relStop2)) {
			my $smallProductCount = ($adjProductCount1 < $adjProductCount2) ? $adjProductCount1 : $adjProductCount2;
			my $bigProductCount = ($adjProductCount1 < $adjProductCount2) ? $adjProductCount2 : $adjProductCount1;
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
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop) = @{$product};
	$productTotal += $adjProductCount;
    }
    my @sortedProducts = sort {$a->[5] <=> $b->[5]} @{$productInfo};  #sorted from least to most abundant
    my $maxOverlapProductAbundance = 0;
    for (my $i = 0; $i < @sortedProducts; $i++) {
	my $product1 = $sortedProducts[$i];
	my($side1,$newType1,$productList1,$productCount1,$maxProductCount1,$adjProductCount1,$adjMaxProductCount1,$relStart1,$relStop1) = @{$product1};
	unless(($newType1 =~ /long/) || ($newType1 =~ /short/)) {
	    for (my $j = $i + 1; $j < @sortedProducts; $j++) {
		my $product2 = $sortedProducts[$j];
		my($side2,$newType2,$productList2,$productCount2,$maxProductCount2,$adjProductCount2,$adjMaxProductCount2,$relStart2,$relStop2) = @{$product2};
		unless (($newType2 =~ /long/) || ($newType2 =~ /short/)) {
		    if (my $overlap = getOverlap($relStart1,$relStop1,$relStart2,$relStop2)) {
			if ($overlap > $minOverlap) {
			    my $overlapProductAbundance = $adjProductCount1 / $productTotal;
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
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop) = @{$product};
	$productTotal += $adjProductCount;
    }
    my @sortedProducts = sort {$a->[5] <=> $b->[5]} @{$productInfo};  #sorted from least to most abundant
    my $overlapProductAbundance = 0;
    for (my $i = 0; $i < @sortedProducts; $i++) {
	my $product1 = $sortedProducts[$i];
	my($side1,$newType1,$productList1,$productCount1,$maxProductCount1,$adjProductCount1,$adjMaxProductCount1,$relStart1,$relStop1) = @{$product1};
	unless(($newType1 =~ /long/) || ($newType1 =~ /short/)) {
	    for (my $j = $i + 1; $j < @sortedProducts; $j++) {
		my $product2 = $sortedProducts[$j];
		my($side2,$newType2,$productList2,$productCount2,$maxProductCount2,$adjProductCount2,$adjMaxProductCount2,$relStart2,$relStop2) = @{$product2};
		unless (($newType2 =~ /long/) || ($newType2 =~ /short/)) {
		    if (my $overlap = getOverlap($relStart1,$relStop1,$relStart2,$relStop2)) {
			if ($overlap > $minOverlap) {
			    $overlapProductAbundance = $overlapProductAbundance + $overlap * $adjProductCount1 * $adjProductCount2;
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
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop) = @{$product};
	$productTotal += $adjProductCount;
    }
    my @sortedProducts = sort {$a->[5] <=> $b->[5]} @{$productInfo};  #sorted from least to most abundant
    my $maxOverlapProductAbundance = 0;
    for (my $i = 0; $i < @sortedProducts; $i++) {
	my $product1 = $sortedProducts[$i];
	my($side1,$newType1,$productList1,$productCount1,$maxProductCount1,$adjProductCount1,$adjMaxProductCount1,$relStart1,$relStop1) = @{$product1};
	for (my $j = $i + 1; $j < @sortedProducts; $j++) {
	    my $product2 = $sortedProducts[$j];
	    my($side2,$newType2,$productList2,$productCount2,$maxProductCount2,$adjProductCount2,$adjMaxProductCount2,$relStart2,$relStop2) = @{$product2};
	    if (getOverlap($relStart1,$relStop1,$relStart2,$relStop2)) {
		my $overlapProductAbundance = $adjProductCount1 / $productTotal;
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
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop) = @{$product};
	unless(($newType =~ /long/) || ($newType =~ /short/)) {
	    $splitProductTotal += $adjProductCount if ($newType eq 'split');
	    $productTotal += $adjProductCount;
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
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop) = @{$product};
	$splitProductTotal += $adjProductCount if ($newType eq 'split');
	$productTotal += $adjProductCount;
    }
    return ($productTotal == 0) ? 0 : $splitProductTotal / $productTotal;
}

sub buildProducts {
    my($location,$distinctReads,$productStrand,$parameters) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $NUMPRODUCTS = 1;
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
	    push(@{$productHash{$NUMPRODUCTS}},$parsedRead);
	    $NUMPRODUCTS++;
	}
    }
    my @productInfo;
    foreach my $id (keys %productHash) {
	my $adjTotal=0;
	my($maxRelStart,$maxRelStop,$maxGStart);
	my $FIRST=1;
	my $productCount = 0;
	my $maxProductCount = 0;
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
	    $productCount += $readCount;
	    $adjProductCount += $readCount * $adjustedSeqCount;
	    if(($relStart == $maxRelStart)&&($relStop == $maxRelStop)) {
		$maxProductCount += $readCount;
		$adjMaxProductCount += $readCount * $adjustedSeqCount;
		$maxGStart = $dStart;
	    }
	}
	push(@productInfo,[$id,$productHash{$id},$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand]);
    }
    my @sortedProducts = sort {$b->[4] <=> $a->[4]} @productInfo;
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
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart) = @{$product};
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
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
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
	 my($side1,$newType1,$productList1,$productCount1,$maxProductCount1,$adjProductCount1,$adjMaxProductCount1,$relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	 unless(($newType1 =~ /long/) || ($newType1 =~ /short/)) {
	     if($adjProductCount1 >= $countThreshold) {	
		 if($side1 eq "5p") {
		     print "adjProdCount1 = $adjProductCount1 >= $countThreshold = countThreshold for $side1-$newType1\n" if ($test); 
		     for(my $j=0;$j<@{$productInfo};$j++) {
			 if($i != $j) {
			     my($side2,$newType2,$productList2,$productCount2,$maxProductCount2,$adjProductCount2,$adjMaxProductCount2,$relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
			     if (computeOverlapRatio($adjProductCount1,$adjProductCount2) > $parameters->{maxOverlapProductAbundance}) {
				 unless(($newType2 =~ /long/) || ($newType2 =~ /short/)) {
				     if($adjProductCount2 >= $countThreshold) {
					 if($side2 eq "3p") {				
					     print "\tadjProdCount2 = $adjProductCount2 >= $countThreshold = countThreshold for $side2-$newType2\n" if ($test);
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
	my($side1,$newType1,$productList1,$productCount1,$maxProductCount1,$adjProductCount1,$adjMaxProductCount1,$relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	unless(($newType1 =~ /long/) || ($newType1 =~ /short/)) {
	    if($adjProductCount1 >= $countThreshold) {	
		if($side1 eq "5p") {
		    print "adjProdCount1 = $adjProductCount1 >= $countThreshold = countThreshold for $side1-$newType1\n" if ($test); 
		    for(my $j=0;$j<@{$productInfo};$j++) {
			if($i != $j) {
			    my($side2,$newType2,$productList2,$productCount2,$maxProductCount2,$adjProductCount2,$adjMaxProductCount2,$relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
			    unless(($newType2 =~ /long/) || ($newType2 =~ /short/)) {
				if($adjProductCount2 >= $countThreshold) {
				    if($side2 eq "3p") {				
					print "\tadjProdCount2 = $adjProductCount2 >= $countThreshold = countThreshold for $side2-$newType2\n" if ($test);
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
	 my($side1,$newType1,$productList1,$productCount1,$maxProductCount1,$adjProductCount1,$adjMaxProductCount1,
	    $relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	 if($adjProductCount1 >= $countThreshold) {	
	     if($side1 eq "5p") {
		 for(my $j=0;$j<@{$productInfo};$j++) {
		     if($i != $j) {
			 my($side2,$newType2,$productList2,$productCount2,$maxProductCount2,$adjProductCount2,$adjMaxProductCount2,
			    $relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
			 if($adjProductCount2 >= $countThreshold) {
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
	my($side1,$newType1,$productList1,$productCount1,$maxProductCount1,$adjProductCount1,$adjMaxProductCount1,$relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	 unless (($newType1 =~ /long/) || ($newType1 =~ /short/)) {
	     #print "$i: $side1 $newType1 $relStart1 $relStop1 $productCount1 vs $countThreshold\n";
	     if($adjProductCount1 >= $countThreshold) {
		 print "relStart1 = $relStart1\trelStop1 = $relStop1\n" if ($test);
		 print "adjProdCount1 = $adjProductCount1 >= $countThreshold = countthreshold\n" if ($test);
		 for(my $j=0;$j<@{$productInfo};$j++) {
		     if($i != $j) {
			 my($side2,$newType2,$productList2,$productCount2,$maxProductCount2,$adjProductCount2,$adjMaxProductCount2,$relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
			 if (computeOverlapRatio($adjProductCount1,$adjProductCount2) > $parameters->{maxOverlapProductAbundance}) {
			     unless (($newType2 =~ /long/) || ($newType2 =~ /short/)) {
				 if($adjProductCount2 >= $countThreshold) {
				     print "relStart2 = $relStart2\trelStop2 = $relStop2\n" if ($test);
				     print "adjProdCount2 = $adjProductCount2 >= $countThreshold = countThreshold\n" if ($test);
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
	my($side1,$newType1,$productList1,$productCount1,$maxProductCount1,$adjProductCount1,$adjMaxProductCount1,$relStart1,$relStop1,$offset1,$gStart1) = @{$productInfo->[$i]};
	#print "$i: $side1 $newType1 $relStart1 $relStop1 $productCount1 vs $countThreshold\n";
	if($adjProductCount1 >= $countThreshold) {	    
	    for(my $j=0;$j<@{$productInfo};$j++) {
		if($i != $j) {
		    my($side2,$newType2,$productList2,$productCount2,$maxProductCount2,$adjProductCount2,$adjMaxProductCount2,
		       $relStart2,$relStop2,$offset2,$gStart2) = @{$productInfo->[$j]};
		     if($adjProductCount2 >= $countThreshold) {
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
    die "relStart1 > relStop1 ($relStart1 > $relStop1) in getOverlap()\n" if ($relStart1 > $relStop1);
    die "relStart2 > relStop2 ($relStart2 > $relStop2) in getOverlap()\n" if ($relStart2 > $relStop2);
    if(($relStart1 <= $relStart2)&&($relStop2 <= $relStop1)) {
        #relstart2..relstop2 is within relstart1..relstop1
        return $relStop2 - $relStart2 + 1;
    } elsif (($relStart2 <= $relStart1)&&($relStop1 <= $relStop2)) {
        #relstart1..relstop1 is within relstart2..relstop2
        return $relStop1 - $relStart1 + 1;
    } elsif(($relStart1 <= $relStart2)&&($relStart2 <= $relStop1)) {
        #relStart2..relstop2 overlaps on right of relstart1..relstop1 
        return $relStop1 - $relStart2 + 1;
    } elsif(($relStart1 <= $relStop2)&&($relStop2 <= $relStop1)) {
        #relstart2..relStop2 overlaps on left of relstart1..relstop1
        return $relStop2 - $relStart1 + 1;
    } elsif(($relStart2 <= $relStart1)&&($relStart1 <= $relStop2)) {
        #relstart1..relstop1 overlaps on right of relstart2..relstop2
        return $relStop2 - $relStart1 + 1;
    } elsif(($relStart2 <= $relStop1)&&($relStop1<=$relStop2)) {
        #relstart1..relstop1 overlaps on left of relstart2..relstop2
        return $relStop2-$relStop1 + 1;
    }
    return 0;
}

#getOverlap will gradually be replaced with getOverlap during the next testing phase.
#sub getOverlap {
#    my($relStart1,$relStop1,$relStart2,$relStop2) = @_;    
#    if(($relStart1 <= $relStart2)&&($relStart2 <= $relStop1)) {
#	# relStart2 within first window
#	return $relStop1 - $relStart2 + 1;
#    } elsif(($relStart1 <= $relStop2)&&($relStop2 <= $relStop1)) {
#	# if relStop2 within first window
#	return $relStop2 - $relStart1 + 1;	
#    } elsif(($relStart2 <= $relStart1)&&($relStart1 <= $relStop2)) {
#	return $relStop2 - $relStart1 + 1;
#    } elsif(($relStart2 <= $relStop1)&&($relStop1<=$relStop2)) {
#	return $relStop2-$relStop1 + 1;
#    }
#    return 0;
#}

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
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
	   $relStart,$relStop,$offset,$gStart) = @{$productInfo->[$i1]};
	my $start = $relStart;
	my $stop = $relStop;
	for(my $i2=0;$i2<@{$revProductInfo};$i2++) {
	    my($aType,$aNewType,$aProdList,$aProductCount,$aMaxProductCount,$aAdjProdCount,$aAdjMaxProdCount,
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
    # make global variable
    my $minFrac = 0.0005;
    my $sum = 0;
    my $total = 0;
    my %USED;
    my %AS_USED;
    my @distances;
    for(my $i1=0;$i1<@{$productInfo};$i1++) {
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart) = @{$productInfo->[$i1]};
	unless (($newType =~ /long/) || ($newType =~ /short/)) {
	    my $start = $relStart;
	    my $stop = $relStop;
	    for(my $i2=0;$i2<@{$revProductInfo};$i2++) {
		my($aType,$aNewType,$aProdList,$aProductCount,$aMaxProductCount,$aAdjProdCount,$aAdjMaxProdCount,
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
			    if ((computeOverlapRatio($adjProductCount,$aAdjProdCount) > $parameters->{maxOverlapProductAbundance}) &&
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
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart) = @{$productInfo->[$i1]};
	unless (($newType =~ /long/) || ($newType =~ /short/)) {
	    print "side=$side\ttype=$newType\tadjProdCount=$adjProductCount\tadjMaxProdCount=$adjMaxProductCount\n" if ($test);
	    print "relStart=$relStart\trelStop=$relStop\toffset=$offset\tgStart=$gStart\n" if ($test);
	    my $start = $relStart;
	    my $stop = $relStop;
	    my $minDisp = $minDispInit;
	    print "minDisp=$minDisp\n" if ($test);
	    my $minI2;
	    for(my $i2=0;$i2<@{$revProductInfo};$i2++) {
		my($aType,$aNewType,$aProdList,$aProdCount,$aMaxProductCount,$aAdjProdCount,$aAdjMaxProdCount,
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
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
	   $relStart,$relStop,$offset,$gStart) = @{$productInfo->[$i1]};
	print "side=$side\ttype=$newType\tadjProdCount=$adjProductCount\tadjMaxProdCount=$adjMaxProductCount\n" if ($test);
	print "relStart=$relStart\trelStop=$relStop\toffset=$offset\tgStart=$gStart\n" if ($test);
	my $start = $relStart;
	my $stop = $relStop;
	my $minDisp = $minDispInit;
	print "minDisp=$minDisp\n" if ($test);
	my $minI2;
	for(my $i2=0;$i2<@{$revProductInfo};$i2++) {
	    my($aType,$aNewType,$aProdList,$aProdCount,$aMaxProductCount,$aAdjProdCount,$aAdjMaxProdCount,
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
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart) = @{$product};
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
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
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
	 my($side,$newType,$productList,$productCount,$maxProductCount) = @{$product};
	 my $fivePrimeHet = computeFivePrimeHet($productList);
	 #print "prodCount: $productCount, 5' het: $fivePrimeHet\n";
	 if(($fivePrimeHet < $maxFivePrimeHet)&&($productCount > 1)) {
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
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount) = @{$product};
	unless (($newType =~ /long/) || ($newType =~ /short/)) {
	    if($newType eq "miR") {
		my $fivePrimeHet = computeFivePrimeHet($productList);
		if($adjProductCount > $maxCount) {
		    $maxProdFivePrimeHet = $fivePrimeHet;
		    $maxCount = $adjMaxProductCount;
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
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount) = @{$product};
	if($newType eq "miR") {
	    my $fivePrimeHet = computeFivePrimeHet($productList);
	    if($adjProductCount > $maxCount) {
		$maxProdFivePrimeHet = $fivePrimeHet;
		$maxCount = $adjMaxProductCount;
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
	my($oSide,$oNewType,$oProductlist,$oProductCount,$oMaxProductCount,$oAdjustedProductCount,$oAdjustedMaxProdCount,$oRelStart,$oRelStop) = @{$oProduct};
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
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart) = @{$product};
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
		if($adjProductCount > $maxCount) {
		    $maxCount = $adjProductCount;
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
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
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
	     if($adjProductCount > $maxCount) {
		 $maxCount = $adjProductCount;
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
	my($id,$productReads,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$product};
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
	push(@productInfo,[$side,\@newProductReadsData,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
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
	my($side,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount) = @{$products};
	$totalReads += $adjProductCount;
    }   
    my @sortedProductInfo = sort {abs($a->[8]) <=> abs($b->[8])} @{$productInfo};
    foreach my $products (@sortedProductInfo) {
	my($side,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
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
	push(@newProductInfo1,[$side,$newType,$productList,$productCount,
			       $maxProductCount,$adjProductCount,$adjMaxProductCount,
			       $relStart,$relStop,$offset,$gStart,$productStrand]);
    }
    foreach my $products (@sortedProductInfo) {
	my($side,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
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
			       $productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
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
	my($side,$type,$productList,$productCount,$maxProductCount,$adjProductcount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
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
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart) = @{$product};
	unless (($newType =~ /long/) || ($newType =~ /short/)) { 
	    $adjTotalReads += $adjProductCount;
	}
    }
    return $adjTotalReads;
}

#this version of getAdjTotalProductReads does not take into account product size
sub oldGetAdjTotalProductReads {
    my $productInfo = shift;
    my $adjTotalReads = 0;
    foreach my $product (@{$productInfo}) {
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,
	   $relStart,$relStop,$offset,$gStart) = @{$product};
	$adjTotalReads += $adjProductCount;
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

sub getAdjTotalLibraryCounts_noCenter {
    #for a special case when the center isn't present in mirPreprocess
    my $productList = shift;
    my %adjTotalLibraryCounts;
    foreach my $read (@{$productList}) {
	my($rRelStart,$rRelStop,$rGStart,$rCount,$rHitCount,$libraryCounts,$seq,$adjSeqCount,$adjLibraryCounts)=@{$read};
	foreach my $sample (keys %{$adjLibraryCounts}) {
	    $adjTotalLibraryCounts{$sample} += $adjLibraryCounts->{$sample};
	}
    }
    return \%adjTotalLibraryCounts;
}


########################################
#  Annotations Functions               #
########################################

sub annotateProducts {
   my($location,$fold,$products,$annotProducts,$hairpinId) = @_;
   my $minOffsetLimit = 2;  #offset cannot be greater than 2 for an annotation;
   my @productAnnotations;
   my %annotProductCandidates;
   my %annotProductOffsets;
   #finding all products that overlap annotations
#   print $hairpinId . "\n";
   my($foldRelStart,$foldRelStop) = getFoldRelStartStop($fold);
   foreach my $product (@{$products}) {
       my($side,$type,$productList,$productCount,$maxProductCount,$adjProductcount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
       foreach my $annotProduct (@{$annotProducts->{$hairpinId}}) {
	   my($chrom,$start,$stop,$strand,$id,$name) = @{$annotProduct};
	   my($annotRelStart,$annotRelStop) = getRelStartStop($location,$start,$stop,$strand);
	   my $overlap = getOverlap($relStart,$relStop,$annotRelStart,$annotRelStop);
	   if ($overlap) {
	       if ($strand eq $productStrand) {
		   my $offset = $annotRelStart - $relStart;
		   if ($annotProductOffsets{$name}) {
		       if (abs($offset) < abs($annotProductOffsets{$name})) {
			   $annotProductOffsets{$name} = $offset;
			   $annotProductCandidates{$name} = $product;
		       }
		   } else {
		       $annotProductOffsets{$name} = $offset;
		       $annotProductCandidates{$name} = $product;		    
		   }
	       }
	   }
       }
   }
   #filling @productAnnotations array and filtering out all products that aren't within $minOffsetLimit
   foreach my $annotProduct (@{$annotProducts->{$hairpinId}}) {
       my($chrom,$start,$stop,$strand,$id,$name) = @{$annotProduct};
       if ($annotProductCandidates{$name}) {
	   if (abs($annotProductOffsets{$name}) <= $minOffsetLimit) {
	       my $product = $annotProductCandidates{$name};
	       my $annotOffset = $annotProductOffsets{$name};
	       my($side,$type,$productList,$productCount,$maxProductCount,$adjProductcount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand) = @{$product};
	       push(@productAnnotations, [$side,$type,$relStart,$relStop,$strand,$gStart,$id,$name,$annotOffset]);
	   }
       }
   }
   my @sortedProductAnnotations = sort {$a->[2] <=> $b->[2]} @productAnnotations; #sorted based on relstart
   return \@sortedProductAnnotations;
}

sub annotateHairpin {
    my($location,$fold,$annotHairpins,$annotProducts) = @_;
    my $foldLocation = getFoldGenomicLocation($location,$fold);
    my $hpAnnotationId = "-";
    my $hpAnnotationName = "-";
    my $maxOverlap = 0;
    my ($hairpinCoverage,$annotCoverage,$productCoverage,$numAnnotProducts,$annotProdsOverlapped) = (0,0,0,0,0);
    my($foldChrom,$foldStart,$foldStop,$foldStrand) = parseLocation($foldLocation);
    #retrieving greatest overlapping annotation information
    foreach my $hpInfo (@{$annotHairpins->{$foldChrom}}) {
	my($start,$stop,$strand,$id,$name) = @{$hpInfo};
	if ($strand eq $foldStrand) {
	    my $overlap = getOverlap($foldStart,$foldStop,$start,$stop);
	    if ($overlap > $maxOverlap) {
		$maxOverlap = $overlap;
		$hairpinCoverage = $overlap / ($foldStop - $foldStart + 1);
		$annotCoverage = $overlap / ($stop - $start + 1);
		$hpAnnotationId = $id;
		$hpAnnotationName = $name;
#		print "$id foldstart=$foldStart\tfoldStop=$foldStop\tstart=$start\tstop=$stop\noverlap=$overlap\nhpCov=$hairpinCoverage\nannotCov=$annotCoverage\n";
	    }
	}
    }
    my $mirbaseId = $hpAnnotationId;
    if ($hpAnnotationId eq "-") {
	#no overlapping Sense annotation found
	#retrieving greatest overlappign antisense annotation information.
	foreach my $hpInfo (@{$annotHairpins->{$foldChrom}}) {
	    my($start,$stop,$strand,$id,$name) = @{$hpInfo};
	    if ($strand ne $foldStrand) {
		my $overlap = getOverlap($foldStart,$foldStop,$start,$stop);
		if ($overlap > $maxOverlap) {
		    $maxOverlap = $overlap;
		    $hairpinCoverage = $overlap / ($foldStop - $foldStart + 1);
		    $annotCoverage = $overlap / ($stop - $start + 1);
		    $mirbaseId = $id;
		    $hpAnnotationId = "as-" . $id;
		    $hpAnnotationName = "as-" . $name;
		}
	    }
	}
    }
    #determining the percentage of annotated products that are overlapped
    unless ($mirbaseId eq "-") {
	my $overlapTotal = 0;
	my $prodLengthTotal = 0;
#	print "$hpAnnotationId\n";
	foreach my $annotProduct (@{$annotProducts->{$mirbaseId}}) {
	    my($chrom,$start,$stop,$strand,$id,$name) = @{$annotProduct};
#	    print "$name\t$chrom\t$start\t$stop\n";
	    $numAnnotProducts += 1;
	    my $overlap = getOverlap($foldStart,$foldStop,$start,$stop);
	    if ($overlap) {
		$annotProdsOverlapped += 1;
	    }
	    $overlapTotal += $overlap;
	    $prodLengthTotal += ($stop - $start + 1);
	}

#	##TESTING 
#	if ($prodLengthTotal == 0) {
#	    print "miRBaseId = $mirbaseId\n";
#	}

#	print "mirBaseId = $mirbaseId\n";
#	print "prodLengthTotal = $prodLengthTotal\n";
	$productCoverage = ($prodLengthTotal) ? $overlapTotal / $prodLengthTotal : 1;
    }
#    unless ($mirbaseId eq "-") {
#	print "hairpinCoverage = $hairpinCoverage\nannotCoverage = $annotCoverage\n";
#    }
    return($hpAnnotationId,$hpAnnotationName,$hairpinCoverage,$annotCoverage,$productCoverage,$numAnnotProducts,$annotProdsOverlapped);
}

sub annotateHairpin_old {
    my($location,$fold,$annotHairpins,$annotProducts) = @_;
    my $foldLocation = getFoldGenomicLocation($location,$fold);
    my $hpAnnotationId = "-";
    my $hpAnnotationName = "-";
    my $maxOverlap = 0;
    my ($hairpinCoverage,$annotCoverage,$productCoverage,$numAnnotProducts) = (0,0,0,0);
    my($foldChrom,$foldStart,$foldStop,$foldStrand) = parseLocation($foldLocation);
    #retrieving greatest overlapping annotation information
    foreach my $hpInfo (@{$annotHairpins->{$foldChrom}}) {
	my($start,$stop,$strand,$id,$name) = @{$hpInfo};
	if ($strand eq $foldStrand) {
	    my $overlap = getOverlap($foldStart,$foldStop,$start,$stop);
	    if ($overlap > $maxOverlap) {
		$maxOverlap = $overlap;
		$hairpinCoverage = $overlap / ($foldStop - $foldStart + 1);
		$annotCoverage = $overlap / ($stop - $start + 1);
		$hpAnnotationId = $id;
		$hpAnnotationName = $name;
#		print "$id foldstart=$foldStart\tfoldStop=$foldStop\tstart=$start\tstop=$stop\noverlap=$overlap\nhpCov=$hairpinCoverage\nannotCov=$annotCoverage\n";
	    }
	}
    }
    my $mirbaseId = $hpAnnotationId;
    if ($hpAnnotationId eq "-") {
	#no overlapping Sense annotation found
	#retrieving greatest overlappign antisense annotation information.
	foreach my $hpInfo (@{$annotHairpins->{$foldChrom}}) {
	    my($start,$stop,$strand,$id,$name) = @{$hpInfo};
	    if ($strand ne $foldStrand) {
		my $overlap = getOverlap($foldStart,$foldStop,$start,$stop);
		if ($overlap > $maxOverlap) {
		    $maxOverlap = $overlap;
		    $hairpinCoverage = $overlap / ($foldStop - $foldStart + 1);
		    $annotCoverage = $overlap / ($stop - $start + 1);
		    $mirbaseId = $id;
		    $hpAnnotationId = "as-" . $id;
		    $hpAnnotationName = "as-" . $name;
		}
	    }
	}
    }
    #determining the percentage of annotated products that are overlapped
    unless ($mirbaseId eq "-") {
	my $overlapTotal = 0;
	my $prodLengthTotal = 0;
#	print "$hpAnnotationId\n";
	foreach my $annotProduct (@{$annotProducts->{$mirbaseId}}) {
	    my($chrom,$start,$stop,$strand,$id,$name) = @{$annotProduct};
#	    print "$name\t$chrom\t$start\t$stop\n";
	    $numAnnotProducts += 1;
	    $overlapTotal += getOverlap($foldStart,$foldStop,$start,$stop);
	    $prodLengthTotal += ($stop - $start + 1);
	}

#	##TESTING 
#	if ($prodLengthTotal == 0) {
#	    print "miRBaseId = $mirbaseId\n";
#	}

#	print "mirBaseId = $mirbaseId\n";
#	print "prodLengthTotal = $prodLengthTotal\n";
	$productCoverage = $overlapTotal / $prodLengthTotal;
    }
#    unless ($mirbaseId eq "-") {
#	print "hairpinCoverage = $hairpinCoverage\nannotCoverage = $annotCoverage\n";
#    }
    return($hpAnnotationId,$hpAnnotationName,$hairpinCoverage,$annotCoverage,$productCoverage,$numAnnotProducts);
}

sub getOtherAnnotations {
    my($otherAnnotationsBuckets,$location,$bucketSize) = @_;
    my @otherAnnotations;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $nearbyAnnotations = getBucketData($otherAnnotationsBuckets,$location,$bucketSize);
    foreach my $annotInfo (@{$nearbyAnnotations}) {
	my($annoStart,$annoStop,$annoId,$score,$annoStrand) = @{$annotInfo};
	if (($strand eq $annoStrand) && getOverlap($start,$stop,$annoStart,$annoStop)) {
	    push(@otherAnnotations,$annotInfo);
	}
    }
    return \@otherAnnotations;
}

sub loadOtherAnnotations {
    my($otherAnnotFileList,$bucketSize) = @_;
    die "bucketSize entry not loaded (bucketSize = 0)" unless ($bucketSize);
    print "loading $otherAnnotFileList\n";
    print "bucketSize = $bucketSize\n";
    my %otherAnnotationsBuckets;    #to contain buckets of annotations
    open(OAFL,$otherAnnotFileList) or die "failed to open $otherAnnotFileList\n";
    while (<OAFL>) {
	chomp;
	unless ( /^#/ ) {
	    my($annotType,$annotFile) = split(/\t/);
	    #the following function reads $annotFile into %otherAnnotations
	    print "loading $annotFile into buckets\n";
	    readBedFileIntoBuckets(\%otherAnnotationsBuckets,$annotFile,$bucketSize);
	}
    }
    close(OAFL);
    return \%otherAnnotationsBuckets;
}

sub generateAnnotations {
    my($hairpinsFile,$productFile,$parameters) = @_;
    my $hairpinAnnotationFile = $parameters->{hpAnnotFile} or die "Error: hpAnnotFile parameter not set\n";
    my $productAnnotationFile = $parameters->{prodAnnotFile} or die "Error: prodAnnotFile parameter not set\n";
    my $otherAnnotationFile = $parameters->{otherAnnotFile} or die "Error: otherAnnotFile parameter not set\n";
    my $otherAnnotFileList = $parameters->{otherAnnotFileList};
    my $bucketSize = $parameters->{bucketSize} or die "Error: bucketSize parameter not set";
    my($annotHairpins,$annotProducts) = (0,0);
    if ($parameters->{mirbaseGff}) {
	($annotHairpins,$annotProducts) = readMirbaseGff3($parameters->{mirbaseGff});
    }
    my $otherAnnotationsBuckets = 0;
    if (-e $otherAnnotFileList) {
	print "loading annotations from $otherAnnotFileList\n";
	#$otherAnnotationsBuckets is a hash of arrays keyed by chrom and bucket number
	$otherAnnotationsBuckets = loadOtherAnnotations($otherAnnotFileList,$bucketSize) or die "failed to load $otherAnnotFileList\n";
    }
    open(HPA,">$hairpinAnnotationFile") or die "failed to open $hairpinAnnotationFile for writing\n";
    open(PA,">$productAnnotationFile") or die "failed to open $productAnnotationFile for writing\n";
    open(OA,">$otherAnnotationFile") or die "failed to open $otherAnnotationFile for writing\n";
    print HPA "#tag\tlocation\tAnnotId\thpAnnotName\thpCoverage\tannotCoverage\tannotMatCoverage\tannotMatCount\n";
    print PA "#tag\tside\ttype\trelStart\trelStop\tstrand\tannotMatId\tannotMatName\tannotOffset\n";
    print OA "#tag\totherAnnotations\n";
    open(HL,$hairpinsFile) or die "failed to open $hairpinsFile\n";
    open(HPL,$productFile) or die "failed to open $productFile\n";
    my $nextHPProdTag = 0;
    my $nextHPProd;
    while (<HL>) {
	chomp;
	unless ( /^#/ ) {
	    my($tag,$chrom,$start,$stop,$strand,$lc,$rc,$totalSense,$totalAntisense,$mfe,$sequence,$fold) = split(/\t/);
	    my $location = "$chrom:$start..$stop:$strand";
	    my $senseProducts;
	    my $antiSenseProducts;
	    #loading remaining product info from previous loop
	    if ($tag eq $nextHPProdTag) {
		my($prodTag,$sidetype,$type,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$prodStrand) = @{$nextHPProd};
#$side,$type,$productList,$productCount,$maxProductCount,$adjProductcount,$adjMaxProductCount,$relStart,$relStop,$offset,$gStart,$productStrand
#		print $prodTag . "\n";
		if ($prodStrand eq $strand) {
		    push(@{$senseProducts},[$sidetype,$type,0,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,0,0,$prodStrand]);
		} else {
		    push(@{$antiSenseProducts},[$sidetype,$type,0,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,0,0,$prodStrand]);
		}
	    } else {
		print "warning: $tag not equal to $nextHPProdTag in generateAnnotations().  some products may not get annotated.\n" if ($nextHPProdTag);
	    }
	    #loading products from file
	    while (<HPL>) {
		unless ( /^#/ ) {
		    my($prodTag,$sidetype,$type,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$prodStrand) = split(/\t/,$_);
#		    print $prodTag . "\n";
		    if (($prodTag eq $tag) && ($prodStrand eq $strand)) {
			push(@{$senseProducts},[$sidetype,$type,0,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,0,0,$prodStrand]);
		    } elsif ($prodTag eq $tag) {
			push(@{$antiSenseProducts},[$sidetype,$type,0,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,0,0,$prodStrand]);
		    } else {
			#saving product for next loop of the hairpin file
			$nextHPProdTag = $prodTag;
			$nextHPProd = [$prodTag,$sidetype,$type,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop,$prodStrand];
			last;
		    }
		}
	    }
	    my $hpAnnotationId = "-";
	    my $hpAnnotationName = "-";
	    my ($hairpinCoverage,$annotCoverage,$productCoverage,$numAnnotProducts) = (0,0,0,0);
	    my ($senseProductAnnotations,$antiSenseProductAnnotations);
	    my ($otherAnnotations);
	    if ($parameters->{mirbaseGff}) {
		($hpAnnotationId,$hpAnnotationName,$hairpinCoverage,$annotCoverage,$productCoverage,$numAnnotProducts) = annotateHairpin($location,$fold,$annotHairpins,$annotProducts);
		unless ($hpAnnotationId eq "-") {
		    if ($hpAnnotationId =~ /^as-/) {
			my($asHPAnnotId) = $hpAnnotationId =~ /^as-(.*)$/;
			$antiSenseProductAnnotations = annotateProducts($location,$fold,$antiSenseProducts,$annotProducts,$asHPAnnotId);
		    } else {
			$senseProductAnnotations = annotateProducts($location,$fold,$senseProducts,$annotProducts,$hpAnnotationId);
		    }
		}	    
	    }
	    unless ($hpAnnotationId eq "-") {
		print HPA "$tag\t$location\t$hpAnnotationId\t$hpAnnotationName\t$hairpinCoverage\t$annotCoverage\t$productCoverage\t$numAnnotProducts\n";
	    }
	    if ($senseProductAnnotations) {
		foreach my $prodAnnotInfo (@{$senseProductAnnotations}) {
		    my($side,$type,$relStart,$relStop,$strand,$gStart,$id,$name,$annotOffset) = @{$prodAnnotInfo};
		    print PA "$tag\t$side\t$type\t$relStart\t$relStop\t$strand\t$id\t$name\t$annotOffset\n";
		}
	    }
	    if ($antiSenseProductAnnotations) {
		foreach my $prodAnnotInfo (@{$antiSenseProductAnnotations}) {
		    my($side,$type,$relStart,$relStop,$strand,$gStart,$id,$name,$annotOffset) = @{$prodAnnotInfo};
		    print PA "$tag\t$side\t$type\t$relStart\t$relStop\t$strand\t$id\t$name\t$annotOffset\n";
		}
	    }
	    if ((-e $otherAnnotFileList) && ($hpAnnotationId eq "-")) {
		my $otherAnnotations = getOtherAnnotations($otherAnnotationsBuckets,$location,$bucketSize);
		if (@{$otherAnnotations}) {
		    my $otherAnnotLine = "$tag\t";
		    for (my $oaEntry = 0; $oaEntry < (@{$otherAnnotations}); $oaEntry++) { 
			my $otherAnnotInfo = $otherAnnotations->[$oaEntry];
			my($annoStart,$annoStop,$annoId,$score,$annoStrand) = @{$otherAnnotInfo};
			$otherAnnotLine .= "$annoId";
			$otherAnnotLine .= "," unless ($oaEntry == @{$otherAnnotations} - 1);
		    }
		    print OA "$otherAnnotLine\n";
		}
	    }   
	}
    }
    close(HL);
    close(HPL);
    close(HPA);
    close(PA);
    close(OA);
}

sub loadOtherAnnotationsFile {
    my($otherAnnotFile) = @_;
    my %otherAnnotData;
    open(OAF,$otherAnnotFile) or die "failed to open $otherAnnotFile\n";
    while (<OAF>) {
	chomp;
	unless ( /^#/ ) {
	    my($tag,$annotation) = split(/\t/);
	    $otherAnnotData{$tag} = $annotation;
	}
    }
    return \%otherAnnotData;
}

sub loadMatAnnotationsFile {
    my($matAnnotFile) = @_;
    my %matAnnotData;
    open(MAF,$matAnnotFile) or die "failed to open $matAnnotFile\n";
    while (<MAF>) {
	chomp;
	unless ( /^#/ ) {
	    my($tag,$side,$type,$relStart,$relStop,$strand,$annotMatId,$annotMatName,$annotOffset) = split(/\t/);
	    push(@{$matAnnotData{$tag}},[$side,$type,$relStart,$relStop,$strand,$annotMatId,$annotMatName,$annotOffset]);
	}
    }
    return \%matAnnotData;
}

sub loadPrecursorAnnotationsFile {
    my($precAnnotFile) = @_;
    my %precAnnotData;
    open(PAF,$precAnnotFile) or die "failed to open $precAnnotFile\n";
    while (<PAF>) {
	chomp;
	unless ( /^#/ ) {
	    my($tag,$location,$annotId,$annotName,$hpCoverage,$annotCoverage,$annotMatCoverage,$annotMatCount) = split(/\t/);
	    die "more than one annotation for $tag\n" if ($precAnnotData{$tag});
	    $precAnnotData{$tag} = [$location,$annotId,$annotName,$hpCoverage,$annotCoverage,$annotMatCoverage,$annotMatCount];
	}
    }
    return \%precAnnotData;
}

sub annotatePredictions {
    my($predictions,$parameters) = @_;
    my $minHPCoverage = $parameters->{aliasMinHPCoverage};
    my $minPrecCoverage = $parameters->{aliasMinPrecCoverage};
    my $minMatCoverage = $parameters->{aliasMinMatCoverage};
    my $minHPCoverageSingle = $parameters->{aliasMinHPCoverageSingle};
    my $minPrecCoverageSingle = $parameters->{aliasMinPrecCoverageSingle};
    my $minMatCoverageSingle = $parameters->{aliasMinMatCoverageSingle};
    #print "minHPCoverage = $minHPCoverage\tminPrecCoverage = $minPrecCoverage\tminMatCoverage = $minMatCoverage\n";
    #print "minHPCoverageSinge = $minHPCoverageSingle\tminPrecCoverageSingle = $minPrecCoverageSingle\tminMatCoverageSingle = $minMatCoverageSingle\n";
    my $precursorAnnotations = (-e $parameters->{hpAnnotFile}) ? loadPrecursorAnnotationsFile($parameters->{hpAnnotFile}) : {};
    my $otherAnnotations = (-e $parameters->{otherAnnotFile}) ? loadOtherAnnotationsFile($parameters->{otherAnnotFile}) : {};
    my %annotations;
    foreach my $chrom (keys %{$predictions}) {
	foreach my $hpInfo (@{$predictions->{$chrom}}) {
	    my($start,$stop,$strand,$geneId,$name) = @{$hpInfo};
	    if($precursorAnnotations->{$name}) {
		my($location,$annotId,$annotName,$hpCoverage,$annotCoverage,$annotMatCoverage,$annotMatCount) = @{$precursorAnnotations->{$name}};
		###############
		if (($hpCoverage == 1) || ($annotCoverage == 1) || ($annotMatCoverage >= $minMatCoverage)) {
			$annotations{$name} = $annotName;
		} else {
			$annotations{$name} = $annotName . "_overlap";
		}
#		if ($annotMatCount == 1) {
#		    if (($hpCoverage > $minHPCoverageSingle) && ($annotCoverage > $minPrecCoverageSingle) && ($annotMatCoverage > $minMatCoverageSingle)) {
#			$annotations{$name} = $annotName;		    
#		    } else {
#			$annotations{$name} = $annotName . "_overlap";
#		    }
#		} else {
#		    if (($hpCoverage > $minHPCoverage) && ($annotCoverage > $minPrecCoverage) && ($annotMatCoverage > $minMatCoverage)) {
#			$annotations{$name} = $annotName;		    
#		    } else {
#			$annotations{$name} = $annotName . "_overlap";
#		    }
#		}
	    } elsif ($otherAnnotations->{$name}) {
		$annotations{$name} = $otherAnnotations->{$name};
	    }
	}
    }
    return \%annotations;
}

########################################
#  General Machine Learning Functions  #
########################################

sub createRemovedScoresSet {
    #returns an array of scores not to use in training or applying the random forest model
    my($MLParameters) = @_;
    my @scoresToRemove;
    foreach my $score (keys %{$MLParameters}) {
	unless ($MLParameters->{$score}) {
	    push(@scoresToRemove, $score);
	}
    }
    #the model was origially trained with this variable named nonMirNeighborCount.  The variable was later
    #changed to be called neighborCount because miRs aren't used in finding it.  So for now until a new model 
    #is created this line will be used to change the variable name.
    unless ($MLParameters->{"neighborCount"}) {
	push(@scoresToRemove, "nonMirNeighborCount");
    }
    return \@scoresToRemove;
}

sub countNegClasses {
    my($trainFile) = @_;
    my $totalNeg = 0;
    open(TF,$trainFile) or die "failed to open $trainFile\n";
    my $header = <TF>;
    while (<TF>) {
	chomp;
	my($geneId,$class) = split(/\t/);
	if ($class == -1) {
	    $totalNeg++;
	}
    }
    print "$totalNeg positive classes found in $trainFile\n";
    return $totalNeg;
}

sub countPosClasses {
    my($trainFile) = @_;
    my $totalPos = 0;
    open(TF,$trainFile) or die "failed to open $trainFile\n";
    my $header = <TF>;
    while (<TF>) {
	chomp;
	my($geneId,$class) = split(/\t/);
	if ($class == 1) {
	    $totalPos++;
	}
    }
    print "$totalPos positive classes found in $trainFile\n";
    return $totalPos;
}

sub createRFModel_roseOver {
    my($R,$modelFile,$trainFile,$scoresToRemove,$parameters) = @_;
    my $totalNeg = countNegClasses($trainFile);
    my $numSamples = $totalNeg * 2;
    my $out = $R->set('featureVecTrainFile', $trainFile);
    $out .= $R->set('modelFile', $modelFile);
    $out .= $R->run(q`print(paste("reading in ",featureVecTrainFile, " for the training data", sep=""))`);
    # Loading Libraries and Parameters 
    $out .= $R->run(q`library(ROCR)`);
    $out .= $R->run(q`library(randomForest)`);
    $out .= $R->run(q`library(ROSE)`);
    $out .= $R->run(q`classificationtype <-'C-classification'`);
    $out .= $R->run(q`crossvalidation <- 10`);
    $out .= $R->run(q`options(scipen=999)`);
    #loading train data
    $out .= $R->run(q`inputData <- subset(read.table(featureVecTrainFile, header=TRUE, sep="\t"), select=c(-GeneId))`);
    #the model was origially trained with this variable named nonMirNeighborCount.  The variable was later
    #changed to be called neighborCount because miRs aren't used in finding it.  So for now until a new model 
    #is created this line will be used to change the variable name.
    $out .= $R->run(q`colnames(inputData)[colnames(inputData)=="neighborCount"] <- "nonMirNeighborCount"`);
    foreach my $score (@{$scoresToRemove}) {
	my $command = "if (\"". $score . "\" %in% colnames(inputData)) inputData\$" . $score. " <- NULL";
	$out .= $R->run($command);
	$R->run(q`columnNames <- colnames(inputData)`);
#       #printing column names for testing
#	my $names = $R->get( 'columnNames' );
#	foreach my $name (@{$names}) {
#	    print "$name\t";
#	}
#	print "\n";
    }
    $out .= $R->run(q`inputData$y <- as.factor(inputData$y)`);
    $out .= $R->set('numSamples',$numSamples);
    $out .= $R->run(q`trainData <- ovun.sample(y ~ ., data = inputData, method = "over", N = numSamples)$data`);
    $out .= $R->run(q`trainClass <- as.matrix(subset(trainData, select = y))`);
    $out .= $R->run(q`trainClass <- as.factor(trainClass)`);
    $out .= $R->run(q`trainValues <- as.matrix(subset(trainData, select = -y))`);
    #creating random forest from train data
    $out .= $R->run(q`model <-rfcv(trainy=trainClass,trainx=trainValues,scale = "log" ,cv.fold=crossvalidation)`);
    $out .= $R->run(q`iMIN = which.min(model$error.cv)`);
    my $iMin = $R->get( 'iMIN' );
    $out .= $R->run(q`predicted <-model$predicted[[iMIN]]`);
    $out .= $R->run(q`predict.rf <-randomForest(as.factor(trainClass) ~ .,data=trainValues, mtry=iMIN, ntree=1000,keep.forest=TRUE,importance=TRUE,test=data$val)`);
    print "model created\nmTry = $iMin\n";
#    $out .= $R->run(q`predict.rf$predicted = predicted[,"y"]`);
    $out .= $R->run(q`predict.rf$predicted = predicted`);
    print "\n";
    $out .= $R->run(q`print(paste("saving model to ", modelFile, "", sep=""))`);
    $out .= $R->run(q`save(predict.rf,file=modelFile)`);
    print "$out\n";
}

sub createRFModel_roseUnder {
    my($R,$modelFile,$trainFile,$scoresToRemove,$parameters) = @_;
    my $totalPos = countPosClasses($trainFile);
    my $numSamples = $totalPos * 2;
    my $out = $R->set('featureVecTrainFile', $trainFile);
    $out .= $R->set('modelFile', $modelFile);
    $out .= $R->run(q`print(paste("reading in ",featureVecTrainFile, " for the training data", sep=""))`);
    # Loading Libraries and Parameters 
    $out .= $R->run(q`library(ROCR)`);
    $out .= $R->run(q`library(randomForest)`);
    $out .= $R->run(q`library(ROSE)`);
    $out .= $R->run(q`classificationtype <-'C-classification'`);
    $out .= $R->run(q`crossvalidation <- 10`);
    $out .= $R->run(q`options(scipen=999)`);
    #loading train data
    $out .= $R->run(q`inputData <- subset(read.table(featureVecTrainFile, header=TRUE, sep="\t"), select=c(-GeneId))`);
    #the model was origially trained with this variable named nonMirNeighborCount.  The variable was later
    #changed to be called neighborCount because miRs aren't used in finding it.  So for now until a new model 
    #is created this line will be used to change the variable name.
    $out .= $R->run(q`colnames(inputData)[colnames(inputData)=="neighborCount"] <- "nonMirNeighborCount"`);
    foreach my $score (@{$scoresToRemove}) {
	my $command = "if (\"". $score . "\" %in% colnames(inputData)) inputData\$" . $score. " <- NULL";
	$out .= $R->run($command);
	$R->run(q`columnNames <- colnames(inputData)`);
#       #printing column names for testing
#	my $names = $R->get( 'columnNames' );
#	foreach my $name (@{$names}) {
#	    print "$name\t";
#	}
#	print "\n";
    }
    $out .= $R->run(q`inputData$y <- as.factor(inputData$y)`);
    $out .= $R->set('numSamples',$numSamples);
    $out .= $R->run(q`trainData <- ovun.sample(y ~ ., data = inputData, method = "under", N = numSamples, seed = 1)$data`);
    $out .= $R->run(q`trainClass <- as.matrix(subset(trainData, select = y))`);
    $out .= $R->run(q`trainClass <- as.factor(trainClass)`);
    $out .= $R->run(q`trainValues <- as.matrix(subset(trainData, select = -y))`);
    #creating random forest from train data
    $out .= $R->run(q`model <-rfcv(trainy=trainClass,trainx=trainValues,scale = "log" ,cv.fold=crossvalidation)`);
    $out .= $R->run(q`iMIN = which.min(model$error.cv)`);
    my $iMin = $R->get( 'iMIN' );
    $out .= $R->run(q`predicted <-model$predicted[[iMIN]]`);
    $out .= $R->run(q`predict.rf <-randomForest(as.factor(trainClass) ~ .,data=trainValues, mtry=iMIN, ntree=1000,keep.forest=TRUE,importance=TRUE,test=data$val)`);
    print "model created\nmTry = $iMin\n";
#    $out .= $R->run(q`predict.rf$predicted = predicted[,"y"]`);
    $out .= $R->run(q`predict.rf$predicted = predicted`);
    print "\n";
    $out .= $R->run(q`print(paste("saving model to ", modelFile, "", sep=""))`);
    $out .= $R->run(q`save(predict.rf,file=modelFile)`);
    print "$out\n";
}

sub createRFModel_roseSyn {
    my($R,$modelFile,$trainFile,$scoresToRemove,$parameters) = @_;
    my $out = $R->set('featureVecTrainFile', $trainFile);
    $out .= $R->set('modelFile', $modelFile);
    $out .= $R->run(q`print(paste("reading in ",featureVecTrainFile, " for the training data", sep=""))`);
    # Loading Libraries and Parameters 
    $out .= $R->run(q`library(ROCR)`);
    $out .= $R->run(q`library(randomForest)`);
    $out .= $R->run(q`library(ROSE)`);
    $out .= $R->run(q`classificationtype <-'C-classification'`);
    $out .= $R->run(q`crossvalidation <- 10`);
    $out .= $R->run(q`options(scipen=999)`);
    #loading train data
    $out .= $R->run(q`inputData <- subset(read.table(featureVecTrainFile, header=TRUE, sep="\t"), select=c(-GeneId))`);
    #the model was origially trained with this variable named nonMirNeighborCount.  The variable was later
    #changed to be called neighborCount because miRs aren't used in finding it.  So for now until a new model 
    #is created this line will be used to change the variable name.
    $out .= $R->run(q`colnames(inputData)[colnames(inputData)=="neighborCount"] <- "nonMirNeighborCount"`);
    foreach my $score (@{$scoresToRemove}) {
	my $command = "if (\"". $score . "\" %in% colnames(inputData)) inputData\$" . $score. " <- NULL";
	$out .= $R->run($command);
	$R->run(q`columnNames <- colnames(inputData)`);
#       #printing column names for testing
#	my $names = $R->get( 'columnNames' );
#	foreach my $name (@{$names}) {
#	    print "$name\t";
#	}
#	print "\n";
    }
    $out .= $R->run(q`inputData$y <- as.factor(inputData$y)`);
    $out .= $R->run(q`trainData <- ROSE(y ~ ., data = inputData, seed = 1)$data`);
    $out .= $R->run(q`trainClass <- as.matrix(subset(trainData, select = y))`);
    $out .= $R->run(q`trainClass <- as.factor(trainClass)`);
    $out .= $R->run(q`trainValues <- as.matrix(subset(trainData, select = -y))`);
    #creating random forest from train data
    $out .= $R->run(q`model <-rfcv(trainy=trainClass,trainx=trainValues,scale = "log" ,cv.fold=crossvalidation)`);
    $out .= $R->run(q`iMIN = which.min(model$error.cv)`);
    my $iMin = $R->get( 'iMIN' );
    $out .= $R->run(q`predicted <-model$predicted[[iMIN]]`);
    $out .= $R->run(q`predict.rf <-randomForest(as.factor(trainClass) ~ .,data=trainValues, mtry=iMIN, ntree=1000,keep.forest=TRUE,importance=TRUE,test=data$val)`);
    print "model created\nmTry = $iMin\n";
#    $out .= $R->run(q`predict.rf$predicted = predicted[,"y"]`);
    $out .= $R->run(q`predict.rf$predicted = predicted`);
    print "\n";
    $out .= $R->run(q`print(paste("saving model to ", modelFile, "", sep=""))`);
    $out .= $R->run(q`save(predict.rf,file=modelFile)`);
    print "$out\n";
}

sub createRFModel_wgtdStratSampling {
    my($R,$modelFile,$trainFile,$scoresToRemove,$posSampSize,$negSampSize,$posWeight,$negWeight,$parameters) = @_;
    my $out = $R->set('featureVecTrainFile', $trainFile);
    $out .= $R->set('modelFile', $modelFile);
    $out .= $R->set('posSampSize', $posSampSize);
    $out .= $R->set('negSampSize', $negSampSize);
    $out .= $R->set('posWeight', $posWeight);
    $out .= $R->set('negWeight', $negWeight);
    $out .= $R->run(q`print(paste("reading in ",featureVecTrainFile, " for the training data", sep=""))`);
    # Loading Libraries and Parameters 
    $out .= $R->run(q`library(ROCR)`);
    $out .= $R->run(q`library(randomForest)`);
    $out .= $R->run(q`classificationtype <-'C-classification'`);
    $out .= $R->run(q`crossvalidation <- 10`);
    $out .= $R->run(q`options(scipen=999)`);
    #loading train data
    $out .= $R->run(q`trainData <- subset(read.table(featureVecTrainFile, header=TRUE, sep="\t"), select=c(-GeneId))`);
    #the model was origially trained with this variable named nonMirNeighborCount.  The variable was later
    #changed to be called neighborCount because miRs aren't used in finding it.  So for now until a new model 
    #is created this line will be used to change the variable name.
    $out .= $R->run(q`colnames(trainData)[colnames(trainData)=="neighborCount"] <- "nonMirNeighborCount"`);
    foreach my $score (@{$scoresToRemove}) {
	my $command = "if (\"". $score . "\" %in% colnames(trainData)) trainData\$" . $score. " <- NULL";
	$out .= $R->run($command);
	$R->run(q`columnNames <- colnames(trainData)`);
#       #printing column names for testing
#	my $names = $R->get( 'columnNames' );
#	foreach my $name (@{$names}) {
#	    print "$name\t";
#	}
#	print "\n";
    }
    $out .= $R->run(q`trainClass <- as.matrix(subset(trainData, select = y))`);
    $out .= $R->run(q`trainClass <- as.factor(trainClass)`);
    $out .= $R->run(q`trainValues <- as.matrix(subset(trainData, select = -y))`);

#    my $trainFactors = $R->run(q`as.factor(trainClass)`);
#    print "$trainFactors";

#    die "test End\n";

    #creating random forest from train data
    my $iMin;
    if ($parameters->{mTry}) {
	$iMin = $parameters->{mTry};
	$R->set('iMIN',$iMin);
	print "iMin set to $iMin\n";
    } else {
	$out .= $R->run(q`model <-rfcv(trainy=trainClass,trainx=trainValues,scale = "log" ,cv.fold=crossvalidation,sampsize=c('-1'=negSampSize,'1'=posSampSize),classwt=c('-1'=negWeight,'1'=posWeight),strata=as.factor(trainClass))`);
	$out .= $R->run(q`iMIN = which.min(model$error.cv)`);
	$iMin = $R->get( 'iMIN' );
	$out .= $R->run(q`predicted <-model$predicted[[iMIN]]`);
    }
#    my $factor = $R->run(q`as.factor(trainClass)`);
#    print "factor = $factor\n";
    $out .= $R->run(q`predict.rf <-randomForest(as.factor(trainClass) ~ .,data=trainValues, mtry=iMIN, ntree=5000,keep.forest=TRUE,importance=TRUE,test=data$val,sampsize=c('-1'=negSampSize,'1'=posSampSize),classwt=c('-1'=negWeight,'1'=posWeight),strata=as.factor(trainClass))`);
    print "model created\nmTry = $iMin\n";
#    $out .= $R->run(q`predixct.rf$predicted = predicted[,"y"]`);
    unless ($parameters->{mTry}) {
	$out .= $R->run(q`predict.rf$predicted = predicted`);
    }
    print "\n";
    $out .= $R->run(q`print(paste("saving model to ", modelFile, "", sep=""))`);
    $out .= $R->run(q`save(predict.rf,file=modelFile)`);
    print "$out\n";
}

sub createRFModel_sampSize {
    my($R,$modelFile,$trainFile,$scoresToRemove,$posSampSize,$negSampSize,$parameters) = @_;
    my $out = $R->set('featureVecTrainFile', $trainFile);
    $out .= $R->set('modelFile', $modelFile);
    $out .= $R->set('posSampSize', $posSampSize);
    $out .= $R->set('negSampSize', $negSampSize);
    if ($parameters->{nodeSize}) {
	print "in createRFModel_sampSize: using $parameters->{nTree} as the number of trees\n";
	$out .= $R->set('nodeSize',$parameters->{nodeSize});
    } else {
	$R->set('nodeSize', 1);
    }
    if ($parameters->{nTree}) {
	print "in createRFModel_sampSize: using $parameters->{nTree} as the number of trees\n";
	$out .= $R->set('numTrees', $parameters->{nTree});
    } else {
       	$out .= $R->set('numTrees', 1000);
    }
    $out .= $R->run(q`print(paste("reading in ",featureVecTrainFile, " for the training data", sep=""))`);
    # Loading Libraries and Parameters 
    $out .= $R->run(q`library(ROCR)`);
    $out .= $R->run(q`library(randomForest)`);
    $out .= $R->run(q`classificationtype <-'C-classification'`);
    $out .= $R->run(q`crossvalidation <- 10`);
    $out .= $R->run(q`options(scipen=999)`);
    #loading train data
    $out .= $R->run(q`trainData <- subset(read.table(featureVecTrainFile, header=TRUE, sep="\t"), select=c(-GeneId))`);
    #the model was origially trained with this variable named nonMirNeighborCount.  The variable was later
    #changed to be called neighborCount because miRs aren't used in finding it.  So for now until a new model 
    #is created this line will be used to change the variable name.
    $out .= $R->run(q`colnames(trainData)[colnames(trainData)=="neighborCount"] <- "nonMirNeighborCount"`);
    foreach my $score (@{$scoresToRemove}) {
	my $command = "if (\"". $score . "\" %in% colnames(trainData)) trainData\$" . $score. " <- NULL";
	$out .= $R->run($command);
	$R->run(q`columnNames <- colnames(trainData)`);
#       #printing column names for testing
#	my $names = $R->get( 'columnNames' );
#	foreach my $name (@{$names}) {
#	    print "$name\t";
#	}
#	print "\n";
    }
    $out .= $R->run(q`trainClass <- as.matrix(subset(trainData, select = y))`);
    $out .= $R->run(q`trainClass <- as.factor(trainClass)`);
    $out .= $R->run(q`trainValues <- as.matrix(subset(trainData, select = -y))`);

#    my $trainFactors = $R->run(q`as.factor(trainClass)`);
#    print "$trainFactors";

#    die "test End\n";

    #creating random forest from train data
    my $iMin;
    if ($parameters->{mTry}) {
	$iMin = $parameters->{mTry};
	$R->set('iMIN',$iMin);
	print "iMin set to $iMin\n";
    } else {
	$out .= $R->run(q`model <-rfcv(trainy=trainClass,trainx=trainValues,scale = "log" ,cv.fold=crossvalidation,sampsize=c('-1'=negSampSize,'1'=posSampSize),strata=as.factor(trainClass))`);
	$out .= $R->run(q`iMIN = which.min(model$error.cv)`);
	$iMin = $R->get( 'iMIN' );
	$out .= $R->run(q`predicted <-model$predicted[[iMIN]]`);
    }
#    my $factor = $R->run(q`as.factor(trainClass)`);
#    print "factor = $factor\n";
    $out .= $R->run(q`predict.rf <-randomForest(as.factor(trainClass) ~ .,data=trainValues, mtry=iMIN, ntree=numTrees, nodesize=nodeSize, keep.forest=TRUE,importance=TRUE,test=data$val,sampsize=c('-1'=negSampSize,'1'=posSampSize),strata=as.factor(trainClass))`);
    print "model created\nmTry = $iMin\n";
#    $out .= $R->run(q`predixct.rf$predicted = predicted[,"y"]`);
    unless ($parameters->{mTry}) {
	$out .= $R->run(q`predict.rf$predicted = predicted`);
    }
    print "\n";
    $out .= $R->run(q`print(paste("saving model to ", modelFile, "", sep=""))`);
    $out .= $R->run(q`save(predict.rf,file=modelFile)`);
    print "$out\n";
}

sub createRFModel_classwt {
    my($R,$modelFile,$trainFile,$scoresToRemove,$posWeight,$negWeight,$parameters) = @_;
    my $out = $R->set('featureVecTrainFile', $trainFile);
    $out .= $R->set('modelFile', $modelFile);
    $out .= $R->set('posWeight', $posWeight);
    $out .= $R->set('negWeight', $negWeight);
    $out .= $R->run(q`print(paste("reading in ",featureVecTrainFile, " for the training data", sep=""))`);
    # Loading Libraries and Parameters 
    $out .= $R->run(q`library(ROCR)`);
    $out .= $R->run(q`library(randomForest)`);
    $out .= $R->run(q`classificationtype <-'C-classification'`);
    $out .= $R->run(q`crossvalidation <- 10`);
    $out .= $R->run(q`options(scipen=999)`);
    #loading train data
    $out .= $R->run(q`trainData <- subset(read.table(featureVecTrainFile, header=TRUE, sep="\t"), select=c(-GeneId))`);
    #the model was origially trained with this variable named nonMirNeighborCount.  The variable was later
    #changed to be called neighborCount because miRs aren't used in finding it.  So for now until a new model 
    #is created this line will be used to change the variable name.
    $out .= $R->run(q`colnames(trainData)[colnames(trainData)=="neighborCount"] <- "nonMirNeighborCount"`);
    foreach my $score (@{$scoresToRemove}) {
	my $command = "if (\"". $score . "\" %in% colnames(trainData)) trainData\$" . $score. " <- NULL";
	$out .= $R->run($command);
	$R->run(q`columnNames <- colnames(trainData)`);
#       #printing column names for testing
#	my $names = $R->get( 'columnNames' );
#	foreach my $name (@{$names}) {
#	    print "$name\t";
#	}
#	print "\n";
    }
    $out .= $R->run(q`trainClass <- as.matrix(subset(trainData, select = y))`);
    $out .= $R->run(q`trainClass <- as.factor(trainClass)`);
    $out .= $R->run(q`trainValues <- as.matrix(subset(trainData, select = -y))`);
    #creating random forest from train data
    $out .= $R->run(q`model <-rfcv(trainy=trainClass,trainx=trainValues,scale = "log" ,cv.fold=crossvalidation)`);
    $out .= $R->run(q`iMIN = which.min(model$error.cv)`);
    my $iMin = $R->get( 'iMIN' );
    $out .= $R->run(q`predicted <-model$predicted[[iMIN]]`);
#    my $factor = $R->run(q`as.factor(trainClass)`);
#    print "factor = $factor\n";
    $out .= $R->run(q`predict.rf <-randomForest(as.factor(trainClass) ~ .,data=trainValues, mtry=iMIN, ntree=1000,keep.forest=TRUE,importance=TRUE,test=data$val, classwt=c(negWeight,posWeight))`);
    print "model created\nmTry = $iMin\n";
#    $out .= $R->run(q`predict.rf$predicted = predicted[,"y"]`);
    $out .= $R->run(q`predict.rf$predicted = predicted`);
    print "\n";
    $out .= $R->run(q`print(paste("saving model to ", modelFile, "", sep=""))`);
    $out .= $R->run(q`save(predict.rf,file=modelFile)`);
    print "$out\n";
}

sub createRFModel {
    my($R,$modelFile,$trainFile,$scoresToRemove,$parameters) = @_;
    my $out = $R->set('featureVecTrainFile', $trainFile);
    if ($parameters->{nTree}) {
	print "in createRFModel: using $parameters->{nTree} as the number of trees\n";
	$out .= $R->set('numTrees', $parameters->{nTree});
    } else {
       	$out .= $R->set('numTrees', 1000);
    }
    $out .= $R->set('modelFile', $modelFile);
    $out .= $R->run(q`print(paste("reading in ",featureVecTrainFile, " for the training data", sep=""))`);
    # Loading Libraries and Parameters 
    $out .= $R->run(q`library(ROCR)`);
    $out .= $R->run(q`library(randomForest)`);
    $out .= $R->run(q`classificationtype <-'C-classification'`);
    $out .= $R->run(q`crossvalidation <- 10`);
    $out .= $R->run(q`options(scipen=999)`);
    #loading train data
    $out .= $R->run(q`trainData <- subset(read.table(featureVecTrainFile, header=TRUE, sep="\t"), select=c(-GeneId))`);
    #the model was origially trained with this variable named nonMirNeighborCount.  The variable was later
    #changed to be called neighborCount because miRs aren't used in finding it.  So for now until a new model 
    #is created this line will be used to change the variable name.
    $out .= $R->run(q`colnames(trainData)[colnames(trainData)=="neighborCount"] <- "nonMirNeighborCount"`);
    foreach my $score (@{$scoresToRemove}) {
	my $command = "if (\"". $score . "\" %in% colnames(trainData)) trainData\$" . $score. " <- NULL";
	$out .= $R->run($command);
	$R->run(q`columnNames <- colnames(trainData)`);
#       #printing column names for testing
#	my $names = $R->get( 'columnNames' );
#	foreach my $name (@{$names}) {
#	    print "$name\t";
#	}
#	print "\n";
    }
    $out .= $R->run(q`trainClass <- as.matrix(subset(trainData, select = y))`);
    $out .= $R->run(q`trainClass <- as.factor(trainClass)`);
    $out .= $R->run(q`trainValues <- as.matrix(subset(trainData, select = -y))`);
    #creating random forest from train data
    my $iMin;
    if ($parameters->{mTry}) {
	$iMin = $parameters->{mTry};
	$R->set('iMIN',$iMin);
	print "iMin set to $iMin\n";
    } else {
	$out .= $R->run(q`model <-rfcv(trainy=trainClass,trainx=trainValues,scale = "log" ,cv.fold=crossvalidation)`);
	$out .= $R->run(q`iMIN = which.min(model$error.cv)`);
	$iMin = $R->get( 'iMIN' );
	$out .= $R->run(q`predicted <-model$predicted[[iMIN]]`);
    }
    $out .= $R->run(q`predict.rf <-randomForest(as.factor(trainClass) ~ .,data=trainValues, mtry=iMIN, ntree=numTrees,keep.forest=TRUE,importance=TRUE,test=data$val)`);
    print "model created\nmTry = $iMin\n";
#    $out .= $R->run(q`predict.rf$predicted = predicted[,"y"]`);
    unless ($parameters->{mTry}) {
	$out .= $R->run(q`predict.rf$predicted = predicted`);
    }
    print "\n";
    $out .= $R->run(q`print(paste("saving model to ", modelFile, "", sep=""))`);
    $out .= $R->run(q`save(predict.rf,file=modelFile)`);
    print "$out\n";
}

sub RFModelPredict {
    my($R,$modelFile,$featureVectorFile,$outputFile,$scoresToRemove,$parameters) = @_;
    my $out .= $R->set('modelFile', $modelFile);
    $out .= $R->set('featureVectorFile', $featureVectorFile);
    $out .= $R->set('outputFile', $outputFile);
    $out .= $R->run(q`unlink(outputFile)`);
    $out .= $R->run(q`print(paste("reading in ",modelFile, " for the model", sep=""))`);
    $out .= $R->run(q`print(paste("reading in ",featureVectorFile, " for the data", sep=""))`);
    $out .= $R->run(q`print(paste("reading in ",outputFile, " for the predicted scores file", sep=""))`);
    # Loading Libraries and Parameters 
    $out .= $R->run(q`library(ROCR)`);
    $out .= $R->run(q`library(randomForest)`);
    $out .= $R->run(q`classificationtype <-'C-classification'`);
    $out .= $R->run(q`crossvalidation <- 10`);
    $out .= $R->run(q`options(scipen=999)`);
    #loading test data
    $out .= $R->run(q`testRows <- read.table(featureVectorFile,header=TRUE)`);
    $out .= $R->run(q`num <- nrow(testRows)`);
    $out .= $R->run(q`testData <- subset(testRows, select = c(-GeneId))`);
    #the model was origially trained with this variable named nonMirNeighborCount.  The variable was later
    #changed to be called neighborCount because miRs aren't used in finding it.  So for now until a new model 
    #is created this line will be used to change the variable name.
    $out .= $R->run(q`colnames(testData)[colnames(testData)=="neighborCount"] <- "nonMirNeighborCount"`);
    foreach my $score (@{$scoresToRemove}) {
	my $command = "if (\"". $score . "\" %in% colnames(testData)) testData\$" . $score. " <- NULL";
	$out .= $R->run($command);
	$R->run(q`columnNames <- colnames(testData)`);
#       #printing column names for testing
#	my $names = $R->get( 'columnNames' );
#	foreach my $name (@{$names}) {
#	    print "$name\t";
#	}
#	print "\n";
    }
    $out .= $R->run(q`testClass <- as.matrix(subset(testData, select = y))`);
    $out .= $R->run(q`testValues<- as.matrix(subset(testData, select = c(-y)))`);
    #creating random forest from train data
    $out .= $R->run(q`load(file=modelFile)`);
    $out .= $R->run(q`predict.rf.pr = predict(predict.rf,type="prob",newdata=testValues)[,2]`);
    $out .= $R->run(q`classPredictions <- as.matrix(predict(predict.rf, testValues, decision.values=TRUE))`);
    $out .= $R->run(q`testedGeneIds <- subset(testRows[attributes(predict.rf.pr)$names,], select = GeneId)`);
    $out .= $R->run(q`testedProbs <- as.matrix(predict.rf.pr)`);
    $out .= $R->run(q`testOutput <- data.frame(testedGeneIds,classPredictions,testedProbs)`);
    $out .= $R->run(q`write(t(testOutput),outputFile,ncolumns=3,sep="\t")`);
    print "$out\n";
}

sub featureSelection {
    my($R,$trainFile,$scoresToRemove,$parameters) = @_;
    my $out .= $R->set('featurFile', $trainFile);
    my $label = $parameters->{outputPrefix};
    unless ($parameters->{mTry}) {
	die "parameters mTry is not set\n";
    }
    $out .= $R->set('mTryVal',$parameters->{mTry});
    $label .= "\_featSelect\_mTry\_$parameters->{mTry}";
    $out .= $R->set('label',$label);
    $out .= $R->run(q`library(mlbench)`);
    $out .= $R->run(q`library(randomForest)`);
    $out .= $R->run(q`library("Boruta")`);
    $out .= $R->run(q`pdf(file=paste(label,".pdf", sep=""),width=20,height=10)`);
    $out .= $R->run(q`data<-subset(read.table(featurFile, header=TRUE, sep="\t"), select=c(-GeneId))`);
    $out .= $R->run(q`data[is.na(data)] <- 0`);
    #the model was origially trained with this variable named nonMirNeighborCount.  The variable was later
    #changed to be called neighborCount because miRs aren't used in finding it.  So for now until a new model 
    #is created this line will be used to change the variable name.
    $out .= $R->run(q`colnames(data)[colnames(data)=="neighborCount"] <- "nonMirNeighborCount"`);
    foreach my $score (@{$scoresToRemove}) {
	my $command = "if (\"". $score . "\" %in% colnames(data)) data\$" . $score. " <- NULL";
	$out .= $R->run($command);
	$R->run(q`columnNames <- colnames(data)`);
#       #printing column names for testing
#	my $names = $R->get( 'columnNames' );
#	foreach my $name (@{$names}) {
#	    print "$name\t";
#	}
#	print "\n";
    }
    $out .= $R->run(q`dataSet <-as.matrix(data)`);
    $out .= $R->run(q`dataClass<-as.matrix(subset(data, select = y))`);
    $out .= $R->run(q`dataValues<- as.matrix(subset(data, select = c(-y)))`);
    $out .= $R->run(q`features <- Boruta(y ~ ., data = data, doTrace = 2, ntree =1000,maxRuns=500,mtry=mTryVal)`);
    $out .= $R->run(q`par(mar=c(10, 5.1, 5.1, 3.1),las=2)`);
    $out .= $R->run(q`par(las=2)`);
    $out .= $R->run(q`plot(features,cex.axis=1.2,xlab="",cex.lab=2)`);
    $out .= $R->run(q`attStats(features)`);
    $out .= $R->run(q`print(features)`);
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
	#the model was origially trained with this variable named nonMirNeighborCount.  The variable was later
	#changed to be called neighborCount because miRs aren't used in finding it.  So for now until a new model 
	#is created this line will be used to change the variable name.
	$out .= $R->run(q`colnames(data)[colnames(data)=="neighborCount"] <- "nonMirNeighborCount"`);
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
	my($id,$productReads,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$productInfo};
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

sub getProductGenomicLocation {
    my($product,$chrom) = @_;
    my($id,$reads,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$product};
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
    my($id,$productReads,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$product};
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

sub addOffsetToProductReads {
    #this is a function used to add an offset to the product reads.  buildProducts is often followed by getProductInfo.  getProductInfo
    #is the function that normally adds the offset.  However, buildProducts is used without getProductInfo in the extractProductFeatures step.
    #This function adds a dummy offset.  In the future, recoding buildProducts to add a dummy variable for the offset and placing a comment in all
    #functions that require getProductInfo to be called first may be preferable.
    my($productReads) = @_;
    my @newProductReads;
    foreach my $read (@{$productReads}) {
	my($relStart,$relStop,@readData) = @{$read};
	#a nonnumber was added for the offset to cause an error if a function uses it
	push(@newProductReads,[$relStart,$relStop,".",@readData]);
    }
    return \@newProductReads;
}

sub getProductFeatures {
    my($product,$seq,$class,$rrLocation,$genome,$chromLengths,$prodMLParameters,$parameters) = @_;
    #calculating scores for products
    my($id,$productReads,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$product};
    my $fivePrimeHet = computeFivePrimeHet(addOffsetToProductReads($productReads)) if ($prodMLParameters->{fivePrimeHet});
    my $length = length($seq);
    my $medianLength = getProductMedianLength($productReads) if ($prodMLParameters->{medianLength});
    my $gcContent = getGCcontent($seq,$length) if ($prodMLParameters->{gcContent});
    my $WFC = getWFC($seq) if ($prodMLParameters->{WFC});
    my($aa,$ac,$ag,$at,$ca,$cc,$cg,$ct,$ga,$gc,$gg,$gt,$ta,$tc,$tg,$tt) = getDinucleotideFrequencies($seq) unless ($prodMLParameters->{noDNF});
    my($r10,$r9,$r8,$r7,$r6,$r5,$r4,$r3,$r2,$r1,$s0,
       $f1,$f2,$f3,$f4,$f5,$f6,$f7,$f8,$f9,$f10) = getProductRelStartCounts($product) unless ($prodMLParameters->{noPositionalData});
    my $duplexEnergy = getProductDuplexEnergy($genome,$seq,$rrLocation,$chromLengths,$parameters) if ($prodMLParameters->{duplexEnergy});
    my @productFeatures = ($class);
    #pushing product scores into product vector array
    push(@productFeatures,$adjMaxProductCount) if ($prodMLParameters->{adjMaxProductCount});
    push(@productFeatures,$adjProductCount) if ($prodMLParameters->{adjProductCount});
    push(@productFeatures,$fivePrimeHet) if ($prodMLParameters->{fivePrimeHet});
    push(@productFeatures,$length) if ($prodMLParameters->{length});
    push(@productFeatures,$medianLength) if ($prodMLParameters->{medianLength});
    push(@productFeatures,$gcContent) if ($prodMLParameters->{gcContent});
    push(@productFeatures,$aa) if ($prodMLParameters->{aa});
    push(@productFeatures,$ac) if ($prodMLParameters->{ac});
    push(@productFeatures,$ag) if ($prodMLParameters->{ag});
    push(@productFeatures,$at) if ($prodMLParameters->{at});
    push(@productFeatures,$ca) if ($prodMLParameters->{ca});
    push(@productFeatures,$cc) if ($prodMLParameters->{cc});
    push(@productFeatures,$cg) if ($prodMLParameters->{cg});
    push(@productFeatures,$ct) if ($prodMLParameters->{ct});
    push(@productFeatures,$ga) if ($prodMLParameters->{ga});
    push(@productFeatures,$gc) if ($prodMLParameters->{gc});
    push(@productFeatures,$gg) if ($prodMLParameters->{gg});
    push(@productFeatures,$gt) if ($prodMLParameters->{gt});
    push(@productFeatures,$ta) if ($prodMLParameters->{ta});
    push(@productFeatures,$tc) if ($prodMLParameters->{tc});
    push(@productFeatures,$tg) if ($prodMLParameters->{tg});
    push(@productFeatures,$tt) if ($prodMLParameters->{tt});
    push(@productFeatures,$r10) if ($prodMLParameters->{r10});
    push(@productFeatures,$r9) if ($prodMLParameters->{r9});
    push(@productFeatures,$r8) if ($prodMLParameters->{r8});
    push(@productFeatures,$r7) if ($prodMLParameters->{r7});
    push(@productFeatures,$r6) if ($prodMLParameters->{r6});
    push(@productFeatures,$r5) if ($prodMLParameters->{r5});
    push(@productFeatures,$r4) if ($prodMLParameters->{r4});
    push(@productFeatures,$r3) if ($prodMLParameters->{r3});
    push(@productFeatures,$r2) if ($prodMLParameters->{r2});
    push(@productFeatures,$r1) if ($prodMLParameters->{r1});
    push(@productFeatures,$s0) if ($prodMLParameters->{s0});
    push(@productFeatures,$f1) if ($prodMLParameters->{f1});
    push(@productFeatures,$f2) if ($prodMLParameters->{f2});
    push(@productFeatures,$f3) if ($prodMLParameters->{f3});
    push(@productFeatures,$f4) if ($prodMLParameters->{f4});
    push(@productFeatures,$f5) if ($prodMLParameters->{f5});
    push(@productFeatures,$f6) if ($prodMLParameters->{f6});
    push(@productFeatures,$f7) if ($prodMLParameters->{f7});
    push(@productFeatures,$f8) if ($prodMLParameters->{f8});
    push(@productFeatures,$f9) if ($prodMLParameters->{f9});
    push(@productFeatures,$f10) if ($prodMLParameters->{f10});
    push(@productFeatures,$WFC) if ($prodMLParameters->{WFC});
    push(@productFeatures,$duplexEnergy) if ($prodMLParameters->{duplexEnergy});
    return \@productFeatures;
}

sub printProductTrainData {
    my($bamList,$genomeDir,$chromLengths,$mirBaseGff,$class,$productTrainFile,$prodMLParameters,$parameters) = @_;
    my($annotHairpins,$annotProducts) = readMirbaseGff3($mirBaseGff);
    my $minMirOverlap = 0.60;  #require at least a 60% overlap of mirbase annotation for positive cases
    my $extraScanDistance = 15;  #distance added to the scanRegion so that retrieveReadData won't ignore reads that extend beyond region
    open(PVF,">>".$productTrainFile) or die "failed to open $productTrainFile for appending\n";
    foreach my $hairpinChrom (keys %{$annotHairpins}) {
	my $genome = loadGenome("$genomeDir/$hairpinChrom.fa");	
	foreach my $hairpinInfo (@{$annotHairpins->{$hairpinChrom}}) {
	    my($hairpinStart,$hairpinStop,$strand,$hairpinId,$hairpinName) = @{$hairpinInfo};
	    foreach my $productInfo (@{$annotProducts->{$hairpinId}}) {
		my($chrom,$start,$stop,$strand,$id,$name) = @{$productInfo};
		my $location = "$chrom:$start..$stop:$strand";
		my $scanRegion = "$chrom:". ($start - $extraScanDistance) . ".." . ($stop + $extraScanDistance) . ":$strand";
		my($distinctReads,$uniqueReadCount,$totalDistinctReads,$adjustedReadCount) = retrieveReadData($scanRegion,$bamList,$parameters);
		if ($distinctReads->{$strand}) {  
		    my($products) = buildProducts($location,$distinctReads->{$strand},$strand,$parameters);
		    my ($mirProduct,$percentOverlap) = getGreatestOverlappingProduct($products,$location,$name);
		    my $mirProductLocation = getProductGenomicLocation($mirProduct,$hairpinChrom);
		    my $sequence = getSequence($mirProductLocation,$genome);
		    if ($percentOverlap > $minMirOverlap) {
			my $productName = $name . "_" . $mirProductLocation;
			my $productFeatures = getProductFeatures($mirProduct,$sequence,$class,$location,$genome,$chromLengths,$prodMLParameters,$parameters);
			my $productFeatVectorLine = join("\t", @{$productFeatures});
			print PVF "$productName\t$productFeatVectorLine\n";
		    }
		}
	    }
	}
    }
    close(PVF);
}

sub printProductVectorFile {
    my($readRegionsFile,$bamList,$genomeDir,$chromLengths,$prodMLParameters,$parameters) = @_;
    my($posClass,$negClass) = (1,-1);
    my $productFeatVectorFile = $parameters->{productFeatVectorFile};
    open(PVF,">".$productFeatVectorFile) or die "failed to open $productFeatVectorFile for writing\n";
    open(RRF, $readRegionsFile) or die "failed to open $readRegionsFile\n";
    my $productFeaturesFileHeader = getProdFeaturesFileHeader($prodMLParameters);
    print PVF "$productFeaturesFileHeader\n";
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
		    my($id,$productReads,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$maxRelStart,$maxRelStop,$maxGStart,$productStrand) = @{$product};
		    my $productLocation = getProductGenomicLocation($product,$chrom);
		    my($productId) = @{$product};
		    my $productName = $rrId . "-". $productId . "_" . $productLocation;
		    my $sequence = getSequence($productLocation,$genome);
		    my $classNum = $negClass;
		    my $productFeatures = getProductFeatures($product,$sequence,$classNum,$location,$genome,$chromLengths,$prodMLParameters,$parameters);
		    my $productFeatVectorLine = join("\t", @{$productFeatures});
		    print PVF "$productName\t$productFeatVectorLine\n";
		}
	    }
	}
    }
    close(PVF);
    close(RRF);
}

sub getProdFeaturesFileHeader {
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
    my $productFile = $parameters->{predProductFile} or die "error: predProductFile parameter not set\n";
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

sub readProductFile {
    my($productFile) = @_;
    my %products;
    open(PF, $productFile) or die "faiiled to open $productFile for reading\n";
    while (<PF>) {
	chomp;
	unless ( /^#/ ) {
	    my($name,$side,$type,$total,$totalMostAbundant,$adjTotal,$adjTotalMostAbundant,$start,$stop,$strand,$seq,@libraryCounts) = split("\t",$_);
	    push(@{$products{$name}},[$name,$side,$type,$total,$totalMostAbundant,$adjTotal,$adjTotalMostAbundant,$start,$stop,$strand,$seq,@libraryCounts]);
	}
    }
    close(PF);
    return \%products;
}

sub getMajorProductInfo {
    #this function obtains the interval for the major products within the fold of each hairpin
    my($hairpinsFile,$productFile) = @_;
    my %mpInfo;
    my %miRProducts;
    my $products = readProductFile($productFile);
    open(HP, $hairpinsFile) or die "failed to open $hairpinsFile for reading\n";
    while (<HP>) {
	chomp;
	unless ( /^#/ ) {
	    my($name,$chrom,$start,$stop,$strand,$leftCenter,$rightCenter,$totalSense,$totalAntisense,$mfe,$sequence,$fold) = split("\t", $_);
	    my($middle,$left) = $fold =~ /^((\+*).*[^\+])/;
	    my $foldRelStart = length($left);
	    my $foldRelStop = length($middle);
	    my($foldStart,$foldStop);
	    if ($strand eq '+') {
		$foldStart = $start + $foldRelStart;
		$foldStop = $start + $foldRelStop;
	    } else {
		$foldStart = $stop - $foldRelStop;
		$foldStop = $stop - $foldRelStart;
	    }
	    my $majorProductCount = 0;
	    my $miRCount = 0;
	    foreach my $product (@{$products->{$name}}) {
		my($prodName,$side,$type,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relProdStart,$relProdStop,$productStrand,$productSequence) = @{$product};
		if ($type eq /miR/) {
		    $miRCount += 1;  #counts all miRs including antisense
		}
		push(@{$miRProducts{$name}},$product) if ($type eq "miR");
		my($prodStart,$prodStop) = (-1,-2);  #initialized with dummy variables for testing
		if ($strand eq '+') {
		    $prodStart = $start + $relProdStart;
		    $prodStop = $start + $relProdStop;
		} elsif ($strand eq '-') {
		    $prodStop = $stop - $relProdStart;
		    $prodStart = $stop - $relProdStop;
		}
		if ($prodStart > $prodStop) {
		    die "prodStart is greater than prodStop\n";
		}
		if ($productCount > $majorProductCount && getOverlap($prodStart,$prodStop,$foldStart,$foldStop)) {
		    #in the array below $miRCount is the total number of miRs in the hairpin that contains the major product
		    #this is included since mirCount may be used to avoid filtering hairpins based on size later
		    $mpInfo{$name} = [$prodStart,$prodStop,$productCount,$adjProductCount,$miRCount];
		    $majorProductCount = $productCount;
		}
	    }
	}
    }
    close(HP);
    return(\%mpInfo,\%miRProducts);
}

sub getHairpinARPM {
    my($hairpinsFile) = @_;
    my %hpARPM;
    my $ARPMColName = "senseRPM";
    my $ARPMcol;
    open(HP, $hairpinsFile) or die "Failed to open hairpins file for reading\n";
    while (<HP>) {
	chomp;
	unless ( /^#/ ) {
	    my(@entries) = split("\t",$_);
	    my $name = $entries[0];
	    $hpARPM{$name} = $entries[$ARPMcol];
	} else {
	    my $header = $_;
	    $ARPMcol = getColNum($header,$ARPMColName);
	}
    }
    close(HP);
    return \%hpARPM;
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

sub classifyPredictions {
    my($predictions,$RFPredictions,$hpARPM,$dropOverlaps,$mpInfo,$parameters) = @_;
    my $overlapCutoff = 1;  #require at least this many nucleotides to overlap before dropping
                             #this may be changed later to drop based on overlapping major products
    my $RFCutoff = $parameters->{HPDVCutoff} or die "error: HPDVCutoff not set\n";
    my $minARPM = $parameters->{ARPMCutoff};
    my $letThrough = $parameters->{LetThroughCutoff};
    my $minMajProdCount = $parameters->{minMajProdCount};
    my $minAdjMajProdCount = $parameters->{minAdjMajProdCount};
    my $noAbundCheckMiRCutoff = $parameters->{noAbundCheckMiRCutoff};
    #print "minARPM = $minARPM\n";
    #print "letThrough = $letThrough\n";
    #print "RFCutoff = $RFCutoff\n";
    #print "minMajProdCount = $minMajProdCount\n";
    #print "minAdjMajProdCount = $minAdjMajProdCount\n";
    #print "noAbundCheckMiRCutoff = $noAbundCheckMiRCutoff\n";
    my %positivePredictions;
    my %negativePredictions;
    my %sizeFilteredPredictions;
    my %overlaps;
    foreach my $chrom (keys %{$predictions}) {
	my %hpScores;
	#converting predictions hash to a hash keyed by strand so that we determine overlaps on the same strand
	foreach my $predictionInfo (@{$predictions->{$chrom}}) {
	    my($start,$stop,$strand,$geneId,$name) = @{$predictionInfo};
	    push(@{$hpScores{$strand}}, [$start,$stop,$RFPredictions->{$name},$geneId,$name]);
	}
	#now that the scores are keyed by strand
	#we are determining overlapping hairpins and seperating predictions into positives and negatives
	foreach my $strand (keys %hpScores) {
	    my @usedScores;
	    my @sortedHPScores = sort {$b->[2] <=> $a->[2]} @{$hpScores{$strand}};
	    foreach my $hpScoreInfo (@sortedHPScores) {
		my($start,$stop,$hpScore,$geneId,$name) = @{$hpScoreInfo};
		my $ARPM = $hpARPM->{$name};
		my $predictedOverlap = overlapsPrediction($hpScoreInfo,\@usedScores,$overlapCutoff,$mpInfo);
		my($mpStart,$mpStop,$mpProductCount,$mpAdjProductCount,$mirCount) = @{$mpInfo->{$name}};
		if ($dropOverlaps && $predictedOverlap) {
		    push(@{$overlaps{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
		    #print "dropping $name because it overlaps $predictedOverlap\n";
		} elsif ((($mpProductCount < $minMajProdCount) || ($mpAdjProductCount < $minAdjMajProdCount)) 
			 && ($mirCount < $noAbundCheckMiRCutoff)) {
		    push(@{$sizeFilteredPredictions{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
		    #print "$name filtered by major product size \(mpProductCount = $mpProductCount mpAdjProdCount = $mpAdjProductCount mirCount = $mirCount\)\n";
		} elsif (($ARPM < $minARPM) && ($hpScore < $letThrough)) {
#		    push(@{$sizeFilteredPredictions{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
		    push(@{$negativePredictions{$chrom}}, [$start,$stop,$strand,$geneId,$name]);		
		    #print "$name filtered by hairpin size \(ARPM = " .$ARPM. "score = ".$hpScore."\n";
		} elsif ($hpScore >= $RFCutoff) {
		    push(@{$positivePredictions{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
		} else {
		    push(@{$negativePredictions{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
		}
		push(@usedScores, [$start,$stop,$hpScore,$geneId,$name]);
	    }
	}
    }
    return (\%positivePredictions,\%negativePredictions,\%sizeFilteredPredictions,\%overlaps);
}

sub classifyPredictions_strictARPMCutoff {
    my($predictions,$RFPredictions,$hpARPM,$dropOverlaps,$mpInfo,$parameters) = @_;
    my $overlapCutoff = 1;  #require at least this many nucleotides to overlap before dropping
                             #this may be changed later to drop based on overlapping major products
    my $RFCutoff = $parameters->{HPDVCutoff} or die "error: HPDVCutoff not set\n";
    my $minARPM = $parameters->{ARPMCutoff};
    my $letThrough = $parameters->{LetThroughCutoff};
    my $minMajProdCount = $parameters->{minMajProdCount};
    my $minAdjMajProdCount = $parameters->{minAdjMajProdCount};
    my $noAbundCheckMiRCutoff = $parameters->{noAbundCheckMiRCutoff};
    my $verbose = $parameters->{verbose};
    print "minARPM = $minARPM\n";
    print "letThrough = $letThrough\n";
    print "RFCutoff = $RFCutoff\n";
    print "minMajProdCount = $minMajProdCount\n";
    print "minAdjMajProdCount = $minAdjMajProdCount\n";
    print "noAbundCheckMiRCutoff = $noAbundCheckMiRCutoff\n";
    my %positivePredictions;
    my %negativePredictions;
    my %sizeFilteredPredictions;
    my %overlaps;
    foreach my $chrom (keys %{$predictions}) {
	my %hpScores;
	#converting predictions hash to a hash keyed by strand so that we determine overlaps on the same strand
	foreach my $predictionInfo (@{$predictions->{$chrom}}) {
	    my($start,$stop,$strand,$geneId,$name) = @{$predictionInfo};
	    push(@{$hpScores{$strand}}, [$start,$stop,$RFPredictions->{$name},$geneId,$name]);
	}
	#now that the scores are keyed by strand
	#we are determining overlapping hairpins and seperating predictions into positives and negatives
	foreach my $strand (keys %hpScores) {
	    my @usedScores;
	    my @sortedHPScores = sort {$b->[2] <=> $a->[2]} @{$hpScores{$strand}};
	    foreach my $hpScoreInfo (@sortedHPScores) {
		my($start,$stop,$hpScore,$geneId,$name) = @{$hpScoreInfo};
		my $ARPM = $hpARPM->{$name};
		my $predictedOverlap = overlapsPrediction($hpScoreInfo,\@usedScores,$overlapCutoff,$mpInfo);
		my($mpStart,$mpStop,$mpProductCount,$mpAdjProductCount,$mirCount) = @{$mpInfo->{$name}};
		if ($dropOverlaps && $predictedOverlap) {
		    push(@{$overlaps{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
		    #print "dropping $name because it overlaps $predictedOverlap\n";
		} elsif ((($mpProductCount < $minMajProdCount) || ($mpAdjProductCount < $minAdjMajProdCount)) 
			 && ($mirCount < $noAbundCheckMiRCutoff)) {
		    push(@{$sizeFilteredPredictions{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
		    print "$name filtered by major product size \(mpProductCount = $mpProductCount mpAdjProdCount = $mpAdjProductCount mirCount = $mirCount\)\n" if $verbose;
		} elsif ($ARPM < $minARPM) {
		    push(@{$sizeFilteredPredictions{$chrom}}, [$start,$stop,$strand,$geneId,$name]);		
		    print "$name filtered by hairpin size \(ARPM = " .$ARPM. "score = ".$hpScore."\n" if $verbose;
		} elsif ($hpScore >= $RFCutoff) {
		    push(@{$positivePredictions{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
		} else {
		    push(@{$negativePredictions{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
		}
		push(@usedScores, [$start,$stop,$hpScore,$geneId,$name]);
	    }
	}
    }
    return (\%positivePredictions,\%negativePredictions,\%sizeFilteredPredictions,\%overlaps);
}

sub dropRepeatRegions {
    my($predictions,$repeatRegions,$annotatedMirNames,$hairpinFoldRegions) = @_;
    my %newPredictions;
    foreach my $chrom (%{$predictions}) {
	foreach my $hpInfo (@{$predictions->{$chrom}}) {
	    my($start,$stop,$strand,$geneId,$name) = @{$hpInfo};
	    my($foldStart,$foldStop) = @{$hairpinFoldRegions->{$name}};
	    my $inRepeat = 0;
	    foreach my $repeatRegion (@{$repeatRegions->{$chrom}}) {
		my($repeatStart,$repeatStop) = @{$repeatRegion};
		if (getOverlap($foldStart,$foldStop,$repeatStart,$repeatStop)) {
		    $inRepeat = 1;
		}
	    }
	    if ($inRepeat) {
		if ($annotatedMirNames->{$name}) {
		    print "$name \(annotName=". $annotatedMirNames->{$name} . "\)\($geneId\) found in repeat region\n";
		} else {
		    print "$name \($geneId\) found in repeat region\n";
		}
	    } else {
		push(@{$newPredictions{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
	    }
	}
    }
    return \%newPredictions;
}

sub printPredictedMirGff {
    my($predHairpinClassFile,$productFile,$hairpinsFile,$repeatRegions,$parameters) = @_; 
    my $dropOverlaps = 1;
    my $dropHairpinOverlaps = 1;
    my $minOverlapPercent = 0.7;
#    my $matureAnnotations = (-e $parameters->{prodAnnotFile}) ? loadMatAnnotationsFile($parameters->{prodAnnotFile}) : {};
    print "in predicted mir gff\n";
    print "predHairpinClassFile = $predHairpinClassFile\n";
    print "productFile = $productFile\n";
    print "hairpinsFile = $hairpinsFile\n";
    my($posGff,$negGff,$sizeFilteredGff,$overlapGff) = ("positivePredictions.gff","negativePredictions.gff","sizeFiltered.gff","droppedOverlaps.gff");
    my %predictions;
    my %RFPredictions;
    my $hpARPM = getHairpinARPM($hairpinsFile);
    #retrieving predictions and scores from $predHairpinClassFile
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
    #print "\n\nretrieving information about  major products:\n";
    my($mpInfo,$miRProducts) = getMajorProductInfo($parameters->{hairpinsFile},$parameters->{productFile});
    my $annotations = annotatePredictions(\%predictions,$parameters);
    #classifying predictions and dropping overlaping folds with the same major product
    my($posPredictions,$negPredictions,
       $sizeFilteredPred,$overlaps) = classifyPredictions(\%predictions,\%RFPredictions,$hpARPM,$dropOverlaps,$mpInfo,$parameters);
    if ($dropHairpinOverlaps) {
	#dropping hairpins which have overlapping miR products
	($posPredictions,$overlaps) = dropOverlappingHairpins($posPredictions,$overlaps,$miRProducts,\%RFPredictions);
    }
    if ($parameters->{postPredExcludeFile}) {
	print "Dropping Positive Predictions in Repeat Regions\n";
	$posPredictions = dropRepeatRegions($posPredictions,$repeatRegions,$annotations,$mpInfo);
    }
    print "printing gff\'s\n";
    open(PGFF,">$posGff") or die "failed to open $posGff for writing\n";
    foreach my $chrom (%{$posPredictions}) {
	foreach my $hpInfo (@{$posPredictions->{$chrom}}) {
	    my($start,$stop,$strand,$geneId,$name) = @{$hpInfo};
	    my $RFavg = $RFPredictions{$name};
#	    my $readCount = $hpCounts->{$name};
	    my $hpInfo = "ID=$name;Alias=$name;Name=$name";
	    if($annotations->{$name}) {
		my $mirName = $annotations->{$name};
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
#	    my $readCount = $hpCounts->{$name};
	    my $hpInfo = "ID=$name;Alias=$name;Name=$name";
	    if($annotations->{$name}) {
		my $mirName = $annotations->{$name};
		$hpInfo = "ID=$name;Alias=$mirName;Name=$name";
	    }
	    print NGFF "$chrom\tmiRWOODS\tmiRNA_primary_transcript\t$start\t$stop\t$RFavg\t$strand\t.\t$hpInfo\n";
	}
    }
    close(NGFF);
    open(SFGFF,">$sizeFilteredGff") or die "failed to open $sizeFilteredGff for writing\n";
    foreach my $chrom (%{$sizeFilteredPred}) {
	foreach my $hpInfo (@{$sizeFilteredPred->{$chrom}}) {
	    my($start,$stop,$strand,$geneId,$name) = @{$hpInfo};
	    my $RFavg = $RFPredictions{$name};
#	    my $readCount = $hpCounts->{$name};
	    my $hpInfo = "ID=$name;Alias=$name;Name=$name";
	    if($annotations->{$name}) {
		my $mirName = $annotations->{$name};
		$hpInfo = "ID=$name;Alias=$mirName;Name=$name";
	    }
	    print SFGFF "$chrom\tmiRWOODS\tmiRNA_primary_transcript\t$start\t$stop\t$RFavg\t$strand\t.\t$hpInfo\n";
	}
    }
    close(SFGFF);
    #the following is for testing purposes and will eventually be included in the negGFF predicitons
    if ($dropOverlaps) {
	open(OGFF,">$overlapGff") or die "failed to open overlappingPosPredictions.gff for writing\n";
	foreach my $chrom (%{$overlaps}) {
	    foreach my $hpInfo (@{$overlaps->{$chrom}}) {
		my($start,$stop,$strand,$geneId,$name) = @{$hpInfo};
		my $RFavg = $RFPredictions{$name};
#		my $readCount = $hpCounts->{$name};
		my $hpInfo = "ID=$name;Alias=$name;Name=$name";
		if($annotations->{$name}) {
		    my $mirName = $annotations->{$name};
		    $hpInfo = "ID=$name;Alias=$mirName;Name=$name";
		}
		print OGFF "$chrom\tmiRWOODS\tmiRNA_primary_transcript\t$start\t$stop\t$RFavg\t$strand\t.\t$hpInfo\n";
	    }
	}
    close(OGFF);
    }
}

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

sub determineMiROverlap {
    my($hpScoreInfo,$usedScores,$overlapCutoff,$miRProducts) = @_;
    my $greatestOverlapAmt = 0;
    my $greatestOverlappingHP = 0;
    my($start1,$stop1,$hpScoreInfo1,$geneId1,$name1) = @{$hpScoreInfo};
    foreach my $productInfo1 (@{$miRProducts->{$name1}}) {
	my($prodName1,$side1,$type1,$productCount1,$maxProductCount1,$adjProductCount1,$adjMaxProductCount1,
	   $relProdStart1,$relProdStop1,$productStrand1,$productSequence1) = @{$productInfo1};
	my($prodStart1,$prodStop1) = getProductGenomicStartStop($start1,$stop1,$relProdStart1,$relProdStop1,$productStrand1);
	foreach my $scoreInfo (@{$usedScores}) {
	    my($start2,$stop2,$hpScoreInfo2,$geneId2,$name2) = @{$scoreInfo};
	    foreach my $productInfo2 (@{$miRProducts->{$name2}}) {
		my($prodName2,$side2,$type2,$productCount2,$maxProductCount2,$adjProductCount2,$adjMaxProductCount2,
		   $relProdStart2,$relProdStop2,$productStrand2,$productSequence2) = @{$productInfo2};
		my($prodStart2,$prodStop2) = getProductGenomicStartStop($start2,$stop2,$relProdStart2,$relProdStop2,$productStrand2);
		my $overlapAmt = getOverlap($prodStart1,$prodStop1,$prodStart2,$prodStop2);
		if ($overlapAmt > $overlapCutoff && $overlapAmt > $greatestOverlapAmt && $productStrand1 eq $productStrand2) {
		    $greatestOverlapAmt = $overlapAmt;
		    $greatestOverlappingHP = $name2;
		}
	    }	
	}
    }
    return $greatestOverlappingHP;
}

sub dropOverlappingHairpins {
    my($posPredictions,$overlaps,$miRProducts,$RFPredictions) = @_;
    my %newPosPredictions;
    my $overlapCutoff = 1;  #require at least this many nucleotides to overlap before dropping
    foreach my $chrom (keys %{$posPredictions}) {
	my %hpScores;
	#converting predictions hash to a hash keyed by strand so that we determine overlaps on the same strand
	foreach my $predictionInfo (@{$posPredictions->{$chrom}}) {
	    my($start,$stop,$strand,$geneId,$name) = @{$predictionInfo};
	    push(@{$hpScores{$strand}}, [$start,$stop,$RFPredictions->{$name},$geneId,$name]);
	}
	#now that the scores are keyed by strand
	#we are determining overlapping hairpins and seperating predictions into positives and negatives
	foreach my $strand (keys %hpScores) {
	    my @usedScores;
	    my @sortedHPScores = sort {$b->[2] <=> $a->[2]} @{$hpScores{$strand}};
	    foreach my $hpScoreInfo (@sortedHPScores) {
		my($start,$stop,$hpScore,$geneId,$name) = @{$hpScoreInfo};
		my $hairpinMiROverlap = determineMiROverlap($hpScoreInfo,\@usedScores,$overlapCutoff,$miRProducts);
		if ($hairpinMiROverlap) {
		    print "dropping $name because a miR it overlaps another miR on $hairpinMiROverlap\n";
		    push(@{$overlaps->{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
		} else {
		    push(@{$newPosPredictions{$chrom}}, [$start,$stop,$strand,$geneId,$name]);
		}
		push(@usedScores, [$start,$stop,$hpScore,$geneId,$name]);
	    }
	}
    }
    return(\%newPosPredictions,$overlaps);
}

sub getProductGenomicStartStop {
    my($start,$stop,$relProdStart,$relProdStop,$strand) = @_;
    my($prodStart,$prodStop) = (-1,-2);  #initialized with dummy variables for testing
    if ($strand eq '+') {
	$prodStart = $start + $relProdStart;
	$prodStop = $start + $relProdStop;
    } elsif ($strand eq '-') {
	$prodStop = $stop - $relProdStart;
	$prodStart = $stop - $relProdStop;
    }
    if ($prodStart > $prodStop) {
	die "prodStart is greater than prodStop\n";
    }
    return($prodStart,$prodStop);
}

sub getGlobalStartStop {
    my($location,$relStart,$relStop,$relStrand) = @_;
    my($chrom,$start,$stop) = parseLocation($location);
    my($gStart,$gStop) = (-1,-2);  #(-1,-2) to produce error if something goes wrong.
    if ($relStrand eq '+') {
	$gStart = $start + $relStart;
	$gStop = $start + $relStop;
    } elsif ($relStrand eq '-') {
	$gStop = $stop - $relStart;
	$gStart = $stop - $relStop;
    }
    if ($gStart > $gStop) {
	die "gStart is greater than gStop ($gStart > $gStop) in getGlobalStartStop:\nloc=$location\trelStart=$relStart\trelStop=$relStop\trelStrand=$relStrand\n";
    } 
    return($gStart,$gStop);
}

sub getRelStartStop {
    #converts gStart and gStop into a relative start and stop within a region bound by $location
    my($location,$gStart,$gStop,$gStrand) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $relStart = ($strand eq "+") ? $gStart - $start : $stop - $gStop;
    my $relStop = ($strand eq "+") ? $gStop - $start : $stop - $gStart;
    return ($relStart,$relStop);
}


sub getBestAnnot {
    my($precAnnotFile) = @_;
    my %annotations;
    my %bestAnnotated;
    open(PAF,$precAnnotFile) or die "failed to open $precAnnotFile\n";
    while (<PAF>) {
	chomp;
	unless ( /^#/ ) {
	    my($tag,$location,$annotId,$annotName,$hpCoverage,$annotCoverage,$annotMatCoverage,$annotMatCount) = split(/\t/);
	    push(@{$annotations{$annotName}},["$tag\_$location",$hpCoverage,$annotCoverage,$annotMatCoverage,$annotMatCount]);
	}
    }
    foreach my $annotName (keys %annotations) {
	my @sortedAnnots = sort {$b->[3] <=> $a->[3] || $b->[1] <=> $a->[1] || $b->[2] <=> $a->[2]} @{$annotations{$annotName}};
#	print $annotName . "\n";
#	foreach my $annotInfo (@sortedAnnots) {
#	    my($geneId,$hpCoverage,$annotCoverage,$annotMatCoverage,$annotMatCount) = @{$annotInfo};
#	    print "$geneId\t$hpCoverage\t$annotCoverage\t$annotMatCoverage\t$annotMatCount\n";
#	}
#	print "\n";
	my($geneId,$hpCoverage,$annotCoverage,$annotMatCoverage,$annotMatCount) = @{$sortedAnnots[0]};
	$bestAnnotated{$geneId} = [$hpCoverage,$annotCoverage,$annotMatCoverage,$annotMatCount];
    }
    return \%bestAnnotated;
}

sub printHairpinTrainingData {
    my($scoresFile,$gff,$genomeDir,$chromLengths,$class,$trainFile,$productFile,$parameters) = @_;
    my $minHPCoverage = 0.9;
    my $minAnnotCoverage = 0;
    my $minAnnotMatCoverage = 1;
    my $hpAnnotFile = $parameters->{hpAnnotFile} or die "Error: parameter hpAnnotFile not set";   
    unless (-e $hpAnnotFile) {
	die "error:" . $hpAnnotFile . " not found";
    }
    my $bestAnnotated = getBestAnnot($hpAnnotFile);
    open(HVF,">".$trainFile) or die "failed to open $trainFile for writing\n";
    open(HPS,$scoresFile) or die "failed to open $scoresFile for reading\n";
    my $first = 1;
    while (<HPS>) {
	chomp;
	my $line = $_;
	unless($first) {
	    my($geneId,$class,@other) = split(/\t/,$line);
#	    print $geneId . "\n" if ($bestAnnotated->{$geneId});
	    if ($bestAnnotated->{$geneId}) {
		my($hpCoverage,$annotCoverage,$annotMatCoverage,$annotMatCount) = @{$bestAnnotated->{$geneId}};
		if ($annotMatCount > 1 && $annotMatCoverage >= $minAnnotMatCoverage) {
		    print HVF "$geneId\t1\t" . join("\t", @other)  . "\n";
		} elsif (($annotMatCount == 1) && ($hpCoverage >= $minHPCoverage) && 
			 ($annotMatCoverage >= $minAnnotMatCoverage) && ($annotCoverage >= $minAnnotCoverage)) {
		    print HVF "$geneId\t1\t" . join("\t", @other)  . "\n";
		}
 	    }
	} else {
	    $first = 0;
	    print HVF $line . "\n";
	}	
    }
    close(HPS);
    close(HVF);
}

sub addRandomNegHPFeatures {
    #adds random negative features from feature vector file to the training data
    my($hairpinVecTrainFile,$hairpinScoresFile) = @_;
    my %posNames;
    my $posHairpinCount = 0;
    my @negVectorData;
    my $negHairpinMult = 1;
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
    $negHairpinCount *= $negHairpinMult;
    for (my $i = 0; $i < $negHairpinCount; $i++) {
	print HVFA pop(@shuffledNegVectorData) . "\n";
    }
    close(HVFA);
}

sub printHairpinVectorFile {
    my($hairpinsFile,$hairpinVectorFile,$genomeDir,$chromLengths,$class,$neighborCounts,$parameters) = @_;
    my $useHPWithMirProdOnly = ($parameters->{useHPWithMirProdOnly}) ? $parameters->{useHPWithMirProdOnly} : 0;
    my $hpWithMirProducts = getTagsWithMirProducts($parameters);
    open(HVF,">".$hairpinVectorFile) or die "failed to open $hairpinVectorFile for writing\n";
    print HVF "GeneId\ty\ttotalSense\ttotalAntisense\tmfe\taapd\ttapd\turf\tahc\tafh\tpbp\tsameShift\tbothShift\tSPA\tOPA\t";
    print HVF "aa\tac\tag\tat\tca\tcc\tcg\tct\tga\tgc\tgg\tgt\tta\ttc\ttg\ttt\tduplexEnergy\tGCcontent\t";
    #print HVF "zScore\t";
    print HVF "foldDupCmp\texactFoldDupCmp\tdupPBP\tdupLoopLength";
    print HVF "\tAPV\twAPV";
    print HVF "\tARV\twARV";
    print HVF "\tmpCount\tdupCount\tdupOverlap\tmpLoopDistance\tdupLoopDistance\tloopSize\tmpOverlapScore\tdupOverlapScore\twghtMPOverlapIntensity\twghtDupOverlapIntensity\tOPA2\ttotalOverlapAmount\taverageOverlapAmount\ttotalRelativeOverlapAmount\taverageRelativeOverlapAmount\tmaxOverlap\tinnerLoopGapCount\ttotalSenseRPM\ttotalAntisenseRPM";
    print HVF "\toverlapSizeRatioSum\ttotalOverlap\tmaxBulge\tmaxInteriorLoop\tintLoopSideDiff\tmaxUnboundOverhang\tnumOffshoots\tdupSize";
    print HVF "\trel5pOutCount\trel5pMorCount\trel5pMirCount\trelLoopCount\trelSplitCount\trel3pMirCount\trel3pMorCount\trel3pOutCount\textraOutCount";
    print HVF "\tmiRmoR5pOverlap\tmiRmoR3pOverlap\tmiRLoop5pOverlap\tmiRLoop3pOverlap\tloopLoopOverlap\tout5pOverlap\tout3pOverlap\toutOut5pOverlap\toutOut3pOverlap\tinProdOverlap\tmiRSplit5pOverlap\tmiRSplit3pOverlap\tmoRSplit5pOverlap\tmoRSplit3pOverlap\tsplitLoopOverlap\tloopSplitOverlap";
    print HVF "\tneighborCount";
    print HVF "\tRFProductAvg\n";
    open(HPF, $hairpinsFile) or die "failed to open $hairpinsFile\n";
    while (<HPF>) {
        chomp;
	unless(/\#/) {	  
	    my($name,$chrom,$start,$stop,$strand,$leftCenter,$rightCenter,$totalSense,$totalAntisense,$mfe,$sequence,$fold,$aapd,$tapd,$urf,$ahc,$afh,$pbp,$sameShift,$bothShift,$spa,$opa,$aa,$ac,$ag,$at,$ca,$cc,$cg,$ct,$ga,$gc,$gg,$gt,$ta,$tc,$tg,$tt,$duplexEnergy,$GCcontent,$foldDupCmp,$exactFoldDupCmp,$dupPBP,$dupLoopLength,$apv,$wApv,$arv,$wArv,$maxProductCount,$mpDuplexCount,$duplexOverlap,$mpLoopDistance,$dupLoopDistance,$loopSize,$mpOverlapScore,$dupOverlapScore,$wghtMPOverlapIntensity,$wghtDupOverlapIntensity,$opa2,$totalOverlapAmount,$averageOverlapAmount,$totalRelativeOverlapAmount,$averageRelativeOverlapAmount,$maxOverlap,$innerLoopGapCount,$totalSenseRPM,$totalAntisenseRPM,$overlapSizeRatioSum,$totalOverlap,$maxBulge,$maxInteriorLoop,$intLoopSideDiff,$maxUnboundOverhang,$numOffshoots,$dupSize,$rel5pOutCount,$rel5pMorCount,$rel5pMirCount,$relLoopCount,$relSplitCount,$rel3pMirCount,$rel3pMorCount,$rel3pOutCount,$extraOutCount,$miRmoR5pOverlap,$miRmoR3pOverlap,$miRLoop5pOverlap,$miRLoop3pOverlap,$loopLoopOverlap,$out5pOverlap,$out3pOverlap,$outOut5pOverlap,$outOut3pOverlap,$inProdOverlap,$miRSplit5pOverlap,$miRSplit3pOverlap,$moRSplit5pOverlap,$moRSplit3pOverlap,$splitLoopOverlap,$loopSplitOverlap,$RFProductAvg) = split(/\t/,$_);
	    if (($useHPWithMirProdOnly && $hpWithMirProducts->{$name}) || not $useHPWithMirProdOnly) {
		my $neighbors = $neighborCounts->{$name};
		my $location = "$chrom:$start..$stop:$strand";
		print HVF "$name" . "_" . "$location\t$class\t$totalSense\t$totalAntisense\t$mfe\t$aapd\t$tapd\t$urf\t$ahc\t$afh\t$pbp\t$sameShift\t$bothShift\t$spa\t$opa\t";
		print HVF "$aa\t$ac\t$ag\t$at\t$ca\t$cc\t$cg\t$ct\t$ga\t$gc\t$gg\t$gt\t$ta\t$tc\t$tg\t$tt\t";
		print HVF "$duplexEnergy\t$GCcontent";
		print HVF "\t$foldDupCmp\t$exactFoldDupCmp\t$dupPBP\t$dupLoopLength";
		print HVF "\t$apv\t$wApv\t$arv\t$wArv";
		print HVF "\t$maxProductCount\t$mpDuplexCount\t$duplexOverlap\t$mpLoopDistance\t$dupLoopDistance\t$loopSize\t$mpOverlapScore\t$dupOverlapScore\t$wghtMPOverlapIntensity\t$wghtDupOverlapIntensity\t$opa2\t$totalOverlapAmount\t$averageOverlapAmount\t$totalRelativeOverlapAmount\t$averageRelativeOverlapAmount\t$maxOverlap\t$innerLoopGapCount\t$totalSenseRPM\t$totalAntisenseRPM";
		print HVF "\t$overlapSizeRatioSum\t$totalOverlap\t$maxBulge\t$maxInteriorLoop\t$intLoopSideDiff\t$maxUnboundOverhang\t$numOffshoots\t$dupSize";
		print HVF "\t$rel5pOutCount\t$rel5pMorCount\t$rel5pMirCount\t$relLoopCount\t$relSplitCount\t$rel3pMirCount\t$rel3pMorCount\t$rel3pOutCount\t$extraOutCount\t$miRmoR5pOverlap\t$miRmoR3pOverlap\t$miRLoop5pOverlap\t$miRLoop3pOverlap\t$loopLoopOverlap\t$out5pOverlap\t$out3pOverlap\t$outOut5pOverlap\t$outOut3pOverlap\t$inProdOverlap\t$miRSplit5pOverlap\t$miRSplit3pOverlap\t$moRSplit5pOverlap\t$moRSplit3pOverlap\t$splitLoopOverlap\t$loopSplitOverlap";
		print HVF "\t$neighbors";
		print HVF "\t$RFProductAvg";
		print HVF "\n";
	    }
	}
    }
    close(HPF);
    close(HVF);
}

sub printHairpinVectorFile_old {
    my($hairpinsFile,$hairpinVectorFile,$genomeDir,$chromLengths,$class,$neighborCounts,$parameters) = @_;
    my $useHPWithMirProdOnly = ($parameters->{useHPWithMirProdOnly}) ? $parameters->{useHPWithMirProdOnly} : 0;
    my $hpWithMirProducts = getTagsWithMirProducts($parameters);
    open(HVF,">".$hairpinVectorFile) or die "failed to open $hairpinVectorFile for writing\n";
    print HVF "GeneId\ty\ttotalSense\ttotalAntisense\tmfe\taapd\ttapd\turf\tahc\tafh\tpbp\tsameShift\tbothShift\tSPA\tOPA\t";
    print HVF "aa\tac\tag\tat\tca\tcc\tcg\tct\tga\tgc\tgg\tgt\tta\ttc\ttg\ttt\tduplexEnergy\tGCcontent\t";
    #print HVF "zScore\t";
    print HVF "foldDupCmp\texactFoldDupCmp\tdupPBP\tdupLoopLength";
    print HVF "\tAPV\twAPV";
    print HVF "\tARV\twARV";
    print HVF "\tmpCount\tdupCount\tdupOverlap\tmpLoopDistance\tdupLoopDistance\tloopSize\tmpOverlapScore\tdupOverlapScore\twghtMPOverlapIntensity\twghtDupOverlapIntensity\tOPA2\ttotalOverlapAmount\taverageOverlapAmount\ttotalRelativeOverlapAmount\taverageRelativeOverlapAmount\tmaxOverlap\tinnerLoopGapCount\ttotalSenseRPM\ttotalAntisenseRPM";
    print HVF "\toverlapSizeRatioSum\ttotalOverlap\tmaxBulge\tmaxInteriorLoop\tintLoopSideDiff\tmaxUnboundOverhang\tnumOffshoots\tdupSize";
    print HVF "\tneighborCount";
    print HVF "\tRFProductAvg\n";
    open(HPF, $hairpinsFile) or die "failed to open $hairpinsFile\n";
    while (<HPF>) {
        chomp;
	unless(/\#/) {	  
	    my($name,$chrom,$start,$stop,$strand,$leftCenter,$rightCenter,$totalSense,$totalAntisense,$mfe,$sequence,$fold,$aapd,$tapd,$urf,$ahc,$afh,$pbp,$sameShift,$bothShift,$spa,$opa,$aa,$ac,$ag,$at,$ca,$cc,$cg,$ct,$ga,$gc,$gg,$gt,$ta,$tc,$tg,$tt,$duplexEnergy,$GCcontent,$foldDupCmp,$exactFoldDupCmp,$dupPBP,$dupLoopLength,$apv,$wApv,$arv,$wArv,$maxProductCount,$mpDuplexCount,$duplexOverlap,$mpLoopDistance,$dupLoopDistance,$loopSize,$mpOverlapScore,$dupOverlapScore,$wghtMPOverlapIntensity,$wghtDupOverlapIntensity,$opa2,$totalOverlapAmount,$averageOverlapAmount,$totalRelativeOverlapAmount,$averageRelativeOverlapAmount,$maxOverlap,$innerLoopGapCount,$totalSenseRPM,$totalAntisenseRPM,$overlapSizeRatioSum,$totalOverlap,$maxBulge,$maxInteriorLoop,$intLoopSideDiff,$maxUnboundOverhang,$numOffshoots,$dupSize,$RFProductAvg) = split(/\t/,$_);
	    if (($useHPWithMirProdOnly && $hpWithMirProducts->{$name}) || not $useHPWithMirProdOnly) {
		my $neighbors = $neighborCounts->{$name};
		my $location = "$chrom:$start..$stop:$strand";
		print HVF "$name" . "_" . "$location\t$class\t$totalSense\t$totalAntisense\t$mfe\t$aapd\t$tapd\t$urf\t$ahc\t$afh\t$pbp\t$sameShift\t$bothShift\t$spa\t$opa\t";
		print HVF "$aa\t$ac\t$ag\t$at\t$ca\t$cc\t$cg\t$ct\t$ga\t$gc\t$gg\t$gt\t$ta\t$tc\t$tg\t$tt\t";
		print HVF "$duplexEnergy\t$GCcontent";
		print HVF "\t$foldDupCmp\t$exactFoldDupCmp\t$dupPBP\t$dupLoopLength";
		print HVF "\t$apv\t$wApv\t$arv\t$wArv";
		print HVF "\t$maxProductCount\t$mpDuplexCount\t$duplexOverlap\t$mpLoopDistance\t$dupLoopDistance\t$loopSize\t$mpOverlapScore\t$dupOverlapScore\t$wghtMPOverlapIntensity\t$wghtDupOverlapIntensity\t$opa2\t$totalOverlapAmount\t$averageOverlapAmount\t$totalRelativeOverlapAmount\t$averageRelativeOverlapAmount\t$maxOverlap\t$innerLoopGapCount\t$totalSenseRPM\t$totalAntisenseRPM";
		print HVF "\t$overlapSizeRatioSum\t$totalOverlap\t$maxBulge\t$maxInteriorLoop\t$intLoopSideDiff\t$maxUnboundOverhang\t$numOffshoots\t$dupSize";
		print HVF "\t$neighbors";
		print HVF "\t$RFProductAvg";
		print HVF "\n";
	    }
	}
    }
    close(HPF);
    close(HVF);
}

sub getTagsWithMirProducts {
    my($parameters) = @_;
    my %hpWithMirProducts;
    my $productFile = $parameters->{productFile};
    open(PROD,$productFile) or die "failed to open $productFile\n";
    while (<PROD>) {
	unless ( /^#/ ) {
	    my($tag,$side,$type) = split(/\t/);
	    if ($type eq 'miR') {
		$hpWithMirProducts{$tag} = 1;
	    }
	}
    }
    close(PROD);
    return \%hpWithMirProducts;
}

#####################
# INPUT SUBROUTINES #
#####################

sub getColNum {
    #getColNum retrieves the column number for a header entry.  
    my($header,$colName) = @_;
    my @headerNames = split(/\t/,$header);
    for (my $colNum = 0; $colNum < @headerNames; $colNum++) {
	if ($headerNames[$colNum] eq $colName) {
	    return $colNum;
	}
    }
    die "failed to find column number for $colName\nColumns = ". join(",",@headerNames) ."\n";
}

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
    my $fileType = checkFileFormat($repeatRegionsFile,$chromSizes);
    if ($fileType eq "fileList") {
	open(RRFL,$repeatRegionsFile) or die "could not open $repeatRegionsFile\n";
	while(<RRFL>) {
	    chomp;
	    my $repeatRegionsFile = $_;
	    my $listedFileType = checkFileFormat($repeatRegionsFile,$chromSizes);
	    $repeatRegions = readRepeatRegionsFile($repeatRegionsFile,$repeatRegions,$listedFileType);
	}
	close(RRFL);
    } else {
	$repeatRegions = readRepeatRegionsFile($repeatRegionsFile,$repeatRegions,$fileType);
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

sub checkFileFormat {
    my $fileName = shift;
    my $chromSizes = shift;
    my @chroms = keys %{$chromSizes};
    my $fileType;

    if (($fileType) = $fileName =~ /\.(bed|gff)$/) {
	return $fileType;
    }
    open(FP, $fileName) or die "failed to open $fileName in checkFileFormat";
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

######################################
#test functions                      #
######################################

sub printProducts {
    my($products) = @_;
    my @sortedProducts = sort {$a->[7] <=> $b->[7]} @{$products};
    foreach my $product (@sortedProducts) {
	my($side,$newType,$productList,$productCount,$maxProductCount,$adjProductCount,$adjMaxProductCount,$relStart,$relStop) = @{$product};
	my $length = $relStop - $relStart + 1;
	print ' ' x $relStart . 'x' x $length . ' ' x 5 . "\(adjProdCount = $adjProductCount , length = $length\n"; 
    }
}

sub insideOffshoot {
    my($position,$offshoots) = @_;
    foreach my $offshoot (@{$offshoots}) {
	if ($offshoot->[0] <= $position && $offshoot->[1] >= $position) {
	    return 1;
	}
    }
    return 0;
}

sub printMPAndDup {
    my($fold,$maxProductCoords) = @_;
    my $map = pairMap($fold);
    my $offshoots = getOffshoots($map,$maxProductCoords);
    my $duplexCoords = getDuplexCoords($map,$maxProductCoords,$offshoots);
    my $dupOffshoots = getOffshoots($map,$duplexCoords);
    print "$fold\n";
    print ' ' x $maxProductCoords->[0] . 'x' x ($maxProductCoords->[1] - $maxProductCoords->[0] + 1) . "\n";
    if ($duplexCoords->[0] != -1 && $duplexCoords->[1] != -1) {
	print ' ' x $duplexCoords->[0] . 'x' x ($duplexCoords->[1] - $duplexCoords->[0] + 1) . "\n";
    } 
    if (@{$offshoots}) {
	print "Major Product Offshoots:\n";
	foreach my $offshoot (@{$offshoots}) {
	    print ' ' x $offshoot->[0] . 'x' x ($offshoot->[1] - $offshoot->[0] + 1) . "\n";	    
	}
    }
    if (@{$dupOffshoots}) {
	print "Dup Offshoots:\n";
	foreach my $offshoot (@{$dupOffshoots}) {
	    print ' ' x $offshoot->[0] . 'x' x ($offshoot->[1] - $offshoot->[0] + 1) . "\n";	    
	}
    }
}



1;
__END__


=head1 NAME

miRWoods - Perl module (and script) for sensitive detection of microRNAs using Stacked Random Forests

=head1 SYNOPSIS

miRWoods is a software using a stacked random forest strategy for the sensitive detection microRNAs, including those with only one read.

=head1 DESCRIPTION

miRWoods uses a stacked random forest strategy to detect microRNAs.  While most softwares will place limits on the minimum number of reads within potential miR loci, miRWoods does not.  miRWoods uses two random forests as part of it's strategy.  The first random forest referred to as the miR Product Random Forsest (MPRF) assembles read stacks into products and scores them.  Products which don't pass through the MRPF are filtered out, and those which do are combined with surrounding products and folded to produce hairpins.  Each hairpin is scored and then passed through a Hairpin Precursor Random Forest (HPRF) to obtain the final results.  miRWoods will evaluate several possible overlapping precursors for each loci, and pick the one with the best score.

=head1 SEE ALSO

http://hendrixlab.cgrb.oregonstate.edu/

=head1 AUTHORS

Jimmy Bell, bellji@oregonstate.edu

Dr. David Hendrix, david.hendrix@oregonstate.edu 

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2019 by Jimmy Bell

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

=cut

