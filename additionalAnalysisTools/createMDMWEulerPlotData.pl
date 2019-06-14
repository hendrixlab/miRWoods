#!/usr/bin/perl -w
use miRWoods;
use mwPaperCommonCode;
use strict;
$| = 1;

my $usage = "USAGE:\n$0 <miRWoods Pos Predictions> <mirwoods hairpins File > <mirdeep results> <mirbase gff> <expression File> <output Prefix>\n";
my $posPredGFF = $ARGV[0] or die $usage;
my $hairpinsFile = $ARGV[1] or die $usage;
my $mirDeepResultsFile = $ARGV[2] or die $usage;
my $mirbaseGff = $ARGV[3] or die $usage;
my $expressionFile = $ARGV[4] or die $usage;
my $outputPrefix = $ARGV[5] or die $usage;
my $eulerInputFile = $outputPrefix . '_eulerPlotData.txt';
my $overlappingInfoFile = $outputPrefix . '_eulerPlotInfo.txt';
my $warningsOutput = $outputPrefix . '_warnings.txt';

my $warning = "";

my $parameters = mwPaperCommonCode::loadDefaultParameters();

#loading miRWoods Data
my $mwPredictions = mwPaperCommonCode::parseMiRWoodsResults($posPredGFF);
my $mwFolds = mwPaperCommonCode::getFolds($hairpinsFile);
my($annoHairpins,$annoProducts) = miRWoods::readMirbaseGff3($mirbaseGff);
my $mwAnnotPredictions = mwPaperCommonCode::annotateMirWoods($mwPredictions,$mwFolds,$annoHairpins,$annoProducts,\$warning);
my($newMWAnnotPredictions,
   $mwDuplicatePredictions,$mwDuplicatedPredictionCount) = mwPaperCommonCode::checkForDuplicatePredictions($mwAnnotPredictions,\$warning);
my($mwKnown,$mwNovel) = mwPaperCommonCode::getKnownAndNovelPredictions($newMWAnnotPredictions,$parameters);
#loading miRDeep Data
my($mirDeepAccData,$mdPredictions) = mwPaperCommonCode::parseMirDeepResults($mirDeepResultsFile);
$mdPredictions = mwPaperCommonCode::applyMirdeepCutoff($mirDeepAccData,$mdPredictions,$parameters);
my($mdAnnotPredictions) = mwPaperCommonCode::annotateMirDeep($mdPredictions,$annoHairpins,$annoProducts,\$warning);
my($newMDAnnotPredictions,
   $mdDuplicatePredictions,$mdDuplicatedPredictionCount) = mwPaperCommonCode::checkForDuplicatePredictions($mdAnnotPredictions,\$warning);
my($mdKnown,$mdNovel) = mwPaperCommonCode::getKnownAndNovelPredictions($newMDAnnotPredictions,$parameters);
#loading annotation expression Data
my($hpExpressionInfo) = mwPaperCommonCode::readExpressionFile($expressionFile); 
my $hpExpressedCounts = mwPaperCommonCode::countExpressedAnnotations($hpExpressionInfo);
my($miRWoodsAnnotOnly,$miRDeepAnnotOnly,$bothAnnot) = mwPaperCommonCode::getEulerIntersectStats($mwKnown,$mdKnown);
my($miRWoodsNovelOnly,$miRDeepNovelOnly,$bothNovel) = mwPaperCommonCode::getEulerIntersectStats($mwNovel,$mdNovel);
printOverlappingInfo($miRWoodsAnnotOnly,$miRDeepAnnotOnly,$bothAnnot,
		     $miRWoodsNovelOnly,$miRDeepNovelOnly,$bothNovel,$overlappingInfoFile);
printEulerInputFile($miRWoodsAnnotOnly,$miRDeepAnnotOnly,$bothAnnot,$miRWoodsNovelOnly,
		    $miRDeepNovelOnly,$bothNovel,$hpExpressedCounts,$eulerInputFile);

sub printOverlappingInfo {
    my($miRWoodsAnnotOnly,$miRDeepAnnotOnly,$bothAnnot,$miRWoodsNovelOnly,$miRDeepNovelOnly,$bothNovel,$overlappingInfoFile) = @_;
    open(OIF,">$overlappingInfoFile") or die "failed to open $overlappingInfoFile for writing\n";
    print OIF "Overlapping: miRWoods - miRDeep - Annotations\n";
    foreach my $entry (@{$bothAnnot}) {
	my($miRWoodsId,$miRDeepId) = @{$entry};
	print OIF "$miRWoodsId\t$miRDeepId\n"
    }
    print OIF "\nOverlapping: miRWoods - Annotations\n";
    foreach my $name (@{$miRWoodsAnnotOnly}) {
	print OIF "$name\n";
    }
    print OIF "\nOverlapping: miRDeep - Annotations\n";
    foreach my $name (@{$miRDeepAnnotOnly}) {
	print OIF "$name\n";
    }
    print OIF "\nOverlapping: miRWoods - miRDeep\n";
    foreach my $entry (@{$bothNovel}) {
	my($miRWoodsId,$miRDeepId) = @{$entry};
	print OIF "$miRWoodsId\t$miRDeepId\n"
    }
    print OIF "\nOverlapping: miRWoods\n";
    foreach my $name (@{$miRWoodsNovelOnly}) {
	print OIF "$name\n";
    }
    print OIF "\nOverlapping: miRDeep\n";
    foreach my $name (@{$miRDeepNovelOnly}) {
	print OIF "$name\n";
    }    
    close(OIF);
}

sub printEulerInputFile {
    my($miRWoodsAnnotOnly,$miRDeepAnnotOnly,$bothAnnot,$miRWoodsNovelOnly,$miRDeepNovelOnly,$bothNovel,$hpExpressedCount,$eulerInputFile) = @_;
    my $numAnnotPredByBoth = @{$bothAnnot};
    my $numAnnotPredByMirWoods = @{$miRWoodsAnnotOnly};
    my $numAnnotPredByMirDeep = @{$miRDeepAnnotOnly};
    my $numNovelPredByBoth = @{$bothNovel};
    my $numNovelPredByMirWoods = @{$miRWoodsNovelOnly};
    my $numNovelPredByMirDeep = @{$miRDeepNovelOnly};
    my $numNonPredictedAnnotations = $hpExpressedCount - $numAnnotPredByBoth - $numAnnotPredByMirWoods - $numAnnotPredByMirDeep;
    open(EUL,">$eulerInputFile") or die "failed to open $eulerInputFile for writing\n";
    print EUL "miRWoods\tmiRDeep\tannotations\tmiRWoods_miRDeep\tmiRWoods_annotations\tmiRDeep_annotations\tmiRWoods_miRDeep_annotations\n";
    print EUL "$numNovelPredByMirWoods\t$numNovelPredByMirDeep\t$numNonPredictedAnnotations\t$numNovelPredByBoth\t$numAnnotPredByMirWoods\t$numAnnotPredByMirDeep\t$numAnnotPredByBoth\n";
    close(EUL);
}
