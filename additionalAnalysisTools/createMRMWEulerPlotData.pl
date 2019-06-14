#!/usr/bin/perl -w
use miRWoods;
use mwPaperCommonCode;
use strict;
$| = 1;

my $usage = "USAGE:\n$0 <miRWoods Pos Predictions> <mirwoods hairpins File > <miReap results> <mirbase gff> <expression File> <output Prefix>\n";
my $posPredGFF = $ARGV[0] or die $usage;
my $hairpinsFile = $ARGV[1] or die $usage;
my $miReapResultsFile = $ARGV[2] or die $usage;
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
#loading miReap Data
my($mrPredictions) = mwPaperCommonCode::parseMiReapResults($miReapResultsFile);
my $mrAnnotPredictions = mwPaperCommonCode::annotateMiReap($mrPredictions,$annoHairpins,$annoProducts,\$warning);
my($newMRAnnotPredictions,$mrDuplicatePredictions,$mrDuplicatedPredictionCount) = mwPaperCommonCode::checkForDuplicatePredictions($mrAnnotPredictions,\$warning);
my($mdKnown,$mdNovel) = mwPaperCommonCode::getKnownAndNovelPredictions($newMRAnnotPredictions,$parameters);
#loading annotation expression Data
my($hpExpressionInfo) = mwPaperCommonCode::readExpressionFile($expressionFile); 
my $hpExpressedCounts = mwPaperCommonCode::countExpressedAnnotations($hpExpressionInfo);
my($miRWoodsAnnotOnly,$miReapAnnotOnly,$bothAnnot) = mwPaperCommonCode::getEulerIntersectStats($mwKnown,$mdKnown);
my($miRWoodsNovelOnly,$miReapNovelOnly,$bothNovel) = mwPaperCommonCode::getEulerIntersectStats($mwNovel,$mdNovel);
printOverlappingInfo($miRWoodsAnnotOnly,$miReapAnnotOnly,$bothAnnot,
		     $miRWoodsNovelOnly,$miReapNovelOnly,$bothNovel,$overlappingInfoFile);
printEulerInputFile($miRWoodsAnnotOnly,$miReapAnnotOnly,$bothAnnot,$miRWoodsNovelOnly,
		    $miReapNovelOnly,$bothNovel,$hpExpressedCounts,$eulerInputFile);

sub printOverlappingInfo {
    my($miRWoodsAnnotOnly,$miReapAnnotOnly,$bothAnnot,$miRWoodsNovelOnly,$miReapNovelOnly,$bothNovel,$overlappingInfoFile) = @_;
    open(OIF,">$overlappingInfoFile") or die "failed to open $overlappingInfoFile for writing\n";
    print OIF "Overlapping: miRWoods - miReap - Annotations\n";
    foreach my $entry (@{$bothAnnot}) {
	my($miRWoodsId,$miReapId) = @{$entry};
	print OIF "$miRWoodsId\t$miReapId\n"
    }
    print OIF "\nOverlapping: miRWoods - Annotations\n";
    foreach my $name (@{$miRWoodsAnnotOnly}) {
	print OIF "$name\n";
    }
    print OIF "\nOverlapping: miReap - Annotations\n";
    foreach my $name (@{$miReapAnnotOnly}) {
	print OIF "$name\n";
    }
    print OIF "\nOverlapping: miRWoods - miReap\n";
    foreach my $entry (@{$bothNovel}) {
	my($miRWoodsId,$miReapId) = @{$entry};
	print OIF "$miRWoodsId\t$miReapId\n"
    }
    print OIF "\nOverlapping: miRWoods\n";
    foreach my $name (@{$miRWoodsNovelOnly}) {
	print OIF "$name\n";
    }
    print OIF "\nOverlapping: miReap\n";
    foreach my $name (@{$miReapNovelOnly}) {
	print OIF "$name\n";
    }    
    close(OIF);
}

sub printEulerInputFile {
    my($miRWoodsAnnotOnly,$miReapAnnotOnly,$bothAnnot,$miRWoodsNovelOnly,$miReapNovelOnly,$bothNovel,$hpExpressedCount,$eulerInputFile) = @_;
    my $numAnnotPredByBoth = @{$bothAnnot};
    my $numAnnotPredByMirWoods = @{$miRWoodsAnnotOnly};
    my $numAnnotPredByMirDeep = @{$miReapAnnotOnly};
    my $numNovelPredByBoth = @{$bothNovel};
    my $numNovelPredByMirWoods = @{$miRWoodsNovelOnly};
    my $numNovelPredByMirDeep = @{$miReapNovelOnly};
    my $numNonPredictedAnnotations = $hpExpressedCount - $numAnnotPredByBoth - $numAnnotPredByMirWoods - $numAnnotPredByMirDeep;
    open(EUL,">$eulerInputFile") or die "failed to open $eulerInputFile for writing\n";
    print EUL "miRWoods\tmiReap\tannotations\tmiRWoods_miReap\tmiRWoods_annotations\tmiReap_annotations\tmiRWoods_miReap_annotations\n";
    print EUL "$numNovelPredByMirWoods\t$numNovelPredByMirDeep\t$numNonPredictedAnnotations\t$numNovelPredByBoth\t$numAnnotPredByMirWoods\t$numAnnotPredByMirDeep\t$numAnnotPredByBoth\n";
    close(EUL);
}
