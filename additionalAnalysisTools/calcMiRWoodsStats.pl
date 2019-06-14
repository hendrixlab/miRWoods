#!/usr/bin/perl -w
use miRWoods;
use mwPaperCommonCode;
use strict;

my $usage = "USAGE:\n$0 <mirwoods predictions gff> <miRWoods predictions hairpinsFile> <mirbase gff> <expression File> <output File> [<homology file>]\n";
my $posPredGFF = $ARGV[0] or die $usage;
my $hairpinsFile = $ARGV[1] or die $usage;
my $mirbaseGff = $ARGV[2] or die $usage;
my $expressionFile = $ARGV[3] or die $usage;
my $outputPrefix = $ARGV[4] or die $usage;
my $homologyFile = $ARGV[5];                       #can be zero for no homology file

my $outputFile = $outputPrefix . ".txt";
my $infoFile = $outputPrefix . "_info.txt";

my $warning = "";

my $parameters = mwPaperCommonCode::loadDefaultParameters();

my $mwPredictions = mwPaperCommonCode::parseMiRWoodsResults($posPredGFF);
my $mwFolds = mwPaperCommonCode::getFolds($hairpinsFile);
my($annoHairpins,$annoProducts) = miRWoods::readMirbaseGff3($mirbaseGff);
my $mwAnnotPredictions = mwPaperCommonCode::annotateMirWoods($mwPredictions,$mwFolds,$annoHairpins,$annoProducts,\$warning);
if ($homologyFile) {
    $mwAnnotPredictions = mwPaperCommonCode::addMiRWoodsHomology($mwAnnotPredictions,$homologyFile);
}
my($newAnnotPredictions,$duplicatePredictions,$duplicatedPredictionCount) = mwPaperCommonCode::checkForDuplicatePredictions($mwAnnotPredictions,\$warning);
my($numAnnotated,$numNovel,$numAntisense,$numOverlaps,$numHomologs) = mwPaperCommonCode::countAnnotationTypes($newAnnotPredictions);
my($hpExpressionInfo) = mwPaperCommonCode::readExpressionFile($expressionFile); 
my $totalExpressed = mwPaperCommonCode::countExpressedAnnotations($hpExpressionInfo);
my($hpCount,$novelCount,$totalPos,$precision,$recall,$f1Score) = mwPaperCommonCode::calcPrecisionAndRecall($numAnnotated,$numNovel,$numAntisense,$numOverlaps,$numHomologs,$totalExpressed,$parameters);
mwPaperCommonCode::printInfoFile($infoFile,$mwAnnotPredictions,$warning);

print $warning;

open(OTF,">$outputFile") or die "failed to open $outputFile for writing\n";
print OTF "miRs in data = $totalPos\n";
print OTF "True Positives = $hpCount\n";
print OTF "False Positives = $novelCount\n";
if ($homologyFile) {
    print OTF "Homology Count = $numHomologs\n";
}
print OTF "Novel Count = " . ($novelCount - $numHomologs) . "\n";
print OTF "precision = $precision\n";
print OTF "recall = $recall\n";
print OTF "f1Score = $f1Score\n";
print OTF "Annotated Count = $numAnnotated\n";
print OTF "Other Novel = $numNovel\n";
print OTF "Antisense Count = $numAntisense\n";
print OTF "Overlap Count = $numOverlaps\n";
print OTF "Total = " . ($hpCount + $novelCount) . "\n";
close(OTF);
