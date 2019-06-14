#!/usr/bin/perl -w
use miRWoods;
use mwPaperCommonCode;
use strict;


my $usage = "USAGE:\n$0 <mireap predictions gff> <mirbase gff> <expression File> <output File> [<homology File>]\n";
my $posPredGFF = $ARGV[0] or die $usage;
my $mirbaseGff = $ARGV[1] or die $usage;
my $expressionFile = $ARGV[2] or die $usage;
my $outputPrefix = $ARGV[3] or die $usage;
my $homologyFile = $ARGV[4];

my $outputFile = $outputPrefix . ".txt";
my $infoFile = $outputPrefix . "_info.txt";

my $warning = "";

my $parameters = mwPaperCommonCode::loadDefaultParameters();

my($mrPredictions) = mwPaperCommonCode::parseMiReapResults($posPredGFF);
my($annoHairpins,$annoProducts) = miRWoods::readMirbaseGff3($mirbaseGff);
my $mrAnnotPredictions = mwPaperCommonCode::annotateMiReap($mrPredictions,$annoHairpins,$annoProducts,\$warning);
if ($homologyFile) {
    $mrAnnotPredictions = mwPaperCommonCode::addMiReapHomology($mrAnnotPredictions,$homologyFile);
}
my($newAnnotPredictions,$duplicatePredictions,$duplicatedPredictionCount) = mwPaperCommonCode::checkForDuplicatePredictions($mrAnnotPredictions,\$warning);
my($numAnnotated,$numNovel,$numAntisense,$numOverlaps,$numHomologs) = mwPaperCommonCode::countAnnotationTypes($newAnnotPredictions);
my($hpExpressionInfo) = mwPaperCommonCode::readExpressionFile($expressionFile); 
my $totalExpressed = mwPaperCommonCode::countExpressedAnnotations($hpExpressionInfo);
my($hpCount,$novelCount,$totalPos,$precision,$recall,$f1Score) = mwPaperCommonCode::calcPrecisionAndRecall($numAnnotated,$numNovel,$numAntisense,$numOverlaps,$numHomologs,$totalExpressed,$parameters);
mwPaperCommonCode::printInfoFile($infoFile,$mrAnnotPredictions,$warning);

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
