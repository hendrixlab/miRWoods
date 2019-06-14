#!/usr/bin/perl -w
use miRWoods;
use mwPaperCommonCode;
use strict;
$| = 1;

my $usage = "USAGE:\n$0 <miRWoods Pos Predictions> <mirwoods hairpins File > <mirbase gff> <expression File> <output Prefix>\n";
my $posPredGFF = $ARGV[0] or die $usage;
my $hairpinsFile = $ARGV[1] or die $usage;
my $mirbaseGff = $ARGV[2] or die $usage;
my $expressionFile = $ARGV[3] or die $usage;
my $outputPrefix = $ARGV[4] or die $usage;

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
#loading annotation expression Data
my($hpExpressionInfo) = mwPaperCommonCode::readExpressionFile($expressionFile); 
my $hpExpressedCounts = mwPaperCommonCode::countExpressedAnnotations($hpExpressionInfo);
printEulerInputFile_annot($mwKnown,$mwNovel,$hpExpressedCounts,$eulerInputFile);

sub printEulerInputFile_annot {
    my($mwKnown,$mwNovel,$hpExpressedCount,$eulerInputFile) = @_;
    my $numAnnotPredByMirWoods = 0;  
    my $numNovelPredByMirWoods = 0;
    print "Annot:\n";
    foreach my $chrom (keys %{$mwKnown}) {
	foreach my $hpInfo (@{$mwKnown->{$chrom}}) {
	    my($start,$stop,$strand,$id,$annotation,$annotationType) = @{$hpInfo};
	    print "$annotation\t$annotationType\n";
	    $numAnnotPredByMirWoods++;
	}
    }
    print "\n\n\nNov:\n";
    foreach my $chrom (keys %{$mwNovel}) {
	foreach my $hpInfo (@{$mwNovel->{$chrom}}) {
	    my($start,$stop,$strand,$id,$annotation,$annotationType) = @{$hpInfo};
	    print "$annotation\t$annotationType\n";
	    $numNovelPredByMirWoods++;
	}
    }
    my $numNonPredictedAnnotations = $hpExpressedCount - $numAnnotPredByMirWoods;
    open(EUL,">$eulerInputFile") or die "failed to open $eulerInputFile for writing\n";
    print EUL "miRWoods\tannotations\tmiRWoods_annotations\n";
    print EUL "$numNovelPredByMirWoods\t$numNonPredictedAnnotations\t$numAnnotPredByMirWoods\n";
    close(EUL);
}
