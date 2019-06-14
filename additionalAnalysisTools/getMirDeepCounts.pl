#!/usr/bin/perl -w
use miRWoods;
use mwPaperCommonCode;
use strict;


my $usage = "USAGE:\n$0 <miRDeepResultsFile> <mirbase gff> <sample name> <output prefix> <get annotated Only ? 1 : 0>\n";

my $miRDeepResultsFile = $ARGV[0] or die $usage;
my $mirbaseGff = $ARGV[1] or die $usage;
my $sampleName = $ARGV[2] or die $usage;
my $outputPrefix = $ARGV[3] or die $usage;
my $annotOnly = $ARGV[4];

my $outputFile = ($annotOnly) ? "$outputPrefix\_annotProductSizes\_miRDeep\_counts.txt" : "$outputPrefix\_majorProductSizes\_miRDeep\_counts.txt";



my($newMiRDeepAnnotPredictions,$mirDeepProducts,$mdPredictions) = getMiRDeepAnnotations($miRDeepResultsFile,$mirbaseGff);
my $mpInfo = getMajorProductCounts($mdPredictions,$mirDeepProducts);
printResults($newMiRDeepAnnotPredictions,$mpInfo,$annotOnly,$outputFile);

sub printResults {
    my($newMiRDeepAnnotPredictions,$mpInfo,$annotOnly,$outputFile) = @_;
    open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
    print OPTF "#tag\tannotationType\tproduct\t$sampleName\n";
    foreach my $chrom (keys %{$newMiRDeepAnnotPredictions}) {
	foreach my $hpInfo (@{$newMiRDeepAnnotPredictions->{$chrom}}) {
	    my($hpStart,$hpStop,$hpStrand,$hpId,$annoName,$annoType) = @{$hpInfo};
	    if ($annoType eq "Annotated" || !$annotOnly) {
		my($prodName,$matureCount) = @{$mpInfo->{$hpId}};
		print OPTF "$hpId\t$annoType\t$prodName\t$matureCount\n";
	    }
	}
    }
    close(OPTF);
}

sub getMiRDeepAnnotations {
    my($miRDeepResultsFile,$mirbaseGff) = @_;
    my $warning = "";
    my $parameters = mwPaperCommonCode::loadDefaultParameters();
    my($mirDeepAccData,$mdPredictions,$mirDeepProducts) = mwPaperCommonCode::parseMirDeepResults($miRDeepResultsFile);
    $mdPredictions = mwPaperCommonCode::applyMirdeepCutoff($mirDeepAccData,$mdPredictions,$parameters);
    my($annoHairpins,$annoProducts) = miRWoods::readMirbaseGff3($mirbaseGff);
    my($mdAnnotPredictions) = mwPaperCommonCode::annotateMirDeep($mdPredictions,$annoHairpins,$annoProducts,\$warning);
    my($newMiRDeepAnnotPredictions,$duplicatePredictions,$duplicatedPredictionCount) = mwPaperCommonCode::checkForDuplicatePredictions($mdAnnotPredictions,
																       \$warning);
    print $warning . "\n";
    return($newMiRDeepAnnotPredictions,$mirDeepProducts,$mdPredictions);
}

sub getMajorProductCounts {
    my($mdPredictions,$mirDeepProducts) = @_;
    my $mpInfo = {};
    foreach my $chrom (keys %{$mdPredictions}) {
	foreach my $miRDeepResult (@{$mdPredictions->{$chrom}}) {
	    my($id,$score,$TPProb,$rfamAlert,$total,$matureCount,$loopCount,$starCount,$randfoldSignificant,$miRNA,$otherSpecMiRNA,$UCSCbrowser,$NCBIblastn,$mature,$star,$precursor,$location,$mdPredictionClass) = @{$miRDeepResult};
	    foreach my $prodInfo (@{$mirDeepProducts->{$id}}) {
		my($chrom,$start,$stop,$strand,$prodName) = @{$prodInfo};
		if ($prodName =~ /-mat-\dp$/) {
		    @{$mpInfo->{$id}} = ($prodName,$matureCount);
		}
	    }
	}
    }
    return $mpInfo;
}
