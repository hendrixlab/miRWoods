#!/usr/bin/perl -w
use miRWoods;
use mwPaperCommonCode;
use strict;


my $usage = "USAGE:\n$0 <mirwoods predictions gff> <miRWoods predictions hairpinsFile> <miRWoods predictions productFile> <mirbase gff> <sample name> <output prefix> <get annotated Only ? 1 : 0>\n";

my $miRWoodsPredGff = $ARGV[0] or die $usage;
my $miRWoodsHairpinsFile = $ARGV[1] or die $usage;
my $miRWoodsProductFile = $ARGV[2] or die $usage;
my $mirbaseGff = $ARGV[3] or die $usage;
my $sampleName = $ARGV[4] or die $usage;
my $outputPrefix = $ARGV[5] or die $usage;
my $annotOnly = $ARGV[6];

my $outputFile = ($annotOnly) ? "$outputPrefix\_annotProductSizes\_miRWoods\_counts.txt" : "$outputPrefix\_majorProductSizes\_miRWoods\_counts.txt";

my($newMiRWoodsAnnotPredictions,$miRWoodsProducts) = getMiRWoodsAnnotations($miRWoodsPredGff,$miRWoodsHairpinsFile,$miRWoodsProductFile,$mirbaseGff);
my $mpInfo = getMajorProducts($miRWoodsProductFile);
printResults($newMiRWoodsAnnotPredictions,$mpInfo,$annotOnly,$outputFile);


sub printResults {
    my($newMiRWoodsAnnotPredictions,$mpInfo,$annotOnly,$outputFile) = @_;
    open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
    print OPTF "#tag\tannotationType\tproduct\t$sampleName\n";
    foreach my $chrom (keys %{$newMiRWoodsAnnotPredictions}) {
	foreach my $hpInfo (@{$newMiRWoodsAnnotPredictions->{$chrom}}) {
	    my($hpStart,$hpStop,$hpStrand,$hpId,$annoName,$annoType) = @{$hpInfo};
	    if ($annoType eq "Annotated" || !$annotOnly) {
		my($prodName,$matureCount) = @{$mpInfo->{$hpId}};
		print OPTF "$hpId\t$annoType\t$prodName\t$matureCount\n";
	    }
	}
    }
    close(OPTF);
}

sub getMiRWoodsAnnotations {
    my($miRWoodsPredGff,$miRWoodsHairpinsFile,$miRWoodsProductFile,$mirbaseGff) = @_;
    my $warning = "";
    my $parameters = mwPaperCommonCode::loadDefaultParameters();
    my $mwPredictions = mwPaperCommonCode::parseMiRWoodsResults($miRWoodsPredGff);
    my $miRWoodsProducts = mwPaperCommonCode::getMiRWoodsMatProducts($mwPredictions,$miRWoodsProductFile);
    my $mwFolds = mwPaperCommonCode::getFolds($miRWoodsHairpinsFile);
    my($annoHairpins,$annoProducts) = miRWoods::readMirbaseGff3($mirbaseGff);
    my $mwAnnotPredictions = mwPaperCommonCode::annotateMirWoods($mwPredictions,$mwFolds,$annoHairpins,$annoProducts,\$warning);
    my($newMiRWoodsAnnotPredictions,$duplicatePredictions,$duplicatedPredictionCount) = mwPaperCommonCode::checkForDuplicatePredictions($mwAnnotPredictions,
																	\$warning);
    print $warning . "\n";
    return($newMiRWoodsAnnotPredictions,$miRWoodsProducts);
}

sub getMajorProducts {
    my($miRWoodsProductFile) = @_;
    my $mpInfo = {};
    my $mwProducts = readProductFile($miRWoodsProductFile);
    foreach my $tag (keys %{$mwProducts}) {
	my @sortedProducts = sort {$b->[1] <=> $a->[1]} @{$mwProducts->{$tag}};
	my($mpName,$mpSize) = @{$sortedProducts[0]};
	@{$mpInfo->{$tag}} = ($mpName,$mpSize);
    }
    return $mpInfo;
}


sub readProductFile {
    my($miRWoodsProductFile) = @_;
    my $mwProducts = {};
    open(PF,$miRWoodsProductFile) or die "failed to open $miRWoodsProductFile";
    while (<PF>) {
	chomp;
	unless ( /^#/ ) {
	    my($tag,$side,$type,$total) = split;      
	    if (($type =~ /miR$/) && !($type =~ /as-/)) {
		my $product = "$tag-$side-$type";
		push(@{$mwProducts->{$tag}},[$product,$total]);
	    }
	}
    }
    close(PF);
    return $mwProducts;
}
