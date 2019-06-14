#!/usr/bin/perl -w
use miRWoods;
use mwPaperCommonCode;
use strict;

my $usage = "USAGE:\n$0 <miReap results gff> <mirbase gff> <sample name> <output prefix> <get annotated Only ? 1 : 0>\n";

my $miReapPredGff = $ARGV[0] or die $usage;
my $mirbaseGff = $ARGV[1] or die $usage;
my $sampleName = $ARGV[2] or die $usage;
my $outputPrefix = $ARGV[3] or die $usage;
my $annotOnly = $ARGV[4];


my $outputFile = ($annotOnly) ? "$outputPrefix\_annotProductSizes\_miReap\_counts.txt" : "$outputPrefix\_majorProductSizes\_miReap\_counts.txt";


my($newMiReapAnnotPredictions,$miReapProducts) = getMiReapAnnotations($miReapPredGff,$mirbaseGff);
my $mpInfo = getMajorProductCounts($miReapPredGff);
printResults($newMiReapAnnotPredictions,$mpInfo,$annotOnly,$outputFile);

sub printResults {
    my($newMiReapAnnotPredictions,$mpInfo,$annotOnly,$outputFile) = @_;
    open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
    print OPTF "#tag\tannotationType\tproduct\t$sampleName\n";
    foreach my $chrom (keys %{$newMiReapAnnotPredictions}) {
	foreach my $hpInfo (@{$newMiReapAnnotPredictions->{$chrom}}) {
	    my($hpStart,$hpStop,$hpStrand,$hpId,$annoName,$annoType) = @{$hpInfo};
	    if ($annoType eq "Annotated" || !$annotOnly) {
		my($prodName,$matureCount) = @{$mpInfo->{$hpId}};
		print OPTF "$hpId\t$annoType\t$prodName\t$matureCount\n";
	    }
	}
    }
    close(OPTF);
}


sub getMiReapAnnotations {
    my($miReapPredGff,$mirbaseGff) = @_;
    my $warning = "";
    my $parameters = mwPaperCommonCode::loadDefaultParameters();
    my($mrPredictions,$miReapProducts) = mwPaperCommonCode::parseMiReapResults($miReapPredGff);
    my($annoHairpins,$annoProducts) = miRWoods::readMirbaseGff3($mirbaseGff);
    my $mrAnnotPredictions = mwPaperCommonCode::annotateMiReap($mrPredictions,$annoHairpins,$annoProducts,\$warning);
    my($newMiReapAnnotPredictions,$duplicatePredictions,$duplicatedPredictionCount) = mwPaperCommonCode::checkForDuplicatePredictions($mrAnnotPredictions,
																\$warning);
    print $warning . "\n";
    return($newMiReapAnnotPredictions,$miReapProducts);
}


sub getMajorProductCounts {
    my($miReapPredGff) = @_;
    my $mpInfo = {};
    my($hairpins,$products) = readMiReapGff($miReapPredGff);
    foreach my $chrom (keys %{$hairpins}) {
	foreach my $hairpinInfo (@{$hairpins->{$chrom}}) {
	    my($hpStart,$hpStop,$hpStrand,$hpId) = @{$hairpinInfo};
	    my @sortedProducts = sort {$b->[5] <=> $a->[5]} @{$products->{$hpId}};
	    my($mpChrom,$mpStart,$mpStop,$mpStrand,$mpId,$mpCount) = @{$sortedProducts[0]};
	    @{$mpInfo->{$hpId}} = ($mpId,$mpCount);
	}
    }
    return $mpInfo;
}


sub readMiReapGff {
    my($miReapPredGff) = @_;
    my %hairpins;
    my %products;
    open(MBGFF3,$miReapPredGff) or die "could not open $miReapPredGff\n";
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
	    if($type eq "precursor") {
		# hairpin region from the gff file
		my $id = $info{ID} or die "No ID found for the line:\n$_\n"; 
		my $count = $info{Count} or die "No Count found for the line:\n$_\n"; 
		push(@{$hairpins{$chrom}},[$start,$stop,$strand,$id]);
	    }
	    if($type eq "mature-5p" || $type eq "mature-3p") {
		# mature product from the hairpin file
		my $id = $info{ID} or die "No ID found for the line:\n$_\n";
                my $count = $info{Count} or die "No Count found for the line:\n$_\n";
		my $parentId = $info{Parent} or die "No Derives_from found for the line:\n$_\n";
		push(@{$products{$parentId}},[$chrom,$start,$stop,$strand,$id,$count]);
	    }
	}
    }
    close(MBGFF3);
    return(\%hairpins,\%products);
}
