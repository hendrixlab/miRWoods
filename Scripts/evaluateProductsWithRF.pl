#!/usr/bin/perl -w
use Memory::Usage;
use miRWoods;
use Bio::DB::Sam;
use Getopt::Long;
use strict;
$| = 1;

my $USAGE = "USAGE:\n$0 -L <config file>\n";

die $USAGE unless (@ARGV);

my $parameters = miRWoods::loadDefaultParameters();

Getopt::Long::Configure("no_ignore_case");

GetOptions ("lengthMin=i" => \$parameters->{lengthMin},
	    "totalLength=i" => \$parameters->{totalLength},
	    "hairpinShortLength=i" => \$parameters->{hairpinShortLength},
	    "distanceMin=i" => \$parameters->{distanceMin},
	    "reverseMax=i" => \$parameters->{reverseMax},
	    "countMinLocus=i" => \$parameters->{countMinLocus},
	    "fivePrimeHetMax=f" => \$parameters->{fivePrimeHetMax},
	    "shiftMin=i" => \$parameters->{shiftMin},
	    "OverlapMin=i" => \$parameters->{OverlapMin},
	    "InHairpinBuffer=i" => \$parameters->{InHairpinBuffer},
	    "OutHairpinBuffer=i" => \$parameters->{OutHairpinBuffer},
	    "RangeOfHairpin=i" => \$parameters->{RangeOfHairpin},
	    "outputPrefix=s" => \$parameters->{outputPrefix},
	    "bamListFile=s" => \$parameters->{bamListFile},
	    "FileRepeatRegions=s" => \$parameters->{FileRepeatRegions},
	    "genomeDir=s" => \$parameters->{genomeDir},
	    "SizesFile=s" => \$parameters->{SizesFile},
	    "LoadFromConfigFile=s" => \$parameters->{LoadFromConfigFile}
    );

unless ($parameters->{LoadFromConfigFile} eq "") {
    my $configFile = $parameters->{LoadFromConfigFile};
    miRWoods::readConfigFile($configFile,$parameters);
}
miRWoods::createOutputFileParameters($parameters);

my $prodMLParameters = miRWoods::loadDefaultMLProductParameters() or die "failed to load product machine learning parameters\n";

#checking file Parameters
my $productBedFile = $parameters->{predProductFile} or die "error: predProductFile parameter not set\n";
my $productRF = $parameters->{productRF} or die "error: productRF parameter not set\n";
my $productFeatVectorFile = $parameters->{productFeatVectorFile} or die "error: productFeatVectorFile parameter not set\n";
my $predProductClassFile = $parameters->{predProductClasses} or die "error: predProductClasses parameter not set\n";
die "productFeatVectorFile not found\n" unless (-e $productFeatVectorFile);

#scoresToRemove is an array of scores not to include if they are present in the train file or vector file
my $scoresToRemove = miRWoods::createRemovedScoresSet($prodMLParameters);

#starting R
my $R = Statistics::R->new();
$R->start();

#training the product random forest model if trainModels is set to 1
#if ($parameters->{trainModels}) {
#    my $productTrainFile = $parameters->{productTrainFile} or die "error: productTrainFile parameter not set";
#    die "$productTrainFile not found.\n" unless (-e $productTrainFile);
#    miRWoods::createRFModel($R,$productRF,$productTrainFile,$scoresToRemove,$parameters);
#}
print "$parameters->{trainModels} = ". $parameters->{trainModels} . "\n";
if ($parameters->{trainModels}) {
    if ($parameters->{prodMTry}) {
	$parameters->{mTry} = $parameters->{prodMTry}; #the random forest function uses $parameters->{mTry} to train on.
    } else {
	$parameters->{mTry} = 0;  #cross validation will be used to determine mtry
    }
    my $productTrainFile = $parameters->{productTrainFile} or die "error: productTrainFile parameter not set\n";
    die "$productTrainFile not found.\n" unless (-e $productTrainFile);
    my $prodTrainMethod = $parameters->{prodTrainMethod};
    my($numPos,$numNeg) = countClasses($productTrainFile);
    print "\nClass Counts:\n\t\'1\' = $numPos\n\t\'-1\' = $numNeg\n";
    if ($prodTrainMethod eq "weightedRF") {
	print "Train Method = Weighted Random Forest\n";
	my($posWeight,$negWeight) = (0,0);
	if ($parameters->{prodPosWeight} && $parameters->{prodNegWeight}) {
	    $posWeight = $parameters->{prodPosWeight};
	    $negWeight = $parameters->{prodNegWeight};
	} elsif ($parameters->{prodPosWeight}) {
	    $posWeight = $parameters->{prodPosWeight};
	    $negWeight = 1 - $posWeight; 
	} elsif ($parameters->{prodNegWeight}) {
	    $negWeight = $parameters->{prodNegWeight};
	    $posWeight = 1 - $negWeight;
	} else {
	    my $weightShift = ($parameters->{prodWeightShift}) ? $parameters->{prodWeightShift} : 0;
	    if ($weightShift) {
		print "shifting the weights of the pos and neg classes by $weightShift\n";
	    }
	    $negWeight = $numPos / ($numPos + $numNeg) - $weightShift;
	    $posWeight = $numNeg / ($numPos + $numNeg) + $weightShift;
	}
	print "Class Weights:\n\t\'1\' = $posWeight\n\t\'-1\' = $negWeight\n";
	miRWoods::createRFModel_classwt($R,$productRF,$productTrainFile,$scoresToRemove,$posWeight,$negWeight,$parameters);
    } elsif ($prodTrainMethod eq "wghtStratSampling") {
	print "Train Method = Weighted Stratified Sampling\n";
	my($stratSampPos,$stratSampNeg) = (0,0);
	if($parameters->{prodStratSampPos} && $parameters->{prodStratSampNeg}) {
	    $stratSampPos = $parameters->{prodStratSampPos};
	    $stratSampNeg = $parameters->{prodStratSampNeg};
	} else {
	    my $stratSampMultiplier = ($parameters->{prodStratSampMult}) ? $parameters->{prodStratSampMult} : 1;
	    unless ($stratSampMultiplier == 1) {
		print "using $stratSampMultiplier times the number of the positive class for the stratified sample size of the negative class\n";
	    }
	    $stratSampPos = $numPos;
	    $stratSampNeg = $numPos * $stratSampMultiplier;
	}
	my($posWeight,$negWeight) = (0,0);
	if ($parameters->{prodPosWeight} && $parameters->{prodNegWeight}) {
	    $posWeight = $parameters->{prodPosWeight};
	    $negWeight = $parameters->{prodNegWeight};
	} elsif ($parameters->{prodPosWeight}) {
	    $posWeight = $parameters->{prodPosWeight};
	    $negWeight = 1 - $posWeight; 
	} elsif ($parameters->{prodNegWeight}) {
	    $negWeight = $parameters->{prodNegWeight};
	    $posWeight = 1 - $negWeight;
	} else {
	    my $weightShift = ($parameters->{prodWeightShift}) ? $parameters->{prodWeightShift} : 0;
	    if ($weightShift) {
		print "shifting the weights of the pos and neg classes by $weightShift\n";
	    }
	    $negWeight = $stratSampPos / ($stratSampPos + $stratSampNeg) - $weightShift;
	    $posWeight = $stratSampNeg / ($stratSampPos + $stratSampNeg) + $weightShift;
	}
	print "Stratified Sample Count:\n\t\'1\' = $stratSampPos\n\t\'-1\' = $stratSampNeg\n";
	print "Class Weights:\n\t\'1\' = $posWeight\n\t\'-1\' = $negWeight\n";
	miRWoods::createRFModel_wgtdStratSampling($R,$productRF,$productTrainFile,$scoresToRemove,$stratSampPos,$stratSampNeg,$posWeight,$negWeight,$parameters);
    } elsif ($prodTrainMethod eq "stratifiedSampling") {
	print "Train Method = Stratified Sampling\n";
	my($stratSampPos,$stratSampNeg) = (0,0);
	if($parameters->{prodStratSampPos} && $parameters->{prodStratSampNeg}) {
	    $stratSampPos = $parameters->{prodStratSampPos};
	    $stratSampNeg = $parameters->{prodStratSampNeg};
	} else {
	    my $stratSampMultiplier = ($parameters->{prodStratSampMult}) ? $parameters->{prodStratSampMult} : 1;
	    unless ($stratSampMultiplier == 1) {
		print "using $stratSampMultiplier times the number of the positive class for the stratified sample size of the negative class\n";
	    }
	    $stratSampPos = $numPos;
	    $stratSampNeg = $numPos * $stratSampMultiplier;
	}
	print "Stratified Sample Count:\n\t\'1\' = $stratSampPos\n\t\'-1\' = $stratSampNeg\n";
	miRWoods::createRFModel_sampSize($R,$productRF,$productTrainFile,$scoresToRemove,$stratSampPos,$stratSampNeg,$parameters);
    } else {
	print "Train Method = normal\n";      
	miRWoods::createRFModel($R,$productRF,$productTrainFile,$scoresToRemove,$parameters);
    }
}


#using the model to determine scores for each product in the product vector file
miRWoods::RFModelPredict($R,$productRF,$productFeatVectorFile,$predProductClassFile,$scoresToRemove,$parameters);

#miRWoods::evaluateFeatureVectorFile_new($R,$productTrainFile,$productFeatVectorFile,$predProductClassFile,$scoresToRemove,$parameters);
$R->stop();
miRWoods::createProductBedFile($predProductClassFile,$parameters);



sub countClasses {
    my($hairpinTrainFile) = @_;
    my($numPos,$numNeg) = (0,0);
    open(HPTF,$hairpinTrainFile) or die "failed to open $hairpinTrainFile\n";
    my $header = <HPTF>;
    while (<HPTF>) {
	chomp;
	my($geneId,$class) = split(/\t/);
	if ($class == 1) {
	    $numPos++;
	} elsif ($class == -1) {
	    $numNeg++;
	} else {
	    die "unknown class (class = $class) in $hairpinTrainFile\n";
	}
    }
    close(HPTF);
    return($numPos,$numNeg);
}
