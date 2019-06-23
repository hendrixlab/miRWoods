#!/usr/bin/perl -w
use miRWoods;
use Bio::DB::Sam;
use Getopt::Long;
use Statistics::R;
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
	    "stratSampMult=i" => \$parameters->{stratSampMult},
	    "stratSampPos=i" => \$parameters->{stratSampPos},
	    "stratSampNeg=i" => \$parameters->{stratSampNeg},
	    "hairpinRF=s" => \$parameters->{hairpinRF},
	    "LoadFromConfigFile=s" => \$parameters->{LoadFromConfigFile}
    );
unless ($parameters->{LoadFromConfigFile} eq "") {
    my $configFile = $parameters->{LoadFromConfigFile};
    miRWoods::readConfigFile($configFile,$parameters);
}
miRWoods::createOutputFileParameters($parameters);

#my $HPMLParameters = miRWoods::loadDefaultMLHairpinParameters();
my $HPMLParameters = miRWoods::loadDefaultMLHairpinParameters_prodOverlaps();

#checking file Parameters
my $hairpinRF = $parameters->{hairpinRF} or die "error: hairpinRF parameter not set";
my $productFile = $parameters->{productFile} or die "error: productFile parameter not set";
my $hairpinsFile = $parameters->{hairpinsFile} or die "error: hairpinsFile parameter not set\n";
my $distinctReadsFile = $parameters->{distinctReadsFile} or die "error: distinctReadsFile parameter not set\n";
my $hairpinVectorFile = $parameters->{hairpinVectorFile} or die "error: hairpinVectorFile parameter not set\n";
my $predHairpinClassFile = $parameters->{predHairpinClasses} or die "error: predHairpinClasses parameter not set\n";
die "$productFile not found\n" unless (-e $productFile);
die "$hairpinsFile not found\n" unless (-e $hairpinsFile);
die "$distinctReadsFile not found\n" unless (-e $distinctReadsFile);
die "$hairpinVectorFile not found\n" unless (-e $hairpinVectorFile);

#loading repeat regions if a repeat regions file exists
my $repeatRegions;
if ($parameters->{postPredExcludeFile}) {
    my $repeatRegionsFile = $parameters->{postPredExcludeFile};
    my $sizesFile = $parameters->{SizesFile} or die "error: SizesFile parameter not set\n";
    die "$repeatRegionsFile not found\n" unless (-e $repeatRegionsFile);
    die "$sizesFile not found\n" unless (-e $sizesFile);
    my $chromLengths = miRWoods::readChromLengths($parameters->{SizesFile});
    $repeatRegions = miRWoods::loadRepeatRegions($repeatRegionsFile, $chromLengths) or die "Failed to load repeat regions file.\n";
}

#starting R
my $R = Statistics::R->new();
print "starting R\n";
$R->start();
print "R is started\n";

#scoresToRemove is an array of scores not to include if they are present in the train file or vector file
my $scoresToRemove = miRWoods::createRemovedScoresSet($HPMLParameters);

#training the hairpin random forest model if trainModels is set to 1
print "$parameters->{trainModels} = ". $parameters->{trainModels} . "\n";
if ($parameters->{trainModels}) {
    my $hairpinTrainFile = $parameters->{hairpinTrainFile} or die "error: hairpinTrainFile parameter not set\n";
    die "$hairpinTrainFile not found.\n" unless (-e $hairpinTrainFile);
    my $trainMethod = $parameters->{trainMethod};
    my($numPos,$numNeg) = countClasses($hairpinTrainFile);
    print "\nClass Counts:\n\t\'1\' = $numPos\n\t\'-1\' = $numNeg\n";
    if ($trainMethod eq "weightedRF") {
	print "Train Method = Weighted Random Forest\n";
	my($posWeight,$negWeight) = (0,0);
	if ($parameters->{posWeight} && $parameters->{negWeight}) {
	    $posWeight = $parameters->{posWeight};
	    $negWeight = $parameters->{negWeight};
	} elsif ($parameters->{posWeight}) {
	    $posWeight = $parameters->{posWeight};
	    $negWeight = 1 - $posWeight; 
	} elsif ($parameters->{negWeight}) {
	    $negWeight = $parameters->{negWeight};
	    $posWeight = 1 - $negWeight;
	} else {
	    my $weightShift = ($parameters->{weightShift}) ? $parameters->{weightShift} : 0;
	    if ($weightShift) {
		print "shifting the weights of the pos and neg classes by $weightShift\n";
	    }
	    $negWeight = $numPos / ($numPos + $numNeg) - $weightShift;
	    $posWeight = $numNeg / ($numPos + $numNeg) + $weightShift;
	}
	print "Class Weights:\n\t\'1\' = $posWeight\n\t\'-1\' = $negWeight\n";
	miRWoods::createRFModel_classwt($R,$hairpinRF,$hairpinTrainFile,$scoresToRemove,$posWeight,$negWeight,$parameters);
    } elsif ($trainMethod eq "wghtStratSampling") {
	print "Train Method = Weighted Stratified Sampling\n";
	my($stratSampPos,$stratSampNeg) = (0,0);
	if($parameters->{stratSampPos} && $parameters->{stratSampNeg}) {
	    $stratSampPos = $parameters->{stratSampPos};
	    $stratSampNeg = $parameters->{stratSampNeg};
	} else {
	    my $stratSampMultiplier = ($parameters->{stratSampMult}) ? $parameters->{stratSampMult} : 1;
	    unless ($stratSampMultiplier == 1) {
		print "using $stratSampMultiplier times the number of the positive class for the stratified sample size of the negative class\n";
	    }
	    $stratSampPos = $numPos;
	    $stratSampNeg = $numPos * $stratSampMultiplier;
	}
	my($posWeight,$negWeight) = (0,0);
	if ($parameters->{posWeight} && $parameters->{negWeight}) {
	    $posWeight = $parameters->{posWeight};
	    $negWeight = $parameters->{negWeight};
	} elsif ($parameters->{posWeight}) {
	    $posWeight = $parameters->{posWeight};
	    $negWeight = 1 - $posWeight; 
	} elsif ($parameters->{negWeight}) {
	    $negWeight = $parameters->{negWeight};
	    $posWeight = 1 - $negWeight;
	} else {
	    my $weightShift = ($parameters->{weightShift}) ? $parameters->{weightShift} : 0;
	    if ($weightShift) {
		print "shifting the weights of the pos and neg classes by $weightShift\n";
	    }
	    $negWeight = $stratSampPos / ($stratSampPos + $stratSampNeg) - $weightShift;
	    $posWeight = $stratSampNeg / ($stratSampPos + $stratSampNeg) + $weightShift;
	}
	print "Stratified Sample Count:\n\t\'1\' = $stratSampPos\n\t\'-1\' = $stratSampNeg\n";
	print "Class Weights:\n\t\'1\' = $posWeight\n\t\'-1\' = $negWeight\n";
	miRWoods::createRFModel_wgtdStratSampling($R,$hairpinRF,$hairpinTrainFile,$scoresToRemove,$stratSampPos,$stratSampNeg,$posWeight,$negWeight,$parameters);
    } elsif ($trainMethod eq "stratifiedSampling") {
	print "Train Method = Stratified Sampling\n";
	my($stratSampPos,$stratSampNeg) = (0,0);
	if($parameters->{stratSampPos} && $parameters->{stratSampNeg}) {
	    $stratSampPos = $parameters->{stratSampPos};
	    $stratSampNeg = $parameters->{stratSampNeg};
	} else {
	    my $stratSampMultiplier = ($parameters->{stratSampMult}) ? $parameters->{stratSampMult} : 1;
	    unless ($stratSampMultiplier == 1) {
		print "using $stratSampMultiplier times the number of the positive class for the stratified sample size of the negative class\n";
	    }
	    $stratSampPos = $numPos;
	    $stratSampNeg = $numPos * $stratSampMultiplier;
	}
	print "Stratified Sample Count:\n\t\'1\' = $stratSampPos\n\t\'-1\' = $stratSampNeg\n";
	miRWoods::createRFModel_sampSize($R,$hairpinRF,$hairpinTrainFile,$scoresToRemove,$stratSampPos,$stratSampNeg,$parameters);
    } else {
	print "Train Method = Normal\n";
	miRWoods::createRFModel($R,$hairpinRF,$hairpinTrainFile,$scoresToRemove,$parameters);
    }
}
die "$hairpinRF file not found.\n" unless (-e $hairpinRF);
print "model = $hairpinRF\n";

#using the model to determine scores for each hairpin in the hairpin vector file
miRWoods::RFModelPredict($R,$hairpinRF,$hairpinVectorFile,$predHairpinClassFile,$scoresToRemove,$parameters);

$R->stop();
miRWoods::printPredictedMirGff($predHairpinClassFile,$productFile,$hairpinsFile,$repeatRegions,$parameters);
miRWoods::extractGffReadRegions("positivePredictions",$hairpinsFile,$productFile,$distinctReadsFile);
miRWoods::extractGffReadRegions("negativePredictions",$hairpinsFile,$productFile,$distinctReadsFile);
miRWoods::extractGffReadRegions("sizeFiltered",$hairpinsFile,$productFile,$distinctReadsFile);


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
