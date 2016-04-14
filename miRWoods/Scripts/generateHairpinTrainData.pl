#!/usr/bin/perl -w
use Memory::Usage;
use miRWoods;
use Bio::DB::Sam;
use Getopt::Long;
use File::Copy;
use strict;
$| = 1;

my $USAGE = "USAGE:\n$0 -L <config file>";

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

my $mirBaseGff = $parameters->{mirbaseGff} or die "mirBase gff not loaded.\n";
my($posClass,$negClass) = miRWoods::getClassValues($parameters);
my $hairpinVectorFile = $parameters->{hairpinVectorFile};
my $hairpinVecTrainFile = $parameters->{hairpinTrainFile};
my $balancedHairpinVecTrainFile = $parameters->{hairpinTrainFileSB};
my $prodBalancedHPVecTrainFile = "readRegions_hairpinTrainFilePB.txt";
my $productFile = $parameters->{productFile};
my $chromSizesFile = $parameters->{SizesFile};
my $genomeDir = $parameters->{genomeDir};
my $chromLengths = miRWoods::readChromLengths($chromSizesFile);
miRWoods::printHairpinTrainingData($hairpinVectorFile,$mirBaseGff,$genomeDir,$chromLengths,$posClass,$hairpinVecTrainFile,$productFile,$parameters);
copy($hairpinVecTrainFile,$balancedHairpinVecTrainFile) or die "Failed to copy $hairpinVecTrainFile to $balancedHairpinVecTrainFile\n";
copy($hairpinVecTrainFile,$prodBalancedHPVecTrainFile) or die "Failed to copy $hairpinVecTrainFile to $prodBalancedHPVecTrainFile";

if ($parameters->{negGff}) {
    my $curatedNegGff = $parameters->{negGff};
    miRWoods::printHairpinTrainingData($hairpinVectorFile,$curatedNegGff,$genomeDir,$chromLengths,$negClass,$hairpinVecTrainFile,$parameters);
} else {
    #negative products have not been curated so we will assume that most features from the productVectorFile are not true miRs and add from that
    miRWoods::addRandomNegHPFeatures($hairpinVecTrainFile,$hairpinVectorFile);
}

miRWoods::addRandomNegHPFeatures_SizeBalanced($balancedHairpinVecTrainFile,$hairpinVectorFile,$mirBaseGff,$parameters);
miRWoods::addRandomNegHPFeatures_ProdBalanced($prodBalancedHPVecTrainFile,$hairpinVectorFile,$mirBaseGff,$parameters);

