#!/usr/bin/perl -w
use miRWoods;
use Bio::DB::Sam;
use Getopt::Long;
use Statistics::R;
use RNA;
use List::Util;
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

my $scriptDir = $parameters->{scriptDir} ? $parameters->{scriptDir} : "";
my $checkMiRBaseGff = $scriptDir . "/checkMiRBaseGff.pl";
my $getLibrarySizes = $scriptDir . "/countMappedReads.pl";
my $printReadRegions = $scriptDir . "/printReadRegions.pl";
my $extractProductFeatures = $scriptDir . "/extractProductFeatures.pl";
my $generateProductTrainData = $scriptDir . "/generateProductTrainData.pl";
my $evaluateProductsWithRF = $scriptDir . "/evaluateProductsWithRF.pl";
my $processReadRegions = $scriptDir . "/processReadRegions.pl";
my $generateAnnotations = $scriptDir . "/generateAnnotations.pl";
my $extractHairpinFeatures = $scriptDir . "/extractHairpinFeatures.pl";
my $generateHairpinTrainData = $scriptDir . "/generateHairpinTrainData.pl";
my $evaluateHairpinsWithRF = $scriptDir . "/evaluateHairpinsWithRF.pl";
my $printMiRWoodsPredictions = $scriptDir . "/printMiRWoodsPredictions.pl";

die "error: $checkMiRBaseGff not found\n" unless (-e $checkMiRBaseGff);
die "error: $getLibrarySizes not found\n" unless (-e $getLibrarySizes);
die "error: $printReadRegions not found\n" unless (-e $printReadRegions);
die "error: $extractProductFeatures not found\n" unless (-e $extractProductFeatures);
die "error: $generateProductTrainData not found\n" unless (-e $generateProductTrainData);
die "error: $evaluateProductsWithRF not found\n" unless (-e $evaluateProductsWithRF);
die "error: $processReadRegions not found\n" unless (-e $processReadRegions);
die "error: $generateAnnotations not found\n" unless (-e $generateAnnotations);
die "error: $extractHairpinFeatures not found\n" unless (-e $extractHairpinFeatures);
die "error: $generateHairpinTrainData not found\n" unless (-e $generateHairpinTrainData);
die "error: $evaluateHairpinsWithRF not found\n" unless (-e $evaluateHairpinsWithRF);
die "error: $printMiRWoodsPredictions not found\n" unless (-e $printMiRWoodsPredictions);

print "Checking the miRBase Gff File:\n";
runCheckMiRBaseGff($checkMiRBaseGff,$parameters);
print "finished checking miRBaseGffFile\n";
unless (-e $parameters->{librarySizes}) {
    print "warning: library sizes file not found\n";
    print "creating " . $parameters->{librarySizes} . "\n\n";
    runGetLibrarySizes($getLibrarySizes,$parameters);
}
print "running $printReadRegions:\n";
runPrintReadRegions($printReadRegions,$parameters);
print "$printReadRegions finished\n\n";
print "running $extractProductFeatures:\n";
runExtractProductFeatures($extractProductFeatures,$parameters);
print "$extractProductFeatures finished\n\n";
if ($parameters->{trainModels}) {
    print "running $generateProductTrainData:\n";
    runGenerateProductTrainData($generateProductTrainData,$parameters);
    print "$generateProductTrainData finished\n\n"
}
print "running $evaluateProductsWithRF:\n";
runEvaluateProductsWithRF($evaluateProductsWithRF,$parameters);
print "$evaluateProductsWithRF finished\n\n";
print "running $processReadRegions\n";
runProcessReadRegions($processReadRegions,$parameters);
print "$processReadRegions finished\n\n";
if ($parameters->{mirbaseGff} || $parameters->{otherAnnotFileList}) {
    print "running $generateAnnotations:\n";
    runGenerateAnnotations($generateAnnotations,$parameters);
    print "$generateAnnotations finished\n\n";
}
print "running $extractHairpinFeatures:\n";
runExtractHairpinFeatures($extractHairpinFeatures,$parameters);
print "$extractHairpinFeatures finished\n\n";
if ($parameters->{trainModels}) {
    print "running $generateHairpinTrainData\n";
    runGenerateHairpinTrainData($generateHairpinTrainData,$parameters);
    print "$generateHairpinTrainData finished\n\n";
}
print "running $evaluateHairpinsWithRF\n";
runEvaluateHairpinsWithRF($evaluateHairpinsWithRF,$parameters);
print "$evaluateHairpinsWithRF finished\n\n";
print "printing output files\n";
runPrintMiRWoodsPredictions($printMiRWoodsPredictions,$parameters);
print "output files printed\n\n";
print "predictions are found in: predictedMiRs.gff\n";
print "product info is found in: predictedMiRs_productInfo.txt\n";
print "done.\n\n";

sub runCheckMiRBaseGff {
    my($checkMiRBaseGff,$parameters) = @_;
    my $configFile = $parameters->{LoadFromConfigFile} or die "error: LoadFromConfigFile parameter not set";
    my $miRBaseGff = $parameters->{mirbaseGff} or die "error: mirbaseGff parameter not set";
    system("$checkMiRBaseGff $configFile");
}

sub runGetLibrarySizes {
    my($getLibrarySizes,$parameters) = @_;
    my $configFile = $parameters->{LoadFromConfigFile} or die "error: LoadFromConfigFile parameter not set";
    my $bamListFile = $parameters->{bamListFile} or die "error: bamListFile parameter not set\n";
    my $librarySizesFile = $parameters->{librarySizes} or die "error: librarySizes parameter not set";
    die "$bamListFile not found\n" unless (-e $bamListFile);
    system("$getLibrarySizes $bamListFile");    
};

sub runPrintReadRegions {
    my($printReadRegions,$parameters) = @_;
    my $configFile = $parameters->{LoadFromConfigFile} or die "error: LoadFromConfigFile parameter not set";
    my $bamListFile = $parameters->{bamListFile} or die "error: bamListFile parameter not set\n";
    my $chromSizesFile = $parameters->{SizesFile} or die "error: SizesFile parameter not set\n";
    die "$configFile not Found\n" unless (-e $configFile);
    die "$bamListFile not Found\n" unless (-e $bamListFile);
    die "$chromSizesFile not Found\n" unless (-e $chromSizesFile);
    if ($parameters->{RepeatRegionsFile}) {
	my $repeatRegionsFile = $parameters->{RepeatRegionsFile};
	die "$repeatRegionsFile not Found\n" unless (-e $repeatRegionsFile);
    }
    system("$printReadRegions -L $configFile");
}

sub runExtractProductFeatures {
    my($extractProductFeatures,$parameters) = @_;
    my $configFile = $parameters->{LoadFromConfigFile} or die "error: LoadFromConfigFile parameter not set";
    my $bamListFile = $parameters->{bamListFile} or die "error: bamListFile parameter not set\n";
    my $readRegionsFile = $parameters->{readRegionsFile} or die "error: readRegionsFile parameter not set\n";
    my $chromSizesFile = $parameters->{SizesFile} or die "error: SizesFile parameter not set\n";
    my $genomeDir = $parameters->{genomeDir} or die "error: genomeDir parameter not set\n";
    die "error: $bamListFile not found\n" unless (-e $bamListFile);
    die "error: $readRegionsFile not found\n" unless (-e $readRegionsFile);
    die "error: $chromSizesFile not found\n" unless (-e $chromSizesFile);
    system("$extractProductFeatures -L $configFile");
}

sub runGenerateProductTrainData {
    my($generateProductTrainData,$parameters) = @_;
    my $configFile = $parameters->{LoadFromConfigFile} or die "error: LoadFromConfigFile parameter not set";
    my $bamListFile = $parameters->{bamListFile} or die "error: bamListFile parameter not set\n";
    my $readRegionsFile = $parameters->{readRegionsFile} or die "error: readRegionsFile parameter not set\n";
    my $chromSizesFile = $parameters->{SizesFile} or die "error: chromSizesFile parameter not set\n";
    my $genomeDir = $parameters->{genomeDir} or die "error: genomeDir parameter not set\n";
    my $mirBaseGff = $parameters->{mirbaseGff} or die "error: mirbaseGff parameter not set\n";
    my $productTrainFile = $parameters->{productTrainFile} or die "error: productTrainFile parameter not set\n";
    my $productVectorFile = $parameters->{productFeatVectorFile} or die "error: productFeatVectorFile parameter not set\n";
    die "$bamListFile not found\n" unless (-e $bamListFile);
    die "$readRegionsFile not found\n" unless (-e $readRegionsFile);
    die "$chromSizesFile not found\n" unless (-e $chromSizesFile);
    die "$mirBaseGff not found\n" unless (-e $mirBaseGff);
    system("$generateProductTrainData -L $configFile");
}

sub runEvaluateProductsWithRF {
    my($evaluateProductsWithRF,$parameters) = @_;
    my $configFile = $parameters->{LoadFromConfigFile} or die "error: LoadFromConfigFile parameter not set";
    my $productRF = $parameters->{productRF} or die "error: productRF parameter not set";
    my $productFeatVectorFile = $parameters->{productFeatVectorFile} or die "error: productFeatVectorFile parameter not set\n";
    my $predProductClassFile = $parameters->{predProductClasses} or die "error: predProductClasses parameter not set\n";
    die "$productFeatVectorFile not found\n" unless (-e $productFeatVectorFile);
    if ($parameters->{trainModels}) {
	my $productTrainFile = $parameters->{productTrainFile} or die "error: productTrainFile parameter not set";
	die "$productTrainFile not found.\n" unless (-e $productTrainFile);
    }
    system("$evaluateProductsWithRF -L $configFile");
}

sub runProcessReadRegions {
    my($processReadRegions,$parameters) = @_;
    my $configFile = $parameters->{LoadFromConfigFile} or die "error: LoadFromConfigFile parameter not set";
    my $bamListFile = $parameters->{bamListFile} or die "error: bamListFile parameter not set\n";
    my $predProductRF = $parameters->{predProductFile} or die "error: predProductRF parameter not set\n";
    my $chromSizesFile = $parameters->{SizesFile} or die "error: SizesFile parameter not set\n";
    my $genomeDir = $parameters->{genomeDir} or die "error: genomeDir parameter not set\n";
    my $librarySizesFile = $parameters->{librarySizes} or die "no librarySizes entry found in config file\n";
    die "$bamListFile not found\n" unless (-e $bamListFile);
    die "$predProductRF not found\n" unless (-e $predProductRF);
    die "$chromSizesFile not found\n" unless (-e $chromSizesFile);
    die "$librarySizesFile not found\n" unless (-e $librarySizesFile);
    system("$processReadRegions -L $configFile");    
}

sub runGenerateAnnotations {
    my($generateAnnotations,$parameters) = @_;
    my $configFile = $parameters->{LoadFromConfigFile} or die "error: LoadFromConfigFile parameter not set";
    my $hairpinsFile = $parameters->{hairpinsFile} or die "error: hairpinsFile parameter not set\n";
    my $productFile = $parameters->{productFile} or die "error: productFile parameter not set\n";
    die "$hairpinsFile not found\n" unless (-e $hairpinsFile);
    die "$productFile not found\n" unless (-e $productFile);
    system("$generateAnnotations -L $configFile");
}

sub runExtractHairpinFeatures {
    my($extractHairpinFeatures,$parameters) = @_;
    my $configFile = $parameters->{LoadFromConfigFile} or die "error: LoadFromConfigFile parameter not set";
    my $chromSizesFile = $parameters->{SizesFile} or die "error: SizesFile parameter not set\n";
    my $genomeDir = $parameters->{genomeDir} or die "error: genomeDir parameter not set\n";
    my $hairpinsFile = $parameters->{hairpinsFile} or die "error: hairpinsFile parameter not set\n";
    my $hairpinVectorFile = $parameters->{hairpinVectorFile} or die "error: hairpinVectorFile not set\n";
    die "$chromSizesFile not found\n" unless (-e $chromSizesFile);
    die "$hairpinsFile not found\n" unless (-e $hairpinsFile);
    system("$extractHairpinFeatures -L $configFile");
}

sub runGenerateHairpinTrainData {
    my($generateHairpinTrainData,$parameters) = @_;
    my $configFile = $parameters->{LoadFromConfigFile} or die "error: LoadFromConfigFile parameter not set";
    my $bamListFile = $parameters->{bamListFile} or die "error: bamListFile parameter not set\n";
    my $bamList = miRWoods::loadBamList($bamListFile) or die "error: failed to load $bamListFile\n";
    my $mirBaseGff = $parameters->{mirbaseGff} or die "error: mirBaseGFF parameter not set\n";
    my $hairpinVectorFile = $parameters->{hairpinVectorFile} or die "error: hairpinVectorFile parameter not set\n";
    my $hairpinVecTrainFile = $parameters->{hairpinTrainFile} or die "error: hairpinTrainFile parameter not set\n";
    my $productFile = $parameters->{productFile} or die "error: productFile parameter not set\n";
    my $chromSizesFile = $parameters->{SizesFile} or die "error: SizesFile parameter not set\n";
    my $genomeDir = $parameters->{genomeDir} or die "error: genomeDir parameter not set\n";
    die "$bamListFile not found" unless (-e $bamListFile);
    die "$mirBaseGff not found" unless (-e $mirBaseGff);
    die "$hairpinVectorFile not found" unless (-e $hairpinVectorFile);
    die "$productFile not found" unless (-e $productFile);
    die "$chromSizesFile not found" unless (-e $chromSizesFile);
    if ($parameters->{negGff}) {
	my $curatedNegGff = $parameters->{negGff} or die "error: negGff parameter not set\n";
	die "$curatedNegGff not found\n" unless (-e $curatedNegGff);
    }
    system("$generateHairpinTrainData -L $configFile");    
}

sub runEvaluateHairpinsWithRF {
    my($evaluateHairpinsWithRF,$parameters) = @_;
    my $configFile = $parameters->{LoadFromConfigFile} or die "error: LoadFromConfigFile parameter not set";
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
    if ($parameters->{RepeatRegionsFile}) {
	my $repeatRegionsFile = $parameters->{RepeatRegionsFile};
	my $sizesFile = $parameters->{SizesFile} or die "error: SizesFile parameter not set\n";
	die "$repeatRegionsFile not found\n" unless (-e $repeatRegionsFile);
	die "$sizesFile not found\n" unless (-e $sizesFile);
    }
    if ($parameters->{trainModels}) {
	my $hairpinTrainFile = $parameters->{hairpinTrainFile} or die "error: hairpinTrainFile parameter not set\n";
	die "$hairpinTrainFile not found.\n" unless (-e $hairpinTrainFile);
    }
    system("$evaluateHairpinsWithRF -L $configFile");
}

sub runPrintMiRWoodsPredictions {
    my($printMiRWoodsPredictions,$parameters) = @_;
    die "could not find positivePredictions.gff" unless (-e "positivePredictions.gff");
    die "could not find positivePredictions_hairpins.txt" unless (-e "positivePredictions_hairpins.txt");
    die "could not find positivePredictions_products.txt" unless (-e "positivePredictions_products.txt");
    system("$printMiRWoodsPredictions positivePredictions.gff positivePredictions_hairpins.txt positivePredictions_products.txt");
}
