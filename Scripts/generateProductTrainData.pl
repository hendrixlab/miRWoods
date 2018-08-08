#!/usr/bin/perl -w
use Memory::Usage;
use miRWoods;
use Bio::DB::Sam;
use Getopt::Long;
use strict;
$| = 1;

my $USAGE = "USAGE:\n$0 -S <SizesFile> -b <BamListFile> -g <genome directory> -o <outputPrefix>
\t-S\tSizesFile\tchrom.sizes file containing the size of each chromosome.
\t-b\tbamListFile\t\tfile containing a list of sample names and bamfile locations
\t\t\t\tbamListFileExample:
\t\t\t\t\t<sample name 1>\t<bam file 1>
\t\t\t\t\t<sample name 2>\t<bam file 2>
\t\t\t\t\t<sample name 3>\t<bam file 3>
\t\t\t\t\tetc.
\t-g\tgenomeDir\tDirectory containing chromosomes in sperate fasta files
\t-o\toutputPrefix\tprefix for output files;
\t-l\tlengthMin\tminimum length for the arms of the hairpin
\t-t\ttotalLength\tmaximum length of the entire hairpin
\t-d\tdistanceMin\tdistance within which a read can still be associated with a product
\t-f\tfivePrimeHetMax\tMaximum 5` hererogeneity
\t-c\tcountMinLocus\tmaximum number of producs reads neccessary to be a real product (use total of values returned from addNHTags.pl)
\t-r\treverseMax\tmaximum allowed fraction of total antisense product reads to total (sense + antisense) product reads
\t-s\tshiftMin\tminimum shift of products
\t-h\thairpinShortLength\tlength of short arms allowed in the middle of hammer head loops
\t-O\tOverlapMin\tmin amount of product overlap needed for a same shifted and both shifted value to be recorded
\t-I\tInHairpinBuffer\tamount 5p products are allowed to cross into the loop without being considered a loop product
\t-O\tOutHairpinBuffer\tamount 3p products are allowed to cross into the loop without being considered a loop product
\t-R\tRangeOfHairpint\tRange of hairpin arms.  Outside this range are out products.
\t-M\tMaxLength\tThe maximum length for each read region reported in the read regions file\n";


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

my $bamListFile = $parameters->{bamListFile} or die "error: bamListFile parameter not set\n";
my $readRegionsFile = $parameters->{readRegionsFile} or die "error: readRegionsFile parameter not set\n";
my $chromSizesFile = $parameters->{SizesFile} or die "error: SizesFile parameter not set\n";
my $genomeDir = $parameters->{genomeDir} or die "error: genomeDir parameter not set\n";
my $mirBaseGff = $parameters->{mirbaseGff} or die "error: mirbaseGff parameter not set\n";
die "$bamListFile not found\n" unless (-e $bamListFile);
die "$readRegionsFile not found\n" unless (-e $readRegionsFile);
die "$chromSizesFile not found\n" unless (-e $chromSizesFile);
die "$mirBaseGff not found\n" unless (-e $mirBaseGff);
my $chromLengths = miRWoods::readChromLengths($chromSizesFile) or die "failed to load $chromSizesFile.\n";
my $bamList = miRWoods::loadBamList($bamListFile) or die "failed to load read data using $bamListFile.\n";
my $prodMLParameters = miRWoods::loadDefaultMLProductParameters() or die "failed to load default product phase parameters\n";
my($posClass,$negClass) = (1,-1);
my $productFeaturesFileHeader = miRWoods::getProdFeaturesFileHeader($prodMLParameters) or die "failed to generate header for product train file\n";
my $productTrainFile = $parameters->{productTrainFile} or die "error: productTrainFile parameter not set\n";
my $productVectorFile = $parameters->{productFeatVectorFile} or die "error: productFeatVectorFile parameter not set\n";
open(PTF, ">$productTrainFile") or die "failed to open $productTrainFile for writing\n";
print PTF "$productFeaturesFileHeader\n";
close(PTF);
miRWoods::printProductTrainData($bamList,$genomeDir,$chromLengths,$mirBaseGff,$posClass,$productTrainFile,$prodMLParameters,$parameters);
if ($parameters->{negGff}) {
    my $curatedNegGff = $parameters->{negGff};
    die "$curatedNegGff not found\n" unless (-e $curatedNegGff);
    miRWoods::printProductTrainData($bamList,$genomeDir,$chromLengths,$curatedNegGff,$negClass,$productTrainFile,$prodMLParameters,$parameters);
} else {
    #negative products have not been curated so we will assume that most features from the productVectorFile are not true miRs and add from that
    die "$productTrainFile not found\n" unless (-e $productTrainFile);
    die "$productVectorFile not found\n" unless (-e $productVectorFile);
    print "prodVecFile = $productVectorFile\n";
    print "prodTrainFile = $productTrainFile\n";
#    miRWoods::addRandomNegFeatures($productTrainFile,$productVectorFile);
    printNonMiRproductTrainData($productVectorFile,$mirBaseGff,-1,$productTrainFile,$prodMLParameters,$parameters);
}


sub printNonMiRproductTrainData {
    my($productVectorFile,$mirBaseGff,$class,$productTrainFile,$prodMLParameters,$parameters) = @_;
    my($annotHairpins,$annotProducts) = miRWoods::readMirbaseGff3($mirBaseGff);
    open(PTF,">>".$productTrainFile) or die "failed to open $productTrainFile for appending\n";
    open(PVF, $productVectorFile) or die "failed to open $productVectorFile for reading\n";
    my $header = <PVF>;
    while (<PVF>) {
	chomp;
	my $productLine = $_;
	my($geneId,$y,@scores) = split(/\t/,$productLine);
	unless (checkProductAnnotationOverlap($geneId,$annotHairpins,$annotProducts)) {
	    print PTF join("\t",$geneId,$class,@scores) . "\n";
	}
    }
    close(PVF);
    close(PTF);
}

sub checkProductAnnotationOverlap {
    my($geneId,$annotHairpins,$annotProducts) = @_;
    my($tag,$location) = miRWoods::parseGeneId($geneId);
    my($chrom,$start,$stop,$strand) = miRWoods::parseLocation($location);
    foreach my $hairpinInfo (@{$annotHairpins->{$chrom}}) {
	my($hpStart,$hpStop,$hpStrand,$hpId,$hpName) = @{$hairpinInfo};
	foreach my $productInfo (@{$annotProducts->{$hpId}}) {
	    my($prodChrom,$prodStart,$prodStop,$prodStrand,$prodId,$prodName) = @{$productInfo};
	    if (miRWoods::getOverlap($start,$stop,$prodStart,$prodStop)) {
		return 1;
	    }
	}
    }
    return 0;
}
