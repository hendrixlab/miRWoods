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


my $bamListFile = $parameters->{bamListFile};
my $readRegionsFile = $parameters->{readRegionsFile};
my $chromSizesFile = $parameters->{SizesFile};
my $genomeDir = $parameters->{genomeDir};
my $mirbaseGff = $parameters->{mirbaseGff};
my $chromLengths = miRWoods::readChromLengths($chromSizesFile);
my $bamList = miRWoods::loadBamList($bamListFile);
my $prodMLParameters = miRWoods::loadDefaultMLProductParameters();
miRWoods::generateTrainAndTestProdData($bamList,$readRegionsFile,$genomeDir,$chromLengths,$prodMLParameters,$parameters);
