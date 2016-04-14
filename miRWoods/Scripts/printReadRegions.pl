#!/usr/bin/perl -w
use miRWoods;
use Bio::DB::Sam;
use Getopt::Long;
use strict;
$|=1;

my $USAGE = "USAGE:\n$0 -S <SizesFile> -b <BamListFile> -o <outputPrefix> [-R <Repeat regions file>]
\t-S\tSizesFile\tchrom.sizes file containing the size of each chromosome.
\t-b\tbamListFile\t\tfile containing a list of sample names and bamfile locations
\t\t\t\tbamListFileExample:
\t\t\t\t\t<sample name 1>\t<bam file 1>
\t\t\t\t\t<sample name 2>\t<bam file 2>
\t\t\t\t\t<sample name 3>\t<bam file 3>
\t\t\t\t\tetc.
\t-R\tRepeatRegionsFile\tFile containing a list of repeat regions. (this file is optional)
\t-o\toutputPrefix\tprefix for output files;
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


my $bamListFile = $parameters->{bamListFile} or die "FAIL: Bam List file not loaded (not found in parameters).\n";
my $chromSizesFile = $parameters->{SizesFile} or die "FAIL: chrom.sizes file not loaded (not found in parameters).\n";
my $repeatRegionsFile = $parameters->{RepeatRegionsFile};
my $chromLengths = miRWoods::readChromLengths($chromSizesFile) or die "Failed to load chrom.sizes file.\n";
my $bamList = miRWoods::loadBamList($bamListFile) or die "Failed to load bamList file.\n";
my $repeatRegions = {};
if($repeatRegionsFile) {
    $repeatRegions = miRWoods::loadRepeatRegions($repeatRegionsFile, $chromLengths) or die "Failed to load repeat regions file.\n";
}
miRWoods::printReadRegions($bamList, $chromLengths, $repeatRegions, $parameters);
