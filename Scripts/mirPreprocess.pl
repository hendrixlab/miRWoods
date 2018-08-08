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

miRWoods::createOutputFileParameters($parameters);
unless ($parameters->{LoadFromConfigFile} eq "") {
    my $configFile = $parameters->{LoadFromConfigFile};
    miRWoods::readConfigFile($configFile,$parameters);
}

my $bamListFile = $parameters->{bamListFile};
my $librarySizesFile = $parameters->{librarySizes} or die "no librarySizes entry found in config file\n";
my $librarySizes = miRWoods::loadLibrarySizes($librarySizesFile) or die "failed to load $librarySizesFile\n";
my $gffFile = $parameters->{GFF};
my($hairpins) = miRWoods::readMirbaseGff3($gffFile);
my $chromSizesFile = $parameters->{SizesFile};
my $genomeDir = $parameters->{genomeDir};
my $chromLengths = miRWoods::readChromLengths($chromSizesFile);
my $bamList = miRWoods::loadBamList($bamListFile);
miRWoods::mirPreprocess($bamList,$hairpins,$genomeDir,$chromLengths,$librarySizes,$parameters);
