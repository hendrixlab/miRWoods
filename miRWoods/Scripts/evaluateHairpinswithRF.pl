#!/usr/bin/perl -w
use Memory::Usage;
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
	    "LoadFromConfigFile=s" => \$parameters->{LoadFromConfigFile}
    );

unless ($parameters->{LoadFromConfigFile} eq "") {
    my $configFile = $parameters->{LoadFromConfigFile};
    miRWoods::readConfigFile($configFile,$parameters);
}
miRWoods::createOutputFileParameters($parameters);

my $HPMLParameters = miRWoods::loadDefaultMLHairpinParameters();
my $sizeAdjHPMLParameters = miRWoods::loadDefaultMLHairpinParameters_sizeAdj();
my $gffFile = $parameters->{GFF};
my $productFile = $parameters->{productFile};
my $hairpinsFile = $parameters->{hairpinsFile} or die "failed to open hairpins file\n";
my $distinctReadsFile = $parameters->{distinctReadsFile} or die "failed to open distinct reads file\n";
my $hairpinTrainFile = $parameters->{hairpinTrainFile} or die "failed to open hairpin training file\n";
my $hairpinVectorFile = $parameters->{hairpinVectorFile};
my $newHairpinVectorFile = $parameters->{newHairpinVectorFile};
my $newHairpinClassFile = $parameters->{newPredHairpinClasses};
my $balancedHPTrainFile = $parameters->{hairpinTrainFileSB};
my $predHairpinClassFile = $parameters->{predHairpinClasses};


my $R = Statistics::R->new();
print "starting R\n";
$R->start();
print "R is started\n";
my $scoresToRemove = createRemovedScoresSet($sizeAdjHPMLParameters);
#my $scoresToRemove = createRemovedScoresSet($HPMLParameters);
print "scores to remove = $scoresToRemove\n";
print "balancedHPTrainFile = $balancedHPTrainFile\n";
miRWoods::evaluateFeatureVectorFile_new($R,$balancedHPTrainFile,$hairpinVectorFile,$predHairpinClassFile,$scoresToRemove,$parameters);
#miRWoods::evaluateFeatureVectorFile_new($R,$hairpinTrainFile,$hairpinVectorFile,$predHairpinClassFile,$scoresToRemove,$parameters);
print "evaluation complete\n";
createNewScoresFile($hairpinVectorFile,$predHairpinClassFile,$newHairpinVectorFile);
print "scores file created\n";

#if ($parameters->{includeHPSizes}) {
#    #when this option is picked the random forest is run once more with sizes and a small set of other parameters to utilize size 
#    #differences and result in a smaller positive list.
#    print "using hp sizes in another random forest\n";
#    $scoresToRemove = createRemovedScoresSet($HPMLParameters);
#    miRWoods::evaluateFeatureVectorFile_new($R,$hairpinTrainFile,$newHairpinVectorFile,$newHairpinClassFile,$scoresToRemove,$parameters);
#    miRWoods::printPredictedMirGff($newHairpinClassFile,$gffFile,$productFile,$hairpinsFile,$parameters);
#    $predHairpinClassFile = $newHairpinClassFile;
#}

$R->stop();
print "R is now stopped\ngoing into print predicted mirgff\n";
miRWoods::printPredictedMirGff($predHairpinClassFile,$gffFile,$productFile,$hairpinsFile,$parameters);
print "extracting gff read regions\n";
miRWoods::extractGffReadRegions("positivePredictions",$hairpinsFile,$productFile,$distinctReadsFile);
miRWoods::extractGffReadRegions("negativePredictions",$hairpinsFile,$productFile,$distinctReadsFile);
miRWoods::extractGffReadRegions("droppedOverlaps",$hairpinsFile,$productFile,$distinctReadsFile);

sub createRemovedScoresSet {
    my($MLParameters) = @_;
    my @scoresToRemove;
    foreach my $score (keys %{$MLParameters}) {
	unless ($MLParameters->{$score}) {
	    push(@scoresToRemove, $score);
	    #print "removeing $score from scores\n";
	}
    }
    return \@scoresToRemove;
}

sub createNewScoresFile {
    # this prints out subset of hairpins that are positives 
    my($hairpinVectorFile,$predHairpinClassFile,$newHairpinVectorFile) = @_;
    my %predictedHairpins;
    open(PHC,$predHairpinClassFile) or die "failed to open $predHairpinClassFile for reading\n";
    while (<PHC>) {
	chomp;
	my($geneId,$class) = split("\t",$_);
	if ($class == 1) {
	    $predictedHairpins{$geneId} = 1;
	    #print "adding $geneId to new scores file\n";
	}
    }
    close(PHC);
    open(HVF,$hairpinVectorFile) or die "failed to open $hairpinVectorFile for reading\n";
    open(NHVF,">$newHairpinVectorFile") or die "failed to open $newHairpinVectorFile for writing\n";
    my $FIRST = 1;
    while (<HVF>) {
	chomp;
	my $line = $_;
	if ($FIRST) {
	    print NHVF $line . "\n";
	    $FIRST = 0;
	} else {
	    my($geneId) = split("\t",$line);
	    if ($predictedHairpins{$geneId}) {
		print NHVF $line . "\n";
	    }
	}
    }
    close(HVF);
    close(NHVF);
}
