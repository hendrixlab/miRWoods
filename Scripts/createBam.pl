#!/usr/bin/perl -w
use strict;
$| = 1;

#this is a temporary script is used to create bam files from the starting fastq files
#it will eventually be replaced by a more elaborate script when miRWoods is finished.

my $USAGE = "USAGE:\n$0 <fastq> <adapter> <min quality> <bowtie index>\n";

my $scriptDir = 0;

my $fastqAvgQualityFilter = ($scriptDir) ? $scriptDir . "/fastqAvgQualityFilter.pl" : "fastqAvgQualityFilter.pl";
my $trimMirReads = ($scriptDir) ? $scriptDir . "/trimMirReads.pl" : "trimMirReads.pl";
my $addNHTags = ($scriptDir) ? $scriptDir . "/addNHTags.pl" : "addNHTags.pl";


#cutadapt parameters
my $caError = 0.2; #the allowable error rate for cutadapt
my $caQuality = 10;  #trim on 3' end until quality equals $caQuality
my $caMinSize = 17; #throw away reads not equal to $caMinSize

#bowtie parameters
my $seedLength = 18;
my $mismatchInSeed = 1;
my $minErrorOutsideSeed = 50;
my $maxNumMatches = 10;

my $fastq = $ARGV[0] or die $USAGE;
my $adapter = $ARGV[1] or die $USAGE;
my $avgQual = $ARGV[2] or die $USAGE;
my $bowtieIndex = $ARGV[3] or die $USAGE;

my ($fileBase) = $fastq =~ /([^\/]*)\.f.*$/;
#my $outputBase = $fileBase . "_cadpte02u3m".$caMinSize."q".$caQuality."_minAvgQual".$avgQual."_".$bowtieIndex."_bowtie_l".$seedLength."n".$mismatchInSeed."e".$minErrorOutsideSeed."m".$maxNumMatches."aS_NH_trimmed";
my $outputBase = $fileBase . "_cadpte02m".$caMinSize."q".$caQuality."_minAvgQual".$avgQual."_".$bowtieIndex."_bowtie_l".$seedLength."n".$mismatchInSeed."e".$minErrorOutsideSeed."m".$maxNumMatches."aS_NH_trimmed";
my $caFastq = $outputBase."_cadptTemp.fastq";
my $qualFiltFastq = $outputBase.".fastq";
my $caOutputFile = $outputBase."_cutadapt.out";

my $caOutput = `cutadapt -e $caError -m $caMinSize -q $caQuality -a $adapter $fastq -o $caFastq`;
#my $caOutput = `cutadapt -e $caError -u 3 -m $caMinSize -q $caQuality -a $adapter $fastq -o $caFastq`;

open(CA,">$caOutputFile") or die "failed to open $caOutputFile for writing\n";
print CA $caOutput;
close(CA);

my $out = `$fastqAvgQualityFilter $caFastq $qualFiltFastq $avgQual`;
print $out;
$out = `rm $caFastq`;
print $out;

my $samFile = $outputBase . ".sam";
my $trimmedSamFile = $outputBase . "_trimmedTemp.sam";
my $bowtieOutput = $outputBase . "_bowtieOutput.out";
$out = `bowtie -n $mismatchInSeed -e $minErrorOutsideSeed -l $seedLength -a -m $maxNumMatches -S --best --strata $bowtieIndex $qualFiltFastq $samFile >& $bowtieOutput`;
#$out = `bowtie -n $mismatchInSeed -e $minErrorOutsideSeed -l $seedLength -a -m $maxNumMatches -S $bowtieIndex $qualFiltFastq $samFile >& $bowtieOutput`;
print $out;
$out = `$trimMirReads $samFile $trimmedSamFile`;
print $out;
$out = `mv $trimmedSamFile $samFile`;
print $out;

my $NHSamFile = $outputBase . "_tempAddNHTags.sam";
my $NHTagOut = $outputBase . "_addNHTags.out\n";
$out = `$addNHTags $samFile $NHSamFile >& $NHTagOut`;
print $out;
$out = `mv $NHSamFile $samFile`;
print $out;

my $bamFile = $outputBase . ".bam";
my $sortBase = $outputBase . "_sort";
my $sortedBam = $sortBase . ".bam";
my $sortOutputFile = $outputBase . "_sortBam.out";

$out = `samtools view -bS -o /scratch/$bamFile $samFile`;
print $out;
$out = `mv /scratch/$bamFile $bamFile`;
print $out;
$out = `samtools sort $bamFile /scratch/$sortBase`;
print $out;
$out = `mv /scratch/$sortedBam $bamFile`;
print $out;
$out = `samtools index $bamFile`;
print $out;
