#!/usr/bin/perl -w
use miRWoods;
use strict;

my $usage = "USAGE:\n$0 <predictions gff> <Hairpins File> <Product File>\n";

my $predGff = $ARGV[0] or die $usage;
my $hairpinsFile = $ARGV[1] or die $usage;
my $productsFile = $ARGV[2] or die $usage;

my $outputGff = "predictedMiRs.gff";
my $outputProductFile = "predictedMiRs_productInfo.txt";

my($gffEntries,$gffLines,$aliases) = readGffLines($predGff);
my($foldLocations,$hpLocations) = retrieveGffFoldLocations($gffEntries,$hairpinsFile);
my $productLocations = getProductGffInfo($hpLocations,$productsFile);
printGff($gffEntries,$foldLocations,$productLocations,$aliases,$outputGff);
printProductFile($productsFile,$hpLocations,$outputProductFile);

sub printProductFile {
    my($productsFile,$hpLocations,$outputProductFile) = @_;
    open(OPF,">$outputProductFile") or die "failed to open $outputProductFile for writing\n";
    open(PF,$productsFile) or die "failed to open $productsFile\n";
    my $FIRST = 1;
    while(<PF>) {
	chomp;
	unless ($FIRST) {
	    my($tag,$side,$type,$totalReads,$totalMostAbundant,$adjustedReads,$adjTotalMostAbundant,$start,$stop,$strand,$sequence,@data) = split(/\t/);
	    my $prodLocation = getGenomicProductLocation($hpLocations->{$tag},$start,$stop,$strand);
	    $adjustedReads = sprintf("%.3f", $adjustedReads);
	    print OPF "$tag-$side-$type\t$prodLocation\t$totalReads\t$totalMostAbundant\t$adjustedReads\t$sequence\t".join("\t",@data)."\n";
	} else {
	    #tag    side type       type    total reads     total most abundant     adjusted reads  adj total most abundant start   stop    strand  sequence        SRR326279
	    my($tag,$side,$type,$totalReads,$totalMostAbundant,$adjustedReads,$adjTotalMostAbundant,$start,$stop,$strand,$sequence,@headers) = split(/\t/);
	    print OPF "#product\tlocation\t$totalReads\t$totalMostAbundant\t$adjustedReads\tsequence\t".join("\t",@headers)."\n";
	    $FIRST = 0;
	}
    }
    close(PF);
    close(OPF);
}

sub printGff {
    my($gffEntries,$foldLocations,$productLocations,$aliases,$outputGff) = @_;
    open(OGFF,">$outputGff") or die "failed to open $outputGff for writing\n";
    foreach my $gChrom (sort {$a cmp $b} keys %{$gffEntries}) {
	foreach my $tag (keys %{$gffEntries->{$gChrom}}) {
	    my $location = $foldLocations->{$tag};
	    my($chrom,$start,$stop,$strand) = miRWoods::parseLocation($location);
	    my $hpInfo = "ID=$tag;Alias=$aliases->{$chrom}{$tag};Name=$tag";
	    print OGFF "$chrom\t.\tmiRNA_primary_transcript\t$start\t$stop\t.\t$strand\t.\t$hpInfo\n";
	    foreach my $productInfo (@{$productLocations->{$tag}}) {
		my($prodChrom,$prodStart,$prodStop,$prodStrand,$side,$type,$sequence) = @{$productInfo};
		my $gffProdInfo = "ID=$tag\-$side\-$type;Name=$tag\-$side\-$type;Derives_from=$tag";
		print OGFF "$prodChrom\t.\tmiRNA\t$prodStart\t$prodStop\t.\t$prodStrand\t.\t$gffProdInfo\n";
	    }
	}    
    }
    close(OGFF);
}

sub getProductGffInfo {
    my($hpLocations,$productsFile) = @_;
    my %productLocations;
    open(PF,$productsFile) or die "failed to open $productsFile\n";
    while (<PF>) {
	chomp;
	unless ( /^\#/ ) {
	    my($rr,$side,$type,$total,$totalMostAbundant,$adjTotal,$adjTotalMostAbundant,$relStart,$relStop,$strand,$sequence) = split(/\t/);
	    if ($type eq "miR" ||
		$type eq "long-miR" ||
		$type eq "short-miR") {
		my $prodLocation = getGenomicProductLocation($hpLocations->{$rr},$relStart,$relStop,$strand);
		my($prodChrom,$prodStart,$prodStop,$prodStrand) = miRWoods::parseLocation($prodLocation);
		push(@{$productLocations{$rr}},[$prodChrom,$prodStart,$prodStop,$prodStrand,$side,$type,$sequence]);
	    }
	} 
    }
    close(PF);
    return \%productLocations;
}

sub getGenomicProductLocation {
    my($hairpinLocation,$relStart,$relStop,$productStrand) = @_;
    my($chrom,$start,$stop,$strand) = miRWoods::parseLocation($hairpinLocation);
    my $prodStart = ($strand eq '+') ? $start + $relStart : $stop - $relStop;
    my $prodStop = ($strand eq '+') ? $start + $relStop : $stop - $relStart;
    my $productLocation = "$chrom:$prodStart..$prodStop:$strand";
    return $productLocation;
}

sub retrieveGffFoldLocations {
    my($gffEntries,$hairpinsFile) = @_;
    my %foldLocations;
    my %hpLocations;
    open(HP,$hairpinsFile) or die "failed to open $hairpinsFile\n";
    while (<HP>) {
	chomp;
	unless ( /^#/ ) {
	    my($tag,$chrom,$start,$stop,$strand,$lc,$rc,$totalSense,$totalAntisense,$mfe,$sequence,$fold) = split(/\t/);
	    if ($gffEntries->{$chrom}{$tag}) {
		my $location = "$chrom:$start..$stop:$strand";
		my $foldLocation = miRWoods::getFoldGenomicLocation($location,$fold);
		$foldLocations{$tag} = $foldLocation;
		$hpLocations{$tag} = $location;
	    }
	}
    }
    close(HP);
    return(\%foldLocations,\%hpLocations);
}

sub readGffLines {
    my($predGff) = @_;
    my @gffLines;
    my %gffEntries;
    my %aliases;
    open(PGFF,$predGff) or die "failed to open $predGff";
    while (<PGFF>) {
	chomp;
	unless ( /^#/ ) {
	    my $gffLine = $_; 
	    push(@gffLines,$gffLine);
	    my($chrom,$source,$type,$start,$stop,$score,$strand,$phase,$info) = split(/\t/);
	    my %info;
	    my @terms = split(/;/,$info);
	    foreach my $term (@terms) {
		my($key,$value) = $term =~ /(.*)=(.*)/;
		$info{$key} = $value;
	    }
	    $gffEntries{$chrom}{$info{ID}} = 1;
	    $aliases{$chrom}{$info{ID}} = $info{Alias};
	}
    }
    return(\%gffEntries,\@gffLines,\%aliases);
}
