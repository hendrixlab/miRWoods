#!/usr/bin/perl -w
use miRWoods;
use strict;

my $usage = "USAGE:\n$0 <miRWoods config file>";

my $configFile = $ARGV[0] or die $usage;

my($params,$paramArray) = readConfig($configFile);
my $miRBaseGff = $params->{mirbaseGff};

my($gffLines,$createNewFile,$changes) = checkMirbaseGff3($miRBaseGff);

if($createNewFile) {
    my(@filePath) = split(/\//,$miRBaseGff);
    my $fileName = pop(@filePath);
    my($gffPrefix) = $fileName =~ /^(.*).gff\d*$/;
    my $newGffFile = "miRWoodsFixed\_$gffPrefix.gff3";

    my $oldConfigFile = "$configFile.old";
    my $itr = 1;
    while (-e $oldConfigFile) {
	$itr++;
	$oldConfigFile = "$configFile.old$itr";
    } 
    my $out = `cp $configFile $oldConfigFile`;
    print $out;
    unless (-e $oldConfigFile) {
	die "failed to create $oldConfigFile\n";
    }
    printNewGffFile($newGffFile,$gffLines);
    createNewConfigFile($configFile,$oldConfigFile,$newGffFile);
    print "some issues with the gff file may cause problems with miRWoods.\n";
    print "a new config file was created with an updated gff file $newGffFile\n";
    print "the following changes were made in the updaded gff:\n";
    foreach my $change (@{$changes}) {
	print "\t$change\n";
    }
}

sub printNewGffFile {
    my($newGffFile,$gffLines) = @_;
    open(NGFF,">$newGffFile") or die "failed to open $newGffFile for writing\n";
    foreach my $line (@{$gffLines}) {
	print NGFF $line . "\n";
    }
    close(NGFF);
}

sub createNewConfigFile {
    my($configFile,$oldConfigFile,$newGffFile) = @_;
    open(NCF,">$configFile") or die "failed to open $configFile for writing\n";
    open(OCF,$oldConfigFile) or die "failed to open $oldConfigFile\n";
    while (<OCF>) {
	chomp;
	my $line = $_;
	unless ($line eq "" || ($line =~ /^\#/) ) {
	    my($param,$value) = $line =~ /^(.*)\s+=\s+(.*)\s*$/;
	    if ($param eq "mirbaseGff") {
		$line = "$param\t=\t$newGffFile";
	    }
	}
	print NCF $line . "\n";
    }
    close(OCF);
    close(NCF);
}

sub readConfig {
    my($configFile) = @_;
    my @paramArray;
    my %params;
    open(CF,$configFile) or die "failed to open $configFile\n";
    while (<CF>) {
	chomp;
	unless ($_ eq "" || ($_ =~ /^\#/) ) {
	    my($param,$value) = $_ =~ /^(.*)\s+=\s+(.*)\s*$/;
	    push(@paramArray,$param);
	    $params{$param} = $value;
	}
    }
    close(CF);
    return(\%params,\@paramArray);
}

sub checkMirbaseGff3 {
    my($mirbaseGff3) = @_;
    my @gffLines;
    my @changes;
    my $createNewFile = 0;
    my $lastHPId = "";
    open(MBGFF3,$mirbaseGff3) or die "could not open $mirbaseGff3\n";
    while(<MBGFF3>) {
	chomp;
	my $line = $_;
	unless(/^\#/) {
	    my($chrom,$source,$type,$start,$stop,$score,$strand,$phase,$info) = split(/\t/,$line);
	    my %info;
	    my @terms = split(/;/,$info);
	    foreach my $term (@terms) {
		my($key,$value) = $term =~ /(.*)=(.*)/;
		$info{$key} = $value;
	    }
	    if($type eq "miRNA_primary_transcript") {
		# hairpin region from the gff file
		my $id = $info{ID} or die "No ID found for the line:\n$_\n";
		my $name = $info{Name} or die "No Name found for the line:\n$_\n";
		push(@gffLines,$line);
		$lastHPId = $id;
	    }
	    if($type eq "miRNA") {
		# mature product from the hairpin file
		my $id = $info{ID} or die "No ID found for the line:\n$_\n";
                my $name = $info{Name} or die "No Name found for the line:\n$_\n";
		my $parentId = $info{Derives_from} or die "No Derives_from found for the line:\n$_\n";
		if ($parentId ne $lastHPId) {
		    my($testId) = split(/\_/,$lastHPId);
		    if ($parentId eq $testId) {
			push(@changes,  "changed Derives_from=$parentId to Derives_from=$lastHPId for $name");
			$parentId = $lastHPId;
			$createNewFile = 1;
			for (my $i=0; $i < @terms; $i++) {
			    if ($terms[$i] =~ /^Derives_from/) {
				$terms[$i] = "Derives_from=$parentId";
			    }
			}
			$info = join(";",@terms);
		    } else {
			print "Warning: last hairpin id in gff was $lastHPId and this product derives from $parentId.\n";
		    }
		}
		my $line = join("\t",($chrom,$source,$type,$start,$stop,$score,$strand,$phase,$info));
		push(@gffLines,$line);
	    }
	} else {
	    push(@gffLines,$line);
	}
    }
    close(MBGFF3);
    return(\@gffLines,$createNewFile,\@changes);
}
