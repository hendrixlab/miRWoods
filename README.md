# miRWoods

miRWoods is a software for the
    sensitive detection microRNAs, including those with only one read.
    miRWoods uses two random forests. The 
    first random forest referred to as the miR Product Random Forsest (MPRF) 
    assembles read stacks into products and scores them. Products which 
    don't pass through the MRPF are filtered out, and those which do are 
    combined with surrounding products and folded to produce hairpins. 
    Each hairpin is scored and then passed through a Hairpin Precursor 
    Random Forest (HPRF) to obtain the final results. miRWoods
    evaluates several possible overlapping precursors for each loci, and pick
    the one with the best score.
    
# Dependencies

miRWoods Requires the Vienna RNA fold package and Perl Library which may be downloaded here:
https://www.tbi.univie.ac.at/RNA/

Cutadapt is used for adapter trimming and may be downloaded using the following command:
`pip3 install --user --upgrade cutadapt`

Bowtie is used for read mapping and may be downloaded here:
http://bowtie-bio.sourceforge.net/index.shtml

Samtools is needed to handle sam and bam files and may be downloaded here:
http://www.htslib.org/

The following libraries are also needed and can be downloaded from CPAN:
- Bio::DB::Sam
- Statistics::R
- Memory::Usage

miRWoods also uses some R code which may be downloaded here:
https://www.r-project.org/  

R will require the following libraries which may be downloaded from CRAN:
- ROCR
- randomForest
- mlbench

# Installation

1. Add the Scripts directory to your path
2. Add the Perl Modules directory to your PERL5LIB path

# Preparing Bam Files To Be Used By miRWoods

The software takes it's inputs in the form of bam files.  These may be created from fastq files using the following command:
`createBam.pl < fastq file > < 3' adaptor > < min avg quality score > < bowtie index >`

example: `createBam.pl SRR326279.fastq ATCTCGTATGCCGTCTTCTGCTTGT 30 hg19`

Here < min avg quality score > is set to the average quality score requred for a read to be kept.  We set < min avg quality score > to 30 for our experiments but depending on the fastq file you may need to set this differently.

After creating each of the bam files you should place them into a tab delimited bam file list (see example_bamlist.txt for an example.)

# Running miRWoods

The easiest way to pass inputs to miRWoods is through a config file (see example_config.txt)  The config file includes the following information:

- scriptDir - the location of the miRWoods/Scripts directory
- hairpinRF - the location of the HPRF model file
- productRF - the location of the MPRF model file
- SizesFile - A tab chrom.sizes file containing each chrom and it's size 
- genomeDir - the directory of the genome (each chrom must be split into individual files)
- geneModels - a gff file containing the genome annotations
- mirbaseGff - the gff with the mirbase annotations 
- bamListFile - a tab delimited file containing the list of bam files
- minReadLength - the minimum read length for a mapped read to be processed
- maxReadLength - the maximum read length for a mapped read to be processed
- HPDVCutoff - the HPDVCutoff threshold
- ARPMCutoff - the ARPMCutoff threshold

After the config file is setup, miRWoods may be run with the following command:
`miRWoods -L <config file>`

# Outputs

The following two files are output by miRWoods:
- predictedMiRs.gff - a gff file containing the hairpins and products predicted by miRWoods
- predictedMiRs_productInfo.txt - a tab delimited file containing read count information for each product. 
