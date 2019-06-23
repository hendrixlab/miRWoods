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
