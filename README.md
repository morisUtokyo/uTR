## Introduction
Decompose a DNA string into mosaic tandem repeats with different repeat units

## Usage
uTR [-f input fasta file] [-u input representative unit string] [-h input haplotype file] [-i output table file] [-r maximum discrepancy ratio] [-o output fasta file with decomposition] [-tsd]

-f : Feed a fasta file. Input a smaller input fasta file to reduce the computation time. For this purpose, cTR is useful because it clusters reads collected from individuals to unique TR representatives (https://github.com/morisUtokyo/cTR). cTR outputs the following information to the annotation of the representative read of each group, and uTR feeds the annotation. 

> \>GroupSize = N, Diameter = D, RadiusFromCentroid = R, CentroidReadName = ID, CentroidReadLength = L

For example:

> \> GroupSize = 10, Diameter = 4, RadiusFromCentroid = 4, CentroidReadName = sampleID,readID, CentroidReadLength = 166

-u : Feed a single representative unit of a putative tandem repeat in the input read, which can be computed by using mTR (https://github.com/morisUtokyo/mTR). If it is not specified, uTR automatically estimates tandem repeat units; however, mTR is better at predicting longer units than uTR (see the results in https://academic.oup.com/bioinformatics/article/37/5/612/5919583). 

-h : Feed SNV information surrounding a focal TR in a read; namely, a list of tuples of the form sampleID, readID, and a pair of SNV positions closest to the focal TR (e.g., 14882386|14883645, where two positions are separated by the bar "|"). For each read, the pair of nearest SNVs is put into the annotation of the read.

-r : Give a maximum threshold on the mismatch ratio between the representaitve unit and a tandem repeat of the unit. No tandem repeat is output if the mismatch ratio exceeds this threshold. The default parameter is 0.3, which is set by:
> #define MAX_DIS_RATIO_DEFAULT 0.3 in uTR.h

-o : uTR outputs the given input fasta file annotated with decompositions identified in the reads. For example:

> \> (ID,read1,166,10,2,0.01) [2 (0,AAAG,4,116,116) (1,AG,2,50,50)] \<AAAG\>6\<AG\>25\<AAAG\>23
> AAAGAAAGAAAGAAAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAAAGAAAGAAAGGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA

- (ID,read1,166,10,2,0.01): ID is the identifier of a sample, read1 is the name of a read collected from the sample, and the pair of ID and read1 is the centroid of a group of strings. 166 shows the length of the centroid, 10 is the number of elements in the group, 2 is the number of key repeat units in the centroid read, and 0.01 is the mismatch ratio between the read and the decomposition \<AAAG\>6\<AG\>25\<AAAG\>23 of the above input string, which concatenates 6 copies of AAAG, 25 copies of AG, and 23 copies of AAAG.
- [2 (0,AAAG,4,116,116) (1,AG,2,50,50)]: The first 2 shows the number of repeat units. In tuple (0,AAAG,4,116,116), 0 is the identifier of the unit, AAAG is the string of the unit, 4 is the unit length, 116 is the total bases in the unit occurrences, and the last 116 the total bases in the tandem repeats of the unit. Similarly, in tuple (1,AG,2,50,50), 1 is the identifier, AG the unit string, 2 the length of AG, 50 the total bases and the last 50 the total bases in the tandem repeats of the unit.
- \<AAAG\>6\<AG\>25\<AAAG\>23 : Decomposition that concatenates 6 copies of AAAG, 25 copies of AG, and 23 copies of AAAG.
- At the end of the fasta file, all repeat units used in the decompositions of reads are appended. Each repeat unit is annotated with its occurrence frequency in reads.

-i : Output a table file in which each line shows, for example: nID,read1 (166,10,2,0.01) [2 (0,AAAG,4,116,116) (1,AG,2,50,50)]

-t : Output the wall clock time to process each read in the input fasta file.　

-s : Do not print the pair of sample identifier ID and the name of the read. 

> \> (166,10,0.01) [2 (0,AAAG,4,116,116) (1,AG,2,50,50)] \<AAAG\>6\<AG\>25\<AAAG\>23

-d : Print the decomposition only.

> Given -s and -d, 
> 
> \> (166,10,0.01) \<AAAG\>6\<AG\>25\<AAAG\>23
