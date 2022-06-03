## Introduction
Detect key Units in mosaic Tandem Repeats from representative reads at the same locus.

## Usage
uTR [-f input fasta file] [-u input representative unit string] [-h input haplotype file] [-i output table file] [-r maximum discrepancy ratio] [-o output fasta file with decomposition] [-tsd]

-f : Use a smaller input fasta file to reduce the computation time. The users are asked to use cTR (https://github.com/morisUtokyo/cTR) to cluster reads collected from individuals to unique representative reads. Afterwards, use this tool to retrieve decompositions from representative reads. The program cTR outputs the following information to the annotation of the representative read of each group, and uTR feeds the information. See the details of parameters in  https://github.com/morisUtokyo/cTR#readme.
> \>GroupSize = N, Diameter = D, RadiusFromCentroid = R, CentroidReadName = ID, CentroidReadLength = L

-u : Input a single representative unit of a putative tandem repeat in the input read, which can be computed by using mTR (https://github.com/morisUtokyo/mTR). If it is not specified, uTR automatically estimates tandem repeat units; however, mTR is better at predicting longer units than uTR (see the results in https://academic.oup.com/bioinformatics/article/37/5/612/5919583). 

-r : Give a maximum threshold on the mismatch ratio between the representaitve unit and a tandem repeat of the unit. No tandem repeat is output if the mismatch ratio exceeds this threshold. The default parameter is 0.3, which is set by:
> #define MAX_DIS_RATIO_DEFAULT 0.3 in uTR.h

-o : uTR outputs the given input fasta file annotated with decompositions identified in the reads. It also adds information, for example:

> \> (166,10,2,MTR,162),B483,read1, S=(1,AAAG,4,112,112),(2,AG,2,50,50), D=[1,AAAG,4,24],[2,AG,2,50],[1,AAAG,4,12],G,[1,AAAG,4,76],AAA, P=1111111111111111111111112222222222222222222222222222222222222222222222222211111111111101111111111111111111111111111111111111111111111111111111111111111111111111111000

when the input read is:

> AAAGAAAGAAAGAAAGAAAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAAAGAAAGAAAGGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA

(166,10,2,MTR,162),B483,read1: The first value, 166, shows the length of the centroid read "B483,read1" in the group, 10 is the number of elements in the group, 2 is the number of key units in the centroid read, MTR means that uTR detects mosaic tandem repeats (MTR) in the read, and 162 shows that of 166 bases in the read, 162 are found in mosaic tandem repeats.

S=(1,AAAG,4,112,112),(2,AG,2,50,50) : Each tuple shows a unit. In tuple (1,AAAG,4,112,112), 1 is the identifier of the unit, AAAG is the string of the unit, 4 is the unit length, 112 is the total bases in the unit occurrences (namely, 112/4=28 unit occurrences are found), and the last 112 the total bases in the tandem repeats of the unit.  

D=[1,AAAG,4,24],[2,AG,2,50],[1,AAAG,4,12],G,[1,AAAG,4,76],AAA : The decomposition of the input read by the units in S. For example, [1,AAAG,4,24] means that ID 1 (string AAAG of length 4) occurs 24/4 = 6 times, [2,AG,2,50], ID 2 (string AG of length 2) occurs 12/4 = 3 times, and so on.

P=11111111111111111111111122222222222222222222222222222222222222222222222222 : Each base is coded by each unit identifier (1,2). 

-t : Output the wall clock time to process each read in the input fasta file.　

cTR outputs a list of tandem repeat haplotypes for each locus. For examples:

> 10321 chr1 14883297 14883426 #haplotypes=54 (166,10,2,MTR,162) (152,8,2,MTR,148) (158,6,2,MTR,154) (162,6,2,MTR,158) (154,5,2,MTR,150)

The first value, 10321, is the identifier of the locus in the human chromosome 1 (chr1) that ranges from 14883297 to 14883426. #haplotypes=54 means that 54 haplotypes are found, and the first haplotype is characterized by (166,10,2,MTR,162), the meaning of which is described above.
