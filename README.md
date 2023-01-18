## Introduction
The program uTR decomposes a DNA string into mosaic tandem repeats with different repeat units.

## Usage
uTR [-f input fasta file] [-o output fasta file with annotation] [-l locus information] [-u input representative unit string] [-r maximum discrepancy ratio] [-t]

-f : Feed a fasta file. For example:

    > SAND12(control,BAFME)
    AAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAATAAAATAAAATAAAATAAAATAAAATAAAATAAAAATGAACAAAA

-o : uTR outputs the input fasta file annotated with a tandem repeat pattern identified. In the running example, the annotation contains "#Pat \<AAAAT\>23", which shows a tandem repeat pattern, and "#Info (116,1,0.06)", where 116 shows the length of the DNA string, 1 means the number of unit is 1, and 0.06 is the divergence between the tandem repeat pattern and the string:

    > #Info (116,1,0.06) #Pat <AAAAT>23 #Annotation SAND12(control,BAFME)  
    AAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAATAAAATAAAATAAAATAAAATAAAATAAAATAAAAATGAACAAAA

-l : Feed the locus information (e.g., chr8:118366813-118366928), which is added to the annotation of each string in the output fasta file with tag #Locus. For example:

    #Locus chr8:118366813-118366928

-u : Feed a single representative unit of a putative tandem repeat in the input read. If it is not specified, uTR automatically estimates tandem repeat units of length 20 nt or less. Long units (of length > 20nt, for example) can be computed by using mTR: 

    https://github.com/morisUtokyo/mTR 
   
See the details of mTR:
    
    https://academic.oup.com/bioinformatics/article/37/5/612/5919583

-t : Output the wall clock time to process each read in the input fasta file.　

-r : Give a maximum threshold on the mismatch ratio between the representaitve unit and a tandem repeat of the unit. No tandem repeat is output if the mismatch ratio exceeds this threshold. The default parameter is 0.3, which is set in uTR.h by:

    >#define MAX_DIS_RATIO_DEFAULT 0.3 


## Simulated datasets to show the accuracy of the code

In the directory named "test_public", codes for generating simulated datasets and checking the accuracy of uTR and RepeatMasker are available in the following sub-directories:

- gendata: Codes for generating simulated datasets with correct answer mosaic repeat patterns.
- check_uTR: Codes for checking the accuracy of uTR according to the correct answer patterns.
- parse_RepeatMasker: Codes for parsing the output of RepeaseMasker after entering the simulated datasets (named $TR_file) into RepeatMasker by issuing the command: "repeatmasker -e hmmer -noint -pa 4 -div 0 -xsmall $TR_file" 
- check_RepeatMasker: Codes for checking the accuracy of RepeatMasker.

Use "make.sh" to compile all codes in the above subdirectories.

Use "test.sh" to perform all the steps in a batch manner. The script reports the accuracy of uTR and RepeatMasker in "accuracy_uTR.txt" and "accuracy_RepeatMasker.txt". For example, the top six rows of accuracy_uTR.txt are:

    AC AG
    2 20   1000    1000    42
    2 50   1000    1000    104
    2 100  1000    1000    209
    2 200  1000    1000    403
    2 500  998     1000    1011

The first row shows that 1000 instances mosaic tandem repeats of the form (AC)m (AG)n are generated at random. In the second row,  "2 20" implies that values of variables m and n ranged from 2 to 20, "1000 1000" means that 1000 of 1000 mosaic tandem repeat patterns were predicted correctly, the last "42" shows the average length of all patterns. Similarly, in the sixth row, "2 500" shows that m and n ranged from 2 to 500, "998 1000" means 998 of 1000 patterns were correctly estimated, and the last "1011" shows the average length. 

## Real data 

The following read data are found in the directory:

  realdata/realdata.fasta

The fasta file includes the following DNA sequences:

- SAND12(control,BAFME) hg38_dna range=chr8:118366813-118366928 in the 4th intron of SAMD12, pattern=(AAAAT)23

- SAND12(case,BAFME) Tandem repeat in the 4th intron of SAND12 found in Patient II-1 in family F6115 (Supplementary Figure 6, Ishiura, H. et al. Expansions of intronic TTTCA and TTTTA repeats in benign adult familial myoclonic epilepsy. Nat Genet 50, 581–590 (2018))  pattern=(ATTTT)221(ATTTC)221(ATTTT)82

- SAND12(case,BAFME) Tandem repeat in the 4th intron of SAND12 found in (Supplementary Figure 6, Ishiura, H. et al. Expansions of intronic TTTCA and TTTTA repeats in benign adult familial myoclonic epilepsy. Nat Genet 50, 581–590 (2018))  pattern=(ATTTT)613(ATTTC)320(ATTTT)5(ATTTC)130

- RFC1(control,CANVAS) hg38_dna range=chr4:39348425-39348483 pattern=(AAAAG)11

- KAZN(control) hg38_dna range=chr1:14883297-14883426 pattern=(AAAG)6(AG)11(AAAG)20

- ZNF37A(control) hg38_dna range=chr10:38112731-38112826 pattern=(CTTTT)12(CTTGT)3(CTTTT)2


## Handling a fasta file output by cTR

To reduce computation time, it is reasonable to generate smaller input Fasta files. For this purpose, cTR is useful because it clusters reads into groups of similar reads and is available at: 

    https://github.com/morisUtokyo/cTR 

cTR outputs the following information to the annotation of the representative read (with sampleID,readID) of each group, and uTR can feed the annotation. 

    > GroupSize = N, Diameter = D, RadiusFromCentroid = R, CentroidReadName =  sampleID,readID, CentroidReadLength = L

For example:

    > GroupSize = 10, Diameter = 4, RadiusFromCentroid = 4, CentroidReadName = sampleID,readID, CentroidReadLength = 166

uTR parses the above information and outputs the annotation:
    
    > #Info (ID,read1,166,10,2,0.01) #Pat <AAAG>6<AG>25<AAAG>23 #Decomp [2 (0,AAAG,4,116,116) (1,AG,2,50,50)] 

- #Info (ID,read1,166,10,2,0.01): ID is the identifier of a sample, read1 is the name of a read collected from the sample, and the pair of ID and read1 is the centroid of a group of strings. 166 shows the length of the centroid, 10 is the number of elements in the group, 2 is the number of key repeat units in the centroid read, and 0.01 is the mismatch ratio between the read and the decomposition \<AAAG\>6\<AG\>25\<AAAG\>23 of the above input string, which concatenates 6 copies of AAAG, 25 copies of AG, and 23 copies of AAAG.

- #Decomp [2 (0,AAAG,4,116,116) (1,AG,2,50,50)]: The first 2 shows the number of repeat units. In tuple (0,AAAG,4,116,116), 0 is the identifier of the unit, AAAG is the string of the unit, 4 is the unit length, 116 is the total bases in the unit occurrences, and the last 116 the total bases in the tandem repeats of the unit. Similarly, in tuple (1,AG,2,50,50), 1 is the identifier, AG the unit string, 2 the length of AG, 50 the total bases and the last 50 the total bases in the tandem repeats of the unit.

- At the end of the fasta file, all repeat units used in the decompositions of reads are appended. Each repeat unit is annotated with its occurrence frequency in reads.

To process the above information associated with sampleID and readID, several input parameters can be used.

uTR [-sd] [-i output table file] [-h input haplotype file]
    
-s : Do not print the pair of sample identifier ID and the name of the read. 

    > #Info (166,10,0.01) #Pat <AAAG>6<AG>25<AAAG>23 #Decomp [2 (0,AAAG,4,116,116) (1,AG,2,50,50)] 

-d : Print the decomposition only. For example, given -s and -d, we have: 
    
    > #Info (166,10,0.01) #Pat <AAAG>6<AG>25<AAAG>23 

-i : Output a table file in which each line shows, for example: 

    nID,read1 (166,10,2,0.01) [2 (0,AAAG,4,116,116) (1,AG,2,50,50)]

-h : Feed SNV information surrounding a input TR in a read; namely, a list of tuples of the form sampleID, readID, and a pair of SNV positions closest to the focal TR (e.g., 14882386|14883645, where two positions are separated by the bar "|"). For each read, the pair of nearest SNVs is put into the annotation of the read with tag #Hap.

 
