## Introduction

For comparing uTR, RepeatMasker, and TRF using synthetic datasets　in terms of accuracy of predicting mosaic tandem repeat (TR) patterns, first compile all codes by executing:

    bash make.sh

To compare the performance (in terms of accuracy and wall clock time) of uTR, RepeatMasker, and TRF, execute:

    bash test.sh

The above program creates a test data set consisting of widely different units with typical mosaic TRs present in the human genome:

    (AC)i(AG)j (ACC)i(GTT)j (AAG)i(AG)j (AAG)i(AGG)j (AAAG)i(AG)j 
    (AAAG)i(AG)j(AAAG)k (AAAG)i(AG)j(AGGG)k(AG)l(AAAG)m
    (AGGGG)i(AAAAGAAAGAGAGGG)j(AGGGG)k
    
The variable (e.g., i, j, k, l, m) next to each unit in parentheses represents the number of unit occurrences. Mosaic TRs are harder to decompose when different units are more similar, and more units are present. Mosaic TRs are harder to decompose correctly when different units are more similar, and more distinct units are present.
To understand the hardness of the decomposition, we generate a variety of datasets of different average lengths for each mosaic TR pattern; i.e., variables in each pattern are set to random values ranging from 10 to 200.
The program in the following directory generates these datasets.

    gendata/

To see how sequencing errors affect the prediction of the original mosaic TR patterns, letters of strings in each dataset are modified at random by sequencing errors (substitutions, insertions, and deletions) at the rate of 0%, 1%, 3%, 5%, 10%, and 15%.
Precisely, for example, when a mosaic TR has three units (U)i(V)j(W)k, all of the three units need to be predicted nearly correctly, and a series of units (U)i is accurate if the value of i differs by at most X% (e.g., 0%, 1%, 2%, and 3%) of the true value,  where we call X an allowance.
Accuracy increases by setting allowance X to a larger value, and this mitigation is reasonable and necessary when dealing with two homologous units (i.e., AAG and AG) because it becomes ambiguous to correctly determine the boundary between two similar units in the presence of sequencing errors.

The program compares the prediction accuracy of uTR with TRF (Version 4.09) and RepeatMasker (version open-4.0.7).
The program uses RepeatMasker with default parameter settings (-e hmmer -noint -pa 4 -div 0 -xsmall) and TRF with default parameter settings except for lowering the minimum alignment score from 50 to 10 (i.e., 2 7 7 80 10 10 1000 -h -ngs) in order to detect small TRs with 10 or more units in our benchmark datasets.
TRF sometimes returns a single most likely mosaic TR, but it often outputs a number of tandem repeats some of which overlap each other.
To find a mosaic TR, a series of non-overlapping TRs has to be selected, which is actually solved by RepeatMasker.
Therefore, RepeatMasker seems to be better suited to detect mosaic TRs than TRFs.

For each of eight mosaic TR patterns, the program considers six sequencing error rates, and creates a total of 48 ($=8 \times 6$) datasets with 1000 strings.
Programs in the following subdirectories evaluate the accuracy of uTR, RepeatMasker, and TRF in terms of the given allowance.

    check_uTR/ 
    parse_RepeatMasker/ 
    check_RepeatMasker/ 
    checkTRF/

When tool T (uTR, RM, or TRF) is used with allowance X (0%, 1%, 2%, or 3%), the accuracy table named 

    tmp/accuracy_T_allowanceX.txt 

is generated in the directory tmp.
For example, the top five lines of the file

    tmp/accuracy_uTR_allowance0.02.txt 

are:

    AC_AG_10_200_0.0 992 1000 414
    ACC_GTT_10_200_0.0 998 1000 617
    AAG_AG_10_200_0.0 991 1000 524
    AAG_AGG_10_200_0.0 998 1000 617
    AAAG_AG_10_200_0.0 982 1000 624

The first row means that of 1000 mosaic TR pattern (AC)i(AG)j, where the values of i and j are selected from 10 to 200 at random, 992 are predicted correctly by uTR when the allowance is set to 0.02, and the average length is 414. 
In the directory tmp, the Excel table named 

    tmp/accuracy_time.xlsx 

summaries the accuracy and wall clock time.
In most cases, uTR outperformed RepeatMasker and TRF in terms of prediction accuracy, and this is especially true when mosaic TRs have three or more series of units.
Prediction accuracy of uTR, RepeatMasker, and TRF tends to decrease as the sequencing error rate increases because sequencing errors obscure the original unit patterns and make prediction difficult.
