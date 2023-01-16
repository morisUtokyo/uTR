## Introduction

For comparing uTR, RepeatMasker, and TRF using synthetic datasets, first compile all codes by executing:

> bash make.sh

To compare the performance (in terms of accuracy and wall clock time) of uTR, RepeatMasker, and TRF, execute:

> bash test.sh

The above program creates a test data set consisting of widely different units with typical mosaic TRs present in the human genome:

    (AC)i(AG)j (ACC)i(GTT)j (AAG)i(AG)j (AAG)i(AGG)j (AAAG)i(AG)j 
    (AAAG)i(AG)j(AAAG)k (AAAG)i(AG)j(AGGG)k(AG)l(AAAG)m
    (AGGGG)i(AAAAGAAAGAGAGGG)j(AGGGG)k
    
The variable (e.g., i, j, k, l, m) next to each unit in parentheses represents the number of unit occurrences. Mosaic TRs are harder to decompose when different units are more similar, and more units are present. To understand the hardness of the decomposition, we generate a variety of datasets of different average lengths for each mosaic TR pattern; i.e., variables in each pattern were set to 1000 combinations of random values ranging from 10 to 200.
To see how sequencing errors affect the prediction of the original mosaic TR patterns, strings in each dataset were replaced with random sequencing errors (substitutions, insertions, and deletions) at the rate of 0%, 1%, 3%, 5%, 10%, and 15%.

Precisely, for example, when a mosaic TR has three units (U)i(V)j(W)k, all of the three units need to be predicted nearly correctly, and a series of units (U)i is accurate if the value of i differs by at most X% (e.g., 0%, 1%, 2%, and 3%) of the true value, which we call an allowance. 
The mitigation condition (e.g., 1%, 2%, and 3%) is reasonable and necessary because when to handle two homologous units (i.e., AAG and AG), correctly determining the boundary between two similar units becomes ambiguous (especially in the presence of sequencing errors).

The program compares the prediction accuracy of uTR with TRF (Version 4.09) and RepeatMasker (version open-4.0.7) in the presence of sequencing errors.
TRF sometimes returns a single most likely mosaic TR, but it typically outputs a number of tandem repeats, some of which overlap each other and demands us to select a series of non-overlapping TRs, which is a nontrivial task to do.
RepeatMasker partly solves this problem because RepeatMasker uses TRF to generate overlapping TRs first, selects a series of non-overlapping TRs, and treats it as a mosaic TR.
For each of eight mosaic TR patterns, the program considers six sequencing error rates, and creates a total of 48 ($=8 \times 6$) datasets with 1000 strings.
The program uses RepeatMasker with default parameter settings (-e hmmer -noint -pa 4 -div 0 -xsmall) and TRF with default parameter settings except for lowering the minimum alignment score from 50 to 10 (i.e., 2 7 7 80 10 10 1000 -h -ngs) in order to detect small TRs with 10 or more units in our benchmark datasets.

To define the accuracy of prediction by each method, the predicted number of unit occurrences is allowed to differ by at most X% of the true value, where allowance X is set to 0%, 1%, 2%, or 3%.
When tool T (uTR, RM, or TRF) is used with allowance X (0%, 1%, 2%, or 3%), the accuracy table named accuracy_T_allowanceX.txt is generated in the directory tmp.
For example, the top five lines of the accuracy_uTR_allowance0.02.txt are:

    AC_AG_10_200_0.0 992 1000 414
    ACC_GTT_10_200_0.0 998 1000 617
    AAG_AG_10_200_0.0 991 1000 524
    AAG_AGG_10_200_0.0 998 1000 617
    AAAG_AG_10_200_0.0 982 1000 624

The first row means that of 1000 mosaic TR pattern (AC)i(AG)j, where the values of i and j are selected from 10 to 200 at random, 992 are predicted correctly by uTR when the allowance is set to 0.02, and the average length is 414. 
In the directory tmp, the Excel table named 

> accuracy_time.xlsx 

summaries the accuracy and wall clock time.
uTR outperformed RepeatMasker and TRF in terms of prediction accuracy, and this is especially true when mosaic TRs have three or more series of units.
Prediction accuracy of uTR, RepeatMasker, and TRF tends to decrease as the sequencing error rate increases because sequencing errors obscure the original unit patterns and make prediction difficult.


