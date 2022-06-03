#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include "uTR.h"

#define MATCH_GAIN          1
#define MISMATCH_PENALTY    1
#define INDEL_PENALTY       1

//#define DEGUG_string_decomposer1
//#define DEGUG_string_decomposer2

void string_decomposer(Read *currentRead, Unit *keyUnits, int numKeyUnits, int *prio2unit){
    if(numKeyUnits < 1) return;
    // Although currentRead->intString[i] and keyUnits[b].intString[j] are 0-origin indexing, their matrix are 1-origin indexing to device wrap-around DP.
    
#ifdef DEGUG_string_decomposer1
    for(int b=0; b<numKeyUnits; b++){
        fprintf(stderr, "%d %s\n", b, keyUnits[b].string);
    }
    fprintf(stderr, "%s (len=%d)\n", currentRead->string, currentRead->len);
#endif
    
    int i, b, j, lenU;
    int lenR = currentRead->len;
    int bU = (lenR + 1); // begin of the unit
    for(int b=0; b < numKeyUnits; b++){         // Process each unit in keyUnits
        lenU = keyUnits[b].len;   // Length of the unit
        for(j=1; j <= keyUnits[b].len; j++){    // Scan the unit. 1-origin indexing.
            WrapDP[bU + lenU*0 + j] = 0;   // Each unit uses local alignment
            //WrapDP[bU + lenU*0 + j] = (-1) * INDEL_PENALTY * j; // global alignment
        }
        bU += (lenR + 1) * lenU;
    }
    
    int max_column_b[MAX_READ_LENGTH];
    
    max_column_b[0] = -1;
    int max_wrd = 0;
    int max_b = 0;
    int max_j = 0;
    int val_match, val_mismatch, val_insertion, val_deletion;
    WrapDP[0] = 0;
    for(i=1; i <= lenR; i++){   // Scan repeat. 1-origin indexing.
        int min_lenU = 10000;
        bU = (lenR + 1);
        int max_column = (-1) * INDEL_PENALTY * lenR;
        for(b=0; b < numKeyUnits; b++){     // Process each unit in keyUnits
            lenU = keyUnits[b].len;         // Length of the unit
            for(j=1; j <= lenU; j++){       // Scan the unit. 1-origin indexing.
                if( WrapDPsize <= bU + j ){
                    fprintf(stderr, "Increse WrapDPsize.\n"); exit(EXIT_FAILURE);
                }
                if(currentRead->intString[i-1] == keyUnits[b].intString[j-1]){
                    // Switch to 0-origin indexting to access intString !!
                    if(j == 1)
                        WrapDP[bU+lenU*i+j] = WrapDP[i-1] + MATCH_GAIN; // The first row
                    else
                        WrapDP[bU+lenU*i+j] = WrapDP[bU+lenU*(i-1)+(j-1)] + MATCH_GAIN;
                }else{
                    val_insertion = WrapDP[bU + lenU*(i-1) + j] - INDEL_PENALTY;
                    if(j == 1){
                        val_mismatch  = WrapDP[i-1] - MISMATCH_PENALTY; // The first row
                        // Ignore the deletion from WrapDP[i-1]
                        WrapDP[bU+lenU*i+j] = MAX(val_mismatch,val_insertion);
                    }else{
                        val_mismatch = WrapDP[bU + lenU*(i-1) + (j-1)] - MISMATCH_PENALTY;
                        val_deletion = WrapDP[bU + lenU*i     + (j-1)] - INDEL_PENALTY;
                        WrapDP[bU + lenU*i + j] = MAX( MAX( val_mismatch, val_insertion), val_deletion);
                    }
                }
                     
                if( i == lenR && max_wrd <= WrapDP[bU + lenU*i + j]  )
                    // Extend the alignment longer
                //if( i == lenR && max_wrd < WrapDP[bU + lenU*i + j] )
                {   // At the last column (i==lenR), compute the maximum.
                    max_wrd = WrapDP[bU + lenU*i + j];
                    max_b = b;  max_j = j; 
                }
            }
            // Note that j == lenU (= keyUnits[b].len)
            if( max_column <= WrapDP[bU + lenU*i + lenU] ){
            //if( max_column < WrapDP[bU + lenU*i + lenU] ){
                max_column = WrapDP[bU + lenU*i + lenU];
                max_column_b[i] = b;
            }
            bU += (lenR + 1) * lenU;
        }
        WrapDP[i] = max_column;      // wrap around. The first row.
    }
    
#ifdef DEGUG_string_decomposer1
    fprintf(stderr, "max_b=%d\tmax_j=%d\tmax_wrd=%d\n", max_b, max_j, max_wrd);
    for(i=1; i <= lenR; i++)   // Scan repeat. 1-origin indexing.
        fprintf(stderr, "%d", max_column_b[i]);
    fprintf(stderr, "\n");
#endif
    
#ifdef DEGUG_string_decomposer2
    for(i=1; i <= lenR; i++)   // Scan repeat. 1-origin indexing.
        fprintf(stderr, "%d", max_column_b[i]);
    fprintf(stderr, "\n");
    bU = 0;
    for(i=0; i <= lenR; i++)   // Scan repeat. 1-origin indexing.
        fprintf(stderr, "(%d,0,%d,%d,%d) ", i, bU, max_column_b[i], WrapDP[i] );
    fprintf(stderr, "\n\n");
    bU = (lenR + 1);
    for(b=0; b < numKeyUnits; b++){   // Process each unit in keyUnits
        lenU = keyUnits[b].len;   // Length of the unit
        for(j=1; j <= lenU; j++){     // Scan the unit. 1-origin indexing.
            for(i=0; i <= lenR; i++)   // Scan repeat. 1-origin indexing.
                fprintf(stderr, "(%d,%d,%d,%d,%d) ", i, j, bU, b, WrapDP[bU + lenU*i + j] );
            fprintf(stderr, "\n");
        }
        bU += (lenR+1) * lenU;
    }
    fprintf(stderr, "\n");
#endif

    // trace back the optimal alignment while storing it in the data structure "alignment"
    int Num_matches = 0;
    int Num_mismatches = 0;
    int Num_insertions = 0;
    int Num_deletions  = 0;
    int Num_scanned_unit = 0;
    
    i = lenR;  j = max_j;
    for(bU=(lenR+1), b=0; b < max_b; bU += (lenR+1) * keyUnits[b].len, b++){}
    // Note that b == max_b, but we make sure this.
    lenU = keyUnits[b].len;
    
    if(j == 0){ // Wrap around.
        for(bU=(lenR+1), b=0; b < max_column_b[i]; bU += (lenR+1) * keyUnits[b].len, b++){}
        // Note that b == max_column_b[i], here.
        lenU = keyUnits[b].len;
        j = lenU;
    }
    
    int blocks[MAX_READ_LENGTH+1];
    int deleted = 0;    // Print block number if the character is not deleted.
    while(i > 0){
    //while(i > 0 && WrapDP[bU + lenU*i + j] > 0){
#ifdef DEGUG_string_decomposer2
        fprintf(stderr, "(%d,%d,%d,%d,%d)->", i, j, bU, b, WrapDP[bU+lenU*i+j] );
#endif
        if(j == 1){
            val_match   = WrapDP[i-1]  + MATCH_GAIN;
            val_mismatch= WrapDP[i-1]  - MISMATCH_PENALTY;
            val_deletion= (-1) * INDEL_PENALTY * lenR;  // Avoid selecting a deletion.
        }else{
            val_match   = WrapDP[bU + lenU*(i-1) + j-1]  + MATCH_GAIN;
            val_mismatch= WrapDP[bU + lenU*(i-1) + j-1]  - MISMATCH_PENALTY;
            val_deletion= WrapDP[bU + lenU*i     + j-1]  - INDEL_PENALTY;
        }
        val_insertion   = WrapDP[bU + lenU*(i-1) + j]    - INDEL_PENALTY;
        
        if( max_wrd == val_match
            && currentRead->intString[i-1] == keyUnits[b].intString[j-1] ){
            // Switch to 0-origin indexting to access intString !!
            max_wrd -= MATCH_GAIN;
            i--; j--;
            Num_matches++;
            Num_scanned_unit++;
        }else if( max_wrd == val_mismatch
                  && currentRead->intString[i-1] != keyUnits[b].intString[j-1]){
                    // Switch to 0-origin indexting to access intString !!
            max_wrd += MISMATCH_PENALTY;
            i--; j--;
            Num_mismatches++;
            Num_scanned_unit++;
        }else if( max_wrd == val_deletion){     // deletion
            max_wrd += INDEL_PENALTY;
            j--;
            Num_deletions++;    // Num_insertions++;
            Num_scanned_unit++;
            deleted = 1;
        }else if( max_wrd == val_insertion){    // insertion
            max_wrd += INDEL_PENALTY;
            i--;
            Num_insertions++;
        }else if( max_wrd == 0){
            break;
        }else{
            fprintf(stderr, "\nfatal error in wrap-around DP at (%d,%d,%d,%d,%d)\n", i, j, bU, b, WrapDP[bU+lenU*i+j] );
            exit(EXIT_FAILURE);
        }
        if(deleted == 0){
            blocks[i] = b;  // When j is not deleted, we assign b to blocks[i] because i has been decremented and is equal to i in terms of 0-origin indexing.
        }else   deleted = 0;
        // Revise the block if necessary
        if(j == 0){
            for(bU=(lenR+1), b=0; b < max_column_b[i]; bU+=(lenR+1)*keyUnits[b].len, b++){}
            lenU = keyUnits[b].len;
            j = lenU;
        }
        // a bug ...
        //if(deleted == 0){
        //    blocks[i] = b;  // When j is not deleted, we assign b to blocks[i] because i has been decremented and is equal to i in terms of 0-origin indexing.
        //}else deleted = 0;
    }
#ifdef DEGUG_string_decomposer1
    for(int k=0; k < lenR; k++) fprintf(stderr, "%d", blocks[k]);
    fprintf(stderr, "\nmat=%d mis=%d ins=%d del=%d\n", Num_matches, Num_mismatches, Num_insertions, Num_deletions);
#endif

    //----------------------------------------------
    // Represent the decomposition by a regular expression
    //----------------------------------------------
    currentRead->discrepancy_ratio = (float) (Num_mismatches + Num_deletions + Num_insertions) / lenR;
    //currentRead->discrepancy_ratio = 1 - ((float) Num_matches / lenR);
    for(int b=0; b < numKeyUnits; b++){
        // Reset and recalculate these values.
        Units[prio2unit[b]].sumOccurrences = 0;
        Units[prio2unit[b]].sumTandem = 0;
    }

    char Decomp[MAX_READ_LENGTH];
    sprintf(Decomp, "");
    //sprintf(Decomp, "D=");
    int run = 0;
    int prevPrio = blocks[0];
    for(int i=0; i <= lenR; i++){
        if(prevPrio != blocks[i] || i == lenR){
            // Do not break a run of units in the presence of mutations
            Units[prio2unit[prevPrio]].sumOccurrences += run;
            Units[prio2unit[prevPrio]].sumTandem += run;
            // Use a pair of square brackets for tandem repeat units
            sprintf(Decomp, "%s<%s>%d", Decomp,  Units[prio2unit[prevPrio]].string, run/Units[prio2unit[prevPrio]].len);
            run = 1;
            prevPrio = blocks[i];
        }else
            run++;
    }

    // Compute a decomposition or a list of blocks
    char ListPrios[MAX_READ_LENGTH];
    // add the list of prios of all nucleotides
    sprintf(ListPrios, "");
    for(int i=0; i<lenR; i++)
        sprintf(ListPrios, "%s%d", ListPrios, blocks[i]);
    int lenDecomp=0;
    for(; Decomp[lenDecomp] !='\0'; lenDecomp++)
    if(lenDecomp < 100){ // If the description length is 100 letters or less, we assume it is readable.
        currentRead->RegExpressionDecomp = 1;
        sprintf(currentRead->RegExpression, "%s", Decomp);
    }else{
        currentRead->RegExpressionDecomp = 0;
        sprintf(currentRead->RegExpression, "%s", ListPrios);
    }
    
#ifdef DEGUG_string_decomposer1
    fprintf(stderr, "%s\n\n", currentRead->RegExpression);
#endif
}
