#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include "uTR.h"

//#define DEBUG_wrap_around_DP
int wrap_around_DP(char *a_rep_unit, char *a_rep, int *t_mat, int *t_mis, int *t_ins, int *t_del){ // wraparound DP for handling a single unit
    
    int unit_len, rep_len;
    for(unit_len=0; a_rep_unit[unit_len] != '\0'; unit_len++){}
    for(rep_len=0;  a_rep[rep_len] != '\0'; rep_len++){}
    
    // Move 0-origin to 1-origin
    char *rep_unit = (char *) malloc(sizeof(char) * (unit_len+2));
    char *rep      = (char *) malloc(sizeof(char) * (rep_len +2));
    for(int i=0; i<unit_len; i++) rep_unit[i+1] = a_rep_unit[i];
    rep_unit[unit_len+1] = '\0';
    for(int i=0; i< rep_len; i++) rep[i+1] = a_rep[i];
    rep[rep_len+1] = '\0';
    
    int i, j;
    int next = unit_len+1;

    for(j=0; j<=rep_len; j++){   // Scan rep_unit
        WrapDP[next*0 + j] = 0; // local alignment
    }
    
    int max_wrd = 0;
    int max_i = 0;
    int max_j = 0;
    int val_match, val_mismatch, val_insertion, val_deletion;
    for(i=1; i <= rep_len; i++){        // Scan repeat
        for(j=1; j<=unit_len; j++){   // Scan rep_unit
            if( WrapDPsize <= next*i + j ){
                fprintf(stderr, "You need to increse the value of WrapDPsize.\n"); exit(EXIT_FAILURE);
            }
            if(rep[i] == rep_unit[j]){    // *1*-origin index !!!!
                WrapDP[next*i + j] = WrapDP[next*(i-1) + j-1]  + MATCH_GAIN;
            }else{
                val_mismatch    = WrapDP[next*(i-1) + j-1]  - MISMATCH_PENALTY;
                val_insertion   = WrapDP[next*(i-1) + j]    - INDEL_PENALTY;
                if(j > 1){
                    val_deletion = WrapDP[next*i + j-1] - INDEL_PENALTY;
                    WrapDP[next*i + j] = MAX(0, MAX( MAX( val_mismatch, val_insertion), val_deletion));
                }else{
                    WrapDP[next*i + j] = MAX(0, MAX( val_mismatch, val_insertion));
                }
            }
            if(max_wrd < WrapDP[next*i + j])
            {
                max_wrd = WrapDP[next*i + j];
                max_i = i;
                max_j = j;
            }
        }
        // wrap around
        WrapDP[next*i + 0] = WrapDP[next*i + unit_len];
    }
#ifdef DEBUG_wrap_around_DP
    fprintf(stderr, "%s\t%s\n", rep_unit, rep);
    fprintf(stderr, "max_wrd=%d\tmax_i=%d\tmax_j=%d\n", max_wrd, max_i, max_j);
#endif
    
    // trace back the optimal alignment while storing it in the data structure "alignment"
    int Num_matches = 0;
    int Num_mismatches = 0;
    int Num_insertions = 0;
    int Num_deletions  = 0;
    int Num_scanned_unit = 0;
    
    i = max_i;
    j = max_j;
    if(j == 0){ j = unit_len; } // 1-origin index
    int answer;
    //answer = max_j % unit_len;  // When max_j+1 be the begin of the repeat unit according to the 1-based indexing, max_j is the begin in the 0-based indexing.
    
    while(i > 0 && WrapDP[next*i + j] > 0){                 // global alignment
        val_match       = WrapDP[next*(i-1) + j-1]  + MATCH_GAIN;
        val_mismatch    = WrapDP[next*(i-1) + j-1]  - MISMATCH_PENALTY;
        val_insertion   = WrapDP[next*(i-1) + j]    - INDEL_PENALTY;
        val_deletion    = WrapDP[next*i + j-1]      - INDEL_PENALTY;
        
        if( max_wrd == val_match          && rep[i] == rep_unit[j]){
            max_wrd -= MATCH_GAIN;
            i--; j--;
            Num_matches++;
            Num_scanned_unit++;
        }else if( max_wrd == val_mismatch && rep[i] != rep_unit[j]){     // mismatch
            max_wrd += MISMATCH_PENALTY;
            i--; j--;
            Num_mismatches++;
            Num_scanned_unit++;
        }else if( max_wrd == val_deletion){     // deletion
            max_wrd += INDEL_PENALTY;
            j--;
            Num_deletions++;    // Num_insertions++;
            Num_scanned_unit++;
        }else if( max_wrd == val_insertion){    // insertion
            max_wrd += INDEL_PENALTY;
            i--;
            Num_insertions++;
            //Num_scanned_unit++;       // The base of the repeat unit is skipped.
        }else if( max_wrd == 0){
            break;
        }else{
            fprintf(stderr, "fatal error in wrap-around DP max_wrd = %i\n", max_wrd);
            exit(EXIT_FAILURE);
        }
        if(j == 0){
            j = unit_len;
        }
    }
    *t_mat = Num_matches;
    *t_mis = Num_mismatches;
    *t_ins = Num_insertions;
    *t_del = Num_deletions;
    free(rep_unit);
    free(rep);
    answer = j % unit_len;  // When j+1 be the begin of the repeat unit according to the 1-based indexing, j%unit_len is the begin in the 0-based indexing.
    return(answer);
}
