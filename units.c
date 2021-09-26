/*
 Copyright (c) 2021, Shinichi Morishita
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 The views and conclusions contained in the software and documentation are those
 of the authors and should not be interpreted as representing official policies,
 either expressed or implied, of the FreeBSD Project.
 */

#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include "MT.h"     // Use the Mersenne Twister.
#include "uTR.h"
//#include "lzp-sa-lcp_sm.h"



//  Allow mismatches for units of length > TH_LONG_UNITS, (int) 1/MAX_DIS_RATIO

//#define DEBUG_match_bounded_DP
int match_bounded_DP(char *s0, int n0, char *s1, int n1, int **mat){
    // Compute the maximum limit of disagreements
    int max_dis =  ceil(MAX_DIS_RATIO * n0);
        #ifdef DEBUG_match_bounded_DP
        fprintf(stderr, "max_dis=%d\n", max_dis);
        #endif
    int window = MAX(1, ceil(n0 * MAX_DIS_RATIO));
        #ifdef DEBUG_match_bounded_DP
        fprintf(stderr, "window=%d\n", window);
        #endif
    // Initialize the score matrix
    int window_plus = window+1;
    for(int i=1; i < n0+1; i++)
        for(int j = MAX(1, i-window_plus); j < MIN(n1+1, i+window_plus); j++)
            mat[i][j] = -(n0+n1);
    
    // Base cases
    for(int i=0; i < n0+1; i++) mat[i][0] = i * INDEL;
    for(int j=0; j < n1+1; j++) mat[0][j] = j * INDEL;
    // Inductive steps
    int maxScore = 0;
    for(int i=1; i < n0+1; i++){
        for(int j=MAX(1,i-window); j < MIN(n1+1,i+window); j++){ // Bounded
            int oneMatch;
            if( s0[i-1] == s1[j-1] ) oneMatch = MATCH; else oneMatch = MISMATCH;
            mat[i][j] = MAX(  mat[i-1][j-1] + oneMatch,
                        MAX(  mat[i-1][j] + INDEL,
                              mat[i][j-1] + INDEL ));
            maxScore = MAX(maxScore, mat[i][j]);
        }
        if( 2*max_dis < (i - maxScore) ){
            #ifdef DEBUG_match_bounded_DP
            fprintf(stderr, "maxScore=%d\n", maxScore);
            #endif
            return(0);
        }
    }
        #ifdef DEBUG_match_bounded_DP
        fprintf(stderr, "score=%d\n", mat[n0][n1]);
        #endif
    return(1);
}

void match_bounded_DP_traceback(char *s0, int n0, char *s1, int n1, int *covered){
    // Generate an (n0+1) x (n1+1) matrix
    int **mat;
    mat = malloc( sizeof(int *) * (n0+1) );
    for(int i=0; i < n0+1; i++)
        mat[i] = malloc( sizeof(int) * (n1+1) );
    
    // Compute the maximum limit of disagreements
    int max_dis =  ceil(MAX_DIS_RATIO * n0);
    int window = MAX(1, ceil(n0 * MAX_DIS_RATIO) );
    // Initialize the score matrix
    int window_plus = window+1;
    for(int i=1; i < n0+1; i++)
        for(int j = MAX(1, i-window_plus); j < MIN(n1+1, i+window_plus); j++)
            mat[i][j] = -(n0+n1);
    
    // Base cases
    for(int i=0; i < n0+1; i++)
        for(int j=0; j < n1+1; j++) mat[i][j] = (i+j) * INDEL;
    // Inductive steps
    int maxScore = 0; int i_max = 0; int j_max = 0;
    for(int i=1; i < n0+1; i++){
        for(int j=MAX(1,i-window); j < MIN(n1+1,i+window); j++){ // Bounded
            int oneMatch;
            if( s0[i-1] == s1[j-1] ) oneMatch = MATCH;
            else                     oneMatch = MISMATCH;
            mat[i][j] = MAX(  mat[i-1][j-1] + oneMatch,
                        MAX(  mat[i-1][j] + INDEL,
                              mat[i][j-1] + INDEL ));
            if(maxScore < mat[i][j]){
                maxScore = mat[i][j]; i_max = i; j_max = j;
            }
        }
    }
    // Trace back the above DP
    int i = i_max;
    int j = j_max;
    for(;;){
        if(0 < i && 0 < j && maxScore == mat[i-1][j-1] + MATCH){
            covered[i-1] = 1;
            i--; j--; maxScore -= MATCH;             }else
        if(0 < i && 0 < j && maxScore == mat[i-1][j-1] + MISMATCH){
            i--; j--; maxScore -= MISMATCH;          }else
        if(0 < i && 0 <= j && maxScore == mat[i-1][j]   + INDEL){
            i--;      maxScore -= INDEL;             }else
        if(0 <= i && 0 < j && maxScore == mat[i][j-1]   + INDEL){
            j--;      maxScore -= INDEL;
        }else
            break;
    }
    for(int i=0; i<(n0+1); i++) free(mat[i]);
    free(mat);
}


void rotate(int start, char *s1, int n1, char *rotated_s1){
    for(int i=0; i<n1; i++){
        rotated_s1[i] = s1[ (start + i ) % n1 ];
    }
    rotated_s1[n1] = '\0';
}

int rotate_match(char *S0, int n0, char *S1, int n1){
    // Make sure that n0 >= n1
    char *s0, *s1;
    if(n0 < n1){
        char *tmp_s = S1; s1 = S0; s0 = tmp_s;
        int   tmp_n = n1; n1 = n0; n0 = tmp_n;
    }else{
        s0 = S0;  s1 = S1;
    }
    // Generate an (n0+1) x (n1+1) matrix
    int **mat;
    mat = malloc( sizeof(int *) * (n0+1) );
    for(int i=0; i < n0+1; i++)
        mat[i] = malloc( sizeof(int) * (n1+1) );
    //
    char *rotated_s1 = malloc(sizeof(char) * (n1+1));
    
    int match = 0;
    for(int start=0; start<n1; start++){
        rotate(start, s1, n1, rotated_s1);
        match = match_bounded_DP(s0, n0, rotated_s1, n1, mat);
        if(match == 1) break;
    }
    for(int i=0; i<(n1+1); i++) free(mat[i]);
    free(mat);
    free(rotated_s1);
    return(match);
}

int get_threshold_top_k(int *a, int len, int topK){
    // descending order by the insertion sort
    for(int i=0; i<len; i++){
        for(int j=i; 0 <= j-1; j--){
            if(a[j-1] < a[j]){
                int tmp = a[j]; a[j]=a[j-1]; a[j-1]=tmp;
            }
        }
    }
    // Many could be given the same rank.
    int b = 0;
    for(int i=0; i < MIN(len, topK)-1; i++){
        if(a[i] > a[i+1])   b = i+1;
    }
    return(a[b]);
}

void retain_top_k_units(int topK){
    if(unit_cnt <= topK)  return;  // Do not reduce units.
    
    int *array_sumOccurrences = malloc(sizeof(int) * unit_cnt);
    for(int i=0; i<unit_cnt; i++)
        array_sumOccurrences[i] = Units[i].sumOccurrences;
    int threshold_top_k =
        get_threshold_top_k(array_sumOccurrences, unit_cnt, topK);
    free(array_sumOccurrences);
    
    int j=0;
    for(int i=0; i<unit_cnt; i++){
        if(threshold_top_k == 0) threshold_top_k++;
        if(threshold_top_k <= Units[i].sumOccurrences){
            // It is safe to overwrite Units[j] as j <= i.
            Units[j].ID  = Units[i].ID;
            Units[j].len = Units[i].len;
            Units[j].sumOccurrences = Units[i].sumOccurrences;
            for(int l=0; l < Units[i].len; l++)
                Units[j].string[l] = Units[i].string[l];
            Units[j].string[Units[i].len] = '\0';
            j++;
        }
    }
    int prev_unit_cnt = unit_cnt;
    unit_cnt = j;
    // Clear the remining elements
    for(int j=unit_cnt; j<prev_unit_cnt; j++){
        Units[j].ID  = 0;
        Units[j].len = 0;
        Units[j].string[0] = '\0';
        Units[j].sumOccurrences = 0;
    }
}

int char2quadratic(char c){
    int ans;
    switch(c){
        case 'A': ans = 0; break;
        case 'C': ans = 1; break;
        case 'G': ans = 2; break;
        case 'T': ans = 3; break;
        default: fprintf(stderr, "Invalid char: %c in char2quadratic\n", c); exit(EXIT_FAILURE);
    }
    return(ans);
}

int min_quadratic_ID(char *s, int len){
    int n, min_n;
    n = 0;
    if(len <= 15){ // 4^15 = 2^30
        for(int i=0; i<len; i++)
            n = 4*n + char2quadratic(s[i]);
    }else{
        int p = 536870909; // a prime number "p" s.t. 4p+3 < 2^31-1
        for(int i=0; i<len; i++)
            n =  (4*n + char2quadratic(s[i])) % p;
    }
    return(n);
}

int min_quadratic_ID_with_rotation(char *s, int len){
    int n, min_n;
    if(len <= 15){ // 4^15 = 2^30
        n = 0;
        for(int i=0; i<len; i++)
            n = 4*n + char2quadratic(s[i]);
        // Rotate s
        int pow4 = pow(4,len-1);
        min_n = n;
        for(int i=0; i<len; i++){
            n = 4*(n % pow4) + char2quadratic(s[i]);
            min_n = MIN( min_n, n );
        }
        return(min_n);
    }else{
        int p = 536870909; // a prime number "p" s.t. 4p+3 < 2^31-1
        for(int i=0; i<len; i++){
            n = 0;
            for(int j=0; j<len; j++)
                n =  (4*n + char2quadratic(s[(i+j)%len])) % p; // rotate s
            if(i==0) min_n = n; else min_n = MIN(min_n, n);
        }
        return(min_n);
    }
}

void put_repUnit(char *tmpUnit){
    int len;
    for(len=0; tmpUnit[len] != '\0'; len++);
    if(len == 0) return;
    int tmpID = min_quadratic_ID_with_rotation(tmpUnit, len);
    
    // RepUnit is the first unit and is put into the database.
    Units[0].ID  = tmpID;
    Units[0].len = len;
    free(Units[0].string);
    Units[0].string = malloc( sizeof(char) * (len+1) );
    for(int i=0; i<len; i++)
        Units[0].string[i] = tmpUnit[i];
    Units[0].string[len] = '\0';
    unit_cnt=1;
    
    GlobalUnits[0].ID  = tmpID;
    GlobalUnits[0].len = len;
    free(GlobalUnits[0].string);
    GlobalUnits[0].string = malloc( sizeof(char) * (len+1) );
    for(int i=0; i<len; i++)
        GlobalUnits[0].string[i] = tmpUnit[i];
    GlobalUnits[0].string[len] = '\0';
    GlobalUnits[0].sumOccurrences = 0;
    global_unit_cnt = 1;
}

void put_unit(char *tmpUnit){
    // Assume that tmpUnit is non-self-overlapping
    int len;
    for(len=0; tmpUnit[len] != '\0'; len++);
    if(len == 0) return;
    if( MAX_UNIT_LENGTH < len ) return; // Discard if the length exceeds the max.
    // Check if tmpUnit is in the database.
    int tmpID = min_quadratic_ID_with_rotation(tmpUnit, len);
    
    if(len < LONG_UNIT_LEN_TH){
        // Use exact match
        for(int i = 0; i < unit_cnt; i++)
            if(Units[i].len == len && Units[i].ID == tmpID)    return;
    }else{
        // Allow some mismatches for long units
        for(int i = 0; i < unit_cnt; i++){
            int diff = MAX(len, Units[i].len) - MIN(len, Units[i].len);
            if( diff <= (int) (MAX(len, Units[i].len) * MAX_DIS_RATIO) ){
                int match = rotate_match( tmpUnit, len, Units[i].string, Units[i].len);
                if(match == 1)  return;
            }
        }
    }
    // As tmpUnit was NOT found in units, put tmpUnit into units.
    Units[unit_cnt].ID  = tmpID;
    Units[unit_cnt].len = len;
    for(int i=0; i<len; i++)
        Units[unit_cnt].string[i] = tmpUnit[i];
    Units[unit_cnt].string[len] = '\0';
    unit_cnt++;
    
    if(MAX_NUMBER_UNITS < unit_cnt){
        fprintf(stderr, "The number of unit counts %d exceeds %d\n", unit_cnt, MAX_NUMBER_UNITS);
        exit(EXIT_FAILURE);
    }
}

void put_into_GlobalUnits(char *tmpUnit){
    // Assume that tmpUnit is non-self-overlapping
    
    int len;
    for(len=0; tmpUnit[len] != '\0'; len++);
    //if( MAX_UNIT_LENGTH < len || len == 0) return; // Discard if the length exceeds the max.

    // Check if tmpUnit is in the database.
    int tmpID = min_quadratic_ID_with_rotation(tmpUnit, len);
    
    if(len < LONG_UNIT_LEN_TH){
        // Use exact match
        for(int i = 0; i < global_unit_cnt; i++)
            if(GlobalUnits[i].len == len && GlobalUnits[i].ID == tmpID){
                GlobalUnits[i].sumOccurrences++;
                return;
            }
    }else{
        // Allow some mismatches for long units
        for(int i = 0; i < global_unit_cnt; i++){
            int diff = MAX(len, GlobalUnits[i].len) - MIN(len, GlobalUnits[i].len);
            if( diff <= (int) (MAX(len, GlobalUnits[i].len) * MAX_DIS_RATIO) ){
                int match = rotate_match( tmpUnit, len, GlobalUnits[i].string, GlobalUnits[i].len);
                if(match == 1){
                    GlobalUnits[i].sumOccurrences++;
                    return;
                }
            }
        }
    }
    
    // As tmpUnit was NOT found in units, put tmpUnit into units.
    GlobalUnits[global_unit_cnt].ID  = tmpID;
    GlobalUnits[global_unit_cnt].len = len;
    for(int i=0; i<len; i++)
        GlobalUnits[global_unit_cnt].string[i] = tmpUnit[i];
    GlobalUnits[global_unit_cnt].string[len] = '\0';
    GlobalUnits[global_unit_cnt].sumOccurrences = 1;
    
    global_unit_cnt++;
    if(MAX_NUMBER_UNITS < global_unit_cnt){
        fprintf(stderr, "The number of unit counts %d exceeds %d\n", global_unit_cnt, MAX_NUMBER_UNITS);
        exit(EXIT_FAILURE);
    }
}

void print_GlobalUnits(){
    //printf("List of units\n");
    for(int i=0; i<global_unit_cnt; i++)
        printf("> frequent unit. freq. = %d\n%s\n", GlobalUnits[i].sumOccurrences, GlobalUnits[i].string);
    
}
