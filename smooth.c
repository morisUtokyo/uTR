#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include "uTR.h"

int majority(int *blocks, int start, int end, int *unit_freq, int numKeyUnits){
    for(int i=0; i<numKeyUnits; i++)   unit_freq[i]=0;
    for(int i=start; i<= end; i++){
        if(0 <= blocks[i] && blocks[i] < numKeyUnits)
            unit_freq[ blocks[i] ]++;
    }
    int max=0;
    int max_i=0;
    for(int i=0; i<numKeyUnits; i++)
        if(max < unit_freq[i]){
            max = unit_freq[i];
            max_i = i;
        }
    return(max_i);
}

#define minWindow 11
#define maxWindow 41
int window_length(int len){
    int threshold = len/5;
    if( threshold%2 == 0 ) threshold++; // threshold must be an odd number.
    int wlen;
    if(threshold < minWindow)       wlen = minWindow;
    else if(maxWindow < threshold)  wlen = maxWindow;
    else                            wlen = threshold;
    return(wlen);
}

void smooth_sub(int *input_blocks, int len, int numKeyUnits){
    int *unit_freq = (int *) malloc( sizeof(int) * numKeyUnits);
    //int *unit_freq = (int *) malloc( sizeof(int) * MAX_NUMBER_UNITS);
    int *output_blocks = (int *) malloc( sizeof(int) * (len+1) );
    
    int wlen = window_length(len);
    for(int i=0; i < len; i++){
        int start_i = i - wlen/2;
        int end_i = i + wlen/2;
        if( start_i < 0){
            start_i = 0;
            end_i = MIN(len, wlen)-1;
        }
        if( len-1 < end_i){
            end_i = len-1;
            start_i = MAX(0, end_i - wlen +1);
        }
        if(start_i < 0 || len-1 < end_i){
            fprintf(stderr, "Out of range: start_i = %d, end_i = %d\n", start_i, end_i);
            exit(EXIT_FAILURE);
        }
        output_blocks[i] = majority(input_blocks, start_i, end_i, unit_freq, numKeyUnits);
    }
    for(int i=0; i < len; i++)
        input_blocks[i] = output_blocks[i];
    
    free(unit_freq);
    free(output_blocks);
}

void smooth(int *input_blocks, int len, int numKeyUnits){
    // Smooth the blocks twice
    smooth_sub(input_blocks, len, numKeyUnits);
    smooth_sub(input_blocks, len, numKeyUnits);
}

