/*
 Copyright (c) 2019, Shinichi Morishita
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

#include <stdio.h>

// Key default parameters

// The maximum discrepancy ratio
#define MAX_DIS_RATIO_DEFAULT       0.3     // Setting this to 0.05 or 0.1 is too strict to overlook meaningful decompositions
float MAX_DIS_RATIO;

#define MIN_unit_occupancy_ratio  0.8
#define MIN_unit_length      2
#define TOP_k_units 10   // Print top k units for debugging

// Internal variables and data structures
#define MAX_NUMBER_READS    1000
#define MAX_READ_LENGTH     20000
#define MAX_ID_LENGTH       10000
#define BLK 4096
#define MAX_NUMBER_UNITS    10000
#define MAX_UNIT_LENGTH     20 // 100
int repUnitLen;
int read_cnt;
char *nextReadID;

// Parameters for processing long units
#define LONG_UNIT_LEN_TH    10 // 20
#define WINDOW_LEN          3 //5



// Parameters for controling bounded alignment among long units
#define MATCH       1
#define MISMATCH   -1
#define INDEL      -1
#define Mosaic_tandem_repeat 0
#define Mosaic_repeat 1

typedef struct{
    int len;
    char string[MAX_READ_LENGTH];
    char ID[MAX_ID_LENGTH];
    int  numKeyUnits;
    int  mosaic_mode; //Mosaic_tandem_repeat or Mosaic_repeat
    int  sumTandem;   // The sum of bases in mosaic TANDEM repeats
    char RegExpression[MAX_READ_LENGTH*5];
} Read;

typedef struct{
    int len;
    int numReads;
    int numKeyUnits;
    int mosaic_mode; //Mosaic_tandem_repeat or Mosaic_repeat
    int sumTandem;   // The sum of bases in mosaic TANDEM repeats
} QualifiedRead;
QualifiedRead *Qreads;

typedef struct{
    char *string;   // MAX_UNIT_LENGTH
    int ID;         // quadratic number
    int len;
    int sumOccurrences; // Total number of bases that match the unit
    int sumTandem;      // Total number of bps that match tandem repeat of the unit 
    int prio;
    int *covered;   // MAX_READ_LENGTH
} Unit; // The priority 1,2,... in the set of mosaic repeat units. -1 if the unit is not in the set.
Unit *Units;
int  unit_cnt;

Unit *GlobalUnits;
int  global_unit_cnt;

FILE*   init_handle_one_file(char *inputFile);
void    return_one_read(FILE *fp, Read *currentRead);

void    malloc_Units();
void    free_Units();
void    clear_Units_incrementally();

void    malloc_GlobalUnits();
void    free_GlobalUnits();
void    put_into_GlobalUnits(char *tmpUnit);
void    print_GlobalUnits();

char    int2char(int i);
void    put_repUnit(char *repUnit);

void    SA_IS(unsigned char *s, int *SA, int n, int K, int cs);

void    coverage_by_units(char *S, int MIN_number_repetitions);
void    set_cover_greedy(FILE *ofp, Read *currentRead, int MIN_number_repetitions);

// Interface between C and C++ functions
#ifndef __CSUB_H__
#define __CSUB_H__ 1
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    // C function called from a C++ function
    extern void put_unit(char *unit);
    extern void retain_top_k_units(int topK);
    extern int  char2int(char c);
    extern void match_bounded_DP_traceback(char *s0, int n0, char *s1, int n1, int *covered);
    // C++ function called from a C function
    extern void get_non_self_overlapping_prefixes(char *aString);
    extern int lzp(char *aString);
    extern void count_occurrences_long_unit(char *S, int n, int *SA, int *C, int **OCC, char *unit, int unitLen, int *tmpCovered, int MIN_number_repetitions );
#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* __CSUB_H__ */


// External functions
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define DIFF(x, y) ((x) > (y) ? ((x) - (y)) : ((y) - (x)))

// Debug Mode
//#define DEBUG_EDDC_mode  // Show the frequency of each init.


float time_get_non_self_overlapping_prefixes, time_coverage_by_units, time_set_cover_greedy;


