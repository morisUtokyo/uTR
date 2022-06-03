#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include "MT.h"     // Use the Mersenne Twister.

#define NUM_MOSAIC_TRS 1000
#define ERROR_RATE 0 // 0.05

#define MAX_UNIT_OCC 20
#define MIN_UNIT_OCC 10

#define ST_LEN 100000

//#define DEBUG_errors

void print_units(char *st, char *unit, int n){
    for(int i=0; i<n; i++)
        sprintf(st, "%s%s", st, unit); // append the give unit n times.
}

char another_char(char c){
    char c1;
    for(;;){
        switch(genrand_int32()%4){
            case 0: c1='A';
            case 1: c1='C';
            case 2: c1='G';
            case 3: c1='T';
        }
        if(c!=c1) return(c1);
    }
}

float put_errors(char *st, float error_ratio){
    char *tmp_st = (char *) malloc( sizeof(char) * ST_LEN );
    int i, j;
    for(i=0; st[i]!='\0'; i++)
        tmp_st[i] = st[i];
    tmp_st[i]='\0';
    
    int num_errors=0;

    for(i=0, j=0; tmp_st[i]!='\0'; ){
        if( (genrand_int32() % 10000) < 10000*(error_ratio) ){
            switch(genrand_int32()%3){
                case 0: st[++j] = tmp_st[i++]; j++; break; // insertion
                case 1: i++; break; // deletion
                case 2: st[j++] = another_char(tmp_st[i++]); num_errors++; break;
            }
        }else
            st[j++] = tmp_st[i++];
    }
    free(tmp_st);
    float discrepancy = (float)num_errors / i;
    return(discrepancy);
}


int main(int argc, char *argv[])
{
    int num_units, max_unit_occ=MAX_UNIT_OCC, min_unit_occ=MIN_UNIT_OCC;
    int num_TRs = NUM_MOSAIC_TRS;
    float error_ratio = 0;
    int opt;
    while ((opt = getopt(argc, argv, "k:l:m:n:e:")) != -1) {
        switch(opt){
            case 'k':
                sscanf(optarg, "%d", &min_unit_occ);
                //fprintf(stderr, "Minimum number of unit occrrences = %d\n", min_unit_occ);
                break;
            case 'l':
                sscanf(optarg, "%d", &max_unit_occ);
                //fprintf(stderr, "Maximum number of unit occrrences = %d\n", max_unit_occ);
                break;
            case 'm':
                sscanf(optarg, "%d", &num_units);
                //fprintf(stderr, "Number of units = %d\n", num_units);
                break;
            case 'n':
                sscanf(optarg, "%d", &num_TRs);
                //fprintf(stderr, "Number of TRs = %d\n", num_TRs);
                break;
            case 'e':
                sscanf(optarg, "%f", &error_ratio);
                //fprintf(stderr, "Error ratio = %3.2f\n", error_ratio);
                break;
            default:
                fprintf(stderr, "Usage: gen -f (number of units) units -n (number of TRs) -e (error rate) \n");
                exit(EXIT_FAILURE);
        }
    }
    
    if(max_unit_occ - min_unit_occ < 1){
        fprintf(stderr, "max_unit_occ must be greater than min_unit_occ\n");
        exit(EXIT_FAILURE);
    }

    int *randNums  = (int *) malloc( sizeof(int) * num_units );
    char *st = (char *) malloc( sizeof(char) * ST_LEN );

    for(int j=0; j<num_TRs; ){
        for (int i=0; i<num_units; i++){
            randNums[i] = (genrand_int32()%(max_unit_occ - min_unit_occ)) + min_unit_occ;
        }
        
        sprintf(st, "");
        int start_units = argc - num_units; // The first unit is argv[i + start_units]
        for (int i=0; i<num_units; i++)
            print_units(st, argv[i + start_units], randNums[i]);
        float discrepancy = put_errors(st, error_ratio);
        
        printf("> ");
        for (int i=0; i<num_units; i++)
            printf("<%s>%d", argv[i + start_units], randNums[i]);
        printf(" discrepancy=%3.2f\n", discrepancy);
        printf("%s\n", st);
        j++;
    }

    free(randNums);
    free(st);
    return EXIT_SUCCESS;
}
