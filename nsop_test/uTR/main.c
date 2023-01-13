#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include "uTR.h"

void put_qualified_read(Read *currentRead, int i){
    int numReads, diameter, radius;
    // The annotation starts with a space " " !!!
    sscanf(currentRead->ID, " GroupSize = %d, Diameter = %d, RadiusFromCentroid = %d", &numReads, &diameter, &radius);
    Qreads[i].numReads     = numReads;
    Qreads[i].len          = currentRead->len;
    Qreads[i].numKeyUnits  = currentRead->numKeyUnits;
    Qreads[i].mosaic_mode  = currentRead->mosaic_mode;
}

int main(int argc, char *argv[])
{    
    char inputFile[500];    // For the input file name
    char repUnit[1000];     // For storing a representative unit
    char outputFile[500];   // For the output file name
    int inputFile_given = 0;
    int repUnit_given = 0;
    int print_time = 0;
    int print_EDDC = 0;
    int opt;
    while ((opt = getopt(argc, argv, "f:u:o:t")) != -1) {
        switch(opt){
            case 'f':
                strcpy(inputFile,optarg);   inputFile_given = 1;    break;
            case 'u':
                strcpy(repUnit, optarg);    repUnit_given= 1;    break;
            case 'o':
                strcpy(outputFile, optarg); print_EDDC = 1; break;
            case 't':
                print_time = 1; break;
            default:
                fprintf(stderr, "Usage: uTR -f <fasta file> (-u <representative unit string>) -o <output file>\n");
                exit(EXIT_FAILURE);
        }
    }
    if(inputFile_given == 0){
        fprintf(stderr, "Input file is not given.\n");
        exit(EXIT_FAILURE);
    }
    
    struct timeval s, e;    gettimeofday(&s, NULL);
    float time_get_non_self_overlapping_prefixes = time_coverage_by_units = time_set_cover_greedy = 0;
    
    // get non-self-overlapping units
    FILE *fp = init_handle_one_file(inputFile);
    Read *currentRead = malloc(sizeof(Read));
    malloc_Units(); malloc_GlobalUnits();
    for(int i=0;;i++){
        return_one_read(fp, currentRead);
        if(currentRead->len == 0) break;
        // get non-self-overlapping units
        get_non_self_overlapping_prefixes(currentRead->string);
    }
    int nsop_unit_cnt = unit_cnt;
    fclose(fp);
    free_Units();
    free_GlobalUnits();
    
    // get all substrings of length <= MAX_UNIT_LENGTH
    fp = init_handle_one_file(inputFile);
    currentRead = malloc(sizeof(Read));
    malloc_Units(); malloc_GlobalUnits();
    for(int i=0;;i++){
        return_one_read(fp, currentRead);
        if(currentRead->len == 0) break;
        // get all substrings of length <= MAX_UNIT_LENGTH
        put_all_substrings(currentRead->string);
    }
    int all_unit_cnt = unit_cnt;
    fclose(fp);
    free_Units();
    free_GlobalUnits();
    
    printf("\t%f\t%d\t%d\n", (float)nsop_unit_cnt/all_unit_cnt, nsop_unit_cnt, all_unit_cnt);
    
    return EXIT_SUCCESS;
}
