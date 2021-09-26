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
    Qreads[i].sumTandem    = currentRead->sumTandem;
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
    
    MAX_DIS_RATIO = MAX_DIS_RATIO_DEFAULT;
    float ratio;
    
    while ((opt = getopt(argc, argv, "f:u:o:r:t")) != -1) {
        switch(opt){
            case 'f':
                strcpy(inputFile,optarg);   inputFile_given = 1;    break;
            case 'u':
                strcpy(repUnit, optarg);    repUnit_given= 1;    break;
            case 'o':
                strcpy(outputFile, optarg); print_EDDC = 1; break;
            case 't':
                print_time = 1; break;
            case 'r':
                sscanf(optarg, "%f", &ratio);
                MAX_DIS_RATIO = ratio;
                break;
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
    
    FILE *fp = init_handle_one_file(inputFile);
    FILE *ofp;
    if(print_EDDC == 1)
        ofp = fopen(outputFile, "w");
    else
        ofp = stdout;
    
    Read *currentRead = malloc(sizeof(Read));
    malloc_Units(); malloc_GlobalUnits();
    if(repUnit_given == 1)  put_repUnit(repUnit);
    for(int i=0; repUnit[i]!='\0'; i++, repUnitLen=i){}
    
    int numQualifiedReads = 0;
    // Feed each read and make the database of units
    for(int i=0;;i++){
        return_one_read(fp, currentRead);
        if(currentRead->len == 0) break;

        for(int j=0; j<2; j++){
            int MIN_reps = 2-j; // First MTR, then MR
            // mosaic tandem repeats = 2, mosaic repeats = 1
            clear_Units_incrementally();
            get_non_self_overlapping_prefixes(currentRead->string);
            coverage_by_units(currentRead->string, MIN_reps);
            //int numKeyUnits = 0;
            if(0 < unit_cnt)
                set_cover_greedy(ofp, currentRead, MIN_reps);
            int numKeyUnits = currentRead->numKeyUnits;
            if(MIN_reps == 2)
                currentRead->mosaic_mode = Mosaic_tandem_repeat;
            else if(MIN_reps == 1)
                currentRead->mosaic_mode = Mosaic_repeat;
            int numReads, diameter, radius, readlen;
            char readID[1000];
            // The annotation starts with a space " " !!!
            sscanf(currentRead->ID, " GroupSize = %d, Diameter = %d, RadiusFromCentroid = %d, CentroidReadName = %s, CentroidReadLength = %d", &numReads, &diameter, &radius, readID, &readlen);
            
            if(0 < numKeyUnits){ // Print reads with primary units
                if(print_EDDC == 1){ // Print annotation of the focal read
                    fprintf(ofp, "> (%d,%d,%d,", currentRead->len, numReads, numKeyUnits);
                    if(MIN_reps == 2)
                        fprintf(ofp, "MTR,");// mosaic tandem repeats
                    else if(MIN_reps == 1)
                        fprintf(ofp, "MR,"); // mosaic repeats
                    fprintf(ofp, "%d),%s S=", currentRead->sumTandem, readID);
                }
                // 1-origin index for prio !!!
                int prevSum = 0;
                for(int prio=1; prio <= MIN(TOP_k_units, unit_cnt); prio++){
                    for(int j=0; j<unit_cnt; j++)
                        if(prio == Units[j].prio){  // Note that 1 <= Units[j].prio
                            put_into_GlobalUnits(Units[j].string);
                            if(print_EDDC == 1)
                                fprintf(ofp, "(%d,%s,%d,%d,%d),", prio, Units[j].string, Units[j].len,  Units[j].sumOccurrences - prevSum, Units[j].sumTandem);
                            prevSum = Units[j].sumOccurrences;
                        }
                }
                put_qualified_read(currentRead, numQualifiedReads);
                numQualifiedReads++;
                if(print_EDDC == 1){// Print the focal read
                    fprintf(ofp, " D=%s\n", currentRead->RegExpression);
                    fprintf(ofp, "%s\n",   currentRead->string);
                }
                break;
            }
        }
    }
    if(0 < numQualifiedReads){
        printf(" #haplotypes=%d ", numQualifiedReads);
        for(int i=0; i<numQualifiedReads; i++){
            printf("(%d,%d,%d,", Qreads[i].len, Qreads[i].numReads, Qreads[i].numKeyUnits);
            if(Qreads[i].mosaic_mode == Mosaic_tandem_repeat) printf("MTR,");
            else printf("MR,");
            printf("%d) ", Qreads[i].sumTandem);
        }
        printf("#units=%d ", global_unit_cnt);
        for(int i=0; i<global_unit_cnt; i++){
            printf("(%s,%d) ", GlobalUnits[i].string, GlobalUnits[i].sumOccurrences);
            if(print_EDDC == 1) // Print primary units for EDDC algorithm
                fprintf(ofp, "> frequent unit. freq. = %d\n%s\n", GlobalUnits[i].sumOccurrences, GlobalUnits[i].string);
        }
        printf("\n");
    }else{
        printf(" #haplotypes=0\n");
    }
    fclose(fp);
    if(print_EDDC == 1) fclose(ofp);
    free_Units();
    free_GlobalUnits();
    

    if(print_time == 1){
        gettimeofday(&e, NULL);
        float time_all = (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
        fprintf(stderr, "file=%s\ttime=%f\tunit_cnt=%d\tMAX_DIS_RATIO=%.3f\n", inputFile,time_all, unit_cnt, MAX_DIS_RATIO);
    }
    
    return EXIT_SUCCESS;
}
