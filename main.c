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
    // To parse CentroidReadName = B483,read1, use %[^,] in place of %s
    int return_sscanf = sscanf(currentRead->ID, " GroupSize = %d, Diameter = %d, RadiusFromCentroid = %d, CentroidReadName = %[^,],%[^,], CentroidReadLength =", &numReads, &diameter, &radius, Qreads[i].individualID, Qreads[i].readID);
    
    if(return_sscanf == 0){
        numReads = 1;
        strcpy( Qreads[i].individualID, "NA");
        strcpy( Qreads[i].readID, "NA");
    }
    Qreads[i].numReads     = numReads;
    Qreads[i].len          = currentRead->len;
    Qreads[i].numKeyUnits  = currentRead->numKeyUnits;
    //Qreads[i].mosaic_mode  = currentRead->mosaic_mode;
    Qreads[i].sumTandem    = currentRead->sumTandem;
    Qreads[i].discrepancy_ratio = currentRead->discrepancy_ratio;
    strcpy(Qreads[i].RegExpression, currentRead->RegExpression);
    strcpy(Qreads[i].decomposition, currentRead->decomposition);
}

void build_Haps(char *hapFile){
    FILE *hfp;
    char s[BLK+1];
    hap_cnt = 0;

    hfp = fopen(hapFile, "r");
    while (fgets(s, BLK, hfp) != NULL) hap_cnt++;
    fclose(hfp);
    
    Haps = (Hap *) malloc( sizeof(Hap) * hap_cnt );
    hfp = fopen(hapFile, "r");
    int i=0;
    while (fgets(s, BLK, hfp) != NULL){
        sscanf(s, "%s\t%d\t%s", Haps[i].individualID, &Haps[i].readID, Haps[i].pairHaps);
        i++;
    }
    fclose(hfp);
}

void print_Haps(){
    for(int i=0; i<hap_cnt; i++)
        printf("%s\t%d\t%s\n", Haps[i].individualID, Haps[i].readID, Haps[i].pairHaps);
}

char *query_hap(char *individualID, int readID){
    for(int i=0; i<hap_cnt; i++){
        if(strcmp(individualID, Haps[i].individualID) == 0 && readID == Haps[i].readID ){
            return(Haps[i].pairHaps);
        }
    }
    return(NULL);
}

//#define DEBUG_main
int main(int argc, char *argv[])
{    
    char inputFile[500];    // For the input file name
    char repUnit[1000];     // For storing a representative unit
    char hapFile[500];      // For the haplotype file name
    char outputFile[500];   // For the output file name
    char tableFile[500];    // For the table file name (repRead, information)
    char statTRpatFile[500];// For the file name of statistics with tandem repeat patterns
    char locus[500];        // For storing the locus ID (e.g., 1234:chr1:456-567)
    int inputFile_given = 0;
    int repUnit_given = 0;
    int hapFile_given = 0;
    int print_time = 0;
    int print_EDDC = 0;
    int print_table = 0;
    int print_statTRpat = 0;
    int hide_IDs = 1;   // hide IDs by default
    int regular_expression_only=0;  // Donot print decomposition
    int print_locus=0;
    int print_input_annotation_as_it_is=0;
    int print_stdout_decomp = 0;
    int mode_longer_TRs = 0;
    int mode_smooth = 1;

    int opt;
    
    MAX_DIS_RATIO = MAX_DIS_RATIO_DEFAULT;
    float ratio;
    
    while ((opt = getopt(argc, argv, "l:f:u:h:o:i:r:p:stdaxyz")) != -1) {
        switch(opt){
            case 'l': strcpy(locus,optarg);       print_locus = 1;    break;
            case 'f': strcpy(inputFile,optarg);   inputFile_given = 1;    break;
            case 'u': strcpy(repUnit, optarg);    repUnit_given= 1;    break;
            case 'h': strcpy(hapFile,optarg); build_Haps(hapFile); hapFile_given = 1;    break;
            case 'o': strcpy(outputFile, optarg); print_EDDC = 1; break;
            case 'i': strcpy(tableFile, optarg);  print_table = 1; break;
            case 'p': strcpy(statTRpatFile, optarg);  print_statTRpat = 1; break;
            case 'r': sscanf(optarg, "%f", &ratio);   MAX_DIS_RATIO = ratio;  break;
            case 't': print_time = 1; break;
            case 's': hide_IDs = 0; break;  // Print IDs
            case 'd': regular_expression_only = 1; break;
            case 'a': print_input_annotation_as_it_is = 1; break;
            case 'x': print_stdout_decomp = 1; break;
            case 'y': mode_smooth = 0; break;   // mode to avoid smoothing units
            case 'z': mode_longer_TRs = 1; break;   // mode for finding longer TR units
            default:
                fprintf(stderr, "Usage: uTR -l <locus info> -f <fasta file> (-u <representative unit string> -h <haplotype file> -i <table file> -r <maximum discrepancy ratio> -t (print wall clock time) -s (print IDs) -d (do not print decomposition) -a (print the input annotation) -y (Print more complex patterns with fewer mismatches) -z (Print patterns with longer units) -o <EDDC output file> -p <TR patterns>\n");
                exit(EXIT_FAILURE);
        }
    }
    if(inputFile_given == 0){
        fprintf(stderr, "Input file is not given.\n");
        exit(EXIT_FAILURE);
    }
    
    struct timeval s, e;
    gettimeofday(&s, NULL);
    float time_get_non_self_overlapping_prefixes = time_coverage_by_units = time_set_cover_greedy = 0;
    
    int num_reads = count_reads(inputFile);
    FILE *fp = init_handle_one_file(inputFile);
    FILE *ofp;
    if(print_EDDC == 1) ofp = fopen(outputFile, "w"); else ofp = stdout;
    FILE *tfp;
    if(print_table== 1) tfp = fopen(tableFile, "w");  else tfp = stdout;
    FILE *pfp;
    if(print_statTRpat== 1) pfp = fopen(statTRpatFile, "w");  else pfp = stdout;
    
    Read *currentRead = malloc(sizeof(Read));
    if(currentRead==NULL){ fprintf(stderr, "Failure to malloc currentRead\n"); exit(EXIT_FAILURE); }
    malloc_Units(num_reads);
    //malloc_Units();
    malloc_GlobalVars();
    if(repUnit_given == 1)  put_repUnit(repUnit);
    for(int i=0; repUnit[i]!='\0'; i++, repUnitLen=i){}
    
    int numQualifiedReads = 0;
    char *decomposition = (char *) malloc(sizeof(char) * MAX_READ_LENGTH);
    char *stat = (char *) malloc(sizeof(char) * MAX_READ_LENGTH);
    char *stat_table = (char *) malloc(sizeof(char) * MAX_READ_LENGTH);
    char individualID[MAX_ID_LENGTH];
    char readID[MAX_ID_LENGTH];
    // Feed each read and make the database of units
    for(int i=0;;i++){
        return_one_read(fp, currentRead);
#ifdef DEBUG_main
        fprintf(stderr, "Scanned one read\n");
#endif
        if(currentRead->len == 0) break;
        
        int MIN_reps;
        //if(mode_smooth == 1) MIN_reps=2; else MIN_reps=1;
        for(MIN_reps=2; 0 < MIN_reps; MIN_reps--){
            if(MIN_reps == 2)
                currentRead->mosaic_mode = Mosaic_tandem_repeat;
            else if(MIN_reps == 1)
                currentRead->mosaic_mode = Mosaic_repeat;
            
            set_cover_greedy(currentRead, MIN_reps, mode_longer_TRs, mode_smooth); // Select TOP_k_units units in a batch manner and put them into prio2unit, a local array
            // The annotation starts with a space " " !!!
        
            int numReads, diameter, radius, readlen;
            int return_sscanf = sscanf(currentRead->ID, " GroupSize = %d, Diameter = %d, RadiusFromCentroid = %d, CentroidReadName = %[^,],%[^,], CentroidReadLength = %d", &numReads, &diameter, &radius, individualID, readID, &readlen);
            if(return_sscanf == 0){
                strcpy(individualID, "NA"); strcpy(readID, "NA"); numReads=1;
            }
            
            if(0 < currentRead->numKeyUnits){ // Print reads with primary units
                sprintf(stat, "");
                sprintf(stat_table, "");

                put_qualified_read(currentRead, numQualifiedReads);
                numQualifiedReads++;
                
                if(hide_IDs == 0)
                    sprintf(stat, "(%s,%s,%d,%d,%d,%4.3f)", individualID, readID, currentRead->len, numReads, currentRead->numKeyUnits, currentRead->discrepancy_ratio);
                else
                    sprintf(stat, "(%d,%d,%4.3f)", currentRead->len, numReads, currentRead->discrepancy_ratio);
                sprintf(stat_table, "%s,%s (%d,%d,%d,%4.3f) ", individualID, readID, currentRead->len, numReads, currentRead->numKeyUnits, currentRead->discrepancy_ratio);
                
                if(print_table == 1){
                    fprintf(tfp, "%s", stat_table);
                    fprintf(tfp, "%s ", currentRead->decomposition);
                    fprintf(tfp,"\n");
                    fflush(tfp);
                }
                
                // Print annotation of the focal read
                if(print_EDDC == 1){    // Print the focal read
                    fprintf(ofp, "> ");
                    if(print_locus == 1)    fprintf(ofp, "#Locus %s ", locus);
                    fprintf(ofp, "#Info %s ", stat);
                    
                    fprintf(ofp, "#Pat %s ", currentRead->RegExpression);
                    // Print a pair of nearest SNVs
                    int int_readID;
                    sscanf(readID, "read%d", &int_readID);
                    if(hapFile_given == 1)
                        fprintf(ofp, "#Hap <%s> ", query_hap(individualID, int_readID) );
                    if(regular_expression_only  == 0)
                        fprintf(ofp, "#Decomp %s ", currentRead->decomposition);
                    // Print the input annotation
                    if(print_input_annotation_as_it_is == 1)
                        fprintf(ofp, "#Annotation %s ", currentRead->ID);
                    // Print the string with TR
                    fprintf(ofp, "\n%s\n", currentRead->string);
                    fflush(ofp);
                }
                break;
            }
        }
    }
    if(0 < numQualifiedReads){
	if(print_stdout_decomp == 1)
	        printf(" #haplotypes=%d", numQualifiedReads);
        for(int i=0; i<numQualifiedReads; i++){
	    if(print_stdout_decomp == 1)
	            printf(" (%s,%s,%d,%d,%d,%3.2f) %s", Qreads[i].individualID, Qreads[i].readID, Qreads[i].len, Qreads[i].numReads, Qreads[i].numKeyUnits, Qreads[i].discrepancy_ratio, Qreads[i].decomposition);
        }
        // Print the statistics with TR patterns if(print_statTRpat== 1)
        if(print_statTRpat == 1){
            fprintf(pfp, " #haplotypes=%d", numQualifiedReads);
            if(hide_IDs == 1){
                for(int i=0; i<numQualifiedReads; i++){
                    fprintf(pfp, " (%d,%d,%3.2f) %s", Qreads[i].len, Qreads[i].numReads, Qreads[i].discrepancy_ratio, Qreads[i].RegExpression);
                }
            }else{
                for(int i=0; i<numQualifiedReads; i++){
                    fprintf(pfp, " (%s,%s,%d,%d,%3.2f) %s", Qreads[i].individualID, Qreads[i].readID, Qreads[i].len, Qreads[i].numReads, Qreads[i].discrepancy_ratio, Qreads[i].RegExpression);
                }
            }
        }
        if(print_stdout_decomp == 1){
            printf(" #units=%d ", global_unit_cnt);
        }
        for(int i=0; i<global_unit_cnt; i++){
            if(print_stdout_decomp == 1){
                printf("(%s,%d) ", GlobalUnits[i].string, GlobalUnits[i].sumOccurrences);
            }
            if(print_EDDC == 1){ // Print primary units for EDDC algorithm
            //if(print_EDDC == 1 && print_input_annotation_as_it_is == 0){ // Print primary units for EDDC algorithm unless the mode of testing accuracy
                fprintf(ofp, ">");
                if(print_locus == 1)    fprintf(ofp, "#Locus %s", locus);
                fprintf(ofp, " frequent unit. freq. = %d\n%s\n", GlobalUnits[i].sumOccurrences, GlobalUnits[i].string);
                fflush(ofp);
            }
        }
        if(print_stdout_decomp == 1){ printf("\n"); }
    }else{
        if(print_stdout_decomp == 1){ printf(" #haplotypes=0\n"); }
    }
    fclose(fp);
    if(print_EDDC == 1) fclose(ofp);
    if(print_table== 1) fclose(tfp);
    free_Units();
    free_GlobalVars();
    

    if(print_time == 1){
        gettimeofday(&e, NULL);
        float time_all = (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
        fprintf(stderr, "%s\t%f sec\n", inputFile, time_all);
    }
    
    return EXIT_SUCCESS;
}
