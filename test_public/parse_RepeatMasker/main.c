// Retrieve TRs from a RepeatMasker output file

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#define BLK 4096

int main(int argc, char *argv[])
{    
        
    int opt;
    int print_homopolymer = 0;
    char inputFile[BLK], outputFile[BLK];
    while ((opt = getopt(argc, argv, "i:o:h")) != -1) {
        switch(opt){
            case 'i':
                strcpy(inputFile,optarg);
                break;
            case 'o':
                strcpy(outputFile,optarg);
                break;
            case 'h':   // Print homopolymer
                print_homopolymer = 1;
                break;
            default:
                exit(EXIT_FAILURE);
        }
    }
    
    // Count the frequency of each SNV
    FILE *fp_in = fopen(inputFile, "r");
    FILE *fp_out = fopen(outputFile, "w");
    
    char *s = (char *)malloc(sizeof(char)*BLK);
    char a1[BLK], a2[BLK], a3[BLK], a4[BLK], ID[BLK], strand[BLK], a_repeat[BLK], repeat_unit[BLK], repeat_class[BLK], others[BLK];
    int  beginTR, endTR, a_left;
    
    int i=0;
    while (fgets(s, BLK, fp_in) != NULL) {
        if(i++ > 2){ // Skip the first three lines
            //printf("%s", s);
            sscanf(s, "%s %s %s %s %s %d %d (%d) %s %s %s %[^\n]",
                   a1,a2,a3,a4,ID,&beginTR,&endTR,&a_left,strand,a_repeat,repeat_class,others);
            if(strcmp(repeat_class,"Simple_repeat")==0){
                sscanf(a_repeat,"(%[^)])n",repeat_unit);
                strcpy(a_repeat,repeat_unit);
            }
            if(1 < strlen(a_repeat) || print_homopolymer == 1){
                fprintf(fp_out, "%s\t%d\t%d\t%d\t%d\t%s\t%s\n", ID,(endTR+a_left),beginTR,endTR,(endTR-beginTR+1),a_repeat,repeat_class);
            }
        }
    }
    fclose(fp_in);
    fclose(fp_out);

    free(s);

    return EXIT_SUCCESS;
}
