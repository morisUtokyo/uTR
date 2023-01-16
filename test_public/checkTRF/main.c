// Retrieve TRs from a RepeatMasker output file

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#define BLK 4096
#define largeBLK 40000
#define OCC_DIFF 1
#define numUnits 500
#define maxUnitLen 200

void rotate(int start, char *s1, int n1, char *rotated_s1){
    for(int i=0; i<n1; i++){
        rotated_s1[i] = s1[ (start + i ) % n1 ];
    }
    rotated_s1[n1] = '\0';
}

int rotate_match(char *s0, int n0, char *s1, int n1){
    char *rotated_s1 = malloc(sizeof(char) * (n1+1));
    
    int match=0;
    if(n0!=n1)
        match=0;
    else{
        // Find any rotation of s1 that matches s0
        for(int st=0; st<n0; st++){
            rotate(st, s1, n1, rotated_s1);
            if(strcmp(s0, rotated_s1) == 0){
                match=1; break;
            }
        }
    }
    free(rotated_s1);
    return(match);
}

struct unitStruct{
    char string[maxUnitLen+1];
    int  len;
};

int decompose(char *TRpattern, struct unitStruct *orig){
    
    int prev_pat_len = strlen(TRpattern);
    int n;
    for(int i=0; ; i++){
        sscanf( TRpattern, "%[^0-9]%d%s", orig[i].string, &orig[i].len, TRpattern );
        sscanf( orig[i].string, "(%[^)])", orig[i].string);
        if(prev_pat_len == strlen(TRpattern)){
            n=i+1; break;
        }
        else
            prev_pat_len = strlen(TRpattern);
    }
    return(n);
}

int print_orig_pred(FILE *fp_out, struct unitStruct *orig, struct unitStruct *pred, int len_orig, int len_pred, float allowance){

    int match = 1;
    if(len_orig == len_pred){
        for(int i=0; i<len_orig; i++){
            int string_match = rotate_match(
                    orig[i].string, strlen(orig[i].string),
                    pred[i].string, strlen(pred[i].string));
            int diff = abs(orig[i].len - pred[i].len);
            int th = ceil(orig[i].len * allowance)+1;
            
            if(string_match == 1 && diff <= th){} // allow at most "th" differences
            //if(string_match == 1 && diff <= OCC_DIFF){}
            else{
                match=0;
                break;
            }
        }
    }else
        match = 0;
    
    if(match == 1) fprintf(fp_out, "Y\t"); else fprintf(fp_out, "N\t");
        
    for(int i=0; i<len_orig; i++)
        fprintf(fp_out, "%s %d ", orig[i].string, orig[i].len);
    fprintf(fp_out, "\t");
    for(int i=0; i<len_pred; i++)
        fprintf(fp_out, "%s %d ", pred[i].string, pred[i].len);
    fprintf(fp_out, "\n");
    
    fflush(fp_out);
    
    return(match);
}

int main(int argc, char *argv[])
{    
        
    int opt;
    int print_homopolymer = 0;
    char inputFile[BLK], outputFile[BLK];
    float allowance = 0.1; // default
    while ((opt = getopt(argc, argv, "i:o:a:")) != -1) {
        switch(opt){
            case 'i':
                strcpy(inputFile,optarg);
                break;
            case 'a':
                sscanf(optarg, "%f", &allowance); break;
            case 'o':
                strcpy(outputFile,optarg);
                break;
            default:
                exit(EXIT_FAILURE);
        }
    }
    
    // Count the frequency of each SNV
    FILE *fp_in = fopen(inputFile, "r");
    FILE *fp_out = fopen(outputFile, "w");
    
    char *s = (char *) malloc( sizeof(char)*largeBLK );
    char head[BLK], TRpattern[largeBLK], unitString[BLK], trString[largeBLK];
    int  beginTR, endTR, periodSize, consensusSize, percentMatches, d1, d2, d3, d4, d5, d6;
    float copyNum, entropy;
    
    struct unitStruct orig[numUnits];
    struct unitStruct pred[numUnits];
    
    int first=1;
    int len_orig;
    int len_pred=0;
    int total_TRs = 0;
    int total_matches = 0;
    
    while (fgets(s, BLK, fp_in) != NULL) {
        sscanf(s, "%s %[^\n]", head, TRpattern);
        
        if(strcmp("@",head) == 0){
            if(first == 1)
                first = 0;
            else{
                total_matches += print_orig_pred(fp_out, orig, pred, len_orig, len_pred, allowance);
                total_TRs++;
                len_pred=0;
            }
            len_orig = decompose(TRpattern, orig);
        }else{
            sscanf(s, "%d %d %d %f %d %d %d %d %d %d %d %d %f %s %[^\n]", &beginTR, &endTR, &periodSize, &copyNum, &consensusSize, &percentMatches, &d1, &d2, &d3, &d4, &d5, &d6, &entropy, unitString, trString);
            int intCopyNum = (int) copyNum;
            
            //fprintf(stderr, "%d\n", beginTR);
            
            if(consensusSize < maxUnitLen){    // Ignore units of length maxUnitLen or more
                strcpy(pred[len_pred].string, unitString);
                pred[len_pred].len = intCopyNum;
                len_pred++;
            }
        }
    }
    total_matches += print_orig_pred(fp_out, orig, pred, len_orig, len_pred, allowance);
    total_TRs++;
    
    fprintf(fp_out, "%d %d\n", total_matches, total_TRs);
    
    fclose(fp_in);
    fclose(fp_out);

    free(s);

    return EXIT_SUCCESS;
}
