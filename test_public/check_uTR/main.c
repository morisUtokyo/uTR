#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#define BLK 100000
#define OCC_DIFF 1

//#define DEBUG
int pattern_cmp(char *s1_pat, char *s2_pat, float allowance){
//int pattern_cmp(char *s1_pat, char *s2_pat, int occ_diff){
    char unit1[100], unit2[100];
    int occ1, occ2;
    
    #ifdef DEBUG
    printf("%s\t%s\n", s1_pat, s2_pat);
    #endif

    int match = 1;
    int diff1 = strlen(s1_pat);
    int diff2 = strlen(s2_pat);
    while(0 < diff1  && 0 < diff2){
        diff1 = strlen(s1_pat);
        diff2 = strlen(s2_pat);
        sscanf(s1_pat, "%s %d%[^\0]%*c", unit1, &occ1, s1_pat);
        sscanf(s2_pat, "%s %d%[^\0]%*c", unit2, &occ2, s2_pat);
        // The two sets of units are different.
        if(strcmp(unit1,unit2) != 0 ){
            match = 0; break;
        }
        int d = abs(occ1 - occ2);
        int th = ceil(occ2 * allowance)+1;
        
        #ifdef DEBUG
        fprintf(stderr, "%s\t%d\t%s\t%d\t%d\t%d\n", unit1, occ1, unit2, occ2, d, th);
        #endif
        
        if(strcmp(unit1,unit2) == 0 && th < d ){
        //if(strcmp(unit1,unit2) == 0 && occ_diff < d ){
            match = 0; break;
        }
        diff1 -= strlen(s1_pat);
        diff2 -= strlen(s2_pat);
        
        if(diff1 == 0 && diff2 == 0){   // The two sets of units are identical.
            match = 1; break;
        }else if(diff1 == 0 || diff2 == 0){ // Different.
            match = 0; break;
        }
    }
    if(match == 1) return(1); else return(0);
}

void ID2patterns(char *ID, char *s1_pat, char *s2_pat){
    
    char *tmpS      = (char *) malloc(sizeof(char) * BLK);
    char *content   = (char *) malloc(sizeof(char) * BLK);
    char tag[100];
    
    strcpy(tmpS, ID);

    int prev_tmpS_len= strlen(tmpS);
    int detected_a_pattern=0;
    for(;;){
        sscanf(tmpS, "#%s %s %[^\0]", tag, content, tmpS);
        if(strcmp(tag, "Annotation") == 0){
            strcpy(s1_pat, content);
            detected_a_pattern = 1;
        }
        if(strcmp(tag, "Pat") == 0){
            strcpy(s2_pat, content);
            detected_a_pattern = 1;
        }
        if(prev_tmpS_len == (int)strlen(tmpS))  break;
        else    prev_tmpS_len = strlen(tmpS);
    }
    free(tmpS); free(content);
    if(detected_a_pattern == 0){
        strcpy(s1_pat, ""); strcpy(s2_pat, "");
    }
    // Replace < and > with ( and ), respectively.
    for(int j=0; j<strlen(s1_pat); j++)
        if(s1_pat[j]=='<' || s1_pat[j]=='(' || s1_pat[j]=='>' || s1_pat[j]==')' )  s1_pat[j]=' ';
    for(int j=0; j<strlen(s2_pat); j++)
        if(s2_pat[j]=='<' || s2_pat[j]=='(' || s2_pat[j]=='>' || s2_pat[j]==')' )  s2_pat[j]=' ';
}


int main(int argc, char *argv[])
{
    int opt;
    char inputFile[1000];
    float allowance = 0.1; // default
    while ((opt = getopt(argc, argv, "i:a:")) != -1) {
        switch(opt){
            case 'i':
                strcpy(inputFile, optarg); break;
            case 'a':
                sscanf(optarg, "%f", &allowance); break;
            default:
                fprintf(stderr, "Usage: -i (file name) -a allowance (e.g., 0.01 by default)\n");
                exit(EXIT_FAILURE);
        }
    }
    FILE *fp = fopen(inputFile, "r");
    
    // Input fasta file begins with the form:
    // > (60,1,0.00) <AAG>12<AG>12 (AAG)12(AG)12
   
    char *s     = (char *) malloc( sizeof(char) * BLK );
    char *s1_pat= (char *) malloc( sizeof(char) * BLK );
    char *s2_pat= (char *) malloc( sizeof(char) * BLK );

    int i=0, j=0, sum_len=0;
    for(; ;){
        if( fgets(s, BLK, fp) != NULL){
            if(s[0] == '>'){
                sscanf(s, "> %[^\0]", s);
                ID2patterns(s, s1_pat, s2_pat);
                //printf("%s\t%s\t%s\n", s, s1_pat, s2_pat);
                i++;
                if(pattern_cmp(s1_pat,s2_pat,allowance) == 1)  j++;
                
            }else{
                int k;
                for(k=0; s[k]!='\0'; k++);
                sum_len += k;
            }
        }else{
            break;
        }
    }
    printf("%d %d %d\n", j, i, sum_len/i);
    
    free(s); free(s1_pat); free(s2_pat);
    return EXIT_SUCCESS;
}
