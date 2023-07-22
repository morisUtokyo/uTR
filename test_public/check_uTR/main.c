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

void rotate(int start, char *s1, int n1, char *rotated_s1){
    for(int i=0; i<n1; i++){
        rotated_s1[i] = s1[ (start + i ) % n1 ];
    }
    rotated_s1[n1] = '\0';
}

int rotate_match(char *s0, int n0, char *s1, int n1){
    // Return 1 if s0==s2, and 0 otherwise
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


int pattern_cmp(char *s1_pat, char *s2_pat, float allowance){
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
        if(rotate_match(unit1, strlen(unit1), unit2, strlen(unit2)) == 0 ){ match = 0; break; }
        //if(strcmp(unit1,unit2) != 0 ){ match = 0; break; }
        int d = abs(occ1 - occ2);
        int th = ceil(occ2 * allowance)+1;
        
        #ifdef DEBUG
        fprintf(stderr, "%s\t%d\t%s\t%d\t%d\t%d\n", unit1, occ1, unit2, occ2, d, th);
        #endif
        
        if(rotate_match(unit1, strlen(unit1), unit2, strlen(unit2)) == 1 && th < d ){ match = 0; break; }
        //if(strcmp(unit1,unit2) == 0 && th < d ){ match = 0; break; }
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

//#define DEBUG_main
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
                strcpy(s1_pat, ""); strcpy(s2_pat, "");
                ID2patterns(s, s1_pat, s2_pat);
                #ifdef DEBUG_main
                printf("%s\t%s\t%s\n", s, s1_pat, s2_pat);
                #endif
                if(strcmp(s1_pat,"") != 0 && strcmp(s2_pat,"") != 0 ){
                    i++;
                    if(pattern_cmp(s1_pat,s2_pat,allowance) == 1)  j++;
                }
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
