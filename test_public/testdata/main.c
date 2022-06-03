#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#define BLK 100000
//#define DEBUG

int main(int argc, char *argv[])
{
    char file1[100], file2[100];
    FILE *fp1, *fp2;
    
    if(argc==3){
        fp1 = fopen(argv[1], "r");
        fp2 = fopen(argv[2], "r");
    }else{
        fprintf(stderr, "Usage: accuracy (file1 name) (file2 name)\n");
        exit(EXIT_FAILURE);
    }
    
    
    //char s1[BLK+1], s2[BLK+1], s1_pat[BLK], s2_pat[BLK], s1_d[BLK], s2_d[BLK];
    
    char *s1 = (char *) malloc( sizeof(char) * BLK );
    char *s2 = (char *) malloc( sizeof(char) * BLK );
    char *s1_pat = (char *) malloc( sizeof(char) * BLK );
    char *s2_pat = (char *) malloc( sizeof(char) * BLK );
    char *s1_d = (char *) malloc( sizeof(char) * BLK );
    char *s2_d = (char *) malloc( sizeof(char) * BLK );
    
    int i=0, j=0, sum_len=0;
    for(; ;){
        if( fgets(s1, BLK, fp1) != NULL && fgets(s2, BLK, fp2) != NULL){
            if(s1[0] == '>' && s2[0] == '>'){
                sscanf(s1, "> %s %s", s1_pat, s1_d);
                sscanf(s2, "> %s %s", s2_pat, s2_d);
                #ifdef DEBUG
                printf("%s\t%s\t%s\t%s\n", s1_pat, s2_pat, s1_d, s2_d);
                #endif
                i++;
                if( strcmp(s1_pat, s2_pat) == 0 )
                    j++;
            }else{
                int k;
                for(k=0; s1[k]!='\0'; k++);
                sum_len += k;
            }
        }else{
            break;
        }
    }
    printf("%d %d %d\n", j, i, sum_len/i);
    
    free(s1); free(s2); free(s1_pat); free(s2_pat); free(s1_d); free(s2_d);
    return EXIT_SUCCESS;
}
