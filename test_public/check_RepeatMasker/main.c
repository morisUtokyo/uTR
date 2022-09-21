// rTR: Retrieve TRs from a RepeatMasker output file

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#define BLK 4096
#define DEBUG

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
    char string[100];
    int  len;
};

void trim_dup(char *pat){
    int n=strlen(pat);
    for(int i=0; i<n; i++){
        if(pat[i]=='_'){
            pat[i]='\0';
            break;
        }
    }
}

int process_one_TR(char *prev_pat, struct unitStruct *unit_list_RM, int j, int print_each_TR){
    
    // For example, remove "_1" in (AC)174(AG)440_1 that is appended to distinguish its duplicate (AC)174(AG)440
    trim_dup(prev_pat);
    
    struct unitStruct unit_list[10];
    
    int prev_pat_len = strlen(prev_pat);
    int n;
    for(int i=0; ; i++){
        sscanf(prev_pat, "%[^0-9]%d%s", unit_list[i].string,&unit_list[i].len,prev_pat);
        sscanf(unit_list[i].string, "(%[^)])", unit_list[i].string);
        if(prev_pat_len == strlen(prev_pat)){
            n=i+1; break;
        }
        else
            prev_pat_len = strlen(prev_pat);
    }
    int match=1;
    if(n != j)
        match=0;
    else{
        for(int i=0; i<n; i++){
            int diff = abs(unit_list[i].len - unit_list_RM[i].len);
            int string_match = rotate_match(
                    unit_list[i].string,
                    strlen(unit_list[i].string),
                    unit_list_RM[i].string,
                    strlen(unit_list_RM[i].string));
            if(string_match == 1 && diff <= 1){}
            else{
                match=0; break;
            }
        }
    }
    if(print_each_TR == 1){
        if(match==1)
            printf("Yes ");
        else
            printf("No  ");
        for(int i=0; i<n; i++)
            printf("%s %d ", unit_list[i].string,unit_list[i].len);
        printf("\t");
        for(int i=0; i<j; i++)
            printf("%s %d ",
                unit_list_RM[i].string,unit_list_RM[i].len);
        printf("\n");
    }
    return(match);
}

int main(int argc, char *argv[])
{    
        
    int opt;
    char inputFile[BLK], outputFile[BLK];
    int print_each_TR=0;
    while ((opt = getopt(argc, argv, "i:o:p")) != -1) {
        switch(opt){
            case 'i':
                strcpy(inputFile,optarg); break;
            //case 'o':
            //    strcpy(outputFile,optarg); break;
            case 'p':
                print_each_TR = 1; break;
            default:
                exit(EXIT_FAILURE);
        }
    }
    
    // Count the frequency of each SNV
    FILE *fp_in = fopen(inputFile, "r");
    
    //char *s = (char *)malloc(sizeof(char)*BLK);
    char s[BLK];
    char pat[BLK], unit[BLK], annotation[BLK];
    int  lenTR, beginUnit, endUnit, lenUnit;
    char prev_pat[BLK];

    //struct unitStruct unit_list[10];
    struct unitStruct unit_list_RM[10];
    
    int count_all = 0;
    int count_match = 0;
    int match;
    int j=0;
    strcpy(prev_pat,"");
    while (fgets(s, BLK, fp_in) != NULL) {
        sscanf(s, "%s %d %d %d %d %s %[^\n]",
               pat,&lenTR,&beginUnit,&endUnit,&lenUnit,unit,annotation);
        if(strcmp(prev_pat,pat)!=0 && strcmp(prev_pat,"")!=0)
        {
            match = process_one_TR(prev_pat, unit_list_RM, j, print_each_TR);
            strcpy(prev_pat, pat);
            j=0;
            count_all++;
            if(match == 1) count_match++;
        }
        strcpy(unit_list_RM[j].string, unit);
        unit_list_RM[j].len = lenUnit/ strlen(unit_list_RM[j].string);
        strcpy(prev_pat, pat);
        j++;
    }
    match = process_one_TR(prev_pat, unit_list_RM, j, print_each_TR);
    count_all++;
    if(match == 1) count_match++;
    
    printf("%d\t%d\n", count_match, count_all);
    
    fclose(fp_in);

    return EXIT_SUCCESS;
}
