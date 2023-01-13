//-------------------------------------------------------------
// Induced Sorting
// Reference: Nong, Ge, Sen Zhang, and Wai Hong Chan. "Two efficient algorithms for linear time suffix array construction.” IEEE Transactions on Computers, 60.10 (2011): 1471-1484.
// To compile, perform: gcc -std=c99 SAIS.c -o SAIS
// To run, SALS　<length of the string>. e.g., SAIS 100
//-------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include "uTR.h"

char int2char(int i){
    char ans;
    switch(i){
        case 1: ans = 'A'; break;
        case 2: ans = 'C'; break;
        case 3: ans = 'G'; break;
        case 4: ans = 'T'; break;
        default: fprintf(stderr, "Invalid integer: %d in int2char\n", i); exit(EXIT_FAILURE);
    }
    return(ans);
}

int char2int(char c){
    int ans;
    switch(c){
        case 'A': ans = 1; break;
        case 'C': ans = 2; break;
        case 'G': ans = 3; break;
        case 'T': ans = 4; break;
        default: fprintf(stderr, "Invalid char: %c in char2int\n", c); exit(EXIT_FAILURE);
    }
    return(ans);
}

int count_occurrences(char *S, int n, int *SA, int *C, int **OCC, char *unit, int unitLen, int *covered, int MIN_number_repetitions )
{
    for(int i=0; i<n; i++) covered[i]=0;
    if(unitLen < LONG_UNIT_LEN_TH){
        // Do not rotate the unit
        int lb = 0;
        int ub = n-1;
        // A tandem unit must occur more than once.
        for(int i=0; i < (int)(MIN_number_repetitions * unitLen); i++ ){
            // parse from the end to the begin
            int x = char2int( unit[ unitLen - 1 - (i % unitLen) ] );
            if(0 < lb)
                lb = C[x] + OCC[x][lb-1];
            else
                lb = C[x];
            ub = C[x] + OCC[x][ub] - 1;
            if(lb > ub || lb < 0 || n-1 < ub)
                return(0);
        }
        for(int k=lb; k<=ub; k++){
            // Mark all letters in the unit occurrences.
            for(int j=0; j < (MIN_number_repetitions * unitLen) && SA[k]+j < n; j++)
                covered[ SA[k]+j ] = 1;
        }
    }else{
        // Divide a unit into shorter windows to escape from mutations and errors
        count_occurrences_long_unit(S, n, SA, C, OCC, unit, unitLen, covered, MIN_number_repetitions );
    }
    int total_cnt = 0;
    for(int i=0; i<n; i++){
        if(covered[i] == 1) total_cnt++;
    }
    return(total_cnt);
}

//#define PRINT_covered
void compute_sumOccurrences(char *S, int n, int *SA, int *C, int **OCC, int MIN_number_repetitions ){
    
    for(int j=0; j<unit_cnt; j++){
        int cnt = count_occurrences(S, n, SA, C, OCC, Units[j].string, Units[j].len, Units[j].covered, MIN_number_repetitions);
        Units[j].sumOccurrences += cnt;
        
        #ifdef PRINT_covered
        for(int i=0; i<Units[j].len; i++)    printf("%c", Units[j].string[i]);
        printf("\tlen = %d\tNum Occs = \t%d\n", Units[j].len, Units[j].sumOccurrences);
        #endif
    }
}

void coverage_by_units(char *S, int MIN_number_repetitions){
    struct timeval s, e;    gettimeofday(&s, NULL);
    
    int n;  for(n=0; S[n]!='\0'; n++){}
    //int Slen = n;     // One day to detect this bug ...
    n++;    // For building SA. The last letter of S is terminal symbol $
    
    unsigned char *Sint  = (unsigned char *)malloc(sizeof(unsigned char) * n);
    if(Sint == NULL)  exit(EXIT_FAILURE);
    
    // We encode $,A,C,G,T by 0,1,2,3,4 so that $=0<A=1<C=2<G=3<T=4
    for(int i=0; i<=(n-2); i++)
        Sint[i] = char2int(S[i]);
    Sint[n-1] = 0; // The last character must be $.
    
    int *SA = (int *)malloc(sizeof(int)*n);
    if(SA == NULL){ exit(EXIT_FAILURE); }
    int K = 4;
    SA_IS(Sint, SA, n, K, sizeof(unsigned char));
    
    int *BWT = (int *)malloc(sizeof(int) * n);
    if(BWT == NULL){ exit(EXIT_FAILURE); }
    
    int **OCC = malloc(sizeof(int *) * (K+1));
    if(OCC == NULL){ exit(EXIT_FAILURE); }
    
    int C0[5] = {0,0,0,0,0};
    for(int j=1; j <= K; j++){
        OCC[j] = malloc(sizeof(int) * n);
        if(OCC[j] == NULL){ exit(EXIT_FAILURE); }
        OCC[j][0] = 0;
    }
    for(int i=0; i<n; i++){
        BWT[i] = Sint[SA[i]-1];
        for(int j=1; j<= K; j++)
            if(0 < i){ OCC[j][i] = OCC[j][i-1]; }
        if(1 <= BWT[i] && BWT[i] <= 4){ // Skip updates when BWT[i] has $
            OCC[ BWT[i] ][i]++;
            C0[ BWT[i] ]++;
        }
    }
    int C[5]; C[0] = 0; C[1] = 1;
    for(int j=2; j<=K; j++){ C[j] = C[j-1] + C0[j-1]; }
    
    compute_sumOccurrences(S, n, SA, C, OCC, MIN_number_repetitions );
    
    free(Sint);
    free(SA);
    free(BWT);
    for(int j=1; j <= K; j++){ free(OCC[j]);}
    free(OCC);
    
    gettimeofday(&e, NULL);
    time_coverage_by_units += (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
}

#define UNDEFINED -1
int set_cover_greedy(FILE *ofp, char *S, int MIN_number_repetitions){
    // Solve the set cover problem in a greedy manner
    struct timeval s, e;    gettimeofday(&s, NULL);
    
    int n;  for(n=0; S[n]!='\0'; n++){} // n is the length of S.
    
    int *tmpCovered = (int *)malloc( sizeof(int) * n );
    if(tmpCovered == NULL){ exit(EXIT_FAILURE); }
    for(int i=0; i<n; i++) tmpCovered[i] = 0;
    int numKeyUnits = -1;
    int fixed = 0;
    int max_cum_cnt = 0; // max of cumulative count
    
    for(int j=0; j<unit_cnt; j++) Units[j].prio = UNDEFINED;
    
    for(int prio=0; prio < TOP_k_units; prio++){
        int max_unit = UNDEFINED;
        for(int j=0; j<unit_cnt; j++){
            if(Units[j].prio == UNDEFINED){
                int cum_cnt = 0;
                for(int i=0; i<n; i++)
                    if(tmpCovered[i] == 1 || Units[j].covered[i] == 1) cum_cnt++;
                // cumulative count of the j-th unit
                if(max_cum_cnt < cum_cnt){
                    max_cum_cnt = cum_cnt;
                    max_unit = j;
                }
            }
        }
        // max unit was computed.
        if(max_unit != UNDEFINED){
            Units[max_unit].prio = prio;
            Units[max_unit].sumOccurrences = max_cum_cnt;
            for(int i=0; i<n; i++)
                if(Units[max_unit].covered[i] == 1)
                    tmpCovered[i] = 1;
            /*
            #ifdef DEBUG_EDDC_mode 
                fprintf(ofp, "%d\t%s\n", max_cum_cnt,  Units[max_unit].string);
                for(int i=0; i<n; i++)  printf("%d", tmpCovered[i]);
                fprintf(ofp, "\n");
            #endif
             */
            if(MIN_unit_occupancy_ratio * n <= max_cum_cnt && fixed == 0){
                numKeyUnits = prio+1;  // i is the 0-origin indexing !
                fixed = 1;
            }
        }else
            break;
    }
    free(tmpCovered);
    gettimeofday(&e, NULL);
    time_set_cover_greedy += (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
    
    return(numKeyUnits);
}
