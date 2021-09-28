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
#include <math.h>
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
            // does not occur more than once.
            if(lb > ub || lb < 0 || n-1 < ub) return(0);
        }
        // Mark all letters in the unit occurrences.
        for(int k=lb; k<=ub; k++){
            for(int j=0; j < (int)(MIN_number_repetitions * unitLen) && SA[k]+j < n; j++)
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

void compute_sumOccurrences(char *S, int n, int *SA, int *C, int **OCC, int MIN_number_repetitions ){
    
    for(int j=0; j<unit_cnt; j++){
        int cnt = count_occurrences(S, n, SA, C, OCC, Units[j].string, Units[j].len, Units[j].covered, MIN_number_repetitions);
        Units[j].sumOccurrences += cnt;
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


int cumulative_count(int n, int *tmpCovered, int unitIndex, int MIN_number_repetitions, int prio)
{
    // Enter chech_state=1 when the focal unit stops covering bases; namely,
    // case 1) we reach the end of reads
    // case 2) the focal base is suddenly covered by the other key units, but the previous is not.
    // case 3) The unit covers the focal base but does not match it and the number of mismatches in the current run, mismatch_run exceeds a threshold.
    // Otherwise, the unit covers the focal base, and we tries to extend the coverage.
    
    int run, check_state, mismatches, cum_cnt, mismatch_run, tandem_cnt;
    run = check_state = mismatches = cum_cnt = mismatch_run = tandem_cnt = 0;
    for(int i=0; i<=n; i++){
        check_state = 0;
        if(i == n)  // case 1: Terminate and process the remaining run
            check_state = 1;
        else if(tmpCovered[i] > 0){ // Already covered by other units
            cum_cnt++;
            if(0 < run) check_state = 1;    // case 2
            else check_state = 0;   // The base has been covered by the other units.
        }else if(Units[unitIndex].covered[i] == 0){
            // case 3: The unit covers the focal base but does not match it.
            run++;  mismatches++; mismatch_run++;
            // Check when the number of mistamtches, mismatch_run,
            // exceeds ceil(Units[unitIndex].len * MAX_DIS_RATIO)
            if(floor(run * MAX_DIS_RATIO) < mismatches ||
               ceil(Units[unitIndex].len * MAX_DIS_RATIO) <= mismatch_run){
                check_state = 1; mismatch_run = 0;
            }else check_state = 0;  // Continue the unit coverage extension.
        }else // The unit covers and matches the unit. Continue the unit coverage extension.
            run++;
        if(check_state == 1){
            // MIN_number_repetitions is either 2(MTR) or 1(MR), and 2(MTR) is examined first.
            if((Units[unitIndex].len * MIN_number_repetitions) <= run){
                // If the run of matches is qualified,
                // add the number of matches to the cumulative count, cum_cnt, and
                // memorize that matches are covered by the unit with priority "prio."
                cum_cnt += (run - mismatches);
                if(0 < prio){
                // The unit is an optimal one if prio is positive.
                // prio is set to 0 when we test whether the unit is optimal one or not.
                    for(int k=1; k <= run; k++)
                        if(Units[unitIndex].covered[i-k] == 1)
                            tmpCovered[i-k] = prio;
                }
                // Calculate the total length of tandem units, 2(MTR), and add the number of matches in the tandem units to tandem_cnt.
                if( (Units[unitIndex].len * 2) <= run)
                    tandem_cnt += (run - mismatches);
            }
            run = mismatches = 0;   // Reset the run.
        }
    }
    if(0 < prio){   // The unit is an optimal one if prio is positive.
        Units[unitIndex].sumOccurrences = cum_cnt;
        Units[unitIndex].sumTandem = tandem_cnt;
    }
    return(cum_cnt);
}


#define UNDEFINED -1
void set_cover_greedy(FILE *ofp, Read *currentRead, int MIN_number_repetitions){
    // Solve the set cover problem in a greedy manner
    struct timeval s, e;    gettimeofday(&s, NULL);
    
    char *S = currentRead->string;
    int n;  for(n=0; S[n]!='\0'; n++){} // n is the length of S.
    
    int *tmpCovered = (int *)malloc( sizeof(int) * (n+1) );
    if(tmpCovered == NULL){ exit(EXIT_FAILURE); }
    for(int i=0; i<n; i++) tmpCovered[i] = 0;
    int numKeyUnits = -1;
    int fixed = 0;
    int prev_cum_cnt = 0; // max of cumulative count
    int min_penalty = n;
    int sumTandem = 0;
    
    for(int j=0; j<unit_cnt; j++) Units[j].prio = UNDEFINED;
    int prio2unit[TOP_k_units+1];
    // 1-origin index. Filled from prio2unit[1]
    for(int prio=1; prio <= TOP_k_units; prio++){
        // Compute an optimal unit that minimizes the penalty.
        int min_penalty_unit = UNDEFINED;
        for(int j=0; j<unit_cnt; j++){
            if(Units[j].prio == UNDEFINED){
                // Setting the last arg to 0 untouches tmpCovered
                int tmp_cum_cnt = cumulative_count(n, tmpCovered, j, MIN_number_repetitions, 0);
                int tmp_penalty = (min_penalty - (n - prev_cum_cnt)) + (Units[j].len + (tmp_cum_cnt - prev_cum_cnt) / Units[j].len) + (n - tmp_cum_cnt);
                // cumulative count of the j-th unit
                if(tmp_penalty < min_penalty){
                    min_penalty = tmp_penalty;
                    prev_cum_cnt = tmp_cum_cnt;
                    min_penalty_unit = j; }
            }
        }
        // After computing an optimal unit that minimizes the penalty.
        if(min_penalty_unit != UNDEFINED){
            Units[min_penalty_unit].prio = prio;
            prio2unit[prio] = min_penalty_unit;
            // Setting the last arg to prio(>0) updates tmpCovered and sumOccurrences of the optimal unit
            int cum_min_penalty_unit = cumulative_count(n, tmpCovered, min_penalty_unit, MIN_number_repetitions, prio);
            sumTandem += Units[min_penalty_unit].sumTandem;
            // The last argument prio is set to 0 when we test whether the unit is optimal one or not.
            if(MIN_unit_occupancy_ratio * n <= cum_min_penalty_unit && fixed == 0){
                numKeyUnits = prio;  // i is the 1-origin indexing !
                fixed = 1;
                currentRead->sumTandem = sumTandem;
            }
        }else
            break;
    }
    sprintf(currentRead->RegExpression, "");
    int run = 0;
    int prevPrio = tmpCovered[0];
    for(int i=0; i<=n; i++){
        if(prevPrio != tmpCovered[i] || i == n){
            if(prevPrio == 0){
                for(int j=i-run; j<i; j++)
                    sprintf(currentRead->RegExpression, "%s%c", currentRead->RegExpression, S[j]);
            }else{
                if(Units[prio2unit[prevPrio]].len < run){
                    // Use a pair of square brackets for representing tandem repeat units
                    sprintf(currentRead->RegExpression, "%s[%d,%s,%d,%d]", currentRead->RegExpression, prevPrio, Units[prio2unit[prevPrio]].string, Units[prio2unit[prevPrio]].len, run);
                }else{
                    // Use a pair of parentheses for representing non-tandem repeat units
                    sprintf(currentRead->RegExpression, "%s(%d,", currentRead->RegExpression, prevPrio);
                    for(int j=i-run; j<i; j++)
                        sprintf(currentRead->RegExpression, "%s%c", currentRead->RegExpression, S[j]);
                    sprintf(currentRead->RegExpression, "%s,%d,%d)",  currentRead->RegExpression, Units[prio2unit[prevPrio]].len, run);
                }
            }
            if(i < n)
                sprintf(currentRead->RegExpression, "%s,", currentRead->RegExpression);
            run = 1;
            prevPrio = tmpCovered[i];
        }else
            run++;
    }
    // add the list of prios of all nucleotides
    sprintf(currentRead->RegExpression, "%s, P=", currentRead->RegExpression);
    for(int i=0; i<n; i++){
        int validPrio = tmpCovered[i];
        // Print all prios
        //if(numKeyUnits < tmpCovered[i]) validPrio = 0;
        sprintf(currentRead->RegExpression, "%s%d", currentRead->RegExpression, validPrio);
    }

    free(tmpCovered);
    gettimeofday(&e, NULL);
    time_set_cover_greedy += (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
    
    currentRead->numKeyUnits = numKeyUnits;
    //return(numKeyUnits);
}
