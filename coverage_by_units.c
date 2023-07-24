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
#include <string.h>
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
    if(Sint == NULL) { fprintf(stderr, "Failure to malloc Sint\n"); exit(EXIT_FAILURE); }
    
    // We encode $,A,C,G,T by 0,1,2,3,4 so that $=0<A=1<C=2<G=3<T=4
    for(int i=0; i<=(n-2); i++)
        Sint[i] = char2int(S[i]);
    Sint[n-1] = 0; // The last character must be $.
    
    int *SA = (int *)malloc(sizeof(int)*n);
    if(SA == NULL){ fprintf(stderr, "Failure to malloc SA\n"); exit(EXIT_FAILURE); }
    int K = 4;
    SA_IS(Sint, SA, n, K, sizeof(unsigned char));
    
    int *BWT = (int *)malloc(sizeof(int) * n);
    if(BWT == NULL){ fprintf(stderr, "Failure to malloc BWA\n"); exit(EXIT_FAILURE); }
    
    int **OCC = malloc(sizeof(int *) * (K+1));
    if(OCC == NULL){ fprintf(stderr, "Failure to malloc OCC\n"); exit(EXIT_FAILURE); }
    
    int C0[5] = {0,0,0,0,0};
    for(int j=1; j <= K; j++){
        OCC[j] = malloc(sizeof(int) * n);
        if(OCC[j] == NULL){ fprintf(stderr, "Failure to malloc OCC[j]\n"); exit(EXIT_FAILURE); }
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
                //if(0 <= prio){
                // The unit is an optimal one if prio is non-negative.
                // prio is set to -1 when we test whether the unit is optimal one or not.
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
    //if(0 <= prio){   // The unit is an optimal one if prio is non-negative.
        Units[unitIndex].sumOccurrences = cum_cnt;
        Units[unitIndex].sumTandem = tandem_cnt;
    }
    return(cum_cnt);
}

int partition(int* target, int* Pos, int left, int right) {
    // Generate a random number, genrand_int32(), using the the Mersenne Twister.
    // int random = left + genrand_int32() % (right - left + 1);
    int random = left + rand() % (right - left + 1);
    // Exchange target[right] and target[random].
    int tmp;
    tmp = target[right]; target[right] = target[random]; target[random] = tmp;
    tmp = Pos[right];    Pos[right] = Pos[random];       Pos[random] = tmp;
    int pivot = target[right];
    int i = left-1; // i scans the array from the left.
    int j = right;  // j scans the array from the right.
    for (;;) {
        // Move from the left until hitting a value no less than the pivot.
        for(i++; target[i] < pivot; i++){}
        // Move from the right until hitting a value no greater than the pivot.
        for(j--; pivot < target[j] && i < j; j--){}
        if (i >= j)  break;
        // Exchange target[i] and target[j].
        tmp = target[i];  target[i] = target[j];  target[j] = tmp;
        tmp = Pos[i];     Pos[i]    = Pos[j];     Pos[j] = tmp;
    }
    // Exchange target[i] and target[right].
    tmp = target[i];  target[i] = target[right];  target[right] = tmp;
    tmp = Pos[i];     Pos[i]    = Pos[right];     Pos[right] = tmp;
    return i;
}

void randomQuickSort3(int* target, int* Pos, int aLeft, int aRight) {
    int left = aLeft; int right = aRight;
    while (left < right) {
        int i = partition(target, Pos, left, right);
        if( i - left <= right - i ){ // The left interval is shorter.
            randomQuickSort3(target, Pos, left, i-1);
            left=i+1;
        }else{                       // The right interval is shorter.
            randomQuickSort3(target, Pos, i+1, right);
            right=i-1;
        }
    }
}

//#define DEBUG_batch_selection  // Print units selected
int unit_selection(Read *currentRead, int *prio2unit, int MIN_number_repetitions, int longer_unit_mode){
    
    for(int j=0; j<unit_cnt; j++) Units[j].prio = UNDEFINED;
    
    int string_len;     // string_len is the length of the string.
    for(string_len=0; currentRead->string[string_len]!='\0'; string_len++){}
    
    //MIN_number_repetitions = 1;
    int *tmpCovered = (int *) malloc(sizeof(int)*(string_len+1));
    for(int i=0; i<string_len+1; i++) tmpCovered[i] = 0;
    int *listPenalty = (int *) malloc(sizeof(int)*unit_cnt);
    int *listUnitInd = (int *) malloc(sizeof(int)*unit_cnt);
    for(int j=0; j<unit_cnt; j++){
        //Setting the last arg to 0 untouches tmpCovered
        int cum_cnt = cumulative_count(string_len, tmpCovered, j, MIN_number_repetitions, 0);
        int penalty = Units[j].len + (int)(cum_cnt / Units[j].len) + (string_len - cum_cnt);
        Units[j].penalty = penalty;
        listPenalty[j]   = penalty;
        listUnitInd[j]   = j;
    }
    #ifdef  DEBUG_batch_selection
    fprintf(stderr, "unit_cnt=%d\n", unit_cnt);
    #endif
    
    randomQuickSort3(listPenalty, listUnitInd, 0, unit_cnt-1); // Sort units according to tentative penalty
    
    int numKeyUnits;
    if(longer_unit_mode == 1){
        // Use top k units as they are
        numKeyUnits = MIN(TOP_k_units, unit_cnt);
        for(int p=0; p<numKeyUnits; p++)
            prio2unit[p] = listUnitInd[p];
        
        #ifdef  DEBUG_batch_selection
        fprintf(stderr, "\nCandidate of units (numKeyUnits=%d)\n", numKeyUnits);
        for(int p=0; p<numKeyUnits; p++)
            fprintf(stderr, "%d\t%s\n", listPenalty[p], Units[listUnitInd[p]].string);
        #endif
    }else{
        // greedy selection of units
        int k = MIN(TOP_k_units, unit_cnt);
        for(int i=0; i<string_len+1; i++) tmpCovered[i] = 0; // Initialize tmpCovered
        int prev_cum_cnt = 0; // max of cumulative count
        int min_penalty = string_len;   // Set min_penalty a large value
        numKeyUnits=0;
        for(int prio=0; prio<k; prio++){ // 0-origin
            // Compute an optimal unit that minimizes the penalty.
            int min_penalty_unit = UNDEFINED;
            for(int j=0; j<k; j++){
                int unitIndex = listUnitInd[j];
                if(Units[unitIndex].prio == UNDEFINED){
                    //Setting the last arg to 0 untouches tmpCovered
                    int tmp_cum_cnt = cumulative_count(string_len, tmpCovered, unitIndex, MIN_number_repetitions, 0);
                    int tmp_penalty = min_penalty + (Units[unitIndex].len + (tmp_cum_cnt - prev_cum_cnt) / Units[unitIndex].len) - (tmp_cum_cnt - prev_cum_cnt); // The second is the penalty of the tentatively added unit, and the last is the decrease of uncovered bases by the added unit
                    #ifdef  DEBUG_batch_selection
                    fprintf(stderr, "%s\tunitIndex=%d\ttmp_cum_cnt=%d\ttmp_penalty=%d\tmin_penalty=%d\n", Units[unitIndex].string, unitIndex, tmp_cum_cnt, tmp_penalty, min_penalty);
                    //for(int x=0; x<string_len; x++) fprintf(stderr, "%d", tmpCovered[x]);
                    //fprintf(stderr, "\n");
                    #endif
                    if(tmp_penalty <= min_penalty){ // Include two occurrences of a 2bp unit
                        min_penalty = tmp_penalty;
                        prev_cum_cnt = tmp_cum_cnt;
                        min_penalty_unit = unitIndex;
                    }
                }
            }
            
            // After computing an optimal unit that minimizes the penalty.
            if(min_penalty_unit == UNDEFINED) // No units with the min penalty
                break;
            else{
                #ifdef  DEBUG_batch_selection
                fprintf(stderr, "---min_penalty\t\t%d\t%s\n", min_penalty, Units[min_penalty_unit].string);
                #endif
                Units[min_penalty_unit].prio = prio;
                prio2unit[prio] = min_penalty_unit; // Setting the last arg to prio(>=0) updates tmpCovered and sumOccurrences of the optimal unit. The last argument prio is set to -1 when we test whether the unit is optimal one or not.
                cumulative_count(string_len, tmpCovered, min_penalty_unit, MIN_number_repetitions, prio+1);
                numKeyUnits = prio+1;  // i is the 0-origin indexing
            }
        }
    }
    free(tmpCovered);
    free(listPenalty);
    free(listUnitInd);
    currentRead->numKeyUnits = numKeyUnits;
    return(numKeyUnits);
}


void copy_Read(Read *a, Read *b){
    a->len = b->len;
    strcpy(a->string, b->string);
    for(int i=0; i < b->len; i++)
        a->intString[i] = b->intString[i];
    strcpy(a->ID, b->ID);
    a->numKeyUnits = b->numKeyUnits;
    a->mosaic_mode = b->mosaic_mode;
    a->sumTandem = b->sumTandem;
    a->discrepancy_ratio = b->discrepancy_ratio;
    a->mismatch_ratio  = b->mismatch_ratio;
    a->deletion_ratio  = b->deletion_ratio;
    a->insertion_ratio = b->insertion_ratio;
    strcpy(a->RegExpression, b->RegExpression);
    a->RegExpressionDecomp = b->RegExpressionDecomp;
    strcpy(a->decomposition, b->decomposition);
}

void set_cover_greedy_prep(Read *currentRead, int MIN_number_repetitions)
{
    clear_Units_incrementally();
    for(int j=0; j<unit_cnt; j++) Units[j].prio = UNDEFINED;
    // Set the set of units to the empty set
    get_non_self_overlapping_prefixes(currentRead->string);
    //for(int j=0; j<unit_cnt; j++) fprintf(stderr, "%s\n", Units[j].string);
    // Put non-self-overlapping units into Units by calling nsop and put_unit.
    coverage_by_units(currentRead->string, MIN_number_repetitions);
    // Set Units[j].sumOccurrences to the tentative number of occurrences of each candidate unit by calling compute_sumOccurrences, count_occurrences (SA for short motifs), and count_occurrences_long_unit (DP for long motifs).
}


//#define DEBUG_set_cover_greedy
void set_cover_greedy(Read *currentRead, int MIN_number_repetitions, int mode_longer_TRs, int mode_smooth){
    // Solve the set cover problem in a greedy manner
    struct timeval s, e;    gettimeofday(&s, NULL);
    
    int prio2unit[TOP_k_units+1];
    int numKeyUnits;
    //int mode_smooth = 1;
    
    Read *shortUnitRead = malloc(sizeof(Read));
    Read *longUnitRead = malloc(sizeof(Read));
        
    #ifdef DEBUG_set_cover_greedy
    fprintf(stderr, "\nID:  \t%s\n", currentRead->ID);
    #endif
    if(mode_longer_TRs == 1){
        // Mode of computing longer units
        copy_Read(longUnitRead, currentRead);
        set_cover_greedy_prep(longUnitRead, MIN_number_repetitions);
        numKeyUnits = unit_selection(longUnitRead, prio2unit, MIN_number_repetitions, 1); // 1=long read mode!!
        if(numKeyUnits == 0){
            longUnitRead->discrepancy_ratio = 1;
            #ifdef DEBUG_set_cover_greedy
            fprintf(stderr, "long:\tnumKeyUnits == 0\n");
            #endif
        }else{
            string_decomposer(longUnitRead, numKeyUnits, prio2unit, MIN_number_repetitions, mode_smooth, 1); // 1=long read mode!!
            #ifdef DEBUG_set_cover_greedy
            fprintf(stderr, "long:\t%s\t%f\t%f\t%f\t%f\n", longUnitRead->RegExpression, longUnitRead->discrepancy_ratio, longUnitRead->mismatch_ratio, longUnitRead->deletion_ratio, longUnitRead->insertion_ratio);
            #endif
        }
    }

    // Mode of computing shorter units
    copy_Read(shortUnitRead, currentRead);
    set_cover_greedy_prep(shortUnitRead, MIN_number_repetitions);
    numKeyUnits = unit_selection(shortUnitRead, prio2unit, MIN_number_repetitions, 0); // 0 = shorter read mode!!
    if(numKeyUnits == 0){
        shortUnitRead->discrepancy_ratio = 1;
        #ifdef DEBUG_set_cover_greedy
        fprintf(stderr, "short:\tnumKeyUnits == 0\n");
        #endif
        //free(shortUnitRead); free(longUnitRead); return;
    }else{
        string_decomposer(shortUnitRead, numKeyUnits, prio2unit, MIN_number_repetitions, mode_smooth, 0); // 0 = shorter read mode!!
        // Compute major key units that occupy MIN_COVERAGE
        int sumLen=0;
        int revised_numKeyUnits;
        for(int p=numKeyUnits-1; 0<=p; p--){
            sumLen += Units[prio2unit[p]].sumOccurrences;
            if( shortUnitRead->len * MIN_COVERAGE < sumLen){
                revised_numKeyUnits = numKeyUnits - p;
                break;
            }
        }
        // Revise prio2unit and redo the string decomposer
        for(int j=0, i=numKeyUnits-revised_numKeyUnits; j<revised_numKeyUnits; j++, i++)
            prio2unit[j] = prio2unit[i];
        numKeyUnits = revised_numKeyUnits;
        string_decomposer(shortUnitRead, numKeyUnits, prio2unit, MIN_number_repetitions, mode_smooth, 0);
        #ifdef DEBUG_set_cover_greedy
        fprintf(stderr, "short:\t%s\t%f\t%f\t%f\t%f\n", shortUnitRead->RegExpression, shortUnitRead->discrepancy_ratio, shortUnitRead->mismatch_ratio, shortUnitRead->deletion_ratio, shortUnitRead->insertion_ratio);
        #endif
    }

    if(mode_longer_TRs == 1){
        if(shortUnitRead->discrepancy_ratio <= longUnitRead->discrepancy_ratio){
            copy_Read(currentRead, shortUnitRead);
            #ifdef DEBUG_set_cover_greedy
            fprintf(stderr, "shorter was selected\n");
            #endif
        }else{
            copy_Read(currentRead, longUnitRead);
            #ifdef DEBUG_set_cover_greedy
            fprintf(stderr, "longer was selected\n");
            #endif
        }
    }else{
        copy_Read(currentRead, shortUnitRead);
    }
    
    free(shortUnitRead);
    free(longUnitRead);
    
    gettimeofday(&e, NULL);
    time_set_cover_greedy += (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
}
