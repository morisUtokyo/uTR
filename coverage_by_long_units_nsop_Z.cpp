#include <iostream>
//#include <string>
#include <set>
#include <map>
#include <ctime>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <numeric>
#include <stack>
#include <vector>
#include <string>
#include <sys/time.h>
#include "uTR.h"
using namespace std;

class Alignment{
public:
	int start_x, start_y, end_x, end_y;
	float initial_score, score;
	string name;
	Alignment* predecessor;

	Alignment(int x0, int y0, int x1, int y1, float s, string n1){
		start_x = x0;
		start_y	= y0;
		end_x	= x1;
		end_y	= y1;
		initial_score	= s;
		score	= s;
		predecessor	= 0;
	}
	void print()const{
		//if(name.length() > 0) cout << name << "\t";
		cout << "(" <<	start_x << "," << start_y << ")\t-> (" << 
						end_x << ","  << end_y << ")\t score = " << score << endl;  
	}
    bool extend(char *S, int n, char *unit, int unitLen, int *covered, int MIN_number_repetitions ){
        int x, y;
        int numMismatches = 0;
        int maxMismatches = ceil(unitLen * MAX_DIS_RATIO);
        int numCopies = 0;
        
        for(x=end_x, y=end_y; x<n; ){
            if(S[x] == unit[ y % unitLen ]){
                x++, y++, score++;  // Half-open intervals
            }else{
                numMismatches++;
                if(numMismatches > maxMismatches) // If the number of mismatches exceed the threshold, stop extension.
                    break;
                else{x++, y++;}
            }
            // Update the ends if the length is unitLen or more.
            if(unitLen <= x - end_x){
                end_x = x; end_y = y;
                numMismatches = 0;
                numCopies++;
            }
        }
        // The unit must be duplicated.
        if(MIN_number_repetitions <= numCopies){
            // Use the last portion if the number of mismatches <= threshold
            if(numMismatches <= maxMismatches){
                end_x = x; end_y = y;
            }
            for(x=start_x, y=start_y; x<end_x; x++, y++) // Removal of the if statement fills mismatches in a partial mismatch of the unit, which is not implemented, but we show mismatches between the unit and read explicitly.
                if( S[x] == unit[ y % unitLen] )
                    covered[x] = 1;
            return(true);
        }else
            return(false);
    }
};

void count_occurrences_long_unit(char *S, int n, int *SA, int *C, int **OCC, char *unit, int unitLen, int *covered, int MIN_number_repetitions ){
    
    set<Alignment*> set_of_alignments;
    for(int l=0; l<unitLen; l++){
        int lb = 0;
        int ub = n-1;
        for(int i=0; i < WINDOW_LEN; i++ ){
            // parse from the end to the begin
            int x = char2int( unit[ (l+WINDOW_LEN-1-i)%unitLen ] );
            if(0 < lb)  lb = C[x] + OCC[x][lb-1];
            else        lb = C[x];
            ub = C[x] + OCC[x][ub] - 1;
            if(lb > ub || lb < 0 || n-1 < ub)   break;
        }
        for(int k=lb; k<=ub; k++)
            // Insert an alignment staring from (SA[k],l) of length WINDOW_LEN.
            set_of_alignments.insert( new Alignment( SA[k], l, SA[k] + WINDOW_LEN, l + WINDOW_LEN, WINDOW_LEN, "") );  // Half-open intervals are used
    }
    // Seed extension in a greedy manner
    multimap<int, Alignment*> sorted_by_X;
    for(set<Alignment*>::iterator A = set_of_alignments.begin();
        A != set_of_alignments.end();A++){
        sorted_by_X.insert(make_pair((*A)->start_x,*A));
    }
    int aln_start_x;    // The start position of a new alignment
    int aln_next_y = 0; // The y coord of the new alignment must be closest to this.
    // Scan sorted_by_X in the ascending order
    for(multimap<int, Alignment*>::iterator P = sorted_by_X.begin();
        P != sorted_by_X.end(); )
    {
        aln_start_x = P->second->start_x; // The start position of a new alignment
        // Find the seed alignment closest to aln_next_y
        auto seedP = P;  // Candidate of the seed
        int min_distance_P = (seedP->second->start_y - aln_next_y) % unitLen;
        auto nextP = P; nextP++;    // Iterator
        while(nextP != sorted_by_X.end() && nextP->second->start_x == aln_start_x){
            int distance_nextP = (nextP->second->start_y - aln_next_y) % unitLen;
            if(distance_nextP < min_distance_P){
                min_distance_P = distance_nextP;
                seedP = nextP;
            }
            nextP++;
        }
        // Extend the seed alignment
        if( seedP->second->extend(S, n, unit, unitLen, covered, MIN_number_repetitions ) ){
            aln_next_y = seedP->second->end_y;
            while( P != sorted_by_X.end() ){
                if( P->second->start_x <= seedP->second->end_x )
                    P++;
                else
                    break;
            }
        }else
            P++; // Move to the next
    }
    // Delete all alignments from set_of_alignments
    for(set<Alignment*>::iterator iter = set_of_alignments.begin(); iter != set_of_alignments.end(); iter++){
        delete *iter;
    }
}


//#define DUMP_Kawahara_nsop_Z
//-------------------------------------------------------------------
//
// Riki Kawahara has the copyright of the following program.
// Shinichi Morishita added an API on June 28, 2021.
//
//-------------------------------------------------------------------

void induced_sort(vector<int> &vec, int val_range, vector<int> &sa,
                  vector<bool> &sl, vector<int> &lms_idx) {
    vector<int> l(val_range, 0), r(val_range, 0);
    for (int c : vec) {
        if (c + 1 < val_range) ++l[c + 1];
        ++r[c];
    }
    partial_sum(l.begin(), l.end(), l.begin());
    partial_sum(r.begin(), r.end(), r.begin());

    fill(sa.begin(), sa.end(), -1);

    for (int i = lms_idx.size() - 1; i >= 0; --i) {
        sa[--r[vec[lms_idx[i]]]] = lms_idx[i];
    }

    for (int i : sa)
        if (i >= 1 && sl[i - 1]) {
            sa[l[vec[i - 1]]++] = i - 1;
        }

    fill(r.begin(), r.end(), 0);
    for (int c : vec) ++r[c];
    partial_sum(r.begin(), r.end(), r.begin());
    for (int k = sa.size() - 1, i = sa[k]; k >= 1; --k, i = sa[k])
        if (i >= 1 && !sl[i - 1]) {
            sa[--r[vec[i - 1]]] = i - 1;
        }
}

vector<int> sa_is(vector<int> &vec, int val_range) {
    const int n = vec.size();
    vector<int> sa(n), lms_idx;
    vector<bool> sl(n);

    sl[n - 1] = false;
    for (int i = n - 2; i >= 0; --i) {
        sl[i] = (vec[i] > vec[i + 1] || (vec[i] == vec[i + 1] && sl[i + 1]));
        if (sl[i] && !sl[i + 1]) lms_idx.push_back(i + 1);
    }
    reverse(lms_idx.begin(), lms_idx.end());

    induced_sort(vec, val_range, sa, sl, lms_idx);

    vector<int> new_lms_idx(lms_idx.size()), lms_vec(lms_idx.size());
    for (int i = 0, k = 0; i < n; ++i)
        if (!sl[sa[i]] && sa[i] >= 1 && sl[sa[i] - 1]) {
            new_lms_idx[k++] = sa[i];
        }

    int cur = 0;
    sa[n - 1] = cur;
    for (size_t k = 1; k < new_lms_idx.size(); ++k) {
        int i = new_lms_idx[k - 1], j = new_lms_idx[k];
        if (vec[i] != vec[j]) {
            sa[j] = ++cur;
            continue;
        }
        bool flag = false;
        for (int a = i + 1, b = j + 1;; ++a, ++b) {
            if (vec[a] != vec[b]) {
                flag = true;
                break;
            }
            if ((!sl[a] && sl[a - 1]) || (!sl[b] && sl[b - 1])) {
                flag = !((!sl[a] && sl[a - 1]) && (!sl[b] && sl[b - 1]));
                break;
            }
        }
        sa[j] = (flag ? ++cur : cur);
    }
    for (size_t i = 0; i < lms_idx.size(); ++i) {
        lms_vec[i] = sa[lms_idx[i]];
    }

    if (cur + 1 < (int)lms_idx.size()) {
        auto lms_sa = sa_is(lms_vec, cur + 1);

        for (size_t i = 0; i < lms_idx.size(); ++i) {
            new_lms_idx[i] = lms_idx[lms_sa[i]];
        }
    }

    induced_sort(vec, val_range, sa, sl, new_lms_idx);

    return sa;
}

vector<int> suffix_array(string s) {
    s += '$';
    vector<int> vec(s.size());
    for (int i = 0; i < (int)s.size(); ++i) vec[i] = s[i];
    auto sa = sa_is(vec, 128);
    sa.erase(sa.begin());
    return sa;
}

vector<int> lcp_array(string &s, vector<int> &sa) {
    int n = s.size();
    vector<int> rnk(n);
    for (int i = 0; i < n; ++i) rnk[sa[i]] = i;
    vector<int> lcp(n); //lcp(n - 1);
    int h = 0;
    for (int i = 0; i < n; ++i) {
        if (h > 0) --h;
        if (rnk[i] == 0) continue;
        for (int j = sa[rnk[i] - 1]; j + h < n && i + h < n; ++h) {
            if (s[i + h] != s[j + h]) break;
        }
        lcp[rnk[i] - 1] = h;
    }
    return lcp;
}

void dump_int_array(int *a, int len, string name){
    cout << name << "\t";
    for(int i=0; i<len; i++)    cout << a[i] << " ";
    cout << "\n";
}

// Riki Kawahara's algorithm that lists all non-self-overlapping substrings for a given string S in O(n^2)-time
void nsop(string S, int start){
        #ifdef DUMP_Kawahara_nsop_Z
        cout << "String\t";
        for(int i=0; i<S.size(); i++)
            cout << S[i] << " ";
        cout << "\t[" << start << "," << start+S.size()-1 << "]\n";
        cout << "i\t";
        for(int i=0; i<S.size(); i++)
            cout << start+i << " ";
        cout << "\n";
        #endif
    int *Z = (int*) malloc(sizeof(int) * S.size() );
    // Z-algorithm
    // Copied from https://qiita.com/Pro_ktmr/items/16904c9570aa0953bf05
    Z[0] = S.size();
    int i = 1, j = 0;
    while(i < S.size()){
        while(i + j < S.size() && S[j] == S[i + j]) j++;
        Z[i] = j;
        if(j == 0){ i++; continue; }
        int k = 1;
        while(k < j && k + Z[k] < j){
            Z[i + k] = Z[k]; k++;
        }
        i += k; j -= k;
    }
    // Kawahara's algorithm
    int *OL = (int*) malloc(sizeof(int) * (S.size()+1) );
    for(int i=0; i<S.size()+1; i++)
        OL[i]=0;
    for(int j=1; j<S.size(); j++){
        OL[j]++;
        OL[j+Z[j]]--;
    }
        #ifdef DUMP_Kawahara_nsop_Z
        dump_int_array(Z,  S.size(), "Z");
        dump_int_array(OL, S.size(), "OL");
        #endif
    for(int j=1; j<S.size(); j++){
        OL[j] += OL[j-1];
    }
        #ifdef DUMP_Kawahara_nsop_Z
        dump_int_array(OL, S.size(), "OL");
        for(int j=1; j<S.size(); j++){
            if(OL[j] == 0){  // non-self-overlapping
                //string tmp = S.substr(0,j+1);
                cout << "nsop\t" << &(S.substr(0,j+1))[0] << "\t[" << start << "," << start+j <<  "]\n";
            }
        }
        #endif
    #ifndef DUMP_Kawahara_nsop_Z
    for(int j=1; j<S.size(); j++){
        if(OL[j] == 0){
            put_unit(&(S.substr(0,j+1))[0]);
        }
    }
    #endif
    free(Z);
    free(OL);
}

void get_non_self_overlapping_prefixes(char *aS){
    std::string s(aS);
    int n = s.size();
    auto sa = suffix_array(s);
    auto lcp = lcp_array(s, sa);    // lcp of the i-th and (i+1)-th suffixes
    
    #ifdef DUMP_Kawahara_nsop_Z
    cout << "i\t";
    for(int i=0; i<sa.size(); i++)  cout << i << "\t";
    cout << "\nSA\t";
    for(int i=0; i<sa.size(); i++)  cout << sa[i] << "\t";
    cout << "\nlcp\t";
    for(int i=0; i<lcp.size(); i++) cout << lcp[i] << "\t";
    cout << "\n";
    #endif
    // Process repetitive non-self-overlapping substrings using lcp
        for(int i=0, max_lcp=0, max_i=0, up_state=1; i<sa.size(); i++){
            if(max_lcp < lcp[i]){
                max_lcp = lcp[i];
                max_i   = i;
                up_state = 1;
            }else if(max_lcp > lcp[i]){
                if(up_state == 1){
                    // Position max_i is locally maximum.
                    // Search non-self-overlapping substrings.
                    nsop( s.substr(sa[max_i], lcp[max_i]), sa[max_i] );
                    up_state = 0;
                }
                max_lcp = lcp[i];
                max_i   = i;
            }else{
                // max_lcp == lcp[i]
                // Untouch up_state, max_lcp, and max_i
        }
    }
}

#ifdef DUMP_Kawahara_nsop_Z
int main()
{
    //char *S = "AAAGAAAAG";
    //char *S = "momomosumomomomo";
    char S[] = "ATAATACGATAATAA";
    get_non_self_overlapping_prefixes(S);
    return(0);
}
#endif
