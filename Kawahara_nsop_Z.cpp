// Riki Kawahara's algorithm that lists all non-self-overlapping substrings for a given string S in O(n^2)-time

//#define DUMP_Kawahara_nsop_Z

#ifndef DUMP_Kawahara_nsop_Z
#include "uTR.h"
#endif

#include <iostream>
#include <string>
#include <ctime>
using namespace std;

void dump_int_array(int *a, int len, string name){
    cout << name << "\t";
    for(int i=0; i<len; i++)    cout << a[i] << " ";
    cout << "\n";
}

// Kawahara's algorithm
void nsopSub(string S, int start){
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
    for(int j=1; j<S.size(); j++)
        put_unit(&(S.substr(0,j+1))[0]);
    #endif
    free(Z);
    free(OL);
}

void get_non_self_overlapping_prefixes(char *aS){
    std::string S(aS);
    for(int i=0; i<S.size(); i++)
        nsopSub( S.substr(i, S.size()-i), i );
}

#ifdef DUMP_Kawahara_nsop_Z
int main()
{
    //char *S = "AAAGAAAAG";
    //char *S = "momomosumomomomo";
    char S[] = "AGAGAGAGAGAGAGAGAGAGAGA";
    get_non_self_overlapping_prefixes(S);
    return(0);
}
#endif
