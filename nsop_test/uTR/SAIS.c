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

unsigned char mask[]={0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01};
// S と L から構成される列を８文字毎に区切り、tyope S を 1, L を 0 でビット表現して 8ビット（1バイト）とし、配列 t に入れる. たとえば LSLSSLSS は 01011011 にコードする（主記憶を節約する本格的なプログラミング技術）
// 配列 t から i 番目の文字の type b を取り出すマクロ
// ビット演算　or演算は |, and 演算は &, ビットの反転は ~ である
#define tget(i) ( (t[(i)/8] & mask[(i)%8]) ? 1 : 0 )
// i 番目の文字の type が b であることを配列 t に入れるマクロ
#define tset(i, b) t[(i)/8] = (b) ? (mask[(i)%8] | t[(i)/8]) : ((~mask[(i)%8]) & t[(i)/8])
// 配列 s からのデータの取り出しを cs が sizeof(int) なら整数として s[i] を取り出し、そうでなければ unsigned char として取り出すマクロ
#define chr(i) ( cs==sizeof(int) ? ((int*)s)[i] : ((unsigned char *)s)[i] )
// i 番目の文字がleft-most S-type (LMS) ならば true、そうでなければ false を返すマクロ
#define isLMS(i) (i>0 && tget(i) && !tget(i-1))

// find the start or end of each bucket
void getBuckets(unsigned char *s, int *bkt, int n, int K, int cs, bool end) {
    int i, sum=0;
    // clear all buckets. K は $ を除いた文字の数
    for(i=0; i<=K; i++) bkt[i]=0;
    // compute the size of each bucket
    for(i=0; i<n; i++) bkt[chr(i)]++;
    for(i=0; i<=K; i++){ sum+=bkt[i]; bkt[i]=end ? sum : sum-bkt[i]; }
    // end が true ならば bucket の最後の位置を bkt[i] に入れ、false ならば最初の位置を入れる
}
// compute SAl　L-type の文字をソート
void induceSAl(unsigned char *t, int *SA, unsigned char *s, int *bkt, int n, int K, int cs, bool end) {
    int i, j;
    // find starts of buckets　配列 bkt に最初の位置を入れる
    getBuckets(s, bkt, n, K, cs, end);
    for(i=0; i<n; i++) {    // i は @ の置かれた位置を表現
        j=SA[i]-1;
        if(j>=0 && !tget(j)) SA[bkt[chr(j)]++]=j;
        // j 番目の文字が L-type (0) ならば、 ^ が置かれた位置　bkt[chr(j)] に j を入れて ^ を右へ1つ移動
    }
}
// compute SAs  S-type の文字をソート
void induceSAs(unsigned char *t, int *SA, unsigned char *s, int *bkt, int n, int K, int cs, bool end) {
    int i, j;
    // find ends of buckets　配列 bkt に最後の位置を入れる
    getBuckets(s, bkt, n, K, cs, end);
    for(i=n-1; i>=0; i--) {     // i は @ の置かれた位置を表現
        j=SA[i]-1;
        if(j>=0 && tget(j)) SA[--bkt[chr(j)]]=j;
        // j 番目の文字が S-type (1) ならば、 ^ が置かれた位置　bkt[chr(j)] を１つ左に移動して j を入れる
    }
}
// find the suffix array SA of s[0..n-1] in {1..K}n require s[n-1]=0 (the sentinel!), n>=2
// use a working space (excluding s and SA) of at most 2.25n+O(1) for a constant alphabet
void SA_IS(unsigned char *s, int *SA, int n, int K, int cs) {
    // LS-type array in bits　タイプ S もしくは L をビット表現して記載する配列を t とする
    unsigned char *t=(unsigned char *)malloc(n/8+1);
    int i, j;
    // classify the type of each character the sentinel must be in s1, important!!!
    // 各文字が S-type か L-type のどちらであるかを計算して配列 t に入れる
    tset(n-2, 0); tset(n-1, 1); //　最後から2番めの文字は 0 (L-type) で、最後の文字 $ は 1 (S-type)
    for(i=n-3; i>=0; i--)       // 後ろの文字から埋めてゆく
        tset(i, (chr(i)<chr(i+1) || (chr(i)==chr(i+1) && tget(i+1)==1)) ? 1 : 0);
    // reduce the problem by at least 1/2 sort all the S-substrings bucket array
    // 各 bucket の ^ の位置を格納する配列 bkt を定義
    int *bkt = (int *)malloc(sizeof(int)*(K+1));
    // find ends of buckets
    getBuckets(s, bkt, n, K, cs, true);     // true は各 bucket の最後の位置に ^ を置くことを意味する
    for(i=0; i<n; i++) SA[i]=-1;
    for(i=1; i<n; i++) if(isLMS(i)) SA[--bkt[chr(i)]]=i;
        // LMS ならば bucket の最後から詰める. LMS suffix の順序は未知なので等しいと考えて indude sort する.
    induceSAl(t, SA, s, bkt, n, K, cs, false);
        //  L-type の suffix をソートする. false は各 bucket の最初の位置に ^ を置くことを意味する.
    induceSAs(t, SA, s, bkt, n, K, cs, true);
        //  S-type の suffix をソートする. true は各 bucket の最後の位置に ^ を置くことを意味する.
    free(bkt);
    // compact all the sorted substrings into the first n1 items of SA
    // 2*n1 must be not larger than n (proveable)
    // LMS prefix の順番を保持しながら SA の前半部分に移動して、SAの後半部分を空ける
    // n1 は LMS prefix を入れる位置を表現しており、計算が終わると SA[n1] 以降は空けられる
    int n1=0;
    for(i=0; i<n; i++)
        if(isLMS(SA[i])) SA[n1++]=SA[i];

    // find the lexicographic names of substrings
    // init the name array buffer
    // SA の n1 以降の場所を、LMS prefix をコードした文字列を作成する場所として使う
    for(i=n1; i<n; i++) SA[i]=-1;   // SA の n1 以降の場所を -1 で初期化
    int name=0, prev=-1;            // name は LMS prefix に付与する順番
    for(i=0; i<n1; i++) {
        int pos=SA[i]; bool diff=false; // 1つ前の LMS prefix と一致しないとき diff を true とする
        for(int d=0; d<n; d++)      // 現在の LMS prefix が1つ前の LMS prefix と一致するか調べる
            if(prev==-1 || chr(pos+d)!=chr(prev+d) || tget(pos+d)!=tget(prev+d)){
                diff=true; break;   // 一致しなければ diff を true にして抜ける
            }else if(d>0 && (isLMS(pos+d) || isLMS(prev+d)))
                break;    // 次の LMS prefix まで到達し、一致したので diff = falseのまま抜ける
        if(diff) { name++; prev=pos; }  // LMS prefix が一致しなければ新しい順番を name に設定
        pos=(pos%2==0) ? pos/2 : (pos-1)/2; // LMS は少なくとも1つおきに出現するので場所をつめる
        SA[n1+pos]=name-1;  // 更新した name の1つ前の順番を付与する
    }
    for(i=n-1, j=n-1; i>=n1; i--)   // -1 を除きながら LMS prefix をコードした文字列を SA の後半に詰める
        if(SA[i]>=0) SA[j--]=SA[i];

    // solve the reduced problem
    // recurse if names are not yet unique
    // s1 は LMS prefix をコードした文字列で SA の後半を使うのに対して、その suffix array SA1 は SA の前半に計算することで主記憶を節約する(主記憶を節約する本格的なプログラミング技術)
    int *SA1=SA, *s1=SA+n-n1;
    // 同一順序の LMS prefix が存在するならば、LMS prefix の総数 n1 より name が小さくなるので SA_IS を再帰呼び出しして s1 の suffix array SA1 を決定
    if(name<n1)
        SA_IS((unsigned char*)s1, SA1, n1, name-1, sizeof(int));
    else // generate the suffix array of s1 directly　存在しないなら s1 そのものが suffix array
        for(i=0; i<n1; i++) SA1[s1[i]] = i;

    // induce the result for the original problem
    // LMS prefix をコード化した文字列の suffix array が計算できたので、コード化する前の文字列の位置にまず戻す
    bkt = (int *)malloc(sizeof(int)*(K+1));
    // put all the LMS characters into their buckets
    // find ends of buckets
    getBuckets(s, bkt, n, K, cs, true);
    for(i=1, j=0; i<n; i++)
        if(isLMS(i)) s1[j++]=i;
        // コード化する前の文字列での LMS prefix の位置 i を、使用済みの s1 の場所に暫定的に入れる
    // get index in s   コード化する前の文字列の位置で、LMS suffix の順番を SA1 に入れる
    for(i=0; i<n1; i++) SA1[i]=s1[SA1[i]];
    // init SA[n1..n-1]
    for(i=n1; i<n; i++) SA[i]=-1;
    for(i=n1-1; i>=0; i--) {
        j=SA[i]; SA[i]=-1;
        SA[--bkt[chr(j)]]=j;    // LMS suffix の順番を bucket に入れる
    }
    // LMS suffix の順番が bucket 内に入ったので L-type をソートして後に S-type をソートする
    induceSAl(t, SA, s, bkt, n, K, cs, false);
    induceSAs(t, SA, s, bkt, n, K, cs, true);
    free(bkt); free(t);
}
