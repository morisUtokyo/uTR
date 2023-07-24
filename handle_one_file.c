/*
 Copyright (c) 2021, Shinichi Morishita
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 The views and conclusions contained in the software and documentation are those
 of the authors and should not be interpreted as representing official policies,
 either expressed or implied, of the FreeBSD Project.
 */

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "uTR.h"

//#define DEBUG_feed

void clear_Units_incrementally(){
    for(int i=0; i<unit_cnt; i++){
        Units[i].sumOccurrences = 0;
        Units[i].prio = -1;
        // Leave other values untouched
    }
}

void malloc_Units(int num_reads){
//void malloc_Units(){
    unit_cnt = 0;
    
    Units = (Unit *) malloc( sizeof(Unit) * MAX_NUMBER_UNITS);
    if(Units == NULL) exit(EXIT_FAILURE);
    for(int i=0; i<MAX_NUMBER_UNITS; i++)
        Units[i].sumOccurrences = 0;
    
    keyUnits = (Unit *) malloc( sizeof(Unit) * TOP_k_units);
    if(keyUnits == NULL){
        fprintf(stderr, "Failure to malloc keyUnits(=%d)\n", TOP_k_units);
        exit(EXIT_FAILURE);
    }
    
    Qreads = (QualifiedRead *) malloc( sizeof(QualifiedRead) * num_reads);
    //Qreads = (QualifiedRead *) malloc( sizeof(QualifiedRead) * MAX_NUMBER_READS);
    if(Qreads == NULL){
        fprintf(stderr, "The number of reads in the input fasta file, %d, is too large.\n", num_reads);
        exit(EXIT_FAILURE);
    }
    #ifdef DEBUG_feed
    fprintf(stderr, "Succeeded in malooc (handle_one_file.c)\n");
    #endif
}
void free_Units(){
    if(Units != NULL) free(Units);
    if(keyUnits != NULL) free(keyUnits);
    if(Qreads != NULL) free(Qreads);
}

void malloc_GlobalVars(){
    WrapDP = malloc( sizeof(int) * WrapDPsize );
    if(WrapDP == NULL) exit(EXIT_FAILURE);
    
    nextReadID = malloc( sizeof(char) * (MAX_ID_LENGTH+1) );
    if(nextReadID == NULL) exit(EXIT_FAILURE);
    
    global_unit_cnt = 0;
    GlobalUnits =
        (Unit *) malloc( sizeof(Unit) * MAX_NUMBER_GLOBAL_UNITS);
    if(GlobalUnits == NULL) exit(EXIT_FAILURE);
    for(int i=0; i<MAX_NUMBER_GLOBAL_UNITS; i++)
        GlobalUnits[i].sumOccurrences = 0;
}
void free_GlobalVars(){
    free(WrapDP);
    free(nextReadID);
    if(GlobalUnits != NULL) free(GlobalUnits);
}


char capitalize(char c){
    char charCode;
    switch(c){
        case 'A':
        case 'a':
            charCode = 'A'; break;
        case 'C':
        case 'c':
            charCode = 'C'; break;
        case 'G':
        case 'g':
            charCode = 'G'; break;
        case 'T':
        case 't':
            charCode = 'T'; break;
        default:
            fprintf(stderr, "Invalid character: %c in capitalize\n", c); exit(EXIT_FAILURE);
    }
    return(charCode);
}

int count_reads(char *inputFile){
    FILE *fp = fopen(inputFile, "r");
    if(fp == NULL){
        fprintf(stderr, "fatal error: cannot open %s\n", inputFile);
        fflush(stderr);
        exit(EXIT_FAILURE);
    }
    char s[BLK+1];
    int cnt = 0;
    while (fgets(s, BLK, fp) != NULL)
        if(s[0] == '>') cnt++;
    fclose(fp);
    return(cnt);
}

FILE* init_handle_one_file(char *inputFile){
    FILE *fp = fopen(inputFile, "r");
    if(fp == NULL){
        fprintf(stderr, "fatal error: cannot open %s\n", inputFile);
        fflush(stderr);
        exit(EXIT_FAILURE);
    }
    read_cnt = -1;
    return(fp);
}

void return_one_read(FILE *fp, Read *currentRead){
    char s[BLK+1];
    int i;
    char charCode;
    int cnt=0;
    int no_read = 1;
    
    while (fgets(s, BLK, fp) != NULL) { // Feed a string of size BLK from fp into string s
        no_read = 0;
        fflush(fp);
        
        if(s[0] == '>'){
            if(read_cnt != -1){ // This is NOT the first read
                // Set the ID of currentRead to the ID of nextRead
                int j;
                for(j=0;  nextReadID[j] != '\0'; j++)
                    currentRead->ID[j] = nextReadID[j];
                currentRead->ID[j] = '\0';
            }
            // Feed the ID of the current read into the ID of nextRead
            int shift;
            if(s[1]==' ') shift=2; else shift=1;  // Skip the space at the head
            for(i=0; s[shift+i]!='\0' && s[shift+i]!='\n' && s[shift+i]!='\r' && i<BLK; i++)
                nextReadID[i] = s[shift+i];
            nextReadID[i] = '\0';
            
            if(read_cnt == -1){ // This is the first read
                read_cnt = 0;
            }else{
                read_cnt++;
                // Finalize the currentRead string by appending '\0'
                currentRead->string[cnt] = '\0';
                currentRead->len = cnt;
                return;
            }
        }else{
            // Feed the string            
            for(i=0; s[i]!='\0' && s[i]!='\n' && s[i]!='\r' && i<BLK; i++){
                currentRead->string[cnt] = capitalize(s[i]);
                currentRead->intString[cnt] = char2int(capitalize(s[i]));
                cnt++;
                if( MAX_READ_LENGTH <= cnt){
                    fprintf(stderr, "fatal error: The length %d is tentatively at most %i.\nread ID = %s\nSet MAX_READ_LENGTH to a larger value", cnt, MAX_READ_LENGTH, currentRead->ID);
                    free_Units();
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    if(no_read == 1){  // No reads
        currentRead->len = 0;
    }else{
        // Process the last read.
        // Set the ID of currentRead to the ID of nextRead
        int j;
        for(j=0;  nextReadID[j] != '\0'; j++)
            currentRead->ID[j] = nextReadID[j];
        currentRead->ID[j] = '\0';
        // Finalize the currentRead string by appending '\0'
        currentRead->string[cnt] = '\0';
        currentRead->len = cnt;
        read_cnt++;
    }
}
