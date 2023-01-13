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

void free_Units(){
    for(int i=0; i<MAX_NUMBER_UNITS; i++){
        if(Units[i].string != NULL) free(Units[i].string);
        if(Units[i].covered!= NULL) free(Units[i].covered);
    }
    if(Units != NULL) free(Units);
    #ifdef DEBUG_feed
    fprintf(stderr, "Succeeded in free (handle_one_file.c)\n");
    #endif
}

void malloc_Units(){
    unit_cnt = 0;
    Units = (Unit *) malloc( sizeof(Unit) * MAX_NUMBER_UNITS);
    if(Units == NULL) exit(EXIT_FAILURE);
    for(int i=0; i<MAX_NUMBER_UNITS; i++)
        Units[i].sumOccurrences = 0;
    for(int i=0; i<MAX_NUMBER_UNITS; i++){
        Units[i].string = malloc( sizeof(char) * MAX_UNIT_LENGTH );
        if(Units[i].string == NULL){ free_Units(); exit(EXIT_FAILURE); }
        Units[i].covered = malloc( sizeof(int) * MAX_READ_LENGTH );
        if(Units[i].covered == NULL){ free_Units(); exit(EXIT_FAILURE); }
    }
    Qreads = (QualifiedRead *) malloc( sizeof(QualifiedRead) * MAX_NUMBER_READS);
    if(Qreads == NULL) exit(EXIT_FAILURE);
    #ifdef DEBUG_feed
    fprintf(stderr, "Succeeded in malooc (handle_one_file.c)\n");
    #endif
}

void free_GlobalUnits(){
    free(nextReadID);
    for(int i=0; i<MAX_NUMBER_UNITS; i++){
        if(GlobalUnits[i].string != NULL) free(GlobalUnits[i].string);
    }
    if(GlobalUnits != NULL) free(GlobalUnits);
    #ifdef DEBUG_feed
    fprintf(stderr, "Succeeded in free (handle_one_file.c)\n");
    #endif
}

void malloc_GlobalUnits(){
    nextReadID = malloc( sizeof(char) * (MAX_ID_LENGTH+1) );
    global_unit_cnt = 0;
    GlobalUnits = (Unit *) malloc( sizeof(Unit) * MAX_NUMBER_UNITS);
    if(GlobalUnits == NULL) exit(EXIT_FAILURE);
    for(int i=0; i<MAX_NUMBER_UNITS; i++)
        GlobalUnits[i].sumOccurrences = 0;
    for(int i=0; i<MAX_NUMBER_UNITS; i++){
        GlobalUnits[i].string = malloc( sizeof(char) * MAX_UNIT_LENGTH );
        if(GlobalUnits[i].string == NULL){
            free_Units(); exit(EXIT_FAILURE);
        }
    }
    #ifdef DEBUG_feed
    fprintf(stderr, "Succeeded in malooc (handle_one_file.c)\n");
    #endif
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
    char *s = (char *)malloc(sizeof(char)*BLK);
    int i;
    char charCode;
    int cnt=0;
    int no_read = 1;
    
    while (fgets(s, BLK, fp) != NULL) { // Feed a string of size BLK from fp into string s
        no_read = 0;
        
        if(s[0] == '>'){
            if(read_cnt == -1){ // This is the first read
                read_cnt = 0;
                // Feed the ID of the current read into the ID of nextRead
                for(i=1; s[i]!='\0' && s[i]!='\n' && i<MAX_ID_LENGTH; i++)
                    nextReadID[i-1] = s[i];
                nextReadID[i-1] = '\0';
                // Move on to feed the current string
            }else{
                // Set the ID of currentRead to the ID of nextRead
                int j;
                for(j=0;  nextReadID[j] != '\0'; j++)
                    currentRead->ID[j] = nextReadID[j];
                currentRead->ID[j] = '\0';
                // Feed the ID of the current read into the ID of nextRead
                for(i=1; s[i]!='\0' && s[i]!='\n' && i<BLK; i++)
                    nextReadID[i-1] = s[i];
                nextReadID[i-1] = '\0';
                // Finalize the currentRead string by appending '\0'
                currentRead->string[cnt] = '\0';
                currentRead->len = cnt;
                read_cnt++;
                free(s);
                return;
            }
        }else{
            // Feed the string            
            for(i=0; s[i]!='\0' && s[i]!='\n' && s[i]!='\r'; i++){
                currentRead->string[cnt++] = capitalize(s[i]);
                if( MAX_READ_LENGTH <= cnt){
                    fprintf(stderr, "fatal error: The length %d is tentatively at most %i.\nread ID = %s\nSet MAX_READ_LENGTH to a larger value", cnt, MAX_READ_LENGTH, currentRead->ID);
                    free_Units();
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    free(s);
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
