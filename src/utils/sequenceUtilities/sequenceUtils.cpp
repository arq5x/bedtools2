//
//  sequenceUtils.cpp
//  BEDTools
//
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//  Summary:  Contains common functions for manipulating DNA sequences.
//
//  Acknowledgment: I am grateful to Michael Stromberg for the code below to
//                  reverse complement a sequence.

#include "sequenceUtils.h"

// Performs an in-place sequence reversal
void reverseSequence(string &seq) {
    std::reverse(seq.begin(), seq.end());
}

// Performs an in-place reverse complement conversion
void reverseComplement(string &seq) {

    // reverse the sequence
    reverseSequence(seq);

    // swap the bases
    for(unsigned int i = 0; i < seq.length(); i++) {
        switch(seq[i]) {
            case 'A':
                seq[i] = 'T';
                break;
            case 'C':
                seq[i] = 'G';
                break;
            case 'G':
                seq[i] = 'C';
                break;
            case 'T':
                seq[i] = 'A';
                break;
            case 'a':
                seq[i] = 't';
                break;
            case 'c':
                seq[i] = 'g';
                break;
            case 'g':
                seq[i] = 'c';
                break;
            case 't':
                seq[i] = 'a';
                break;
            default:
                break;
        }
    }
}


void toLowerCase(std::string &seq)
{
    const int length = seq.length();
    for(int i=0; i < length; ++i)
    {
        seq[i] = std::tolower(seq[i]);
    }
}


void toUpperCase(std::string &seq)
{
    const int length = seq.length();
    for(int i=0; i < length; ++i)
    {
        seq[i] = std::toupper(seq[i]);
    }
}


void getDnaContent(const string &seq, int &a, int &c, int &g, int &t, int &n, int &other)
{
    // swap the bases
    for(unsigned int i = 0; i < seq.length(); i++) {
        switch(seq[i]) {
            case 'A':
            case 'a':
                a++;
                break;
            case 'C':
            case 'c':
                c++;
                break;
            case 'G':
            case 'g':
                g++;
                break;
            case 'T':
            case 't':
                t++;
                break;
            case 'N':
            case 'n':
                n++;
                break;
            default:
                other++;
                break;
        }
    }    
}


int countPattern(const string &seq, const string &pattern, bool ignoreCase)
{
    string s = seq;
    string p = pattern;
    // standardize the seq and the pattern 
    // if case should be ignored
    if (ignoreCase) {
        toUpperCase(s);
        toUpperCase(p);
    }
    int patternLength = p.size();
    int patternCount = 0;
    for(unsigned int i = 0; i < s.length(); i++) {
        if (s.substr(i,patternLength) == p) {
            patternCount++;
        }
    }
    return patternCount;
}


