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
void reverseComplement(string &seq, bool isRNA) {

    // reverse the sequence
    reverseSequence(seq);

    // swap the bases
    for(unsigned int i = 0; i < seq.length(); i++) {
        switch(seq[i]) {
            // nucleotides
            case 'A':
                if (!isRNA) seq[i] = 'T';
                else seq[i] = 'U';
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
            // issue 682
            case 'U':  
                seq[i] = 'A';
                break;
            case 'a':
                if (!isRNA) seq[i] = 't';
                else seq[i] = 'u';
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
            // issue 682
            case 'u':  
                seq[i] = 'a';
                break;

            // unknown (N)
            case 'N':
                seq[i] = 'N';
                break;
            case 'n':
                seq[i] = 'n';
                break;

            // unknown (X)
            case 'X':
                seq[i] = 'X';
                break;
            case 'x':
                seq[i] = 'x';
                break;

            // IUPAC
            case 'Y':
                seq[i] = 'R';
                break;
            case 'y':
                seq[i] = 'r';
                break;

            case 'R':
                seq[i] = 'Y';
                break;
            case 'r':
                seq[i] = 'y';
                break;

            case 'S':
                seq[i] = 'S';
                break;
            case 's':
                seq[i] = 's';
                break;

            case 'W':
                seq[i] = 'W';
                break;
            case 'w':
                seq[i] = 'w';
                break;

            case 'K':
                seq[i] = 'M';
                break;
            case 'k':
                seq[i] = 'm';
                break;

            case 'M':
                seq[i] = 'K';
                break;
            case 'm':
                seq[i] = 'k';
                break;

            case 'B':
                seq[i] = 'V';
                break;
            case 'b':
                seq[i] = 'v';
                break;

            case 'V':
                seq[i] = 'B';
                break;
            case 'v':
                seq[i] = 'b';
                break;

            case 'D':
                seq[i] = 'H';
                break;
            case 'd':
                seq[i] = 'h';
                break;

            case 'H':
                seq[i] = 'D';
                break;
            case 'h':
                seq[i] = 'd';
                break;

            default:
                break;
        }
    }
}


void toLowerCase(std::string &seq)
{
    const CHRPOS length = seq.length();
    for(CHRPOS i=0; i < length; ++i)
    {
        seq[i] = std::tolower(seq[i]);
    }
}


void toUpperCase(std::string &seq)
{
    const CHRPOS length = seq.length();
    for(CHRPOS i=0; i < length; ++i)
    {
        seq[i] = std::toupper(seq[i]);
    }
}


void getDnaContent(const string &seq, CHRPOS &a, CHRPOS &c, CHRPOS &g, CHRPOS &t, CHRPOS &n, CHRPOS &other)
{
    // swap the bases
    for(CHRPOS i = 0; i < (CHRPOS)seq.length(); i++) {
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
            case 'U':
            case 'u':
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


