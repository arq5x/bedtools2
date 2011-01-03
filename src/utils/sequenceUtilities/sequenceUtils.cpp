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
void reverseSequence(string &sequence) {
    std::reverse(sequence.begin(), sequence.end());
}

// Performs an in-place reverse complement conversion
void reverseComplement(string &sequence) {

    // reverse the sequence
    reverseSequence(sequence);

    // swap the bases
    for(unsigned int i = 0; i < sequence.length(); i++) {
        switch(sequence[i]) {
            case 'A':
                sequence[i] = 'T';
                break;
            case 'C':
                sequence[i] = 'G';
                break;
            case 'G':
                sequence[i] = 'C';
                break;
            case 'T':
                sequence[i] = 'A';
                break;
            case 'a':
                sequence[i] = 't';
                break;
            case 'c':
                sequence[i] = 'g';
                break;
            case 'g':
                sequence[i] = 'c';
                break;
            case 't':
                sequence[i] = 'a';
                break;
            default:
                break;
        }
    }
}


void toLowerCase(std::string &str)
{

    const int length = str.length();
    for(int i=0; i < length; ++i)
    {
        str[i] = std::tolower(str[i]);
    }

    // alternate, C++ style.
    //transform(str.start(), str.end(), str.start(), std::tolower);
}


void toUpperCase(std::string &str)
{

    const int length = str.length();
    for(int i=0; i < length; ++i)
    {
        str[i] = std::toupper(str[i]);
    }

    // alternate, C++ style.
    //transform(str.start(), str.end(), str.start(), std::toupper);
}
