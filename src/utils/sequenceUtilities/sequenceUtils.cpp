// 
//  sequenceUtils.cpp
//  BEDTools
//  
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//  Summary:  Contains common functions for finding BED overlaps.
//
//  Acknowledgment: I am grateful to Michael Stromberg for the code below to
//                  reverse complement a sequence.

#include "sequenceUtils.h"

// Performs an in-place sequence reversal
void ReverseSequence(string &sequence) {
	std::reverse(sequence.begin(), sequence.end());
}

// Performs an in-place reverse complement conversion
void ReverseComplement(string &sequence) {

	// reverse the sequence
	ReverseSequence(sequence);

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

