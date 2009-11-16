// 
//  maskFastaFromBed.cpp
//  BEDTools
//  
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//  Summary:  Mask fasta sequences based on BED intervals.
//
#include "lineFileUtilities.h"
#include "maskFastaFromBed.h"


MaskFastaFromBed::MaskFastaFromBed(string &fastaInFile, string &bedFile, string &fastaOutFile, bool &softMask) {

	this->softMask = false;
	if (softMask) {
		this->softMask = true;
	}
	
	this->fastaInFile = fastaInFile;
	this->bedFile = bedFile;
	this->fastaOutFile = fastaOutFile;
	
	this->bed = new BedFile(this->bedFile);
	bed->loadBedFileIntoMapNoBin();
}


MaskFastaFromBed::~MaskFastaFromBed(void) {
}


//******************************************************************************
// Mask the Fasta file based on the coordinates in the BED file.
//******************************************************************************
void MaskFastaFromBed::MaskFasta() {

	/* Make sure that we can open all of the files successfully*/
	
	// open the fasta database for reading
	ifstream fa(this->fastaInFile.c_str(), ios::in);
	if ( !fa ) {
		cerr << "Error: The requested fasta file (" << this->fastaInFile << ") could not be opened. Exiting!" << endl;
		exit (1);
	}
	
	// open the fasta database for reading
	ofstream faOut(this->fastaOutFile.c_str(), ios::out);
	if ( !faOut ) {
		cerr << "Error: The requested fasta output file (" << this->fastaOutFile << ") could not be opened. Exiting!" << endl;
		exit (1);
	}	
	

	/* Read the fastaDb chromosome by chromosome*/
	
	string fastaInLine;
	string currChrom;
	string currDNA = "";
	currDNA.reserve(500000000);
	int fastaWidth = -1;
	bool widthSet  = false;
	int start, end, length;
	string replacement;
	
	while (getline(fa,fastaInLine)) {
		
		if (fastaInLine.find(">",0) != 0 ) {
			if (widthSet == false) {
				fastaWidth = fastaInLine.size();
				widthSet = true;
			}
			currDNA += fastaInLine;
		}
		else {
			if (currDNA.size() > 0) {

				vector<BED> bedList = bed->bedMapNoBin[currChrom];

				/*
					loop through each BED entry for this chrom and 
					mask the requested sequence in the FASTA file.
				*/
				for (unsigned int i = 0; i < bedList.size(); i++) {
					start = bedList[i].start;
					end = bedList[i].end;
					length = end - start;
					
					/* 
					   (1) if soft masking, extract the sequence, lowercase it,
					       then put it back
					    (2) otherwise replace with Ns
					*/
					if (this->softMask) {
						replacement = currDNA.substr(start, length);
						toLowerCase(replacement);
						currDNA.replace(start, length, replacement);
					}
					else {
						string hardmask(length, 'N');
						currDNA.replace(start, length, hardmask);
					}
				}
				// write the masked chrom to the output file
				PrettyPrintChrom(faOut, currChrom, currDNA, fastaWidth);
			}
			
			// reset for the next chromosome.
			currChrom = fastaInLine.substr(1, fastaInLine.find_first_of(" ")-1);
			currDNA = "";
		}
	}
	
	// process the last chromosome.
	// exact same logic as in the main loop.
	if (currDNA.size() > 0) {

		vector<BED> bedList = bed->bedMapNoBin[currChrom];

		for (unsigned int i = 0; i < bedList.size(); i++) {
			start = bedList[i].start;
			end = bedList[i].end;
			length = end - start;
			
			if (this->softMask) {
				replacement = currDNA.substr(start, length);
				toLowerCase(replacement);
				currDNA.replace(start, length, replacement);
			}
			else {
				string hardmask(length, 'N');
				currDNA.replace(start, length, hardmask);
			}
		}
		PrettyPrintChrom(faOut, currChrom, currDNA, fastaWidth);
	}
	
	fa.close();
	faOut.close();
}

void MaskFastaFromBed::PrettyPrintChrom(ofstream &out, string chrom, const string &sequence, int width) {
	
	int seqLength = sequence.size();
	 
	out << ">" << chrom << endl;
	for(int i=0; i < seqLength; i+=width)
	{
		if (i+width < seqLength) out << sequence.substr(i, width) << endl;
		else out << sequence.substr(i, seqLength-i) << endl;
	}
}

