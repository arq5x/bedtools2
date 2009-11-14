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
	
	while (getline(fa,fastaInLine)) {
		
		if (fastaInLine.find(">",0) != 0 ) {
			currDNA += fastaDbLine;
		}
		else {
			if (currDNA.size() > 0) {

				vector<BED> bedList = bed->bedMapNoBin[currChrom];

				// loop through each BED entry for this chrom and 
				// mask the requested sequence in the FASTA file.
				for (unsigned int i = 0; i < bedList.size(); i++) {
					if (this->softMask) {
						//currDNA.replace(bedList[i].start, ((bedList[i].end - bedList[i].start)));
					}
					else {
						//currDNA.replace(bedList[i].start, ((bedList[i].end - bedList[i].start)));
					}
				}
				faOut << ">" << currChrom << endl << currDna << endl;
			}
			currChrom = fastaInLine.substr(1, fastaInLine.find_first_of(" ")-1);
			currDNA = "";
		}
	}
	
	// process the last chromosome.
	if (currDNA.size() > 0) {

		vector<BED> bedList = bed->bedMapNoBin[currChrom];

		// loop through each BED entry for this chrom and 
		// mask the requested sequence in the FASTA file.
		for (unsigned int i = 0; i < bedList.size(); i++) {
			if (this->softMask) {
				currDNA.replace(bedList[i].start, ((bedList[i].end - bedList[i].start)));
			}
			else {
				currDNA.replace(bedList[i].start, ((bedList[i].end - bedList[i].start)));
			}
		}
		faOut << ">" << currChrom << endl << currDna << endl;
	}
}

