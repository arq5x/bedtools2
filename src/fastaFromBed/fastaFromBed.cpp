// 
//  fastaFromBed.cpp
//  BEDTools
//  
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//  Summary:  Extract fasta sequences based on BED intervals.
//
#include "lineFileUtilities.h"
#include "fastaFromBed.h"

Bed2Fa::Bed2Fa(bool &useName, string &dbFile, string &bedFile, string &fastaOutFile) {

	if (useName) {
		this->useName = true;
	}
	
	this->dbFile = dbFile;
	this->bedFile = bedFile;
	this->fastaOutFile = fastaOutFile;
	
	this->bed = new BedFile(this->bedFile);
	bed->loadBedFileIntoMapNoBin();
}


Bed2Fa::~Bed2Fa(void) {
}


//******************************************************************************
// ExtractDNA
//******************************************************************************
void Bed2Fa::ExtractDNA() {

	/* Make sure that we can oen all of the files successfully*/
	
	// open the fasta database for reading
	ifstream faDb(this->dbFile.c_str(), ios::in);
	if ( !faDb ) {
		cerr << "Error: The requested fasta database file (" << this->dbFile << ") could not be opened. Exiting!" << endl;
		exit (1);
	}
	
	// open the fasta database for reading
	ofstream faOut(this->fastaOutFile.c_str(), ios::out);
	if ( !faOut ) {
		cerr << "Error: The requested fasta output file (" << this->fastaOutFile << ") could not be opened. Exiting!" << endl;
		exit (1);
	}	
	

	/* Read the fastaDb chromosome by chromosome*/
	
	string fastaDbLine;
	string currChrom;
	string currDNA = "";
	currDNA.reserve(500000000);
	
	while (getline(faDb,fastaDbLine)) {
		
		if (fastaDbLine.find(">chr",0) != 0 ) {
			currDNA += fastaDbLine;
		}
		else {
			if (currDNA.size() > 0) {

				// retive the vecvtor of BED elements that come from the current 
				// chromosome that we are processing in the fasta "database".
						
				vector<BED> bedList = bed->bedMapNoBin[currChrom];

				// loop through each BED entry for this chrom and 
				// print the sequence
				for (unsigned int i = 0; i < bedList.size(); i++) {
						
					string dna = currDNA.substr(bedList[i].start, ((bedList[i].end - bedList[i].start)));
				
					if (!(this->useName)) {
				    	faOut << ">" << currChrom << ":"
					  		<< bedList[i].start << "-" << bedList[i].end
					  		<< "_" << bedList[i].name << endl
					  		<< dna << endl;
				  	}
				  	else {
				    	faOut << ">" << bedList[i].name << endl << dna << endl;
				  	}
				}
				currDNA = "";	
			}
			currChrom = fastaDbLine.substr(1, fastaDbLine.find_first_of(" ")-1);
		}
	}
	
	// process the last chromosome in the fasta file.
	if (currDNA.size() > 0) {

		// retrieve the vector of BED elements that come from the current 
		// chromosome that we are processing in the fasta "database".
		
		vector<BED> bedList = bed->bedMapNoBin[currChrom];

		// loop through each BED entry for this chrom and 
		// print the sequence.
		for (unsigned int i = 0; i < bedList.size(); i++) {
				
			string dna = currDNA.substr(bedList[i].start, ((bedList[i].end - bedList[i].start)));
		
			if (bedList[i].strand == "-")  {
				ReverseComplement(dna);
			}
			if (!(this->useName)) {
		    	faOut << ">" << currChrom << ":"
			  		<< bedList[i].start << "-" << bedList[i].end
			  		<< "_" << bedList[i].name << endl
			  		<< dna << endl;
		  	}
		  	else {
		    	faOut << ">" << bedList[i].name << endl << dna << endl;
		  	}
		}
		currDNA = "";	
	}
}

