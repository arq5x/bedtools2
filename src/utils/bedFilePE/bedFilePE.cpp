// 
//  bedFilePE.cpp
//  BEDTools
//  
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//  Summary:  Contains common functions for finding BED overlaps.
//
//  Acknowledgments: Much of the code herein is taken from Jim Kent's
//                   BED processing code.  I am grateful for his elegant
//					 genome binning algorithm and therefore use it extensively.


#include "bedFilePE.h"

//***********************************************
// Common functions
//***********************************************

// Constructor
BedFilePE::BedFilePE(string &bedFile) {
	this->bedFile = bedFile;
}

// Destructorc
BedFilePE::~BedFilePE(void) {
}

/*
	reportBedPE
	
	Writes the _original_ BED entry for A.
	Works for BEDPE only.
*/
void BedFilePE::reportBedPE(const BEDPE &a) {

	if (this->bedType == 6) {
		cout << a.chrom1 << "\t" << a.start1 << "\t" << a.end1 << "\t"
			 << a.chrom2 << "\t" << a.start2 << "\t" << a.end2;
	}
	else if (this->bedType == 7) {
		cout << a.chrom1 << "\t" << a.start1 << "\t" << a.end1 << "\t"
			 << a.chrom2 << "\t" << a.start2 << "\t" << a.end2 << "\t" 
			 << a.name;
	}	
	else if (this->bedType == 8) {
		cout << a.chrom1 << "\t" << a.start1 << "\t" << a.end1 << "\t"
			 << a.chrom2 << "\t" << a.start2 << "\t" << a.end2 << "\t" 
			 << a.name << "\t" << a.score;
	}
	else if (this->bedType == 10) {
		cout << a.chrom1 << "\t" << a.start1 << "\t" << a.end1 << "\t" 
			 << a.chrom2 << "\t" << a.start2 << "\t" << a.end2 << "\t"
			 << a.name << "\t" << a.score << "\t" << a.strand1 << "\t" << a.strand2;
	}
}

bool BedFilePE::parseBedPELine (BEDPE &bed, const vector<string> &lineVector, const int &lineNum) {

	if ((lineNum == 1) && (lineVector.size() >= 3)) {

		this->bedType = lineVector.size();

		if (this->bedType == 6) {
			bed.chrom1 = lineVector[0];
			bed.start1 = atoi(lineVector[1].c_str());
			bed.end1 = atoi(lineVector[2].c_str());
						
			bed.chrom2 = lineVector[3];
			bed.start2 = atoi(lineVector[4].c_str());
			bed.end2 = atoi(lineVector[5].c_str());

			return true;
		}
		else if (this->bedType == 7) {
			bed.chrom1 = lineVector[0];
			bed.start1 = atoi(lineVector[1].c_str());
			bed.end1 = atoi(lineVector[2].c_str());
						
			bed.chrom2 = lineVector[3];
			bed.start2 = atoi(lineVector[4].c_str());
			bed.end2 = atoi(lineVector[5].c_str());

			bed.name = lineVector[6];
			return true;
		}	
		else if (this->bedType == 8) {
			bed.chrom1 = lineVector[0];
			bed.start1 = atoi(lineVector[1].c_str());
			bed.end1 = atoi(lineVector[2].c_str());
						
			bed.chrom2 = lineVector[3];
			bed.start2 = atoi(lineVector[4].c_str());
			bed.end2 = atoi(lineVector[5].c_str());

			bed.name = lineVector[6];
			bed.score = atoi(lineVector[7].c_str());
			return true;
		}	
		else if (this->bedType == 10) {
			bed.chrom1 = lineVector[0];
			bed.start1 = atoi(lineVector[1].c_str());
			bed.end1 = atoi(lineVector[2].c_str());
						
			bed.chrom2 = lineVector[3];
			bed.start2 = atoi(lineVector[4].c_str());
			bed.end2 = atoi(lineVector[5].c_str());

			bed.name = lineVector[6];
			bed.score = atoi(lineVector[7].c_str());

			bed.strand1 = lineVector[8];
			bed.strand2 = lineVector[9];

			return true;
		}
		else {
			cerr << "Unexpected number of fields: " << lineNum << ".  Verify that your files are TAB-delimited and that your BEDPE file has 6,7,8 or 10 fields.  Exiting..." << endl;
			exit(1);
		}
		
		if (bed.start1 > bed.end1) {
			cerr << "Error: malformed BED entry at line " << lineNum << ". Start1 was greater than End1. Ignoring it and moving on." << endl;
			return false;
		}
		else if (bed.start2 > bed.end2) {
			cerr << "Error: malformed BED entry at line " << lineNum << ". Start2 was greater than End2. Ignoring it and moving on." << endl;
			return false;
		}
		else if ( (bed.start1 < 0) || (bed.end1 < 0) || (bed.start2 < 0) || (bed.end2 < 0) ) {
			cerr << "Error: malformed BED entry at line " << lineNum << ". Coordinate <= 0. Ignoring it and moving on." << endl;
			return false;
		}
	}
	else if ( (lineNum > 1) && (lineVector.size() == this->bedType)) {

		if (this->bedType == 6) {
			bed.chrom1 = lineVector[0];
			bed.start1 = atoi(lineVector[1].c_str());
			bed.end1 = atoi(lineVector[2].c_str());
						
			bed.chrom2 = lineVector[3];
			bed.start2 = atoi(lineVector[4].c_str());
			bed.end2 = atoi(lineVector[5].c_str());

			return true;
		}
		else if (this->bedType == 7) {
			bed.chrom1 = lineVector[0];
			bed.start1 = atoi(lineVector[1].c_str());
			bed.end1 = atoi(lineVector[2].c_str());
						
			bed.chrom2 = lineVector[3];
			bed.start2 = atoi(lineVector[4].c_str());
			bed.end2 = atoi(lineVector[5].c_str());

			bed.name = lineVector[6];
			return true;
		}	
		else if (this->bedType == 8) {
			bed.chrom1 = lineVector[0];
			bed.start1 = atoi(lineVector[1].c_str());
			bed.end1 = atoi(lineVector[2].c_str());
						
			bed.chrom2 = lineVector[3];
			bed.start2 = atoi(lineVector[4].c_str());
			bed.end2 = atoi(lineVector[5].c_str());

			bed.name = lineVector[6];
			bed.score = atoi(lineVector[7].c_str());
			return true;
		}	
		else if (this->bedType == 10) {
			bed.chrom1 = lineVector[0];
			bed.start1 = atoi(lineVector[1].c_str());
			bed.end1 = atoi(lineVector[2].c_str());
						
			bed.chrom2 = lineVector[3];
			bed.start2 = atoi(lineVector[4].c_str());
			bed.end2 = atoi(lineVector[5].c_str());

			bed.name = lineVector[6];
			bed.score = atoi(lineVector[7].c_str());

			bed.strand1 = lineVector[8];
			bed.strand2 = lineVector[9];

			return true;
		}
		else {
			cerr << "Unexpected number of fields: " << lineNum << ".  Verify that your files are TAB-delimited and that your BEDPE file has 6,7,8 or 10 fields.  Exiting..." << endl;
			exit(1);
		}

		if (bed.start1 > bed.end1) {
			cerr << "Error: malformed BED entry at line " << lineNum << ". Start1 was greater than End1. Ignoring it and moving on." << endl;
			return false;
		}
		else if (bed.start2 > bed.end2) {
			cerr << "Error: malformed BED entry at line " << lineNum << ". Start2 was greater than End2. Ignoring it and moving on." << endl;
			return false;
		}
		else if ( (bed.start1 < 0) || (bed.end1 < 0) || (bed.start2 < 0) || (bed.end2 < 0) ) {
			cerr << "Error: malformed BED entry at line " << lineNum << ". Coordinate <= 0. Ignoring it and moving on." << endl;
			return false;
		}
	}
	else if (lineVector.size() == 1) {
		cerr << "Only one BED field detected: " << lineNum << ".  Verify that your files are TAB-delimited.  Exiting..." << endl;
		exit(1);		
	}
	else if (lineVector.size() != this->bedType) {
		cerr << "Differing number of BED fields encountered at line: " << lineNum << ".  Exiting..." << endl;
		exit(1);
	}
	else if (lineVector.size() < 8) {
		cerr << "TAB delimited BED file with 10 fields (chrom1, start1, end1, strand1, chrom2, start2, end2, strand2) is required.  Exiting..." << endl;
		exit(1);
	}
	return false;
}


