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


void BedFilePE::binKeeperFind(map<int, vector<BED>, std::less<int> > &bk, const int start, const int end, vector<BED> &hits) {

	int startBin, endBin;
	startBin = (start >>_binFirstShift);
	endBin = ((end-1) >>_binFirstShift);
	
	for (int i = 0; i < 6; ++i) {
		int offset = binOffsetsExtended[i];

		for (int j = (startBin+offset); j <= (endBin+offset); ++j)  {
			for (vector<BED>::const_iterator el = bk[j].begin(); el != bk[j].end(); ++el) {
				
				if (overlaps(el->start, el->end, start, end) > 0) {
					hits.push_back(*el);
				}
				
			}
		}
		startBin >>= _binNextShift;
		endBin >>= _binNextShift;
	}
}


//***********************************************
// Common functions
//***********************************************

// Constructor
BedFilePE::BedFilePE(string &bedFile) {
	this->bedFile = bedFile;
}

// Destructor
BedFilePE::~BedFilePE(void) {
}


/*
	reportBedPETab
	
	Writes the _original_ BED entry for A.
	Works for BEDPE only.
*/
void BedFilePE::reportBedPETab(const BEDPE &a) {

	if (this->bedType == 6) {
		printf("%s\t%d\t%d\t%s\t%d\t%d\t", a.chrom1.c_str(), a.start1, a.end1,
		 									a.chrom2.c_str(), a.start2, a.end2);
	}
	else if (this->bedType == 7) {
		printf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t", a.chrom1.c_str(), a.start1, a.end1,
		 									a.chrom2.c_str(), a.start2, a.end2,
											a.name.c_str());
	}	
	else if (this->bedType == 8) {
		printf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t", a.chrom1.c_str(), a.start1, a.end1,
		 									a.chrom2.c_str(), a.start2, a.end2,
											a.name.c_str(), a.score.c_str());
	}
	else if (this->bedType == 10) {
		printf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t", a.chrom1.c_str(), a.start1, a.end1,
		 									a.chrom2.c_str(), a.start2, a.end2,
											a.name.c_str(), a.score.c_str(), a.strand1.c_str(), a.strand2.c_str());
	}
}



/*
	reportBedPENewLine
	
	Writes the _original_ BED entry for A.
	Works for BEDPE only.
*/
void BedFilePE::reportBedPENewLine(const BEDPE &a) {

	if (this->bedType == 6) {
		printf("%s\t%d\t%d\t%s\t%d\t%d\n", a.chrom1.c_str(), a.start1, a.end1,
		 									a.chrom2.c_str(), a.start2, a.end2);
	}
	else if (this->bedType == 7) {
		printf("%s\t%d\t%d\t%s\t%d\t%d\t%s\n", a.chrom1.c_str(), a.start1, a.end1,
		 									a.chrom2.c_str(), a.start2, a.end2,
											a.name.c_str());
	}	
	else if (this->bedType == 8) {
		printf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n", a.chrom1.c_str(), a.start1, a.end1,
		 									a.chrom2.c_str(), a.start2, a.end2,
											a.name.c_str(), a.score.c_str());
	}
	else if (this->bedType == 10) {
		printf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n", a.chrom1.c_str(), a.start1, a.end1,
		 									a.chrom2.c_str(), a.start2, a.end2,
											a.name.c_str(), a.score.c_str(), a.strand1.c_str(), a.strand2.c_str());
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
			bed.score = lineVector[7].c_str();
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
			bed.score = lineVector[7].c_str();

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
			bed.score = lineVector[7].c_str();
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
			bed.score = lineVector[7].c_str();

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
	else if ((lineVector.size() != this->bedType) && (lineVector.size() != 0)) {
		cerr << "Differing number of BEDPE fields encountered at line: " << lineNum << ".  Exiting..." << endl;
		exit(1);
	}
	else if ((lineVector.size() < 6) && (lineVector.size() != 0)) {
		cerr << "TAB delimited BEDPE file with at least 6 fields (chrom1, start1, end1, chrom2, start2, end2) is required at line: "<< lineNum << ".  Exiting..." << endl;
		exit(1);
	}
	return false;
}


void BedFilePE::loadBedPEFileIntoMap() {

	// open the BEDPE file for reading                                                                                                                                      
	ifstream bed(bedFile.c_str(), ios::in);
	if ( !bed ) {
		cerr << "Error: The requested bed file (" <<bedFile << ") could not be opened. Exiting!" << endl;
		exit (1);
	}

	string bedLine;
	BEDPE bedpeEntry;
	BED bedEntry1, bedEntry2;
                                                                                                                        
	unsigned int lineNum = 0;
	int bin1, bin2;
	vector<string> bedFields;	// vector of strings for each column in BED file.
	bedFields.reserve(10);		// reserve space for worst case (BED 6)

	while (getline(bed, bedLine)) {
		lineNum++;
		
		if (lineNum > 1) {

			Tokenize(bedLine,bedFields);

			if (parseBedPELine(bedpeEntry, bedFields, lineNum)) {
				
				// separate the BEDPE entry into separate
				// BED entries
				splitBedPEIntoBeds(bedpeEntry, lineNum, bedEntry1, bedEntry2);

				// load end1 into a UCSC bin map
				bin1 = getBin(bedEntry1.start, bedEntry1.end);
				this->bedMapEnd1[bedEntry1.chrom][bin1].push_back(bedEntry1);	

				// load end2 into a UCSC bin map
				bin2 = getBin(bedEntry2.start, bedEntry2.end);				
				this->bedMapEnd2[bedEntry2.chrom][bin2].push_back(bedEntry2);		
			}
			bedFields.clear();

		}
		else {
			if ((bedLine.find("track") != string::npos) || (bedLine.find("browser") != string::npos)) {
				continue;
			}
			else {
				Tokenize(bedLine,bedFields);

				if (parseBedPELine(bedpeEntry, bedFields, lineNum)) {
					// separate the BEDPE entry into separate
					// BED entries
					splitBedPEIntoBeds(bedpeEntry, lineNum, bedEntry1, bedEntry2);

					// load end1 into a UCSC bin map
					bin1 = getBin(bedEntry1.start, bedEntry1.end);
					this->bedMapEnd1[bedEntry1.chrom][bin1].push_back(bedEntry1);	

					// load end2 into a UCSC bin map
					bin2 = getBin(bedEntry2.start, bedEntry2.end);				
					this->bedMapEnd2[bedEntry2.chrom][bin2].push_back(bedEntry2);
				}
				bedFields.clear();			
			}
		}
	}
}


void BedFilePE::splitBedPEIntoBeds(const BEDPE &bedpeEntry, unsigned int lineNum, BED &bedEntry1, BED &bedEntry2) {
	
	/* 
	   Split the BEDPE entry into separate BED entries
	
	   NOTE: I am using a trick here where I store
	   the lineNum of the BEDPE from the original file
	   in the "count" column.  This allows me to later
	   resolve whether the hits found on both ends of BEDPE A
	   came from the same entry in BEDPE B.  Tracking by "name"
	   alone with fail when there are multiple mappings for a given
	   read-pair.
	*/
	
	bedEntry1.chrom = bedpeEntry.chrom1;
	bedEntry1.start = bedpeEntry.start1;
	bedEntry1.end = bedpeEntry.end1;
	bedEntry1.name = bedpeEntry.name;
	bedEntry1.score = bedpeEntry.score;
	bedEntry1.strand = bedpeEntry.strand1;
	bedEntry1.count = lineNum;
	bedEntry1.minOverlapStart = INT_MAX;
	
	bedEntry2.chrom = bedpeEntry.chrom2;
	bedEntry2.start = bedpeEntry.start2;
	bedEntry2.end = bedpeEntry.end2;
	bedEntry2.name = bedpeEntry.name;
	bedEntry2.score = bedpeEntry.score;
	bedEntry2.strand = bedpeEntry.strand2;		
	bedEntry2.count = lineNum;
	bedEntry2.minOverlapStart = INT_MAX;
}



