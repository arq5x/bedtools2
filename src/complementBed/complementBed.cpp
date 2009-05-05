// 
//  complementBed.cpp
//  BEDTools
//  
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//  Summary:  Return the intervals NOT spanned by a BED file.
//

#include "lineFileUtilities.h"
#include "complementBed.h"

//==========================
// Constructor
//
BedComplement::BedComplement(string &bedFile, string &genomeFile) {

	this->bedFile = bedFile;
	this->genomeFile = genomeFile;
	this->bed = new BedFile(bedFile);
}


//
// Destructor
//
BedComplement::~BedComplement(void) {
}



//
// Merge overlapping BED entries into a single entry 
//
void BedComplement::ComplementBed() {

	// open the GENOME file for reading
	ifstream genome(this->genomeFile.c_str(), ios::in);
	if ( !genome ) {
		cerr << "Error: The requested genome file (" <<this->genomeFile << ") could not be opened. Exiting!" << endl;
		exit (1);
	}

	string chrom;
	unsigned int size;

	map<string, int, less<string> > chromSizes; 

	while (genome >> chrom >> size) {
		chromSizes[chrom] = size;
	}

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bed->loadBedFileIntoMapNoBin();

	string currChrom;
	
	// loop through each chromosome and merge their BED entries
	for (masterBedMapNoBin::iterator m = bed->bedMapNoBin.begin(); m != bed->bedMapNoBin.end(); ++m) {

		currChrom = m->first;
		// bedList is already sorted by start position.
		vector<BED> bedList = m->second; 

		int minStart = INT_MAX;
		int maxEnd = 0;
		bool OIP = false;       // OIP = Overlap In Progress.  Lame, I realize.
		unsigned int i;
		int mergeCount = 1;
		bool firstSegment = true;
		// loop through the BED entries for this chromosome
		// and look for overlaps
		for (i = 1; i < bedList.size(); ++i) {
			
			// Are the current and previous entries within the maximum distance allowed by the user?
			// By default, maxDistance is 0, which allows for the following to be reported:
			//          OVERLAP                        BOOK-END
			// 		********                ***************
			//      	   *******		                   ***************
			// Ans. **************          ******************************
			if ( overlaps(bedList[i-1].start, bedList[i-1].end, 
						  bedList[i].start, bedList[i].end) 
						>= 0 ) {

				OIP = true;
				mergeCount++;
				minStart = min(bedList[i-1].start, min(minStart, bedList[i].start));
				maxEnd = max(bedList[i-1].end, max(maxEnd, bedList[i].end));

			}
			else if ( overlaps(minStart, maxEnd, 
							   bedList[i].start, bedList[i].end) 
							>= 0 ) {

				mergeCount++;
				minStart = min(minStart, bedList[i].start);
				maxEnd = max(maxEnd, bedList[i].end);

			}
			// either an overlap was broken or we have an 
			// interstice between two non-overlapping features
			else {

				// was there an overlap before the current entry broke it?
				// if so, report the complement as the position after the end of the
				// merged segment to the position before the start of the current segment.
				// For example:
				// *****************                  *******************
				//        ******************
				// Complement:              !!!!!!!!!!
				if (OIP) {
					if (firstSegment) {
						cout << currChrom << "\t" << 0 << "\t" << minStart << endl;
						firstSegment = false;
					}
					cout << currChrom << "\t" << maxEnd << "\t" << bedList[i].start << endl;
				}
				// otherwise report the interstice between the current and previous feature.
				// For example:
				// *****************                  *******************
				// Complement:      !!!!!!!!!!!!!!!!!!
				else {
					if (firstSegment) {
						cout << currChrom << "\t" << 0 << "\t" << bedList[i-1].start << endl;
						firstSegment = false;
					}
					cout << currChrom << "\t" << bedList[i-1].end << "\t" << bedList[i].start << endl;
				}

				// reset things for the next overlapping "block"
				OIP = false;
				mergeCount = 1;
				minStart = INT_MAX;
				maxEnd = 0;

			}
		}

		// clean up based on the last entry for the current chromosome
		// same as above code...
		if (OIP) {
			if (firstSegment) {
				cout << currChrom << "\t" << 0 << "\t" << minStart << endl;
				firstSegment = false;
			}
			cout << currChrom << "\t" << maxEnd << "\t" << chromSizes[currChrom] << endl;
		}
		else {
			if (firstSegment) {
				cout << currChrom << "\t" << 0 << "\t" << bedList[i-1].start-1 << endl;
				firstSegment = false;
			}
			cout << currChrom << "\t" << bedList[i-1].end << "\t" << chromSizes[currChrom] << endl;
		}
	}
}

