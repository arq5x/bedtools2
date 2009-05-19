// 
//  mergeBed.cpp
//  BEDTools
//  
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//  Summary:  Combines overlapping BED entries into a single entry.
//
#include "lineFileUtilities.h"
#include "mergeBed.h"

// ===============
// = Constructor =
// ===============
BedMerge::BedMerge(string &bedFile, bool &numEntries, int &maxDistance, bool &forceStrand) {

	this->bedFile = bedFile;
	this->numEntries = numEntries;
	this->maxDistance = -1 * maxDistance;
	this->forceStrand = forceStrand;
	
	this->bed = new BedFile(bedFile);

}


// =================
// =   Destructor  =
// =================
BedMerge::~BedMerge(void) {
}


// =====================================================
// = Merge overlapping BED entries into a single entry =
// =====================================================
void BedMerge::MergeBed() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bed->loadBedFileIntoMapNoBin();

	// loop through each chromosome and merge their BED entries
	for (masterBedMapNoBin::iterator m = bed->bedMapNoBin.begin(); m != bed->bedMapNoBin.end(); ++m) {

		// bedList is already sorted by start position.
		vector<BED> bedList = m->second; 

		int minStart = INT_MAX;
		int maxEnd = 0;
		bool OIP = false;       // OIP = Overlap In Progress.  Lame, I realize.
		unsigned int prev = 0;
		unsigned int curr = 0;
		int mergeCount = 1;


		// loop through the BED entries for this chromosome
		// and look for overlaps
		for (curr = 1; curr < bedList.size(); ++curr) {
			
			// Is there an overlap between the current and previous entries?		
			if ( overlaps(bedList[prev].start, bedList[prev].end, 
			 			bedList[curr].start, bedList[curr].end) >= this->maxDistance) {
				
				OIP = true;
				mergeCount++;
				minStart = min(bedList[prev].start, min(minStart, bedList[curr].start));
				maxEnd = max(bedList[prev].end, max(maxEnd, bedList[curr].end));

			}
			else if ( overlaps(minStart, maxEnd, 
							bedList[curr].start, bedList[curr].end) >= this->maxDistance) {
				mergeCount++;
				minStart = min(minStart, bedList[curr].start);
				maxEnd = max(maxEnd, bedList[curr].end);
			}
			else {

				// was there an overlap befor the current entry broke it?
				if (OIP) {
					if (this->numEntries) {
						cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << mergeCount << endl;
					}
					else {
						cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << endl;
					}
				}
				else {
					if (this->numEntries) {
						cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << mergeCount << endl;
					}
					else {
						cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << endl;
					}
				}

				// reset things for the next overlapping "block"
				OIP = false;
				mergeCount = 1;			
				minStart = INT_MAX;
				maxEnd = 0;

			}
			prev = curr;
		}

		// clean up based on the last entry for the current chromosome
		if (OIP) {
			if (this->numEntries) {
				cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << mergeCount << endl;
			}
			else {
				cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << endl;
			}
		}
		else {
			if (this->numEntries) {
				cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end  << "\t" << mergeCount << endl;
			}
			else {
				cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << endl;	
			}
		}
	}
}


// ==================================================================================
// = Merge overlapping BED entries into a single entry, accounting for strandedness =
// ==================================================================================
void BedMerge::MergeBedStranded() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bed->loadBedFileIntoMapNoBin();

	// loop through each chromosome and merge their BED entries
	for (masterBedMapNoBin::iterator m = bed->bedMapNoBin.begin(); m != bed->bedMapNoBin.end(); ++m) {

		// bedList is already sorted by start position.
		vector<BED> bedList = m->second; 

		// make a list of the two strands to merge separately.
		vector<string> strands(2);
		strands[0] = "+";
		strands[1] = "-";

		// do two passes, one for each strand.
		for (unsigned int s = 0; s < strands.size(); s++) {

			int minStart = INT_MAX;
			int maxEnd = 0;
			bool OIP = false;       // OIP = Overlap In Progress.  Lame, I realize.
			int prev = -1;
			unsigned int curr = 0;
			int mergeCount = 1;
			int numOnStrand = 0;
				
			// loop through the BED entries for this chromosome
			// and look for overlaps
			for (curr = 0; curr < bedList.size(); ++curr) {

				// if forcing strandedness, move on if the hit
				// is not on the current strand.
				
				if (bedList[curr].strand != strands[s]) {
					continue;		// continue force the next iteration of the for loop.
				}
				else {
					numOnStrand++;
				}

				// make sure prev points to an actual element on the
				// current strand
				if (prev < 0) {
					if (bedList[curr].strand == strands[s]) {
						prev = curr;
					}
					continue;
				}

				
				if ( overlaps(bedList[prev].start, bedList[prev].end, 
				 			bedList[curr].start, bedList[curr].end) >= this->maxDistance) {
					
					OIP = true;
					mergeCount++;
					minStart = min(bedList[prev].start, min(minStart, bedList[curr].start));
					maxEnd = max(bedList[prev].end, max(maxEnd, bedList[curr].end));

				}
				else if ( overlaps(minStart, maxEnd, 
								bedList[curr].start, bedList[curr].end) >= this->maxDistance) {

					mergeCount++;
					minStart = min(minStart, bedList[curr].start);
					maxEnd = max(maxEnd, bedList[curr].end);
				}
				else {

					// was there an overlap befor the current entry broke it?
					if (OIP) {
						if (this->numEntries) {
							cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << mergeCount << "\t" << strands[s] << endl;
						}
						else {
							cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << strands[s] << endl;
						}
					}
					else {
						if ((this->numEntries) && (numOnStrand > 0)) {
							cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << mergeCount << "\t" << strands[s] << endl;
						}
						else if (numOnStrand > 0) {
							cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << strands[s] << endl;
						}
					}

					// reset things for the next overlapping "block"
					OIP = false;
					mergeCount = 1;			
					minStart = INT_MAX;
					maxEnd = 0;
				}
				prev = curr;
			}

			// clean up based on the last entry for the current chromosome
			if (OIP) {
				if (this->numEntries) {
					cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << mergeCount << "\t" << strands[s] << endl;
				}
				else {
					cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << strands[s] << endl;
				}
			}
			else {
				if ((this->numEntries) && (numOnStrand > 0)) {
					cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << mergeCount << "\t" << strands[s] << endl;
				}
				else if (numOnStrand > 0) {
					cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << strands[s] << endl;
				}
			}
		}
	}
}
