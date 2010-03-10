/*****************************************************************************
  mergeBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "mergeBed.h"

// ===============
// = Constructor =
// ===============
BedMerge::BedMerge(string &bedFile, bool &numEntries, int &maxDistance, bool &forceStrand, bool &reportNames) {

	this->bedFile = bedFile;
	this->numEntries = numEntries;
	this->maxDistance = -1 * maxDistance;
	this->forceStrand = forceStrand;
	this->reportNames = reportNames;
	
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
		vector<string> names;

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
				
				//names.push_back(bedList[prev].name);
				names.push_back(bedList[curr].name);
			}
			else if ( overlaps(minStart, maxEnd, 
							bedList[curr].start, bedList[curr].end) >= this->maxDistance) {
				mergeCount++;
				minStart = min(minStart, bedList[curr].start);
				maxEnd = max(maxEnd, bedList[curr].end);
				
				names.push_back(bedList[curr].name);
			}
			else {

				// was there an overlap befor the current entry broke it?
				if (OIP) {
					if (this->numEntries) {
						cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << mergeCount << endl;
					}
					else if (this->reportNames) {
						cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t";
						for (unsigned int n = 0; n < names.size(); ++n) {
							if (n < (names.size() - 1)) {cout << names[n] << ";";}
							else {cout << names[n];}
						}
						cout << endl;
					}
					else {
						cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << endl;
					}
				}
				else {
					if (this->numEntries) {
						cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << mergeCount << endl;
					}
					else if (this->reportNames) {
						cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << bedList[prev].name << endl;
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
				names.clear();
				names.push_back(bedList[prev].name);
			}
			prev = curr;
		}

		// clean up based on the last entry for the current chromosome
		if (OIP) {
			if (this->numEntries) {
				cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << mergeCount << endl;
			}
			else if (this->reportNames) {
				cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t";
				for (unsigned int n = 0; n < names.size(); ++n) {
					if (n < (names.size() - 1)) {cout << names[n] << ";";}
					else {cout << names[n];}
				}
				cout << endl;
			}
			else {
				cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << endl;
			}
		}
		else {
			if (this->numEntries) {
				cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end  << "\t" << mergeCount << endl;
			}
			else if (this->reportNames) {
				cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << bedList[prev].name << endl;
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
			vector<string> names;	
			
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

					names.push_back(bedList[curr].name);
				}
				else if ( overlaps(minStart, maxEnd, 
								bedList[curr].start, bedList[curr].end) >= this->maxDistance) {

					mergeCount++;
					minStart = min(minStart, bedList[curr].start);
					maxEnd = max(maxEnd, bedList[curr].end);
					
					names.push_back(bedList[curr].name);
				}
				else {

					// was there an overlap before the current entry broke it?
					if (OIP) {
						if (this->numEntries) {
							cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << mergeCount << "\t" << strands[s] << endl;
						}
						else if (this->reportNames) {
							cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t";
							for (unsigned int n = 0; n < names.size(); ++n) {
								if (n < (names.size() - 1)) {cout << names[n] << ";";}
								else {cout << names[n];}
							}
							cout << "\t" << strands[s] << endl;
						}
						else {
							cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << strands[s] << endl;
						}
					}
					else {
						if ((this->numEntries) && (numOnStrand > 0)) {
							cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << mergeCount << "\t" << strands[s] << endl;
						}
						else if (this->reportNames) {
							cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << bedList[prev].name << "\t" << strands[s] << endl;
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
					names.clear();
					
					// add the name of the 
					names.push_back(bedList[curr].name);
				}
				prev = curr;
			}

			// clean up based on the last entry for the current chromosome
			if (OIP) {
				if (this->numEntries) {
					cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << mergeCount << "\t" << strands[s] << endl;
				}
				else if (this->reportNames) {
					cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t";
					for (unsigned int n = 0; n < names.size(); ++n) {
						if (n < (names.size() - 1)) {cout << names[n] << ";";}
						else {cout << names[n];}
					}
					cout << "\t" << strands[s] << endl;
				}
				else {
					cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << strands[s] << endl;
				}
			}
			else {
				if ((this->numEntries) && (numOnStrand > 0)) {
					cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << mergeCount << "\t" << strands[s] << endl;
				}
				else if ((this->reportNames) && (numOnStrand > 0)) {
					cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << bedList[prev].name << "\t" << strands[s] << endl;
				}
				else if (numOnStrand > 0) {
					cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << strands[s] << endl;
				}
			}
		}
	}
}
