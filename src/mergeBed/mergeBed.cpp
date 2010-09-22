/*****************************************************************************
  mergeBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "mergeBed.h"


void ReportMergedNames(const map<string, bool> &names) {
	unsigned int n = 0;
	map<string, bool>::const_iterator nameItr = names.begin();
	map<string, bool>::const_iterator nameEnd = names.end();
	for (; nameItr != nameEnd; ++nameItr) {
		if (n < (names.size() - 1)) {cout << nameItr->first << ";";}
		else {cout << nameItr->first;}
		n++;
	}
}

// ===============
// = Constructor =
// ===============
BedMerge::BedMerge(string &bedFile, bool &numEntries, int &maxDistance, bool &forceStrand, bool &reportNames) {

	_bedFile = bedFile;
	_numEntries = numEntries;
	_maxDistance = -1 * maxDistance;
	_forceStrand = forceStrand;
	_reportNames = reportNames;
	
	_bed = new BedFile(bedFile);
	
	if (_forceStrand == false)
		MergeBed();
	else
		MergeBedStranded();			
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
	_bed->loadBedFileIntoMapNoBin();

	// loop through each chromosome and merge their BED entries
	for (masterBedMapNoBin::iterator m = _bed->bedMapNoBin.begin(); m != _bed->bedMapNoBin.end(); ++m) {

		// bedList is already sorted by start position.
		vector<BED> bedList = m->second; 

		CHRPOS minStart = INT_MAX;
		CHRPOS maxEnd = 0;
		bool OIP = false;       // OIP = Overlap In Progress.  Lame, I realize.
		int prev = -1;
		unsigned int curr = 0;
		int mergeCount = 1;
		map<string, bool> names;

		// loop through the BED entries for this chromosome
		// and look for overlaps
		for (curr = 0; curr < bedList.size(); ++curr) {
			
			// make sure prev points to an actual element
			if (prev < 0) {
				prev = curr;
				continue;
			}

			// Is there an overlap between the current and previous entries?		
			if ( overlaps(bedList[prev].start, bedList[prev].end, 
			 			bedList[curr].start, bedList[curr].end) >= _maxDistance) {
				OIP = true;
				mergeCount++;
				minStart = min(bedList[prev].start, min(minStart, bedList[curr].start));
				maxEnd = max(bedList[prev].end, max(maxEnd, bedList[curr].end));

				names[bedList[prev].name] = true;
				names[bedList[curr].name] = true;
			}
			else if ( overlaps(minStart, maxEnd, 
							bedList[curr].start, bedList[curr].end) >= _maxDistance) {
				mergeCount++;
				minStart = min(minStart, bedList[curr].start);
				maxEnd = max(maxEnd, bedList[curr].end);
				names[bedList[curr].name] = true;
			}
			else {
				// was there an overlap befor the current entry broke it?
				if (OIP) {
					if (_numEntries) {
						cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << mergeCount << endl;
					}
					else if (_reportNames) {
						cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t";
						ReportMergedNames(names);
						cout << endl;
					}
					else {
						cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << endl;
					}
				}
				else {
					if (_numEntries) {
						cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << mergeCount << endl;
					}
					else if (_reportNames) {
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
				names[bedList[curr].name] = true;
			}
			prev = curr;
		}

		// clean up based on the last entry for the current chromosome
		if (OIP) {
			if (_numEntries) {
				cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << mergeCount << endl;
			}
			else if (_reportNames) {
				cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t";
				ReportMergedNames(names);
				cout << endl;
			}
			else {
				cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << endl;
			}
		}
		else {
			if (_numEntries) {
				cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end  << "\t" << mergeCount << endl;
			}
			else if (_reportNames) {
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
	_bed->loadBedFileIntoMapNoBin();

	// loop through each chromosome and merge their BED entries
	masterBedMapNoBin::const_iterator m    = _bed->bedMapNoBin.begin(); 
	masterBedMapNoBin::const_iterator mEnd = _bed->bedMapNoBin.end(); 
    for (; m != mEnd; ++m) {
		// bedList is already sorted by start position.
		vector<BED> bedList = m->second; 

		// make a list of the two strands to merge separately.
		vector<string> strands(2);
		strands[0] = "+";
		strands[1] = "-";

		// do two passes, one for each strand.
		for (unsigned int s = 0; s < strands.size(); s++) {

			CHRPOS minStart = INT_MAX;
			CHRPOS maxEnd = 0;
			bool OIP = false;       // OIP = Overlap In Progress.  Lame, I realize.
			int prev = -1;
			unsigned int curr = 0;
			int mergeCount = 1;
			int numOnStrand = 0;
			map<string, bool> names;	
			
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
				 			bedList[curr].start, bedList[curr].end) >= _maxDistance) {					
					OIP = true;
					mergeCount++;
					minStart = min(bedList[prev].start, min(minStart, bedList[curr].start));
					maxEnd = max(bedList[prev].end, max(maxEnd, bedList[curr].end));

					names[bedList[prev].name] = true;
					names[bedList[curr].name] = true;
				}
				else if ( overlaps(minStart, maxEnd, 
								bedList[curr].start, bedList[curr].end) >= _maxDistance) {
					mergeCount++;
					minStart = min(minStart, bedList[curr].start);
					maxEnd = max(maxEnd, bedList[curr].end);
					names[bedList[curr].name] = true;
				}
				else {

					// was there an overlap before the current entry broke it?
					if (OIP) {
						if (_numEntries) {
							cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << mergeCount << "\t" << strands[s] << endl;
						}
						else if (_reportNames) {
							cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t";
							ReportMergedNames(names);
							cout << "\t" << strands[s] << endl;
						}
						else {
							cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << strands[s] << endl;
						}
					}
					else {
						if ((_numEntries) && (numOnStrand > 0)) {
							cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << mergeCount << "\t" << strands[s] << endl;
						}
						else if (_reportNames) {
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
					
					// add the name of the current element in prep for the next block
					names[bedList[curr].name] = true;
				}
				prev = curr;
			}

			// clean up based on the last entry for the current chromosome
			if (OIP) {
				if (_numEntries) {
					cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << mergeCount << "\t" << strands[s] << endl;
				}
				else if (_reportNames) {
					cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t";
					ReportMergedNames(names);
					cout << "\t" << strands[s] << endl;
				}
				else {
					cout << bedList[prev].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << strands[s] << endl;
				}
			}
			else {
				if ((_numEntries) && (numOnStrand > 0)) {
					cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << mergeCount << "\t" << strands[s] << endl;
				}
				else if ((_reportNames) && (numOnStrand > 0)) {
					cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << bedList[prev].name << "\t" << strands[s] << endl;
				}
				else if (numOnStrand > 0) {
					cout << bedList[prev].chrom << "\t" << bedList[prev].start << "\t" << bedList[prev].end << "\t" << strands[s] << endl;
				}
			}
		}
	}
}
