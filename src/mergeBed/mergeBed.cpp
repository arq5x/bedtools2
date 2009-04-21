#include "mergeBed.h"

//::::::::::::::::::::::::::::::::::
// Constructor
//::::::::::::::::::::::::::::::::::

BedMerge::BedMerge(string &bedFile, bool &numEntries, int &maxDistance) {

	this->bedFile = bedFile;
	this->numEntries = numEntries;
	this->maxDistance = -1 * maxDistance;
	this->bed = new BedFile(bedFile);

}


//::::::::::::::::::::::::::::::::::
// Destructor
//::::::::::::::::::::::::::::::::::

BedMerge::~BedMerge(void) {
}



//::::::::::::::::::::::::::::::::::::::::::::::::::
// Merge overlapping BED entries into a single entry 
//::::::::::::::::::::::::::::::::::::::::::::::::::

void BedMerge::MergeBed() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bed->loadBedFileIntoMapNoBin();

	// loop through each chromosome and merge their BED entries
	for (masterBedMapNoBin::iterator m = bed->bedMapNoBin.begin(); m != bed->bedMapNoBin.end(); ++m) {

		// bedList is already sorted by start position.
		vector<BED> bedList = m->second; 

		int minStart = 999999999;
		int maxEnd = 0;
		bool OIP = false;       // OIP = Overlap In Progress.  Lame, I realize.
		int i;
		int mergeCount = 1;

		// loop through the BED entries for this chromosome
		// and look for overlaps
		for (i = 1; i < bedList.size(); ++i) {
			
			// Is there an overlap between the current and previous entries?		
			if ( overlaps(bedList[i-1].start, bedList[i-1].end, 
			bedList[i].start, bedList[i].end) >= this->maxDistance) {

				OIP = true;
				mergeCount++;
				minStart = min(bedList[i-1].start, min(minStart, bedList[i].start));
				maxEnd = max(bedList[i-1].end, max(maxEnd, bedList[i].end));

			}
			else if ( overlaps(minStart, maxEnd, 
			bedList[i].start, bedList[i].end) >= this->maxDistance) {

				mergeCount++;
				minStart = min(minStart, bedList[i].start);
				maxEnd = max(maxEnd, bedList[i].end);

			}
			else {

				// was there an overlap befor the current entry broke it?
				if (OIP) {
					if (this->numEntries) {
						cout << bedList[i-1].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << mergeCount << endl;
					}
					else {
						cout << bedList[i-1].chrom << "\t" << minStart << "\t" << maxEnd << endl;
					}
				}
				else {
					if (this->numEntries) {
						cout << bedList[i-1].chrom << "\t" << bedList[i-1].start << "\t" << bedList[i-1].end << "\t" << mergeCount << endl;
					}
					else {
						cout << bedList[i-1].chrom << "\t" << bedList[i-1].start << "\t" << bedList[i-1].end << endl;
					}
				}

				// reset things for the next overlapping "block"
				OIP = false;
				mergeCount = 1;
				minStart = 999999999;
				maxEnd = 0;

			}
		}

		// clean up based on the last entry for the current chromosome
		if (OIP) {
			if (this->numEntries) {
				cout << bedList[i-1].chrom << "\t" << minStart << "\t" << maxEnd << "\t" << mergeCount << endl;
			}
			else {
				cout << bedList[i-1].chrom << "\t" << minStart << "\t" << maxEnd << endl;	
			}
		}
		else {
			if (this->numEntries) {
				cout << bedList[i-1].chrom << "\t" << bedList[i-1].start << "\t" << bedList[i-1].end  << "\t" << mergeCount << endl;
			}
			else {
				cout << bedList[i-1].chrom << "\t" << bedList[i-1].start << "\t" << bedList[i-1].end << endl;	
			}
		}
	}
}

