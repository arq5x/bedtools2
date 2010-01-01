/*****************************************************************************
  coverageBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "coverageBed.h"


BedCoverage::BedCoverage(string &bedAFile, string &bedBFile, bool &forceStrand) {

	this->bedAFile = bedAFile;
	this->bedBFile = bedBFile;
	
	this->bedA = new BedFile(bedAFile);
	this->bedB = new BedFile(bedBFile);
	
	this->forceStrand = forceStrand;
}



BedCoverage::~BedCoverage(void) {
}


 
void BedCoverage::GetCoverage(istream &bedInput) {
	
	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bedB->loadBedFileIntoMap();

	string bedLine;                                                                                                                    
	int lineNum = 0;					// current input line number
	vector<string> bedFields;			// vector for a BED entry
	
	bedFields.reserve(12);	
		
	// process each entry in A
	while (getline(bedInput, bedLine)) {

		lineNum++;
		Tokenize(bedLine,bedFields);
		BED a;
		
		// find the overlaps with B if it's a valid BED entry. 
		if (bedA->parseLine(a, bedFields, lineNum)) {
	
			// increment the count of overlaps for each feature in B that
			// overlaps the current A interval
			bedB->countHits(bedB->bedMap[a.chrom], a, this->forceStrand);
		}	
		// reset for the next input line
		bedFields.clear();
	}
	
	// now, report the count of hist for each feature in B.
	for (masterBedMap::iterator c = bedB->bedMap.begin(); c != bedB->bedMap.end(); ++c) {
		map<int, vector<BED> > bin2Beds = c->second;

		for (map<int, vector<BED> >::iterator b = bin2Beds.begin(); b != bin2Beds.end(); ++b) {

			vector<BED> beds = b->second;
			for (unsigned int i = 0; i < beds.size(); i++) {
							
				int zeroDepthCount = 0;
				int depth = 0;
				
				int start = min(beds[i].minOverlapStart, beds[i].start);
				
				for (int pos = start+1; pos <= beds[i].end; pos++) {

					if (beds[i].depthMap.find(pos) != beds[i].depthMap.end()) {
						depth += beds[i].depthMap[pos].starts;
																				
						if ((depth == 0) && (pos > beds[i].start) && (pos <= beds[i].end)) {
							zeroDepthCount++;
						}
						
						depth = depth - beds[i].depthMap[pos].ends;
					}
					else {
						if ((depth == 0) && (pos > beds[i].start) && (pos <= beds[i].end)) {
							zeroDepthCount++;
						}
					}
				}

				// Report the coverage for the current interval.
				int length = beds[i].end - beds[i].start;

				bedB->reportBedTab(beds[i]);
				printf("%d\t%d\t%d\t%0.7f\n", beds[i].count, (length-zeroDepthCount), length, (float) (length-zeroDepthCount)/length);
			}
		}
	}
}


void BedCoverage::DetermineBedInput() {
	if (bedA->bedFile != "stdin") {   // process a file
		ifstream beds(bedA->bedFile.c_str(), ios::in);
		if ( !beds ) {
			cerr << "Error: The requested bed file (" << bedA->bedFile << ") could not be opened. Exiting!" << endl;
			exit (1);
		}
		GetCoverage(beds);
	}
	else {   						// process stdin
		GetCoverage(cin);		
	}
}


