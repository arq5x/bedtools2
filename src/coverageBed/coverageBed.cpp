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

		if (bedA->parseLine(a, bedFields, lineNum)) {	
			// count a as a hit with all the relevant features in B
			bedB->countHits(a, this->forceStrand);
		}	
		// reset for the next input line
		bedFields.clear();
	}
	
	//vector<int> depths;	// track the discrete depths for each base in B
						// used to calculate, min, max, median, etc.
						
	// now, report the count of hits for each feature in B.
	masterBedMap::const_iterator chromItr = bedB->bedMap.begin();
	masterBedMap::const_iterator chromEnd = bedB->bedMap.end();
	for (; chromItr != chromEnd; ++chromItr) {
	
		binsToBeds::const_iterator binItr = chromItr->second.begin();
		binsToBeds::const_iterator binEnd = chromItr->second.end();
		for (; binItr != binEnd; ++binItr) {

			vector<BED>::const_iterator bedItr = binItr->second.begin();
			vector<BED>::const_iterator bedEnd = binItr->second.end();
			for (; bedItr != bedEnd; ++bedItr) {
									
				int zeroDepthCount = 0;
				int depth = 0;
				int start = min(bedItr->minOverlapStart, bedItr->start);
				
				for (int pos = start+1; pos <= bedItr->end; pos++) {
					
					if (bedItr->depthMap.find(pos) != bedItr->depthMap.end()) {
						
						map<unsigned int, DEPTH> dMap = bedItr->depthMap;
						depth += dMap[pos].starts;
						//depths.push_back(depth);
																				
						if ((depth == 0) && (pos > bedItr->start) && (pos <= bedItr->end)) {
							zeroDepthCount++;
						}
						
						depth = depth - dMap[pos].ends;
					}
					else {
						if ((depth == 0) && (pos > bedItr->start) && (pos <= bedItr->end)) {
							zeroDepthCount++;
						}
						//depths.push_back(depth);
					}
				}

				// Report the coverage for the current interval.
				int length = bedItr->end - bedItr->start;
				int nonZeroBases =  (length-zeroDepthCount);
				float fractCovered = (float) nonZeroBases /length;
				//sort(depths.begin(), depths.end());
				
				bedB->reportBedTab(*bedItr);
				printf("%d\t%d\t%d\t%0.7f\n", bedItr->count, nonZeroBases, length, fractCovered);
				/*
				cout << bedItr->count << "\t";
				cout << (length-zeroDepthCount) << "\t";
				cout << length << "\t";
				cout << setw (10);
				cout << fractCovered << "\t";							
				cout << *min_element(depths.begin(), depths.end()) << "\t";
				cout << *max_element(depths.begin(), depths.end()) << "\t";
				cout << *(depths.begin()+depths.size()/2) << "\n";
				depths.clear();
				*/				
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


