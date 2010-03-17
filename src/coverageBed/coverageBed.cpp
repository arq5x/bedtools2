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


BedCoverage::BedCoverage(string &bedAFile, string &bedBFile, bool &forceStrand, bool &writeHistogram) {
	
	this->bedAFile       = bedAFile;
	this->bedBFile       = bedBFile;
	
	this->bedA           = new BedFile(bedAFile);
	this->bedB           = new BedFile(bedBFile);
	
	this->forceStrand    = forceStrand;
	this->writeHistogram = writeHistogram;
}



BedCoverage::~BedCoverage(void) {
}


 
void BedCoverage::CollectCoverage(istream &bedInput) {
	
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
	
	// report the coverage (summary or histogram) for BED B.
	ReportCoverage();					
}


void BedCoverage::ReportCoverage() {

	map<unsigned int, unsigned int> allDepthHist;
	unsigned int totalLength = 0;

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
				
				// track the numnber of bases in the feature covered by
				// 0, 1, 2, ... n features in A
				map<unsigned int, unsigned int> depthHist;
				map<unsigned int, DEPTH>::const_iterator depthItr;
				
				for (int pos = start+1; pos <= bedItr->end; pos++) {
					
					depthItr = bedItr->depthMap.find(pos);
					
					if (depthItr != bedItr->depthMap.end()) {
						depth += depthItr->second.starts;
						if ((pos > bedItr->start) && (pos <= bedItr->end)) {	
							if (depth == 0) zeroDepthCount++;
							depthHist[depth]++;
							allDepthHist[depth]++;
						}
						depth = depth - depthItr->second.ends;
					}
					else {
						if ((pos > bedItr->start) && (pos <= bedItr->end)) {	
							if (depth == 0) zeroDepthCount++;
							depthHist[depth]++;
							allDepthHist[depth]++;
						}
					}
				}

				// Report the coverage for the current interval.
				int length = bedItr->end - bedItr->start;
				totalLength += length;
				
				int nonZeroBases =  (length-zeroDepthCount);
				float fractCovered = (float) nonZeroBases /length;
				
				if (this->writeHistogram == false) {
					bedB->reportBedTab(*bedItr);
					printf("%d\t%d\t%d\t%0.7f\n", bedItr->count, nonZeroBases, length, fractCovered);
				}
				else {
					map<unsigned int, unsigned int>::const_iterator histItr = depthHist.begin();
					map<unsigned int, unsigned int>::const_iterator histEnd = depthHist.end();
					for (; histItr != histEnd; ++histItr) {
						float fractAtThisDepth = (float) histItr->second / length;
						bedB->reportBedTab(*bedItr);
						printf("%d\t%d\t%d\t%0.7f\n", histItr->first, histItr->second, length, fractAtThisDepth);
					}
				}
			}
		}
	}
	if (this->writeHistogram == true) {
		map<unsigned int, unsigned int>::const_iterator histItr = allDepthHist.begin();
		map<unsigned int, unsigned int>::const_iterator histEnd = allDepthHist.end();
		for (; histItr != histEnd; ++histItr) {
			float fractAtThisDepth = (float) histItr->second / totalLength;
			printf("all\t%d\t%d\t%d\t%0.7f\n", histItr->first, histItr->second, totalLength, fractAtThisDepth);
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
		CollectCoverage(beds);
	}
	else {   						// process stdin
		CollectCoverage(cin);		
	}
}


