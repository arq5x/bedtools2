/*****************************************************************************
  annotateBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "annotateBed.h"


BedAnnotate::BedAnnotate (const string &bedFile, const vector<string> &annotationFiles, bool &forceStrand) {

	this->bedFile = bedFile;
	this->annotationFiles = annotationFiles;
	
	this->bed = new BedFile(bedFile);
	this->forceStrand = forceStrand;
}



BedAnnotate::~BedAnnotate (void) {
}


void BedAnnotate::ProcessAnnotations (istream &bedInput) {
	
	// loop through each of the annotation files and compute the coverage
	// of each with respect to the input BED file.
	vector<string>::const_iterator annotItr = this->annotationFiles.begin();
	vector<string>::const_iterator annotEnd = this->annotationFiles.end();
	for (; annotItr != annotEnd; ++annotItr) {
		// create a BED file of the current annotation file.
		BedFile annotation = new BedFile(*annotItr);
		// compute the coverage of the annotation file with respect to the input file.
		GetCoverage(annotation, bedInput);
	}	
}

 
void BedAnnotate::GetCoverage (const BedFile &annotation, istream &bedInput) {
	
	// load the annotation bed file into a map so
	// that we can easily compare "A" to it for overlaps
	annotation->loadBedFileIntoMap();

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
																				
						if ((depth == 0) && (pos > bedItr->start) && (pos <= bedItr->end)) {
							zeroDepthCount++;
						}
						
						depth = depth - dMap[pos].ends;
					}
					else {
						if ((depth == 0) && (pos > bedItr->start) && (pos <= bedItr->end)) {
							zeroDepthCount++;
						}
					}
				}

				// Report the coverage for the current interval.
				int length = bedItr->end - bedItr->start;
				int nonZeroBases =  (length-zeroDepthCount);
				float fractCovered = (float) nonZeroBases /length;
				
				bedB->reportBedTab(*bedItr);
				printf("%d\t%d\t%d\t%0.7f\n", bedItr->count, nonZeroBases, length, fractCovered);			
			}
		}
	}
}


void BedAnnotate::DetermineBedInput () {
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


