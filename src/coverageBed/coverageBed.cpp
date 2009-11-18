// 
//  coverageBed.cpp
//  BEDTools
//  
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//
#include "lineFileUtilities.h"
#include "coverageBed.h"


BedGraph::BedGraph(string &bedAFile, string &bedBFile, bool &forceStrand) {

	this->bedAFile = bedAFile;
	this->bedBFile = bedBFile;
	
	this->bedA = new BedFile(bedAFile);
	this->bedB = new BedFile(bedBFile);
	
	this->forceStrand = forceStrand;
}



BedGraph::~BedGraph(void) {
}


 
void BedGraph::GraphBed() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bedB->loadBedFileIntoMap();


	string bedLine;
	BED bedEntry;                                                                                                                        
	int lineNum = 0;

	// are we dealing with a file?
	if (bedA->bedFile != "stdin") {
		
		ifstream bed(bedA->bedFile.c_str(), ios::in);
		if ( !bed ) {
			cerr << "Error: The requested bed file (" <<bedA->bedFile << ") could not be opened. Exiting!" << endl;
			exit (1);
		}
	
		BED a;
		while (getline(bed, bedLine)) {
		
			if ((bedLine.find("track") != string::npos) || (bedLine.find("browser") != string::npos)) {
				continue;
			}
			else {
				vector<string> bedFields;
				Tokenize(bedLine,bedFields);
				lineNum++;
				// process the feature in A IFF it's valid.
				if (bedA->parseBedLine(a, bedFields, lineNum)) {
			
					// increment the count of overlaps for each feature in B that
					// overlaps the current A interval
					bedB->countHits(bedB->bedMap[a.chrom], a, this->forceStrand);
				}
			}
		}
	}
	// A is being passed via stdin.
	else {
	
		BED a;
		while (getline(cin, bedLine)) {
		
			if ((bedLine.find("track") != string::npos) || (bedLine.find("browser") != string::npos)) {
				continue;
			}
			else {
				vector<string> bedFields;
				Tokenize(bedLine,bedFields);
				lineNum++;
				// process the feature in A IFF it's valid.
				if (bedA->parseBedLine(a, bedFields, lineNum)) {
			
					// increment the count of overlaps for each feature in B that
					// overlaps the current A interval
					bedB->countHits(bedB->bedMap[a.chrom], a, this->forceStrand);
				}
			}
		}
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


