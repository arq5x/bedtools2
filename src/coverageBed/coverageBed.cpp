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


BedGraph::BedGraph(string &bedAFile, string &bedBFile) {

	this->bedAFile = bedAFile;
	this->bedBFile = bedBFile;
	
	this->bedA = new BedFile(bedAFile);
	this->bedB = new BedFile(bedBFile);
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

	ifstream bed(bedA->bedFile.c_str(), ios::in);
	if ( !bed ) {
		cerr << "Error: The requested bed file (" <<bedA->bedFile << ") could not be opened. Exiting!" << endl;
		exit (1);
	}
	
	BED a;
	while (getline(bed, bedLine)) {
		
		if ((bedLine.find_first_of("track") == 1) || (bedLine.find_first_of("browser") == 1)) {
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
				bedB->countHits(bedB->bedMap[a.chrom], a.start, a.end);
			}
		}
	}
	
	// now, report the count of hist for each feature in B.
	for (masterBedMap::iterator c = bedB->bedMap.begin(); c != bedB->bedMap.end(); ++c) {
		map<int, vector<BED> > bin2Beds = c->second;
		for (map<int, vector<BED> >::iterator b = bin2Beds.begin(); b != bin2Beds.end(); ++b) {
			vector<BED> beds = b->second;
			for (unsigned int i = 0; i < beds.size(); i++) {
				cout << c->first << "\t" << beds[i].start << "\t" << beds[i].end << "\t" << beds[i].count << endl;
			}
		}
	}
}


