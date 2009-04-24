// 
//  sortBed.cpp
//  BEDTools
//  
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//  Summary:  Sorts a BED file in ascending order by chrom then by start position.
//

#include "sortBed.h"

//
// Constructor
//
BedSort::BedSort(string &bedFile) {
	this->bedFile = bedFile;
	this->bed = new BedFile(bedFile);
}

//
// Destructor
//
BedSort::~BedSort(void) {
}

//
// Merge overlapping BED entries into a single entry 
//
void BedSort::SortBed() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bed->loadBedFileIntoMapNoBin();

	// loop through each chromosome and merge their BED entries
	for (masterBedMapNoBin::iterator m = bed->bedMapNoBin.begin(); m != bed->bedMapNoBin.end(); ++m) {

		// bedList is already sorted by start position.
		vector<BED> bedList = m->second; 

		for (int i = 0; i < bedList.size(); ++i) {
			cout << bedList[i].chrom << "\t" << bedList[i].start << "\t" << bedList[i].end << endl;
		}
	}
}

