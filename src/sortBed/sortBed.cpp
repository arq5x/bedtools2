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


/*
	reportBed
	
	Writes the _original_ BED entry for A.
	Works for BED3 - BED6.
*/
void BedSort::reportBed(const BED &a) {
	
	if (bed->bedType == 3) {
		cout << a.chrom << "\t" << a.start << "\t" << a.end;
	}
	else if (bed->bedType == 4) {
		cout << a.chrom << "\t" << a.start << "\t" << a.end << "\t"
		<< a.name;
	}
	else if (bed->bedType == 5) {
		cout << a.chrom << "\t" << a.start << "\t" << a.end << "\t"
		<< a.name << "\t" << a.score;
	}
	else if (bed->bedType == 6) {
		cout << a.chrom << "\t" << a.start << "\t" << a.end << "\t" 
		<< a.name << "\t" << a.score << "\t" << a.strand;
	}
}

void BedSort::SortBed() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bed->loadBedFileIntoMapNoBin();

	// loop through each chromosome and merge their BED entries
	for (masterBedMapNoBin::iterator m = bed->bedMapNoBin.begin(); m != bed->bedMapNoBin.end(); ++m) {

		// bedList is already sorted by start position.
		vector<BED> bedList = m->second; 

		for (unsigned int i = 0; i < bedList.size(); ++i) {
			reportBed(bedList[i]); cout << "\n";
			///cout << bedList[i].chrom << "\t" << bedList[i].start << "\t" << bedList[i].end << endl;
		}
	}
}


void BedSort::SortBedBySizeAsc() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bed->loadBedFileIntoMapNoBin();

	vector<BED> masterList;
	masterList.reserve(1000000);
	
	// loop through each chromosome and merge their BED entries
	for (masterBedMapNoBin::iterator m = bed->bedMapNoBin.begin(); m != bed->bedMapNoBin.end(); ++m) {

		// bedList is already sorted by start position.
		vector<BED> bedList = m->second; 

		// add the entries from this chromosome to the current list
		for (unsigned int i = 0; i < m->second.size(); ++i) {
			masterList.push_back(m->second[i]);
		}
	}
	
	// sort the master list by size (asc.)
	sort(masterList.begin(), masterList.end(), sortBySizeAsc);
	
	// report the entries in ascending order
	for (unsigned int i = 0; i < masterList.size(); ++i) {
		reportBed(masterList[i]); cout << "\n";
		//cout << masterList[i].chrom << "\t" << masterList[i].start << "\t" << masterList[i].end << endl;
	}
}


void BedSort::SortBedBySizeDesc() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bed->loadBedFileIntoMapNoBin();

	vector<BED> masterList;
	masterList.reserve(1000000);
	
	// loop through each chromosome and merge their BED entries
	for (masterBedMapNoBin::iterator m = bed->bedMapNoBin.begin(); m != bed->bedMapNoBin.end(); ++m) {

		// bedList is already sorted by start position.
		vector<BED> bedList = m->second; 

		// add the entries from this chromosome to the current list
		for (unsigned int i = 0; i < m->second.size(); ++i) {
			masterList.push_back(m->second[i]);
		}
	}
	
	// sort the master list by size (asc.)
	sort(masterList.begin(), masterList.end(), sortBySizeDesc);
	
	// report the entries in ascending order
	for (unsigned int i = 0; i < masterList.size(); ++i) {
		reportBed(masterList[i]); cout << "\n";
		//cout << masterList[i].chrom << "\t" << masterList[i].start << "\t" << masterList[i].end << endl;
	}
}

void BedSort::SortBedByChromThenSizeAsc() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bed->loadBedFileIntoMapNoBin();

	// loop through each chromosome and merge their BED entries
	for (masterBedMapNoBin::iterator m = bed->bedMapNoBin.begin(); m != bed->bedMapNoBin.end(); ++m) {

		// bedList is already sorted by start position.
		vector<BED> bedList = m->second; 
		sort(bedList.begin(), bedList.end(), sortBySizeAsc);
		
		for (unsigned int i = 0; i < bedList.size(); ++i) {
			reportBed(bedList[i]); cout << "\n";
			//cout << bedList[i].chrom << "\t" << bedList[i].start << "\t" << bedList[i].end << endl;
		}
	}
}


void BedSort::SortBedByChromThenSizeDesc() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bed->loadBedFileIntoMapNoBin();

	// loop through each chromosome and merge their BED entries
	for (masterBedMapNoBin::iterator m = bed->bedMapNoBin.begin(); m != bed->bedMapNoBin.end(); ++m) {

		// bedList is already sorted by start position.
		vector<BED> bedList = m->second; 

		sort(bedList.begin(), bedList.end(), sortBySizeDesc);
		
		for (unsigned int i = 0; i < bedList.size(); ++i) {
			reportBed(bedList[i]); cout << "\n";
			//cout << bedList[i].chrom << "\t" << bedList[i].start << "\t" << bedList[i].end << endl;
		}
	}
}


void BedSort::SortBedByChromThenScoreAsc() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bed->loadBedFileIntoMapNoBin();

	if (bed->bedType >= 5) {
		// loop through each chromosome and merge their BED entries
		for (masterBedMapNoBin::iterator m = bed->bedMapNoBin.begin(); m != bed->bedMapNoBin.end(); ++m) {

			// bedList is already sorted by start position.
			vector<BED> bedList = m->second; 
			sort(bedList.begin(), bedList.end(), sortByScoreAsc);
			
			for (unsigned int i = 0; i < bedList.size(); ++i) {
				reportBed(bedList[i]); cout << "\n";
				//cout << bedList[i].chrom << "\t" << bedList[i].start << "\t" << bedList[i].end << endl;
			}
		}
	}
	else {
		cerr << "Error: Requested a sort by score, but your BED file does not appear to be in BED 5 format or greater.  Exiting." << endl;
		exit(1);
	}
}


void BedSort::SortBedByChromThenScoreDesc() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bed->loadBedFileIntoMapNoBin();

	if (bed->bedType >= 5) {
		// loop through each chromosome and merge their BED entries
		for (masterBedMapNoBin::iterator m = bed->bedMapNoBin.begin(); m != bed->bedMapNoBin.end(); ++m) {

			// bedList is already sorted by start position.
			vector<BED> bedList = m->second; 
			sort(bedList.begin(), bedList.end(), sortByScoreDesc);
		
			for (unsigned int i = 0; i < bedList.size(); ++i) {
				reportBed(bedList[i]); cout << "\n";
				//cout << bedList[i].chrom << "\t" << bedList[i].start << "\t" << bedList[i].end << endl;
			}
		}
	}
	else {
		cerr << "Error: Requested a sort by score, but your BED file does not appear to be in BED 5 format or greater.  Exiting." << endl;
		exit(1);
	}
}

