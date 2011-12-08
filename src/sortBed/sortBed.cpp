/*****************************************************************************
  sortBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "sortBed.h"

//
// Constructor
//
BedSort::BedSort(string &bedFile, bool printHeader) {
    _bedFile = bedFile;
    _bed = new BedFile(bedFile);
    
    _bed->loadBedFileIntoMapNoBin();
    // report the header first if asked.
    if (printHeader == true) {
        _bed->PrintHeader();
    }
}

//
// Destructor
//
BedSort::~BedSort(void) {
}


void BedSort::SortBed() {

    // loop through each chromosome and merge their BED entries
    for (masterBedMapNoBin::iterator m = _bed->bedMapNoBin.begin(); m != _bed->bedMapNoBin.end(); ++m) {

        // bedList is already sorted by start position.
        vector<BED> bedList = m->second;

        for (unsigned int i = 0; i < bedList.size(); ++i) {
            _bed->reportBedNewLine(bedList[i]);
        }
    }
}


void BedSort::SortBedBySizeAsc() {

    vector<BED> masterList;
    masterList.reserve(1000000);

    // loop through each chromosome and merge their BED entries
    for (masterBedMapNoBin::iterator m = _bed->bedMapNoBin.begin(); m != _bed->bedMapNoBin.end(); ++m) {

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
        _bed->reportBedNewLine(masterList[i]);
    }
}


void BedSort::SortBedBySizeDesc() {

    vector<BED> masterList;
    masterList.reserve(1000000);

    // loop through each chromosome and merge their BED entries
    for (masterBedMapNoBin::iterator m = _bed->bedMapNoBin.begin(); m != _bed->bedMapNoBin.end(); ++m) {

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
        _bed->reportBedNewLine(masterList[i]);
    }
}

void BedSort::SortBedByChromThenSizeAsc() {

    // loop through each chromosome and merge their BED entries
    for (masterBedMapNoBin::iterator m = _bed->bedMapNoBin.begin(); m != _bed->bedMapNoBin.end(); ++m) {

        // bedList is already sorted by start position.
        vector<BED> bedList = m->second;
        sort(bedList.begin(), bedList.end(), sortBySizeAsc);

        for (unsigned int i = 0; i < bedList.size(); ++i) {
            _bed->reportBedNewLine(bedList[i]);
        }
    }
}


void BedSort::SortBedByChromThenSizeDesc() {

    // loop through each chromosome and merge their BED entries
    for (masterBedMapNoBin::iterator m = _bed->bedMapNoBin.begin(); m != _bed->bedMapNoBin.end(); ++m) {

        // bedList is already sorted by start position.
        vector<BED> bedList = m->second;

        sort(bedList.begin(), bedList.end(), sortBySizeDesc);

        for (unsigned int i = 0; i < bedList.size(); ++i) {
            _bed->reportBedNewLine(bedList[i]);
        }
    }
}


void BedSort::SortBedByChromThenScoreAsc() {

    if (_bed->bedType >= 5) {
        // loop through each chromosome and merge their BED entries
        for (masterBedMapNoBin::iterator m = _bed->bedMapNoBin.begin(); m != _bed->bedMapNoBin.end(); ++m) {

            // bedList is already sorted by start position.
            vector<BED> bedList = m->second;
            sort(bedList.begin(), bedList.end(), sortByScoreAsc);

            for (unsigned int i = 0; i < bedList.size(); ++i) {
                _bed->reportBedNewLine(bedList[i]);
            }
        }
    }
    else {
        cerr << "Error: Requested a sort by score, but your BED file does not appear to be in BED 5 format or greater.  Exiting." << endl;
        exit(1);
    }
}


void BedSort::SortBedByChromThenScoreDesc() {

    if (_bed->bedType >= 5) {
        // loop through each chromosome and merge their BED entries
        for (masterBedMapNoBin::iterator m = _bed->bedMapNoBin.begin(); m != _bed->bedMapNoBin.end(); ++m) {

            // bedList is already sorted by start position.
            vector<BED> bedList = m->second;
            sort(bedList.begin(), bedList.end(), sortByScoreDesc);

            for (unsigned int i = 0; i < bedList.size(); ++i) {
                _bed->reportBedNewLine(bedList[i]);
            }
        }
    }
    else {
        cerr << "Error: Requested a sort by score, but your BED file does not appear to be in BED 5 format or greater.  Exiting." << endl;
        exit(1);
    }
}

