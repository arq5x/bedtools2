#ifndef BEDFILEPE_H
#define BEDFILEPE_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <algorithm>
#include "bedFile.h"
#include "lineFileUtilities.h"

using namespace std;


/*
    Structure for paired-end records
*/
struct BEDPE {

    // UCSC BED fields
    string chrom1;
    CHRPOS start1;
    CHRPOS end1;

    string chrom2;
    CHRPOS start2;
    CHRPOS end2;

    string name;
    string score;

    string strand1;
    string strand2;

    // all of the original fields in the record
    vector<string> fields;
    // indices of the "other" fields
    vector<uint16_t> other_idxs;
};




//************************************************
// BedFile Class methods and elements
//************************************************
class BedFilePE {

public:

    // Constructor
    BedFilePE(string &);

    // Destructor
    ~BedFilePE(void);

    // Open a BEDPE file for reading (creates an istream pointer)
    void Open(void);

    // Close an opened BEDPE file.
    void Close(void);

    // Get the next BED entry in an opened BED file.
    BedLineStatus GetNextBedPE (BEDPE &bedpe, int &lineNum);


    // Methods

    void reportBedPETab(const BEDPE &a);
    void reportBedPENewLine(const BEDPE &a);
    void loadBedPEFileIntoMap();
    void splitBedPEIntoBeds(const BEDPE &a, const int &lineNum, MATE *bedEntry1, MATE *bedEntry2);


    void FindOverlapsPerBin(int bEnd, string chrom, CHRPOS start, CHRPOS end, string name, string strand,
        vector<MATE> &hits, float overlapFraction, bool forceStrand, bool enforceDiffNames);


    string bedFile;
    unsigned int bedType;

    masterMateMap bedMapEnd1;
    masterMateMap bedMapEnd2;

private:
    istream *_bedStream;

    // methods
    BedLineStatus parseLine (BEDPE &bedpe, const vector<string> &lineVector, int &lineNum);
    bool parseBedPELine (BEDPE &bed, const vector<string> &lineVector, const int &lineNum);
};

#endif /* BEDFILEPE_H */
