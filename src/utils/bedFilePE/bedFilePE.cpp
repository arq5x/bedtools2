//
//  bedFilePE.cpp
//  BEDTools
//
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//  Summary:  Contains common functions for finding BED overlaps.
//
//  Acknowledgments: Much of the code herein is taken from Jim Kent's
//                   BED processing code.  I am grateful for his elegant
//                   genome binning algorithm and therefore use it extensively.


#include "bedFilePE.h"

using namespace std;

// Constructor
BedFilePE::BedFilePE(string &bedFile) {
    this->bedFile = bedFile;
}

// Destructor
BedFilePE::~BedFilePE(void) {
}

void BedFilePE::Open(void) {
    if (bedFile == "stdin" || bedFile == "-") {
        _bedStream = &cin;
    }
    else {
        _bedStream = new ifstream(bedFile.c_str(), ios::in);

        if (isGzipFile(_bedStream) == true) {
            delete _bedStream;
            _bedStream = new igzstream(bedFile.c_str(), ios::in);
        }
        // can we open the file?
        if ( _bedStream->fail() ) {
            cerr << "Error: The requested BEDPE file (" 
                 << bedFile
                 << ") " 
                 << "could not be opened. "
                 << "Error message: ("
                 << strerror(errno)
                 << "). Exiting!" << endl;
            exit (1);
        }
    }
}



// Close the BEDPE file
void BedFilePE::Close(void) {
    if (bedFile != "stdin" && bedFile != "-") delete _bedStream;
}


BedLineStatus BedFilePE::GetNextBedPE (BEDPE &bedpe, int &lineNum) {

    // make sure there are still lines to process.
    // if so, tokenize, validate and return the BEDPE entry.
    if (_bedStream->good()) {
        string bedPELine;
        vector<string> bedPEFields;
        bedPEFields.reserve(10);

        // parse the bedStream pointer
        getline(*_bedStream, bedPELine);
        lineNum++;

        // split into a string vector.
        Tokenize(bedPELine,bedPEFields);

        // load the BEDPE struct as long as it's a valid BEDPE entry.
        return parseLine(bedpe, bedPEFields, lineNum);
    }
    // default if file is closed or EOF
    return BED_INVALID;
}


/*
    reportBedPETab

    Writes the _original_ BED entry for A.
    Works for BEDPE only.
*/
void BedFilePE::reportBedPETab(const BEDPE &a) {

    if (this->bedType == 6) {
        printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t",
                                            a.chrom1.c_str(), a.start1, a.end1,
                                            a.chrom2.c_str(), a.start2, a.end2);
    }
    else if (this->bedType == 7) {
        printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t",
                                            a.chrom1.c_str(), a.start1, a.end1,
                                            a.chrom2.c_str(), a.start2, a.end2,
                                            a.name.c_str());
    }
    else if (this->bedType == 8) {
        printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t",
                                            a.chrom1.c_str(), a.start1, a.end1,
                                            a.chrom2.c_str(), a.start2, a.end2,
                                            a.name.c_str(), a.score.c_str());
    }
    else if (this->bedType == 10) {
        printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\t%s\t",
                                            a.chrom1.c_str(), a.start1, a.end1,
                                            a.chrom2.c_str(), a.start2, a.end2,
                                            a.name.c_str(), a.score.c_str(), a.strand1.c_str(), a.strand2.c_str());
    }
    else if (this->bedType > 10) {
        printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\t%s",
                                            a.chrom1.c_str(), a.start1, a.end1,
                                            a.chrom2.c_str(), a.start2, a.end2,
                                            a.name.c_str(), a.score.c_str(), a.strand1.c_str(), a.strand2.c_str());
        vector<uint16_t>::const_iterator othIt  = a.other_idxs.begin();
        vector<uint16_t>::const_iterator othEnd = a.other_idxs.end();
        for ( ; othIt != othEnd; ++othIt) {
            printf("\t%s", a.fields[*othIt].c_str());
        }
        printf("\t");
    }
}



/*
    reportBedPENewLine

    Writes the _original_ BED entry for A.
    Works for BEDPE only.
*/
void BedFilePE::reportBedPENewLine(const BEDPE &a) {

    if (this->bedType == 6) {
        printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\n",
                                            a.chrom1.c_str(), a.start1, a.end1,
                                            a.chrom2.c_str(), a.start2, a.end2);
    }
    else if (this->bedType == 7) {
        printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\n",
                                            a.chrom1.c_str(), a.start1, a.end1,
                                            a.chrom2.c_str(), a.start2, a.end2,
                                            a.name.c_str());
    }
    else if (this->bedType == 8) {
        printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\n",
                                            a.chrom1.c_str(), a.start1, a.end1,
                                            a.chrom2.c_str(), a.start2, a.end2,
                                            a.name.c_str(), a.score.c_str());
    }
    else if (this->bedType == 10) {
        printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\t%s\n",
                                            a.chrom1.c_str(), a.start1, a.end1,
                                            a.chrom2.c_str(), a.start2, a.end2,
                                            a.name.c_str(), a.score.c_str(), a.strand1.c_str(), a.strand2.c_str());
    }
    else if (this->bedType > 10) {
        printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\t%s",
                                            a.chrom1.c_str(), a.start1, a.end1,
                                            a.chrom2.c_str(), a.start2, a.end2,
                                            a.name.c_str(), a.score.c_str(), a.strand1.c_str(), a.strand2.c_str());
        vector<uint16_t>::const_iterator othIt  = a.other_idxs.begin();
        vector<uint16_t>::const_iterator othEnd = a.other_idxs.end();
        for ( ; othIt != othEnd; ++othIt) {
            printf("\t%s", a.fields[*othIt].c_str());
        }
        printf("\n");
    }
}


BedLineStatus BedFilePE::parseLine (BEDPE &bedpe, const vector<string> &lineVector, int &lineNum) {

    // bail out if we have a blank line
    if (lineVector.empty())
        return BED_BLANK;

    if ((lineVector[0].find("track") == string::npos) && (lineVector[0].find("browser") == string::npos) && (lineVector[0].find("#") == string::npos) ) {
        // we need at least 6 columns
        if (lineVector.size() >= 6) {
            if (parseBedPELine(bedpe, lineVector, lineNum) == true)
                return BED_VALID;
            else return BED_INVALID;
        }
        else  {
            cerr << "It looks as though you have less than 6 columns.  Are you sure your files are tab-delimited?" << endl;
            exit(1);
        }
    }
    else {
        lineNum--;
        return BED_HEADER;
    }

    // default
    return BED_INVALID;
}


bool BedFilePE::parseBedPELine (BEDPE &bed, const vector<string> &lineVector, const int &lineNum) {

    if ((lineNum == 1) && (lineVector.size() >= 6)) {

        this->bedType = lineVector.size();
        bed.fields = lineVector;
        if (this->bedType == 6) {
            bed.chrom1 = lineVector[0];
            bed.start1 = stoll(lineVector[1].c_str());
            bed.end1 = stoll(lineVector[2].c_str());

            bed.chrom2 = lineVector[3];
            bed.start2 = stoll(lineVector[4].c_str());
            bed.end2 = stoll(lineVector[5].c_str());

            return true;
        }
        else if (this->bedType == 7) {
            bed.chrom1 = lineVector[0];
            bed.start1 = stoll(lineVector[1].c_str());
            bed.end1 = stoll(lineVector[2].c_str());

            bed.chrom2 = lineVector[3];
            bed.start2 = stoll(lineVector[4].c_str());
            bed.end2 = stoll(lineVector[5].c_str());

            bed.name = lineVector[6];
            return true;
        }
        else if (this->bedType == 8) {
            bed.chrom1 = lineVector[0];
            bed.start1 = stoll(lineVector[1].c_str());
            bed.end1 = stoll(lineVector[2].c_str());

            bed.chrom2 = lineVector[3];
            bed.start2 = stoll(lineVector[4].c_str());
            bed.end2 = stoll(lineVector[5].c_str());

            bed.name = lineVector[6];
            bed.score = lineVector[7].c_str();
            return true;
        }
        else if (this->bedType == 10) {
            bed.chrom1 = lineVector[0];
            bed.start1 = stoll(lineVector[1].c_str());
            bed.end1 = stoll(lineVector[2].c_str());

            bed.chrom2 = lineVector[3];
            bed.start2 = stoll(lineVector[4].c_str());
            bed.end2 = stoll(lineVector[5].c_str());

            bed.name = lineVector[6];
            bed.score = lineVector[7].c_str();

            bed.strand1 = lineVector[8];
            bed.strand2 = lineVector[9];

            return true;
        }
        else if (this->bedType > 10) {
            bed.chrom1 = lineVector[0];
            bed.start1 = stoll(lineVector[1].c_str());
            bed.end1 = stoll(lineVector[2].c_str());

            bed.chrom2 = lineVector[3];
            bed.start2 = stoll(lineVector[4].c_str());
            bed.end2 = stoll(lineVector[5].c_str());

            bed.name = lineVector[6];
            bed.score = lineVector[7].c_str();

            bed.strand1 = lineVector[8];
            bed.strand2 = lineVector[9];

            for (unsigned int i = 10; i < lineVector.size(); ++i) {
                bed.other_idxs.push_back(i);
            }
            return true;
        }
        else {
            cerr << "Unexpected number of fields: " << lineNum << ".  Verify that your files are TAB-delimited and that your BEDPE file has 6,7,8 or 10 fields.  Exiting..." << endl;
            exit(1);
        }

        if (bed.start1 > bed.end1) {
            cerr << "Error: malformed BEDPE entry at line " << lineNum << ". Start1 was greater than End1. Ignoring it and moving on." << endl;
            return false;
        }
        else if (bed.start2 > bed.end2) {
            cerr << "Error: malformed BEDPE entry at line " << lineNum << ". Start2 was greater than End2. Ignoring it and moving on." << endl;
            return false;
        }
    }
    else if ( (lineNum > 1) && (lineVector.size() == this->bedType)) {

        bed.fields = lineVector;
        
        if (this->bedType == 6) {
            bed.chrom1 = lineVector[0];
            bed.start1 = stoll(lineVector[1].c_str());
            bed.end1 = stoll(lineVector[2].c_str());

            bed.chrom2 = lineVector[3];
            bed.start2 = stoll(lineVector[4].c_str());
            bed.end2 = stoll(lineVector[5].c_str());

            return true;
        }
        else if (this->bedType == 7) {
            bed.chrom1 = lineVector[0];
            bed.start1 = stoll(lineVector[1].c_str());
            bed.end1 = stoll(lineVector[2].c_str());

            bed.chrom2 = lineVector[3];
            bed.start2 = stoll(lineVector[4].c_str());
            bed.end2 = stoll(lineVector[5].c_str());

            bed.name = lineVector[6];
            return true;
        }
        else if (this->bedType == 8) {
            bed.chrom1 = lineVector[0];
            bed.start1 = stoll(lineVector[1].c_str());
            bed.end1 = stoll(lineVector[2].c_str());

            bed.chrom2 = lineVector[3];
            bed.start2 = stoll(lineVector[4].c_str());
            bed.end2 = stoll(lineVector[5].c_str());

            bed.name = lineVector[6];
            bed.score = lineVector[7].c_str();
            return true;
        }
        else if (this->bedType == 10) {
            bed.chrom1 = lineVector[0];
            bed.start1 = stoll(lineVector[1].c_str());
            bed.end1 = stoll(lineVector[2].c_str());

            bed.chrom2 = lineVector[3];
            bed.start2 = stoll(lineVector[4].c_str());
            bed.end2 = stoll(lineVector[5].c_str());

            bed.name = lineVector[6];
            bed.score = lineVector[7].c_str();

            bed.strand1 = lineVector[8];
            bed.strand2 = lineVector[9];

            return true;
        }
        else if (this->bedType > 10) {
            bed.chrom1 = lineVector[0];
            bed.start1 = stoll(lineVector[1].c_str());
            bed.end1 = stoll(lineVector[2].c_str());

            bed.chrom2 = lineVector[3];
            bed.start2 = stoll(lineVector[4].c_str());
            bed.end2 = stoll(lineVector[5].c_str());

            bed.name = lineVector[6];
            bed.score = lineVector[7].c_str();

            bed.strand1 = lineVector[8];
            bed.strand2 = lineVector[9];
            for (unsigned int i = 10; i < lineVector.size(); ++i) {
                bed.other_idxs.push_back(i);
            }
            return true;
        }
        else {
            cerr << "Unexpected number of fields: " << lineNum << ".  Verify that your files are TAB-delimited and that your BEDPE file has 6,7,8 or 10 fields.  Exiting..." << endl;
            exit(1);
        }

        if (bed.start1 > bed.end1) {
            cerr << "Error: malformed BED entry at line " << lineNum << ". Start1 was greater than End1. Ignoring it and moving on." << endl;
            return false;
        }
        else if (bed.start2 > bed.end2) {
            cerr << "Error: malformed BED entry at line " << lineNum << ". Start2 was greater than End2. Ignoring it and moving on." << endl;
            return false;
        }
    }
    else if (lineVector.size() == 1) {
        cerr << "Only one BED field detected: " << lineNum << ".  Verify that your files are TAB-delimited.  Exiting..." << endl;
        exit(1);
    }
    else if ((lineVector.size() != this->bedType) && (lineVector.size() != 0)) {
        cerr << "Differing number of BEDPE fields encountered at line: " << lineNum << ".  Exiting..." << endl;
        exit(1);
    }
    else if ((lineVector.size() < 6) && (lineVector.size() != 0)) {
        cerr << "TAB delimited BEDPE file with at least 6 fields (chrom1, start1, end1, chrom2, start2, end2) is required at line: "<< lineNum << ".  Exiting..." << endl;
        exit(1);
    }
    return false;
}


/*
    Adapted from kent source "binKeeperFind"
*/
void BedFilePE::FindOverlapsPerBin(int bEnd, string chrom, CHRPOS start, CHRPOS end, string name, string strand,
                                   vector<MATE> &hits, float overlapFraction, bool forceStrand, bool enforceDiffNames) {

    CHRPOS startBin, endBin;
    startBin = (start >> _binFirstShift);
    endBin = ((end-1) >> _binFirstShift);

    // loop through each bin "level" in the binning hierarchy
    for (int i = 0; i < _binLevels; ++i) {

        // loop through each bin at this level of the hierarchy
        CHRPOS offset = _binOffsetsExtended[i];
        for (CHRPOS j = (startBin+offset); j <= (endBin+offset); ++j)  {

            // loop through each feature in this chrom/bin and see if it overlaps
            // with the feature that was passed in.  if so, add the feature to
            // the list of hits.
            vector<MATE>::const_iterator bedItr;
            vector<MATE>::const_iterator bedEnd;
            if (bEnd == 1) {
                bedItr = bedMapEnd1[chrom][j].begin();
                bedEnd = bedMapEnd1[chrom][j].end();
            }
            else if (bEnd == 2) {
                bedItr = bedMapEnd2[chrom][j].begin();
                bedEnd = bedMapEnd2[chrom][j].end();
            }
            else {
                cerr << "Unexpected end of B requested" << endl;
            }
            for (; bedItr != bedEnd; ++bedItr) {
                float overlap = overlaps(bedItr->bed.start, bedItr->bed.end, start, end);
                float size    = (float)(end - start);

                if ( (overlap / size) >= overlapFraction ) {

                    // skip the hit if not on the same strand (and we care)
                    if ((forceStrand == false) && (enforceDiffNames == false)) {
                        hits.push_back(*bedItr);    // it's a hit, add it.
                    }
                    else if ((forceStrand == true) && (enforceDiffNames == false)) {
                        if (strand == bedItr->bed.strand)
                            hits.push_back(*bedItr);    // it's a hit, add it.
                    }
                    else if ((forceStrand == true) && (enforceDiffNames == true)) {
                        if ((strand == bedItr->bed.strand) && (name != bedItr->bed.name))
                            hits.push_back(*bedItr);    // it's a hit, add it.
                    }
                    else if ((forceStrand == false) && (enforceDiffNames == true)) {
                        if (name != bedItr->bed.name)
                            hits.push_back(*bedItr);    // it's a hit, add it.
                    }
                }

            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
}


void BedFilePE::loadBedPEFileIntoMap() {

    int lineNum = 0;
    int bin1, bin2;
    BedLineStatus bedStatus;
    BEDPE bedpeEntry, nullBedPE;

    Open();
    bedStatus = this->GetNextBedPE(bedpeEntry, lineNum);
    while (bedStatus != BED_INVALID) {

        if (bedStatus == BED_VALID) {
            MATE *bedEntry1 = new MATE();
            MATE *bedEntry2 = new MATE();
            // separate the BEDPE entry into separate
            // BED entries
            splitBedPEIntoBeds(bedpeEntry, lineNum, bedEntry1, bedEntry2);

            // load end1 into a UCSC bin map
            bin1 = getBin(bedEntry1->bed.start, bedEntry1->bed.end);
            this->bedMapEnd1[bedEntry1->bed.chrom][bin1].push_back(*bedEntry1);

            // load end2 into a UCSC bin map
            bin2 = getBin(bedEntry2->bed.start, bedEntry2->bed.end);
            this->bedMapEnd2[bedEntry2->bed.chrom][bin2].push_back(*bedEntry2);

            bedpeEntry = nullBedPE;
        }
        bedStatus = this->GetNextBedPE(bedpeEntry, lineNum);
    }
    Close();
}


void BedFilePE::splitBedPEIntoBeds(const BEDPE &bedpeEntry, const int &lineNum, MATE *bedEntry1, MATE *bedEntry2) {

    /*
       Split the BEDPE entry into separate BED entries

       NOTE: I am using a trick here where I store
       the lineNum of the BEDPE from the original file
       in the "count" column.  This allows me to later
       resolve whether the hits found on both ends of BEDPE A
       came from the same entry in BEDPE B.  Tracking by "name"
       alone with fail when there are multiple mappings for a given
       read-pair.
    */

    bedEntry1->bed.chrom           = bedpeEntry.chrom1;
    bedEntry1->bed.start           = bedpeEntry.start1;
    bedEntry1->bed.end             = bedpeEntry.end1;
    bedEntry1->bed.name            = bedpeEntry.name;
    bedEntry1->bed.score           = bedpeEntry.score;        // only store the score in end1 to save memory
    bedEntry1->bed.strand          = bedpeEntry.strand1;
    bedEntry1->bed.fields          = bedpeEntry.fields;       // only store the fields in end1 to save memory
    bedEntry1->bed.other_idxs      = bedpeEntry.other_idxs;   // only store the other_idxs in end1 to save memory
    bedEntry1->lineNum             = lineNum;
    bedEntry1->mate                = bedEntry2;               // keep a pointer to end2

    bedEntry2->bed.chrom           = bedpeEntry.chrom2;
    bedEntry2->bed.start           = bedpeEntry.start2;
    bedEntry2->bed.end             = bedpeEntry.end2;
    bedEntry2->bed.name            = bedpeEntry.name;
    bedEntry2->bed.strand          = bedpeEntry.strand2;
    bedEntry2->lineNum             = lineNum;
    bedEntry2->mate                = bedEntry1;               // keep a pointer to end1
}



