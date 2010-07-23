// ***************************************************************************
// BamMultiReader.cpp (c) 2010 Erik Garrison
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 20 July 2010 (DB)
// ---------------------------------------------------------------------------
// Uses BGZF routines were adapted from the bgzf.c code developed at the Broad
// Institute.
// ---------------------------------------------------------------------------
// Functionality for simultaneously reading multiple BAM files.
//
// This functionality allows applications to work on very large sets of files
// without requiring intermediate merge, sort, and index steps for each file
// subset.  It also improves the performance of our merge system as it
// precludes the need to sort merged files.
// ***************************************************************************

// C++ includes
#include <algorithm>
#include <iterator>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

// BamTools includes
#include "BGZF.h"
#include "BamMultiReader.h"
using namespace BamTools;
using namespace std;

// -----------------------------------------------------
// BamMultiReader implementation
// -----------------------------------------------------

// constructor
BamMultiReader::BamMultiReader(void)
    : CurrentRefID(0)
    , CurrentLeft(0)
{ }

// destructor
BamMultiReader::~BamMultiReader(void) {
    Close(); // close the bam files
    // clean up reader objects
    for (vector<pair<BamReader*, BamAlignment*> >::iterator it = readers.begin(); it != readers.end(); ++it) {
        delete it->first;
        delete it->second;
    }
}

// close the BAM files
void BamMultiReader::Close(void) {
    for (vector<pair<BamReader*, BamAlignment*> >::iterator it = readers.begin(); it != readers.end(); ++it) {
        BamReader* reader = it->first;
        reader->Close();  // close the reader
    }
}

// updates the reference id stored in the BamMultiReader
// to reflect the current state of the readers
void BamMultiReader::UpdateReferenceID(void) {
    // the alignments are sorted by position, so the first alignment will always have the lowest reference ID
    if (alignments.begin()->second.second->RefID != CurrentRefID) {
        // get the next reference id
        // while there aren't any readers at the next ref id
        // increment the ref id
        int nextRefID = CurrentRefID;
        while (alignments.begin()->second.second->RefID != nextRefID) {
            ++nextRefID;
        }
        //cerr << "updating reference id from " << CurrentRefID << " to " << nextRefID << endl;
        CurrentRefID = nextRefID;
    }
}

// checks if any readers still have alignments
bool BamMultiReader::HasOpenReaders() {
    return alignments.size() > 0;
}

// get next alignment among all files
bool BamMultiReader::GetNextAlignment(BamAlignment& nextAlignment) {

    // bail out if we are at EOF in all files, means no more alignments to process
    if (!HasOpenReaders())
        return false;

    // when all alignments have stepped into a new target sequence, update our
    // current reference sequence id
    UpdateReferenceID();

    // our lowest alignment and reader will be at the front of our alignment index
    BamAlignment* alignment = alignments.begin()->second.second;
    BamReader* reader = alignments.begin()->second.first;

    // now that we have the lowest alignment in the set, save it by copy to our argument
    nextAlignment = BamAlignment(*alignment);

    // remove this alignment index entry from our alignment index
    alignments.erase(alignments.begin());

    // and add another entry if we can get another alignment from the reader
    if (reader->GetNextAlignment(*alignment)) {
        alignments.insert(make_pair(make_pair(alignment->RefID, alignment->Position),
                                    make_pair(reader, alignment)));
    } else { // do nothing
        //cerr << "reached end of file " << lowestReader->GetFilename() << endl;
    }

    return true;

}

// get next alignment among all files without parsing character data from alignments
bool BamMultiReader::GetNextAlignmentCore(BamAlignment& nextAlignment) {

    // bail out if we are at EOF in all files, means no more alignments to process
    if (!HasOpenReaders())
        return false;

    // when all alignments have stepped into a new target sequence, update our
    // current reference sequence id
    UpdateReferenceID();

    // our lowest alignment and reader will be at the front of our alignment index
    BamAlignment* alignment = alignments.begin()->second.second;
    BamReader* reader = alignments.begin()->second.first;

    // now that we have the lowest alignment in the set, save it by copy to our argument
    nextAlignment = BamAlignment(*alignment);
    //memcpy(&nextAlignment, alignment, sizeof(BamAlignment));

    // remove this alignment index entry from our alignment index
    alignments.erase(alignments.begin());

    // and add another entry if we can get another alignment from the reader
    if (reader->GetNextAlignmentCore(*alignment)) {
        alignments.insert(make_pair(make_pair(alignment->RefID, alignment->Position), 
                                    make_pair(reader, alignment)));
    } else { // do nothing
        //cerr << "reached end of file " << lowestReader->GetFilename() << endl;
    }

    return true;

}

// jumps to specified region(refID, leftBound) in BAM files, returns success/fail
bool BamMultiReader::Jump(int refID, int position) {

    //if ( References.at(refID).RefHasAlignments && (position <= References.at(refID).RefLength) ) {
    CurrentRefID = refID;
    CurrentLeft  = position;

    bool result = true;
    for (vector<pair<BamReader*, BamAlignment*> >::iterator it = readers.begin(); it != readers.end(); ++it) {
        BamReader* reader = it->first;
        result &= reader->Jump(refID, position);
        if (!result) {
            cerr << "ERROR: could not jump " << reader->GetFilename() << " to " << refID << ":" << position << endl;
            exit(1);
        }
    }
    if (result) UpdateAlignments();
    return result;
}

bool BamMultiReader::SetRegion(const int& leftRefID, const int& leftPosition, const int& rightRefID, const int& rightPosition) {

    BamRegion region(leftRefID, leftPosition, rightRefID, rightPosition);

    return SetRegion(region);

}

bool BamMultiReader::SetRegion(const BamRegion& region) {

    Region = region;

    // NB: While it may make sense to track readers in which we can
    // successfully SetRegion, In practice a failure of SetRegion means "no
    // alignments here."  It makes sense to simply accept the failure,
    // UpdateAlignments(), and continue.

    for (vector<pair<BamReader*, BamAlignment*> >::iterator it = readers.begin(); it != readers.end(); ++it) {
        it->first->SetRegion(region);
    }

    UpdateAlignments();

    return true;

}

void BamMultiReader::UpdateAlignments(void) {
    // Update Alignments
    alignments.clear();
    for (vector<pair<BamReader*, BamAlignment*> >::iterator it = readers.begin(); it != readers.end(); ++it) {
        BamReader* br = it->first;
        BamAlignment* ba = it->second;
        if (br->GetNextAlignment(*ba)) {
            alignments.insert(make_pair(make_pair(ba->RefID, ba->Position), 
                                        make_pair(br, ba)));
        } else {
            // assume BamReader end of region / EOF
        }
    }
}

// opens BAM files
bool BamMultiReader::Open(const vector<string> filenames, bool openIndexes, bool coreMode, bool useDefaultIndex) {
    
    // for filename in filenames
    fileNames = filenames; // save filenames in our multireader
    for (vector<string>::const_iterator it = filenames.begin(); it != filenames.end(); ++it) {
        string filename = *it;
        BamReader* reader = new BamReader;

        bool openedOK = true;
        if (openIndexes) {
            if (useDefaultIndex)
                openedOK = reader->Open(filename, filename + ".bai");
            else 
                openedOK = reader->Open(filename, filename + ".bti");
        } else {
            openedOK = reader->Open(filename); // for merging, jumping is disallowed
        }
        
        // if file opened ok, check that it can be read
        if ( openedOK ) {
           
            bool fileOK = true;
            BamAlignment* alignment = new BamAlignment;
            if (coreMode) {
                fileOK &= reader->GetNextAlignmentCore(*alignment);
            } else {
                fileOK &= reader->GetNextAlignment(*alignment);
            }
            
            if (fileOK) {
                readers.push_back(make_pair(reader, alignment)); // store pointers to our readers for cleanup
                alignments.insert(make_pair(make_pair(alignment->RefID, alignment->Position),
                                            make_pair(reader, alignment)));
            } else {
                cerr << "WARNING: could not read first alignment in " << filename << ", ignoring file" << endl;
                // if only file available & could not be read, return failure
                if ( filenames.size() == 1 ) return false;
            }
        
        } 
       
        // TODO; any more error handling on openedOK ??
        else 
            return false;
    }

    // files opened ok, at least one alignment could be read,
    // now need to check that all files use same reference data
    ValidateReaders();
    return true;
}

void BamMultiReader::PrintFilenames(void) {
    for (vector<pair<BamReader*, BamAlignment*> >::iterator it = readers.begin(); it != readers.end(); ++it) {
        BamReader* reader = it->first;
        cout << reader->GetFilename() << endl;
    }
}

// for debugging
void BamMultiReader::DumpAlignmentIndex(void) {
    for (AlignmentIndex::const_iterator it = alignments.begin(); it != alignments.end(); ++it) {
        cerr << it->first.first << ":" << it->first.second << " " << it->second.first->GetFilename() << endl;
    }
}

// returns BAM file pointers to beginning of alignment data
bool BamMultiReader::Rewind(void) { 
    bool result = true;
    for (vector<pair<BamReader*, BamAlignment*> >::iterator it = readers.begin(); it != readers.end(); ++it) {
        BamReader* reader = it->first;
        result &= reader->Rewind();
    }
    return result;
}

// saves index data to BAM index files (".bai") where necessary, returns success/fail
bool BamMultiReader::CreateIndexes(bool useDefaultIndex) {
    bool result = true;
    for (vector<pair<BamReader*, BamAlignment*> >::iterator it = readers.begin(); it != readers.end(); ++it) {
        BamReader* reader = it->first;
        result &= reader->CreateIndex(useDefaultIndex);
    }
    return result;
}

// makes a virtual, unified header for all the bam files in the multireader
const string BamMultiReader::GetHeaderText(void) const {

    string mergedHeader = "";
    map<string, bool> readGroups;

    // foreach extraction entry (each BAM file)
    for (vector<pair<BamReader*, BamAlignment*> >::const_iterator rs = readers.begin(); rs != readers.end(); ++rs) {

        map<string, bool> currentFileReadGroups;

        BamReader* reader = rs->first;

        stringstream header(reader->GetHeaderText());
        vector<string> lines;
        string item;
        while (getline(header, item))
            lines.push_back(item);

        for (vector<string>::const_iterator it = lines.begin(); it != lines.end(); ++it) {

            // get next line from header, skip if empty
            string headerLine = *it;
            if ( headerLine.empty() ) { continue; }

            // if first file, save HD & SQ entries
            if ( rs == readers.begin() ) {
                if ( headerLine.find("@HD") == 0 || headerLine.find("@SQ") == 0) {
                    mergedHeader.append(headerLine.c_str());
                    mergedHeader.append(1, '\n');
                }
            }

            // (for all files) append RG entries if they are unique
            if ( headerLine.find("@RG") == 0 ) {
                stringstream headerLineSs(headerLine);
                string part, readGroupPart, readGroup;
                while(std::getline(headerLineSs, part, '\t')) {
                    stringstream partSs(part);
                    string subtag;
                    std::getline(partSs, subtag, ':');
                    if (subtag == "ID") {
                        std::getline(partSs, readGroup, ':');
                        break;
                    }
                }
                if (readGroups.find(readGroup) == readGroups.end()) { // prevents duplicate @RG entries
                    mergedHeader.append(headerLine.c_str() );
                    mergedHeader.append(1, '\n');
                    readGroups[readGroup] = true;
                    currentFileReadGroups[readGroup] = true;
                } else {
                    // warn iff we are reading one file and discover duplicated @RG tags in the header
                    // otherwise, we emit no warning, as we might be merging multiple BAM files with identical @RG tags
                    if (currentFileReadGroups.find(readGroup) != currentFileReadGroups.end()) {
                        cerr << "WARNING: duplicate @RG tag " << readGroup 
                            << " entry in header of " << reader->GetFilename() << endl;
                    }
                }
            }
        }
    }

    // return merged header text
    return mergedHeader;
}

// ValidateReaders checks that all the readers point to BAM files representing
// alignments against the same set of reference sequences, and that the
// sequences are identically ordered.  If these checks fail the operation of
// the multireader is undefined, so we force program exit.
void BamMultiReader::ValidateReaders(void) const {
    int firstRefCount = readers.front().first->GetReferenceCount();
    BamTools::RefVector firstRefData = readers.front().first->GetReferenceData();
    for (vector<pair<BamReader*, BamAlignment*> >::const_iterator it = readers.begin(); it != readers.end(); ++it) {
        BamReader* reader = it->first;
        BamTools::RefVector currentRefData = reader->GetReferenceData();
        BamTools::RefVector::const_iterator f = firstRefData.begin();
        BamTools::RefVector::const_iterator c = currentRefData.begin();
        if (reader->GetReferenceCount() != firstRefCount || firstRefData.size() != currentRefData.size()) {
            cerr << "ERROR: mismatched number of references in " << reader->GetFilename()
                      << " expected " << firstRefCount 
                      << " reference sequences but only found " << reader->GetReferenceCount() << endl;
            exit(1);
        }
        // this will be ok; we just checked above that we have identically-sized sets of references
        // here we simply check if they are all, in fact, equal in content
        while (f != firstRefData.end()) {
            if (f->RefName != c->RefName || f->RefLength != c->RefLength) {
                cerr << "ERROR: mismatched references found in " << reader->GetFilename()
                          << " expected: " << endl;
                for (BamTools::RefVector::const_iterator a = firstRefData.begin(); a != firstRefData.end(); ++a)
                    cerr << a->RefName << " " << a->RefLength << endl;
                cerr << "but found: " << endl;
                for (BamTools::RefVector::const_iterator a = currentRefData.begin(); a != currentRefData.end(); ++a)
                    cerr << a->RefName << " " << a->RefLength << endl;
                exit(1);
            }
            ++f; ++c;
        }
    }
}

// NB: The following functions assume that we have identical references for all
// BAM files.  We enforce this by invoking the above validation function
// (ValidateReaders) to verify that our reference data is the same across all
// files on Open, so we will not encounter a situation in which there is a
// mismatch and we are still live.

// returns the number of reference sequences
const int BamMultiReader::GetReferenceCount(void) const {
    return readers.front().first->GetReferenceCount();
}

// returns vector of reference objects
const BamTools::RefVector BamMultiReader::GetReferenceData(void) const {
    return readers.front().first->GetReferenceData();
}

const int BamMultiReader::GetReferenceID(const string& refName) const { 
    return readers.front().first->GetReferenceID(refName);
}
