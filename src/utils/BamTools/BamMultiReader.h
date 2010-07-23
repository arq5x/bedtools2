// ***************************************************************************
// BamMultiReader.h (c) 2010 Erik Garrison
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 20 July 2010 (DB)
// ---------------------------------------------------------------------------
// Functionality for simultaneously reading multiple BAM files
// ***************************************************************************

#ifndef BAMMULTIREADER_H
#define BAMMULTIREADER_H

// C++ includes
#include <string>
#include <map>
#include <utility> // for pair
#include <sstream>

using namespace std;

// BamTools includes
#include "BamAux.h"
#include "BamReader.h"

namespace BamTools {

// index mapping reference/position pairings to bamreaders and their alignments
typedef multimap<pair<int, int>, pair<BamReader*, BamAlignment*> > AlignmentIndex;


class BamMultiReader {

    // constructor / destructor
    public:
        BamMultiReader(void);
        ~BamMultiReader(void);

    // public interface
    public:

        // positioning
        int CurrentRefID;
        int CurrentLeft;

        // region under analysis, specified using SetRegion
        BamRegion Region;

        // ----------------------
        // BAM file operations
        // ----------------------

        // close BAM files
        void Close(void);

        // opens BAM files (and optional BAM index files, if provided)
        // @openIndexes - triggers index opening, useful for suppressing
        // error messages during merging of files in which we may not have
        // indexes.
        // @coreMode - setup our first alignments using GetNextAlignmentCore();
        // also useful for merging
        bool Open(const vector<string> filenames, bool openIndexes = true, bool coreMode = false, bool useDefaultIndex = true);

        // performs random-access jump to reference, position
        bool Jump(int refID, int position = 0);

        // sets the target region
        bool SetRegion(const BamRegion& region);
        bool SetRegion(const int&, const int&, const int&, const int&); // convenience function to above

        // returns file pointers to beginning of alignments
        bool Rewind(void);

        // ----------------------
        // access alignment data
        // ----------------------
        // updates the reference id marker to match the lower limit of our readers
        void UpdateReferenceID(void);

        // retrieves next available alignment (returns success/fail) from all files
        bool GetNextAlignment(BamAlignment&);
        // retrieves next available alignment (returns success/fail) from all files
        // and populates the support data with information about the alignment
        // *** BUT DOES NOT PARSE CHARACTER DATA FROM THE ALIGNMENT
        bool GetNextAlignmentCore(BamAlignment&);
        // ... should this be private?
        bool HasOpenReaders(void);

        // ----------------------
        // access auxiliary data
        // ----------------------

        // returns unified SAM header text for all files
        const string GetHeaderText(void) const;
        // returns number of reference sequences
        const int GetReferenceCount(void) const;
        // returns vector of reference objects
        const BamTools::RefVector GetReferenceData(void) const;
        // returns reference id (used for BamMultiReader::Jump()) for the given reference name
        const int GetReferenceID(const std::string& refName) const;
        // validates that we have a congruent set of BAM files that are aligned against the same reference sequences
        void ValidateReaders() const;

        // ----------------------
        // BAM index operations
        // ----------------------

        // creates index for BAM files which lack them, saves to files (default = bamFilename + ".bai")
        bool CreateIndexes(bool useDefaultIndex = true);

        //const int GetReferenceID(const string& refName) const;

        // utility
        void PrintFilenames(void);
        void DumpAlignmentIndex(void);
        void UpdateAlignments(void); // updates our alignment cache

    // private implementation
    private:

        // the set of readers and alignments which we operate on, maintained throughout the life of this class
        vector<pair<BamReader*, BamAlignment*> > readers;

        // readers and alignments sorted by reference id and position, to keep track of the lowest (next) alignment
        // when a reader reaches EOF, its entry is removed from this index
        AlignmentIndex alignments;

        vector<string> fileNames;
};

} // namespace BamTools

#endif // BAMMULTIREADER_H
