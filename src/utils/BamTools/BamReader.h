// ***************************************************************************
// BamReader.h (c) 2009 Derek Barnett, Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 9 July 2010 (DB)
// ---------------------------------------------------------------------------
// Uses BGZF routines were adapted from the bgzf.c code developed at the Broad
// Institute.
// ---------------------------------------------------------------------------
// Provides the basic functionality for reading BAM files
// ***************************************************************************

#ifndef BAMREADER_H
#define BAMREADER_H

// C++ includes
#include <string>

// BamTools includes
#include "BamAux.h"

namespace BamTools {
  
class BamReader {

    // constructor / destructor
    public:
        BamReader(void);
        ~BamReader(void);

    // public interface
    public:

        // ----------------------
        // BAM file operations
        // ----------------------

        // close BAM file
        void Close(void);
        // returns whether reader is open for reading or not
        bool IsOpen(void) const;
        // performs random-access jump to reference, position
        bool Jump(int refID, int position = 0);
        // opens BAM file (and optional BAM index file, if provided)
        bool Open(const std::string& filename, const std::string& indexFilename = "");
        // returns file pointer to beginning of alignments
        bool Rewind(void);
        // sets a region of interest (with left & right bound reference/position)
        // attempts a Jump() to left bound as well
        // returns success/failure of Jump()
        bool SetRegion(const BamRegion& region);
        bool SetRegion(const int& leftRefID, const int& leftBound, const int& rightRefID, const int& rightBound);

        // ----------------------
        // access alignment data
        // ----------------------

        // retrieves next available alignment (returns success/fail)
        bool GetNextAlignment(BamAlignment& bAlignment);
        
        // retrieves next available alignment core data (returns success/fail)
        // ** DOES NOT parse any character data (read name, bases, qualities, tag data)
        //    these can be accessed, if necessary, from the supportData 
        // useful for operations requiring ONLY positional or other alignment-related information
        bool GetNextAlignmentCore(BamAlignment& bAlignment);

        // ----------------------
        // access auxiliary data
        // ----------------------

        // returns SAM header text
        const std::string GetHeaderText(void) const;
        // returns number of reference sequences
        int GetReferenceCount(void) const;
        // returns vector of reference objects
        const BamTools::RefVector& GetReferenceData(void) const;
        // returns reference id (used for BamReader::Jump()) for the given reference name
        int GetReferenceID(const std::string& refName) const;
        // returns the name of the file associated with this BamReader
        const std::string GetFilename(void) const;

        // ----------------------
        // BAM index operations
        // ----------------------

        // creates index for BAM file, saves to file (default = bamFilename + ".bai")
        bool CreateIndex(bool useDefaultIndex = true);
        
    // private implementation
    private:
        struct BamReaderPrivate;
        BamReaderPrivate* d;
};

} // namespace BamTools

#endif // BAMREADER_H
