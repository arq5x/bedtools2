// ***************************************************************************
// BamAux.h (c) 2009 Derek Barnett, Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 8 December 2009 (DB)
// ---------------------------------------------------------------------------
// Provides the basic constants, data structures, etc. for using BAM files
// ***************************************************************************

#ifndef BAMAUX_H
#define BAMAUX_H

// C inclues
#include <cstdlib>
#include <cstring>

// C++ includes
#include <exception>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace BamTools {

// BAM constants
const int BAM_CORE_SIZE   = 32;
const int BAM_CMATCH      = 0;
const int BAM_CINS        = 1;
const int BAM_CDEL        = 2;
const int BAM_CREF_SKIP   = 3;
const int BAM_CSOFT_CLIP  = 4;
const int BAM_CHARD_CLIP  = 5;
const int BAM_CPAD        = 6;
const int BAM_CIGAR_SHIFT = 4;
const int BAM_CIGAR_MASK  = ((1 << BAM_CIGAR_SHIFT) - 1);

// BAM index constants
const int MAX_BIN           = 37450;	// =(8^6-1)/7+1
const int BAM_MIN_CHUNK_GAP = 32768;
const int BAM_LIDX_SHIFT    = 14;

// Explicit variable sizes
const int BT_SIZEOF_INT = 4;

struct CigarOp;

struct BamAlignment {

    // Queries against alignment flag
    public:
        // Returns true if this read is a PCR duplicate (determined by external app)
        bool IsDuplicate(void) const { return ( (AlignmentFlag & DUPLICATE) != 0 ); }
        // Returns true if this read failed quality control (determined by external app)
        bool IsFailedQC(void) const { return ( (AlignmentFlag & QC_FAILED) != 0 ); }
        // Returns true if alignment is first mate on read
        bool IsFirstMate(void) const { return ( (AlignmentFlag & READ_1) != 0 ); }
        // Returns true if alignment is mapped
        bool IsMapped(void) const { return ( (AlignmentFlag & UNMAPPED) == 0 ); }
        // Returns true if alignment's mate is mapped
        bool IsMateMapped(void) const { return ( (AlignmentFlag & MATE_UNMAPPED) == 0 ); }
        // Returns true if alignment's mate mapped to reverse strand
        bool IsMateReverseStrand(void) const { return ( (AlignmentFlag & MATE_REVERSE)  != 0 ); }
        // Returns true if alignment part of paired-end read
        bool IsPaired(void) const { return ( (AlignmentFlag & PAIRED) != 0 ); }
        // Returns true if this position is primary alignment (determined by external app)
        bool IsPrimaryAlignment(void) const  { return ( (AlignmentFlag & SECONDARY) == 0 ); }
        // Returns true if alignment is part of read that satisfied paired-end resolution (determined by external app)
        bool IsProperPair(void) const { return ( (AlignmentFlag & PROPER_PAIR) != 0 ); }
        // Returns true if alignment mapped to reverse strand
        bool IsReverseStrand(void) const { return ( (AlignmentFlag & REVERSE) != 0 ); }
        // Returns true if alignment is second mate on read
        bool IsSecondMate(void) const { return ( (AlignmentFlag & READ_2) != 0 ); }


		/*
		  Aaron's Additions
		*/
		
		
		// Sets "PCR duplicate" bit to TRUE 
        void SetAsDuplicate(void) { AlignmentFlag |= DUPLICATE; }
        // Sets "failed quality control" bit to TRUE
        void SetAsFailedQC(void) { AlignmentFlag |= QC_FAILED; }
        // Sets "alignment is first mate" bit to TRUE
        void SetAsFirstMate(void) { AlignmentFlag |= READ_1; }
        // Sets "alignment is mapped" bit to TRUE
        void SetAsUnMapped(void) { AlignmentFlag |= UNMAPPED; }
        // Sets "alignment's mate is mapped" bit to TRUE
        void SetAsMateUnMapped(void) { AlignmentFlag |= MATE_UNMAPPED; }
        // Sets "alignment's mate mapped to reverse strand" bit to TRUE
        void SetAsMateReverseStrand(void) { AlignmentFlag |= MATE_REVERSE; }
        // Sets "alignment part of paired-end read" bit to TRUE
        void SetAsPaired(void) { AlignmentFlag |= PAIRED; }
        // Sets "position is primary alignment (determined by external app)" to true.
        void SetAsSecondaryAlignment(void)  { AlignmentFlag |= SECONDARY; }
        // Sets "alignment is part of read that satisfied paired-end resolution" bit to TRUE
        void SetAsProperPair(void) { AlignmentFlag |= PROPER_PAIR; }
        // Sets "alignment mapped to reverse strand" bit to TRUE
        void SetAsReverseStrand(void) { AlignmentFlag |= REVERSE; }
        // Sets "alignment is second mate on read" bit to TRUE
        void SetAsSecondMate(void) { AlignmentFlag |= READ_2; }

		
		void SetInsertSize(int32_t size) { InsertSize = size; }
		/*
		  END: Aaron's Additions
		*/

    public:

        // get "RG" tag data
        bool GetReadGroup(std::string& readGroup) const {

            if ( TagData.empty() ) { return false; }

            // localize the tag data
            char* pTagData = (char*)TagData.data();
            const unsigned int tagDataLen = TagData.size();
            unsigned int numBytesParsed = 0;

            bool foundReadGroupTag = false;
            while( numBytesParsed < tagDataLen ) {

                const char* pTagType = pTagData;
                const char* pTagStorageType = pTagData + 2;
                pTagData       += 3;
                numBytesParsed += 3;

                // check the current tag
                if ( std::strncmp(pTagType, "RG", 2) == 0 ) {
                    foundReadGroupTag = true;
                    break;
                }

                // get the storage class and find the next tag
                SkipToNextTag( *pTagStorageType, pTagData, numBytesParsed );
            }

            // return if the read group tag was not present
            if ( !foundReadGroupTag ) { return false; }

            // assign the read group
            const unsigned int readGroupLen = std::strlen(pTagData);
            readGroup.resize(readGroupLen);
            std::memcpy( (char*)readGroup.data(), pTagData, readGroupLen );
            return true;
        }

        // get "NM" tag data - contributed by Aaron Quinlan
        bool GetEditDistance(uint8_t& editDistance) const {

            if ( TagData.empty() ) { return false; }

            // localize the tag data
            char* pTagData = (char*)TagData.data();
            const unsigned int tagDataLen = TagData.size();
            unsigned int numBytesParsed = 0;

            bool foundEditDistanceTag = false;
            while( numBytesParsed < tagDataLen ) {

                const char* pTagType = pTagData;
                const char* pTagStorageType = pTagData + 2;
                pTagData       += 3;
                numBytesParsed += 3;

                // check the current tag
                if ( strncmp(pTagType, "NM", 2) == 0 ) {
                    foundEditDistanceTag = true;
                    break;
                }

                // get the storage class and find the next tag
                SkipToNextTag( *pTagStorageType, pTagData, numBytesParsed );
            }
            // return if the edit distance tag was not present
            if ( !foundEditDistanceTag ) { return false; }

            // assign the editDistance value
            memcpy(&editDistance, pTagData, 1);
            return true;
        }

    private:
        static void SkipToNextTag(const char storageType, char* &pTagData, unsigned int& numBytesParsed) {
            switch(storageType) {

                case 'A':
                case 'c':
                case 'C':
                        ++numBytesParsed;
                        ++pTagData;
                        break;

                case 's':
                case 'S':
                case 'f':
                        numBytesParsed += 2;
                        pTagData       += 2;
                        break;

                case 'i':
                case 'I':
                        numBytesParsed += 4;
                        pTagData       += 4;
                        break;

                case 'Z':
                case 'H':
                        while(*pTagData) {
                            ++numBytesParsed;
                            ++pTagData;
                        }
                        break;

                default:
                        printf("ERROR: Unknown tag storage class encountered: [%c]\n", *pTagData);
                        exit(1);
            }
        }

    // Data members
    public:
        std::string  Name;              // Read name
        int32_t      Length;            // Query length
        std::string  QueryBases;        // 'Original' sequence (as reported from sequencing machine)
        std::string  AlignedBases;      // 'Aligned' sequence (includes any indels, padding, clipping)
        std::string  Qualities;         // FASTQ qualities (ASCII characters, not numeric values)
        std::string  TagData;           // Tag data (accessor methods will pull the requested information out)
        int32_t      RefID;             // ID number for reference sequence
        int32_t      Position;          // Position (0-based) where alignment starts
        uint16_t     Bin;               // Bin in BAM file where this alignment resides
        uint16_t     MapQuality;        // Mapping quality score
        uint32_t     AlignmentFlag;     // Alignment bit-flag - see Is<something>() methods for available queries
        std::vector<CigarOp> CigarData; // CIGAR operations for this alignment
        int32_t      MateRefID;         // ID number for reference sequence where alignment's mate was aligned
        int32_t      MatePosition;      // Position (0-based) where alignment's mate starts
        int32_t      InsertSize;        // Mate-pair insert size

    // Alignment flag query constants
    private:
        enum { PAIRED        = 1,
               PROPER_PAIR   = 2,
               UNMAPPED      = 4,
               MATE_UNMAPPED = 8,
               REVERSE       = 16,
               MATE_REVERSE  = 32,
               READ_1        = 64,
               READ_2        = 128,
               SECONDARY     = 256,
               QC_FAILED     = 512,
               DUPLICATE     = 1024
             };
};

// ----------------------------------------------------------------
// Auxiliary data structs & typedefs

struct CigarOp {
    char     Type;   // Operation type (MIDNSHP)
    uint32_t Length; // Operation length (number of bases)
};

struct RefData {
    // data members
    std::string RefName;          // Name of reference sequence
    int         RefLength;        // Length of reference sequence
    bool        RefHasAlignments; // True if BAM file contains alignments mapped to reference sequence
    // constructor
    RefData(void)
        : RefLength(0)
        , RefHasAlignments(false)
    { }
};

typedef std::vector<RefData> RefVector;
typedef std::vector<BamAlignment> BamAlignmentVector;

// ----------------------------------------------------------------
// Indexing structs & typedefs

struct Chunk {
    // data members
    uint64_t Start;
    uint64_t Stop;
    // constructor
    Chunk(const uint64_t& start = 0, const uint64_t& stop = 0)
        : Start(start)
        , Stop(stop)
    { }
};

inline
bool ChunkLessThan(const Chunk& lhs, const Chunk& rhs) {
    return lhs.Start < rhs.Start;
}

typedef std::vector<Chunk> ChunkVector;
typedef std::map<uint32_t, ChunkVector> BamBinMap;
typedef std::vector<uint64_t> LinearOffsetVector;

struct ReferenceIndex {
    // data members
    BamBinMap Bins;
    LinearOffsetVector Offsets;
    // constructor
    ReferenceIndex(const BamBinMap& binMap = BamBinMap(),
                   const LinearOffsetVector& offsets = LinearOffsetVector())
        : Bins(binMap)
        , Offsets(offsets)
    { }
};

typedef std::vector<ReferenceIndex> BamIndex;

} // namespace BamTools

#endif // BAMAUX_H
