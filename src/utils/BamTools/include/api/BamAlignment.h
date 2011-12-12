// ***************************************************************************
// BamAlignment.h (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 22 April 2011 (DB)
// ---------------------------------------------------------------------------
// Provides the BamAlignment data structure
// ***************************************************************************

#ifndef BAMALIGNMENT_H
#define BAMALIGNMENT_H

#include <api/api_global.h>
#include <api/BamAux.h>
#include <string>
#include <vector>

namespace BamTools {

// forward declaration of BamAlignment's friend classes
namespace Internal {
    class BamReaderPrivate;
    class BamWriterPrivate;
} // namespace Internal

// BamAlignment data structure
struct API_EXPORT BamAlignment {

    // constructors & destructor
    public:
        BamAlignment(void);
        BamAlignment(const BamAlignment& other);
        ~BamAlignment(void);

    // queries against alignment flags
    public:        
        bool IsDuplicate(void) const;           // returns true if this read is a PCR duplicate
        bool IsFailedQC(void) const;            // returns true if this read failed quality control
        bool IsFirstMate(void) const;           // returns true if alignment is first mate on read
        bool IsMapped(void) const;              // returns true if alignment is mapped
        bool IsMateMapped(void) const;          // returns true if alignment's mate is mapped
        bool IsMateReverseStrand(void) const;   // returns true if alignment's mate mapped to reverse strand
        bool IsPaired(void) const;              // returns true if alignment part of paired-end read
        bool IsPrimaryAlignment(void) const;    // returns true if reported position is primary alignment
        bool IsProperPair(void) const;          // returns true if alignment is part of read that satisfied paired-end resolution
        bool IsReverseStrand(void) const;       // returns true if alignment mapped to reverse strand
        bool IsSecondMate(void) const;          // returns true if alignment is second mate on read

    // manipulate alignment flags
    public:        
        void SetIsDuplicate(bool ok);           // sets value of "PCR duplicate" flag
        void SetIsFailedQC(bool ok);            // sets value of "failed quality control" flag
        void SetIsFirstMate(bool ok);           // sets value of "alignment is first mate" flag
        void SetIsMapped(bool ok);              // sets value of "alignment is mapped" flag
        void SetIsMateMapped(bool ok);          // sets value of "alignment's mate is mapped" flag
        void SetIsMateReverseStrand(bool ok);   // sets value of "alignment's mate mapped to reverse strand" flag
        void SetIsPaired(bool ok);              // sets value of "alignment part of paired-end read" flag
        void SetIsPrimaryAlignment(bool ok);    // sets value of "position is primary alignment" flag
        void SetIsProperPair(bool ok);          // sets value of "alignment is part of read that satisfied paired-end resolution" flag
        void SetIsReverseStrand(bool ok);       // sets value of "alignment mapped to reverse strand" flag
        void SetIsSecondMate(bool ok);          // sets value of "alignment is second mate on read" flag

        // legacy methods (consider deprecated, but still available)
        void SetIsMateUnmapped(bool ok);        // complement of using SetIsMateMapped()
        void SetIsSecondaryAlignment(bool ok);  // complement of using SetIsPrimaryAlignment()
        void SetIsUnmapped(bool ok);            // complement of using SetIsMapped()

    // tag data access methods
    public:

        // -------------------------------------------------------------------------------------
        // N.B. - The following tag access methods may not be used on BamAlignments fetched
        // using BamReader::GetNextAlignmentCore().  Attempting to use them will not result in 
        // error message (to keep output clean) but will ALWAYS return false.  Only user-created
        // BamAlignments or those retrieved using BamReader::GetNextAlignment() are valid here.
        //
        // You can call BuildCharData() on such an alignment retrieved by GetNextAlignmentCore().
        // This populates all the character data, and will enable subsequent queries on tag data.
        // -------------------------------------------------------------------------------------

        // adds a tag
        bool AddTag(const std::string& tag, const std::string& type, const std::string& value);
        bool AddTag(const std::string& tag, const std::string& type, const uint32_t& value);
        bool AddTag(const std::string& tag, const std::string& type, const int32_t& value);
        bool AddTag(const std::string& tag, const std::string& type, const float& value);

        // adds a "binary array" tag
        bool AddTag(const std::string& tag, const std::vector<uint8_t>& values);
        bool AddTag(const std::string& tag, const std::vector<int8_t>& values);
        bool AddTag(const std::string& tag, const std::vector<uint16_t>& values);
        bool AddTag(const std::string& tag, const std::vector<int16_t>& values);
        bool AddTag(const std::string& tag, const std::vector<uint32_t>& values);
        bool AddTag(const std::string& tag, const std::vector<int32_t>& values);
        bool AddTag(const std::string& tag, const std::vector<float>& values);

        // edits a tag
        bool EditTag(const std::string& tag, const std::string& type, const std::string& value);
        bool EditTag(const std::string& tag, const std::string& type, const uint32_t& value);
        bool EditTag(const std::string& tag, const std::string& type, const int32_t& value);
        bool EditTag(const std::string& tag, const std::string& type, const float& value);

        // edits a "binary array" tag
        bool EditTag(const std::string& tag, const std::vector<uint8_t>& values);
        bool EditTag(const std::string& tag, const std::vector<int8_t>& values);
        bool EditTag(const std::string& tag, const std::vector<uint16_t>& values);
        bool EditTag(const std::string& tag, const std::vector<int16_t>& values);
        bool EditTag(const std::string& tag, const std::vector<uint32_t>& values);
        bool EditTag(const std::string& tag, const std::vector<int32_t>& values);
        bool EditTag(const std::string& tag, const std::vector<float>& values);

        // retrieves data for a tag
        bool GetTag(const std::string& tag, std::string& destination) const;
        bool GetTag(const std::string& tag, uint32_t& destination) const;
        bool GetTag(const std::string& tag, int32_t& destination) const;
        bool GetTag(const std::string& tag, float& destination) const;

        // retrieves data for a "binary array" tag
        bool GetTag(const std::string& tag, std::vector<uint32_t>& destination) const;
        bool GetTag(const std::string& tag, std::vector<int32_t>& destination) const;
        bool GetTag(const std::string& tag, std::vector<float>& destination) const;

        // retrieves the BAM tag-type character for a tag
        bool GetTagType(const std::string& tag, char& type) const;

        // legacy methods (consider deprecated, but still available)
        bool GetEditDistance(uint32_t& editDistance) const;         // retrieves value of "NM" tag
        bool GetReadGroup(std::string& readGroup) const;            // retrieves value of "RG" tag
        
        // returns true if alignment has a record for this tag name
        bool HasTag(const std::string& tag) const;

        // removes a tag
        bool RemoveTag(const std::string& tag);

    // additional methods
    public:
        // populates alignment string fields
        bool BuildCharData(void);
        // calculates alignment end position
        int GetEndPosition(bool usePadded = false, bool zeroBased = true) const;  

    // public data fields
    public:
        std::string Name;               // read name
        int32_t     Length;             // length of query sequence
        std::string QueryBases;         // 'original' sequence (as reported from sequencing machine)
        std::string AlignedBases;       // 'aligned' sequence (includes any indels, padding, clipping)
        std::string Qualities;          // FASTQ qualities (ASCII characters, not numeric values)
        std::string TagData;            // tag data (use provided methods to query/modify)
        int32_t     RefID;              // ID number for reference sequence
        int32_t     Position;           // position (0-based) where alignment starts
        uint16_t    Bin;                // BAM (standard) index bin number for this alignment
        uint16_t    MapQuality;         // mapping quality score
        uint32_t    AlignmentFlag;      // alignment bit-flag (use provided methods to query/modify)
        std::vector<CigarOp> CigarData; // CIGAR operations for this alignment
        int32_t     MateRefID;          // ID number for reference sequence where alignment's mate was aligned
        int32_t     MatePosition;       // position (0-based) where alignment's mate starts
        int32_t     InsertSize;         // mate-pair insert size
        std::string Filename;           // name of BAM file which this alignment comes from

    //! \cond
    // internal utility methods
    private:
        bool FindTag(const std::string& tag,
                     char*& pTagData,
                     const unsigned int& tagDataLength,
                     unsigned int& numBytesParsed) const;
        bool IsValidSize(const std::string& tag,
                         const std::string& type) const;
        bool SkipToNextTag(const char storageType,
                           char*& pTagData,
                           unsigned int& numBytesParsed) const;

    // internal data
    private:

        struct BamAlignmentSupportData {
      
            // data members
            std::string AllCharData;
            uint32_t    BlockLength;
            uint32_t    NumCigarOperations;
            uint32_t    QueryNameLength;
            uint32_t    QuerySequenceLength;
            bool        HasCoreOnly;
            
            // constructor
            BamAlignmentSupportData(void)
                : BlockLength(0)
                , NumCigarOperations(0)
                , QueryNameLength(0)
                , QuerySequenceLength(0)
                , HasCoreOnly(false)
            { }
        };
        BamAlignmentSupportData SupportData;
        friend class Internal::BamReaderPrivate;
        friend class Internal::BamWriterPrivate;
    //! \endcond
};

typedef std::vector<BamAlignment> BamAlignmentVector;

} // namespace BamTools

#endif // BAMALIGNMENT_H
