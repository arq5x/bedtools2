// ***************************************************************************
// BamAux.h (c) 2009 Derek Barnett, Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 27 July 2010 (DB)
// ---------------------------------------------------------------------------
// Provides the basic constants, data structures, etc. for using BAM files
// ***************************************************************************

#ifndef BAMAUX_H
#define BAMAUX_H

// C inclues
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>

// C++ includes
#include <exception>
#include <map>
#include <string>
#include <utility>
#include <vector>

// Platform-specific type definitions
#ifndef BAMTOOLS_TYPES
#define BAMTOOLS_TYPES
    #ifdef _MSC_VER
        typedef char                 int8_t;
        typedef unsigned char       uint8_t;
        typedef short               int16_t;
        typedef unsigned short     uint16_t;
        typedef int                 int32_t;
        typedef unsigned int       uint32_t;
        typedef long long           int64_t;
        typedef unsigned long long uint64_t;
    #else
        #include <stdint.h>
    #endif
#endif // BAMTOOLS_TYPES

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

    // constructors & destructor
    public:
        BamAlignment(void);
        BamAlignment(const BamAlignment& other);
        ~BamAlignment(void);

    // Queries against alignment flags
    public:        
        bool IsDuplicate(void) const;		// Returns true if this read is a PCR duplicate       
        bool IsFailedQC(void) const;		// Returns true if this read failed quality control      
        bool IsFirstMate(void) const;		// Returns true if alignment is first mate on read        
        bool IsMapped(void) const;		// Returns true if alignment is mapped        
        bool IsMateMapped(void) const;		// Returns true if alignment's mate is mapped        
        bool IsMateReverseStrand(void) const;	// Returns true if alignment's mate mapped to reverse strand        
        bool IsPaired(void) const;		// Returns true if alignment part of paired-end read        
        bool IsPrimaryAlignment(void) const;	// Returns true if reported position is primary alignment       
        bool IsProperPair(void) const;		// Returns true if alignment is part of read that satisfied paired-end resolution     
        bool IsReverseStrand(void) const;	// Returns true if alignment mapped to reverse strand
        bool IsSecondMate(void) const;		// Returns true if alignment is second mate on read

    // Manipulate alignment flags
    public:        
        void SetIsDuplicate(bool ok);		// Sets "PCR duplicate" flag        
        void SetIsFailedQC(bool ok);		// Sets "failed quality control" flag        
        void SetIsFirstMate(bool ok);		// Sets "alignment is first mate" flag        
        void SetIsMateUnmapped(bool ok);	// Sets "alignment's mate is mapped" flag        
        void SetIsMateReverseStrand(bool ok);	// Sets "alignment's mate mapped to reverse strand" flag        
        void SetIsPaired(bool ok);		// Sets "alignment part of paired-end read" flag        
	void SetIsProperPair(bool ok);		// Sets "alignment is part of read that satisfied paired-end resolution" flag        
        void SetIsReverseStrand(bool ok);	// Sets "alignment mapped to reverse strand" flag        
        void SetIsSecondaryAlignment(bool ok);	// Sets "position is primary alignment" flag        
        void SetIsSecondMate(bool ok);		// Sets "alignment is second mate on read" flag        
        void SetIsUnmapped(bool ok);		// Sets "alignment is mapped" flag

    // Tag data access methods
    public:
        // -------------------------------------------------------------------------------------
        // N.B. - The following tag-modifying methods may not be used on BamAlignments fetched
        // using BamReader::GetNextAlignmentCore().  Attempting to use them will not result in 
        // error message (to keep output clean) but will ALWAYS return false.  Only user-
        // generated BamAlignments or those retrieved using BamReader::GetNextAlignment() are valid.

        // add tag data (create new TAG entry with TYPE and VALUE)
        // TYPE is one of {A, i, f, Z, H} depending on VALUE - see SAM/BAM spec for details
        // returns true if new data added, false if error or TAG already exists
        // N.B. - will NOT modify existing tag. Use EditTag() instead
        bool AddTag(const std::string& tag, const std::string& type, const std::string& value); // type must be Z or H
        bool AddTag(const std::string& tag, const std::string& type, const uint32_t& value);    // type must be A or i
        bool AddTag(const std::string& tag, const std::string& type, const int32_t& value);     // type must be A or i
        bool AddTag(const std::string& tag, const std::string& type, const float& value);       // type must be A, i, or f
        
        // edit tag data (sets existing TAG with TYPE to VALUE or adds new TAG if not already present)
        // TYPE is one of {A, i, f, Z, H} depending on VALUE - see SAM/BAM spec for details
        // returns true if edit was successfaul, false if error
        bool EditTag(const std::string& tag, const std::string& type, const std::string& value); // type must be Z or H
        bool EditTag(const std::string& tag, const std::string& type, const uint32_t& value);    // type must be A or i
        bool EditTag(const std::string& tag, const std::string& type, const int32_t& value);     // type must be A or i
        bool EditTag(const std::string& tag, const std::string& type, const float& value);       // type must be A, i, or f

        // specific tag data access methods - these only remain for legacy support
        bool GetEditDistance(uint32_t& editDistance) const; // get "NM" tag data (implemented as GetTag("NM", editDistance))
        bool GetReadGroup(std::string& readGroup) const;    // get "RG" tag data (implemented as GetTag("RG", readGroup)) 
        
        // generic tag data access methods 
        bool GetTag(const std::string& tag, std::string& destination) const;    // access variable-length char or hex strings 
        bool GetTag(const std::string& tag, uint32_t& destination) const;       // access unsigned integer data
        bool GetTag(const std::string& tag, int32_t& destination) const;        // access signed integer data
        bool GetTag(const std::string& tag, float& destination) const;          // access floating point data
        
        // remove tag data
        // returns true if removal was successful, false if error
        // N.B. - returns false if TAG does not exist (no removal can occur)
        bool RemoveTag(const std::string& tag);

    // Additional data access methods
    public:
	int GetEndPosition(bool usePadded = false) const;       // calculates alignment end position, based on starting position and CIGAR operations

    // 'internal' utility methods 
    private:
        static bool FindTag(const std::string& tag, char* &pTagData, const unsigned int& tagDataLength, unsigned int& numBytesParsed);
        static bool SkipToNextTag(const char storageType, char* &pTagData, unsigned int& numBytesParsed);

    // Data members
    public:
        std::string Name;              // Read name
        int32_t     Length;            // Query length
        std::string QueryBases;        // 'Original' sequence (as reported from sequencing machine)
        std::string AlignedBases;      // 'Aligned' sequence (includes any indels, padding, clipping)
        std::string Qualities;         // FASTQ qualities (ASCII characters, not numeric values)
        std::string TagData;           // Tag data (accessor methods will pull the requested information out)
        int32_t     RefID;             // ID number for reference sequence
        int32_t     Position;          // Position (0-based) where alignment starts
        uint16_t    Bin;               // Bin in BAM file where this alignment resides
        uint16_t    MapQuality;        // Mapping quality score
        uint32_t    AlignmentFlag;     // Alignment bit-flag - see Is<something>() methods to query this value, SetIs<something>() methods to manipulate 
        std::vector<CigarOp> CigarData; // CIGAR operations for this alignment
        int32_t     MateRefID;         // ID number for reference sequence where alignment's mate was aligned
        int32_t     MatePosition;      // Position (0-based) where alignment's mate starts
        int32_t     InsertSize;        // Mate-pair insert size
          
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
        
        // contains raw character data & lengths
        BamAlignmentSupportData SupportData;   
        
        // allow these classes access to BamAlignment private members (SupportData)
        // but client code should not need to touch this data
        friend class BamReader;
        friend class BamWriter;

    // Alignment flag query constants
    // Use the get/set methods above instead
    private:
        enum { PAIRED        = 1
             , PROPER_PAIR   = 2
             , UNMAPPED      = 4
             , MATE_UNMAPPED = 8
             , REVERSE       = 16
             , MATE_REVERSE  = 32
             , READ_1        = 64
             , READ_2        = 128
             , SECONDARY     = 256
             , QC_FAILED     = 512
             , DUPLICATE     = 1024 
             };
};

// ----------------------------------------------------------------
// Auxiliary data structs & typedefs

struct CigarOp {
  
    // data members
    char     Type;   // Operation type (MIDNSHP)
    uint32_t Length; // Operation length (number of bases)
    
    // constructor
    CigarOp(const char type = '\0', 
            const uint32_t length = 0) 
        : Type(type)
        , Length(length) 
    { }
};

struct RefData {
   
    // data members
    std::string RefName;          // Name of reference sequence
    int32_t     RefLength;        // Length of reference sequence
    bool        RefHasAlignments; // True if BAM file contains alignments mapped to reference sequence
    
    // constructor
    RefData(const int32_t& length = 0, 
            bool ok = false)
        : RefLength(length)
        , RefHasAlignments(ok)
    { }
};

typedef std::vector<RefData>      RefVector;
typedef std::vector<BamAlignment> BamAlignmentVector;

struct BamRegion {
  
    // data members
    int LeftRefID;
    int LeftPosition;
    int RightRefID;
    int RightPosition;
    
    // constructor
    BamRegion(const int& leftID   = -1, 
              const int& leftPos  = -1,
              const int& rightID  = -1,
              const int& rightPos = -1)
        : LeftRefID(leftID)
        , LeftPosition(leftPos)
        , RightRefID(rightID)
        , RightPosition(rightPos)
    { }
};

// ----------------------------------------------------------------
// Added: 3-35-2010 DWB
// Fixed: Routines to provide endian-correctness
// ----------------------------------------------------------------

// returns true if system is big endian
inline bool SystemIsBigEndian(void) {
   const uint16_t one = 0x0001;
   return ((*(char*) &one) == 0 );
}

// swaps endianness of 16-bit value 'in place'
inline void SwapEndian_16(int16_t& x) {
    x = ((x >> 8) | (x << 8));
}

inline void SwapEndian_16(uint16_t& x) {
    x = ((x >> 8) | (x << 8));
}

// swaps endianness of 32-bit value 'in-place'
inline void SwapEndian_32(int32_t& x) {
    x = ( (x >> 24) | 
         ((x << 8) & 0x00FF0000) | 
         ((x >> 8) & 0x0000FF00) | 
          (x << 24)
        );
}

inline void SwapEndian_32(uint32_t& x) {
    x = ( (x >> 24) | 
         ((x << 8) & 0x00FF0000) | 
         ((x >> 8) & 0x0000FF00) | 
          (x << 24)
        );
}

// swaps endianness of 64-bit value 'in-place'
inline void SwapEndian_64(int64_t& x) {
    x = ( (x >> 56) | 
         ((x << 40) & 0x00FF000000000000ll) |
         ((x << 24) & 0x0000FF0000000000ll) |
         ((x << 8)  & 0x000000FF00000000ll) |
         ((x >> 8)  & 0x00000000FF000000ll) |
         ((x >> 24) & 0x0000000000FF0000ll) |
         ((x >> 40) & 0x000000000000FF00ll) |
          (x << 56)
        );
}

inline void SwapEndian_64(uint64_t& x) {
    x = ( (x >> 56) | 
         ((x << 40) & 0x00FF000000000000ll) |
         ((x << 24) & 0x0000FF0000000000ll) |
         ((x << 8)  & 0x000000FF00000000ll) |
         ((x >> 8)  & 0x00000000FF000000ll) |
         ((x >> 24) & 0x0000000000FF0000ll) |
         ((x >> 40) & 0x000000000000FF00ll) |
          (x << 56)
        );
}

// swaps endianness of 'next 2 bytes' in a char buffer (in-place)
inline void SwapEndian_16p(char* data) {
    uint16_t& value = (uint16_t&)*data; 
    SwapEndian_16(value);
}

// swaps endianness of 'next 4 bytes' in a char buffer (in-place)
inline void SwapEndian_32p(char* data) {
    uint32_t& value = (uint32_t&)*data; 
    SwapEndian_32(value);
}

// swaps endianness of 'next 8 bytes' in a char buffer (in-place)
inline void SwapEndian_64p(char* data) {
    uint64_t& value = (uint64_t&)*data; 
    SwapEndian_64(value);
}

// ----------------------------------------------------------------
// BamAlignment member methods

// constructors & destructor
inline BamAlignment::BamAlignment(void) { }

inline BamAlignment::BamAlignment(const BamAlignment& other)
    : Name(other.Name)
    , Length(other.Length)
    , QueryBases(other.QueryBases)
    , AlignedBases(other.AlignedBases)
    , Qualities(other.Qualities)
    , TagData(other.TagData)
    , RefID(other.RefID)
    , Position(other.Position)
    , Bin(other.Bin)
    , MapQuality(other.MapQuality)
    , AlignmentFlag(other.AlignmentFlag)
    , CigarData(other.CigarData)
    , MateRefID(other.MateRefID)
    , MatePosition(other.MatePosition)
    , InsertSize(other.InsertSize)
    , SupportData(other.SupportData)
{ }

inline BamAlignment::~BamAlignment(void) { }

// Queries against alignment flags
inline bool BamAlignment::IsDuplicate(void) const         { return ( (AlignmentFlag & DUPLICATE)     != 0 ); }
inline bool BamAlignment::IsFailedQC(void) const          { return ( (AlignmentFlag & QC_FAILED)     != 0 ); }
inline bool BamAlignment::IsFirstMate(void) const         { return ( (AlignmentFlag & READ_1)        != 0 ); }
inline bool BamAlignment::IsMapped(void) const            { return ( (AlignmentFlag & UNMAPPED)      == 0 ); }
inline bool BamAlignment::IsMateMapped(void) const        { return ( (AlignmentFlag & MATE_UNMAPPED) == 0 ); }
inline bool BamAlignment::IsMateReverseStrand(void) const { return ( (AlignmentFlag & MATE_REVERSE)  != 0 ); }
inline bool BamAlignment::IsPaired(void) const            { return ( (AlignmentFlag & PAIRED)        != 0 ); }
inline bool BamAlignment::IsPrimaryAlignment(void) const  { return ( (AlignmentFlag & SECONDARY)     == 0 ); }
inline bool BamAlignment::IsProperPair(void) const        { return ( (AlignmentFlag & PROPER_PAIR)   != 0 ); }
inline bool BamAlignment::IsReverseStrand(void) const     { return ( (AlignmentFlag & REVERSE)       != 0 ); }
inline bool BamAlignment::IsSecondMate(void) const        { return ( (AlignmentFlag & READ_2)        != 0 ); }

// Manipulate alignment flags 
inline void BamAlignment::SetIsDuplicate(bool ok)          { if (ok) AlignmentFlag |= DUPLICATE;     else AlignmentFlag &= ~DUPLICATE; }
inline void BamAlignment::SetIsFailedQC(bool ok)           { if (ok) AlignmentFlag |= QC_FAILED;     else AlignmentFlag &= ~QC_FAILED; }
inline void BamAlignment::SetIsFirstMate(bool ok)          { if (ok) AlignmentFlag |= READ_1;        else AlignmentFlag &= ~READ_1; }
inline void BamAlignment::SetIsMateUnmapped(bool ok)       { if (ok) AlignmentFlag |= MATE_UNMAPPED; else AlignmentFlag &= ~MATE_UNMAPPED; }
inline void BamAlignment::SetIsMateReverseStrand(bool ok)  { if (ok) AlignmentFlag |= MATE_REVERSE;  else AlignmentFlag &= ~MATE_REVERSE; }
inline void BamAlignment::SetIsPaired(bool ok)             { if (ok) AlignmentFlag |= PAIRED;        else AlignmentFlag &= ~PAIRED; }
inline void BamAlignment::SetIsProperPair(bool ok)         { if (ok) AlignmentFlag |= PROPER_PAIR;   else AlignmentFlag &= ~PROPER_PAIR; }
inline void BamAlignment::SetIsReverseStrand(bool ok)      { if (ok) AlignmentFlag |= REVERSE;       else AlignmentFlag &= ~REVERSE; }
inline void BamAlignment::SetIsSecondaryAlignment(bool ok) { if (ok) AlignmentFlag |= SECONDARY;     else AlignmentFlag &= ~SECONDARY; }
inline void BamAlignment::SetIsSecondMate(bool ok)         { if (ok) AlignmentFlag |= READ_2;        else AlignmentFlag &= ~READ_2; }
inline void BamAlignment::SetIsUnmapped(bool ok)           { if (ok) AlignmentFlag |= UNMAPPED;      else AlignmentFlag &= ~UNMAPPED; }

// calculates alignment end position, based on starting position and CIGAR operations
inline 
int BamAlignment::GetEndPosition(bool usePadded) const {

    // initialize alignment end to starting position
    int alignEnd = Position;

    // iterate over cigar operations
    std::vector<CigarOp>::const_iterator cigarIter = CigarData.begin();
    std::vector<CigarOp>::const_iterator cigarEnd  = CigarData.end();
    for ( ; cigarIter != cigarEnd; ++cigarIter) {
	const char cigarType = (*cigarIter).Type;
	if ( cigarType == 'M' || cigarType == 'D' || cigarType == 'N' ) {
	    alignEnd += (*cigarIter).Length;
	} 
        else if ( usePadded && cigarType == 'I' ) {
            alignEnd += (*cigarIter).Length;
        }
    }
    return alignEnd;
}

inline
bool BamAlignment::AddTag(const std::string& tag, const std::string& type, const std::string& value) {
  
    if ( SupportData.HasCoreOnly ) return false;
    if ( tag.size() != 2 || type.size() != 1 ) return false;
    if ( type != "Z" && type != "H" ) return false;
  
    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) return false;
  
    // otherwise, copy tag data to temp buffer
    std::string newTag = tag + type + value;
    const int newTagDataLength = tagDataLength + newTag.size() + 1; // leave room for null-term
    char originalTagData[newTagDataLength];
    memcpy(originalTagData, TagData.c_str(), tagDataLength + 1);    // '+1' for TagData null-term
    
    // append newTag
    strcat(originalTagData + tagDataLength, newTag.data());  // removes original null-term, appends newTag + null-term
    
    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData;
    TagData.assign(newTagData, newTagDataLength);
    
    // return success
    return true;
}

inline
bool BamAlignment::AddTag(const std::string& tag, const std::string& type, const uint32_t& value) {
  
    if ( SupportData.HasCoreOnly ) return false;
    if ( tag.size() != 2 || type.size() != 1 ) return false;
    if ( type == "f" || type == "Z" || type == "H" ) return false;
  
    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) return false;
  
    // otherwise, convert value to string
    union { unsigned int value; char valueBuffer[sizeof(unsigned int)]; } un;
    un.value = value;

    // copy original tag data to temp buffer
    std::string newTag = tag + type;
    const int newTagDataLength = tagDataLength + newTag.size() + 4; // leave room for new integer
    char originalTagData[newTagDataLength];
    memcpy(originalTagData, TagData.c_str(), tagDataLength + 1);    // '+1' for TagData null-term
    
    // append newTag
    strcat(originalTagData + tagDataLength, newTag.data());
    memcpy(originalTagData + tagDataLength + newTag.size(), un.valueBuffer, sizeof(unsigned int));
    
    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData;
    TagData.assign(newTagData, newTagDataLength);
    
    // return success
    return true;
}

inline
bool BamAlignment::AddTag(const std::string& tag, const std::string& type, const int32_t& value) {
    return AddTag(tag, type, (const uint32_t&)value);
}

inline
bool BamAlignment::AddTag(const std::string& tag, const std::string& type, const float& value) {
  
    if ( SupportData.HasCoreOnly ) return false;
    if ( tag.size() != 2 || type.size() != 1 ) return false;
    if ( type == "Z" || type == "H" ) return false;
  
    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) return false;
  
    // otherwise, convert value to string
    union { float value; char valueBuffer[sizeof(float)]; } un;
    un.value = value;

    // copy original tag data to temp buffer
    std::string newTag = tag + type;
    const int newTagDataLength = tagDataLength + newTag.size() + 4; // leave room for new float
    char originalTagData[newTagDataLength];
    memcpy(originalTagData, TagData.c_str(), tagDataLength + 1);    // '+1' for TagData null-term
    
    // append newTag
    strcat(originalTagData + tagDataLength, newTag.data());
    memcpy(originalTagData + tagDataLength + newTag.size(), un.valueBuffer, sizeof(float));
    
    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData;
    TagData.assign(newTagData, newTagDataLength);
    
    // return success
    return true;
}

inline
bool BamAlignment::EditTag(const std::string& tag, const std::string& type, const std::string& value) {
  
    if ( SupportData.HasCoreOnly ) return false;
    if ( tag.size() != 2 || type.size() != 1 ) return false;
    if ( type != "Z" && type != "H" ) return false;
  
    // localize the tag data
    char* pOriginalTagData = (char*)TagData.data();
    char* pTagData = pOriginalTagData;
    const unsigned int originalTagDataLength = TagData.size();
    
    unsigned int newTagDataLength = 0;
    unsigned int numBytesParsed = 0;
    
    // if tag found, store data in readGroup, return success
    if ( FindTag(tag, pTagData, originalTagDataLength, numBytesParsed) ) {
        
        // make sure array is more than big enough
        char newTagData[originalTagDataLength + value.size()];  

        // copy original tag data up til desired tag
        const unsigned int beginningTagDataLength = numBytesParsed;
        newTagDataLength += beginningTagDataLength;
        memcpy(newTagData, pOriginalTagData, numBytesParsed);
      
        // copy new VALUE in place of current tag data
        const unsigned int dataLength = strlen(value.c_str());
        memcpy(newTagData + beginningTagDataLength, (char*)value.c_str(), dataLength+1 );
        
        // skip to next tag (if tag for removal is last, return true) 
        const char* pTagStorageType = pTagData - 1;
        if ( !SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) ) return true;
         
        // copy everything from current tag (the next one after tag for removal) to end
        const unsigned int skippedDataLength = (numBytesParsed - beginningTagDataLength);
        const unsigned int endTagOffset      = beginningTagDataLength + dataLength + 1;
        const unsigned int endTagDataLength  = originalTagDataLength - beginningTagDataLength - skippedDataLength;
        memcpy(newTagData + endTagOffset, pTagData, endTagDataLength);
        
        // ensure null-terminator
        newTagData[ endTagOffset + endTagDataLength + 1 ] = 0;
        
        // save new tag data
        TagData.assign(newTagData, endTagOffset + endTagDataLength);
        return true;
    }
    
    // tag not found, attempt AddTag
    else return AddTag(tag, type, value);
}

inline
bool BamAlignment::EditTag(const std::string& tag, const std::string& type, const uint32_t& value) {
  
    if ( SupportData.HasCoreOnly ) return false;
    if ( tag.size() != 2 || type.size() != 1 ) return false;
    if ( type == "f" || type == "Z" || type == "H" ) return false;
    
     // localize the tag data
    char* pOriginalTagData = (char*)TagData.data();
    char* pTagData = pOriginalTagData;
    const unsigned int originalTagDataLength = TagData.size();
    
    unsigned int newTagDataLength = 0;
    unsigned int numBytesParsed = 0;
    
    // if tag found, store data in readGroup, return success
    if ( FindTag(tag, pTagData, originalTagDataLength, numBytesParsed) ) {
        
        // make sure array is more than big enough
        char newTagData[originalTagDataLength + sizeof(value)];  

        // copy original tag data up til desired tag
        const unsigned int beginningTagDataLength = numBytesParsed;
        newTagDataLength += beginningTagDataLength;
        memcpy(newTagData, pOriginalTagData, numBytesParsed);
      
        // copy new VALUE in place of current tag data
        union { unsigned int value; char valueBuffer[sizeof(unsigned int)]; } un;
        un.value = value;
        memcpy(newTagData + beginningTagDataLength, un.valueBuffer, sizeof(unsigned int));
        
        // skip to next tag (if tag for removal is last, return true) 
        const char* pTagStorageType = pTagData - 1;
        if ( !SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) ) return true;
         
        // copy everything from current tag (the next one after tag for removal) to end
        const unsigned int skippedDataLength = (numBytesParsed - beginningTagDataLength);
        const unsigned int endTagOffset      = beginningTagDataLength + sizeof(unsigned int);
        const unsigned int endTagDataLength  = originalTagDataLength - beginningTagDataLength - skippedDataLength;
        memcpy(newTagData + endTagOffset, pTagData, endTagDataLength);
        
        // ensure null-terminator
        newTagData[ endTagOffset + endTagDataLength + 1 ] = 0;
        
        // save new tag data
        TagData.assign(newTagData, endTagOffset + endTagDataLength);
        return true;
    }
    
    // tag not found, attempt AddTag
    else return AddTag(tag, type, value);
}

inline
bool BamAlignment::EditTag(const std::string& tag, const std::string& type, const int32_t& value) {
    return EditTag(tag, type, (const uint32_t&)value);
}

inline
bool BamAlignment::EditTag(const std::string& tag, const std::string& type, const float& value) {
  
    if ( SupportData.HasCoreOnly ) return false;
    if ( tag.size() != 2 || type.size() != 1 ) return false;
    if ( type == "Z" || type == "H" ) return false;
    
     // localize the tag data
    char* pOriginalTagData = (char*)TagData.data();
    char* pTagData = pOriginalTagData;
    const unsigned int originalTagDataLength = TagData.size();
    
    unsigned int newTagDataLength = 0;
    unsigned int numBytesParsed = 0;
    
    // if tag found, store data in readGroup, return success
    if ( FindTag(tag, pTagData, originalTagDataLength, numBytesParsed) ) {
        
        // make sure array is more than big enough
        char newTagData[originalTagDataLength + sizeof(value)];  

        // copy original tag data up til desired tag
        const unsigned int beginningTagDataLength = numBytesParsed;
        newTagDataLength += beginningTagDataLength;
        memcpy(newTagData, pOriginalTagData, numBytesParsed);
      
        // copy new VALUE in place of current tag data
        union { float value; char valueBuffer[sizeof(float)]; } un;
        un.value = value;
        memcpy(newTagData + beginningTagDataLength, un.valueBuffer, sizeof(float));
        
        // skip to next tag (if tag for removal is last, return true) 
        const char* pTagStorageType = pTagData - 1;
        if ( !SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) ) return true;
         
        // copy everything from current tag (the next one after tag for removal) to end
        const unsigned int skippedDataLength = (numBytesParsed - beginningTagDataLength);
        const unsigned int endTagOffset      = beginningTagDataLength + sizeof(float);
        const unsigned int endTagDataLength  = originalTagDataLength - beginningTagDataLength - skippedDataLength;
        memcpy(newTagData + endTagOffset, pTagData, endTagDataLength);
        
        // ensure null-terminator
        newTagData[ endTagOffset + endTagDataLength + 1 ] = 0;
        
        // save new tag data
        TagData.assign(newTagData, endTagOffset + endTagDataLength);
        return true;
    }
    
    // tag not found, attempt AddTag
    else return AddTag(tag, type, value);
}

// get "NM" tag data - originally contributed by Aaron Quinlan
// stores data in 'editDistance', returns success/fail
inline 
bool BamAlignment::GetEditDistance(uint32_t& editDistance) const { 
    return GetTag("NM", (uint32_t&)editDistance);
}

// get "RG" tag data
// stores data in 'readGroup', returns success/fail
inline 
bool BamAlignment::GetReadGroup(std::string& readGroup) const {
    return GetTag("RG", readGroup);
}

inline
bool BamAlignment::GetTag(const std::string& tag, std::string& destination) const {

    // make sure tag data exists
    if ( SupportData.HasCoreOnly || TagData.empty() ) 
        return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag found, store data in readGroup, return success
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        const unsigned int dataLength = strlen(pTagData);
        destination.clear();
        destination.resize(dataLength);
        memcpy( (char*)destination.data(), pTagData, dataLength );
        return true;
    }
    
    // tag not found, return failure
    return false;
}

inline
bool BamAlignment::GetTag(const std::string& tag, uint32_t& destination) const {
  
    // make sure tag data exists
    if ( SupportData.HasCoreOnly || TagData.empty() ) 
        return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag found, determine data byte-length, store data in readGroup, return success
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        
        // determine data byte-length
        const char type = *(pTagData - 1);
        int destinationLength = 0;
        switch (type) {
            // 1 byte data
            case 'A':
            case 'c':
            case 'C':
                destinationLength = 1;
                break;

            // 2 byte data
            case 's':
            case 'S':
                destinationLength = 2;
                break;

            // 4 byte data
            case 'i':
            case 'I':
                destinationLength = 4;
                break;

            // unsupported type for integer destination (float or var-length strings)
            case 'f':
            case 'Z':
            case 'H':
                printf("ERROR: Cannot store tag of type %c in integer destination\n", type);
                return false;

            // unknown tag type
            default:
                printf("ERROR: Unknown tag storage class encountered: [%c]\n", type);
                return false;
        }
          
        // store in destination
        destination = 0;
        memcpy(&destination, pTagData, destinationLength);
        return true;
    }
    
    // tag not found, return failure
    return false;
}

inline
bool BamAlignment::GetTag(const std::string& tag, int32_t& destination) const {
    return GetTag(tag, (uint32_t&)destination);
}

inline
bool BamAlignment::GetTag(const std::string& tag, float& destination) const {
  
    // make sure tag data exists
    if ( SupportData.HasCoreOnly || TagData.empty() ) 
        return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag found, determine data byte-length, store data in readGroup, return success
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        //pTagData += numBytesParsed;
        
        // determine data byte-length
        const char type = *(pTagData - 1);
        int destinationLength = 0;
        switch(type) {

            // 1 byte data
            case 'A':
            case 'c':
            case 'C':
                destinationLength = 1;
                break;

            // 2 byte data
            case 's':
            case 'S':
                destinationLength = 2;
                break;

            // 4 byte data
            case 'f':
            case 'i':
            case 'I':
                destinationLength = 4;
                break;
            
            // unsupported type (var-length strings)
            case 'Z':
            case 'H':
                printf("ERROR: Cannot store tag of type %c in integer destination\n", type);
                return false;

            // unknown tag type
            default:
                printf("ERROR: Unknown tag storage class encountered: [%c]\n", type);
                return false;
        }
          
        // store in destination
        destination = 0.0;
        memcpy(&destination, pTagData, destinationLength);
        return true;
    }
    
    // tag not found, return failure
    return false;
}

inline
bool BamAlignment::RemoveTag(const std::string& tag) {
  
    // BamAlignments fetched using BamReader::GetNextAlignmentCore() are not allowed
    // also, return false if no data present to remove
    if ( SupportData.HasCoreOnly || TagData.empty() ) return false;
  
    // localize the tag data
    char* pOriginalTagData = (char*)TagData.data();
    char* pTagData = pOriginalTagData;
    const unsigned int originalTagDataLength = TagData.size();
    unsigned int newTagDataLength = 0;
    unsigned int numBytesParsed = 0;
    
    // if tag found, store data in readGroup, return success
    if ( FindTag(tag, pTagData, originalTagDataLength, numBytesParsed) ) {
        
        char newTagData[originalTagDataLength];

        // copy original tag data up til desired tag
        pTagData -= 3;
        numBytesParsed -= 3;
        const unsigned int beginningTagDataLength = numBytesParsed;
        newTagDataLength += beginningTagDataLength;
        memcpy(newTagData, pOriginalTagData, numBytesParsed);
        
        // skip to next tag (if tag for removal is last, return true) 
        const char* pTagStorageType = pTagData + 2;
        pTagData       += 3;
        numBytesParsed += 3;
        if ( !SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) ) return true;
         
        // copy everything from current tag (the next one after tag for removal) to end
        const unsigned int skippedDataLength = (numBytesParsed - beginningTagDataLength);
        const unsigned int endTagDataLength = originalTagDataLength - beginningTagDataLength - skippedDataLength;
        memcpy(newTagData + beginningTagDataLength, pTagData, endTagDataLength );
        
        // save new tag data
        TagData.assign(newTagData, beginningTagDataLength + endTagDataLength);
        return true;
    }
    
    // tag not found, no removal - return failure
    return false;
}

inline
bool BamAlignment::FindTag(const std::string& tag, char* &pTagData, const unsigned int& tagDataLength, unsigned int& numBytesParsed) {

    while ( numBytesParsed < tagDataLength ) {

        const char* pTagType        = pTagData;
        const char* pTagStorageType = pTagData + 2;
        pTagData       += 3;
        numBytesParsed += 3;

        // check the current tag, return true on match
        if ( std::strncmp(pTagType, tag.c_str(), 2) == 0 ) 
            return true;

        // get the storage class and find the next tag
        if ( *pTagStorageType == '\0' ) return false; 
        if ( !SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) ) return false;
        if ( *pTagData == '\0' ) return false;
    }
  
    // checked all tags, none match
    return false;
}

inline
bool BamAlignment::SkipToNextTag(const char storageType, char* &pTagData, unsigned int& numBytesParsed) {
    
    switch(storageType) {

        case 'A':
        case 'c':
        case 'C':
            ++numBytesParsed;
            ++pTagData;
            break;

        case 's':
        case 'S':
            numBytesParsed += 2;
            pTagData       += 2;
            break;

        case 'f':
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
            // increment for null-terminator
            ++numBytesParsed;
            ++pTagData;
            break;

        default: 
            // error case
            printf("ERROR: Unknown tag storage class encountered: [%c]\n", storageType);
            return false;
    }
    
    // return success
    return true;
}

} // namespace BamTools

#endif // BAMAUX_H
