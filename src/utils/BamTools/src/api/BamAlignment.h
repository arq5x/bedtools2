// ***************************************************************************
// BamAlignment.h (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 16 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides the BamAlignment data structure
// ***************************************************************************

#ifndef BAMALIGNMENT_H
#define BAMALIGNMENT_H

#include "api/api_global.h"
#include "api/BamAux.h"
#include "api/BamConstants.h"
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace BamTools {

//! \cond
// forward declaration of BamAlignment's "friends"
namespace Internal {
    class BamReaderPrivate;
    class BamWriterPrivate;
} // namespace Internal
//! \endcond

// BamAlignment data structure
struct API_EXPORT BamAlignment {

    // constructors & destructor
    public:
        BamAlignment(void);
        BamAlignment(const BamAlignment& other);
        ~BamAlignment(void);

    // queries against alignment flags
    public:        
        bool IsDuplicate(void) const;         // returns true if this read is a PCR duplicate
        bool IsFailedQC(void) const;          // returns true if this read failed quality control
        bool IsFirstMate(void) const;         // returns true if alignment is first mate on read
        bool IsMapped(void) const;            // returns true if alignment is mapped
        bool IsMateMapped(void) const;        // returns true if alignment's mate is mapped
        bool IsMateReverseStrand(void) const; // returns true if alignment's mate mapped to reverse strand
        bool IsPaired(void) const;            // returns true if alignment part of paired-end read
        bool IsPrimaryAlignment(void) const;  // returns true if reported position is primary alignment
        bool IsProperPair(void) const;        // returns true if alignment is part of read that satisfied paired-end resolution
        bool IsReverseStrand(void) const;     // returns true if alignment mapped to reverse strand
        bool IsSecondMate(void) const;        // returns true if alignment is second mate on read

    // manipulate alignment flags
    public:        
        void SetIsDuplicate(bool ok);         // sets value of "PCR duplicate" flag
        void SetIsFailedQC(bool ok);          // sets value of "failed quality control" flag
        void SetIsFirstMate(bool ok);         // sets value of "alignment is first mate" flag
        void SetIsMapped(bool ok);            // sets value of "alignment is mapped" flag
        void SetIsMateMapped(bool ok);        // sets value of "alignment's mate is mapped" flag
        void SetIsMateReverseStrand(bool ok); // sets value of "alignment's mate mapped to reverse strand" flag
        void SetIsPaired(bool ok);            // sets value of "alignment part of paired-end read" flag
        void SetIsPrimaryAlignment(bool ok);  // sets value of "position is primary alignment" flag
        void SetIsProperPair(bool ok);        // sets value of "alignment is part of read that satisfied paired-end resolution" flag
        void SetIsReverseStrand(bool ok);     // sets value of "alignment mapped to reverse strand" flag
        void SetIsSecondMate(bool ok);        // sets value of "alignment is second mate on read" flag

    // tag data access methods
    public:

        // add a new tag
        template<typename T> bool AddTag(const std::string& tag, const std::string& type, const T& value);
        template<typename T> bool AddTag(const std::string& tag, const std::vector<T>& values);

        // edit (or append) tag
        template<typename T> bool EditTag(const std::string& tag, const std::string& type, const T& value);
        template<typename T> bool EditTag(const std::string& tag, const std::vector<T>& values);

        // retrieves tag data
        template<typename T> bool GetTag(const std::string& tag, T& destination) const;
        template<typename T> bool GetTag(const std::string& tag, std::vector<T>& destination) const;

        // retrieves the SAM/BAM type-code for requested tag name
        bool GetTagType(const std::string& tag, char& type) const;

        // returns true if alignment has a record for this tag name
        bool HasTag(const std::string& tag) const;

        // removes a tag
        void RemoveTag(const std::string& tag);

    // additional methods
    public:
        // populates alignment string fields
        bool BuildCharData(void);

        // calculates alignment end position
        int GetEndPosition(bool usePadded = false, bool closedInterval = false) const;

        // returns a description of the last error that occurred
        std::string GetErrorString(void) const;

        // retrieves the size, read locations and reference locations of soft-clip operations
        bool GetSoftClips(std::vector<int>& clipSizes,
                          std::vector<int>& readPositions,
                          std::vector<int>& genomePositions,
                          bool usePadded = false) const;

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

    //! \internal
    // internal utility methods
    private:
        bool FindTag(const std::string& tag,
                     char*& pTagData,
                     const unsigned int& tagDataLength,
                     unsigned int& numBytesParsed) const;
        bool IsValidSize(const std::string& tag, const std::string& type) const;
        void SetErrorString(const std::string& where, const std::string& what) const;
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

        mutable std::string ErrorString; // mutable to allow updates even in logically const methods
    //! \endinternal
};

// ---------------------------------------------------------
// BamAlignment tag access methods

/*! \fn bool AddTag(const std::string& tag, const std::string& type, const T& value)
    \brief Adds a field to the BAM tags.

    Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.

    \param[in] tag   2-character tag name
    \param[in] type  1-character tag type
    \param[in] value data to store
    \return \c true if the \b new tag was added successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
template<typename T>
inline bool BamAlignment::AddTag(const std::string& tag, const std::string& type, const T& value) {

    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();

    // check tag/type size
    if ( !IsValidSize(tag, type) ) {
        // TODO: set error string?
        return false;
    }

    // check that storage type code is OK for T
    if ( !TagTypeHelper<T>::CanConvertTo(type.at(0)) ) {
        // TODO: set error string?
        return false;
    }

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        // TODO: set error string?
        return false;
    }

    // otherwise, convert value to string
    union { T value; char valueBuffer[sizeof(T)]; } un;
    un.value = value;

    // copy original tag data to temp buffer
    const std::string newTag = tag + type;
    const size_t newTagDataLength = tagDataLength + newTag.size() + sizeof(T); // leave room for new T
    RaiiBuffer originalTagData(newTagDataLength);
    memcpy(originalTagData.Buffer, TagData.c_str(), tagDataLength + 1);    // '+1' for TagData null-term

    // append newTag
    strcat(originalTagData.Buffer + tagDataLength, newTag.data());
    memcpy(originalTagData.Buffer + tagDataLength + newTag.size(), un.valueBuffer, sizeof(T));

    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData.Buffer;
    TagData.assign(newTagData, newTagDataLength);
    return true;
}

template<>
inline bool BamAlignment::AddTag<std::string>(const std::string& tag,
                                              const std::string& type,
                                              const std::string& value)
{
    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();

    // check tag/type size
    if ( !IsValidSize(tag, type) ) {
        // TODO: set error string?
        return false;
    }

    // check that storage type code is OK for string
    if ( !TagTypeHelper<std::string>::CanConvertTo(type.at(0)) ) {
        // TODO: set error string?
        return false;
    }

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        // TODO: set error string?
        return false;
    }

    // otherwise, copy tag data to temp buffer
    const std::string newTag = tag + type + value;
    const size_t newTagDataLength = tagDataLength + newTag.size() + 1; // leave room for null-term
    RaiiBuffer originalTagData(newTagDataLength);
    memcpy(originalTagData.Buffer, TagData.c_str(), tagDataLength + 1);    // '+1' for TagData null-term

    // append newTag (removes original null-term, then appends newTag + null-term)
    strcat(originalTagData.Buffer + tagDataLength, newTag.data());

    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData.Buffer;
    TagData.assign(newTagData, newTagDataLength);
    return true;
}

/*! \fn template<typename T> bool AddTag(const std::string& tag, const std::vector<T>& values)
    \brief Adds a numeric array field to the BAM tags.

    Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.

    \param[in] tag    2-character tag name
    \param[in] values vector of data values to store
    \return \c true if the \b new tag was added successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
template<typename T>
inline bool BamAlignment::AddTag(const std::string& tag, const std::vector<T>& values) {

    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();

    // check for valid tag name length
    if ( tag.size() != Constants::BAM_TAG_TAGSIZE )
        return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        // TODO: set error string?
        return false;
    }

    // build new tag's base information
    char newTagBase[Constants::BAM_TAG_ARRAYBASE_SIZE];
    memcpy( newTagBase, tag.c_str(), Constants::BAM_TAG_TAGSIZE );
    newTagBase[2] = Constants::BAM_TAG_TYPE_ARRAY;
    newTagBase[3] = TagTypeHelper<T>::TypeCode();

    // add number of array elements to newTagBase
    const int32_t numElements  = values.size();
    memcpy(newTagBase + 4, &numElements, sizeof(int32_t));

    // copy current TagData string to temp buffer, leaving room for new tag's contents
    const size_t newTagDataLength = tagDataLength +
                                    Constants::BAM_TAG_ARRAYBASE_SIZE +
                                    numElements*sizeof(T);
    RaiiBuffer originalTagData(newTagDataLength);
    memcpy(originalTagData.Buffer, TagData.c_str(), tagDataLength+1); // '+1' for TagData's null-term

    // write newTagBase (removes old null term)
    strcat(originalTagData.Buffer + tagDataLength, (const char*)newTagBase);

    // add vector elements to tag
    int elementsBeginOffset = tagDataLength + Constants::BAM_TAG_ARRAYBASE_SIZE;
    for ( int i = 0 ; i < numElements; ++i ) {
        const T& value = values.at(i);
        memcpy(originalTagData.Buffer + elementsBeginOffset + i*sizeof(T), &value, sizeof(T));
    }

    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData.Buffer;
    TagData.assign(newTagData, newTagDataLength);
    return true;
}

/*! \fn template<typename T> bool EditTag(const std::string& tag, const std::string& type, const T& value)
    \brief Edits a BAM tag field.

    If \a tag does not exist, a new entry is created.

    \param tag[in]   2-character tag name
    \param type[in]  1-character tag type (must be "Z" or "H")
    \param value[in] new data value

    \return \c true if the tag was modified/created successfully

    \sa BamAlignment::RemoveTag()
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
template<typename T>
inline bool BamAlignment::EditTag(const std::string& tag, const std::string& type, const T& value) {

    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();

    // remove existing tag if present, then append tag with new value
    if ( HasTag(tag) )
        RemoveTag(tag);
    return AddTag(tag, type, value);
}

/*! \fn template<typename T> bool EditTag(const std::string& tag, const std::vector<T>& values)
    \brief Edits a BAM tag field containing a numeric array.

    If \a tag does not exist, a new entry is created.

    \param tag[in]   2-character tag name
    \param value[in] vector of data values

    \return \c true if the tag was modified/created successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
template<typename T>
inline bool BamAlignment::EditTag(const std::string& tag, const std::vector<T>& values) {

    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();

    // remove existing tag if present, then append tag with new values
    if ( HasTag(tag) )
        RemoveTag(tag);
    return AddTag(tag, values);
}


/*! \fn template<typename T> bool GetTag(const std::string& tag, T& destination) const
    \brief Retrieves the value associated with a BAM tag.

    \param tag[in]          2-character tag name
    \param destination[out] retrieved value
    \return \c true if found
*/
template<typename T>
inline bool BamAlignment::GetTag(const std::string& tag, T& destination) const {

    // skip if alignment is core-only
    if ( SupportData.HasCoreOnly ) {
        // TODO: set error string?
        return false;
    }

    // skip if no tags present
    if ( TagData.empty() ) {
        // TODO: set error string?
        return false;
    }

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // return failure if tag not found
    if ( !FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        // TODO: set error string?
        return false;
    }

    // fetch data type
    const char type = *(pTagData - 1);
    if ( !TagTypeHelper<T>::CanConvertFrom(type) ) {
        // TODO: set error string ?
        return false;
    }

    // determine data length
    int destinationLength = 0;
    switch ( type ) {

        // 1 byte data
        case (Constants::BAM_TAG_TYPE_ASCII) :
        case (Constants::BAM_TAG_TYPE_INT8)  :
        case (Constants::BAM_TAG_TYPE_UINT8) :
            destinationLength = 1;
            break;

        // 2 byte data
        case (Constants::BAM_TAG_TYPE_INT16)  :
        case (Constants::BAM_TAG_TYPE_UINT16) :
            destinationLength = 2;
            break;

        // 4 byte data
        case (Constants::BAM_TAG_TYPE_INT32)  :
        case (Constants::BAM_TAG_TYPE_UINT32) :
        case (Constants::BAM_TAG_TYPE_FLOAT)  :
            destinationLength = 4;
            break;

        // var-length types not supported for numeric destination
        case (Constants::BAM_TAG_TYPE_STRING) :
        case (Constants::BAM_TAG_TYPE_HEX)    :
        case (Constants::BAM_TAG_TYPE_ARRAY)  :
            SetErrorString("BamAlignment::GetTag",
                           "cannot store variable length tag data into a numeric destination");
            return false;

        // unrecognized tag type
        default:
            const std::string message = std::string("invalid tag type: ") + type;
            SetErrorString("BamAlignment::GetTag", message);
            return false;
    }

    // store data in destination
    destination = 0;
    memcpy(&destination, pTagData, destinationLength);

    // return success
    return true;
}

template<>
inline bool BamAlignment::GetTag<std::string>(const std::string& tag,
                                              std::string& destination) const
{
    // skip if alignment is core-only
    if ( SupportData.HasCoreOnly ) {
        // TODO: set error string?
        return false;
    }

    // skip if no tags present
    if ( TagData.empty() ) {
        // TODO: set error string?
        return false;
    }

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // return failure if tag not found
    if ( !FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        // TODO: set error string?
        return false;
    }

    // otherwise copy data into destination
    const unsigned int dataLength = strlen(pTagData);
    destination.clear();
    destination.resize(dataLength);
    memcpy( (char*)destination.data(), pTagData, dataLength );

    // return success
    return true;
}

/*! \fn template<typename T> bool GetTag(const std::string& tag, std::vector<T>& destination) const
    \brief Retrieves the numeric array associated with a BAM tag.

    \param tag[in]          2-character tag name
    \param destination[out] retrieved values
    \return \c true if found
*/
template<typename T>
inline bool BamAlignment::GetTag(const std::string& tag, std::vector<T>& destination) const {

    // skip if alignment is core-only
    if ( SupportData.HasCoreOnly ) {
        // TODO: set error string?
        return false;
    }

    // skip if no tags present
    if ( TagData.empty() ) {
        // TODO: set error string?
        return false;
    }

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // return false if tag not found
    if ( !FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        // TODO: set error string?
        return false;
    }

    // check that tag is array type
    const char tagType = *(pTagData - 1);
    if ( tagType != Constants::BAM_TAG_TYPE_ARRAY ) {
        SetErrorString("BamAlignment::GetTag", "cannot store a non-array tag in array destination");
        return false;
    }

    // fetch element type
    const char elementType = *pTagData;
    if ( !TagTypeHelper<T>::CanConvertFrom(elementType) ) {
        // TODO: set error string ?
        return false;
    }
    ++pTagData;

    // calculate length of each element in tag's array
    int elementLength = 0;
    switch ( elementType ) {
        case (Constants::BAM_TAG_TYPE_ASCII) :
        case (Constants::BAM_TAG_TYPE_INT8)  :
        case (Constants::BAM_TAG_TYPE_UINT8) :
            elementLength = sizeof(uint8_t);
            break;

        case (Constants::BAM_TAG_TYPE_INT16)  :
        case (Constants::BAM_TAG_TYPE_UINT16) :
            elementLength = sizeof(uint16_t);
            break;

        case (Constants::BAM_TAG_TYPE_INT32)  :
        case (Constants::BAM_TAG_TYPE_UINT32) :
        case (Constants::BAM_TAG_TYPE_FLOAT)  :
            elementLength = sizeof(uint32_t);
            break;

        // var-length types not supported for numeric destination
        case (Constants::BAM_TAG_TYPE_STRING) :
        case (Constants::BAM_TAG_TYPE_HEX)    :
        case (Constants::BAM_TAG_TYPE_ARRAY)  :
            SetErrorString("BamAlignment::GetTag",
                           "invalid array data, variable-length elements are not allowed");
            return false;

        // unknown tag type
        default:
            const std::string message = std::string("invalid array element type: ") + elementType;
            SetErrorString("BamAlignment::GetTag", message);
            return false;
    }

    // get number of elements
    int32_t numElements;
    memcpy(&numElements, pTagData, sizeof(int32_t));
    pTagData += 4;
    destination.clear();
    destination.reserve(numElements);

    // read in elements
    T value;
    for ( int i = 0 ; i < numElements; ++i ) {
        memcpy(&value, pTagData, sizeof(T));
        pTagData += sizeof(T);
        destination.push_back(value);
    }

    // return success
    return true;
}

typedef std::vector<BamAlignment> BamAlignmentVector;

} // namespace BamTools

#endif // BAMALIGNMENT_H
