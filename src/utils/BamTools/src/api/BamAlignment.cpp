// ***************************************************************************
// BamAlignment.cpp (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 22 April 2011 (DB)
// ---------------------------------------------------------------------------
// Provides the BamAlignment data structure
// ***************************************************************************

#include <api/BamAlignment.h>
#include <api/BamConstants.h>
using namespace BamTools;

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iostream>
#include <map>
#include <utility>
using namespace std;

/*! \class BamTools::BamAlignment
    \brief The main BAM alignment data structure.

    Provides methods to query/modify BAM alignment data fields.
*/
/*! \var BamAlignment::Name
    \brief read name
*/
/*! \var BamAlignment::Length
    \brief length of query sequence
*/
/*! \var BamAlignment::QueryBases
    \brief 'original' sequence (as reported from sequencing machine)
*/
/*! \var BamAlignment::AlignedBases
    \brief 'aligned' sequence (includes any indels, padding, clipping)
*/
/*! \var BamAlignment::Qualities
    \brief FASTQ qualities (ASCII characters, not numeric values)
*/
/*! \var BamAlignment::TagData
    \brief tag data (use the provided methods to query/modify)
*/
/*! \var BamAlignment::RefID
    \brief ID number for reference sequence
*/
/*! \var BamAlignment::Position
    \brief position (0-based) where alignment starts
*/
/*! \var BamAlignment::Bin
    \brief BAM (standard) index bin number for this alignment
*/
/*! \var BamAlignment::MapQuality
    \brief mapping quality score
*/
/*! \var BamAlignment::AlignmentFlag
    \brief alignment bit-flag (use the provided methods to query/modify)
*/
/*! \var BamAlignment::CigarData
    \brief CIGAR operations for this alignment
*/
/*! \var BamAlignment::MateRefID
    \brief ID number for reference sequence where alignment's mate was aligned
*/
/*! \var BamAlignment::MatePosition
    \brief position (0-based) where alignment's mate starts
*/
/*! \var BamAlignment::InsertSize
    \brief mate-pair insert size
*/
/*! \var BamAlignment::Filename
    \brief name of BAM file which this alignment comes from
*/

/*! \fn BamAlignment::BamAlignment(void)
    \brief constructor
*/
BamAlignment::BamAlignment(void)
    : RefID(-1)
    , Position(-1)
    , MateRefID(-1)
    , MatePosition(-1)
    , InsertSize(0)
{ }

/*! \fn BamAlignment::BamAlignment(const BamAlignment& other)
    \brief copy constructor
*/
BamAlignment::BamAlignment(const BamAlignment& other)
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
    , Filename(other.Filename)
    , SupportData(other.SupportData)
{ }

/*! \fn BamAlignment::~BamAlignment(void)
    \brief destructor
*/
BamAlignment::~BamAlignment(void) { }

/*! \fn bool BamAlignment::AddTag(const std::string& tag, const std::string& type, const std::string& value)
    \brief Adds a field with string data to the BAM tags.

    Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.

    \param tag   2-character tag name
    \param type  1-character tag type (must be "Z" or "H")
    \param value string data to store

    \return \c true if the \b new tag was added successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::AddTag(const std::string& tag, const std::string& type, const std::string& value) {
  
    // skip if core data not parsed
    if ( SupportData.HasCoreOnly ) return false;

    // validate tag/type size & that type is OK for string value
    if ( !IsValidSize(tag, type) ) return false;
    if ( type.at(0) != Constants::BAM_TAG_TYPE_STRING &&
         type.at(0) != Constants::BAM_TAG_TYPE_HEX
       )
    {
        return false;
    }
  
    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) )
        return false;
  
    // otherwise, copy tag data to temp buffer
    string newTag = tag + type + value;
    const int newTagDataLength = tagDataLength + newTag.size() + 1; // leave room for null-term
    char* originalTagData = new char[newTagDataLength];
    memcpy(originalTagData, TagData.c_str(), tagDataLength + 1);    // '+1' for TagData null-term

    // append newTag
    strcat(originalTagData + tagDataLength, newTag.data());  // removes original null-term, appends newTag + null-term

    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData;
    TagData.assign(newTagData, newTagDataLength);
    
    delete[] originalTagData;

    // return success
    return true;
}

/*! \fn bool BamAlignment::AddTag(const std::string& tag, const std::string& type, const uint32_t& value)
    \brief Adds a field with unsigned integer data to the BAM tags.

    Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.

    \param tag   2-character tag name
    \param type  1-character tag type (must NOT be "f", "Z", "H", or "B")
    \param value unsigned int data to store

    \return \c true if the \b new tag was added successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::AddTag(const std::string& tag, const std::string& type, const uint32_t& value) {
  
    // skip if core data not parsed
    if ( SupportData.HasCoreOnly ) return false;

    // validate tag/type size & that type is OK for uint32_t value
    if ( !IsValidSize(tag, type) ) return false;
    if ( type.at(0) == Constants::BAM_TAG_TYPE_FLOAT  ||
         type.at(0) == Constants::BAM_TAG_TYPE_STRING ||
         type.at(0) == Constants::BAM_TAG_TYPE_HEX    ||
         type.at(0) == Constants::BAM_TAG_TYPE_ARRAY
       )
    {
        return false;
    }
  
    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) )
        return false;
  
    // otherwise, convert value to string
    union { uint32_t value; char valueBuffer[sizeof(uint32_t)]; } un;
    un.value = value;

    // copy original tag data to temp buffer
    string newTag = tag + type;
    const int newTagDataLength = tagDataLength + newTag.size() + 4; // leave room for new integer
    char* originalTagData = new char[newTagDataLength];
    memcpy(originalTagData, TagData.c_str(), tagDataLength + 1);    // '+1' for TagData null-term
    
    // append newTag
    strcat(originalTagData + tagDataLength, newTag.data());
    memcpy(originalTagData + tagDataLength + newTag.size(), un.valueBuffer, sizeof(uint32_t));
    
    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData;
    TagData.assign(newTagData, newTagDataLength);
    delete[] originalTagData;
    
    // return success
    return true;
}

/*! \fn bool BamAlignment::AddTag(const std::string& tag, const std::string& type, const int32_t& value)
    \brief Adds a field with signed integer data to the BAM tags.

    Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.

    \param tag   2-character tag name
    \param type  1-character tag type (must NOT be "f", "Z", "H", or "B")
    \param value signed int data to store

    \return \c true if the \b new tag was added successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::AddTag(const std::string& tag, const std::string& type, const int32_t& value) {
    return AddTag(tag, type, (const uint32_t&)value);
}

/*! \fn bool BamAlignment::AddTag(const std::string& tag, const std::string& type, const float& value)
    \brief Adds a field with floating-point data to the BAM tags.

    Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.

    \param tag  2-character tag name
    \param type  1-character tag type (must NOT be "Z", "H", or "B")
    \param value float data to store

    \return \c true if the \b new tag was added successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::AddTag(const std::string& tag, const std::string& type, const float& value) {
  
    // skip if core data not parsed
    if ( SupportData.HasCoreOnly ) return false;

    // validate tag/type size & that type is OK for float value
    if ( !IsValidSize(tag, type) ) return false;
    if ( type.at(0) == Constants::BAM_TAG_TYPE_STRING ||
         type.at(0) == Constants::BAM_TAG_TYPE_HEX    ||
         type.at(0) == Constants::BAM_TAG_TYPE_ARRAY
       )
    {
        return false;
    }

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) )
        return false;
  
    // otherwise, convert value to string
    union { float value; char valueBuffer[sizeof(float)]; } un;
    un.value = value;

    // copy original tag data to temp buffer
    string newTag = tag + type;
    const int newTagDataLength = tagDataLength + newTag.size() + 4; // leave room for new float
    char* originalTagData = new char[newTagDataLength];
    memcpy(originalTagData, TagData.c_str(), tagDataLength + 1);    // '+1' for TagData null-term
    
    // append newTag
    strcat(originalTagData + tagDataLength, newTag.data());
    memcpy(originalTagData + tagDataLength + newTag.size(), un.valueBuffer, sizeof(float));
    
    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData;
    TagData.assign(newTagData, newTagDataLength);
    
    delete[] originalTagData;

    // return success
    return true;
}

/*! \fn bool AddTag(const std::string& tag, const std::vector<uint8_t>& values);
    \brief Adds a numeric array field to the BAM tags.

    Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.

    \param tag    2-character tag name
    \param values vector of uint8_t values to store

    \return \c true if the \b new tag was added successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::AddTag(const std::string& tag, const std::vector<uint8_t>& values) {

    // skip if core data not parsed
    if ( SupportData.HasCoreOnly ) return false;

    // check for valid tag length
    if ( tag.size() != Constants::BAM_TAG_TAGSIZE ) return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) )
        return false;

    // build new tag's base information
    char newTagBase[Constants::BAM_TAG_ARRAYBASE_SIZE];
    memcpy( newTagBase, tag.c_str(), Constants::BAM_TAG_TAGSIZE );
    newTagBase[2] = Constants::BAM_TAG_TYPE_ARRAY;
    newTagBase[3] = Constants::BAM_TAG_TYPE_UINT8;

    // add number of array elements to newTagBase
    const int32_t numElements  = values.size();
    memcpy(newTagBase + 4, &numElements, sizeof(int32_t));

    // copy current TagData string to temp buffer, leaving room for new tag's contents
    const int newTagDataLength = tagDataLength +
                                 Constants::BAM_TAG_ARRAYBASE_SIZE +
                                 numElements*sizeof(uint8_t);
    char* originalTagData = new char[newTagDataLength];
    memcpy(originalTagData, TagData.c_str(), tagDataLength+1); // '+1' for TagData's null-term

    // write newTagBase (removes old null term)
    strcat(originalTagData + tagDataLength, (const char*)newTagBase);

    // add vector elements to tag
    int elementsBeginOffset = tagDataLength + Constants::BAM_TAG_ARRAYBASE_SIZE;
    for ( int i = 0 ; i < numElements; ++i ) {
        const uint8_t value = values.at(i);
        memcpy(originalTagData + elementsBeginOffset + i*sizeof(uint8_t),
               &value, sizeof(uint8_t));
    }

    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData;
    TagData.assign(newTagData, newTagDataLength);

    delete[] originalTagData;

    // return success
    return true;
}

/*! \fn bool AddTag(const std::string& tag, const std::vector<int8_t>& values);
    \brief Adds a numeric array field to the BAM tags.

    Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.

    \param tag    2-character tag name
    \param values vector of int8_t values to store

    \return \c true if the \b new tag was added successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::AddTag(const std::string& tag, const std::vector<int8_t>& values) {

    // skip if core data not parsed
    if ( SupportData.HasCoreOnly ) return false;

    // check for valid tag length
    if ( tag.size() != Constants::BAM_TAG_TAGSIZE ) return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) )
        return false;

    // build new tag's base information
    char newTagBase[Constants::BAM_TAG_ARRAYBASE_SIZE];
    memcpy( newTagBase, tag.c_str(), Constants::BAM_TAG_TAGSIZE );
    newTagBase[2] = Constants::BAM_TAG_TYPE_ARRAY;
    newTagBase[3] = Constants::BAM_TAG_TYPE_INT8;

    // add number of array elements to newTagBase
    const int32_t numElements  = values.size();
    memcpy(newTagBase + 4, &numElements, sizeof(int32_t));

    // copy current TagData string to temp buffer, leaving room for new tag's contents
    const int newTagDataLength = tagDataLength +
                                 Constants::BAM_TAG_ARRAYBASE_SIZE +
                                 numElements*sizeof(int8_t);
    char* originalTagData = new char[newTagDataLength];
    memcpy(originalTagData, TagData.c_str(), tagDataLength+1); // '+1' for TagData's null-term

    // write newTagBase (removes old null term)
    strcat(originalTagData + tagDataLength, (const char*)newTagBase);

    // add vector elements to tag
    int elementsBeginOffset = tagDataLength + Constants::BAM_TAG_ARRAYBASE_SIZE;
    for ( int i = 0 ; i < numElements; ++i ) {
        const int8_t value = values.at(i);
        memcpy(originalTagData + elementsBeginOffset + i*sizeof(int8_t),
               &value, sizeof(int8_t));
    }

    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData;
    TagData.assign(newTagData, newTagDataLength);

    delete[] originalTagData;

    // return success
    return true;
}

/*! \fn bool AddTag(const std::string& tag, const std::vector<uint16_t>& values);
    \brief Adds a numeric array field to the BAM tags.

    Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.

    \param tag    2-character tag name
    \param values vector of uint16_t values to store

    \return \c true if the \b new tag was added successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::AddTag(const std::string& tag, const std::vector<uint16_t>& values) {

    // skip if core data not parsed
    if ( SupportData.HasCoreOnly ) return false;

    // check for valid tag length
    if ( tag.size() != Constants::BAM_TAG_TAGSIZE ) return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) )
        return false;

    // build new tag's base information
    char newTagBase[Constants::BAM_TAG_ARRAYBASE_SIZE];
    memcpy( newTagBase, tag.c_str(), Constants::BAM_TAG_TAGSIZE );
    newTagBase[2] = Constants::BAM_TAG_TYPE_ARRAY;
    newTagBase[3] = Constants::BAM_TAG_TYPE_UINT16;

    // add number of array elements to newTagBase
    const int32_t numElements  = values.size();
    memcpy(newTagBase + 4, &numElements, sizeof(int32_t));

    // copy current TagData string to temp buffer, leaving room for new tag's contents
    const int newTagDataLength = tagDataLength +
                                 Constants::BAM_TAG_ARRAYBASE_SIZE +
                                 numElements*sizeof(uint16_t);
    char* originalTagData = new char[newTagDataLength];
    memcpy(originalTagData, TagData.c_str(), tagDataLength+1); // '+1' for TagData's null-term

    // write newTagBase (removes old null term)
    strcat(originalTagData + tagDataLength, (const char*)newTagBase);

    // add vector elements to tag
    int elementsBeginOffset = tagDataLength + Constants::BAM_TAG_ARRAYBASE_SIZE;
    for ( int i = 0 ; i < numElements; ++i ) {
        const uint16_t value = values.at(i);
        memcpy(originalTagData + elementsBeginOffset + i*sizeof(uint16_t),
               &value, sizeof(uint16_t));
    }

    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData;
    TagData.assign(newTagData, newTagDataLength);

    delete[] originalTagData;

    // return success
    return true;
}

/*! \fn bool AddTag(const std::string& tag, const std::vector<int16_t>& values);
    \brief Adds a numeric array field to the BAM tags.

    Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.

    \param tag    2-character tag name
    \param values vector of int16_t values to store

    \return \c true if the \b new tag was added successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::AddTag(const std::string& tag, const std::vector<int16_t>& values) {

    // skip if core data not parsed
    if ( SupportData.HasCoreOnly ) return false;

    // check for valid tag length
    if ( tag.size() != Constants::BAM_TAG_TAGSIZE ) return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) )
        return false;

    // build new tag's base information
    char newTagBase[Constants::BAM_TAG_ARRAYBASE_SIZE];
    memcpy( newTagBase, tag.c_str(), Constants::BAM_TAG_TAGSIZE );
    newTagBase[2] = Constants::BAM_TAG_TYPE_ARRAY;
    newTagBase[3] = Constants::BAM_TAG_TYPE_INT16;

    // add number of array elements to newTagBase
    const int32_t numElements  = values.size();
    memcpy(newTagBase + 4, &numElements, sizeof(int32_t));

    // copy current TagData string to temp buffer, leaving room for new tag's contents
    const int newTagDataLength = tagDataLength +
                                 Constants::BAM_TAG_ARRAYBASE_SIZE +
                                 numElements*sizeof(int16_t);
    char* originalTagData = new char[newTagDataLength];
    memcpy(originalTagData, TagData.c_str(), tagDataLength+1); // '+1' for TagData's null-term

    // write newTagBase (removes old null term)
    strcat(originalTagData + tagDataLength, (const char*)newTagBase);

    // add vector elements to tag
    int elementsBeginOffset = tagDataLength + Constants::BAM_TAG_ARRAYBASE_SIZE;
    for ( int i = 0 ; i < numElements; ++i ) {
        const int16_t value = values.at(i);
        memcpy(originalTagData + elementsBeginOffset + i*sizeof(int16_t),
               &value, sizeof(int16_t));
    }

    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData;
    TagData.assign(newTagData, newTagDataLength);

    delete[] originalTagData;

    // return success
    return true;
}

/*! \fn bool AddTag(const std::string& tag, const std::vector<uint32_t>& values);
    \brief Adds a numeric array field to the BAM tags.

    Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.

    \param tag    2-character tag name
    \param values vector of uint32_t values to store

    \return \c true if the \b new tag was added successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::AddTag(const std::string& tag, const std::vector<uint32_t>& values) {

    // skip if core data not parsed
    if ( SupportData.HasCoreOnly ) return false;

    // check for valid tag length
    if ( tag.size() != Constants::BAM_TAG_TAGSIZE ) return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) )
        return false;

    // build new tag's base information
    char newTagBase[Constants::BAM_TAG_ARRAYBASE_SIZE];
    memcpy( newTagBase, tag.c_str(), Constants::BAM_TAG_TAGSIZE );
    newTagBase[2] = Constants::BAM_TAG_TYPE_ARRAY;
    newTagBase[3] = Constants::BAM_TAG_TYPE_UINT32;

    // add number of array elements to newTagBase
    const int32_t numElements  = values.size();
    memcpy(newTagBase + 4, &numElements, sizeof(int32_t));

    // copy current TagData string to temp buffer, leaving room for new tag's contents
    const int newTagDataLength = tagDataLength +
                                 Constants::BAM_TAG_ARRAYBASE_SIZE +
                                 numElements*sizeof(uint32_t);
    char* originalTagData = new char[newTagDataLength];
    memcpy(originalTagData, TagData.c_str(), tagDataLength+1); // '+1' for TagData's null-term

    // write newTagBase (removes old null term)
    strcat(originalTagData + tagDataLength, (const char*)newTagBase);

    // add vector elements to tag
    int elementsBeginOffset = tagDataLength + Constants::BAM_TAG_ARRAYBASE_SIZE;
    for ( int i = 0 ; i < numElements; ++i ) {
        const uint32_t value = values.at(i);
        memcpy(originalTagData + elementsBeginOffset + i*sizeof(uint32_t),
               &value, sizeof(uint32_t));
    }

    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData;
    TagData.assign(newTagData, newTagDataLength);

    delete[] originalTagData;

    // return success
    return true;
}

/*! \fn bool AddTag(const std::string& tag, const std::vector<int32_t>& values);
    \brief Adds a numeric array field to the BAM tags.

    Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.

    \param tag    2-character tag name
    \param values vector of int32_t values to store

    \return \c true if the \b new tag was added successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::AddTag(const std::string& tag, const std::vector<int32_t>& values) {

    // skip if core data not parsed
    if ( SupportData.HasCoreOnly ) return false;

    // check for valid tag length
    if ( tag.size() != Constants::BAM_TAG_TAGSIZE ) return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) )
        return false;

    // build new tag's base information
    char newTagBase[Constants::BAM_TAG_ARRAYBASE_SIZE];
    memcpy( newTagBase, tag.c_str(), Constants::BAM_TAG_TAGSIZE );
    newTagBase[2] = Constants::BAM_TAG_TYPE_ARRAY;
    newTagBase[3] = Constants::BAM_TAG_TYPE_INT32;

    // add number of array elements to newTagBase
    const int32_t numElements  = values.size();
    memcpy(newTagBase + 4, &numElements, sizeof(int32_t));

    // copy current TagData string to temp buffer, leaving room for new tag's contents
    const int newTagDataLength = tagDataLength +
                                 Constants::BAM_TAG_ARRAYBASE_SIZE +
                                 numElements*sizeof(int32_t);
    char* originalTagData = new char[newTagDataLength];
    memcpy(originalTagData, TagData.c_str(), tagDataLength+1); // '+1' for TagData's null-term

    // write newTagBase (removes old null term)
    strcat(originalTagData + tagDataLength, (const char*)newTagBase);

    // add vector elements to tag
    int elementsBeginOffset = tagDataLength + Constants::BAM_TAG_ARRAYBASE_SIZE;
    for ( int i = 0 ; i < numElements; ++i ) {
        const int32_t value = values.at(i);
        memcpy(originalTagData + elementsBeginOffset + i*sizeof(int32_t),
               &value, sizeof(int32_t));
    }

    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData;
    TagData.assign(newTagData, newTagDataLength);

    delete[] originalTagData;

    // return success
    return true;
}

/*! \fn bool AddTag(const std::string& tag, const std::vector<float>& values);
    \brief Adds a numeric array field to the BAM tags.

    Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.

    \param tag    2-character tag name
    \param values vector of float values to store

    \return \c true if the \b new tag was added successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::AddTag(const std::string& tag, const std::vector<float>& values) {

    // skip if core data not parsed
    if ( SupportData.HasCoreOnly ) return false;

    // check for valid tag length
    if ( tag.size() != Constants::BAM_TAG_TAGSIZE ) return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) )
        return false;

    // build new tag's base information
    char newTagBase[Constants::BAM_TAG_ARRAYBASE_SIZE];
    memcpy( newTagBase, tag.c_str(), Constants::BAM_TAG_TAGSIZE );
    newTagBase[2] = Constants::BAM_TAG_TYPE_ARRAY;
    newTagBase[3] = Constants::BAM_TAG_TYPE_FLOAT;

    // add number of array elements to newTagBase
    const int32_t numElements  = values.size();
    memcpy(newTagBase + 4, &numElements, sizeof(int32_t));

    // copy current TagData string to temp buffer, leaving room for new tag's contents
    const int newTagDataLength = tagDataLength +
                                 Constants::BAM_TAG_ARRAYBASE_SIZE +
                                 numElements*sizeof(float);
    char* originalTagData = new char[newTagDataLength];
    memcpy(originalTagData, TagData.c_str(), tagDataLength+1); // '+1' for TagData's null-term

    // write newTagBase (removes old null term)
    strcat(originalTagData + tagDataLength, (const char*)newTagBase);

    // add vector elements to tag
    int elementsBeginOffset = tagDataLength + Constants::BAM_TAG_ARRAYBASE_SIZE;
    for ( int i = 0 ; i < numElements; ++i ) {
        const float value = values.at(i);
        memcpy(originalTagData + elementsBeginOffset + i*sizeof(float),
               &value, sizeof(float));
    }

    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData;
    TagData.assign(newTagData, newTagDataLength);

    delete[] originalTagData;

    // return success
    return true;
}

/*! \fn bool BamAlignment::BuildCharData(void)
    \brief Populates alignment string fields (read name, bases, qualities, tag data).

    An alignment retrieved using BamReader::GetNextAlignmentCore() lacks this data.
    Using that method makes parsing much quicker when only positional data is required.

    However, if you later want to access the character data fields from such an alignment,
    use this method to populate those fields. Provides ability to do 'lazy evaluation' of
    alignment parsing.

    \return \c true if character data populated successfully (or was already available to begin with)
*/
bool BamAlignment::BuildCharData(void) {

    // skip if char data already parsed
    if ( !SupportData.HasCoreOnly )
        return true;

    // check system endianness
    bool IsBigEndian = BamTools::SystemIsBigEndian();

    // calculate character lengths/offsets
    const unsigned int dataLength     = SupportData.BlockLength - Constants::BAM_CORE_SIZE;
    const unsigned int seqDataOffset  = SupportData.QueryNameLength + (SupportData.NumCigarOperations * 4);
    const unsigned int qualDataOffset = seqDataOffset + (SupportData.QuerySequenceLength+1)/2;
    const unsigned int tagDataOffset  = qualDataOffset + SupportData.QuerySequenceLength;
    const unsigned int tagDataLength  = dataLength - tagDataOffset;

    // check offsets to see what char data exists
    const bool hasSeqData  = ( seqDataOffset  < dataLength );
    const bool hasQualData = ( qualDataOffset < dataLength );
    const bool hasTagData  = ( tagDataOffset  < dataLength );

    // set up char buffers
    const char* allCharData = SupportData.AllCharData.data();
    const char* seqData     = ( hasSeqData  ? (((const char*)allCharData) + seqDataOffset)  : (const char*)0 );
    const char* qualData    = ( hasQualData ? (((const char*)allCharData) + qualDataOffset) : (const char*)0 );
          char* tagData     = ( hasTagData  ? (((char*)allCharData) + tagDataOffset)        : (char*)0 );

    // store alignment name (relies on null char in name as terminator)
    Name.assign((const char*)(allCharData));

    // save query sequence
    QueryBases.clear();
    if ( hasSeqData ) {
        QueryBases.reserve(SupportData.QuerySequenceLength);
        for (unsigned int i = 0; i < SupportData.QuerySequenceLength; ++i) {
            char singleBase = Constants::BAM_DNA_LOOKUP[ ( (seqData[(i/2)] >> (4*(1-(i%2)))) & 0xf ) ];
            QueryBases.append(1, singleBase);
        }
    }

    // save qualities, converting from numeric QV to 'FASTQ-style' ASCII character
    Qualities.clear();
    if ( hasQualData ) {
        Qualities.reserve(SupportData.QuerySequenceLength);
        for (unsigned int i = 0; i < SupportData.QuerySequenceLength; ++i) {
            char singleQuality = (char)(qualData[i]+33);
            Qualities.append(1, singleQuality);
        }
    }

    // clear previous AlignedBases
    AlignedBases.clear();

    // if QueryBases has data, build AlignedBases using CIGAR data
    // otherwise, AlignedBases will remain empty (this case IS allowed)
    if ( !QueryBases.empty() ) {

        // resize AlignedBases
        AlignedBases.reserve(SupportData.QuerySequenceLength);

        // iterate over CigarOps
        int k = 0;
        vector<CigarOp>::const_iterator cigarIter = CigarData.begin();
        vector<CigarOp>::const_iterator cigarEnd  = CigarData.end();
        for ( ; cigarIter != cigarEnd; ++cigarIter ) {
            const CigarOp& op = (*cigarIter);

            switch (op.Type) {

                // for 'M', 'I', '=', 'X' - write bases
                case (Constants::BAM_CIGAR_MATCH_CHAR)    :
                case (Constants::BAM_CIGAR_INS_CHAR)      :
                case (Constants::BAM_CIGAR_SEQMATCH_CHAR) :
                case (Constants::BAM_CIGAR_MISMATCH_CHAR) :
                    AlignedBases.append(QueryBases.substr(k, op.Length));
                    // fall through

                // for 'S' - soft clip, do not write bases
                // but increment placeholder 'k'
                case (Constants::BAM_CIGAR_SOFTCLIP_CHAR) :
                    k += op.Length;
                    break;

                // for 'D' - write gap character
                case (Constants::BAM_CIGAR_DEL_CHAR) :
                    AlignedBases.append(op.Length, Constants::BAM_DNA_DEL);
                    break;

                // for 'P' - write padding character
                case (Constants::BAM_CIGAR_PAD_CHAR) :
                    AlignedBases.append( op.Length, Constants::BAM_DNA_PAD );
                    break;

                // for 'N' - write N's, skip bases in original query sequence
                case (Constants::BAM_CIGAR_REFSKIP_CHAR) :
                    AlignedBases.append( op.Length, Constants::BAM_DNA_N );
                    break;

                // for 'H' - hard clip, do nothing to AlignedBases, move to next op
                case (Constants::BAM_CIGAR_HARDCLIP_CHAR) :
                    break;

                // shouldn't get here
                default:
                    cerr << "BamAlignment ERROR: invalid CIGAR operation type: "
                         << op.Type << endl;
                    exit(1);
            }
        }
    }

    // save tag data
    TagData.clear();
    if ( hasTagData ) {
        if ( IsBigEndian ) {
            int i = 0;
            while ( (unsigned int)i < tagDataLength ) {

                i += Constants::BAM_TAG_TAGSIZE;  // skip tag chars (e.g. "RG", "NM", etc.)
                const char type = tagData[i];     // get tag type at position i
                ++i;                              // move i past tag type

                switch (type) {

                    case(Constants::BAM_TAG_TYPE_ASCII) :
                    case(Constants::BAM_TAG_TYPE_INT8)  :
                    case(Constants::BAM_TAG_TYPE_UINT8) :
                        // no endian swapping necessary for single-byte data
                        ++i;
                        break;

                    case(Constants::BAM_TAG_TYPE_INT16)  :
                    case(Constants::BAM_TAG_TYPE_UINT16) :
                        BamTools::SwapEndian_16p(&tagData[i]);
                        i += sizeof(uint16_t);
                        break;

                    case(Constants::BAM_TAG_TYPE_FLOAT)  :
                    case(Constants::BAM_TAG_TYPE_INT32)  :
                    case(Constants::BAM_TAG_TYPE_UINT32) :
                        BamTools::SwapEndian_32p(&tagData[i]);
                        i += sizeof(uint32_t);
                        break;

                    case(Constants::BAM_TAG_TYPE_HEX) :
                    case(Constants::BAM_TAG_TYPE_STRING) :
                        // no endian swapping necessary for hex-string/string data
                        while ( tagData[i] )
                            ++i;
                        // increment one more for null terminator
                        ++i;
                        break;

                    case(Constants::BAM_TAG_TYPE_ARRAY) :

                    {
                        // read array type
                        const char arrayType = tagData[i];
                        ++i;

                        // swap endian-ness of number of elements in place, then retrieve for loop
                        BamTools::SwapEndian_32p(&tagData[i]);
                        int32_t numElements;
                        memcpy(&numElements, &tagData[i], sizeof(uint32_t));
                        i += sizeof(uint32_t);

                        // swap endian-ness of array elements
                        for ( int j = 0; j < numElements; ++j ) {
                            switch (arrayType) {
                                case (Constants::BAM_TAG_TYPE_INT8)  :
                                case (Constants::BAM_TAG_TYPE_UINT8) :
                                    // no endian-swapping necessary
                                    ++i;
                                    break;
                                case (Constants::BAM_TAG_TYPE_INT16)  :
                                case (Constants::BAM_TAG_TYPE_UINT16) :
                                    BamTools::SwapEndian_16p(&tagData[i]);
                                    i += sizeof(uint16_t);
                                    break;
                                case (Constants::BAM_TAG_TYPE_FLOAT)  :
                                case (Constants::BAM_TAG_TYPE_INT32)  :
                                case (Constants::BAM_TAG_TYPE_UINT32) :
                                    BamTools::SwapEndian_32p(&tagData[i]);
                                    i += sizeof(uint32_t);
                                    break;
                                default:
                                    // error case
                                    cerr << "BamAlignment ERROR: unknown binary array type encountered: "
                                         << arrayType << endl;
                                    return false;
                            }
                        }

                        break;
                    }

                    // shouldn't get here
                    default :
                        cerr << "BamAlignment ERROR: invalid tag value type: "
                             << type << endl;
                        exit(1);
                }
            }
        }

        // store tagData in alignment
        TagData.resize(tagDataLength);
        memcpy((char*)TagData.data(), tagData, tagDataLength);
    }

    // clear the core-only flag
    SupportData.HasCoreOnly = false;

    // return success
    return true;
}

/*! \fn bool BamAlignment::EditTag(const std::string& tag, const std::string& type, const std::string& value)
    \brief Edits a BAM tag field containing string data.

    If \a tag does not exist, a new entry is created.

    \param tag   2-character tag name
    \param type  1-character tag type (must be "Z" or "H")
    \param value string data to store

    \return \c true if the tag was modified/created successfully

    \sa BamAlignment::RemoveTag()
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::EditTag(const std::string& tag, const std::string& type, const std::string& value) {
  
    // skip if core data not parsed
    if ( SupportData.HasCoreOnly ) return false;

    // validate tag/type size & that type is OK for string value
    if ( !IsValidSize(tag, type) ) return false;
    if ( type.at(0) != Constants::BAM_TAG_TYPE_STRING &&
         type.at(0) != Constants::BAM_TAG_TYPE_HEX )
        return false;
  
    // localize the tag data
    char* pOriginalTagData = (char*)TagData.data();
    char* pTagData = pOriginalTagData;
    const unsigned int originalTagDataLength = TagData.size();
    
    unsigned int newTagDataLength = 0;
    unsigned int numBytesParsed = 0;
    
    // if tag found
    if ( FindTag(tag, pTagData, originalTagDataLength, numBytesParsed) ) {
        
        // make sure array is more than big enough
        char* newTagData = new char[originalTagDataLength + value.size()];  

        // copy original tag data up til desired tag
        const unsigned int beginningTagDataLength = numBytesParsed;
        newTagDataLength += beginningTagDataLength;
        memcpy(newTagData, pOriginalTagData, numBytesParsed);
      
        // copy new @value in place of current tag data
        const unsigned int dataLength = strlen(value.c_str());
        memcpy(newTagData + beginningTagDataLength, (char*)value.c_str(), dataLength+1 );
        
        // skip to next tag (if tag for removal is last, return true) 
        const char* pTagStorageType = pTagData - 1;
        if ( !SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) )
            return true;
         
        // copy everything from current tag (the next one after tag for removal) to end
        const unsigned int skippedDataLength = (numBytesParsed - beginningTagDataLength);
        const unsigned int endTagOffset      = beginningTagDataLength + dataLength + 1;
        const unsigned int endTagDataLength  = originalTagDataLength - beginningTagDataLength - skippedDataLength;
        memcpy(newTagData + endTagOffset, pTagData, endTagDataLength);
        
        // ensure null-terminator
        newTagData[ endTagOffset + endTagDataLength + 1 ] = 0;
        
        // save new tag data
        TagData.assign(newTagData, endTagOffset + endTagDataLength);

        delete[] newTagData;

        return true;
    }
    
    // tag not found, attempt AddTag
    else return AddTag(tag, type, value);
}

/*! \fn bool BamAlignment::EditTag(const std::string& tag, const std::string& type, const uint32_t& value)
    \brief Edits a BAM tag field containing unsigned integer data.

    If \a tag does not exist, a new entry is created.

    \param tag   2-character tag name
    \param type  1-character tag type (must NOT be "f", "Z", "H", or "B")
    \param value unsigned integer data to store

    \return \c true if the tag was modified/created successfully

    \sa BamAlignment::RemoveTag()
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::EditTag(const std::string& tag, const std::string& type, const uint32_t& value) {
  
    // skip if core data not parsed
    if ( SupportData.HasCoreOnly ) return false;

    // validate tag/type size & that type is OK for uint32_t value
    if ( !IsValidSize(tag, type) ) return false;
    if ( type.at(0) == Constants::BAM_TAG_TYPE_FLOAT  ||
         type.at(0) == Constants::BAM_TAG_TYPE_STRING ||
         type.at(0) == Constants::BAM_TAG_TYPE_HEX    ||
         type.at(0) == Constants::BAM_TAG_TYPE_ARRAY
       )
    {
        return false;
    }

    // localize the tag data
    char* pOriginalTagData = (char*)TagData.data();
    char* pTagData = pOriginalTagData;
    const unsigned int originalTagDataLength = TagData.size();
    
    unsigned int newTagDataLength = 0;
    unsigned int numBytesParsed = 0;
    
    // if tag found
    if ( FindTag(tag, pTagData, originalTagDataLength, numBytesParsed) ) {
        
        // make sure array is more than big enough
        char* newTagData = new char[originalTagDataLength + sizeof(value)];  

        // copy original tag data up til desired tag
        const unsigned int beginningTagDataLength = numBytesParsed;
        newTagDataLength += beginningTagDataLength;
        memcpy(newTagData, pOriginalTagData, numBytesParsed);
      
        // copy new @value in place of current tag data
        union { uint32_t value; char valueBuffer[sizeof(uint32_t)]; } un;
        un.value = value;
        memcpy(newTagData + beginningTagDataLength, un.valueBuffer, sizeof(uint32_t));
        
        // skip to next tag (if tag for removal is last, return true) 
        const char* pTagStorageType = pTagData - 1;
        if ( !SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) )
            return true;
         
        // copy everything from current tag (the next one after tag for removal) to end
        const unsigned int skippedDataLength = (numBytesParsed - beginningTagDataLength);
        const unsigned int endTagOffset      = beginningTagDataLength + sizeof(uint32_t);
        const unsigned int endTagDataLength  = originalTagDataLength - beginningTagDataLength - skippedDataLength;
        memcpy(newTagData + endTagOffset, pTagData, endTagDataLength);
        
        // ensure null-terminator
        newTagData[ endTagOffset + endTagDataLength + 1 ] = 0;
        
        // save new tag data
        TagData.assign(newTagData, endTagOffset + endTagDataLength);

        delete[] newTagData;

        return true;
    }
    
    // tag not found, attempt AddTag
    else return AddTag(tag, type, value);
}

/*! \fn bool BamAlignment::EditTag(const std::string& tag, const std::string& type, const int32_t& value)
    \brief Edits a BAM tag field containing signed integer data.

    If \a tag does not exist, a new entry is created.

    \param tag   2-character tag name
    \param type  1-character tag type (must NOT be "f", "Z", "H", or "B")
    \param value signed integer data to store

    \return \c true if the tag was modified/created successfully

    \sa BamAlignment::RemoveTag()
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::EditTag(const std::string& tag, const std::string& type, const int32_t& value) {
    return EditTag(tag, type, (const uint32_t&)value);
}

/*! \fn bool BamAlignment::EditTag(const std::string& tag, const std::string& type, const float& value)
    \brief Edits a BAM tag field containing floating-point data.

    If \a tag does not exist, a new entry is created.

    \param tag   2-character tag name
    \param type  1-character tag type (must NOT be "Z", "H", or "B")
    \param value float data to store

    \return \c true if the tag was modified/created successfully

    \sa BamAlignment::RemoveTag()
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::EditTag(const std::string& tag, const std::string& type, const float& value) {
  
    // skip if core data not parsed
    if ( SupportData.HasCoreOnly ) return false;

    // validate tag/type size & that type is OK for float value
    if ( !IsValidSize(tag, type) ) return false;
    if ( type.at(0) == Constants::BAM_TAG_TYPE_STRING ||
         type.at(0) == Constants::BAM_TAG_TYPE_HEX    ||
         type.at(0) == Constants::BAM_TAG_TYPE_ARRAY
       )
    {
        return false;
    }

     // localize the tag data
    char* pOriginalTagData = (char*)TagData.data();
    char* pTagData = pOriginalTagData;
    const unsigned int originalTagDataLength = TagData.size();
    
    unsigned int newTagDataLength = 0;
    unsigned int numBytesParsed = 0;
    
    // if tag found
    if ( FindTag(tag, pTagData, originalTagDataLength, numBytesParsed) ) {
        
        // make sure array is more than big enough
        char* newTagData = new char[originalTagDataLength + sizeof(value)];  

        // copy original tag data up til desired tag
        const unsigned int beginningTagDataLength = numBytesParsed;
        newTagDataLength += beginningTagDataLength;
        memcpy(newTagData, pOriginalTagData, numBytesParsed);
      
        // copy new @value in place of current tag data
        union { float value; char valueBuffer[sizeof(float)]; } un;
        un.value = value;
        memcpy(newTagData + beginningTagDataLength, un.valueBuffer, sizeof(float));
        
        // skip to next tag (if tag for removal is last, return true) 
        const char* pTagStorageType = pTagData - 1;
        if ( !SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) )
            return true;
         
        // copy everything from current tag (the next one after tag for removal) to end
        const unsigned int skippedDataLength = (numBytesParsed - beginningTagDataLength);
        const unsigned int endTagOffset      = beginningTagDataLength + sizeof(float);
        const unsigned int endTagDataLength  = originalTagDataLength - beginningTagDataLength - skippedDataLength;
        memcpy(newTagData + endTagOffset, pTagData, endTagDataLength);
        
        // ensure null-terminator
        newTagData[ endTagOffset + endTagDataLength + 1 ] = 0;
        
        // save new tag data
        TagData.assign(newTagData, endTagOffset + endTagDataLength);

        delete[] newTagData;

        return true;
    }
    
    // tag not found, attempt AddTag
    else return AddTag(tag, type, value);
}

/*! \fn bool EditTag(const std::string& tag, const std::vector<uint8_t>& values);
    \brief Edits a BAM tag field containing a numeric array.

    If \a tag does not exist, a new entry is created.

    \param tag   2-character tag name
    \param value vector of uint8_t values to store

    \return \c true if the tag was modified/created successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::EditTag(const std::string& tag, const std::vector<uint8_t>& values) {

    // can't do anything if TagData not parsed
    if ( SupportData.HasCoreOnly )
        return false;

    // remove existing tag if present
    if ( HasTag(tag) )
        RemoveTag(tag);

    // add tag record with new values
    return AddTag(tag, values);
}

/*! \fn bool EditTag(const std::string& tag, const std::vector<int8_t>& values);
    \brief Edits a BAM tag field containing a numeric array.

    If \a tag does not exist, a new entry is created.

    \param tag   2-character tag name
    \param value vector of int8_t values to store

    \return \c true if the tag was modified/created successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::EditTag(const std::string& tag, const std::vector<int8_t>& values) {

    // can't do anything if TagData not parsed
    if ( SupportData.HasCoreOnly )
        return false;

    // remove existing tag if present
    if ( HasTag(tag) )
        RemoveTag(tag);

    // add tag record with new values
    return AddTag(tag, values);
}

/*! \fn bool EditTag(const std::string& tag, const std::vector<uint16_t>& values);
    \brief Edits a BAM tag field containing a numeric array.

    If \a tag does not exist, a new entry is created.

    \param tag   2-character tag name
    \param value vector of uint16_t values to store

    \return \c true if the tag was modified/created successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::EditTag(const std::string& tag, const std::vector<uint16_t>& values) {

    // can't do anything if TagData not parsed
    if ( SupportData.HasCoreOnly )
        return false;

    // remove existing tag if present
    if ( HasTag(tag) )
        RemoveTag(tag);

    // add tag record with new values
    return AddTag(tag, values);
}

/*! \fn bool EditTag(const std::string& tag, const std::vector<int16_t>& values);
    \brief Edits a BAM tag field containing a numeric array.

    If \a tag does not exist, a new entry is created.

    \param tag   2-character tag name
    \param value vector of int16_t values to store

    \return \c true if the tag was modified/created successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::EditTag(const std::string& tag, const std::vector<int16_t>& values) {

    // can't do anything if TagData not parsed
    if ( SupportData.HasCoreOnly )
        return false;

    // remove existing tag if present
    if ( HasTag(tag) )
        RemoveTag(tag);

    // add tag record with new values
    return AddTag(tag, values);
}

/*! \fn bool EditTag(const std::string& tag, const std::vector<uint32_t>& values);
    \brief Edits a BAM tag field containing a numeric array.

    If \a tag does not exist, a new entry is created.

    \param tag   2-character tag name
    \param value vector of uint32_t values to store

    \return \c true if the tag was modified/created successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::EditTag(const std::string& tag, const std::vector<uint32_t>& values) {

    // can't do anything if TagData not parsed
    if ( SupportData.HasCoreOnly )
        return false;

    // remove existing tag if present
    if ( HasTag(tag) )
        RemoveTag(tag);

    // add tag record with new values
    return AddTag(tag, values);
}

/*! \fn bool EditTag(const std::string& tag, const std::vector<int32_t>& values);
    \brief Edits a BAM tag field containing a numeric array.

    If \a tag does not exist, a new entry is created.

    \param tag   2-character tag name
    \param value vector of int32_t values to store

    \return \c true if the tag was modified/created successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::EditTag(const std::string& tag, const std::vector<int32_t>& values) {

    // can't do anything if TagData not parsed
    if ( SupportData.HasCoreOnly )
        return false;

    // remove existing tag if present
    if ( HasTag(tag) )
        RemoveTag(tag);

    // add tag record with new values
    return AddTag(tag, values);
}

/*! \fn bool EditTag(const std::string& tag, const std::vector<float>& values);
    \brief Edits a BAM tag field containing a numeric array.

    If \a tag does not exist, a new entry is created.

    \param tag   2-character tag name
    \param value vector of float values to store

    \return \c true if the tag was modified/created successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::EditTag(const std::string& tag, const std::vector<float>& values) {

    // can't do anything if TagData not parsed
    if ( SupportData.HasCoreOnly )
        return false;

    // remove existing tag if present
    if ( HasTag(tag) )
        RemoveTag(tag);

    // add tag record with new values
    return AddTag(tag, values);
}

/*! \fn bool BamAlignment::FindTag(const std::string& tag, char*& pTagData, const unsigned int& tagDataLength, unsigned int& numBytesParsed)
    \internal

    Searches for requested tag in BAM tag data.

    \param tag            requested 2-character tag name
    \param pTagData       pointer to current position in BamAlignment::TagData
    \param tagDataLength  length of BamAlignment::TagData
    \param numBytesParsed number of bytes parsed so far

    \return \c true if found

    \post If \a tag is found, \a pTagData will point to the byte where the tag data begins.
          \a numBytesParsed will correspond to the position in the full TagData string.

*/
bool BamAlignment::FindTag(const std::string& tag,
                           char*& pTagData,
                           const unsigned int& tagDataLength,
                           unsigned int& numBytesParsed) const
{

    while ( numBytesParsed < tagDataLength ) {

        const char* pTagType        = pTagData;
        const char* pTagStorageType = pTagData + 2;
        pTagData       += 3;
        numBytesParsed += 3;

        // check the current tag, return true on match
        if ( strncmp(pTagType, tag.c_str(), 2) == 0 )
            return true;

        // get the storage class and find the next tag
        if ( *pTagStorageType == '\0' ) return false;
        if ( !SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) ) return false;
        if ( *pTagData == '\0' ) return false;
    }

    // checked all tags, none match
    return false;
}

/*! \fn bool BamAlignment::GetEditDistance(uint32_t& editDistance) const
    \brief Retrieves value of edit distance tag ("NM").

    \deprecated Instead use BamAlignment::GetTag()
        \code
            BamAlignment::GetTag("NM", editDistance);
        \endcode

    \param editDistance destination for retrieved value

    \return \c true if found
*/
bool BamAlignment::GetEditDistance(uint32_t& editDistance) const {
    return GetTag("NM", (uint32_t&)editDistance);
}

/*! \fn int BamAlignment::GetEndPosition(bool usePadded = false, bool zeroBased = true) const
    \brief Calculates alignment end position, based on starting position and CIGAR data.

    \param usePadded Inserted bases affect reported position. Default is false, so that reported
                     position stays 'sync-ed' with reference coordinates.
    \param zeroBased Return (BAM standard) 0-based coordinate. Setting this to false can be useful
                     when using BAM data with half-open formats (e.g. BED).

    \return alignment end position
*/
int BamAlignment::GetEndPosition(bool usePadded, bool zeroBased) const {

    // initialize alignment end to starting position
    int alignEnd = Position;

    // iterate over cigar operations
    vector<CigarOp>::const_iterator cigarIter = CigarData.begin();
    vector<CigarOp>::const_iterator cigarEnd  = CigarData.end();
    for ( ; cigarIter != cigarEnd; ++cigarIter) {
        const char cigarType = (*cigarIter).Type;
        const uint32_t& cigarLength = (*cigarIter).Length;

        if ( cigarType == Constants::BAM_CIGAR_MATCH_CHAR ||
             cigarType == Constants::BAM_CIGAR_DEL_CHAR ||
             cigarType == Constants::BAM_CIGAR_REFSKIP_CHAR )
            alignEnd += cigarLength;
        else if ( usePadded && cigarType == Constants::BAM_CIGAR_INS_CHAR )
            alignEnd += cigarLength;
    }

    // adjust for zero-based coordinates, if requested
    if ( zeroBased ) alignEnd -= 1;

    // return result
    return alignEnd;
}

/*! \fn bool BamAlignment::GetReadGroup(std::string& readGroup) const
    \brief Retrieves value of read group tag ("RG").

    \deprecated Instead use BamAlignment::GetTag()
        \code
            BamAlignment::GetTag("RG", readGroup);
        \endcode

    \param readGroup destination for retrieved value

    \return \c true if found
*/
bool BamAlignment::GetReadGroup(std::string& readGroup) const {
    return GetTag("RG", readGroup);
}

/*! \fn bool BamAlignment::GetTag(const std::string& tag, std::string& destination) const
    \brief Retrieves the string value associated with a BAM tag.

    \param tag         2-character tag name
    \param destination destination for retrieved value

    \return \c true if found
*/
bool BamAlignment::GetTag(const std::string& tag, std::string& destination) const {

    // make sure tag data exists
    if ( SupportData.HasCoreOnly || TagData.empty() ) 
        return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag found
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

/*! \fn bool BamAlignment::GetTag(const std::string& tag, uint32_t& destination) const
    \brief Retrieves the unsigned integer value associated with a BAM tag.

    \param tag         2-character tag name
    \param destination destination for retrieved value

    \return \c true if found
*/
bool BamAlignment::GetTag(const std::string& tag, uint32_t& destination) const {
  
    // make sure tag data exists
    if ( SupportData.HasCoreOnly || TagData.empty() ) 
        return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag found
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        
        // determine data byte-length
        const char type = *(pTagData - 1);
        int destinationLength = 0;
        switch (type) {

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
                destinationLength = 4;
                break;

            // unsupported type for integer destination (float or var-length strings)
            case (Constants::BAM_TAG_TYPE_FLOAT)  :
            case (Constants::BAM_TAG_TYPE_STRING) :
            case (Constants::BAM_TAG_TYPE_HEX)    :
            case (Constants::BAM_TAG_TYPE_ARRAY)  :
                cerr << "BamAlignment ERROR: cannot store tag of type " << type
                     << " in integer destination" << endl;
                return false;

            // unknown tag type
            default:
                cerr << "BamAlignment ERROR: unknown tag type encountered: "
                     << type << endl;
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

/*! \fn bool BamAlignment::GetTag(const std::string& tag, int32_t& destination) const
    \brief Retrieves the signed integer value associated with a BAM tag.

    \param tag         2-character tag name
    \param destination destination for retrieved value

    \return \c true if found
*/
bool BamAlignment::GetTag(const std::string& tag, int32_t& destination) const {
    return GetTag(tag, (uint32_t&)destination);
}

/*! \fn bool BamAlignment::GetTag(const std::string& tag, float& destination) const
    \brief Retrieves the floating-point value associated with a BAM tag.

    \param tag         2-character tag name
    \param destination destination for retrieved value

    \return \c true if found
*/
bool BamAlignment::GetTag(const std::string& tag, float& destination) const {
  
    // make sure tag data exists
    if ( SupportData.HasCoreOnly || TagData.empty() ) 
        return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag found
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        
        // determine data byte-length
        const char type = *(pTagData - 1);
        int destinationLength = 0;
        switch (type) {

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
            case (Constants::BAM_TAG_TYPE_FLOAT)  :
            case (Constants::BAM_TAG_TYPE_INT32)  :
            case (Constants::BAM_TAG_TYPE_UINT32) :
                destinationLength = 4;
                break;
            
            // unsupported type (var-length strings)
            case (Constants::BAM_TAG_TYPE_STRING) :
            case (Constants::BAM_TAG_TYPE_HEX)    :
            case (Constants::BAM_TAG_TYPE_ARRAY)  :
                cerr << "BamAlignment ERROR: cannot store tag of type " << type
                     << " in float destination" << endl;
                return false;

            // unknown tag type
            default:
                cerr << "BamAlignment ERROR: unknown tag type encountered: "
                     << type << endl;
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

/*! \fn bool BamAlignment::GetTag(const std::string& tag, std::vector<uint32_t>& destination) const
    \brief Retrieves the numeric array data associated with a BAM tag

    \param tag         2-character tag name
    \param destination destination for retrieved data

    \return \c true if found
*/
bool BamAlignment::GetTag(const std::string& tag, std::vector<uint32_t>& destination) const {

    // make sure tag data exists
    if ( SupportData.HasCoreOnly || TagData.empty() )
        return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // return false if tag not found
    if ( !FindTag(tag, pTagData, tagDataLength, numBytesParsed) )
        return false;

    // check that tag is array type
    const char tagType = *(pTagData - 1);
    if ( tagType != Constants::BAM_TAG_TYPE_ARRAY ) {
        cerr << "BamAlignment ERROR: Cannot store non-array data from tag: "
             << tag << " in array destination" << endl;
        return false;
    }

    // calculate length of each element in tag's array
    const char elementType = *pTagData;
    ++pTagData;
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
            elementLength = sizeof(uint32_t);
            break;

        // unsupported type for integer destination (float or var-length data)
        case (Constants::BAM_TAG_TYPE_FLOAT)  :
        case (Constants::BAM_TAG_TYPE_STRING) :
        case (Constants::BAM_TAG_TYPE_HEX)    :
        case (Constants::BAM_TAG_TYPE_ARRAY)  :
            cerr << "BamAlignment ERROR: array element type: " << elementType
                 << " cannot be stored in integer value" << endl;
            return false;

        // unknown tag type
        default:
            cerr << "BamAlignment ERROR: unknown element type encountered: "
                 << elementType << endl;
            return false;
    }

    // get number of elements
    int32_t numElements;
    memcpy(&numElements, pTagData, sizeof(int32_t));
    pTagData += 4;
    destination.clear();
    destination.reserve(numElements);

    // read in elements
    uint32_t value;
    for ( int i = 0 ; i < numElements; ++i ) {
        memcpy(&value, pTagData, sizeof(uint32_t));
        pTagData += sizeof(uint32_t);
        destination.push_back(value);
    }

    // return success
    return false;
}

/*! \fn bool BamAlignment::GetTag(const std::string& tag, std::vector<int32_t>& destination) const
    \brief Retrieves the numeric array data associated with a BAM tag

    \param tag         2-character tag name
    \param destination destination for retrieved data

    \return \c true if found
*/
bool BamAlignment::GetTag(const std::string& tag, std::vector<int32_t>& destination) const {

    // make sure tag data exists
    if ( SupportData.HasCoreOnly || TagData.empty() )
        return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // return false if tag not found
    if ( !FindTag(tag, pTagData, tagDataLength, numBytesParsed) )
        return false;

    // check that tag is array type
    const char tagType = *(pTagData - 1);
    if ( tagType != Constants::BAM_TAG_TYPE_ARRAY ) {
        cerr << "BamAlignment ERROR: Cannot store non-array data from tag: "
             << tag << " in array destination" << endl;
        return false;
    }

    // calculate length of each element in tag's array
    const char elementType = *pTagData;
    ++pTagData;
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
            elementLength = sizeof(uint32_t);
            break;

        // unsupported type for integer destination (float or var-length data)
        case (Constants::BAM_TAG_TYPE_FLOAT)  :
        case (Constants::BAM_TAG_TYPE_STRING) :
        case (Constants::BAM_TAG_TYPE_HEX)    :
        case (Constants::BAM_TAG_TYPE_ARRAY)  :
            cerr << "BamAlignment ERROR: array element type: " << elementType
                 << " cannot be stored in integer value" << endl;
            return false;

        // unknown tag type
        default:
            cerr << "BamAlignment ERROR: unknown element type encountered: "
                 << elementType << endl;
            return false;
    }

    // get number of elements
    int32_t numElements;
    memcpy(&numElements, pTagData, sizeof(int32_t));
    pTagData += 4;
    destination.clear();
    destination.reserve(numElements);

    // read in elements
    int32_t value;
    for ( int i = 0 ; i < numElements; ++i ) {
        memcpy(&value, pTagData, sizeof(int32_t));
        pTagData += sizeof(int32_t);
        destination.push_back(value);
    }

    // return success
    return false;

}

/*! \fn bool BamAlignment::GetTag(const std::string& tag, std::vector<float>& destination) const
    \brief Retrieves the numeric array data associated with a BAM tag

    \param tag         2-character tag name
    \param destination destination for retrieved data

    \return \c true if found
*/
bool BamAlignment::GetTag(const std::string& tag, std::vector<float>& destination) const {

    // make sure tag data exists
    if ( SupportData.HasCoreOnly || TagData.empty() )
        return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // return false if tag not found
    if ( !FindTag(tag, pTagData, tagDataLength, numBytesParsed) )
        return false;

    // check that tag is array type
    const char tagType = *(pTagData - 1);
    if ( tagType != Constants::BAM_TAG_TYPE_ARRAY ) {
        cerr << "BamAlignment ERROR: Cannot store non-array data from tag: "
             << tag << " in array destination" << endl;
        return false;
    }

    // calculate length of each element in tag's array
    const char elementType = *pTagData;
    ++pTagData;
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

        // unsupported type for float destination (var-length data)
        case (Constants::BAM_TAG_TYPE_STRING) :
        case (Constants::BAM_TAG_TYPE_HEX)    :
        case (Constants::BAM_TAG_TYPE_ARRAY)  :
            cerr << "BamAlignment ERROR: array element type: " << elementType
                 << " cannot be stored in float value" << endl;
            return false;

        // unknown tag type
        default:
            cerr << "BamAlignment ERROR: unknown element type encountered: "
                 << elementType << endl;
            return false;
    }

    // get number of elements
    int32_t numElements;
    memcpy(&numElements, pTagData, sizeof(int32_t));
    pTagData += 4;
    destination.clear();
    destination.reserve(numElements);

    // read in elements
    float value;
    for ( int i = 0 ; i < numElements; ++i ) {
        memcpy(&value, pTagData, sizeof(float));
        pTagData += sizeof(float);
        destination.push_back(value);
    }

    // return success
    return false;
}

/*! \fn bool BamAlignment::GetTagType(const std::string& tag, char& type) const
    \brief Retrieves the BAM tag type-code associated with requested tag name.

    \param tag  2-character tag name
    \param type destination for the retrieved (1-character) tag type

    \return \c true if found
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::GetTagType(const std::string& tag, char& type) const {
  
    // make sure tag data exists
    if ( SupportData.HasCoreOnly || TagData.empty() ) 
        return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // lookup tag
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        
        // retrieve tag type code
        type = *(pTagData - 1);
        
        // validate that type is a proper BAM tag type
        switch (type) {
            case (Constants::BAM_TAG_TYPE_ASCII)  :
            case (Constants::BAM_TAG_TYPE_INT8)   :
            case (Constants::BAM_TAG_TYPE_UINT8)  :
            case (Constants::BAM_TAG_TYPE_INT16)  :
            case (Constants::BAM_TAG_TYPE_UINT16) :
            case (Constants::BAM_TAG_TYPE_INT32)  :
            case (Constants::BAM_TAG_TYPE_UINT32) :
            case (Constants::BAM_TAG_TYPE_FLOAT)  :
            case (Constants::BAM_TAG_TYPE_STRING) :
            case (Constants::BAM_TAG_TYPE_HEX)    :
            case (Constants::BAM_TAG_TYPE_ARRAY)  :
                return true;

            // unknown tag type
            default:
                cerr << "BamAlignment ERROR: unknown tag type encountered: "
                     << type << endl;
                return false;
        }
    }
    
    // tag not found, return failure
    return false;
}

/*! \fn bool BamAlignment::HasTag(const std::string& tag) const
    \brief Returns true if alignment has a record for requested tag.
    \param tag 2-character tag name
    \return \c true if alignment has a record for tag
*/
bool BamAlignment::HasTag(const std::string& tag) const {

    // return false if no tag data present
    if ( SupportData.HasCoreOnly || TagData.empty() )
        return false;

    // localize the tag data for lookup
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if result of tag lookup
    return FindTag(tag, pTagData, tagDataLength, numBytesParsed);
}

/*! \fn bool BamAlignment::IsDuplicate(void) const
    \return \c true if this read is a PCR duplicate
*/
bool BamAlignment::IsDuplicate(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_DUPLICATE) != 0 );
}

/*! \fn bool BamAlignment::IsFailedQC(void) const
    \return \c true if this read failed quality control
*/
bool BamAlignment::IsFailedQC(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_QC_FAILED) != 0 );
}

/*! \fn bool BamAlignment::IsFirstMate(void) const
    \return \c true if alignment is first mate on paired-end read
*/
bool BamAlignment::IsFirstMate(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_READ_1) != 0 );
}

/*! \fn bool BamAlignment::IsMapped(void) const
    \return \c true if alignment is mapped
*/
bool BamAlignment::IsMapped(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_UNMAPPED) == 0 );
}

/*! \fn bool BamAlignment::IsMateMapped(void) const
    \return \c true if alignment's mate is mapped
*/
bool BamAlignment::IsMateMapped(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_MATE_UNMAPPED) == 0 );
}

/*! \fn bool BamAlignment::IsMateReverseStrand(void) const
    \return \c true if alignment's mate mapped to reverse strand
*/
bool BamAlignment::IsMateReverseStrand(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_MATE_REVERSE_STRAND) != 0 );
}

/*! \fn bool BamAlignment::IsPaired(void) const
    \return \c true if alignment part of paired-end read
*/
bool BamAlignment::IsPaired(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_PAIRED) != 0 );
}

/*! \fn bool BamAlignment::IsPrimaryAlignment(void) const
    \return \c true if reported position is primary alignment
*/
bool BamAlignment::IsPrimaryAlignment(void) const  {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_SECONDARY) == 0 );
}

/*! \fn bool BamAlignment::IsProperPair(void) const
    \return \c true if alignment is part of read that satisfied paired-end resolution
*/
bool BamAlignment::IsProperPair(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_PROPER_PAIR) != 0 );
}

/*! \fn bool BamAlignment::IsReverseStrand(void) const
    \return \c true if alignment mapped to reverse strand
*/
bool BamAlignment::IsReverseStrand(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_REVERSE_STRAND) != 0 );
}

/*! \fn bool BamAlignment::IsSecondMate(void) const
    \return \c true if alignment is second mate on read
*/
bool BamAlignment::IsSecondMate(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_READ_2) != 0 );
}

/*! \fn bool BamAlignment::IsValidSize(const string& tag, const string& type) const
    \internal

    Checks that tag name & type strings are expected sizes.
    \a tag  should have length
    \a type should have length 1

    \param tag  BAM tag name
    \param type BAM tag type-code

    \return \c true if both \a tag and \a type are correct sizes
*/
bool BamAlignment::IsValidSize(const string& tag, const string& type) const {
    return (tag.size()  == Constants::BAM_TAG_TAGSIZE) &&
           (type.size() == Constants::BAM_TAG_TYPESIZE);
}

/*! \fn bool BamAlignment::RemoveTag(const std::string& tag)
    \brief Removes field from BAM tags.

    \return \c true if tag was removed successfully (or didn't exist before)
*/
bool BamAlignment::RemoveTag(const std::string& tag) {
  
    // skip if no tag data available
    if ( SupportData.HasCoreOnly || TagData.empty() )
        return false;
  
    // localize the tag data
    char* pOriginalTagData = (char*)TagData.data();
    char* pTagData = pOriginalTagData;
    const unsigned int originalTagDataLength = TagData.size();
    unsigned int newTagDataLength = 0;
    unsigned int numBytesParsed = 0;
    
    // if tag found
    if ( FindTag(tag, pTagData, originalTagDataLength, numBytesParsed) ) {
        
        char* newTagData = new char[originalTagDataLength];

        // copy original tag data up til desired tag
        pTagData       -= 3;
        numBytesParsed -= 3;
        const unsigned int beginningTagDataLength = numBytesParsed;
        newTagDataLength += beginningTagDataLength;
        memcpy(newTagData, pOriginalTagData, numBytesParsed);
        
        // skip to next tag (if tag for removal is last, return true) 
        const char* pTagStorageType = pTagData + 2;
        pTagData       += 3;
        numBytesParsed += 3;
        if ( !SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) )
            return true;
         
        // copy everything from current tag (the next one after tag for removal) to end
        const unsigned int skippedDataLength = (numBytesParsed - beginningTagDataLength);
        const unsigned int endTagDataLength = originalTagDataLength - beginningTagDataLength - skippedDataLength;
        memcpy(newTagData + beginningTagDataLength, pTagData, endTagDataLength );
        
        // save new tag data
        TagData.assign(newTagData, beginningTagDataLength + endTagDataLength);

        delete[] newTagData;

        return true;
    }
    
    // tag not found, no removal - return failure
    return false;
}

/*! \fn void BamAlignment::SetIsDuplicate(bool ok)
    \brief Sets value of "PCR duplicate" flag to \a ok.
*/
void BamAlignment::SetIsDuplicate(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_DUPLICATE;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_DUPLICATE;
}

/*! \fn void BamAlignment::SetIsFailedQC(bool ok)
    \brief Sets "failed quality control" flag to \a ok.
*/
void BamAlignment::SetIsFailedQC(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_QC_FAILED;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_QC_FAILED;
}

/*! \fn void BamAlignment::SetIsFirstMate(bool ok)
    \brief Sets "alignment is first mate" flag to \a ok.
*/
void BamAlignment::SetIsFirstMate(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_READ_1;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_READ_1;
}

/*! \fn void BamAlignment::SetIsMapped(bool ok)
    \brief Sets "alignment is mapped" flag to \a ok.
*/
void BamAlignment::SetIsMapped(bool ok) {
    if (ok) AlignmentFlag &= ~Constants::BAM_ALIGNMENT_UNMAPPED;
    else    AlignmentFlag |=  Constants::BAM_ALIGNMENT_UNMAPPED;
}

/*! \fn void BamAlignment::SetIsMateMapped(bool ok)
    \brief Sets "alignment's mate is mapped" flag to \a ok.
*/
void BamAlignment::SetIsMateMapped(bool ok) {
    if (ok) AlignmentFlag &= ~Constants::BAM_ALIGNMENT_MATE_UNMAPPED;
    else    AlignmentFlag |=  Constants::BAM_ALIGNMENT_MATE_UNMAPPED;
}

/*! \fn void BamAlignment::SetIsMateUnmapped(bool ok)
    \brief Complement of using SetIsMateMapped().
    \deprecated For sake of symmetry with the query methods
    \sa IsMateMapped(), SetIsMateMapped()
*/
void BamAlignment::SetIsMateUnmapped(bool ok) {
    SetIsMateMapped(!ok);
}

/*! \fn void BamAlignment::SetIsMateReverseStrand(bool ok)
    \brief Sets "alignment's mate mapped to reverse strand" flag to \a ok.
*/
void BamAlignment::SetIsMateReverseStrand(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_MATE_REVERSE_STRAND;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_MATE_REVERSE_STRAND;
}

/*! \fn void BamAlignment::SetIsPaired(bool ok)
    \brief Sets "alignment part of paired-end read" flag to \a ok.
*/
void BamAlignment::SetIsPaired(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_PAIRED;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_PAIRED;
}

/*! \fn void BamAlignment::SetIsPrimaryAlignment(bool ok)
    \brief Sets "position is primary alignment" flag to \a ok.
*/
void BamAlignment::SetIsPrimaryAlignment(bool ok) {
    if (ok) AlignmentFlag &= ~Constants::BAM_ALIGNMENT_SECONDARY;
    else    AlignmentFlag |=  Constants::BAM_ALIGNMENT_SECONDARY;
}

/*! \fn void BamAlignment::SetIsProperPair(bool ok)
    \brief Sets "alignment is part of read that satisfied paired-end resolution" flag to \a ok.
*/
void BamAlignment::SetIsProperPair(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_PROPER_PAIR;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_PROPER_PAIR;
}

/*! \fn void BamAlignment::SetIsReverseStrand(bool ok)
    \brief Sets "alignment mapped to reverse strand" flag to \a ok.
*/
void BamAlignment::SetIsReverseStrand(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_REVERSE_STRAND;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_REVERSE_STRAND;
}

/*! \fn void BamAlignment::SetIsSecondaryAlignment(bool ok)
    \brief Complement of using SetIsPrimaryAlignment().
    \deprecated For sake of symmetry with the query methods
    \sa IsPrimaryAlignment(), SetIsPrimaryAlignment()
*/
void BamAlignment::SetIsSecondaryAlignment(bool ok) {
    SetIsPrimaryAlignment(!ok);
}

/*! \fn void BamAlignment::SetIsSecondMate(bool ok)
    \brief Sets "alignment is second mate on read" flag to \a ok.
*/
void BamAlignment::SetIsSecondMate(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_READ_2;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_READ_2;
}

/*! \fn void BamAlignment::SetIsUnmapped(bool ok)
    \brief Complement of using SetIsMapped().
    \deprecated For sake of symmetry with the query methods
    \sa IsMapped(), SetIsMapped()
*/
void BamAlignment::SetIsUnmapped(bool ok) {
    SetIsMapped(!ok);
}

/*! \fn bool BamAlignment::SkipToNextTag(const char storageType, char*& pTagData, unsigned int& numBytesParsed)
    \internal

    Moves to next available tag in tag data string

    \param storageType    BAM tag type-code that determines how far to move cursor
    \param pTagData       pointer to current position (cursor) in tag string
    \param numBytesParsed report of how many bytes were parsed (cumulatively)

    \return \c if storageType was a recognized BAM tag type
    \post \a pTagData will point to the byte where the next tag data begins.
          \a numBytesParsed will correspond to the cursor's position in the full TagData string.
*/
bool BamAlignment::SkipToNextTag(const char storageType,
                                 char*& pTagData,
                                 unsigned int& numBytesParsed) const
{
    switch (storageType) {

        case (Constants::BAM_TAG_TYPE_ASCII) :
        case (Constants::BAM_TAG_TYPE_INT8)  :
        case (Constants::BAM_TAG_TYPE_UINT8) :
            ++numBytesParsed;
            ++pTagData;
            break;

        case (Constants::BAM_TAG_TYPE_INT16)  :
        case (Constants::BAM_TAG_TYPE_UINT16) :
            numBytesParsed += sizeof(uint16_t);
            pTagData       += sizeof(uint16_t);
            break;

        case (Constants::BAM_TAG_TYPE_FLOAT)  :
        case (Constants::BAM_TAG_TYPE_INT32)  :
        case (Constants::BAM_TAG_TYPE_UINT32) :
            numBytesParsed += sizeof(uint32_t);
            pTagData       += sizeof(uint32_t);
            break;

        case (Constants::BAM_TAG_TYPE_STRING) :
        case (Constants::BAM_TAG_TYPE_HEX)    :
            while( *pTagData ) {
                ++numBytesParsed;
                ++pTagData;
            }
            // increment for null-terminator
            ++numBytesParsed;
            ++pTagData;
            break;

        case (Constants::BAM_TAG_TYPE_ARRAY) :

        {
            // read array type
            const char arrayType = *pTagData;
            ++numBytesParsed;
            ++pTagData;

            // read number of elements
            int32_t numElements;
            memcpy(&numElements, pTagData, sizeof(uint32_t)); // already endian-swapped if necessary
            numBytesParsed += sizeof(uint32_t);
            pTagData       += sizeof(uint32_t);

            // calculate number of bytes to skip
            int bytesToSkip = 0;
            switch (arrayType) {
                case (Constants::BAM_TAG_TYPE_INT8)  :
                case (Constants::BAM_TAG_TYPE_UINT8) :
                    bytesToSkip = numElements;
                    break;
                case (Constants::BAM_TAG_TYPE_INT16)  :
                case (Constants::BAM_TAG_TYPE_UINT16) :
                    bytesToSkip = numElements*sizeof(uint16_t);
                    break;
                case (Constants::BAM_TAG_TYPE_FLOAT)  :
                case (Constants::BAM_TAG_TYPE_INT32)  :
                case (Constants::BAM_TAG_TYPE_UINT32) :
                    bytesToSkip = numElements*sizeof(uint32_t);
                    break;
                default:
                    cerr << "BamAlignment ERROR: unknown binary array type encountered: "
                         << arrayType << endl;
                    return false;
            }

            // skip binary array contents
            numBytesParsed += bytesToSkip;
            pTagData       += bytesToSkip;
            break;
        }

        default:
            cerr << "BamAlignment ERROR: unknown tag type encountered"
                 << storageType << endl;
            return false;
    }

    // return success
    return true;
}
