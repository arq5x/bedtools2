// ***************************************************************************
// BamAlignment.cpp (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 17 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides the BamAlignment data structure
// ***************************************************************************

#include "api/BamAlignment.h"
#include "api/BamConstants.h"
using namespace BamTools;
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
    const unsigned int seqDataOffset  = SupportData.QueryNameLength + (SupportData.NumCigarOperations*4);
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
        for ( size_t i = 0; i < SupportData.QuerySequenceLength; ++i ) {
            const char singleBase = Constants::BAM_DNA_LOOKUP[ ( (seqData[(i/2)] >> (4*(1-(i%2)))) & 0xf ) ];
            QueryBases.append(1, singleBase);
        }
    }

    // save qualities, converting from numeric QV to 'FASTQ-style' ASCII character
    Qualities.clear();
    if ( hasQualData ) {
        Qualities.reserve(SupportData.QuerySequenceLength);
        for ( size_t i = 0; i < SupportData.QuerySequenceLength; ++i ) {
            const char singleQuality = static_cast<const char>(qualData[i]+33);
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

            switch ( op.Type ) {

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

                // invalid CIGAR op-code
                default:
                    const string message = string("invalid CIGAR operation type: ") + op.Type;
                    SetErrorString("BamAlignment::BuildCharData", message);
                    return false;
            }
        }
    }

    // save tag data
    TagData.clear();
    if ( hasTagData ) {
        if ( IsBigEndian ) {
            size_t i = 0;
            while ( i < tagDataLength ) {

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
                        uint32_t numElements;
                        memcpy(&numElements, &tagData[i], sizeof(uint32_t));
                        i += sizeof(uint32_t);

                        // swap endian-ness of array elements
                        for ( size_t j = 0; j < numElements; ++j ) {
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
                                    const string message = string("invalid binary array type: ") + arrayType;
                                    SetErrorString("BamAlignment::BuildCharData", message);
                                    return false;
                            }
                        }

                        break;
                    }

                    // invalid tag type-code
                    default :
                        const string message = string("invalid tag type: ") + type;
                        SetErrorString("BamAlignment::BuildCharData", message);
                        return false;
                }
            }
        }

        // store tagData in alignment
        TagData.resize(tagDataLength);
        memcpy((char*)(TagData.data()), tagData, tagDataLength);
    }

    // clear core-only flag & return success
    SupportData.HasCoreOnly = false;
    return true;
}

/*! \fn bool BamAlignment::FindTag(const std::string& tag, char*& pTagData, const unsigned int& tagDataLength, unsigned int& numBytesParsed) const
    \internal

    Searches for requested tag in BAM tag data.

    \param[in]     tag            requested 2-character tag name
    \param[in,out] pTagData       pointer to current position in BamAlignment::TagData
    \param[in]     tagDataLength  length of BamAlignment::TagData
    \param[in,out] numBytesParsed number of bytes parsed so far

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

/*! \fn int BamAlignment::GetEndPosition(bool usePadded = false, bool closedInterval = false) const
    \brief Calculates alignment end position, based on its starting position and CIGAR data.

    \warning The position returned now represents a zero-based, HALF-OPEN interval.
    In previous versions of BamTools (0.x & 1.x) all intervals were treated
    as zero-based, CLOSED.

    \param[in] usePadded      Allow inserted bases to affect the reported position. Default is
                              false, so that reported position stays synced with reference
                              coordinates.
    \param[in] closedInterval Setting this to true will return a 0-based end coordinate. Default is
                              false, so that his value represents a standard, half-open interval.

    \return alignment end position
*/
int BamAlignment::GetEndPosition(bool usePadded, bool closedInterval) const {

    // initialize alignment end to starting position
    int alignEnd = Position;

    // iterate over cigar operations
    vector<CigarOp>::const_iterator cigarIter = CigarData.begin();
    vector<CigarOp>::const_iterator cigarEnd  = CigarData.end();
    for ( ; cigarIter != cigarEnd; ++cigarIter) {
        const CigarOp& op = (*cigarIter);

        switch ( op.Type ) {

            // increase end position on CIGAR chars [DMXN=]
            case Constants::BAM_CIGAR_DEL_CHAR      :
            case Constants::BAM_CIGAR_MATCH_CHAR    :
            case Constants::BAM_CIGAR_MISMATCH_CHAR :
            case Constants::BAM_CIGAR_REFSKIP_CHAR  :
            case Constants::BAM_CIGAR_SEQMATCH_CHAR :
                alignEnd += op.Length;
                break;

            // increase end position on insertion, only if @usePadded is true
            case Constants::BAM_CIGAR_INS_CHAR :
                if ( usePadded )
                    alignEnd += op.Length;
                break;

            // all other CIGAR chars do not affect end position
            default :
                break;
        }
    }

    // adjust for closedInterval, if requested
    if ( closedInterval )
        alignEnd -= 1;

    // return result
    return alignEnd;
}

/*! \fn std::string BamAlignment::GetErrorString(void) const
    \brief Returns a human-readable description of the last error that occurred

    This method allows elimination of STDERR pollution. Developers of client code
    may choose how the messages are displayed to the user, if at all.

    \return error description
*/
std::string BamAlignment::GetErrorString(void) const {
    return ErrorString;
}

/*! \fn bool BamAlignment::GetSoftClips(std::vector<int>& clipSizes, std::vector<int>& readPositions, std::vector<int>& genomePositions, bool usePadded = false) const
    \brief Identifies if an alignment has a soft clip. If so, identifies the
           sizes of the soft clips, as well as their positions in the read and reference.

    \param[out] clipSizes       vector of the sizes of each soft clip in the alignment
    \param[out] readPositions   vector of the 0-based read locations of each soft clip in the alignment.
                                These positions are basically indexes within the read, not genomic positions.
    \param[out] genomePositions vector of the 0-based genome locations of each soft clip in the alignment
    \param[in]  usePadded       inserted bases affect reported position. Default is false, so that
                                reported position stays 'sync-ed' with reference coordinates.

    \return \c true if any soft clips were found in the alignment
*/
bool BamAlignment::GetSoftClips(vector<int>& clipSizes,
                                vector<int>& readPositions,
                                vector<int>& genomePositions,
                                bool usePadded) const
{
    // initialize positions & flags
    int refPosition  = Position;
    int readPosition = 0;
    bool softClipFound = false;
    bool firstCigarOp  = true;

    // iterate over cigar operations
    vector<CigarOp>::const_iterator cigarIter = CigarData.begin();
    vector<CigarOp>::const_iterator cigarEnd  = CigarData.end();
    for ( ; cigarIter != cigarEnd; ++cigarIter) {
        const CigarOp& op = (*cigarIter);

        switch ( op.Type ) {

            // increase both read & genome positions on CIGAR chars [DMXN=]
            case Constants::BAM_CIGAR_DEL_CHAR      :
            case Constants::BAM_CIGAR_MATCH_CHAR    :
            case Constants::BAM_CIGAR_MISMATCH_CHAR :
            case Constants::BAM_CIGAR_REFSKIP_CHAR  :
            case Constants::BAM_CIGAR_SEQMATCH_CHAR :
                refPosition  += op.Length;
                readPosition += op.Length;
                break;

            // increase read position on insertion, genome position only if @usePadded is true
            case Constants::BAM_CIGAR_INS_CHAR :
                readPosition += op.Length;
                if ( usePadded )
                    refPosition += op.Length;
                break;

            case Constants::BAM_CIGAR_SOFTCLIP_CHAR :

                softClipFound = true;

                //////////////////////////////////////////////////////////////////////////////
                // if we are dealing with the *first* CIGAR operation
                // for this alignment, we increment the read position so that
                // the read and genome position of the clip are referring to the same base.
                // For example, in the alignment below, the ref position would be 4, yet
                //              the read position would be 0. Thus, to "sync" the two,
                //              we need to increment the read position by the length of the
                //              soft clip.
                // Read:  ATCGTTTCGTCCCTGC
                // Ref:   GGGATTTCGTCCCTGC
                // Cigar: SSSSMMMMMMMMMMMM
                //
                // NOTE: This only needs to be done if the soft clip is the _first_ CIGAR op.
                //////////////////////////////////////////////////////////////////////////////
                if ( firstCigarOp )
                    readPosition += op.Length;

                // track the soft clip's size, read position, and genome position
                clipSizes.push_back(op.Length);
                readPositions.push_back(readPosition);
                genomePositions.push_back(refPosition);

            // any other CIGAR operations have no effect
            default :
                break;
        }

        // clear our "first pass" flag
        firstCigarOp = false;
    }

    // return whether any soft clips found
    return softClipFound;
}

/*! \fn bool BamAlignment::GetTagType(const std::string& tag, char& type) const
    \brief Retrieves the BAM tag type-code associated with requested tag name.

    \param[in]  tag  2-character tag name
    \param[out] type retrieved (1-character) type-code

    \return \c true if found
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::GetTagType(const std::string& tag, char& type) const {
  
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
    
    // if tag not found, return failure
    if ( !FindTag(tag, pTagData, tagDataLength, numBytesParsed) ){
        // TODO: set error string?
        return false;
    }

    // otherwise, retrieve & validate tag type code
    type = *(pTagData - 1);
    switch ( type ) {
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
            const string message = string("invalid tag type: ") + type;
            SetErrorString("BamAlignment::GetTagType", message);
            return false;
    }
}

/*! \fn bool BamAlignment::HasTag(const std::string& tag) const
    \brief Returns true if alignment has a record for requested tag.

    \param[in] tag 2-character tag name
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

/*! \fn bool BamAlignment::IsValidSize(const std::string& tag, const std::string& type) const
    \internal

    Checks that tag name & type strings are expected sizes.

    \param tag[in]  BAM tag name
    \param type[in] BAM tag type-code
    \return \c true if both input strings are valid sizes
*/
bool BamAlignment::IsValidSize(const std::string& tag, const std::string& type) const {
    return (tag.size()  == Constants::BAM_TAG_TAGSIZE) &&
           (type.size() == Constants::BAM_TAG_TYPESIZE);
}

/*! \fn void BamAlignment::RemoveTag(const std::string& tag)
    \brief Removes field from BAM tags.

    \param[in] tag 2-character name of field to remove
*/
void BamAlignment::RemoveTag(const std::string& tag) {
  
    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();

    // skip if no tags available
    if ( TagData.empty() )
        return;
  
    // localize the tag data
    char* pOriginalTagData = (char*)TagData.data();
    char* pTagData = pOriginalTagData;
    const unsigned int originalTagDataLength = TagData.size();
    unsigned int newTagDataLength = 0;
    unsigned int numBytesParsed = 0;

    // skip if tag not found
    if  ( !FindTag(tag, pTagData, originalTagDataLength, numBytesParsed) )
        return;

    // otherwise, remove it
    RaiiBuffer newTagData(originalTagDataLength);

    // copy original tag data up til desired tag
    pTagData       -= 3;
    numBytesParsed -= 3;
    const unsigned int beginningTagDataLength = numBytesParsed;
    newTagDataLength += beginningTagDataLength;
    memcpy(newTagData.Buffer, pOriginalTagData, numBytesParsed);

    // attemp to skip to next tag
    const char* pTagStorageType = pTagData + 2;
    pTagData       += 3;
    numBytesParsed += 3;
    if ( SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) ) {

        // squeeze remaining tag data
        const unsigned int skippedDataLength = (numBytesParsed - beginningTagDataLength);
        const unsigned int endTagDataLength = originalTagDataLength - beginningTagDataLength - skippedDataLength;
        memcpy(newTagData.Buffer + beginningTagDataLength, pTagData, endTagDataLength );

        // save modified tag data in alignment
        TagData.assign(newTagData.Buffer, beginningTagDataLength + endTagDataLength);
    }
}

/*! \fn void BamAlignment::SetErrorString(const std::string& where, const std::string& what) const
    \internal

    Sets a formatted error string for this alignment.

    \param[in] where class/method where error occurred
    \param[in] what  description of error
*/
void BamAlignment::SetErrorString(const std::string& where, const std::string& what) const {
    static const string SEPARATOR = ": ";
    ErrorString = where + SEPARATOR + what;
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

/*! \fn void BamAlignment::SetIsSecondMate(bool ok)
    \brief Sets "alignment is second mate on read" flag to \a ok.
*/
void BamAlignment::SetIsSecondMate(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_READ_2;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_READ_2;
}

/*! \fn bool BamAlignment::SkipToNextTag(const char storageType, char*& pTagData, unsigned int& numBytesParsed) const
    \internal

    Moves to next available tag in tag data string

    \param[in]     storageType    BAM tag type-code that determines how far to move cursor
    \param[in,out] pTagData       pointer to current position (cursor) in tag string
    \param[in,out] numBytesParsed report of how many bytes were parsed (cumulatively)

    \return \c if storageType was a recognized BAM tag type

    \post \a pTagData       will point to the byte where the next tag data begins.
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
            memcpy(&numElements, pTagData, sizeof(uint32_t)); // already endian-swapped, if needed
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
                    const string message = string("invalid binary array type: ") + arrayType;
                    SetErrorString("BamAlignment::SkipToNextTag", message);
                    return false;
            }

            // skip binary array contents
            numBytesParsed += bytesToSkip;
            pTagData       += bytesToSkip;
            break;
        }

        default:
            const string message = string("invalid tag type: ") + storageType;
            SetErrorString("BamAlignment::SkipToNextTag", message);
            return false;
    }

    // if we get here, tag skipped OK - return success
    return true;
}
