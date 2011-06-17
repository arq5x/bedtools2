// ***************************************************************************
// BamWriter_p.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 16 June 2011 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for producing BAM files
// ***************************************************************************

#include <api/BamAlignment.h>
#include <api/BamConstants.h>
#include <api/internal/BamWriter_p.h>
using namespace BamTools;
using namespace BamTools::Internal;

#include <cstdio>
#include <cstdlib>
#include <cstring>
using namespace std;

// ctor
BamWriterPrivate::BamWriterPrivate(void)
    : m_isBigEndian( BamTools::SystemIsBigEndian() )
{ }

// dtor
BamWriterPrivate::~BamWriterPrivate(void) {
    m_stream.Close();
}

// calculates minimum bin for a BAM alignment interval
unsigned int BamWriterPrivate::CalculateMinimumBin(const int begin, int end) const {
    --end;
    if ( (begin >> 14) == (end >> 14) ) return 4681 + (begin >> 14);
    if ( (begin >> 17) == (end >> 17) ) return  585 + (begin >> 17);
    if ( (begin >> 20) == (end >> 20) ) return   73 + (begin >> 20);
    if ( (begin >> 23) == (end >> 23) ) return    9 + (begin >> 23);
    if ( (begin >> 26) == (end >> 26) ) return    1 + (begin >> 26);
    return 0;
}

// closes the alignment archive
void BamWriterPrivate::Close(void) {
    m_stream.Close();
}

// creates a cigar string from the supplied alignment
void BamWriterPrivate::CreatePackedCigar(const vector<CigarOp>& cigarOperations, string& packedCigar) {

    // initialize
    const unsigned int numCigarOperations = cigarOperations.size();
    packedCigar.resize(numCigarOperations * Constants::BAM_SIZEOF_INT);

    // pack the cigar data into the string
    unsigned int* pPackedCigar = (unsigned int*)packedCigar.data();

    // iterate over cigar operations
    vector<CigarOp>::const_iterator coIter = cigarOperations.begin();
    vector<CigarOp>::const_iterator coEnd  = cigarOperations.end();
    for ( ; coIter != coEnd; ++coIter ) {

        // store op in packedCigar
        unsigned int cigarOp;
        switch ( coIter->Type ) {
            case (Constants::BAM_CIGAR_MATCH_CHAR)    : cigarOp = Constants::BAM_CIGAR_MATCH;    break;
            case (Constants::BAM_CIGAR_INS_CHAR)      : cigarOp = Constants::BAM_CIGAR_INS;      break;
            case (Constants::BAM_CIGAR_DEL_CHAR)      : cigarOp = Constants::BAM_CIGAR_DEL;      break;
            case (Constants::BAM_CIGAR_REFSKIP_CHAR)  : cigarOp = Constants::BAM_CIGAR_REFSKIP;  break;
            case (Constants::BAM_CIGAR_SOFTCLIP_CHAR) : cigarOp = Constants::BAM_CIGAR_SOFTCLIP; break;
            case (Constants::BAM_CIGAR_HARDCLIP_CHAR) : cigarOp = Constants::BAM_CIGAR_HARDCLIP; break;
            case (Constants::BAM_CIGAR_PAD_CHAR)      : cigarOp = Constants::BAM_CIGAR_PAD;      break;
            case (Constants::BAM_CIGAR_SEQMATCH_CHAR) : cigarOp = Constants::BAM_CIGAR_SEQMATCH; break;
            case (Constants::BAM_CIGAR_MISMATCH_CHAR) : cigarOp = Constants::BAM_CIGAR_MISMATCH; break;
            default:
              fprintf(stderr, "BamWriter ERROR: unknown cigar operation found: %c\n", coIter->Type);
              exit(1);
        }

        *pPackedCigar = coIter->Length << Constants::BAM_CIGAR_SHIFT | cigarOp;
        pPackedCigar++;
    }
}

// encodes the supplied query sequence into 4-bit notation
void BamWriterPrivate::EncodeQuerySequence(const string& query, string& encodedQuery) {

    // prepare the encoded query string
    const unsigned int queryLen = query.size();
    const unsigned int encodedQueryLen = (unsigned int)((queryLen / 2.0) + 0.5);
    encodedQuery.resize(encodedQueryLen);
    char* pEncodedQuery = (char*)encodedQuery.data();
    const char* pQuery = (const char*)query.data();

    unsigned char nucleotideCode;
    bool useHighWord = true;

    while ( *pQuery ) {
        switch ( *pQuery ) {
            case (Constants::BAM_DNA_EQUAL) : nucleotideCode = Constants::BAM_BASECODE_EQUAL; break;
            case (Constants::BAM_DNA_A)     : nucleotideCode = Constants::BAM_BASECODE_A;     break;
            case (Constants::BAM_DNA_C)     : nucleotideCode = Constants::BAM_BASECODE_C;     break;
            case (Constants::BAM_DNA_G)     : nucleotideCode = Constants::BAM_BASECODE_G;     break;
            case (Constants::BAM_DNA_T)     : nucleotideCode = Constants::BAM_BASECODE_T;     break;
            case (Constants::BAM_DNA_N)     : nucleotideCode = Constants::BAM_BASECODE_N;     break;
            default:
                fprintf(stderr, "BamWriter ERROR: only the following bases are supported in the BAM format: {=, A, C, G, T, N}. Found [%c]\n", *pQuery);
                exit(1);
        }

        // pack the nucleotide code
        if ( useHighWord ) {
            *pEncodedQuery = nucleotideCode << 4;
            useHighWord = false;
        } else {
            *pEncodedQuery |= nucleotideCode;
            ++pEncodedQuery;
            useHighWord = true;
        }

        // increment the query position
        ++pQuery;
    }
}

// returns whether BAM file is open for writing or not
bool BamWriterPrivate::IsOpen(void) const {
    return m_stream.IsOpen;
}

// opens the alignment archive
bool BamWriterPrivate::Open(const string& filename,
                            const string& samHeaderText,
                            const RefVector& referenceSequences)
{
    // open the BGZF file for writing, return failure if error
    if ( !m_stream.Open(filename, "wb") )
        return false;

    // write BAM file 'metadata' components
    WriteMagicNumber();
    WriteSamHeaderText(samHeaderText);
    WriteReferences(referenceSequences);
    return true;
}

// saves the alignment to the alignment archive
void BamWriterPrivate::SaveAlignment(const BamAlignment& al) {

    // if BamAlignment contains only the core data and a raw char data buffer
    // (as a result of BamReader::GetNextAlignmentCore())
    if ( al.SupportData.HasCoreOnly ) {

        // write the block size
        unsigned int blockSize = al.SupportData.BlockLength;
        if ( m_isBigEndian ) BamTools::SwapEndian_32(blockSize);
        m_stream.Write((char*)&blockSize, Constants::BAM_SIZEOF_INT);

        // re-calculate bin (in case BamAlignment's position has been previously modified)
        const uint32_t alignmentBin = CalculateMinimumBin(al.Position, al.GetEndPosition());

        // assign the BAM core data
        uint32_t buffer[Constants::BAM_CORE_BUFFER_SIZE];
        buffer[0] = al.RefID;
        buffer[1] = al.Position;
        buffer[2] = (alignmentBin << 16) | (al.MapQuality << 8) | al.SupportData.QueryNameLength;
        buffer[3] = (al.AlignmentFlag << 16) | al.SupportData.NumCigarOperations;
        buffer[4] = al.SupportData.QuerySequenceLength;
        buffer[5] = al.MateRefID;
        buffer[6] = al.MatePosition;
        buffer[7] = al.InsertSize;

        // swap BAM core endian-ness, if necessary
        if ( m_isBigEndian ) {
            for ( int i = 0; i < 8; ++i )
                BamTools::SwapEndian_32(buffer[i]);
        }

        // write the BAM core
        m_stream.Write((char*)&buffer, Constants::BAM_CORE_SIZE);

        // write the raw char data
        m_stream.Write((char*)al.SupportData.AllCharData.data(),
                       al.SupportData.BlockLength-Constants::BAM_CORE_SIZE);
    }

    // otherwise, BamAlignment should contain character in the standard fields: Name, QueryBases, etc
    // ( resulting from BamReader::GetNextAlignment() *OR* being generated directly by client code )
    else {

        // calculate char lengths
        const unsigned int nameLength         = al.Name.size() + 1;
        const unsigned int numCigarOperations = al.CigarData.size();
        const unsigned int queryLength        = al.QueryBases.size();
        const unsigned int tagDataLength      = al.TagData.size();

        // no way to tell if BamAlignment.Bin is already defined (no default, invalid value)
        // force calculation of Bin before storing
        const int endPosition = al.GetEndPosition();
        const unsigned int alignmentBin = CalculateMinimumBin(al.Position, endPosition);

        // create our packed cigar string
        string packedCigar;
        CreatePackedCigar(al.CigarData, packedCigar);
        const unsigned int packedCigarLength = packedCigar.size();

        // encode the query
        string encodedQuery;
        EncodeQuerySequence(al.QueryBases, encodedQuery);
        const unsigned int encodedQueryLength = encodedQuery.size();

        // write the block size
        const unsigned int dataBlockSize = nameLength +
                                           packedCigarLength +
                                           encodedQueryLength +
                                           queryLength +
                                           tagDataLength;
        unsigned int blockSize = Constants::BAM_CORE_SIZE + dataBlockSize;
        if ( m_isBigEndian ) BamTools::SwapEndian_32(blockSize);
        m_stream.Write((char*)&blockSize, Constants::BAM_SIZEOF_INT);

        // assign the BAM core data
        uint32_t buffer[Constants::BAM_CORE_BUFFER_SIZE];
        buffer[0] = al.RefID;
        buffer[1] = al.Position;
        buffer[2] = (alignmentBin << 16) | (al.MapQuality << 8) | nameLength;
        buffer[3] = (al.AlignmentFlag << 16) | numCigarOperations;
        buffer[4] = queryLength;
        buffer[5] = al.MateRefID;
        buffer[6] = al.MatePosition;
        buffer[7] = al.InsertSize;

        // swap BAM core endian-ness, if necessary
        if ( m_isBigEndian ) {
            for ( int i = 0; i < 8; ++i )
                BamTools::SwapEndian_32(buffer[i]);
        }

        // write the BAM core
        m_stream.Write((char*)&buffer, Constants::BAM_CORE_SIZE);

        // write the query name
        m_stream.Write(al.Name.c_str(), nameLength);

        // write the packed cigar
        if ( m_isBigEndian ) {
            char* cigarData = (char*)calloc(sizeof(char), packedCigarLength);
            memcpy(cigarData, packedCigar.data(), packedCigarLength);
            if ( m_isBigEndian ) {
                for ( unsigned int i = 0; i < packedCigarLength; ++i )
                    BamTools::SwapEndian_32p(&cigarData[i]);
            }
            m_stream.Write(cigarData, packedCigarLength);
            free(cigarData);
        }
        else
            m_stream.Write(packedCigar.data(), packedCigarLength);

        // write the encoded query sequence
        m_stream.Write(encodedQuery.data(), encodedQueryLength);

        // write the base qualities
        char* pBaseQualities = (char*)al.Qualities.data();
        for ( unsigned int i = 0; i < queryLength; ++i )
            pBaseQualities[i] -= 33; // FASTQ conversion
        m_stream.Write(pBaseQualities, queryLength);

        // write the read group tag
        if ( m_isBigEndian ) {

            char* tagData = (char*)calloc(sizeof(char), tagDataLength);
            memcpy(tagData, al.TagData.data(), tagDataLength);

            int i = 0;
            while ( (unsigned int)i < tagDataLength ) {

                i += Constants::BAM_TAG_TAGSIZE;  // skip tag chars (e.g. "RG", "NM", etc.)
                const char type = tagData[i];     // get tag type at position i
                ++i;

                switch ( type ) {

                    case(Constants::BAM_TAG_TYPE_ASCII) :
                    case(Constants::BAM_TAG_TYPE_INT8)  :
                    case(Constants::BAM_TAG_TYPE_UINT8) :
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
                                    fprintf(stderr,
                                            "BamWriter ERROR: unknown binary array type encountered: [%c]\n",
                                            arrayType);
                                    exit(1);
                            }
                        }

                        break;
                    }

                    default :
                        fprintf(stderr, "BamWriter ERROR: invalid tag value type\n"); // shouldn't get here
                        free(tagData);
                        exit(1);
                }
            }
            m_stream.Write(tagData, tagDataLength);
            free(tagData);
        }
        else
            m_stream.Write(al.TagData.data(), tagDataLength);
    }
}

void BamWriterPrivate::SetWriteCompressed(bool ok) {

    // warn if BAM file is already open
    // modifying compression is not allowed in this case
    if ( IsOpen() ) {
        cerr << "BamWriter WARNING: attempting to change compression mode on an open BAM file is not allowed. "
             << "Ignoring request." << endl;
        return;
    }

    // set BgzfStream compression mode
    m_stream.SetWriteCompressed(ok);
}

void BamWriterPrivate::WriteMagicNumber(void) {
    // write BAM file 'magic number'
    m_stream.Write(Constants::BAM_HEADER_MAGIC, Constants::BAM_HEADER_MAGIC_LENGTH);
}

void BamWriterPrivate::WriteReferences(const BamTools::RefVector& referenceSequences) {

    // write the number of reference sequences
    uint32_t numReferenceSequences = referenceSequences.size();
    if ( m_isBigEndian ) BamTools::SwapEndian_32(numReferenceSequences);
    m_stream.Write((char*)&numReferenceSequences, Constants::BAM_SIZEOF_INT);

    // foreach reference sequence
    RefVector::const_iterator rsIter = referenceSequences.begin();
    RefVector::const_iterator rsEnd  = referenceSequences.end();
    for ( ; rsIter != rsEnd; ++rsIter ) {

        // write the reference sequence name length
        uint32_t referenceSequenceNameLen = rsIter->RefName.size() + 1;
        if ( m_isBigEndian ) BamTools::SwapEndian_32(referenceSequenceNameLen);
        m_stream.Write((char*)&referenceSequenceNameLen, Constants::BAM_SIZEOF_INT);

        // write the reference sequence name
        m_stream.Write(rsIter->RefName.c_str(), referenceSequenceNameLen);

        // write the reference sequence length
        int32_t referenceLength = rsIter->RefLength;
        if ( m_isBigEndian ) BamTools::SwapEndian_32(referenceLength);
        m_stream.Write((char*)&referenceLength, Constants::BAM_SIZEOF_INT);
    }
}

void BamWriterPrivate::WriteSamHeaderText(const std::string& samHeaderText) {

    // write the SAM header  text length
    uint32_t samHeaderLen = samHeaderText.size();
    if ( m_isBigEndian ) BamTools::SwapEndian_32(samHeaderLen);
    m_stream.Write((char*)&samHeaderLen, Constants::BAM_SIZEOF_INT);

    // write the SAM header text
    if ( samHeaderLen > 0 )
        m_stream.Write(samHeaderText.data(), samHeaderLen);
}
