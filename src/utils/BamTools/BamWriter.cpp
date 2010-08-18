// ***************************************************************************
// BamWriter.cpp (c) 2009 Michael Str�mberg, Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 17 August 2010 (DB)
// ---------------------------------------------------------------------------
// Uses BGZF routines were adapted from the bgzf.c code developed at the Broad
// Institute.
// ---------------------------------------------------------------------------
// Provides the basic functionality for producing BAM files
// ***************************************************************************

#include <iostream>

#include "BGZF.h"
#include "BamWriter.h"
using namespace BamTools;
using namespace std;

struct BamWriter::BamWriterPrivate {

    // data members
    BgzfData mBGZF;
    bool IsBigEndian;
    
    // constructor / destructor
    BamWriterPrivate(void) { 
      IsBigEndian = SystemIsBigEndian();  
    }
    
    ~BamWriterPrivate(void) {
        mBGZF.Close();
    }

    // "public" interface
    void Close(void);
    bool Open(const string& filename, const string& samHeader, const RefVector& referenceSequences, bool isWriteUncompressed);
    void SaveAlignment(const BamAlignment& al);

    // internal methods
    const unsigned int CalculateMinimumBin(const int begin, int end) const;
    void CreatePackedCigar(const vector<CigarOp>& cigarOperations, string& packedCigar);
    void EncodeQuerySequence(const string& query, string& encodedQuery);
};

// -----------------------------------------------------
// BamWriter implementation
// -----------------------------------------------------

// constructor
BamWriter::BamWriter(void) {
    d = new BamWriterPrivate;
}

// destructor
BamWriter::~BamWriter(void) {
    delete d;
    d = 0;
}

// closes the alignment archive
void BamWriter::Close(void) { 
  d->Close(); 
}

// opens the alignment archive
bool BamWriter::Open(const string& filename, const string& samHeader, const RefVector& referenceSequences, bool isWriteUncompressed) {
    return d->Open(filename, samHeader, referenceSequences, isWriteUncompressed);
}

// saves the alignment to the alignment archive
void BamWriter::SaveAlignment(const BamAlignment& al) { 
    d->SaveAlignment(al);
}

// -----------------------------------------------------
// BamWriterPrivate implementation
// -----------------------------------------------------

// closes the alignment archive
void BamWriter::BamWriterPrivate::Close(void) {
    mBGZF.Close();
}

// calculates minimum bin for a BAM alignment interval
const unsigned int BamWriter::BamWriterPrivate::CalculateMinimumBin(const int begin, int end) const {  
    --end;
    if( (begin >> 14) == (end >> 14) ) return 4681 + (begin >> 14);
    if( (begin >> 17) == (end >> 17) ) return  585 + (begin >> 17);
    if( (begin >> 20) == (end >> 20) ) return   73 + (begin >> 20);
    if( (begin >> 23) == (end >> 23) ) return    9 + (begin >> 23);
    if( (begin >> 26) == (end >> 26) ) return    1 + (begin >> 26);
    return 0;
}

// creates a cigar string from the supplied alignment
void BamWriter::BamWriterPrivate::CreatePackedCigar(const vector<CigarOp>& cigarOperations, string& packedCigar) {

    // initialize
    const unsigned int numCigarOperations = cigarOperations.size();
    packedCigar.resize(numCigarOperations * BT_SIZEOF_INT);

    // pack the cigar data into the string
    unsigned int* pPackedCigar = (unsigned int*)packedCigar.data();

    unsigned int cigarOp;
    vector<CigarOp>::const_iterator coIter;
    for(coIter = cigarOperations.begin(); coIter != cigarOperations.end(); ++coIter) {

        switch(coIter->Type) {
            case 'M':
                  cigarOp = BAM_CMATCH;
                  break;
            case 'I':
                  cigarOp = BAM_CINS;
                  break;
            case 'D':
                  cigarOp = BAM_CDEL;
                  break;
            case 'N':
                  cigarOp = BAM_CREF_SKIP;
                  break;
            case 'S':
                  cigarOp = BAM_CSOFT_CLIP;
                  break;
            case 'H':
                  cigarOp = BAM_CHARD_CLIP;
                  break;
            case 'P':
                  cigarOp = BAM_CPAD;
                  break;
            default:
                  printf("ERROR: Unknown cigar operation found: %c\n", coIter->Type);
                  exit(1);
        }

        *pPackedCigar = coIter->Length << BAM_CIGAR_SHIFT | cigarOp;
        pPackedCigar++;
    }
}

// encodes the supplied query sequence into 4-bit notation
void BamWriter::BamWriterPrivate::EncodeQuerySequence(const string& query, string& encodedQuery) {

    // prepare the encoded query string
    const unsigned int queryLen = query.size();
    const unsigned int encodedQueryLen = (unsigned int)((queryLen / 2.0) + 0.5);
    encodedQuery.resize(encodedQueryLen);
    char* pEncodedQuery = (char*)encodedQuery.data();
    const char* pQuery = (const char*)query.data();

    unsigned char nucleotideCode;
    bool useHighWord = true;

    while(*pQuery) {

        switch(*pQuery) {
            
            case '=':
                nucleotideCode = 0;
                break;
                
            case 'A':
                nucleotideCode = 1;
                break;
            
            case 'C':
                nucleotideCode = 2;
                break;
            
            case 'G':
                nucleotideCode = 4;
                break;
            
            case 'T':
                nucleotideCode = 8;
                break;
            
            case 'N':
                nucleotideCode = 15;
                break;
            
            default:
                printf("ERROR: Only the following bases are supported in the BAM format: {=, A, C, G, T, N}. Found [%c]\n", *pQuery);
                exit(1);
        }

        // pack the nucleotide code
        if(useHighWord) {
            *pEncodedQuery = nucleotideCode << 4;
            useHighWord = false;
        } else {
            *pEncodedQuery |= nucleotideCode;
            pEncodedQuery++;
            useHighWord = true;
        }

        // increment the query position
        pQuery++;
    }
}

// opens the alignment archive
bool BamWriter::BamWriterPrivate::Open(const string& filename, const string& samHeader, const RefVector& referenceSequences, bool isWriteUncompressed) {

    // open the BGZF file for writing, return failure if error
    if ( !mBGZF.Open(filename, "wb", isWriteUncompressed) )
        return false;

    // ================
    // write the header
    // ================

    // write the BAM signature
    const unsigned char SIGNATURE_LENGTH = 4;
    const char* BAM_SIGNATURE = "BAM\1";
    mBGZF.Write(BAM_SIGNATURE, SIGNATURE_LENGTH);

    // write the SAM header text length
    uint32_t samHeaderLen = samHeader.size();
    if (IsBigEndian) SwapEndian_32(samHeaderLen);
    mBGZF.Write((char*)&samHeaderLen, BT_SIZEOF_INT);

    // write the SAM header text
    if(samHeaderLen > 0) 
        mBGZF.Write(samHeader.data(), samHeaderLen);

    // write the number of reference sequences
    uint32_t numReferenceSequences = referenceSequences.size();
    if (IsBigEndian) SwapEndian_32(numReferenceSequences);
    mBGZF.Write((char*)&numReferenceSequences, BT_SIZEOF_INT);

    // =============================
    // write the sequence dictionary
    // =============================

    RefVector::const_iterator rsIter;
    for(rsIter = referenceSequences.begin(); rsIter != referenceSequences.end(); rsIter++) {

        // write the reference sequence name length
        uint32_t referenceSequenceNameLen = rsIter->RefName.size() + 1;
        if (IsBigEndian) SwapEndian_32(referenceSequenceNameLen);
        mBGZF.Write((char*)&referenceSequenceNameLen, BT_SIZEOF_INT);

        // write the reference sequence name
        mBGZF.Write(rsIter->RefName.c_str(), referenceSequenceNameLen);

        // write the reference sequence length
        int32_t referenceLength = rsIter->RefLength;
        if (IsBigEndian) SwapEndian_32(referenceLength);
        mBGZF.Write((char*)&referenceLength, BT_SIZEOF_INT);
    }
    
    // return success
    return true;
}

// saves the alignment to the alignment archive
void BamWriter::BamWriterPrivate::SaveAlignment(const BamAlignment& al) {

    // if BamAlignment contains only the core data and a raw char data buffer
    // (as a result of BamReader::GetNextAlignmentCore())
    if ( al.SupportData.HasCoreOnly ) {
      
        // write the block size
        unsigned int blockSize = al.SupportData.BlockLength;
        if (IsBigEndian) SwapEndian_32(blockSize);
        mBGZF.Write((char*)&blockSize, BT_SIZEOF_INT);

        // assign the BAM core data
        uint32_t buffer[8];
        buffer[0] = al.RefID;
        buffer[1] = al.Position;
        buffer[2] = (al.Bin << 16) | (al.MapQuality << 8) | al.SupportData.QueryNameLength;
        buffer[3] = (al.AlignmentFlag << 16) | al.SupportData.NumCigarOperations;
        buffer[4] = al.SupportData.QuerySequenceLength;
        buffer[5] = al.MateRefID;
        buffer[6] = al.MatePosition;
        buffer[7] = al.InsertSize;
        
        // swap BAM core endian-ness, if necessary
        if ( IsBigEndian ) { 
            for ( int i = 0; i < 8; ++i )
                SwapEndian_32(buffer[i]); 
        }
        
        // write the BAM core
        mBGZF.Write((char*)&buffer, BAM_CORE_SIZE);
        
        // write the raw char data
        mBGZF.Write((char*)al.SupportData.AllCharData.data(), al.SupportData.BlockLength-BAM_CORE_SIZE); 
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
        const unsigned int dataBlockSize = nameLength + packedCigarLength + encodedQueryLength + queryLength + tagDataLength;
        unsigned int blockSize = BAM_CORE_SIZE + dataBlockSize;
        if (IsBigEndian) SwapEndian_32(blockSize);
        mBGZF.Write((char*)&blockSize, BT_SIZEOF_INT);

        // assign the BAM core data
        uint32_t buffer[8];
        buffer[0] = al.RefID;
        buffer[1] = al.Position;
        buffer[2] = (alignmentBin << 16) | (al.MapQuality << 8) | nameLength;
        buffer[3] = (al.AlignmentFlag << 16) | numCigarOperations;
        buffer[4] = queryLength;
        buffer[5] = al.MateRefID;
        buffer[6] = al.MatePosition;
        buffer[7] = al.InsertSize;
        
        // swap BAM core endian-ness, if necessary
        if ( IsBigEndian ) { 
            for ( int i = 0; i < 8; ++i )
                SwapEndian_32(buffer[i]); 
        }
        
        // write the BAM core
        mBGZF.Write((char*)&buffer, BAM_CORE_SIZE);
        
        // write the query name
        mBGZF.Write(al.Name.c_str(), nameLength);

        // write the packed cigar
        if ( IsBigEndian ) {
          
            char* cigarData = (char*)calloc(sizeof(char), packedCigarLength);
            memcpy(cigarData, packedCigar.data(), packedCigarLength);
            
            for (unsigned int i = 0; i < packedCigarLength; ++i) {
                if ( IsBigEndian )
                  SwapEndian_32p(&cigarData[i]); 
            }
            
            mBGZF.Write(cigarData, packedCigarLength);
            free(cigarData);    
        } 
        else 
            mBGZF.Write(packedCigar.data(), packedCigarLength);

        // write the encoded query sequence
        mBGZF.Write(encodedQuery.data(), encodedQueryLength);

        // write the base qualities
        string baseQualities(al.Qualities);
        char* pBaseQualities = (char*)al.Qualities.data();
        for(unsigned int i = 0; i < queryLength; i++) { 
            pBaseQualities[i] -= 33; 
        }
        mBGZF.Write(pBaseQualities, queryLength);

        // write the read group tag
        if ( IsBigEndian ) {
          
            char* tagData = (char*)calloc(sizeof(char), tagDataLength);
            memcpy(tagData, al.TagData.data(), tagDataLength);
          
            int i = 0;
            while ( (unsigned int)i < tagDataLength ) {
                
                i += 2;                                 // skip tag type (e.g. "RG", "NM", etc)
                uint8_t type = toupper(tagData[i]);     // lower & upper case letters have same meaning 
                ++i;                                    // skip value type
                
                switch (type) {
                  
                    case('A') :
                    case('C') : 
                        ++i;
                        break;
                                
                    case('S') : 
                        SwapEndian_16p(&tagData[i]); 
                        i+=2; // sizeof(uint16_t)
                        break;
                                
                    case('F') :
                    case('I') : 
                        SwapEndian_32p(&tagData[i]);
                        i+=4; // sizeof(uint32_t)
                        break;
                                
                    case('D') : 
                        SwapEndian_64p(&tagData[i]);
                        i+=8; // sizeof(uint64_t)
                        break;
                                
                    case('H') :
                    case('Z') : 
                        while (tagData[i]) { ++i; }
                        ++i; // increment one more for null terminator
                        break;
                                
                    default : 
                        printf("ERROR: Invalid tag value type\n"); // shouldn't get here
                        free(tagData);
                        exit(1); 
                }
            }
            
            mBGZF.Write(tagData, tagDataLength);
            free(tagData);
        } 
        else 
            mBGZF.Write(al.TagData.data(), tagDataLength);      
    }
}
