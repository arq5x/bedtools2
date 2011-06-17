// ***************************************************************************
// BamHeader_p.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 21 March 2011 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for handling BAM headers.
// ***************************************************************************

#include <api/BamAux.h>
#include <api/BamConstants.h>
#include <api/internal/BamHeader_p.h>
#include <api/internal/BgzfStream_p.h>
using namespace BamTools;
using namespace BamTools::Internal;

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
using namespace std;

// ctor
BamHeader::BamHeader(void) { }

// dtor
BamHeader::~BamHeader(void) { }

// reads magic number from BGZF stream, returns true if valid
bool BamHeader::CheckMagicNumber(BgzfStream* stream) {

    // try to read magic number
    char buffer[Constants::BAM_HEADER_MAGIC_LENGTH];
    if ( stream->Read(buffer, Constants::BAM_HEADER_MAGIC_LENGTH) != (int)Constants::BAM_HEADER_MAGIC_LENGTH ) {
        fprintf(stderr, "BamHeader ERROR: could not read magic number\n");
        return false;
    }

    // validate magic number
    if ( strncmp(buffer, Constants::BAM_HEADER_MAGIC, Constants::BAM_HEADER_MAGIC_LENGTH) != 0 ) {
        fprintf(stderr, "BamHeader ERROR: invalid magic number\n");
        return false;
    }

    // all checks out
    return true;
}

// clear SamHeader data
void BamHeader::Clear(void) {
    m_header.Clear();
}

// return true if SamHeader data is valid
bool BamHeader::IsValid(void) const {
    return m_header.IsValid();
}

// load BAM header ('magic number' and SAM header text) from BGZF stream
// returns true if all OK
bool BamHeader::Load(BgzfStream* stream) {

    // cannot load if invalid stream
    if ( stream == 0 )
        return false;

    // cannot load if magic number is invalid
    if ( !CheckMagicNumber(stream) )
        return false;

    // cannot load header if cannot read header length
    uint32_t length(0);
    if ( !ReadHeaderLength(stream, length) )
        return false;

    // cannot load header if cannot read header text
    if ( !ReadHeaderText(stream, length) )
        return false;

    // otherwise, everything OK
    return true;
}

// reads SAM header text length from BGZF stream, stores it in @length
// returns read success/fail status
bool BamHeader::ReadHeaderLength(BgzfStream* stream, uint32_t& length) {

    // attempt to read BAM header text length
    char buffer[sizeof(uint32_t)];
    if ( stream->Read(buffer, sizeof(uint32_t)) != sizeof(uint32_t) ) {
        fprintf(stderr, "BamHeader ERROR: could not read header length\n");
        return false;
    }

    // convert char buffer to length, return success
    length = BamTools::UnpackUnsignedInt(buffer);
    if ( BamTools::SystemIsBigEndian() )
        BamTools::SwapEndian_32(length);
    return true;
}

// reads SAM header text from BGZF stream, stores in SamHeader object
// returns read success/fail status
bool BamHeader::ReadHeaderText(BgzfStream* stream, const uint32_t& length) {

    // set up destination buffer
    char* headerText = (char*)calloc(length + 1, 1);

    // attempt to read header text
    const unsigned bytesRead = stream->Read(headerText, length);
    const bool readOk = ( bytesRead == length );
    if ( readOk )
        m_header.SetHeaderText( (string)((const char*)headerText) );
    else
        fprintf(stderr, "BamHeader ERROR: could not read header text\n");

    // clean up calloc-ed temp variable (on success or fail)
    free(headerText);

    // return read success
    return readOk;
}

// returns *copy* of SamHeader data object
SamHeader BamHeader::ToSamHeader(void) const {
    return m_header;
}

// returns SAM-formatted string of header data
string BamHeader::ToString(void) const {
    return m_header.ToString();
}
