// ***************************************************************************
// BamFile_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 November 2011 (DB)
// ---------------------------------------------------------------------------
// Provides BAM file-specific IO behavior
// ***************************************************************************

#include "api/internal/io/BamFile_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <cstdio>
#include <iostream>
using namespace std;

BamFile::BamFile(const string& filename)
    : ILocalIODevice()
    , m_filename(filename)
{ }

BamFile::~BamFile(void) { }

void BamFile::Close(void) {
    if ( IsOpen() ) {
        m_filename.clear();
        ILocalIODevice::Close();
    }
}

bool BamFile::IsRandomAccess(void) const {
    return true;
}

bool BamFile::Open(const IBamIODevice::OpenMode mode) {

    // make sure we're starting with a fresh file stream
    Close();

    // attempt to open FILE* depending on requested openmode
    if ( mode == IBamIODevice::ReadOnly )
        m_stream = fopen(m_filename.c_str(), "rb");
    else if ( mode == IBamIODevice::WriteOnly )
        m_stream = fopen(m_filename.c_str(), "wb");
    else if ( mode == IBamIODevice::ReadWrite )
        m_stream = fopen(m_filename.c_str(), "w+b");
    else {
        SetErrorString("BamFile::Open", "unknown open mode requested");
        return false;
    }

    // check that we obtained a valid FILE*
    if ( m_stream == 0 ) {
        const string message_base = string("could not open file handle for ");
        const string message = message_base + ( (m_filename.empty()) ? "empty filename" : m_filename );
        SetErrorString("BamFile::Open", message);
        return false;
    }

    // store current IO mode & return success
    m_mode = mode;
    return true;
}

bool BamFile::Seek(const int64_t& position, const int origin) {
    BT_ASSERT_X( m_stream, "BamFile::Seek() - null stream" );
    return ( fseek64(m_stream, position, origin) == 0 );
}
