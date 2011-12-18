// ***************************************************************************
// ILocalIODevice_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides shared behavior for files & pipes
// ***************************************************************************

#include "api/internal/io/ILocalIODevice_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <cstdio>
using namespace std;

ILocalIODevice::ILocalIODevice(void)
    : IBamIODevice()
    , m_stream(0)
{ }

ILocalIODevice::~ILocalIODevice(void) {
    Close();
}

void ILocalIODevice::Close(void) {

    // skip if not open
    if ( !IsOpen() )
        return;

    // flush & close FILE*
    fflush(m_stream);
    fclose(m_stream);
    m_stream = 0;

    // reset other device state
    m_mode = IBamIODevice::NotOpen;
}

int64_t ILocalIODevice::Read(char* data, const unsigned int numBytes) {
    BT_ASSERT_X( m_stream, "ILocalIODevice::Read: trying to read from null stream" );
    BT_ASSERT_X( (m_mode == IBamIODevice::ReadOnly), "ILocalIODevice::Read: device not in read-only mode");
    return static_cast<int64_t>( fread(data, sizeof(char), numBytes, m_stream) );
}

int64_t ILocalIODevice::Tell(void) const {
    BT_ASSERT_X( m_stream, "ILocalIODevice::Tell: trying to get file position fromnull stream" );
    return ftell64(m_stream);
}

int64_t ILocalIODevice::Write(const char* data, const unsigned int numBytes) {
    BT_ASSERT_X( m_stream, "ILocalIODevice::Write: tryint to write to null stream" );
    BT_ASSERT_X( (m_mode == IBamIODevice::WriteOnly), "ILocalIODevice::Write: device not in write-only mode" );
    return static_cast<int64_t>( fwrite(data, sizeof(char), numBytes, m_stream) );
}
