// ***************************************************************************
// BamHttp_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 8 December 2011 (DB)
// ---------------------------------------------------------------------------
// Provides reading/writing of BAM files on HTTP server
// ***************************************************************************

#include "api/BamAux.h"
#include "api/internal/io/BamHttp_p.h"
#include "api/internal/io/HttpHeader_p.h"
#include "api/internal/io/TcpSocket_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <cassert>
#include <cctype>
#include <algorithm>
#include <sstream>
using namespace std;

namespace BamTools {
namespace Internal {

// -----------
// constants
// -----------

static const string HTTP_PORT   = "80";
static const string HTTP_PREFIX = "http://";
static const size_t HTTP_PREFIX_LENGTH = 7;

static const string DOUBLE_NEWLINE = "\n\n";

static const string GET_METHOD   = "GET";
static const string HOST_HEADER  = "Host";
static const string RANGE_HEADER = "Range";
static const string BYTES_PREFIX = "bytes=";

static const char HOST_SEPARATOR  = '/';
static const char PROXY_SEPARATOR = ':';

// -----------------
// utility methods
// -----------------

static inline
bool endsWith(const string& source, const string& pattern) {
    return ( source.find(pattern) == (source.length() - pattern.length()) );
}

static inline
string toLower(const string& s) {
    string out;
    const size_t sSize = s.size();
    out.reserve(sSize);
    for ( size_t i = 0; i < sSize; ++i )
        out[i] = tolower(s[i]);
    return out;
}

} // namespace Internal
} // namespace BamTools

// ------------------------
// BamHttp implementation
// ------------------------

BamHttp::BamHttp(const string& url)
    : IBamIODevice()
    , m_socket(new TcpSocket)
    , m_port(HTTP_PORT)
    , m_request(0)
    , m_response(0)
    , m_isUrlParsed(false)
    , m_filePosition(-1)
    , m_endRangeFilePosition(-1)
{
    ParseUrl(url);
}

BamHttp::~BamHttp(void) {

    // close connection & clean up
    Close();
    if ( m_socket )
        delete m_socket;
}

void BamHttp::Close(void) {

    // disconnect socket
    m_socket->DisconnectFromHost();

    // clean up request & response
    if ( m_request )  {
        delete m_request;
        m_request = 0;
    }
    if ( m_response ) {
        delete m_response;
        m_response = 0;
    }

    // reset state - necessary??
    m_isUrlParsed = false;
    m_filePosition = -1;
    m_endRangeFilePosition = -1;
}

bool BamHttp::ConnectSocket(void) {

    BT_ASSERT_X(m_socket, "null socket?");

    // any state checks, etc?
    if ( !m_socket->ConnectToHost(m_hostname, m_port, m_mode) ) {
        // TODO: set error string
        return false;
    }

    // attempt initial request
    m_filePosition = 0;
    m_endRangeFilePosition = -1;
    if ( !SendRequest() ) {
        // TODO: set error string
        Close();
        return false;
    }

    // wait for response from server
    if ( !ReceiveResponse() ) {
        // TODO: set error string
        Close();
        return false;
    }

    // return success
    return true;
}

bool BamHttp::EnsureSocketConnection(void) {
    if ( m_socket->IsConnected() )
        return true;
    else return ConnectSocket();
}

bool BamHttp::IsOpen(void) const {
    return IBamIODevice::IsOpen() && m_isUrlParsed;
}

bool BamHttp::IsRandomAccess(void) const {
    return true;
}

bool BamHttp::Open(const IBamIODevice::OpenMode mode) {

    // BamHttp only supports read-only access
    if ( mode != IBamIODevice::ReadOnly ) {
        SetErrorString("BamHttp::Open", "writing on this device is not supported");
        return false;
    }
    m_mode = mode;

    // attempt connection to socket
    if ( !ConnectSocket() ) {
        SetErrorString("BamHttp::Open", m_socket->GetErrorString());
        return false;
    }

    // return success
    return true;
}

void BamHttp::ParseUrl(const string& url) {

    // clear flag to start
    m_isUrlParsed = false;

    // make sure url starts with "http://", case-insensitive
    string tempUrl(url);
    toLower(tempUrl);
    const size_t prefixFound = tempUrl.find(HTTP_PREFIX);
    if ( prefixFound == string::npos )
        return;

    // find end of host name portion (first '/' hit after the prefix)
    const size_t firstSlashFound = tempUrl.find(HOST_SEPARATOR, HTTP_PREFIX_LENGTH);
    if ( firstSlashFound == string::npos ) {
        ;  // no slash found... no filename given along with host?
    }

    // fetch hostname (check for proxy port)
    string hostname = tempUrl.substr(HTTP_PREFIX_LENGTH, (firstSlashFound - HTTP_PREFIX_LENGTH));
    const size_t colonFound = hostname.find(PROXY_SEPARATOR);
    if ( colonFound != string::npos ) {
        ; // TODO: handle proxy port (later, just skip for now)
    } else {
        m_hostname = hostname;
        m_port = HTTP_PORT;
    }

    // store remainder of URL as filename (must be non-empty)
    string filename = tempUrl.substr(firstSlashFound);
    if ( filename.empty() )
        return;
    m_filename = filename;

    // set parsed OK flag
    m_isUrlParsed = true;
}

int64_t BamHttp::Read(char* data, const unsigned int numBytes) {

    // if BamHttp not in a valid state
    if ( !IsOpen() )
        return -1;

    // read until hit desired @numBytes
    int64_t bytesReadSoFar = 0;
    while ( bytesReadSoFar < numBytes ) {

        // calculate number of bytes we're going to try to read this iteration
        const size_t remainingBytes = ( numBytes - bytesReadSoFar );

        // if socket has access to entire file contents
        // i.e. we received response with full data (status code == 200)
        if ( m_endRangeFilePosition < 0 ) {

            // try to read 'remainingBytes' from socket
            const int64_t socketBytesRead = ReadFromSocket(data+bytesReadSoFar, remainingBytes);
            if ( socketBytesRead < 0 ) // error
                return -1;
            else if ( socketBytesRead == 0 ) // EOF
                return bytesReadSoFar;
            bytesReadSoFar += socketBytesRead;
            m_filePosition += socketBytesRead;
        }

        // socket has access to a range of data (might already be in buffer)
        // i.e. we received response with partial data (status code == 206)
        else {

            // there is data left from last request
            if ( m_endRangeFilePosition > m_filePosition ) {

                // try to read either the total 'remainingBytes' or
                // whatever we have remaining from last request range
                const size_t rangeRemainingBytes = m_endRangeFilePosition - m_filePosition;
                const size_t bytesToRead = std::min(remainingBytes, rangeRemainingBytes);
                const int64_t socketBytesRead = ReadFromSocket(data+bytesReadSoFar, bytesToRead);
                if ( socketBytesRead < 0 ) // error
                    return -1;
                else if ( socketBytesRead == 0 ) // EOF
                    return bytesReadSoFar;
                bytesReadSoFar += socketBytesRead;
                m_filePosition += socketBytesRead;
            }

            // otherwise, this is a 1st-time read or
            // we already read everything from the last GET request
            else {

                // request for next range
                if ( !SendRequest(remainingBytes) || !ReceiveResponse() ) {
                    Close();
                    return -1;
                }
            }
        }
    }

    // return actual number bytes successfully read
    return bytesReadSoFar;
}

int64_t BamHttp::ReadFromSocket(char* data, const unsigned int maxNumBytes) {
    return m_socket->Read(data, maxNumBytes);
}

bool BamHttp::ReceiveResponse(void) {

    // clear any prior response
    if ( m_response )
        delete m_response;

    // make sure we're connected
    if ( !EnsureSocketConnection() )
        return false;

    // fetch header, up until double new line
    string responseHeader;
    do {
        // read line & append to full header
        const string headerLine = m_socket->ReadLine();
        responseHeader += headerLine;

    } while ( !endsWith(responseHeader, DOUBLE_NEWLINE) );

    // sanity check
    if ( responseHeader.empty() ) {
        // TODO: set error string
        Close();
        return false;
    }

    // create response from header text
    m_response = new HttpResponseHeader(responseHeader);
    if ( !m_response->IsValid() ) {
        // TODO: set error string
        Close();
        return false;
    }

    // if we got range response as requested
    if ( m_response->GetStatusCode() == 206 )
        return true;

    // if we got the full file contents instead of range
    else if ( m_response->GetStatusCode() == 200 ) {

        // skip up to current file position
        RaiiBuffer tmp(0x8000);
        int64_t numBytesRead = 0;
        while ( numBytesRead < m_filePosition ) {

            const int64_t remaining = m_filePosition - numBytesRead;
            const size_t bytesToRead = static_cast<size_t>( (remaining > 0x8000) ? 0x8000 : remaining );
            const int64_t socketBytesRead = ReadFromSocket(tmp.Buffer, bytesToRead);
            if ( socketBytesRead < 0 ) { // error
                Close();
                return false;
            }
            else if ( socketBytesRead == 0 ) // EOF
                break;

            numBytesRead += socketBytesRead;
        }

        // return success
        return ( numBytesRead == m_filePosition);
    }

    // on any other reponse status
    // TODO: set error string
    Close();
    return false;
}

bool BamHttp::Seek(const int64_t& position, const int origin) {

    // if HTTP device not in a valid state
    if ( !IsOpen() ) {
        // TODO: set error string
        return false;
    }

    // discard socket's buffer contents, update positions, & return success
    m_socket->ClearBuffer();

    if ( origin == SEEK_CUR )
        m_filePosition += position;
    else if ( origin == SEEK_SET )
        m_filePosition = position;
    else {
        // TODO: set error string
        return false;
    }
    m_endRangeFilePosition = m_filePosition;
    return true;
}

bool BamHttp::SendRequest(const size_t numBytes) {

    // remove any currently active request
    if ( m_request )
        delete m_request;

    // create range string
    m_endRangeFilePosition = m_filePosition + numBytes;
    stringstream range("");
    range << BYTES_PREFIX << m_filePosition << '-' << m_endRangeFilePosition;

    // make sure we're connected
    if ( !EnsureSocketConnection() )
        return false;

    // create request
    m_request = new HttpRequestHeader(GET_METHOD, m_filename);
    m_request->SetField(HOST_HEADER,  m_hostname);
    m_request->SetField(RANGE_HEADER, range.str());

    // write request to socket
    const string requestHeader = m_request->ToString();
    const size_t headerSize    = requestHeader.size();
    return ( WriteToSocket(requestHeader.c_str(), headerSize) == headerSize );
}

int64_t BamHttp::Tell(void) const {
    return ( IsOpen() ? m_filePosition : -1 );
}

int64_t BamHttp::Write(const char* data, const unsigned int numBytes) {
    (void)data;
    (void)numBytes;
    BT_ASSERT_X(false, "BamHttp::Write : write-mode not supported on this device");
    SetErrorString("BamHttp::Write", "write-mode not supported on this device");
    return -1;
}

int64_t BamHttp::WriteToSocket(const char* data, const unsigned int numBytes) {
    if ( !m_socket->IsConnected() )
        return -1;
    m_socket->ClearBuffer();
    return m_socket->Write(data, numBytes);
}
