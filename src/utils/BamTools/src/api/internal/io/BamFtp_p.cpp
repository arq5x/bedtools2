// ***************************************************************************
// BamFtp_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 8 December 2011 (DB)
// ---------------------------------------------------------------------------
// Provides reading/writing of BAM files on FTP server
// ***************************************************************************

#include "api/BamAux.h"
#include "api/internal/io/BamFtp_p.h"
#include "api/internal/io/TcpSocket_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <cctype>
#include <cstdlib>
#include <sstream>
#include <vector>
using namespace std;

namespace BamTools {
namespace Internal {

// -----------
// constants
// -----------

static const uint16_t FTP_PORT          = 21;
static const string   FTP_PREFIX        = "ftp://";
static const size_t   FTP_PREFIX_LENGTH = 6;
static const string   FTP_NEWLINE       = "\r\n";

static const string DEFAULT_USER = "anonymous";
static const string DEFAULT_PASS = "anonymous@";

static const string ABOR_CMD = "ABOR";
static const string USER_CMD = "USER";
static const string PASS_CMD = "PASS";
static const string PASV_CMD = "PASV";
static const string REIN_CMD = "REIN";
static const string REST_CMD = "REST";
static const string RETR_CMD = "RETR";
static const string TYPE_CMD = "TYPE";

static const char CMD_SEPARATOR  = ' ';
static const char HOST_SEPARATOR = '/';
static const char IP_SEPARATOR   = '.';

static const char MULTILINE_CONTINUE = '-';

static const char PASV_REPLY_PREFIX    = '(';
static const char PASV_REPLY_SEPARATOR = ',';
static const char PASV_REPLY_SUFFIX    = ')';

// -----------------
// utility methods
// -----------------

static inline
vector<string> split(const string& source, const char delim) {

    stringstream ss(source);
    string field;
    vector<string> fields;

    while ( getline(ss, field, delim) )
        fields.push_back(field);
    return fields;
}

static inline
bool startsWith(const string& source, const string& pattern) {
    return ( source.find(pattern) == 0 );
}

static inline
string toLower(const string& s) {
    string out;
    const size_t sSize = s.size();
    out.resize(sSize);
    for ( size_t i = 0; i < sSize; ++i )
        out[i] = tolower(s[i]);
    return out;
}

} // namespace Internal
} // namespace BamTools

// -----------------------
// BamFtp implementation
// -----------------------

BamFtp::BamFtp(const string& url)
    : IBamIODevice()
    , m_commandSocket(new TcpSocket)
    , m_dataSocket(new TcpSocket)
    , m_port(FTP_PORT)
    , m_dataPort(0)
    , m_username(DEFAULT_USER)
    , m_password(DEFAULT_PASS)
    , m_isUrlParsed(false)
    , m_filePosition(-1)
{
    ParseUrl(url);
}

BamFtp::~BamFtp(void) {

    // close connection & clean up
    Close();
    if ( m_commandSocket )
        delete m_commandSocket;
    if ( m_dataSocket )
        delete m_dataSocket;
}

void BamFtp::Close(void) {

    // disconnect socket
    m_commandSocket->DisconnectFromHost();
    m_dataSocket->DisconnectFromHost();

    // reset state - necessary??
    m_isUrlParsed = false;
    m_filePosition = -1;
    m_username = DEFAULT_USER;
    m_password = DEFAULT_PASS;
    m_dataHostname.clear();
    m_dataPort = 0;
}

bool BamFtp::ConnectCommandSocket(void) {

    BT_ASSERT_X(m_commandSocket, "null command socket?");

    // connect to FTP server
    if ( !m_commandSocket->ConnectToHost(m_hostname, m_port, m_mode) ) {
        SetErrorString("BamFtp::ConnectCommandSocket", "could not connect to host - ");
        return false;
    }

    // receive initial reply from host
    if ( !ReceiveReply() ) {
        Close();
        return false;
    }

    // send USER command
    string userCommand = USER_CMD + CMD_SEPARATOR + m_username + FTP_NEWLINE;
    if ( !SendCommand(userCommand, true) ) {
        Close();
        return false;
    }

    // send PASS command
    string passwordCommand = PASS_CMD + CMD_SEPARATOR + m_password + FTP_NEWLINE;
    if ( !SendCommand(passwordCommand, true) ) {
        Close();
        return false;
    }

    // send TYPE command
    string typeCommand = TYPE_CMD + CMD_SEPARATOR + 'I' + FTP_NEWLINE;
    if ( !SendCommand(typeCommand, true) ) {
        Close();
        return false;
    }

    // return success
    return true;
}

bool BamFtp::ConnectDataSocket(void) {

    // failure if can't connect to command socket first
    if ( !m_commandSocket->IsConnected() ) {
        if ( !ConnectCommandSocket() )
            return false;
    }

    // make sure we're starting with a fresh data channel
    if ( m_dataSocket->IsConnected() )
        m_dataSocket->DisconnectFromHost();

    // send passive connection command
    const string passiveCommand = PASV_CMD + FTP_NEWLINE;
    if ( !SendCommand(passiveCommand, true) ) {
        // TODO: set error string
        return false;
    }

    // retrieve passive connection port
    if ( !ParsePassiveResponse() ) {
        // TODO: set error string
        return false;
    }

    // set up restart command (tell server where to start fetching bytes from)
    if ( m_filePosition >= 0 ) {

        stringstream fpStream("");
        fpStream << m_filePosition;
        string restartCommand = REST_CMD + CMD_SEPARATOR + fpStream.str() + FTP_NEWLINE;
        if ( !SendCommand(restartCommand, true) ) {
            // TODO: set error string
            return false;
        }
    }

    // main file retrieval request
    string retrieveCommand = RETR_CMD + CMD_SEPARATOR + m_filename + FTP_NEWLINE;
    if ( !SendCommand(retrieveCommand, false) ) {
        // TODO: set error string
        return false;
    }

    // make data channel connection
    if ( !m_dataSocket->ConnectToHost(m_dataHostname, m_dataPort) ) {
        // TODO: set error string
        return false;
    }

    // fetch intial reply from server
    if ( !ReceiveReply() ) {
        // TODO: set error string
        m_dataSocket->DisconnectFromHost();
        return false;
    }

    // make sure we have reply code 150 (all good)
    if ( !startsWith(m_response, "150") ) {
        // TODO: set error string
        m_dataSocket->DisconnectFromHost();
        return false;
    }

    // return success
    return true;
}

bool BamFtp::IsOpen(void) const {
    return IBamIODevice::IsOpen() && m_isUrlParsed;
}

bool BamFtp::IsRandomAccess(void) const {
    return true;
}

bool BamFtp::Open(const IBamIODevice::OpenMode mode) {

    // BamFtp only supports read-only access
    if ( mode != IBamIODevice::ReadOnly ) {
        SetErrorString("BamFtp::Open", "writing on this device is not supported");
        return false;
    }

    // initialize basic valid state
    m_mode = mode;
    m_filePosition = 0;

    // attempt connection to command & data sockets
    return ( ConnectCommandSocket() && ConnectDataSocket() );
}

bool BamFtp::ParsePassiveResponse(void) {

    // fail if empty
    if ( m_response.empty() )
        return false;

    // find parentheses
    const size_t leftParenFound  = m_response.find(PASV_REPLY_PREFIX);
    const size_t rightParenFound = m_response.find(PASV_REPLY_SUFFIX);
    if ( leftParenFound == string::npos || rightParenFound == string::npos )
        return false;

    // grab everything between ( should be "h1,h2,h3,h4,p1,p2" )
    string::const_iterator responseBegin = m_response.begin();
    const string hostAndPort(responseBegin+leftParenFound+1, responseBegin+rightParenFound);

    // parse into string fields
    vector<string> fields = split(hostAndPort, PASV_REPLY_SEPARATOR);
    if ( fields.size() != 6 )
        return false;

    // fetch passive connection IP
    m_dataHostname = fields[0] + IP_SEPARATOR +
                     fields[1] + IP_SEPARATOR +
                     fields[2] + IP_SEPARATOR +
                     fields[3];

    // fetch passive connection port
    const uint8_t portUpper = static_cast<uint8_t>(atoi(fields[4].c_str()));
    const uint8_t portLower = static_cast<uint8_t>(atoi(fields[5].c_str()));
    m_dataPort = ( portUpper<<8 ) + portLower;

    // return success
    return true;
}

void BamFtp::ParseUrl(const string& url) {

    // clear flag to start
    m_isUrlParsed = false;

    // make sure url starts with "ftp://", case-insensitive
    string tempUrl(url);
    toLower(tempUrl);
    const size_t prefixFound = tempUrl.find(FTP_PREFIX);
    if ( prefixFound == string::npos )
        return;

    // find end of host name portion (first '/' hit after the prefix)
    const size_t firstSlashFound = tempUrl.find(HOST_SEPARATOR, FTP_PREFIX_LENGTH);
    if ( firstSlashFound == string::npos ) {
        ;  // no slash found... no filename given along with host?
    }

    // fetch hostname
    string hostname = tempUrl.substr(FTP_PREFIX_LENGTH, (firstSlashFound - FTP_PREFIX_LENGTH));
    m_hostname = hostname;
    m_port = FTP_PORT;

    // store remainder of URL as filename (must be non-empty)
    string filename = tempUrl.substr(firstSlashFound);
    if ( filename.empty() )
        return;
    m_filename = filename;

    // set parsed OK flag
    m_isUrlParsed = true;
}

int64_t BamFtp::Read(char* data, const unsigned int numBytes) {

    // if BamHttp not in a valid state
    if ( !IsOpen() )
        return -1;

    // read until hit desired @numBytes
    int64_t bytesReadSoFar = 0;
    while ( bytesReadSoFar < numBytes ) {

        // calculate number of bytes we're going to try to read this iteration
        const size_t remainingBytes = ( numBytes - bytesReadSoFar );

        // if either disconnected somehow, or (more likely) we have seeked since last read
        if ( !m_dataSocket->IsConnected() ) {
            if ( !ConnectDataSocket() ) {
                // TODO: set error string
                return -1;
            }
        }

        // read bytes from data socket
        const int64_t socketBytesRead = ReadDataSocket(data+bytesReadSoFar, remainingBytes);
        if ( socketBytesRead < 0 ) // error
            return -1;
        else if ( socketBytesRead == 0 ) // EOF
            return bytesReadSoFar;
        bytesReadSoFar += socketBytesRead;
        m_filePosition += socketBytesRead;
    }

    // return actual number bytes successfully read
    return bytesReadSoFar;
}

int64_t BamFtp::ReadCommandSocket(char* data, const unsigned int maxNumBytes) {
    return m_commandSocket->Read(data, maxNumBytes);
}

int64_t BamFtp::ReadDataSocket(char* data, const unsigned int maxNumBytes) {
    return m_dataSocket->Read(data, maxNumBytes);
}

bool BamFtp::ReceiveReply(void) {

    // failure if not connected
    if ( !m_commandSocket->IsConnected() ) {
        SetErrorString("BamFtp::ReceiveReply()", "command socket not connected");
        return false;
    }

    m_response.clear();

    // read header data (& discard for now)
    bool headerEnd = false;
    while ( !headerEnd ) {

        const string headerLine = m_commandSocket->ReadLine();
        m_response += headerLine;

        // if line is of form 'xyz ', quit reading lines
        if ( (headerLine.length() >= 4 ) &&
             isdigit(headerLine[0]) &&
             isdigit(headerLine[1]) &&
             isdigit(headerLine[2]) &&
             ( headerLine[3] != MULTILINE_CONTINUE )
           )
        {
            headerEnd = true;
        }
    }

    // return success, depending on response
    if ( m_response.empty() ) {
        SetErrorString("BamFtp::ReceiveReply", "error reading server reply");
        return false;
    }
    return true;
}

bool BamFtp::Seek(const int64_t& position, const int origin) {

    // if FTP device not in a valid state
    if ( !IsOpen() ) {
        // TODO: set error string
        return false;
    }

    // ----------------------
    // UGLY !! but works??
    // ----------------------
    // disconnect from server
    m_dataSocket->DisconnectFromHost();
    m_commandSocket->DisconnectFromHost();

    // update file position & return success
    if ( origin == SEEK_CUR )
        m_filePosition += position;
    else if ( origin == SEEK_SET)
        m_filePosition = position;
    else {
        // TODO: set error string
        return false;
    }
    return true;
}

bool BamFtp::SendCommand(const string& command, bool waitForReply) {

    // failure if not connected
    if ( !m_commandSocket->IsConnected() ) {
        SetErrorString("BamFtp::SendCommand", "command socket not connected");
        return false;
    }

    // write command to 'command socket'
    if ( WriteCommandSocket(command.c_str(), command.length()) == -1 ) {
        SetErrorString("BamFtp::SendCommand", "error writing to socket");
        // get actual error from command socket??
        return false;
    }

    // if we sent a command that receives a response
    if ( waitForReply )
        return ReceiveReply();

    // return success
    return true;
}

int64_t BamFtp::Tell(void) const {
    return ( IsOpen() ? m_filePosition : -1 );
}

int64_t BamFtp::Write(const char* data, const unsigned int numBytes) {
    (void)data;
    (void)numBytes;
    BT_ASSERT_X(false, "BamFtp::Write : write-mode not supported on this device");
    SetErrorString("BamFtp::Write", "write-mode not supported on this device");
    return -1;
}

int64_t BamFtp::WriteCommandSocket(const char* data, const unsigned int numBytes) {
    if ( !m_commandSocket->IsConnected() )
        return -1;
    m_commandSocket->ClearBuffer();
    return m_commandSocket->Write(data, numBytes);
}

int64_t BamFtp::WriteDataSocket(const char* data, const unsigned int numBytes) {
    (void)data;
    (void)numBytes;
    BT_ASSERT_X(false, "BamFtp::WriteDataSocket: write-mode not supported on this device");
    SetErrorString("BamFtp::Write", "write-mode not supported on this device");
    return -1;
}
