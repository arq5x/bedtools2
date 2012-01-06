// ***************************************************************************
// TcpSocket_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 5 January 2012 (DB)
// ---------------------------------------------------------------------------
// Provides basic TCP I/O interface
// ***************************************************************************

#include "api/internal/io/ByteArray_p.h"
#include "api/internal/io/TcpSocket_p.h"
#include "api/internal/io/TcpSocketEngine_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <algorithm>
#include <climits>
#include <sstream>
#include <vector>
using namespace std;

// ------------------------------------
// static utility methods & constants
// ------------------------------------

namespace BamTools {
namespace Internal {

// constants
static const size_t DEFAULT_BUFFER_SIZE = 0x4000;

} // namespace Internal
} // namespace BamTools

// --------------------------
// TcpSocket implementation
// --------------------------

TcpSocket::TcpSocket(void)
    : m_mode(IBamIODevice::NotOpen)
// , m_localPort(0)
    , m_remotePort(0)
    , m_engine(0)
    , m_cachedSocketDescriptor(-1)
    , m_readBuffer(DEFAULT_BUFFER_SIZE)
    , m_error(TcpSocket::UnknownSocketError)
    , m_state(TcpSocket::UnconnectedState)
{ }

TcpSocket::~TcpSocket(void) {
    if ( m_state == TcpSocket::ConnectedState )
        DisconnectFromHost();
}

size_t TcpSocket::BufferBytesAvailable(void) const {
    return m_readBuffer.Size();
}

bool TcpSocket::CanReadLine(void) const {
    return m_readBuffer.CanReadLine();
}

void TcpSocket::ClearBuffer(void) {
    m_readBuffer.Clear();
}

bool TcpSocket::ConnectImpl(const HostInfo& hostInfo,
                            const std::string& port,
                            IBamIODevice::OpenMode mode)
{
    // skip if we're already connected
    if ( m_state == TcpSocket::ConnectedState ) {
        m_error = TcpSocket::SocketResourceError;
        m_errorString = "socket already connected";
        return false;
    }

    // reset socket state
    m_hostName = hostInfo.HostName();
    m_mode = mode;
    m_state = TcpSocket::UnconnectedState;
    m_error = TcpSocket::UnknownSocketError;
// m_localPort = 0;
    m_remotePort = 0;
// m_localAddress.Clear();
    m_remoteAddress.Clear();
    m_readBuffer.Clear();

    // fetch candidate addresses for requested host
    vector<HostAddress> addresses = hostInfo.Addresses();
    if ( addresses.empty() ) {
        m_error = TcpSocket::HostNotFoundError;
        m_errorString = "no IP addresses found for host";
        return false;
    }

    // convert port string to integer
    stringstream ss(port);
    uint16_t portNumber(0);
    ss >> portNumber;

    // iterate through adddresses
    vector<HostAddress>::const_iterator addrIter = addresses.begin();
    vector<HostAddress>::const_iterator addrEnd = addresses.end();
    for ( ; addrIter != addrEnd; ++addrIter) {
        const HostAddress& addr = (*addrIter);

        // try to initialize socket engine with this address
        if ( !InitializeSocketEngine(addr.GetProtocol()) ) {
            // failure to initialize is OK here
            // we'll just try the next available address
            continue;
        }

        // attempt actual connection
        if ( m_engine->Connect(addr, portNumber) ) {

            // if connection successful, update our state & return true
            m_mode = mode;
// m_localAddress = m_engine->GetLocalAddress();
// m_localPort = m_engine->GetLocalPort();
            m_remoteAddress = m_engine->GetRemoteAddress();
            m_remotePort = m_engine->GetRemotePort();
            m_cachedSocketDescriptor = m_engine->GetSocketDescriptor();
            m_state = TcpSocket::ConnectedState;
            return true;
        }
    }

    // if we get here, no connection could be made
    m_error = TcpSocket::HostNotFoundError;
    m_errorString = "could not connect to any host addresses";
    return false;
}

bool TcpSocket::ConnectToHost(const string& hostName,
                              uint16_t port,
                              IBamIODevice::OpenMode mode)
{
    stringstream ss("");
    ss << port;
    return ConnectToHost(hostName, ss.str(), mode);

}

bool TcpSocket::ConnectToHost(const string& hostName,
                              const string& port,
                              IBamIODevice::OpenMode mode)
{
    // create new address object with requested host name
    HostAddress hostAddress;
    hostAddress.SetAddress(hostName);

    HostInfo info;
    // if host name was IP address ("x.x.x.x" or IPv6 format)
    // otherwise host name was 'plain-text' ("www.foo.bar")
    // we need to look up IP address(es)
    if ( hostAddress.HasIPAddress() )
        info.SetAddresses( vector<HostAddress>(1, hostAddress) );
    else
        info = HostInfo::Lookup(hostName, port);

    // attempt connection on requested port
    return ConnectImpl(info, port, mode);
}

void TcpSocket::DisconnectFromHost(void) {

    // close socket engine & delete
    if ( m_state == TcpSocket::ConnectedState )
        ResetSocketEngine();

    // reset connection state
// m_localPort = 0;
    m_remotePort = 0;
// m_localAddress.Clear();
    m_remoteAddress.Clear();
    m_hostName.clear();
    m_cachedSocketDescriptor = -1;

    // for future, make sure there's outgoing data that needs to be flushed
    m_readBuffer.Clear();
}

TcpSocket::SocketError TcpSocket::GetError(void) const {
    return m_error;
}

std::string TcpSocket::GetErrorString(void) const {
    return m_errorString;
}

std::string TcpSocket::GetHostName(void) const {
    return m_hostName;
}

//HostAddress TcpSocket::GetLocalAddress(void) const {
// return m_localAddress;
//}

//uint16_t TcpSocket::GetLocalPort(void) const {
// return m_localPort;
//}

HostAddress TcpSocket::GetRemoteAddress(void) const {
    return m_remoteAddress;
}

uint16_t TcpSocket::GetRemotePort(void) const {
    return m_remotePort;
}

TcpSocket::SocketState TcpSocket::GetState(void) const {
    return m_state;
}

bool TcpSocket::InitializeSocketEngine(HostAddress::NetworkProtocol protocol) {
    ResetSocketEngine();
    m_engine = new TcpSocketEngine;
    return m_engine->Initialize(protocol);
}

bool TcpSocket::IsConnected(void) const {
    if ( m_engine == 0 )
        return false;
    return ( m_engine->IsValid() && (m_state == TcpSocket::ConnectedState) );
}

// may be read in a look until desired data amount has been read
// returns: number of bytes read, or -1 if error
int64_t TcpSocket::Read(char* data, const unsigned int numBytes) {

    // if we have data in buffer, just return it
    if ( !m_readBuffer.IsEmpty() ) {
        const size_t bytesRead = m_readBuffer.Read(data, numBytes);
        return static_cast<int64_t>(bytesRead);
    }

    // otherwise, we'll need to fetch data from socket
    // first make sure we have a valid socket engine
    if ( m_engine == 0 ) {
        // TODO: set error string/state?
        return -1;
    }

    // fetch data from socket, return 0 for success, -1 for failure
    // since this should be called in a loop,
    // we'll pull the actual bytes from the buffer on next iteration
    const int64_t socketBytesRead = ReadFromSocket();
    if ( socketBytesRead < 0 ) {
        // TODO: set error string/state ?
        return -1;
    }

    // we should have data now in buffer, try to fetch requested amount
    // if nothing in buffer, we will return 0 bytes read (signals EOF reached)
    const size_t numBytesRead = m_readBuffer.Read(data, numBytes);
    return static_cast<int64_t>(numBytesRead);
}

int64_t TcpSocket::ReadFromSocket(void) {

    // check for any socket engine errors
    if ( !m_engine->IsValid() ) {
        m_errorString = "TcpSocket::ReadFromSocket - socket disconnected";
        ResetSocketEngine();
        return -1;
    }

    // wait for ready read
    bool timedOut;
    const bool isReadyRead = m_engine->WaitForRead(5000, &timedOut);

    // if not ready
    if ( !isReadyRead ) {

        // if we simply timed out
        if ( timedOut ) {
            // TODO: get add'l error info from engine ?
            m_errorString = "TcpSocket::ReadFromSocket - timed out waiting for ready read";
        }

        // otherwise, there was some other error
        else {
            // TODO: get add'l error info from engine ?
            m_errorString = "TcpSocket::ReadFromSocket - encountered error while waiting for ready read";
        }

        // return failure
        return -1;
    }

    // get number of bytes available from socket
    const int64_t bytesToRead = m_engine->NumBytesAvailable();
    if ( bytesToRead < 0 ) {
        // TODO: get add'l error info from engine ?
        m_errorString = "TcpSocket::ReadFromSocket - encountered error while determining numBytesAvailable";
        return -1;
    }

    // make space in buffer & read from socket
    char* buffer = m_readBuffer.Reserve(bytesToRead);
    const int64_t numBytesRead = m_engine->Read(buffer, bytesToRead);
    if ( numBytesRead == -1 ) {
        // TODO: get add'l error info from engine ?
        m_errorString = "TcpSocket::ReadFromSocket - encountered error while reading bytes";
    }

    // return number of bytes actually read
    return numBytesRead;
}

string TcpSocket::ReadLine(int64_t max) {

    // prep result byte buffer
    ByteArray result;
    size_t bufferMax = ((max > static_cast<int64_t>(UINT_MAX))
                        ? UINT_MAX : static_cast<size_t>(max));
    result.Resize(bufferMax);

    // read data
    int64_t readBytes(0);
    if ( result.Size() == 0 ) {

        if ( bufferMax == 0 )
            bufferMax = UINT_MAX;

        result.Resize(1);

        int64_t readResult;
        do {
            result.Resize( static_cast<size_t>(min(bufferMax, result.Size() + DEFAULT_BUFFER_SIZE)) );
            readResult = ReadLine(result.Data()+readBytes, result.Size()-readBytes);
            if ( readResult > 0 || readBytes == 0 )
                readBytes += readResult;
        } while ( readResult == DEFAULT_BUFFER_SIZE && result[static_cast<size_t>(readBytes-1)] != '\n' );

    } else
        readBytes = ReadLine(result.Data(), result.Size());

    // clean up byte buffer
    if ( readBytes <= 0 )
        result.Clear();
    else
        result.Resize(static_cast<size_t>(readBytes));

    // return byte buffer as string
    return string( result.ConstData(), result.Size() );
}

int64_t TcpSocket::ReadLine(char* dest, size_t max) {

    // wait for buffer to contain line contents
    if ( !WaitForReadLine() ) {
        m_errorString = "TcpSocket::ReadLine - error waiting for read line";
        return -1;
    }

    // leave room for null term
    if ( max < 2 )
        return -1;
    --max;

    // read from buffer, handle newlines
    int64_t readSoFar = m_readBuffer.ReadLine(dest, max);
    if ( readSoFar && dest[readSoFar-1] == '\n' ) {

        // adjust for windows-style '\r\n'
        if ( readSoFar > 1 && dest[readSoFar-2] == '\r') {
            --readSoFar;
            dest[readSoFar-1] = '\n';
        }
    }

    // null terminate & return number of bytes read
    dest[readSoFar] = '\0';
    return readSoFar;
}

void TcpSocket::ResetSocketEngine(void) {

    // shut down socket engine
    if ( m_engine ) {
        m_engine->Close();
        delete m_engine;
        m_engine = 0;
    }

    // reset our state & cached socket handle
    m_state = TcpSocket::UnconnectedState;
    m_cachedSocketDescriptor = -1;
}

bool TcpSocket::WaitForReadLine(void) {

    // wait until we can read a line (will return immediately if already capable)
    while ( !CanReadLine() ) {
        if ( !ReadFromSocket() )
            return false;
    }

    // if we get here, success
    return true;
}

int64_t TcpSocket::Write(const char* data, const unsigned int numBytes) {

    // single-shot attempt at write (not buffered, just try to shove the data through socket)
    // this method purely exists to send 'small' HTTP requests/FTP commands from client to server

    // wait for our socket to be write-able
    bool timedOut;
    const bool isReadyWrite = m_engine->WaitForWrite(3000, &timedOut);

    // if ready, return number of bytes written
    if ( isReadyWrite )
        return m_engine->Write(data, numBytes);

    // otherwise, socket not ready for writing
    // set error string depending on reason & return failure
    if ( !timedOut ) {
        // TODO: get add'l error info from engine ??
        m_errorString = "TcpSocket::Write - timed out waiting for ready-write";
    }
    else {
        // TODO: get add'l error info from engine ??
        m_errorString = "TcpSocket::Write - error encountered while waiting for ready-write";
    }
    return -1;
}

