// ***************************************************************************
// TcpSocketEngine_unix_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 15 November 2011 (DB)
// ---------------------------------------------------------------------------
// Provides low-level implementation of TCP I/O for all UNIX-like systems
// ***************************************************************************

#include "api/internal/io/TcpSocketEngine_p.h"
#include "api/internal/io/NetUnix_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <cerrno>
#include <ctime>
#include <iostream>
using namespace std;

// ------------------------
// static utility methods
// ------------------------

namespace BamTools {
namespace Internal {

} // namespace Internal
} // namespace BamTools

// --------------------------------
// TcpSocketEngine implementation
// --------------------------------

void TcpSocketEngine::nativeClose(void) {
    close(m_socketDescriptor);
}

bool TcpSocketEngine::nativeConnect(const HostAddress& address, const uint16_t port) {

    // setup connection parameters from address/port
    sockaddr_in  sockAddrIPv4;
    sockaddr_in6 sockAddrIPv6;
    sockaddr*    sockAddrPtr  = 0;
    BT_SOCKLEN_T sockAddrSize = 0;

    // IPv6
    if ( address.GetProtocol() == HostAddress::IPv6Protocol ) {

        memset(&sockAddrIPv6, 0, sizeof(sockAddrIPv6));
        sockAddrIPv6.sin6_family = AF_INET6;
        sockAddrIPv6.sin6_port   = htons(port);

        IPv6Address ip6 = address.GetIPv6Address();
        memcpy(&sockAddrIPv6.sin6_addr.s6_addr, &ip6, sizeof(ip6));

        sockAddrSize = sizeof(sockAddrIPv6);
        sockAddrPtr  = (sockaddr*)&sockAddrIPv6;
    }

    // IPv4
    else if ( address.GetProtocol() == HostAddress::IPv4Protocol ) {

        memset(&sockAddrIPv4, 0, sizeof(sockAddrIPv4));
        sockAddrIPv4.sin_family      = AF_INET;
        sockAddrIPv4.sin_port        = htons(port);
        sockAddrIPv4.sin_addr.s_addr = htonl(address.GetIPv4Address());

        sockAddrSize = sizeof(sockAddrIPv4);
        sockAddrPtr  = (sockaddr*)&sockAddrIPv4;
    }

    // unknown (should be unreachable)
    else BT_ASSERT_X(false, "TcpSocketEngine::nativeConnect() : unknown network protocol");

    // attempt connection
    int connectResult = connect(m_socketDescriptor, sockAddrPtr, sockAddrSize);

    // if failed, handle error
    if ( connectResult == -1 ) {

        // ensure state is set before checking errno
        m_socketState = TcpSocket::UnconnectedState;

        // set error type/message depending on errno
        switch ( errno ) { // <-- potential thread issues later? but can't get error type from connectResult

            case EISCONN:
                m_socketState = TcpSocket::ConnectedState; // socket was already connected
                break;
            case ECONNREFUSED:
            case EINVAL:
                m_socketError = TcpSocket::ConnectionRefusedError;
                m_errorString = "connection refused";
                break;
            case ETIMEDOUT:
                m_socketError = TcpSocket::NetworkError;
                m_errorString = "connection timed out";
                break;
            case EHOSTUNREACH:
                m_socketError = TcpSocket::NetworkError;
                m_errorString = "host unreachable";
                break;
            case ENETUNREACH:
                m_socketError = TcpSocket::NetworkError;
                m_errorString = "network unreachable";
                break;
            case EADDRINUSE:
                m_socketError = TcpSocket::SocketResourceError;
                m_errorString = "address already in use";
                break;
            case EACCES:
            case EPERM:
                m_socketError = TcpSocket::SocketAccessError;
                m_errorString = "permission denied";
                break;
            default:
                break;
        }

        // double check that we're not in 'connected' state; if so, return failure
        if ( m_socketState != TcpSocket::ConnectedState )
            return false;
    }

    // otherwise, we should be good
    // update state & return success
    m_socketState = TcpSocket::ConnectedState;
    return true;
}

bool TcpSocketEngine::nativeCreateSocket(HostAddress::NetworkProtocol protocol) {

    // get protocol value for requested protocol type
    const int protocolNum = ( (protocol == HostAddress::IPv6Protocol) ? AF_INET6
                                                                      : AF_INET );

    // attempt to create socket
    int socketFd = socket(protocolNum, SOCK_STREAM, IPPROTO_TCP);

    // if we fetched an invalid socket descriptor
    if ( socketFd <= 0 ) {

        // see what error we got
        switch ( errno ) {
            case EPROTONOSUPPORT:
            case EAFNOSUPPORT:
            case EINVAL:
                m_socketError = TcpSocket::UnsupportedSocketOperationError;
                m_errorString = "protocol not supported";
                break;
            case ENFILE:
            case EMFILE:
            case ENOBUFS:
            case ENOMEM:
                m_socketError = TcpSocket::SocketResourceError;
                m_errorString = "out of resources";
                break;
            case EACCES:
                m_socketError = TcpSocket::SocketAccessError;
                m_errorString = "permission denied";
                break;
            default:
                break;
        }

        // return failure
        return false;
    }

    // otherwise, store our socket FD & return success
    m_socketDescriptor = socketFd;
    return true;
}

int64_t TcpSocketEngine::nativeNumBytesAvailable(void) const {

    // fetch number of bytes, return 0 on error
    int numBytes(0);
    if ( ioctl(m_socketDescriptor, FIONREAD, (char*)&numBytes) < 0 )
        return -1;
    return static_cast<int64_t>(numBytes);
}

int64_t TcpSocketEngine::nativeRead(char* dest, size_t max) {
    const ssize_t ret = read(m_socketDescriptor, dest, max);
    return static_cast<int64_t>(ret);
}

// negative value for msecs will block (forever) until ready
int TcpSocketEngine::nativeSelect(int msecs, bool isRead) const {

    // set up FD set
    fd_set fds;
    FD_ZERO(&fds);
    FD_SET(m_socketDescriptor, &fds);

    // setup our timeout
    timeval tv;
    tv.tv_sec  = msecs / 1000;
    tv.tv_usec = (msecs % 1000) * 1000;

    // do 'select'
    if ( isRead )
        return select(m_socketDescriptor + 1, &fds, 0, 0, (msecs < 0 ? 0 : &tv));
    else
        return select(m_socketDescriptor + 1, 0, &fds, 0, (msecs < 0 ? 0 : &tv));
}

int64_t TcpSocketEngine::nativeWrite(const char* data, size_t length) {
    const ssize_t writtenBytes = write(m_socketDescriptor, data, length);
    return static_cast<int64_t>(writtenBytes);
}
