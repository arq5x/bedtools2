// ***************************************************************************
// TcpSocketEngine_win_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 8 December 2011 (DB)
// ---------------------------------------------------------------------------
// Provides low-level implementation of TCP I/O for all Windows systems
// ***************************************************************************

#include "api/internal/io/TcpSocketEngine_p.h"
#include "api/internal/io/NetWin_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <cstring>
#include <iostream>
#include <sstream>
using namespace std;

// --------------------------------
// TcpSocketEngine implementation
// --------------------------------

void TcpSocketEngine::nativeClose(void) {
    closesocket(m_socketDescriptor);
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

    // attempt conenction
    const int connectResult = WSAConnect(m_socketDescriptor, sockAddrPtr, sockAddrSize, 0, 0, 0, 0);

    // if failed, handle error
    if ( connectResult == SOCKET_ERROR ) {

        // ensure state is set before checking error code
        m_socketState = TcpSocket::UnconnectedState;

        // set error type/message depending on errorCode
        const int errorCode = WSAGetLastError();
        switch ( errorCode ) {
            case WSANOTINITIALISED:
                m_socketError = TcpSocket::UnknownSocketError;
                m_errorString = "Windows socket functionality not properly initialized";
                break;
            case WSAEISCONN:
                m_socketState = TcpSocket::ConnectedState; // socket already connected
                break;
            case WSAECONNREFUSED:
            case WSAEINVAL:
                m_socketError = TcpSocket::ConnectionRefusedError;
                m_errorString = "connection refused";
                break;
            case WSAETIMEDOUT:
                m_socketError = TcpSocket::NetworkError;
                m_errorString = "connection timed out";
                break;
            case WSAEHOSTUNREACH:
                m_socketError = TcpSocket::NetworkError;
                m_errorString = "host unreachable";
                break;
            case WSAENETUNREACH:
                m_socketError = TcpSocket::NetworkError;
                m_errorString = "network unreachable";
                break;
            case WSAEADDRINUSE:
                m_socketError = TcpSocket::SocketResourceError;
                m_errorString = "address already in use";
                break;
            case WSAEACCES:
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
    const int protocolNum = ( (protocol == HostAddress::IPv6Protocol) ? AF_INET6 : AF_INET );

    // attempt to create socket
    SOCKET socketFd = WSASocket(protocolNum, SOCK_STREAM, IPPROTO_TCP, 0, 0, WSA_FLAG_OVERLAPPED);

    // if we fetched an invalid socket descriptor
    if ( socketFd == INVALID_SOCKET ) {

        // set error type/message depending on error code
        const int errorCode = WSAGetLastError();
        switch ( errorCode ) {
            case WSANOTINITIALISED:
                m_socketError = TcpSocket::UnknownSocketError;
                m_errorString = "Windows socket functionality not properly initialized";
                break;
            case WSAEAFNOSUPPORT:
            case WSAESOCKTNOSUPPORT:
            case WSAEPROTOTYPE:
            case WSAEINVAL:
                m_socketError = TcpSocket::UnsupportedSocketOperationError;
                m_errorString = "protocol not supported";
                break;
            case WSAEMFILE:
            case WSAENOBUFS:
                m_socketError = TcpSocket::SocketResourceError;
                m_errorString = "out of resources";
                break;
            default:
                m_socketError = TcpSocket::UnknownSocketError;
                stringstream errStream("");
                errStream << "WSA ErrorCode: " << errorCode;
                m_errorString = errStream.str();
                break;
        }

        // return failure
        return false;
    }

    // otherwise, store our socket FD & return success
    m_socketDescriptor = static_cast<int>(socketFd);
    return true;
}

int64_t TcpSocketEngine::nativeNumBytesAvailable(void) const {

    int64_t numBytes(0);
    int64_t dummy(0);
    DWORD bytesWritten(0);

    const int ioctlResult = WSAIoctl( m_socketDescriptor, FIONREAD
                                    , &dummy, sizeof(dummy)
                                    , &numBytes, sizeof(numBytes)
                                    , &bytesWritten, 0, 0
                                    );
    return ( ioctlResult == SOCKET_ERROR ? -1 : numBytes );
}

int64_t TcpSocketEngine::nativeRead(char* dest, size_t max) {

    // skip if invalid socket
    if ( !IsValid() )
        return -1;

    // set up our WSA output buffer
    WSABUF buf;
    buf.buf = dest;
    buf.len = max;

    // attempt to read bytes
    DWORD flags = 0;
    DWORD bytesRead = 0;
    const int readResult = WSARecv(m_socketDescriptor, &buf, 1, &bytesRead, &flags, 0, 0);
    if ( readResult == SOCKET_ERROR )
        return -1;

    // return number of bytes read
    return static_cast<int64_t>(bytesRead);
}

// negative value for msecs will block (forever) until
int TcpSocketEngine::nativeSelect(int msecs, bool isRead) const {

    fd_set fds;
    FD_ZERO(&fds);
    FD_SET(m_socketDescriptor, &fds);

    timeval tv;
    tv.tv_sec  = msecs / 1000;
    tv.tv_usec = (msecs % 1000) * 1000;

    // do 'select'
    if ( isRead )
        return select(0, &fds, 0, 0, (msecs < 0 ? 0 : &tv));
    else
        return select(0, 0, &fds, 0, (msecs < 0 ? 0 : &tv));
}

int64_t TcpSocketEngine::nativeWrite(const char* data, size_t length) {

    // setup our WSA write buffer
    WSABUF buf;
    buf.buf = (char*)data;
    buf.len = length;

    // attempt to write bytes
    DWORD flags = 0;
    DWORD bytesWritten = 0;
    const int writeResult = WSASend(m_socketDescriptor, &buf, 1, &bytesWritten, flags, 0, 0);
    if ( writeResult == SOCKET_ERROR )
        return -1;

    // return number of bytes written
    return static_cast<int64_t>(bytesWritten);
}
