// ***************************************************************************
// TcpSocketEngine_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 8 December 2011 (DB)
// ---------------------------------------------------------------------------
// Provides low-level implementation of TCP I/O
// ***************************************************************************

#ifndef TCPSOCKETENGINE_P_H
#define TCPSOCKETENGINE_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include "api/internal/io/HostAddress_p.h"
#include "api/internal/io/TcpSocket_p.h"

#ifdef _WIN32
#  include "api/internal/io/NetWin_p.h"
#endif

namespace BamTools {
namespace Internal {

struct TcpSocketEngine {

    // ctors & dtor
    public:
        TcpSocketEngine(void);
        TcpSocketEngine(const TcpSocketEngine& other);
        ~TcpSocketEngine(void);

    // TcpSocketEngine interface
    public:

        // connection-related methods
        void Close(void);
        bool Connect(const HostAddress& address, const uint16_t port);
        bool Initialize(HostAddress::NetworkProtocol protocol);
        bool IsValid(void) const;

        // IO-related methods
        int64_t NumBytesAvailable(void) const;
        int64_t Read(char* dest, size_t max);
        int64_t Write(const char* data, size_t length);

        bool WaitForRead(int msec, bool* timedOut);
        bool WaitForWrite(int msec, bool* timedOut);

        // query connection state
//        HostAddress GetLocalAddress(void) const;
//        uint16_t GetLocalPort(void) const;
        HostAddress GetRemoteAddress(void) const;
        uint16_t    GetRemotePort(void) const;

        int GetSocketDescriptor(void) const;
        TcpSocket::SocketError GetSocketError(void);
        TcpSocket::SocketState GetSocketState(void);

        std::string GetErrorString(void) const;

    // platform-dependent internal methods
    // provided in the corresponding TcpSocketEngine_<OS>_p.cpp
    private:
        void    nativeClose(void);
        bool    nativeConnect(const HostAddress& address, const uint16_t port);
        bool    nativeCreateSocket(HostAddress::NetworkProtocol protocol);
        void    nativeDisconnect(void);
        int64_t nativeNumBytesAvailable(void) const;
        int64_t nativeRead(char* dest, size_t max);
        int     nativeSelect(int msecs, bool isRead) const;
        int64_t nativeWrite(const char* data, size_t length);

    // data members
    private:
        int m_socketDescriptor;

//        HostAddress m_localAddress;
        HostAddress m_remoteAddress;
//        uint16_t m_localPort;
        uint16_t m_remotePort;

        TcpSocket::SocketError m_socketError;
        TcpSocket::SocketState m_socketState;
        std::string m_errorString;

#ifdef _WIN32
        WindowsSockInit m_win;
#endif
};

} // namespace Internal
} // namespace BamTools

#endif // TCPSOCKETENGINE_P_H
