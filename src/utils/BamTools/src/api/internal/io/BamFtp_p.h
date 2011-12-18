// ***************************************************************************
// BamFtp_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 November 2011 (DB)
// ---------------------------------------------------------------------------
// Provides reading/writing of BAM files on FTP server
// ***************************************************************************

#ifndef BAMFTP_P_H
#define BAMFTP_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include "api/IBamIODevice.h"
#include <string>

namespace BamTools {
namespace Internal {

class TcpSocket;

class BamFtp : public IBamIODevice {

    // ctor & dtor
    public:
        BamFtp(const std::string& url);
        ~BamFtp(void);

    // IBamIODevice implementation
    public:
        void Close(void);
        bool IsOpen(void) const;
        bool IsRandomAccess(void) const;
        bool Open(const IBamIODevice::OpenMode mode);
        int64_t Read(char* data, const unsigned int numBytes);
        bool Seek(const int64_t& position, const int origin = SEEK_SET);
        int64_t Tell(void) const;
        int64_t Write(const char* data, const unsigned int numBytes);

    // internal methods
    private:
        bool ConnectCommandSocket(void);
        bool ConnectDataSocket(void);        
        bool ParsePassiveResponse(void);
        void ParseUrl(const std::string& url);
        int64_t ReadCommandSocket(char* data, const unsigned int numBytes);
        int64_t ReadDataSocket(char* data, const unsigned int numBytes);
        bool ReceiveReply(void);
        bool SendCommand(const std::string& command, bool waitForReply);
        int64_t WriteCommandSocket(const char* data, const unsigned int numBytes);
        int64_t WriteDataSocket(const char* data, const unsigned int numBytes);

    // data members
    private:

        // our main sockets
        TcpSocket* m_commandSocket;
        TcpSocket* m_dataSocket;

        // our connection data
        std::string m_hostname;
        uint16_t    m_port;
        std::string m_dataHostname;
        uint16_t    m_dataPort;
        std::string m_filename;

        std::string m_username;
        std::string m_password;

        std::string m_response;

        // internal state flags
        bool m_isUrlParsed;

        // file position
        int64_t m_filePosition;
};

} // namespace Internal
} // namespace BamTools

#endif // BAMFTP_P_H
