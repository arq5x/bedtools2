// ***************************************************************************
// BamHttp_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 November 2011 (DB)
// ---------------------------------------------------------------------------
// Provides reading/writing of BAM files on HTTP server
// ***************************************************************************

#ifndef BAMHTTP_P_H
#define BAMHTTP_P_H

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

class HttpRequestHeader;
class HttpResponseHeader;
class TcpSocket;

class BamHttp : public IBamIODevice {

    // ctor & dtor
    public:
        BamHttp(const std::string& url);
        ~BamHttp(void);

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
        bool ConnectSocket(void);
        bool EnsureSocketConnection(void);
        void ParseUrl(const std::string& url);
        int64_t ReadFromSocket(char* data, const unsigned int numBytes);
        bool ReceiveResponse(void);
        bool SendRequest(const size_t numBytes = 0);
        int64_t WriteToSocket(const char* data, const unsigned int numBytes);

    // data members
    private:

        // our main socket
        TcpSocket* m_socket;

        // our connection data
        std::string m_hostname;
        std::string m_port;
        std::string m_filename;

        // our last (active) request & response info
        HttpRequestHeader*  m_request;
        HttpResponseHeader* m_response;

        // internal state flags
        bool m_isUrlParsed;

        // file position
        int64_t m_filePosition;
        int64_t m_endRangeFilePosition;
};

} // namespace Internal
} // namespace BamTools

#endif // BAMHTTP_P_H
