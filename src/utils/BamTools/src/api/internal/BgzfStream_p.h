// ***************************************************************************
// BgzfStream_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 5 April 2011(DB)
// ---------------------------------------------------------------------------
// Based on BGZF routines developed at the Broad Institute.
// Provides the basic functionality for reading & writing BGZF files
// Replaces the old BGZF.* files to avoid clashing with other toolkits
// ***************************************************************************

#ifndef BGZFSTREAM_P_H
#define BGZFSTREAM_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include <api/BamAux.h>
#include <api/BamConstants.h>
#include "zlib.h"
#include <cstdio>
#include <string>

namespace BamTools {
namespace Internal {

class BgzfStream {

    // constructor & destructor
    public:
        BgzfStream(void);
        ~BgzfStream(void);

    // main interface methods
    public:
        // closes BGZF file
        void Close(void);
        // opens the BGZF file (mode is either "rb" for reading, or "wb" for writing)
        bool Open(const std::string& filename, const char* mode);
        // reads BGZF data into a byte buffer
        int Read(char* data, const unsigned int dataLength);
        // seek to position in BGZF file
        bool Seek(const int64_t& position);
        // enable/disable compressed output
        void SetWriteCompressed(bool ok);
        // get file position in BGZF file
        int64_t Tell(void) const;
        // writes the supplied data into the BGZF buffer
        unsigned int Write(const char* data, const unsigned int dataLen);

    // internal methods
    private:
        // compresses the current block
        int DeflateBlock(void);
        // flushes the data in the BGZF block
        void FlushBlock(void);
        // de-compresses the current block
        int InflateBlock(const int& blockLength);
        // reads a BGZF block
        bool ReadBlock(void);

    // static 'utility' methods
    public:
        // checks BGZF block header
        static inline bool CheckBlockHeader(char* header);

    // data members
    public:
        unsigned int UncompressedBlockSize;
        unsigned int CompressedBlockSize;
        unsigned int BlockLength;
        unsigned int BlockOffset;
        uint64_t BlockAddress;
        bool IsOpen;
        bool IsWriteOnly;
        bool IsWriteCompressed;
        FILE* Stream;
        char* UncompressedBlock;
        char* CompressedBlock;
};

// -------------------------------------------------------------
// static 'utility' method implementations

// checks BGZF block header
inline
bool BgzfStream::CheckBlockHeader(char* header) {
    return (header[0] == Constants::GZIP_ID1 &&
            header[1] == (char)Constants::GZIP_ID2 &&
            header[2] == Z_DEFLATED &&
            (header[3] & Constants::FLG_FEXTRA) != 0 &&
            BamTools::UnpackUnsignedShort(&header[10]) == Constants::BGZF_XLEN &&
            header[12] == Constants::BGZF_ID1 &&
            header[13] == Constants::BGZF_ID2 &&
            BamTools::UnpackUnsignedShort(&header[14]) == Constants::BGZF_LEN );
}

} // namespace Internal
} // namespace BamTools

#endif // BGZFSTREAM_P_H
