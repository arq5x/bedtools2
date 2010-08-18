// ***************************************************************************
// BGZF.h (c) 2009 Derek Barnett, Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 16 August 2010 (DB)
// ---------------------------------------------------------------------------
// BGZF routines were adapted from the bgzf.c code developed at the Broad
// Institute.
// ---------------------------------------------------------------------------
// Provides the basic functionality for reading & writing BGZF files
// ***************************************************************************

#ifndef BGZF_H
#define BGZF_H

// 'C' includes
#include <cstdio>
#include <cstdlib>
#include <cstring>

// C++ includes
#include <string>

// zlib includes
#include "zlib.h"

// Platform-specific large-file support
#ifndef BAMTOOLS_LFS
#define BAMTOOLS_LFS
    #ifdef WIN32
        #define ftell64(a)     _ftelli64(a)
        #define fseek64(a,b,c) _fseeki64(a,b,c)
    #else
        #define ftell64(a)     ftello(a)
        #define fseek64(a,b,c) fseeko(a,b,c) 
    #endif
#endif // BAMTOOLS_LFS

// Platform-specific type definitions
#ifndef BAMTOOLS_TYPES
#define BAMTOOLS_TYPES
    #ifdef _MSC_VER
        typedef char                 int8_t;
        typedef unsigned char       uint8_t;
        typedef short               int16_t;
        typedef unsigned short     uint16_t;
        typedef int                 int32_t;
        typedef unsigned int       uint32_t;
        typedef long long           int64_t;
        typedef unsigned long long uint64_t;
    #else    
        #include <stdint.h>
    #endif
#endif // BAMTOOLS_TYPES

namespace BamTools {

// zlib constants
const int GZIP_ID1   = 31;
const int GZIP_ID2   = 139;
const int CM_DEFLATE = 8;
const int FLG_FEXTRA = 4;
const int OS_UNKNOWN = 255;
const int BGZF_XLEN  = 6;
const int BGZF_ID1   = 66;
const int BGZF_ID2   = 67;
const int BGZF_LEN   = 2;
const int GZIP_WINDOW_BITS    = -15;
const int Z_DEFAULT_MEM_LEVEL = 8;

// BZGF constants
const int BLOCK_HEADER_LENGTH = 18;
const int BLOCK_FOOTER_LENGTH = 8;
const int MAX_BLOCK_SIZE      = 65536;
const int DEFAULT_BLOCK_SIZE  = 65536;

struct BgzfData {

    // data members
    public:
        unsigned int UncompressedBlockSize;
        unsigned int CompressedBlockSize;
        unsigned int BlockLength;
        unsigned int BlockOffset;
        uint64_t BlockAddress;
        bool     IsOpen;
        bool     IsWriteOnly;
        bool     IsWriteUncompressed;
        FILE*    Stream;
        char*    UncompressedBlock;
        char*    CompressedBlock;

    // constructor & destructor
    public:
        BgzfData(void);
        ~BgzfData(void);

    // main interface methods
    public:       
        // closes BGZF file
        void Close(void);
        // opens the BGZF file (mode is either "rb" for reading, or "wb" for writing)
        bool Open(const std::string& filename, const char* mode, bool isWriteUncompressed = false);
        // reads BGZF data into a byte buffer
        int Read(char* data, const unsigned int dataLength);
        // seek to position in BGZF file
        bool Seek(int64_t position);
        // get file position in BGZF file
        int64_t Tell(void);
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
        // packs an unsigned integer into the specified buffer
        static inline void PackUnsignedInt(char* buffer, unsigned int value);
        // packs an unsigned short into the specified buffer
        static inline void PackUnsignedShort(char* buffer, unsigned short value);
        // unpacks a buffer into a double
        static inline double UnpackDouble(char* buffer);
        static inline double UnpackDouble(const char* buffer);
        // unpacks a buffer into a float
        static inline float UnpackFloat(char* buffer);
        static inline float UnpackFloat(const char* buffer);
        // unpacks a buffer into a signed int
        static inline signed int UnpackSignedInt(char* buffer);
        static inline signed int UnpackSignedInt(const char* buffer);
        // unpacks a buffer into a signed short
        static inline signed short UnpackSignedShort(char* buffer);
        static inline signed short UnpackSignedShort(const char* buffer);
        // unpacks a buffer into an unsigned int
        static inline unsigned int UnpackUnsignedInt(char* buffer);
        static inline unsigned int UnpackUnsignedInt(const char* buffer);
        // unpacks a buffer into an unsigned short
        static inline unsigned short UnpackUnsignedShort(char* buffer);
        static inline unsigned short UnpackUnsignedShort(const char* buffer);
};

// -------------------------------------------------------------
// static 'utility' method implementations

// checks BGZF block header
inline
bool BgzfData::CheckBlockHeader(char* header) {
    return (header[0] == GZIP_ID1 &&
            header[1] == (char)GZIP_ID2 &&
            header[2] == Z_DEFLATED &&
            (header[3] & FLG_FEXTRA) != 0 &&
            BgzfData::UnpackUnsignedShort(&header[10]) == BGZF_XLEN &&
            header[12] == BGZF_ID1 &&
            header[13] == BGZF_ID2 &&
            BgzfData::UnpackUnsignedShort(&header[14]) == BGZF_LEN );
}

// 'packs' an unsigned integer into the specified buffer
inline
void BgzfData::PackUnsignedInt(char* buffer, unsigned int value) {
    buffer[0] = (char)value;
    buffer[1] = (char)(value >> 8);
    buffer[2] = (char)(value >> 16);
    buffer[3] = (char)(value >> 24);
}

// 'packs' an unsigned short into the specified buffer
inline
void BgzfData::PackUnsignedShort(char* buffer, unsigned short value) {
    buffer[0] = (char)value;
    buffer[1] = (char)(value >> 8);
}

// 'unpacks' a buffer into a double (includes both non-const & const char* flavors)
inline
double BgzfData::UnpackDouble(char* buffer) {
    union { double value; unsigned char valueBuffer[sizeof(double)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    un.valueBuffer[4] = buffer[4];
    un.valueBuffer[5] = buffer[5];
    un.valueBuffer[6] = buffer[6];
    un.valueBuffer[7] = buffer[7];
    return un.value;
}

inline
double BgzfData::UnpackDouble(const char* buffer) {
    union { double value; unsigned char valueBuffer[sizeof(double)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    un.valueBuffer[4] = buffer[4];
    un.valueBuffer[5] = buffer[5];
    un.valueBuffer[6] = buffer[6];
    un.valueBuffer[7] = buffer[7];
    return un.value;
}

// 'unpacks' a buffer into a float (includes both non-const & const char* flavors)
inline
float BgzfData::UnpackFloat(char* buffer) {
    union { float value; unsigned char valueBuffer[sizeof(float)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    return un.value;
}

inline
float BgzfData::UnpackFloat(const char* buffer) {
    union { float value; unsigned char valueBuffer[sizeof(float)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    return un.value;
}

// 'unpacks' a buffer into a signed int (includes both non-const & const char* flavors)
inline
signed int BgzfData::UnpackSignedInt(char* buffer) {
    union { signed int value; unsigned char valueBuffer[sizeof(signed int)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    return un.value;
}

inline
signed int BgzfData::UnpackSignedInt(const char* buffer) {
    union { signed int value; unsigned char valueBuffer[sizeof(signed int)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    return un.value;
}

// 'unpacks' a buffer into a signed short (includes both non-const & const char* flavors)
inline
signed short BgzfData::UnpackSignedShort(char* buffer) {
    union { signed short value; unsigned char valueBuffer[sizeof(signed short)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    return un.value;
}

inline
signed short BgzfData::UnpackSignedShort(const char* buffer) {
    union { signed short value; unsigned char valueBuffer[sizeof(signed short)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    return un.value;
}

// 'unpacks' a buffer into an unsigned int (includes both non-const & const char* flavors)
inline
unsigned int BgzfData::UnpackUnsignedInt(char* buffer) {
    union { unsigned int value; unsigned char valueBuffer[sizeof(unsigned int)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    return un.value;
}

inline
unsigned int BgzfData::UnpackUnsignedInt(const char* buffer) {
    union { unsigned int value; unsigned char valueBuffer[sizeof(unsigned int)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    return un.value;
}

// 'unpacks' a buffer into an unsigned short (includes both non-const & const char* flavors)
inline
unsigned short BgzfData::UnpackUnsignedShort(char* buffer) {
    union { unsigned short value; unsigned char valueBuffer[sizeof(unsigned short)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    return un.value;
}

inline
unsigned short BgzfData::UnpackUnsignedShort(const char* buffer) {
    union { unsigned short value; unsigned char valueBuffer[sizeof(unsigned short)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    return un.value;
}

} // namespace BamTools

#endif // BGZF_H
