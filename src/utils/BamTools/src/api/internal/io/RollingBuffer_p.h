// ***************************************************************************
// RollingBuffer_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 7 December 2011 (DB)
// ---------------------------------------------------------------------------
// Provides a dynamic I/O FIFO byte queue, which removes bytes as they are
// read from the front of the buffer and grows to accept bytes being written
// to buffer end.
//
// implementation note: basically a 'smart' wrapper around 1..* ByteArrays
// ***************************************************************************

#ifndef ROLLINGBUFFER_P_H
#define ROLLINGBUFFER_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include "api/api_global.h"
#include "api/internal/io/ByteArray_p.h"
#include <deque>
#include <string>

namespace BamTools {
namespace Internal {

class RollingBuffer {

    // ctors & dtor
    public:
        RollingBuffer(size_t growth);
        ~RollingBuffer(void);

    // RollingBuffer interface
    public:

        // returns current buffer size
        size_t BlockSize(void) const;
        // checks buffer for new line
        bool CanReadLine(void) const;
        // frees @n bytes from end of buffer
        void Chop(size_t n);
        // clears entire buffer structure
        void Clear(void);
        // frees @n bytes from front of buffer
        void Free(size_t n);
        // checks buffer for @c
        size_t IndexOf(char c) const;
        // returns whether buffer contains data
        bool IsEmpty(void) const;
        // reads up to @maxLen bytes into @dest
        // returns exactly how many bytes were read from buffer
        size_t Read(char* dest, size_t max);
        // reads until newline (or up to @maxLen bytes)
        // returns exactly how many bytes were read from buffer
        size_t ReadLine(char* dest, size_t max);
        // returns a C-fxn compatible char* to byte data
        const char* ReadPointer(void) const;
        // ensures that buffer contains space for @n incoming bytes, returns write-able char*
        char* Reserve(size_t n);
        // returns current number of bytes stored in buffer
        size_t Size(void) const;
        // reserves space for @n bytes, then appends contents of @src to buffer
        void Write(const char* src, size_t n);

    // data members
    private:
        size_t m_head;                // index into current data (next char)
        size_t m_tail;                // index into last data position
        size_t m_tailBufferIndex;     // m_data::size() - 1
        size_t m_totalBufferSize;     // total buffer size
        size_t m_bufferGrowth;        // new buffers are typically initialized with this size
        std::deque<ByteArray> m_data; // basic 'buffer of buffers'
};

} // namespace Internal
} // namespace BamTools

#endif // ROLLINGBUFFER_P_H
