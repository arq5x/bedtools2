// ***************************************************************************
// RollingBuffer_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 8 December 2011 (DB)
// ---------------------------------------------------------------------------
// Provides a dynamic I/O FIFO byte queue, which removes bytes as they are
// read from the front of the buffer and grows to accept bytes being written
// to buffer end.
//
// implementation note: basically a 'smart' wrapper around 1..* ByteArrays
// ***************************************************************************

#include "api/internal/io/RollingBuffer_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <climits>
#include <cstring>
#include <algorithm>
#include <string>
using namespace std;

// ------------------------------
// RollingBuffer implementation
// ------------------------------

RollingBuffer::RollingBuffer(size_t growth)
    : m_bufferGrowth(growth)
{
    // buffer always contains at least 1 (maybe empty) byte array
    m_data.push_back( ByteArray() );

    // set cleared state
    Clear();
}

RollingBuffer::~RollingBuffer(void) { }

size_t RollingBuffer::BlockSize(void) const {

    // if only one byte array in buffer <- needed?
    if ( m_tailBufferIndex == 0 )
        return m_tail - m_head;

    // otherwise return remaining num bytes in first array
    const ByteArray& first = m_data.front();
    return first.Size() - m_head;
}

bool RollingBuffer::CanReadLine(void) const {
    return IndexOf('\n') != string::npos;
}

void RollingBuffer::Chop(size_t n) {

    // update buffer size
    if ( n > m_totalBufferSize )
        m_totalBufferSize = 0;
    else
        m_totalBufferSize -= n;

    // loop until target case hit
    for ( ; ; ) {

        // if only one array, decrement tail
        if ( m_tailBufferIndex == 0 ) {
            m_tail -= n;

            // if all data chopped
            if ( m_tail <= m_head ) {
                m_head = 0;
                m_tail = 0;
            }
            return;
        }

        // if there's room in last byte array to 'chop', just decrement tail
        if ( n <= m_tail ) {
            m_tail -= n;
            return;
        }

        // otherwise we're going to overlap our internal byte arrays
        // reduce our chop amount by the amount of data in the last byte array
        n -= m_tail;

        // remove last byte array & set tail to it's end
        m_data.pop_back();
        --m_tailBufferIndex;
        m_tail = m_data.at(m_tailBufferIndex).Size();
    }

    // if buffer is now empty, reset state & clear up memory
    if ( IsEmpty() )
        Clear();
}

void RollingBuffer::Clear(void) {

    // remove all byte arrays (except first)
    m_data.erase( m_data.begin()+1, m_data.end() );

    // clear out first byte array
    m_data[0].Resize(0);
    m_data[0].Squeeze();

    // reset index & size markers
    m_head = 0;
    m_tail = 0;
    m_tailBufferIndex = 0;
    m_totalBufferSize = 0;
}

void RollingBuffer::Free(size_t n) {

    // update buffer size
    if ( n > m_totalBufferSize )
        m_totalBufferSize = 0;
    else
        m_totalBufferSize -= n;

    // loop until target case hit
    for ( ; ; ) {

        const size_t blockSize = BlockSize();

        // if there's room in current array
        if ( n < blockSize ) {

            // shift 'head' over @n bytes
            m_head += n;

            // check for emptied, single byte array
            if ( m_head == m_tail && m_tailBufferIndex == 0 ) {
                m_head = 0;
                m_tail = 0;
            }

            break;
        }

        // otherwise we need to check next byte array
        // first update amount to remove
        n -= blockSize;

        // special case - there was only 1 array
        if ( m_data.size() == 1 ) {
            if ( m_data.at(0).Size() != m_bufferGrowth )
                m_data[0].Resize(m_bufferGrowth);
            m_head = 0;
            m_tail = 0;
            m_tailBufferIndex = 0;
            break;
        }

        // otherwise, remove first array and move to next iteration
        m_data.pop_front();
        --m_tailBufferIndex;
        m_head = 0;
    }

    // if buffer is now empty, reset state & clear up memory
    if ( IsEmpty() )
        Clear();
}

size_t RollingBuffer::IndexOf(char c) const {

    // skip processing if empty buffer
    if ( IsEmpty() )
        return string::npos;

    size_t index(0);

    // iterate over byte arrays
    const size_t numBuffers = m_data.size();
    for ( size_t i = 0; i < numBuffers; ++i ) {
        const ByteArray& current = m_data.at(i);

        // if on first array, use head; else 0
        const size_t start = ( (i==0) ? m_head : 0 );

        // if on last array, set end; else use current byte array size
        const size_t end   = ( (i==m_tailBufferIndex) ? m_tail : current.Size());

        // look through this iteration's byte array for @c
        const char* p = current.ConstData()+start;
        for ( size_t j = start; j < end; ++j ) {
            if ( *p++ == c )
                return index;
            ++index;
        }
    }

    // no match found
    return string::npos;
}

bool RollingBuffer::IsEmpty(void) const {
    return (m_tailBufferIndex == 0) && (m_tail == 0);
}

size_t RollingBuffer::Read(char* dest, size_t max) {

    size_t bytesToRead    = std::min(Size(), max);
    size_t bytesReadSoFar = 0;

    while ( bytesReadSoFar < bytesToRead ) {
        const char* readPtr = ReadPointer();
        size_t blockBytes = std::min( (bytesToRead - bytesReadSoFar), BlockSize() );
        if ( dest )
            memcpy(dest+bytesReadSoFar, readPtr, blockBytes);
        bytesReadSoFar += blockBytes;
        Free(blockBytes);
    }

    return bytesReadSoFar;
}

size_t RollingBuffer::ReadLine(char* dest, size_t max) {

    // if we can't read line or if max is 0
    if ( !CanReadLine() || max == 0 )
        return 0;

    // otherwise, read until we hit newline
    size_t bytesReadSoFar = 0;
    bool finished = false;
    while ( !finished ) {

        const size_t index = IndexOf('\n');
        const char* readPtr = ReadPointer();
        size_t bytesToRead = std::min( (index+1)-bytesReadSoFar, BlockSize() );
        bytesToRead = std::min( bytesToRead, (max-1)-bytesReadSoFar );
        memcpy(dest+bytesReadSoFar, readPtr, bytesToRead);
        bytesReadSoFar += bytesToRead;
        Free(bytesToRead);

        if ( !((bytesReadSoFar < index+1)&&(bytesReadSoFar < max-1)) )
            finished = true;
    }

    // null terminate 'dest' & return numBytesRead
    dest[bytesReadSoFar] = '\0';
    return bytesReadSoFar;
}

const char* RollingBuffer::ReadPointer(void) const {

    // return null if empty buffer
    if ( m_data.empty() )
        return 0;

    // otherwise return pointer to current position
    const ByteArray& first = m_data.front();
    return first.ConstData() + m_head;
}

char* RollingBuffer::Reserve(size_t n) {

    // if empty buffer
    if ( m_totalBufferSize == 0 ) {
        m_data[0].Resize( std::max(m_bufferGrowth, n) );
        m_totalBufferSize += n;
        m_tail = n;
        return m_data[m_tailBufferIndex].Data();
    }

    // increment buffer's byte count
    m_totalBufferSize += n;

    // if buffer already contains enough space to fit @n more bytes
    if ( (m_tail + n) <= m_data.at(m_tailBufferIndex).Size() ) {

        // fetch write pointer at current 'tail', increment tail by @n & return
        char* ptr = m_data[m_tailBufferIndex].Data() + m_tail;
        m_tail += n;
        return ptr;
    }

    // if last byte array isn't half full
    if ( m_tail < m_data.at(m_tailBufferIndex).Size()/2 ) {

        // we'll allow simple resize
        m_data[m_tailBufferIndex].Resize(m_tail + n);

        // fetch write pointer at current 'tail', increment tail by @n & return
        char* ptr = m_data[m_tailBufferIndex].Data() + m_tail;
        m_tail += n;
        return ptr;
    }

    // otherwise, shrink last byte array to current used size
    m_data[m_tailBufferIndex].Resize(m_tail);

    // then append new byte array
    m_data.push_back( ByteArray() );
    ++m_tailBufferIndex;
    m_data[m_tailBufferIndex].Resize( std::max(m_bufferGrowth, n) );
    m_tail = n;

    // return write-able pointer on new array
    return m_data[m_tailBufferIndex].Data();
}

size_t RollingBuffer::Size(void) const {
    return m_totalBufferSize;
}

void RollingBuffer::Write(const char* src, size_t n) {
    char* writePtr = Reserve(n);
    memcpy(writePtr, src, n);
}
