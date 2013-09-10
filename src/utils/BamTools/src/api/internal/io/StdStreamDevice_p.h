/*
 * StdStreamDevice_p.h
 *
 *  Created on: Apr 24, 2013
 *      Author: nek3d
 */

#ifndef STDSTREAMDEVICE_P_H_
#define STDSTREAMDEVICE_P_H_

//NOTE: The stream management methods in this class are placeholders. The class is not intended to manage a
// a stream. Rather, it is intended to accept a stream that is managed elsewhere, and allow that stream to
//be used with bamtools.

#include "api/IBamIODevice.h"
#include <cstring> //for strlen
#include <istream>

namespace BamTools {
namespace Internal {

class StdStreamDevice : public IBamIODevice {

    // ctor & dtor
    public:
        StdStreamDevice(std::istream *stream) : m_stream(stream) {
        	IBamIODevice::m_mode = IBamIODevice::ReadOnly;
        }
        ~StdStreamDevice() {}

    // IBamIODevice implementation
    public:
        void Close(void) {}
        bool IsOpen(void) const { return true; } // assume stream passed in is always open.
        bool IsRandomAccess(void) const { return false; }
        bool Open(const IBamIODevice::OpenMode mode) { return true; } //stream is assumed opened by caller.
        int64_t Read(char* data, const unsigned int numBytes) {
        	m_stream->read(data, numBytes);
        	return (int64_t)(numBytes);
        }
        bool Seek(const int64_t& position, const int origin = SEEK_SET) { return false; }

        //Not sure what Tell and Write should do for streams.
        int64_t Tell(void) const { return 0; }
        int64_t Write(const char* data, const unsigned int numBytes) { return 0; }

    private:
        std::istream* m_stream;
};

} // namespace Internal
} // namespace BamTools

#endif /* STDSTREAMDEVICE_P_H_ */
