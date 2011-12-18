// ***************************************************************************
// ILocalIODevice_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides shared behavior for files & pipes
// ***************************************************************************

#ifndef ILOCALIODEVICE_P_H
#define ILOCALIODEVICE_P_H

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

namespace BamTools {
namespace Internal {

class ILocalIODevice : public IBamIODevice {

    // ctor & dtor
    public:
        ILocalIODevice(void);
        virtual ~ILocalIODevice(void);

    // IBamIODevice implementation
    public:
        virtual void Close(void);
        virtual int64_t Read(char* data, const unsigned int numBytes);
        virtual int64_t Tell(void) const;
        virtual int64_t Write(const char* data, const unsigned int numBytes);

    // data members
    protected:
        FILE* m_stream;
};

} // namespace Internal
} // namespace BamTools

#endif // ILOCALIODEVICE_P_H
