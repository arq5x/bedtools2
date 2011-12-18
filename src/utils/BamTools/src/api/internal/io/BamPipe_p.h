// ***************************************************************************
// BamPipe_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 November 2011 (DB)
// ---------------------------------------------------------------------------
// Provides BAM pipe-specific IO behavior
// ***************************************************************************

#ifndef BAMPIPE_P_H
#define BAMPIPE_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include "api/internal/io/ILocalIODevice_p.h"
#include <string>

namespace BamTools {
namespace Internal {

class BamPipe : public ILocalIODevice {

    // ctor & dtor
    public:
        BamPipe(void);
        ~BamPipe(void);

    // IBamIODevice implementation
    public:
        bool IsRandomAccess(void) const;
        bool Open(const IBamIODevice::OpenMode mode);
        bool Seek(const int64_t& position, const int origin = SEEK_SET);
};

} // namespace Internal
} // namespace BamTools

#endif // BAMPIPE_P_H
