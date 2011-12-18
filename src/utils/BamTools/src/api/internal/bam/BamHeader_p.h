// ***************************************************************************
// BamHeader_p.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for handling BAM headers.
// ***************************************************************************

#ifndef BAMHEADER_P_H
#define BAMHEADER_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include "api/SamHeader.h"
#include <string>

namespace BamTools {
namespace Internal {

class BgzfStream;

class BamHeader {

    // ctor & dtor
    public:
        BamHeader(void);
        ~BamHeader(void);

    // BamHeader interface
    public:
        // clear SamHeader data
        void Clear(void);
        // return true if SamHeader data is valid
        bool IsValid(void) const;
        // load BAM header ('magic number' and SAM header text) from BGZF stream
        // returns true if all OK
        void Load(BgzfStream* stream);
        // returns (editable) copy of SamHeader data object
        SamHeader ToSamHeader(void) const;
        // returns SAM-formatted string of header data
        std::string ToString(void) const;

    // internal methods
    private:
        // reads magic number from BGZF stream
        void CheckMagicNumber(BgzfStream* stream);
        // reads SAM header length from BGZF stream, stores it in @length
        void ReadHeaderLength(BgzfStream* stream, uint32_t& length);
        // reads SAM header text from BGZF stream, stores in SamHeader object
        void ReadHeaderText(BgzfStream* stream, const uint32_t& length);

    // data members
    private:
        SamHeader m_header;
};

} // namespace Internal
} // namespace BamTools

#endif // BAMHEADER_P_H
