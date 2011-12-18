// ***************************************************************************
// SamFormatPrinter.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 6 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides functionality for printing formatted SAM header to string
// ***************************************************************************

#ifndef SAM_FORMAT_PRINTER_H
#define SAM_FORMAT_PRINTER_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include <sstream>
#include <string>

namespace BamTools {

class SamHeader;

namespace Internal {

class SamFormatPrinter {

    // ctor & dtor
    public:
        SamFormatPrinter(const BamTools::SamHeader& header);
        ~SamFormatPrinter(void);

    // generates SAM-formatted string from header data
    public:
        const std::string ToString(void) const;

    // internal methods
    private:
        void PrintHD(std::stringstream& out) const;
        void PrintSQ(std::stringstream& out) const;
        void PrintRG(std::stringstream& out) const;
        void PrintPG(std::stringstream& out) const;
        void PrintCO(std::stringstream& out) const;

    // data members
    private:
        const SamHeader& m_header;
};

} // namespace Internal
} // namespace BamTools

#endif // SAM_FORMAT_PRINTER_H
