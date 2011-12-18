// ***************************************************************************
// SamFormatParser.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 23 December 2010 (DB)
// ---------------------------------------------------------------------------
// Provides functionality for parsing SAM header text into SamHeader object
// ***************************************************************************

#ifndef SAM_FORMAT_PARSER_H
#define SAM_FORMAT_PARSER_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include <string>
#include <vector>

namespace BamTools {

class SamHeader;

namespace Internal {

class SamFormatParser {

    // ctor & dtor
    public:
        SamFormatParser(BamTools::SamHeader& header);
        ~SamFormatParser(void);

    // parse text & populate header data
    public:
        void Parse(const std::string& headerText);

    // internal methods
    private:
        void ParseSamLine(const std::string& line);
        void ParseHDLine(const std::string& line);
        void ParseSQLine(const std::string& line);
        void ParseRGLine(const std::string& line);
        void ParsePGLine(const std::string& line);
        void ParseCOLine(const std::string& line);
        const std::vector<std::string> Split(const std::string& line, const char delim);

    // data members
    private:
        SamHeader& m_header;
};

} // namespace Internal
} // namespace BamTools

#endif // SAM_FORMAT_PARSER_H
