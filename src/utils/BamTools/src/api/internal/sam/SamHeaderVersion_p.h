// ***************************************************************************
// SamHeaderVersion.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides functionality for comparing SAM header versions
// *************************************************************************

#ifndef SAM_HEADERVERSION_P_H
#define SAM_HEADERVERSION_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include "api/SamConstants.h"
#include <sstream>
#include <string>

namespace BamTools {
namespace Internal {

class SamHeaderVersion {

    // ctors & dtor
    public:
        SamHeaderVersion(void)
            : m_majorVersion(0)
            , m_minorVersion(0)
        { }

        explicit SamHeaderVersion(const std::string& version)
            : m_majorVersion(0)
            , m_minorVersion(0)
        {
            SetVersion(version);
        }

        SamHeaderVersion(const unsigned int& major, const unsigned int& minor)
            : m_majorVersion(major)
            , m_minorVersion(minor)
        { }

        ~SamHeaderVersion(void) {
            m_majorVersion = 0;
            m_minorVersion = 0;
        }
    
    // acess data
    public:
        unsigned int MajorVersion(void) const { return m_majorVersion; }
        unsigned int MinorVersion(void) const { return m_minorVersion; }

        void SetVersion(const std::string& version);
        std::string ToString(void) const;

    // data members
    private:
        unsigned int m_majorVersion;
        unsigned int m_minorVersion;
};

inline
void SamHeaderVersion::SetVersion(const std::string& version) {

    // do nothing if version is empty
    if ( !version.empty() ) {

        std::stringstream versionStream("");

        // do nothing if period not found
        const size_t periodFound = version.find(Constants::SAM_PERIOD);
        if ( periodFound != std::string::npos ) {

            // store major version if non-empty and contains only digits
            const std::string& majorVersion = version.substr(0, periodFound);
            versionStream.str(majorVersion);
            if ( !majorVersion.empty() ) {
                const size_t nonDigitFound = majorVersion.find_first_not_of(Constants::SAM_DIGITS);
                if ( nonDigitFound == std::string::npos )
                    versionStream >> m_majorVersion;
            }

            // store minor version if non-empty and contains only digits
            const std::string& minorVersion = version.substr(periodFound + 1);
            versionStream.str(minorVersion);
            if ( !minorVersion.empty() ) {
                const size_t nonDigitFound = minorVersion.find_first_not_of(Constants::SAM_DIGITS);
                if ( nonDigitFound == std::string::npos )
                    versionStream >> m_minorVersion;
            }
        }
    }
}

// -----------------------------------------------------
// printing

inline std::string SamHeaderVersion::ToString(void) const {
    std::stringstream version;
    version << m_majorVersion << Constants::SAM_PERIOD << m_minorVersion;
    return version.str();
}

// -----------------------------------------------------
// comparison operators

inline bool operator==(const SamHeaderVersion& lhs, const SamHeaderVersion& rhs) {
    return (lhs.MajorVersion() == rhs.MajorVersion()) &&
           (lhs.MinorVersion() == rhs.MinorVersion());
}

inline bool operator<(const SamHeaderVersion& lhs, const SamHeaderVersion& rhs) {
    if ( lhs.MajorVersion() == rhs.MajorVersion() )
        return lhs.MinorVersion() < rhs.MinorVersion();
    else 
        return lhs.MajorVersion() < rhs.MajorVersion();
}

inline bool operator> (const SamHeaderVersion& lhs, const SamHeaderVersion& rhs) { return rhs < lhs;  }
inline bool operator<=(const SamHeaderVersion& lhs, const SamHeaderVersion& rhs) { return !(lhs>rhs); }
inline bool operator>=(const SamHeaderVersion& lhs, const SamHeaderVersion& rhs) { return !(lhs<rhs); }

} // namespace Internal 
} // namespace BamTools

#endif // SAM_HEADERVERSION_P_H
